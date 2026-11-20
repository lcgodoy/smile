#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

// Planar (R^2) Hausdorff distance between spatial units.
//
// Supports three methods:
// - "discrete": discrete vertex-to-edge Hausdorff distance (matches GEOS).
// - "boundary": continuous Hausdorff distance between boundaries/linework.
// - "polygon":  continuous Hausdorff distance between closed spatial sets
//               (polygons include 2D interiors excluding holes).
//
// Every distance in this file is in the units of the input coordinates.

namespace {

const double INF = std::numeric_limits<double>::infinity();
const double EPS = std::numeric_limits<double>::epsilon();
const double GOLDEN_RATIO = 0.6180339887498949;

inline double dist2_2d(double x1, double y1, double x2, double y2) {
  const double dx = x1 - x2;
  const double dy = y1 - y2;
  return dx * dx + dy * dy;
}

struct PartR2 {
  int v_begin;
  int v_end;
  int e_begin;
  int e_end;
  double xmin, ymin, xmax, ymax;
  bool is_poly;
};

// A spatial unit, prepared once for all its distances.
//
// Input matrices have columns x, y, optionally a part id (see .to_list_mat()
// in R/dist-areal.R), and optionally an is_poly flag (1 for polygon rings, 0
// for points/linestrings; see .to_polygon_list()). Without a part id, all rows
// form one ring, closed if its first and last rows differ.
struct UnitR2 {
  std::vector<double> v;          // 2 * nr coordinates: (x, y)
  std::vector<int> order;         // Golden-ratio vertex traversal order
  std::vector<int> a, b;          // Edge endpoint indices (b == -1 for points)
  std::vector<double> dx, dy;     // Edge vector B - A
  std::vector<double> inv_len2;   // 1 / ||B - A||^2
  std::vector<double> ex0, ey0;   // Edge AABB min
  std::vector<double> ex1, ey1;   // Edge AABB max
  std::vector<uint8_t> e_is_poly; // 1 if edge belongs to a polygon ring
  std::vector<PartR2> parts;      // Per-part / ring metadata
  double xmin, ymin, xmax, ymax;  // Unit AABB
  double diam;                    // Unit AABB diagonal
  double max_abs;                 // Max(|x|, |y|) for ULP floor
  bool has_poly;
  int nv() const { return static_cast<int>(order.size()); }
  int ne() const { return static_cast<int>(a.size()); }
};

void add_edge_r2(UnitR2& u, int i, int j, bool is_poly) {
  const double ax = u.v[2 * i];
  const double ay = u.v[2 * i + 1];
  double dx = 0.0;
  double dy = 0.0;
  double inv_len2 = 0.0;
  double x0 = ax, y0 = ay, x1 = ax, y1 = ay;
  if (j >= 0) {
    const double bx = u.v[2 * j];
    const double by = u.v[2 * j + 1];
    dx = bx - ax;
    dy = by - ay;
    const double len2 = dx * dx + dy * dy;
    if (len2 == 0.0) {
      j = -1;
    } else {
      inv_len2 = 1.0 / len2;
      x0 = std::min(ax, bx);
      y0 = std::min(ay, by);
      x1 = std::max(ax, bx);
      y1 = std::max(ay, by);
    }
  }
  u.a.push_back(i);
  u.b.push_back(j);
  u.dx.push_back(dx);
  u.dy.push_back(dy);
  u.inv_len2.push_back(inv_len2);
  u.ex0.push_back(x0);
  u.ey0.push_back(y0);
  u.ex1.push_back(x1);
  u.ey1.push_back(y1);
  u.e_is_poly.push_back(is_poly ? 1 : 0);
}

// A unit with a non-finite coordinate (for example, the NaN row that
// wk::wk_coords() gives for an EMPTY member of a MULTIPOINT) is empty, so its
// distances are NA, as in src/cr-dist.cpp.
UnitR2 make_unit_r2(const Eigen::Ref<const Eigen::MatrixXd>& poly) {
  UnitR2 u;
  u.xmin = u.ymin = INF;
  u.xmax = u.ymax = -INF;
  u.diam = 0.0;
  u.max_abs = 0.0;
  u.has_poly = false;

  const int nr = static_cast<int>(poly.rows());
  if (nr == 0) return u;
  for (int i = 0; i < nr; i++) {
    if (!std::isfinite(poly(i, 0)) || !std::isfinite(poly(i, 1))) return u;
  }

  u.v.resize(2 * nr);
  for (int i = 0; i < nr; i++) {
    const double x = poly(i, 0);
    const double y = poly(i, 1);
    u.v[2 * i] = x;
    u.v[2 * i + 1] = y;
    u.xmin = std::min(u.xmin, x);
    u.ymin = std::min(u.ymin, y);
    u.xmax = std::max(u.xmax, x);
    u.ymax = std::max(u.ymax, y);
    u.max_abs = std::max(u.max_abs, std::max(std::fabs(x), std::fabs(y)));
  }
  const double bw = u.xmax - u.xmin;
  const double bh = u.ymax - u.ymin;
  u.diam = std::sqrt(bw * bw + bh * bh);

  const bool has_part = poly.cols() > 2;
  const bool has_poly_col = poly.cols() > 3;

  int p_start = 0;
  while (p_start < nr) {
    int p_end = p_start + 1;
    if (!has_part) {
      p_end = nr;
    } else {
      const double pid = poly(p_start, 2);
      while (p_end < nr && poly(p_end, 2) == pid) p_end++;
    }
    const int plen = p_end - p_start;
    const bool part_poly = has_poly_col ? (poly(p_start, 3) != 0.0) : (plen >= 3);
    if (part_poly) u.has_poly = true;

    PartR2 part;
    part.v_begin = p_start;
    part.v_end = p_end;
    part.e_begin = u.ne();
    part.xmin = part.ymin = INF;
    part.xmax = part.ymax = -INF;
    part.is_poly = part_poly;

    for (int i = p_start; i < p_end; i++) {
      part.xmin = std::min(part.xmin, u.v[2 * i]);
      part.ymin = std::min(part.ymin, u.v[2 * i + 1]);
      part.xmax = std::max(part.xmax, u.v[2 * i]);
      part.ymax = std::max(part.ymax, u.v[2 * i + 1]);
    }

    if (plen == 1) {
      add_edge_r2(u, p_start, -1, false);
    } else {
      for (int i = p_start; i + 1 < p_end; i++) {
        add_edge_r2(u, i, i + 1, part_poly);
      }
      // Closing edge if needed (only for input without a part id)
      if (!has_part &&
          (poly(0, 0) != poly(nr - 1, 0) || poly(0, 1) != poly(nr - 1, 1))) {
        add_edge_r2(u, nr - 1, 0, part_poly);
      }
    }
    part.e_end = u.ne();
    u.parts.push_back(part);
    p_start = p_end;
  }

  // A fixed scrambled order of the vertices (golden-ratio sequence). With it,
  // the largest distance grows fast, and the early break in
  // directed_haus_r2() skips more vertices.
  std::vector<std::pair<double, int> > key(nr);
  for (int i = 0; i < nr; i++) {
    key[i] = std::make_pair(std::fmod(static_cast<double>(i + 1) * GOLDEN_RATIO, 1.0), i);
  }
  std::sort(key.begin(), key.end());
  u.order.resize(nr);
  for (int i = 0; i < nr; i++) u.order[i] = key[i].second;
  return u;
}

// Converts an R object to a numeric matrix with at least 2 columns, or stops.
// Rcpp::NumericMatrix stops if `x` is not a matrix and converts integer
// storage to double.
Rcpp::NumericMatrix as_unit_matrix(SEXP x) {
  Rcpp::NumericMatrix m(x);
  if (m.nrow() > 0 && m.ncol() < 2) {
    Rcpp::stop("Each unit must be a matrix with at least 2 columns.");
  }
  return m;
}

// A view of `m` as an Eigen matrix, with no copy.
inline Eigen::Map<const Eigen::MatrixXd> as_eigen(const Rcpp::NumericMatrix& m) {
  return Eigen::Map<const Eigen::MatrixXd>(m.begin(), m.nrow(), m.ncol());
}

UnitR2 unit_from_r(SEXP x) {
  const Rcpp::NumericMatrix m = as_unit_matrix(x);
  return make_unit_r2(as_eigen(m));
}

std::vector<UnitR2> make_units_r2(const Rcpp::List& polys) {
  std::vector<UnitR2> out;
  out.reserve(polys.size());
  for (R_xlen_t i = 0; i < polys.size(); i++) {
    out.push_back(unit_from_r(polys[i]));
  }
  return out;
}

// Squared Euclidean distance from point (px, py) to edge k of unit y.
inline double seg_dist2(double px, double py, const UnitR2& y, int k) {
  const int ia = y.a[k];
  const double ax = y.v[2 * ia];
  const double ay = y.v[2 * ia + 1];
  if (y.b[k] < 0) {
    return dist2_2d(px, py, ax, ay);
  }
  const double t = ((px - ax) * y.dx[k] + (py - ay) * y.dy[k]) * y.inv_len2[k];
  const double tc = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
  const double qx = ax + tc * y.dx[k];
  const double qy = ay + tc * y.dy[k];
  return dist2_2d(px, py, qx, qy);
}

// Squared distance from (px, py) to the bounding box of edge k of y.
inline double edge_box_dist2(double px, double py, const UnitR2& y, int k) {
  double dx = 0.0;
  if (px < y.ex0[k]) dx = y.ex0[k] - px;
  else if (px > y.ex1[k]) dx = px - y.ex1[k];
  double dy = 0.0;
  if (py < y.ey0[k]) dy = y.ey0[k] - py;
  else if (py > y.ey1[k]) dy = py - y.ey1[k];
  return dx * dx + dy * dy;
}

// Exact minimum squared distance from (px, py) to edge set E(y).
double point_to_edges_dist2(double px, double py, const UnitR2& y,
                            const std::vector<int>* active_edges = NULL) {
  double best = INF;
  if (active_edges == NULL) {
    const int ne = y.ne();
    for (int k = 0; k < ne; k++) {
      if (edge_box_dist2(px, py, y, k) >= best) continue;
      const double d2 = seg_dist2(px, py, y, k);
      if (d2 < best) best = d2;
    }
  } else {
    const int m = static_cast<int>(active_edges->size());
    for (int idx = 0; idx < m; idx++) {
      const int k = (*active_edges)[idx];
      if (edge_box_dist2(px, py, y, k) >= best) continue;
      const double d2 = seg_dist2(px, py, y, k);
      if (d2 < best) best = d2;
    }
  }
  return best;
}

inline bool units_aabb_overlap(const UnitR2& x, const UnitR2& y) {
  return x.xmax >= y.xmin && x.xmin <= y.xmax &&
         x.ymax >= y.ymin && x.ymin <= y.ymax;
}

// Even-odd ray casting over all polygon rings of unit `u`. Coordinates are
// shifted by `(px, py)` first to avoid cancellation on large projected values.
bool point_in_poly_r2(double px, double py, const UnitR2& u) {
  if (!u.has_poly) return false;
  if (px < u.xmin || px > u.xmax || py < u.ymin || py > u.ymax) return false;
  bool inside = false;
  const int np = static_cast<int>(u.parts.size());
  for (int p = 0; p < np; p++) {
    const PartR2& part = u.parts[p];
    if (!part.is_poly) continue;
    if (px < part.xmin || px > part.xmax || py < part.ymin || py > part.ymax) {
      continue;
    }
    for (int k = part.e_begin; k < part.e_end; k++) {
      if (u.b[k] < 0) continue;
      const int ia = u.a[k];
      const int ib = u.b[k];
      const double y1 = u.v[2 * ia + 1] - py;
      const double y2 = u.v[2 * ib + 1] - py;
      if ((y1 > 0.0) != (y2 > 0.0)) {
        const double x1 = u.v[2 * ia] - px;
        const double x2 = u.v[2 * ib] - px;
        const double x_int = x1 - y1 * (x2 - x1) / (y2 - y1);
        if (x_int > 0.0) inside = !inside;
      }
    }
  }
  return inside;
}

inline double point_to_set_dist2(double px, double py, const UnitR2& y,
                                 bool is_polygon, bool check_pip,
                                 const std::vector<int>* active_edges = NULL) {
  if (is_polygon && check_pip && point_in_poly_r2(px, py, y)) {
    return 0.0;
  }
  return point_to_edges_dist2(px, py, y, active_edges);
}

// Directed discrete Hausdorff distance (vertices of x to edges of y), with the
// early break of Taha and Hanbury (2015, IEEE TPAMI 37(11)). It matches GEOS
// and sf::st_distance(which = "Hausdorff").
void directed_haus_r2(const UnitR2& x, const UnitR2& y, double& max_d2) {
  const int nv = x.nv();
  const int ne = y.ne();
  if (nv == 0 || ne == 0) {
    max_d2 = NA_REAL;
    return;
  }
  int start = 0;
  for (int t = 0; t < nv; t++) {
    const int idx = x.order[t];
    const double px = x.v[2 * idx];
    const double py = x.v[2 * idx + 1];
    double best = INF;
    int best_k = start;
    bool closer = false;
    for (int s = 0; s < ne; s++) {
      int k = start + s;
      if (k >= ne) k -= ne;
      const double d2 = seg_dist2(px, py, y, k);
      if (d2 <= max_d2) {
        start = k;
        closer = true;
        break;
      }
      if (d2 < best) {
        best = d2;
        best_k = k;
      }
    }
    if (closer) continue;
    start = best_k;
    max_d2 = best;
  }
}

double haus_r2_units(const UnitR2& x, const UnitR2& y) {
  if (x.nv() == 0 || y.nv() == 0) return NA_REAL;
  double max_d2 = -1.0;
  directed_haus_r2(x, y, max_d2);
  if (R_IsNA(max_d2) || std::isnan(max_d2)) return NA_REAL;
  directed_haus_r2(y, x, max_d2);
  if (R_IsNA(max_d2) || std::isnan(max_d2)) return NA_REAL;
  return max_d2 <= 0.0 ? 0.0 : std::sqrt(max_d2);
}

// ---- Continuous 1D and 2D Hausdorff helpers --------------------------------

// Liang-Barsky test whether segment AB intersects box [x0, x1] x [y0, y1].
// Shifts endpoints relative to box center before multiplying.
bool seg_intersects_box(double ax, double ay, double bx, double by,
                        double x0, double y0, double x1, double y1) {
  const double cx = 0.5 * (x0 + x1);
  const double cy = 0.5 * (y0 + y1);
  const double hx = 0.5 * (x1 - x0);
  const double hy = 0.5 * (y1 - y0);
  const double px = ax - cx;
  const double py = ay - cy;
  const double dx = bx - ax;
  const double dy = by - ay;
  double t0 = 0.0;
  double t1 = 1.0;

  auto clip = [&](double p, double q) -> bool {
    if (p == 0.0) return q >= 0.0;
    const double r = q / p;
    if (p < 0.0) {
      if (r > t1) return false;
      if (r > t0) t0 = r;
    } else {
      if (r < t0) return false;
      if (r < t1) t1 = r;
    }
    return true;
  };

  if (!clip(-dx, px + hx)) return false;
  if (!clip( dx, hx - px)) return false;
  if (!clip(-dy, py + hy)) return false;
  if (!clip( dy, hy - py)) return false;
  return t0 <= t1;
}

struct Pt2 {
  double x, y;
};

// Clips convex polygon `poly` to nx*(p.x - ox) + ny*(p.y - oy) + c0 >= 0.
std::vector<Pt2> clip_linear_halfplane(const std::vector<Pt2>& poly,
                                       double ox, double oy,
                                       double nx, double ny, double c0) {
  std::vector<Pt2> out;
  const int n = static_cast<int>(poly.size());
  if (n == 0) return out;

  auto signed_val = [&](const Pt2& p) -> double {
    return (p.x - ox) * nx + (p.y - oy) * ny + c0;
  };

  Pt2 prev = poly[n - 1];
  double d_prev = signed_val(prev);
  for (int i = 0; i < n; i++) {
    const Pt2 cur = poly[i];
    const double d_cur = signed_val(cur);
    if (d_cur >= 0.0) {
      if (d_prev < 0.0) {
        const double t = d_prev / (d_prev - d_cur);
        out.push_back({prev.x + t * (cur.x - prev.x),
                       prev.y + t * (cur.y - prev.y)});
      }
      out.push_back(cur);
    } else if (d_prev >= 0.0) {
      const double t = d_prev / (d_prev - d_cur);
      out.push_back({prev.x + t * (cur.x - prev.x),
                     prev.y + t * (cur.y - prev.y)});
    }
    prev = cur;
    d_prev = d_cur;
  }
  return out;
}

// Clips convex polygon `poly` against half-plane nx*(x - ax) + ny*(y - ay) >= 0.
std::vector<Pt2> clip_poly_halfplane(const std::vector<Pt2>& poly,
                                     double ax, double ay,
                                     double nx, double ny, double sign) {
  return clip_linear_halfplane(poly, ax, ay, sign * nx, sign * ny, 0.0);
}

// Returns true if cell [x0, x1] x [y0, y1] cannot reach `target` anywhere,
// proved by intersecting the cell with linear upper-bound half-planes for
// active edges of `y`.
bool ridge_prune_r2(double x0, double y0, double x1, double y1,
                    double cx, double cy, double r,
                    const UnitR2& y, const std::vector<int>& y_next,
                    double target) {
  if (target <= 0.0 || y_next.size() < 2) return false;
  std::vector<Pt2> poly = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  const double r_lim = r * (1.0 + 1e-12);

  for (int ky : y_next) {
    const int ia = y.a[ky];
    const double ax = y.v[2 * ia];
    const double ay = y.v[2 * ia + 1];
    double qx = ax, qy = ay;
    double d2;
    bool on_interior = false;
    if (y.b[ky] < 0) {
      d2 = dist2_2d(cx, cy, ax, ay);
    } else {
      const double t = ((cx - ax) * y.dx[ky] + (cy - ay) * y.dy[ky]) * y.inv_len2[ky];
      const double tc = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
      qx = ax + tc * y.dx[ky];
      qy = ay + tc * y.dy[ky];
      d2 = dist2_2d(cx, cy, qx, qy);
      on_interior = (t > 0.0 && t < 1.0);
    }
    const double dk = std::sqrt(d2);
    if (dk <= r_lim) continue;

    bool used_line = false;
    if (on_interior) {
      const int ib = y.b[ky];
      const double lim2 = (dk + 2.0 * r + 1e-12 * std::max(1.0, r));
      const double lim_sq = lim2 * lim2;
      if (dist2_2d(cx, cy, ax, ay) > lim_sq &&
          dist2_2d(cx, cy, y.v[2 * ib], y.v[2 * ib + 1]) > lim_sq) {
        const double ux = (cx - qx) / dk;
        const double uy = (cy - qy) / dk;
        poly = clip_linear_halfplane(poly, cx, cy, ux, uy, dk - target);
        used_line = true;
      }
    }
    if (!used_line) {
      const double ux = (cx - qx) / dk;
      const double uy = (cy - qy) / dk;
      const double err = 0.5 * r * r / (dk - r);
      poly = clip_linear_halfplane(poly, cx, cy, ux, uy, dk + err - target);
    }
    if (poly.empty()) return true;
  }
  return false;
}

double poly_area2(const std::vector<Pt2>& poly) {
  const int n = static_cast<int>(poly.size());
  if (n < 3) return 0.0;
  const double ox = poly[0].x;
  const double oy = poly[0].y;
  double acc = 0.0;
  for (int i = 1; i + 1 < n; i++) {
    const double x1 = poly[i].x - ox;
    const double y1 = poly[i].y - oy;
    const double x2 = poly[i + 1].x - ox;
    const double y2 = poly[i + 1].y - oy;
    acc += x1 * y2 - x2 * y1;
  }
  return std::fabs(acc);
}

Pt2 poly_centroid(const std::vector<Pt2>& poly) {
  const int n = static_cast<int>(poly.size());
  double sx = 0.0, sy = 0.0;
  for (int i = 0; i < n; i++) {
    sx += poly[i].x;
    sy += poly[i].y;
  }
  return {sx / n, sy / n};
}

// Initial vertex-to-set pass giving a strong lower bound `best_d` and upper
// bounds `d_v` for every vertex of `x` (used to prune 1D segments).
void vertex_pass_set_r2(const UnitR2& x, const UnitR2& y, bool is_polygon,
                        bool check_pip, double& best_d,
                        std::vector<double>& d_v) {
  const int nv = x.nv();
  const int ne = y.ne();
  d_v.assign(nv, 0.0);
  double max_d2 = best_d * best_d;
  if (best_d <= 0.0) max_d2 = -1.0;

  int start = 0;
  for (int t = 0; t < nv; t++) {
    const int idx = x.order[t];
    const double px = x.v[2 * idx];
    const double py = x.v[2 * idx + 1];
    if (is_polygon && check_pip && point_in_poly_r2(px, py, y)) {
      d_v[idx] = 0.0;
      continue;
    }
    double best2 = INF;
    int best_k = start;
    bool closer = false;
    for (int s = 0; s < ne; s++) {
      int k = start + s;
      if (k >= ne) k -= ne;
      if (edge_box_dist2(px, py, y, k) >= best2) continue;
      const double d2 = seg_dist2(px, py, y, k);
      if (d2 <= max_d2) {
        start = k;
        d_v[idx] = std::sqrt(d2);
        closer = true;
        break;
      }
      if (d2 < best2) {
        best2 = d2;
        best_k = k;
      }
    }
    if (closer) continue;
    start = best_k;
    max_d2 = best2;
    const double d = std::sqrt(best2);
    d_v[idx] = d;
    if (d > best_d) best_d = d;
  }
}

struct SegTask {
  double ax, ay, bx, by;
  double da, db;
  double half_len;
  double upper;
  std::vector<int> active;
  bool operator<(const SegTask& o) const { return upper < o.upper; }
};

// Directed continuous 1D Hausdorff from edges/points of `x` to `y`
// (either E(y) if !is_polygon, or filled y \cup E(y) if is_polygon).
void directed_haus_1d_r2(const UnitR2& x, const UnitR2& y, bool is_polygon,
                         bool check_pip, const std::vector<double>& d_v,
                         double tol_eff, double min_h, double& best_d) {
  const int ne_x = x.ne();
  const int ne_y = y.ne();
  if (ne_x == 0 || ne_y == 0) return;

  std::vector<int> all_y(ne_y);
  for (int k = 0; k < ne_y; k++) all_y[k] = k;

  std::priority_queue<SegTask> pq;
  for (int k = 0; k < ne_x; k++) {
    if (x.b[k] < 0) continue;  // Point parts already finished in vertex_pass
    const int ia = x.a[k];
    const int ib = x.b[k];
    const double ax = x.v[2 * ia];
    const double ay = x.v[2 * ia + 1];
    const double bx = x.v[2 * ib];
    const double by = x.v[2 * ib + 1];
    const double da = d_v[ia];
    const double db = d_v[ib];
    const double h = 0.5 * std::sqrt(x.dx[k] * x.dx[k] + x.dy[k] * x.dy[k]);
    const double upper = 0.5 * (da + db) + h;
    if (upper > best_d + tol_eff && h > min_h) {
      SegTask task;
      task.ax = ax; task.ay = ay; task.bx = bx; task.by = by;
      task.da = da; task.db = db; task.half_len = h; task.upper = upper;
      task.active = all_y;
      pq.push(std::move(task));
    }
  }

  int iter = 0;
  while (!pq.empty()) {
    if ((++iter & 1023) == 0) Rcpp::checkUserInterrupt();
    SegTask cur = pq.top();
    pq.pop();
    if (cur.upper <= best_d + tol_eff || cur.half_len <= min_h) {
      if (cur.upper <= best_d + tol_eff) break;
      continue;
    }

    const double mx = 0.5 * (cur.ax + cur.bx);
    const double my = 0.5 * (cur.ay + cur.by);
    const double h = cur.half_len;
    const double min_r2 = std::max(0.0, (best_d - h) * (best_d - h));
    const bool check_early = (best_d > h);

    // Check whether midpoint is inside filled polygon y.
    const bool in_y = is_polygon && check_pip && point_in_poly_r2(mx, my, y);

    double best_edge2 = INF;
    bool pruned = false;
    const int na = static_cast<int>(cur.active.size());
    std::vector<double> d_edge(na, INF);
    for (int i = 0; i < na; i++) {
      const int ky = cur.active[i];
      if (edge_box_dist2(mx, my, y, ky) >= best_edge2 && !check_early) continue;
      const double d2 = seg_dist2(mx, my, y, ky);
      d_edge[i] = d2;
      if (d2 < best_edge2) best_edge2 = d2;
      if (check_early && d2 <= min_r2) {
        pruned = true;
        break;
      }
      // Convex edge-cover prune: if a single edge covers both endpoints within
      // best_d, the entire segment is <= best_d everywhere. By convexity, the
      // midpoint distance must also be <= best_d before checking endpoints.
      if (d2 <= best_d * best_d &&
          seg_dist2(cur.ax, cur.ay, y, ky) <= best_d * best_d &&
          seg_dist2(cur.bx, cur.by, y, ky) <= best_d * best_d) {
        pruned = true;
        break;
      }
    }
    if (pruned) continue;

    // If midpoint lies strictly inside y at distance >= h from all boundary edges,
    // the entire segment lies inside y (d = 0).
    if (in_y && best_edge2 >= h * h * (1.0 + 1e-12)) continue;

    const double d_edge_min = std::sqrt(best_edge2);
    const double dm = in_y ? 0.0 : d_edge_min;
    if (dm > best_d) best_d = dm;

    // Filter active edges for children via triangle inequality: any edge
    // farther than d(m, E(y)) + 2*h cannot be closest to any point on S.
    const double keep_lim = d_edge_min + 2.0 * h + 1e-12 * std::max(1.0, h);
    const double keep_lim2 = keep_lim * keep_lim;
    std::vector<int> child_active;
    child_active.reserve(na);
    for (int i = 0; i < na; i++) {
      const int ky = cur.active[i];
      double d2 = d_edge[i];
      if (d2 == INF) d2 = seg_dist2(mx, my, y, ky);
      if (d2 <= keep_lim2) child_active.push_back(ky);
    }

    // 1-edge convexity short-circuit.
    if (child_active.size() == 1) {
      const int ky = child_active[0];
      bool safe_1edge = !is_polygon || !y.has_poly;
      if (!safe_1edge) {
        if (y.b[ky] < 0) {
          safe_1edge = true;
        } else {
          const int va = y.a[ky], vb = y.b[ky];
          safe_1edge = (dist2_2d(mx, my, y.v[2 * va], y.v[2 * va + 1]) > keep_lim2) &&
                       (dist2_2d(mx, my, y.v[2 * vb], y.v[2 * vb + 1]) > keep_lim2);
        }
      }
      if (safe_1edge) {
        const double da_exact = std::sqrt(point_to_set_dist2(cur.ax, cur.ay, y, is_polygon, check_pip, &child_active));
        const double db_exact = std::sqrt(point_to_set_dist2(cur.bx, cur.by, y, is_polygon, check_pip, &child_active));
        best_d = std::max(best_d, std::max(da_exact, db_exact));
        continue;
      }
    }

    const double child_h = 0.5 * h;
    if (child_h <= min_h) continue;
    const double u_left = 0.5 * (cur.da + dm + h);
    if (u_left > best_d + tol_eff) {
      SegTask left;
      left.ax = cur.ax; left.ay = cur.ay; left.bx = mx; left.by = my;
      left.da = cur.da; left.db = dm; left.half_len = child_h;
      left.upper = u_left; left.active = child_active;
      pq.push(std::move(left));
    }
    const double u_right = 0.5 * (dm + cur.db + h);
    if (u_right > best_d + tol_eff) {
      SegTask right;
      right.ax = mx; right.ay = my; right.bx = cur.bx; right.by = cur.by;
      right.da = dm; right.db = cur.db; right.half_len = child_h;
      right.upper = u_right; right.active = std::move(child_active);
      pq.push(std::move(right));
    }
  }
}

struct CellTask {
  double x0, y0, x1, y1;
  double radius;
  double upper;
  std::vector<int> x_poly_edges;
  std::vector<int> y_edges;
  bool operator<(const CellTask& o) const { return upper < o.upper; }
};

// Directed 2D interior Hausdorff over polygon parts of `x` against filled `y`.
void directed_haus_poly_2d_r2(const UnitR2& x, const UnitR2& y,
                              bool check_pip_y, double tol_eff,
                              double min_r, double& best_d) {
  if (!x.has_poly || y.ne() == 0) return;

  std::vector<int> x_poly_all;
  x_poly_all.reserve(x.ne());
  for (int k = 0; k < x.ne(); k++) {
    if (x.e_is_poly[k] && x.b[k] >= 0) x_poly_all.push_back(k);
  }
  if (x_poly_all.empty()) return;

  std::vector<int> y_all(y.ne());
  for (int k = 0; k < y.ne(); k++) y_all[k] = k;

  std::vector<int> poly_parts;
  const int np = static_cast<int>(x.parts.size());
  for (int p = 0; p < np; p++) {
    if (x.parts[p].is_poly) poly_parts.push_back(p);
  }
  std::sort(poly_parts.begin(), poly_parts.end(), [&](int p1, int p2) {
    const PartR2& a = x.parts[p1];
    const PartR2& b = x.parts[p2];
    const double area1 = (a.xmax - a.xmin) * (a.ymax - a.ymin);
    const double area2 = (b.xmax - b.xmin) * (b.ymax - b.ymin);
    return area1 > area2;
  });

  std::vector<PartR2> roots;
  std::priority_queue<CellTask> pq;
  for (int p : poly_parts) {
    const PartR2& part = x.parts[p];
    bool covered = false;
    for (const PartR2& rt : roots) {
      if (part.xmin >= rt.xmin && part.xmax <= rt.xmax &&
          part.ymin >= rt.ymin && part.ymax <= rt.ymax) {
        covered = true;
        break;
      }
    }
    if (covered) continue;  // Holes and enclosed islands are covered by outer root
    roots.push_back(part);

    CellTask root;
    root.x0 = part.xmin; root.y0 = part.ymin;
    root.x1 = part.xmax; root.y1 = part.ymax;
    const double hx = 0.5 * (root.x1 - root.x0);
    const double hy = 0.5 * (root.y1 - root.y0);
    root.radius = std::sqrt(hx * hx + hy * hy);
    if (root.radius <= min_r) continue;
    root.upper = INF;
    root.x_poly_edges = x_poly_all;
    root.y_edges = y_all;
    pq.push(std::move(root));
  }

  int iter = 0;
  while (!pq.empty()) {
    if ((++iter & 1023) == 0) Rcpp::checkUserInterrupt();
    CellTask cur = pq.top();
    pq.pop();
    if (cur.upper <= best_d + tol_eff) break;
    if (cur.radius <= min_r) continue;

    const double cx = 0.5 * (cur.x0 + cur.x1);
    const double cy = 0.5 * (cur.y0 + cur.y1);
    const double r = cur.radius;

    // 1. Filter x polygon edges intersecting this cell.
    std::vector<int> x_cur;
    x_cur.reserve(cur.x_poly_edges.size());
    for (int k : cur.x_poly_edges) {
      const int ia = x.a[k], ib = x.b[k];
      if (seg_intersects_box(x.v[2 * ia], x.v[2 * ia + 1],
                             x.v[2 * ib], x.v[2 * ib + 1],
                             cur.x0, cur.y0, cur.x1, cur.y1)) {
        x_cur.push_back(k);
      }
    }
    const bool in_x = point_in_poly_r2(cx, cy, x);
    if (x_cur.empty() && !in_x) continue;  // Cell strictly outside x

    // 2. Scan active y edges from cell center with early break.
    const bool in_y = check_pip_y && point_in_poly_r2(cx, cy, y);
    const double min_r2 = std::max(0.0, (best_d - r) * (best_d - r));
    const bool check_early = (best_d > r);

    const int ny = static_cast<int>(cur.y_edges.size());
    std::vector<double> dy_edge(ny, INF);
    double best_edge2 = INF;
    double best_poly_edge2 = INF;
    bool pruned = false;

    for (int i = 0; i < ny; i++) {
      const int ky = cur.y_edges[i];
      if (edge_box_dist2(cx, cy, y, ky) >= best_edge2 && !check_early && !in_y) {
        continue;
      }
      const double d2 = seg_dist2(cx, cy, y, ky);
      dy_edge[i] = d2;
      if (d2 < best_edge2) best_edge2 = d2;
      if (y.e_is_poly[ky] && d2 < best_poly_edge2) best_poly_edge2 = d2;
      if (check_early && d2 <= min_r2) {
        pruned = true;
        break;
      }
    }
    if (pruned) continue;

    // Inside-y prune: entire cell is inside filled polygon y.
    const double r_safe = r * (1.0 + 1e-12);
    if (in_y && best_poly_edge2 >= r_safe * r_safe) continue;

    const double d_edge_min = std::sqrt(best_edge2);
    const double dc = in_y ? 0.0 : d_edge_min;
    if (in_x && dc > best_d) best_d = dc;

    const double upper = dc + r;
    if (upper <= best_d + tol_eff) continue;

    // Filter active y edges for sub-cells and 1-site check.
    const double keep_lim = d_edge_min + 2.0 * r + 1e-12 * std::max(1.0, r);
    const double keep_lim2 = keep_lim * keep_lim;
    std::vector<int> y_next;
    y_next.reserve(ny);
    for (int i = 0; i < ny; i++) {
      const int ky = cur.y_edges[i];
      double d2 = dy_edge[i];
      if (d2 == INF) {
        if (edge_box_dist2(cx, cy, y, ky) > keep_lim2) continue;
        d2 = seg_dist2(cx, cy, y, ky);
      }
      if (d2 <= keep_lim2) y_next.push_back(ky);
    }

    // 3. 1-site / convex terminal cell prune (prevents boundary subdivision blowup).
    if (y_next.size() == 1 && x_cur.size() <= 1) {
      const int ky = y_next[0];
      bool safe_y = false;
      if (!y.has_poly || y.e_is_poly[ky] || !in_y) {
        if (y.b[ky] < 0) {
          safe_y = true;
        } else {
          const int va = y.a[ky], vb = y.b[ky];
          safe_y = (dist2_2d(cx, cy, y.v[2 * va], y.v[2 * va + 1]) > keep_lim2) &&
                   (dist2_2d(cx, cy, y.v[2 * vb], y.v[2 * vb + 1]) > keep_lim2);
        }
      }

      if (safe_y) {
        std::vector<Pt2> cand_pts;
        bool safe_x = false;
        const std::vector<Pt2> box = {
          {cur.x0, cur.y0}, {cur.x1, cur.y0}, {cur.x1, cur.y1}, {cur.x0, cur.y1}
        };
        if (x_cur.empty()) {
          if (in_x) {
            cand_pts = box;
            safe_x = true;
          }
        } else {
          const int kx = x_cur[0];
          const int va = x.a[kx], vb = x.b[kx];
          const double ax = x.v[2 * va], ay = x.v[2 * va + 1];
          const double bx = x.v[2 * vb], by = x.v[2 * vb + 1];
          const double pad = 1e-12 * std::max(1.0, r);
          const bool a_in = ax >= cur.x0 - pad && ax <= cur.x1 + pad &&
                            ay >= cur.y0 - pad && ay <= cur.y1 + pad;
          const bool b_in = bx >= cur.x0 - pad && bx <= cur.x1 + pad &&
                            by >= cur.y0 - pad && by <= cur.y1 + pad;
          if (!a_in && !b_in) {
            const double len = std::sqrt(x.dx[kx] * x.dx[kx] + x.dy[kx] * x.dy[kx]);
            if (len > 0.0) {
              const double nx = -x.dy[kx] / len;
              const double ny = x.dx[kx] / len;
              const std::vector<Pt2> pos = clip_poly_halfplane(box, ax, ay, nx, ny, 1.0);
              const std::vector<Pt2> neg = clip_poly_halfplane(box, ax, ay, nx, ny, -1.0);
              const double area_eps = 1e-24 * r * r;
              const bool pos_valid = poly_area2(pos) > area_eps;
              const bool neg_valid = poly_area2(neg) > area_eps;
              const bool pos_in = pos_valid && [&]() {
                const Pt2 c = poly_centroid(pos);
                return point_in_poly_r2(c.x, c.y, x);
              }();
              const bool neg_in = neg_valid && [&]() {
                const Pt2 c = poly_centroid(neg);
                return point_in_poly_r2(c.x, c.y, x);
              }();
              if (pos_in && neg_in) {
                cand_pts = box;
                safe_x = true;
              } else if (pos_in) {
                cand_pts = pos;
                safe_x = true;
              } else if (neg_in) {
                cand_pts = neg;
                safe_x = true;
              } else if (pos_valid && neg_valid) {
                // Degenerate slit outside x on both sides (handled in Stage 1).
                safe_x = true;
              }
            }
          }
        }
        if (safe_x) {
          for (const Pt2& pt : cand_pts) {
            const double d2 = point_to_set_dist2(pt.x, pt.y, y, true, check_pip_y, &y_next);
            const double d = std::sqrt(d2);
            if (d > best_d) best_d = d;
          }
          continue;
        }
      }
    }

    // 4. Ridge-clipping prune: check if linearized edge bounds make
    // best_d + tol_eff unreachable anywhere inside this cell.
    if (best_d > 0.0 && y_next.size() >= 2 &&
        ridge_prune_r2(cur.x0, cur.y0, cur.x1, cur.y1, cx, cy, r,
                       y, y_next, best_d + tol_eff)) {
      continue;
    }

    // 5. Split cell into 4 quadrants.
    const double qx[2][2] = {{cur.x0, cx}, {cx, cur.x1}};
    const double qy[2][2] = {{cur.y0, cy}, {cy, cur.y1}};
    const double child_r = 0.5 * r;
    if (child_r <= min_r) continue;
    for (int ix = 0; ix < 2; ix++) {
      for (int iy = 0; iy < 2; iy++) {
        CellTask child;
        child.x0 = qx[ix][0]; child.x1 = qx[ix][1];
        child.y0 = qy[iy][0]; child.y1 = qy[iy][1];
        child.radius = child_r;
        child.upper = upper;
        child.x_poly_edges = x_cur;
        child.y_edges = y_next;
        pq.push(std::move(child));
      }
    }
  }
}

double haus_r2_method_units(const UnitR2& x_in, const UnitR2& y_in,
                            int method_code, double tol) {
  if (x_in.nv() == 0 || y_in.nv() == 0) return NA_REAL;
  if (method_code == 0) {
    return haus_r2_units(x_in, y_in);
  }
  if (x_in.v == y_in.v && x_in.a == y_in.a && x_in.b == y_in.b &&
      (method_code == 1 || x_in.e_is_poly == y_in.e_is_poly)) {
    return 0.0;
  }

  // Canonicalize unit order for bit-for-bit symmetry between (x, y) and (y, x).
  const bool swap_order =
    (y_in.nv() < x_in.nv()) ||
    (y_in.nv() == x_in.nv() && (y_in.ne() < x_in.ne() ||
      (y_in.ne() == x_in.ne() && std::make_pair(y_in.xmin, y_in.ymin) <
                                 std::make_pair(x_in.xmin, x_in.ymin))));
  const UnitR2& a = swap_order ? y_in : x_in;
  const UnitR2& b = swap_order ? x_in : y_in;

  const bool is_polygon = (method_code == 2);
  const bool overlap = is_polygon && units_aabb_overlap(a, b);
  const bool check_pip_b = overlap && b.has_poly;
  const bool check_pip_a = overlap && a.has_poly;

  const double scale = std::max(1.0, std::max(a.diam, b.diam));
  const double max_coord = std::max(a.max_abs, b.max_abs);
  const double floor_tol = 64.0 * EPS * std::max(1.0, max_coord);
  const double tol_eff = std::max(tol * scale, floor_tol);
  const double min_step = 8.0 * EPS * std::max(1.0, max_coord);

  double best_d = 0.0;
  std::vector<double> dv_a, dv_b;
  vertex_pass_set_r2(a, b, is_polygon, check_pip_b, best_d, dv_a);
  vertex_pass_set_r2(b, a, is_polygon, check_pip_a, best_d, dv_b);

  directed_haus_1d_r2(a, b, is_polygon, check_pip_b, dv_a, tol_eff, min_step, best_d);
  directed_haus_1d_r2(b, a, is_polygon, check_pip_a, dv_b, tol_eff, min_step, best_d);

  if (is_polygon) {
    directed_haus_poly_2d_r2(a, b, check_pip_b, tol_eff, min_step, best_d);
    directed_haus_poly_2d_r2(b, a, check_pip_a, tol_eff, min_step, best_d);
  }
  return best_d;
}

int parse_method_code(const std::string& method) {
  if (method == "discrete") return 0;
  if (method == "boundary") return 1;
  if (method == "polygon") return 2;
  Rcpp::stop("Unknown method: must be \"discrete\", \"boundary\", or \"polygon\".");
}

}  // namespace

//' @title Planar Hausdorff distance between two units (internal)
//' @param x,y matrices with the columns x, y, and optionally a part id.
//' @noRd
// [[Rcpp::export]]
double haus_r2_aux_continuous(SEXP x, SEXP y) {
  return haus_r2_units(unit_from_r(x), unit_from_r(y));
}

//' @title Planar Hausdorff distance matrix (internal)
//' @param poly_list list of matrices from .to_list_mat().
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd dist_haus_r2(const Rcpp::List& poly_list) {
  const int n = poly_list.size();
  const std::vector<UnitR2> u = make_units_r2(poly_list);
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n);

  for (int i = 0; i < n; i++) {
    Rcpp::checkUserInterrupt();
    if (u[i].nv() == 0) out(i, i) = NA_REAL;
    for (int j = i + 1; j < n; j++) {
      const double d = haus_r2_units(u[i], u[j]);
      out(i, j) = d;
      out(j, i) = d;  // The matrix is symmetric
    }
  }
  return out;
}

//' @title Planar Hausdorff distances from units to points (internal)
//' @param poly_list list of matrices from .to_list_mat().
//' @param pred_mat matrix with one point (x, y) for each row.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cdist_haus_r2_mp(const Rcpp::List& poly_list, SEXP pred_mat) {
  const Rcpp::NumericMatrix pm = as_unit_matrix(pred_mat);
  const Eigen::Map<const Eigen::MatrixXd> pred = as_eigen(pm);
  const int n1 = poly_list.size();
  const int n2 = static_cast<int>(pred.rows());
  const std::vector<UnitR2> u = make_units_r2(poly_list);
  // Each row of pred_mat is a unit with one point
  std::vector<UnitR2> p;
  p.reserve(n2);
  for (int j = 0; j < n2; j++) {
    p.push_back(make_unit_r2(Eigen::MatrixXd(pred.row(j))));
  }
  Eigen::MatrixXd out(n1, n2);

  for (int i = 0; i < n1; i++) {
    Rcpp::checkUserInterrupt();
    for (int j = 0; j < n2; j++) {
      out(i, j) = haus_r2_units(u[i], p[j]);
    }
  }
  return out;
}

//' @title Planar Hausdorff cross-distance matrix (internal)
//' @param poly_list,poly_pred lists of matrices from .to_list_mat().
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cdist_haus_r2_lists(const Rcpp::List& poly_list,
                                    const Rcpp::List& poly_pred) {
  const int n1 = poly_list.size();
  const int n2 = poly_pred.size();
  const std::vector<UnitR2> u1 = make_units_r2(poly_list);
  const std::vector<UnitR2> u2 = make_units_r2(poly_pred);
  Eigen::MatrixXd out(n1, n2);

  for (int i = 0; i < n1; i++) {
    Rcpp::checkUserInterrupt();
    for (int j = 0; j < n2; j++) {
      out(i, j) = haus_r2_units(u1[i], u2[j]);
    }
  }
  return out;
}

//' @title Planar Hausdorff distance matrix with method selection (internal)
//' @param poly_list list of matrices from .to_list_mat() or .to_polygon_list().
//' @param method one of "discrete", "boundary", or "polygon".
//' @param tol convergence tolerance for continuous methods.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd dist_haus_r2_method(const Rcpp::List& poly_list,
                                    const std::string& method = "discrete",
                                    double tol = 1e-8) {
  const int method_code = parse_method_code(method);
  if (method_code != 0 && (!std::isfinite(tol) || tol <= 0.0)) {
    Rcpp::stop("tol must be a single finite positive number.");
  }
  if (method_code == 0) return dist_haus_r2(poly_list);

  const int n = poly_list.size();
  const std::vector<UnitR2> u = make_units_r2(poly_list);
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n);

  for (int i = 0; i < n; i++) {
    Rcpp::checkUserInterrupt();
    if (u[i].nv() == 0) out(i, i) = NA_REAL;
    for (int j = i + 1; j < n; j++) {
      const double d = haus_r2_method_units(u[i], u[j], method_code, tol);
      out(i, j) = d;
      out(j, i) = d;
    }
  }
  return out;
}

//' @title Planar Hausdorff cross-distance matrix with method selection (internal)
//' @param poly_list,poly_pred lists of matrices from .to_list_mat() or
//'   .to_polygon_list().
//' @param method one of "discrete", "boundary", or "polygon".
//' @param tol convergence tolerance for continuous methods.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cdist_haus_r2_lists_method(const Rcpp::List& poly_list,
                                           const Rcpp::List& poly_pred,
                                           const std::string& method = "discrete",
                                           double tol = 1e-8) {
  const int method_code = parse_method_code(method);
  if (method_code != 0 && (!std::isfinite(tol) || tol <= 0.0)) {
    Rcpp::stop("tol must be a single finite positive number.");
  }
  if (method_code == 0) return cdist_haus_r2_lists(poly_list, poly_pred);

  const int n1 = poly_list.size();
  const int n2 = poly_pred.size();
  const std::vector<UnitR2> u1 = make_units_r2(poly_list);
  const std::vector<UnitR2> u2 = make_units_r2(poly_pred);
  Eigen::MatrixXd out(n1, n2);

  for (int i = 0; i < n1; i++) {
    Rcpp::checkUserInterrupt();
    for (int j = 0; j < n2; j++) {
      out(i, j) = haus_r2_method_units(u1[i], u2[j], method_code, tol);
    }
  }
  return out;
}
