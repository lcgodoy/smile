#include "haus-common.h"
#include <cstdint>
#include <queue>

// [[Rcpp::depends(RcppEigen)]]

// Planar (R^2) Hausdorff distance between spatial units.
//
// Supports three methods:
// - "discrete": discrete vertex-to-edge Hausdorff distance (matches GEOS).
// - "boundary": continuous Hausdorff distance between boundaries and linework.
// - "polygon":  continuous Hausdorff distance between closed spatial sets
//               (polygons include 2D interiors excluding holes).
//
// Every distance in this file is in the units of the input coordinates.

namespace {

using smile_geom::Pt2;
using smile_geom::cross_dist_matrix;
using smile_geom::dist2_2d;
using smile_geom::golden_order;
using smile_geom::kEps;
using smile_geom::kInf;
using smile_geom::sym_dist_matrix;

using smile_haus::HausMethod;
using smile_haus::PolyPart2D;
using smile_haus::RawPart;
using smile_haus::extract_raw_parts;
using smile_haus::has_valid_coords;
using smile_haus::kBoundary;
using smile_haus::kDiscrete;
using smile_haus::kPolygon;
using smile_haus::make_unit_vec;
using smile_haus::parse_haus_method;
using smile_haus::point_in_poly_2d;
using smile_haus::seg_intersects_box;
using smile_haus::should_swap_units;

// 1. A spatial unit in R^2, prepared once for all its pairwise distances.
struct UnitR2 {
  std::vector<Pt2> v;             // Vertex coordinates
  std::vector<int> order;         // Golden-ratio vertex order
  std::vector<int> a, b;          // Edge endpoints (b == -1 for points)
  std::vector<double> dx, dy;     // Edge vector B - A
  std::vector<double> inv_len2;   // 1 / ||B - A||^2
  std::vector<double> ex0, ey0;   // Edge bounding-box min
  std::vector<double> ex1, ey1;   // Edge bounding-box max
  std::vector<uint8_t> e_is_poly; // 1 if edge belongs to a polygon ring
  std::vector<PolyPart2D> parts;  // Per-part metadata ("polygon")
  double xmin, ymin, xmax, ymax;  // Unit bounding box
  double diam;                    // Unit bounding-box diagonal
  double max_abs;                 // Max(|x|, |y|) for ULP floor
  bool has_poly;

  int nv() const { return static_cast<int>(order.size()); }
  int ne() const { return static_cast<int>(a.size()); }
  bool valid() const { return nv() > 0; }
};

void add_edge_r2(UnitR2& u, int i, int j, bool is_poly) {
  const double ax = u.v[i].x;
  const double ay = u.v[i].y;
  double dx = 0.0;
  double dy = 0.0;
  double inv_len2 = 0.0;
  double x0 = ax, y0 = ay, x1 = ax, y1 = ay;
  if (j >= 0) {
    const double bx = u.v[j].x;
    const double by = u.v[j].y;
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

UnitR2 make_unit_r2(const Eigen::Ref<const Eigen::MatrixXd>& poly,
                    HausMethod method) {
  UnitR2 u;
  u.xmin = u.ymin = kInf;
  u.xmax = u.ymax = -kInf;
  u.diam = 0.0;
  u.max_abs = 0.0;
  u.has_poly = false;

  if (!has_valid_coords(poly)) return u;
  const int nr = static_cast<int>(poly.rows());
  u.v.resize(nr);
  for (int i = 0; i < nr; i++) {
    const double x = poly(i, 0);
    const double y = poly(i, 1);
    u.v[i] = {x, y};
    u.xmin = std::min(u.xmin, x);
    u.ymin = std::min(u.ymin, y);
    u.xmax = std::max(u.xmax, x);
    u.ymax = std::max(u.ymax, y);
    u.max_abs = std::max(u.max_abs, std::max(std::fabs(x), std::fabs(y)));
  }
  const double bw = u.xmax - u.xmin;
  const double bh = u.ymax - u.ymin;
  u.diam = std::sqrt(bw * bw + bh * bh);

  const bool need_poly = (method == kPolygon);
  const std::vector<RawPart> raw_parts = extract_raw_parts(poly);
  for (const RawPart& rp : raw_parts) {
    const bool part_poly = need_poly && rp.is_poly;
    const int e_begin = u.ne();
    if (rp.end - rp.begin == 1) {
      add_edge_r2(u, rp.begin, -1, false);
    } else {
      for (int i = rp.begin; i + 1 < rp.end; i++) {
        add_edge_r2(u, i, i + 1, part_poly);
      }
      if (rp.close_ring) {
        add_edge_r2(u, rp.end - 1, rp.begin, part_poly);
      }
    }

    if (need_poly) {
      if (part_poly) u.has_poly = true;
      PolyPart2D part = {e_begin, u.ne(), part_poly, kInf, kInf, -kInf, -kInf};
      for (int i = rp.begin; i < rp.end; i++) {
        part.xmin = std::min(part.xmin, u.v[i].x);
        part.ymin = std::min(part.ymin, u.v[i].y);
        part.xmax = std::max(part.xmax, u.v[i].x);
        part.ymax = std::max(part.ymax, u.v[i].y);
      }
      u.parts.push_back(part);
    }
  }

  u.order = golden_order(nr);
  return u;
}

// 2. Point-to-edge and point-in-polygon primitives in R^2.

inline double seg_dist2(double px, double py, const UnitR2& y, int k) {
  const Pt2& A = y.v[y.a[k]];
  if (y.b[k] < 0) {
    return dist2_2d(px, py, A.x, A.y);
  }
  const double t =
    ((px - A.x) * y.dx[k] + (py - A.y) * y.dy[k]) * y.inv_len2[k];
  const double tc = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
  return dist2_2d(px, py, A.x + tc * y.dx[k], A.y + tc * y.dy[k]);
}

inline double edge_box_dist2(double px, double py, const UnitR2& y, int k) {
  double dx = 0.0;
  if (px < y.ex0[k]) dx = y.ex0[k] - px;
  else if (px > y.ex1[k]) dx = px - y.ex1[k];
  double dy = 0.0;
  if (py < y.ey0[k]) dy = y.ey0[k] - py;
  else if (py > y.ey1[k]) dy = py - y.ey1[k];
  return dx * dx + dy * dy;
}

inline bool point_in_poly_r2(double px, double py, const UnitR2& u) {
  if (!u.has_poly) return false;
  return point_in_poly_2d(px, py, u.v, u.a, u.b, u.parts,
                          u.xmin, u.ymin, u.xmax, u.ymax);
}

double point_to_set_dist2(double px, double py, const UnitR2& y,
                          bool check_pip,
                          const std::vector<int>& active_edges) {
  if (check_pip && point_in_poly_r2(px, py, y)) return 0.0;
  double best = kInf;
  for (int k : active_edges) {
    if (edge_box_dist2(px, py, y, k) >= best) continue;
    const double d2 = seg_dist2(px, py, y, k);
    if (d2 < best) best = d2;
  }
  return best;
}

inline bool units_aabb_overlap(const UnitR2& x, const UnitR2& y) {
  return x.xmax >= y.xmin && x.xmin <= y.xmax &&
         x.ymax >= y.ymin && x.ymin <= y.ymax;
}

// 3. Directed vertex-to-set pass with the early break of Taha and Hanbury
// (2015, IEEE TPAMI 37(11)). Serves both the exact discrete method (`d_v ==
// nullptr`, matching GEOS) and the initial lower-bound pass of the continuous
// methods.
void vertex_pass_r2(const UnitR2& x, const UnitR2& y, bool check_pip,
                    double& max_d2, double& best_d,
                    std::vector<double>* d_v = nullptr) {
  const int nv = x.nv();
  const int ne = y.ne();
  if (d_v != nullptr) d_v->assign(nv, 0.0);

  int start = 0;
  for (int t = 0; t < nv; t++) {
    const int idx = x.order[t];
    const double px = x.v[idx].x;
    const double py = x.v[idx].y;
    if (check_pip && point_in_poly_r2(px, py, y)) {
      if (d_v != nullptr) (*d_v)[idx] = 0.0;
      continue;
    }
    double best2 = kInf;
    int best_k = start;
    bool closer = false;
    for (int s = 0; s < ne; s++) {
      int k = start + s;
      if (k >= ne) k -= ne;
      if (edge_box_dist2(px, py, y, k) >= best2) continue;
      const double d2 = seg_dist2(px, py, y, k);
      if (d2 <= max_d2) {
        start = k;
        if (d_v != nullptr) (*d_v)[idx] = std::sqrt(d2);
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
    if (d_v != nullptr) {
      const double d = std::sqrt(best2);
      (*d_v)[idx] = d;
      if (d > best_d) best_d = d;
    }
  }
}

double haus_discrete_r2(const UnitR2& x, const UnitR2& y) {
  double max_d2 = -1.0;
  double unused_best = 0.0;
  vertex_pass_r2(x, y, false, max_d2, unused_best);
  vertex_pass_r2(y, x, false, max_d2, unused_best);
  return max_d2 <= 0.0 ? 0.0 : std::sqrt(max_d2);
}

// 4. Continuous 1D (segment) and 2D (cell) branch-and-bound helpers.

bool same_unit_data_r2(const UnitR2& x, const UnitR2& y, bool check_poly) {
  if (x.nv() != y.nv() || x.a != y.a || x.b != y.b) return false;
  if (check_poly && x.e_is_poly != y.e_is_poly) return false;
  for (int i = 0; i < x.nv(); i++) {
    if (x.v[i].x != y.v[i].x || x.v[i].y != y.v[i].y) return false;
  }
  return true;
}

template <typename TaskT>
bool pop_next_task_r2(std::priority_queue<TaskT>& pq, int& iter,
                      double stop_lim, TaskT& out) {
  if (pq.empty()) return false;
  if ((++iter & 1023) == 0) Rcpp::checkUserInterrupt();
  out = pq.top();
  pq.pop();
  return out.upper > stop_lim;
}

inline bool edge_endpoints_outside_r2(double qx, double qy,
                                      const UnitR2& y, int ky,
                                      double keep_lim2) {
  if (y.b[ky] < 0) return true;
  const Pt2& A = y.v[y.a[ky]];
  const Pt2& B = y.v[y.b[ky]];
  return dist2_2d(qx, qy, A.x, A.y) > keep_lim2 &&
         dist2_2d(qx, qy, B.x, B.y) > keep_lim2;
}

// Filters active edges of `y` for child tasks via the triangle inequality.
std::vector<int> filter_active_edges_r2(double qx, double qy,
                                        const UnitR2& y,
                                        const std::vector<int>& active,
                                        const std::vector<double>& d_edge,
                                        double d_edge_min, double radius,
                                        double& keep_lim2_out) {
  const double keep_lim =
    d_edge_min + 2.0 * radius + 1e-12 * std::max(1.0, radius);
  const double keep_lim2 = keep_lim * keep_lim;
  keep_lim2_out = keep_lim2;
  std::vector<int> next;
  next.reserve(active.size());
  for (std::size_t i = 0; i < active.size(); i++) {
    const int ky = active[i];
    double d2 = d_edge[i];
    if (d2 == kInf) {
      if (edge_box_dist2(qx, qy, y, ky) > keep_lim2) continue;
      d2 = seg_dist2(qx, qy, y, ky);
    }
    if (d2 <= keep_lim2) next.push_back(ky);
  }
  return next;
}

// Clips convex polygon `poly` to nx * (p.x - ox) + ny * (p.y - oy) + c0 >= 0.
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

// Returns true when cell [x0, x1] x [y0, y1] cannot reach `target` anywhere,
// proved by intersecting the cell with linear upper-bound half-planes for the
// active edges of `y`.
bool ridge_prune_r2(double x0, double y0, double x1, double y1,
                    double cx, double cy, double r,
                    const UnitR2& y, const std::vector<int>& y_next,
                    double target) {
  if (target <= 0.0 || y_next.size() < 2) return false;
  std::vector<Pt2> poly = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}};
  const double r_lim = r * (1.0 + 1e-12);

  for (int ky : y_next) {
    const Pt2& A = y.v[y.a[ky]];
    double qx = A.x, qy = A.y;
    double d2;
    bool on_interior = false;
    if (y.b[ky] < 0) {
      d2 = dist2_2d(cx, cy, A.x, A.y);
    } else {
      const double t =
        ((cx - A.x) * y.dx[ky] + (cy - A.y) * y.dy[ky]) * y.inv_len2[ky];
      const double tc = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
      qx = A.x + tc * y.dx[ky];
      qy = A.y + tc * y.dy[ky];
      d2 = dist2_2d(cx, cy, qx, qy);
      on_interior = (t > 0.0 && t < 1.0);
    }
    const double dk = std::sqrt(d2);
    if (dk <= r_lim) continue;

    const double ux = (cx - qx) / dk;
    const double uy = (cy - qy) / dk;
    const double lim = dk + 2.0 * r + 1e-12 * std::max(1.0, r);
    if (on_interior && edge_endpoints_outside_r2(cx, cy, y, ky, lim * lim)) {
      poly = clip_linear_halfplane(poly, cx, cy, ux, uy, dk - target);
    } else {
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

struct SegTask {
  double ax, ay, bx, by;
  double da, db;
  double half_len;
  double upper;
  std::vector<int> active;
  bool operator<(const SegTask& o) const { return upper < o.upper; }
};

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
    if (x.b[k] < 0) continue;
    const int ia = x.a[k];
    const int ib = x.b[k];
    const double da = d_v[ia];
    const double db = d_v[ib];
    const double h = 0.5 * std::sqrt(x.dx[k] * x.dx[k] + x.dy[k] * x.dy[k]);
    const double upper = 0.5 * (da + db) + h;
    if (upper > best_d + tol_eff && h > min_h) {
      pq.push({x.v[ia].x, x.v[ia].y, x.v[ib].x, x.v[ib].y,
               da, db, h, upper, all_y});
    }
  }

  int iter = 0;
  SegTask cur;
  while (pop_next_task_r2(pq, iter, best_d + tol_eff, cur)) {
    const double h = cur.half_len;
    if (h <= min_h) continue;

    const double mx = 0.5 * (cur.ax + cur.bx);
    const double my = 0.5 * (cur.ay + cur.by);
    const bool check_early = (best_d > h);
    const double min_r2 = check_early ? (best_d - h) * (best_d - h) : 0.0;
    const double best_d2 = best_d * best_d;
    const bool in_y = check_pip && point_in_poly_r2(mx, my, y);

    double best_edge2 = kInf;
    bool pruned = false;
    const int na = static_cast<int>(cur.active.size());
    std::vector<double> d_edge(na, kInf);
    for (int i = 0; i < na; i++) {
      const int ky = cur.active[i];
      if (!check_early && edge_box_dist2(mx, my, y, ky) >= best_edge2) continue;
      const double d2 = seg_dist2(mx, my, y, ky);
      d_edge[i] = d2;
      if (d2 < best_edge2) best_edge2 = d2;
      if (check_early && d2 <= min_r2) {
        pruned = true;
        break;
      }
      if (d2 <= best_d2 &&
          seg_dist2(cur.ax, cur.ay, y, ky) <= best_d2 &&
          seg_dist2(cur.bx, cur.by, y, ky) <= best_d2) {
        pruned = true;
        break;
      }
    }
    if (pruned) continue;
    if (in_y && best_edge2 >= h * h * (1.0 + 1e-12)) continue;

    const double d_edge_min = std::sqrt(best_edge2);
    const double dm = in_y ? 0.0 : d_edge_min;
    if (dm > best_d) best_d = dm;

    double keep_lim2 = 0.0;
    std::vector<int> child_active = filter_active_edges_r2(
      mx, my, y, cur.active, d_edge, d_edge_min, h, keep_lim2
    );

    if (child_active.size() == 1) {
      const int ky = child_active[0];
      const bool safe_1edge =
        !is_polygon || !y.has_poly ||
        edge_endpoints_outside_r2(mx, my, y, ky, keep_lim2);
      if (safe_1edge) {
        const double da_ex = std::sqrt(
          point_to_set_dist2(cur.ax, cur.ay, y, check_pip, child_active)
        );
        const double db_ex = std::sqrt(
          point_to_set_dist2(cur.bx, cur.by, y, check_pip, child_active)
        );
        best_d = std::max(best_d, std::max(da_ex, db_ex));
        continue;
      }
    }

    const double child_h = 0.5 * h;
    if (child_h <= min_h) continue;
    const double u_left = 0.5 * (cur.da + dm + h);
    const double u_right = 0.5 * (dm + cur.db + h);
    const bool push_left = (u_left > best_d + tol_eff);
    const bool push_right = (u_right > best_d + tol_eff);
    if (push_left) {
      pq.push({cur.ax, cur.ay, mx, my, cur.da, dm, child_h, u_left,
               push_right ? child_active : std::move(child_active)});
    }
    if (push_right) {
      pq.push({mx, my, cur.bx, cur.by, dm, cur.db, child_h, u_right,
               std::move(child_active)});
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

// Evaluates a convex terminal cell against a single active edge of `y` when no
// polygon vertex of `x` lies inside the cell.
bool try_single_site_cell_r2(const CellTask& cur, double cx, double cy,
                             const UnitR2& x, const UnitR2& y,
                             const std::vector<int>& x_cur, bool in_x,
                             bool check_pip_y, bool in_y,
                             const std::vector<int>& y_next,
                             double keep_lim2, double& best_d) {
  if (y_next.size() != 1 || x_cur.size() > 1) return false;
  const int ky = y_next[0];
  if ((y.has_poly && !y.e_is_poly[ky] && in_y) ||
      !edge_endpoints_outside_r2(cx, cy, y, ky, keep_lim2)) {
    return false;
  }

  const std::vector<Pt2> box = {
    {cur.x0, cur.y0}, {cur.x1, cur.y0}, {cur.x1, cur.y1}, {cur.x0, cur.y1}
  };
  std::vector<Pt2> cand_pts;
  bool safe_x = false;

  if (x_cur.empty()) {
    if (in_x) {
      cand_pts = box;
      safe_x = true;
    }
  } else {
    const int kx = x_cur[0];
    const Pt2& A = x.v[x.a[kx]];
    const Pt2& B = x.v[x.b[kx]];
    const double r = cur.radius;
    const double pad = 1e-12 * std::max(1.0, r);
    const bool a_in = A.x >= cur.x0 - pad && A.x <= cur.x1 + pad &&
                      A.y >= cur.y0 - pad && A.y <= cur.y1 + pad;
    const bool b_in = B.x >= cur.x0 - pad && B.x <= cur.x1 + pad &&
                      B.y >= cur.y0 - pad && B.y <= cur.y1 + pad;
    if (!a_in && !b_in) {
      const double len = std::sqrt(x.dx[kx] * x.dx[kx] + x.dy[kx] * x.dy[kx]);
      if (len > 0.0) {
        const double nx = -x.dy[kx] / len;
        const double ny = x.dx[kx] / len;
        const std::vector<Pt2> pos =
          clip_linear_halfplane(box, A.x, A.y, nx, ny, 0.0);
        const std::vector<Pt2> neg =
          clip_linear_halfplane(box, A.x, A.y, -nx, -ny, 0.0);
        const double area_eps = 1e-24 * r * r;
        const bool pos_valid = poly_area2(pos) > area_eps;
        const bool neg_valid = poly_area2(neg) > area_eps;
        auto half_in_x = [&](const std::vector<Pt2>& hp, bool ok) {
          if (!ok) return false;
          const Pt2 c = poly_centroid(hp);
          return point_in_poly_r2(c.x, c.y, x);
        };
        const bool pos_in = half_in_x(pos, pos_valid);
        const bool neg_in = half_in_x(neg, neg_valid);
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
          safe_x = true;
        }
      }
    }
  }

  if (!safe_x) return false;
  for (const Pt2& pt : cand_pts) {
    const double d =
      std::sqrt(point_to_set_dist2(pt.x, pt.y, y, check_pip_y, y_next));
    if (d > best_d) best_d = d;
  }
  return true;
}

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
  for (int p = 0; p < static_cast<int>(x.parts.size()); p++) {
    if (x.parts[p].is_poly) poly_parts.push_back(p);
  }
  std::sort(poly_parts.begin(), poly_parts.end(), [&](int p1, int p2) {
    const PolyPart2D& a = x.parts[p1];
    const PolyPart2D& b = x.parts[p2];
    return (a.xmax - a.xmin) * (a.ymax - a.ymin) >
           (b.xmax - b.xmin) * (b.ymax - b.ymin);
  });

  std::vector<PolyPart2D> roots;
  std::priority_queue<CellTask> pq;
  for (int p : poly_parts) {
    const PolyPart2D& part = x.parts[p];
    bool covered = false;
    for (const PolyPart2D& rt : roots) {
      if (part.xmin >= rt.xmin && part.xmax <= rt.xmax &&
          part.ymin >= rt.ymin && part.ymax <= rt.ymax) {
        covered = true;
        break;
      }
    }
    if (covered) continue;
    roots.push_back(part);

    const double hx = 0.5 * (part.xmax - part.xmin);
    const double hy = 0.5 * (part.ymax - part.ymin);
    const double r0 = std::sqrt(hx * hx + hy * hy);
    if (r0 <= min_r) continue;
    pq.push({part.xmin, part.ymin, part.xmax, part.ymax,
             r0, kInf, x_poly_all, y_all});
  }

  int iter = 0;
  CellTask cur;
  while (pop_next_task_r2(pq, iter, best_d + tol_eff, cur)) {
    const double r = cur.radius;
    if (r <= min_r) continue;
    const double cx = 0.5 * (cur.x0 + cur.x1);
    const double cy = 0.5 * (cur.y0 + cur.y1);

    std::vector<int> x_cur;
    x_cur.reserve(cur.x_poly_edges.size());
    for (int k : cur.x_poly_edges) {
      const Pt2& A = x.v[x.a[k]];
      const Pt2& B = x.v[x.b[k]];
      if (seg_intersects_box(A.x, A.y, B.x, B.y,
                             cur.x0, cur.y0, cur.x1, cur.y1)) {
        x_cur.push_back(k);
      }
    }
    const bool in_x = point_in_poly_r2(cx, cy, x);
    if (x_cur.empty() && !in_x) continue;

    const bool in_y = check_pip_y && point_in_poly_r2(cx, cy, y);
    const bool check_early = (best_d > r);
    const double min_r2 = check_early ? (best_d - r) * (best_d - r) : 0.0;

    const int ny = static_cast<int>(cur.y_edges.size());
    std::vector<double> dy_edge(ny, kInf);
    double best_edge2 = kInf;
    double best_poly_edge2 = kInf;
    bool pruned = false;

    for (int i = 0; i < ny; i++) {
      const int ky = cur.y_edges[i];
      if (!check_early && !in_y &&
          edge_box_dist2(cx, cy, y, ky) >= best_edge2) {
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

    const double r_safe = r * (1.0 + 1e-12);
    if (in_y && best_poly_edge2 >= r_safe * r_safe) continue;

    const double d_edge_min = std::sqrt(best_edge2);
    const double dc = in_y ? 0.0 : d_edge_min;
    if (in_x && dc > best_d) best_d = dc;

    const double upper = dc + r;
    if (upper <= best_d + tol_eff) continue;

    double keep_lim2 = 0.0;
    std::vector<int> y_next = filter_active_edges_r2(
      cx, cy, y, cur.y_edges, dy_edge, d_edge_min, r, keep_lim2
    );

    if (try_single_site_cell_r2(cur, cx, cy, x, y, x_cur, in_x,
                                check_pip_y, in_y, y_next,
                                keep_lim2, best_d)) {
      continue;
    }

    if (best_d > 0.0 && y_next.size() >= 2 &&
        ridge_prune_r2(cur.x0, cur.y0, cur.x1, cur.y1, cx, cy, r,
                       y, y_next, best_d + tol_eff)) {
      continue;
    }

    const double child_r = 0.5 * r;
    if (child_r <= min_r) continue;
    const double qx[2][2] = {{cur.x0, cx}, {cx, cur.x1}};
    const double qy[2][2] = {{cur.y0, cy}, {cy, cur.y1}};
    for (int ix = 0; ix < 2; ix++) {
      for (int iy = 0; iy < 2; iy++) {
        pq.push({qx[ix][0], qy[iy][0], qx[ix][1], qy[iy][1],
                 child_r, upper, x_cur, y_next});
      }
    }
  }
}

double haus_r2_pair(const UnitR2& x_in, const UnitR2& y_in,
                    HausMethod method, double tol) {
  if (method == kDiscrete) {
    return haus_discrete_r2(x_in, y_in);
  }
  const bool is_polygon = (method == kPolygon);
  if (same_unit_data_r2(x_in, y_in, is_polygon)) {
    return 0.0;
  }

  const bool swap_order = should_swap_units(
    x_in.nv(), x_in.ne(), std::make_pair(x_in.xmin, x_in.ymin),
    y_in.nv(), y_in.ne(), std::make_pair(y_in.xmin, y_in.ymin)
  );
  const UnitR2& a = swap_order ? y_in : x_in;
  const UnitR2& b = swap_order ? x_in : y_in;

  const bool overlap = is_polygon && units_aabb_overlap(a, b);
  const bool check_pip_b = overlap && b.has_poly;
  const bool check_pip_a = overlap && a.has_poly;

  const double scale = std::max(1.0, std::max(a.diam, b.diam));
  const double max_coord = std::max(a.max_abs, b.max_abs);
  const double tol_eff =
    std::max(tol * scale, 64.0 * kEps * std::max(1.0, max_coord));
  const double min_step = 8.0 * kEps * std::max(1.0, max_coord);

  double max_d2 = -1.0;
  double best_d = 0.0;
  std::vector<double> dv_a, dv_b;
  vertex_pass_r2(a, b, check_pip_b, max_d2, best_d, &dv_a);
  vertex_pass_r2(b, a, check_pip_a, max_d2, best_d, &dv_b);

  directed_haus_1d_r2(a, b, is_polygon, check_pip_b, dv_a,
                      tol_eff, min_step, best_d);
  directed_haus_1d_r2(b, a, is_polygon, check_pip_a, dv_b,
                      tol_eff, min_step, best_d);

  if (is_polygon) {
    directed_haus_poly_2d_r2(a, b, check_pip_b, tol_eff, min_step, best_d);
    directed_haus_poly_2d_r2(b, a, check_pip_a, tol_eff, min_step, best_d);
  }
  return best_d;
}

}  // namespace

//' @title Planar Hausdorff distance matrix (internal)
//' @param poly_list list of matrices from .to_list_mat() or .to_polygon_list().
//' @param method one of "discrete", "boundary", or "polygon".
//' @param tol convergence tolerance for continuous methods.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd dist_haus_r2(const Rcpp::List& poly_list,
                             const std::string& method = "discrete",
                             double tol = 1e-8) {
  const HausMethod code = parse_haus_method(method, tol);
  const std::vector<UnitR2> u = make_unit_vec<UnitR2>(
    poly_list, [&](const Eigen::Ref<const Eigen::MatrixXd>& m) {
      return make_unit_r2(m, code);
    }
  );
  return sym_dist_matrix(
    static_cast<int>(u.size()),
    [&](int i) { return u[i].valid(); },
    [&](int i, int j) { return haus_r2_pair(u[i], u[j], code, tol); }
  );
}

//' @title Planar Hausdorff cross-distance matrix (internal)
//' @param poly_list,poly_pred lists of matrices from .to_list_mat() or
//'   .to_polygon_list().
//' @param method one of "discrete", "boundary", or "polygon".
//' @param tol convergence tolerance for continuous methods.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cdist_haus_r2(const Rcpp::List& poly_list,
                              const Rcpp::List& poly_pred,
                              const std::string& method = "discrete",
                              double tol = 1e-8) {
  const HausMethod code = parse_haus_method(method, tol);
  const auto builder = [&](const Eigen::Ref<const Eigen::MatrixXd>& m) {
    return make_unit_r2(m, code);
  };
  const std::vector<UnitR2> u1 = make_unit_vec<UnitR2>(poly_list, builder);
  const std::vector<UnitR2> u2 = make_unit_vec<UnitR2>(poly_pred, builder);
  return cross_dist_matrix(
    static_cast<int>(u1.size()),
    static_cast<int>(u2.size()),
    [&](int i) { return u1[i].valid(); },
    [&](int j) { return u2[j].valid(); },
    [&](int i, int j) { return haus_r2_pair(u1[i], u2[j], code, tol); }
  );
}
