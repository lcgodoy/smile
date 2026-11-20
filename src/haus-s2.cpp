#include <RcppEigen.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

// Spherical (S^2) Hausdorff distance between spatial units.
//
// Supports three methods:
// - "discrete": discrete vertex-to-arc Hausdorff distance (matches s2).
// - "boundary": continuous Hausdorff distance between boundaries/linework.
// - "polygon":  continuous Hausdorff distance between closed spatial sets
//               (polygons include 2D interiors excluding holes, following s2
//               smaller-region orientation per ring).
//
// Edges are great-circle arcs, as in s2. Every distance in this file is in
// meters, on a sphere of radius `earth_radius` (the default is
// s2::s2_earth_radius_meters()).
//
// For method = "polygon", each polygon unit must fit within a spherical cap of
// radius < kMaxPolyCapRad (0.45 * pi radians, ~81 degrees). Inside that cap,
// gnomonic projection maps great-circle arcs to exact straight line segments,
// so spherical point-in-polygon and 2D branch-and-bound run in exact planar
// tangent coordinates while distances remain exact S^2 great-circle angles.
// Polygon units exceeding the cap limit return NA with a warning, and
// continuous passes enforce a bounded iteration budget with a warning if tol
// is not reached.

namespace {

const double kPi = 3.14159265358979323846;
const double kMaxPolyCapRad = 0.45 * kPi;
const int kMaxBranchIter = 50000;
const double INF = std::numeric_limits<double>::infinity();
const double EPS = std::numeric_limits<double>::epsilon();
const double GOLDEN_RATIO = 0.6180339887498949;

inline double deg2rad(double degrees) {
  return degrees * kPi / 180.0;
}

inline double dot3(const double* u, const double* v) {
  return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

inline void cross3(const double* u, const double* v, double* out) {
  out[0] = u[1] * v[2] - u[2] * v[1];
  out[1] = u[2] * v[0] - u[0] * v[2];
  out[2] = u[0] * v[1] - u[1] * v[0];
}

inline bool normalize3(double* u) {
  const double n = std::sqrt(dot3(u, u));
  if (n <= 0.0 || !std::isfinite(n)) return false;
  u[0] /= n;
  u[1] /= n;
  u[2] /= n;
  return true;
}

// Squared distance between two 3D vectors.
inline double dist2(const double* u, const double* v) {
  const double d0 = u[0] - v[0];
  const double d1 = u[1] - v[1];
  const double d2 = u[2] - v[2];
  return d0 * d0 + d1 * d1 + d2 * d2;
}

inline double gc_angle_pts(const double* u, const double* v) {
  double c[3];
  cross3(u, v, c);
  return std::atan2(std::sqrt(dot3(c, c)), dot3(u, v));
}

struct PartS2 {
  int e_begin;
  int e_end;
  bool is_poly;
  double min_u, max_u, min_v, max_v;
};

// 1. A spatial unit, prepared once for all its distances.
//
// Input matrices have columns longitude, latitude (degrees), optionally a
// part id (see .to_list_mat() in R/dist-areal.R), and optionally an is_poly
// flag (1 for polygon rings, 0 for points/linestrings; see .to_polygon_list()).
// Without a part id, all rows are one ring, closed if its first and last rows
// differ.
//
// Each edge is a great-circle arc from vertex a to vertex b. It stores the
// unit normal n of its great circle and the tangents ta = n x A and
// tb = B x n. For a point, or for an edge shorter than about 0.6 mm, b is -1,
// and the distance to the edge is the distance to vertex a.
struct Unit {
  std::vector<double> v;          // 3 x (number of vertices), unit vectors
  std::vector<int> order;         // vertex order for the outer loop
  std::vector<int> a, b;          // edge endpoints
  std::vector<double> n, ta, tb;  // 3 x (number of edges)
  std::vector<double> arc_len;    // edge central angle (radians)
  std::vector<uint8_t> e_is_poly; // 1 if edge belongs to a polygon ring
  std::vector<PartS2> parts;      // per-part metadata

  bool has_poly;
  bool has_cap;
  bool poly_cap_exceeded;
  double cap_c[3];
  double cap_rad;
  double cos_cap;

  // Tangent-plane gnomonic basis and projected vertices (populated when
  // has_cap is true).
  double e1[3];
  double e2[3];
  std::vector<double> gu, gv;
  double poly_min_u, poly_max_u, poly_min_v, poly_max_v;

  int nv() const { return static_cast<int>(order.size()); }
  int ne() const { return static_cast<int>(a.size()); }
};

void add_edge(Unit& u, int i, int j, bool is_poly) {
  double nrm[3] = {0.0, 0.0, 0.0};
  double ta[3] = {0.0, 0.0, 0.0};
  double tb[3] = {0.0, 0.0, 0.0};
  double alen = 0.0;
  if (j >= 0) {
    const double* A = &u.v[3 * i];
    const double* B = &u.v[3 * j];
    cross3(A, B, nrm);
    const double len = std::sqrt(dot3(nrm, nrm));
    if (len < 1e-10) {
      j = -1;  // A and B are essentially the same point
    } else {
      for (int c = 0; c < 3; c++) nrm[c] /= len;
      cross3(nrm, A, ta);
      cross3(B, nrm, tb);
      alen = std::atan2(len, dot3(A, B));
    }
  }
  u.a.push_back(i);
  u.b.push_back(j);
  u.n.insert(u.n.end(), nrm, nrm + 3);
  u.ta.insert(u.ta.end(), ta, ta + 3);
  u.tb.insert(u.tb.end(), tb, tb + 3);
  u.arc_len.push_back(alen);
  u.e_is_poly.push_back(is_poly ? 1 : 0);
}

// A unit with a non-finite coordinate (for example, the NaN row that
// wk::wk_coords() gives for an EMPTY member of a MULTIPOINT) is empty, so its
// distances are NA, as in src/cr-dist.cpp.
Unit make_unit(const Eigen::Ref<const Eigen::MatrixXd>& poly) {
  Unit u;
  u.has_poly = false;
  u.has_cap = false;
  u.poly_cap_exceeded = false;
  u.cap_c[0] = 0.0; u.cap_c[1] = 0.0; u.cap_c[2] = 1.0;
  u.cap_rad = kPi;
  u.cos_cap = -1.0;
  u.e1[0] = 1.0; u.e1[1] = 0.0; u.e1[2] = 0.0;
  u.e2[0] = 0.0; u.e2[1] = 1.0; u.e2[2] = 0.0;
  u.poly_min_u = INF; u.poly_max_u = -INF;
  u.poly_min_v = INF; u.poly_max_v = -INF;

  const int nr = static_cast<int>(poly.rows());
  if (nr == 0) return u;
  for (int i = 0; i < nr; i++) {
    if (!std::isfinite(poly(i, 0)) || !std::isfinite(poly(i, 1))) return u;
  }
  u.v.resize(3 * nr);
  for (int i = 0; i < nr; i++) {
    const double lon = deg2rad(poly(i, 0));
    const double lat = deg2rad(poly(i, 1));
    u.v[3 * i] = std::cos(lat) * std::cos(lon);
    u.v[3 * i + 1] = std::cos(lat) * std::sin(lon);
    u.v[3 * i + 2] = std::sin(lat);
  }
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
    bool part_poly = has_poly_col ? (poly(p_start, 3) != 0.0) : (plen >= 3);

    // Check that a polygon ring has at least 3 distinct vertices.
    if (part_poly) {
      int distinct = 0;
      int first_idx = -1, last_idx = -1;
      for (int i = p_start; i < p_end; i++) {
        if (last_idx < 0 || dist2(&u.v[3 * last_idx], &u.v[3 * i]) > 1e-24) {
          if (first_idx < 0) first_idx = i;
          last_idx = i;
          distinct++;
        }
      }
      if (distinct >= 2 && dist2(&u.v[3 * last_idx], &u.v[3 * first_idx]) <= 1e-24) {
        distinct--;
      }
      if (distinct < 3) part_poly = false;
    }

    PartS2 part;
    part.e_begin = u.ne();
    part.is_poly = part_poly;
    part.min_u = INF; part.max_u = -INF;
    part.min_v = INF; part.max_v = -INF;

    if (plen == 1) {
      add_edge(u, p_start, -1, false);
    } else {
      for (int i = p_start; i + 1 < p_end; i++) {
        add_edge(u, i, i + 1, part_poly);
      }
      // Closing edge if needed (only for input without a part id)
      if (!has_part &&
          (poly(0, 0) != poly(nr - 1, 0) || poly(0, 1) != poly(nr - 1, 1))) {
        add_edge(u, nr - 1, 0, part_poly);
      }
    }
    part.e_end = u.ne();
    if (part.is_poly) u.has_poly = true;
    u.parts.push_back(part);
    p_start = p_end;
  }

  // Unit-level spherical bounding cap: test both vertex mean and 3D AABB
  // midpoint and keep whichever gives the tighter cap.
  double mean_c[3] = {0.0, 0.0, 0.0};
  double lo[3] = {INF, INF, INF};
  double hi[3] = {-INF, -INF, -INF};
  for (int i = 0; i < nr; i++) {
    for (int c = 0; c < 3; c++) {
      const double vc = u.v[3 * i + c];
      mean_c[c] += vc;
      if (vc < lo[c]) lo[c] = vc;
      if (vc > hi[c]) hi[c] = vc;
    }
  }
  double mid_c[3] = {
    0.5 * (lo[0] + hi[0]),
    0.5 * (lo[1] + hi[1]),
    0.5 * (lo[2] + hi[2])
  };
  double best_c[3] = {0.0, 0.0, 1.0};
  double best_r = kPi;
  double* cands[2] = {mean_c, mid_c};
  for (int t = 0; t < 2; t++) {
    if (!normalize3(cands[t])) continue;
    double r = 0.0;
    for (int i = 0; i < nr; i++) {
      const double a = gc_angle_pts(cands[t], &u.v[3 * i]);
      if (a > r) r = a;
    }
    if (r < best_r) {
      best_r = r;
      for (int c = 0; c < 3; c++) best_c[c] = cands[t][c];
    }
  }

  if (best_r < kMaxPolyCapRad) {
    u.has_cap = true;
    for (int c = 0; c < 3; c++) u.cap_c[c] = best_c[c];
    u.cap_rad = best_r;
    u.cos_cap = std::cos(best_r);

    // Precompute gnomonic tangent-plane coordinates for polygon units.
    if (u.has_poly) {
      const double ref[3] = {
        std::fabs(u.cap_c[2]) < 0.9 ? 0.0 : 1.0,
        0.0,
        std::fabs(u.cap_c[2]) < 0.9 ? 1.0 : 0.0
      };
      cross3(ref, u.cap_c, u.e1);
      normalize3(u.e1);
      cross3(u.cap_c, u.e1, u.e2);

      u.gu.resize(nr);
      u.gv.resize(nr);
      for (int i = 0; i < nr; i++) {
        const double* V = &u.v[3 * i];
        const double denom = dot3(V, u.cap_c);
        u.gu[i] = dot3(V, u.e1) / denom;
        u.gv[i] = dot3(V, u.e2) / denom;
      }
      for (std::size_t p = 0; p < u.parts.size(); p++) {
        PartS2& part = u.parts[p];
        if (!part.is_poly) continue;
        for (int k = part.e_begin; k < part.e_end; k++) {
          const int ends[2] = {u.a[k], u.b[k]};
          for (int m = 0; m < 2; m++) {
            const int idx = ends[m];
            if (idx < 0) continue;
            part.min_u = std::min(part.min_u, u.gu[idx]);
            part.max_u = std::max(part.max_u, u.gu[idx]);
            part.min_v = std::min(part.min_v, u.gv[idx]);
            part.max_v = std::max(part.max_v, u.gv[idx]);
          }
        }
        u.poly_min_u = std::min(u.poly_min_u, part.min_u);
        u.poly_max_u = std::max(u.poly_max_u, part.max_u);
        u.poly_min_v = std::min(u.poly_min_v, part.min_v);
        u.poly_max_v = std::max(u.poly_max_v, part.max_v);
      }
    }
  } else if (u.has_poly) {
    u.poly_cap_exceeded = true;
  }

  // A fixed scrambled order of the vertices (golden-ratio sequence). With
  // it, the largest distance grows fast, and the early break in
  // directed_haus() skips more vertices.
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
Rcpp::NumericMatrix as_unit_matrix(SEXP x) {
  Rcpp::NumericMatrix m(x);
  if (m.nrow() > 0 && m.ncol() < 2) {
    Rcpp::stop("Each unit must be a matrix with at least 2 columns.");
  }
  return m;
}

inline Eigen::Map<const Eigen::MatrixXd> as_eigen(const Rcpp::NumericMatrix& m) {
  return Eigen::Map<const Eigen::MatrixXd>(m.begin(), m.nrow(), m.ncol());
}

Unit unit_from_r(SEXP x) {
  const Rcpp::NumericMatrix m = as_unit_matrix(x);
  return make_unit(as_eigen(m));
}

std::vector<Unit> make_units(const Rcpp::List& polys) {
  std::vector<Unit> out;
  out.reserve(polys.size());
  for (R_xlen_t i = 0; i < polys.size(); i++) {
    out.push_back(unit_from_r(polys[i]));
  }
  return out;
}

// 2. Squared chord length from the unit vector P to edge k of y. The chord
// grows with the central angle, so it ranks distances without trigonometric
// functions, and it keeps full precision for small angles. `kind` returns the
// closest feature of the edge: 0 (vertex a), 1 (vertex b), or 2 (the inside
// of the arc).
inline double edge_chord2(const double* P, const Unit& y, int k, int& kind) {
  const double dA = dist2(P, &y.v[3 * y.a[k]]);
  if (y.b[k] < 0) {
    kind = 0;
    return dA;
  }
  if (dot3(&y.ta[3 * k], P) >= 0.0 && dot3(&y.tb[3 * k], P) >= 0.0) {
    const double* N = &y.n[3 * k];
    const double s = std::fabs(dot3(P, N));
    double c[3];
    cross3(N, P, c);
    kind = 2;
    return 2.0 * s * s / (1.0 + std::sqrt(dot3(c, c)));
  }
  const double dB = dist2(P, &y.v[3 * y.b[k]]);
  kind = dA <= dB ? 0 : 1;
  return std::min(dA, dB);
}

// 3. Central angle from P to edge k of y, for the closest feature `kind`.
inline double edge_angle(const double* P, const Unit& y, int k, int kind) {
  double c[3];
  if (kind == 2) {
    const double* N = &y.n[3 * k];
    cross3(N, P, c);
    return std::atan2(std::fabs(dot3(P, N)), std::sqrt(dot3(c, c)));
  }
  const double* Q = &y.v[3 * (kind == 0 ? y.a[k] : y.b[k])];
  cross3(P, Q, c);
  return std::atan2(std::sqrt(dot3(c, c)), dot3(P, Q));
}

// 4. Directed Hausdorff distance from the vertices of x to the edges of y,
// with the early break of Taha and Hanbury (2015, IEEE TPAMI 37(11)).
void directed_haus(const Unit& x, const Unit& y,
                   double& max_c2, double& max_angle) {
  const int ne = y.ne();
  int start = 0;
  for (int t = 0; t < x.nv(); t++) {
    const double* P = &x.v[3 * x.order[t]];
    double best = INF;
    int best_k = start;
    int best_kind = 0;
    bool closer = false;
    for (int s = 0; s < ne; s++) {
      int k = start + s;
      if (k >= ne) k -= ne;
      int kind;
      const double c2 = edge_chord2(P, y, k, kind);
      if (c2 <= max_c2) {
        start = k;
        closer = true;
        break;
      }
      if (c2 < best) {
        best = c2;
        best_k = k;
        best_kind = kind;
      }
    }
    if (closer) continue;
    start = best_k;
    max_c2 = best;
    max_angle = std::max(max_angle, edge_angle(P, y, best_k, best_kind));
  }
}

// 5. Discrete vertex-to-arc Hausdorff distance between two units, in meters.
double haus_units(const Unit& x, const Unit& y, double earth_radius) {
  if (x.nv() == 0 || y.nv() == 0) return NA_REAL;
  double max_c2 = -1.0;
  double max_angle = 0.0;
  directed_haus(x, y, max_c2, max_angle);
  directed_haus(y, x, max_c2, max_angle);
  return max_angle * earth_radius;
}

// ---- Continuous spherical 1D and 2D Hausdorff helpers ----------------------

// Exact spherical point-in-polygon in gnomonic tangent coordinates for any
// polygon unit fitting strictly inside its < 0.45 * pi cap.
bool point_in_poly_gnomonic(double pu, double pv, const Unit& u) {
  if (!u.has_poly || !u.has_cap) return false;
  if (pu < u.poly_min_u - 1e-12 || pu > u.poly_max_u + 1e-12 ||
      pv < u.poly_min_v - 1e-12 || pv > u.poly_max_v + 1e-12) {
    return false;
  }
  bool inside = false;
  const int np = static_cast<int>(u.parts.size());
  for (int p = 0; p < np; p++) {
    const PartS2& part = u.parts[p];
    if (!part.is_poly) continue;
    if (pu < part.min_u - 1e-12 || pu > part.max_u + 1e-12 ||
        pv < part.min_v - 1e-12 || pv > part.max_v + 1e-12) {
      continue;
    }
    for (int k = part.e_begin; k < part.e_end; k++) {
      if (u.b[k] < 0) continue;
      const double x1 = u.gu[u.a[k]], y1 = u.gv[u.a[k]];
      const double x2 = u.gu[u.b[k]], y2 = u.gv[u.b[k]];
      if ((y1 > pv) != (y2 > pv)) {
        const double x_int = x1 + (pv - y1) * (x2 - x1) / (y2 - y1);
        if (pu < x_int) inside = !inside;
      }
    }
  }
  return inside;
}

inline bool point_in_poly_s2(const double* P, const Unit& u) {
  if (!u.has_poly || !u.has_cap) return false;
  const double denom = dot3(P, u.cap_c);
  if (denom < u.cos_cap - 1e-12 || denom <= 0.0) return false;
  const double pu = dot3(P, u.e1) / denom;
  const double pv = dot3(P, u.e2) / denom;
  return point_in_poly_gnomonic(pu, pv, u);
}

inline bool units_caps_overlap_s2(const Unit& x, const Unit& y) {
  if (!x.has_cap || !y.has_cap) return true;
  return gc_angle_pts(x.cap_c, y.cap_c) <= x.cap_rad + y.cap_rad + 1e-12;
}

inline double angle_to_chord2(double ang) {
  if (ang <= 0.0) return 0.0;
  if (ang >= kPi) return 4.0;
  const double s = std::sin(0.5 * ang);
  return 4.0 * s * s;
}

// Exact 2D segment vs axis-aligned box intersection (Liang-Barsky), exact for
// great-circle arcs inside a unit's gnomonic tangent plane.
bool seg_intersects_box(double x1, double y1, double x2, double y2,
                        double bx0, double by0, double bx1, double by1) {
  if (std::max(x1, x2) < bx0 || std::min(x1, x2) > bx1 ||
      std::max(y1, y2) < by0 || std::min(y1, y2) > by1) {
    return false;
  }
  const double dx = x2 - x1;
  const double dy = y2 - y1;
  double t0 = 0.0, t1 = 1.0;
  const double p[4] = {-dx, dx, -dy, dy};
  const double q[4] = {x1 - bx0, bx1 - x1, y1 - by0, by1 - y1};
  for (int i = 0; i < 4; i++) {
    if (p[i] == 0.0) {
      if (q[i] < 0.0) return false;
    } else {
      const double r = q[i] / p[i];
      if (p[i] < 0.0) {
        if (r > t1) return false;
        if (r > t0) t0 = r;
      } else {
        if (r < t0) return false;
        if (r < t1) t1 = r;
      }
    }
  }
  return t0 <= t1;
}

// Fast exact check for identical undirected edge sets (handles reversed ring
// orientation, shifted ring starting vertices, and reordered parts).
bool same_undirected_edges_s2(const Unit& x, const Unit& y, bool check_poly) {
  const int ne = x.ne();
  if (ne != y.ne() || ne == 0 || x.has_poly != y.has_poly) return false;
  typedef std::array<double, 7> EdgeKey;
  auto build_keys = [&](const Unit& u) {
    std::vector<EdgeKey> keys(ne);
    for (int k = 0; k < ne; k++) {
      const double* A = &u.v[3 * u.a[k]];
      const double* B = (u.b[k] >= 0) ? &u.v[3 * u.b[k]] : A;
      const bool swap_ab =
        (B[0] < A[0]) ||
        (B[0] == A[0] && (B[1] < A[1] || (B[1] == A[1] && B[2] < A[2])));
      const double* P = swap_ab ? B : A;
      const double* Q = swap_ab ? A : B;
      keys[k] = {{P[0], P[1], P[2], Q[0], Q[1], Q[2],
                  check_poly ? static_cast<double>(u.e_is_poly[k]) : 0.0}};
    }
    std::sort(keys.begin(), keys.end());
    return keys;
  };
  return build_keys(x) == build_keys(y);
}

void vertex_pass_set_s2(const Unit& x, const Unit& y, bool is_polygon,
                        bool check_pip, double& best_ang,
                        std::vector<double>& ang_v) {
  const int nv = x.nv();
  const int ne = y.ne();
  ang_v.assign(nv, 0.0);
  double max_c2 = best_ang <= 0.0 ? -1.0 : angle_to_chord2(best_ang);

  int start = 0;
  for (int t = 0; t < nv; t++) {
    const int idx = x.order[t];
    const double* P = &x.v[3 * idx];
    if (is_polygon && check_pip && point_in_poly_s2(P, y)) {
      ang_v[idx] = 0.0;
      continue;
    }
    double best = INF;
    int best_k = start;
    int best_kind = 0;
    bool closer = false;
    for (int s = 0; s < ne; s++) {
      int k = start + s;
      if (k >= ne) k -= ne;
      int kind;
      const double c2 = edge_chord2(P, y, k, kind);
      if (c2 <= max_c2) {
        start = k;
        ang_v[idx] = edge_angle(P, y, k, kind);
        closer = true;
        break;
      }
      if (c2 < best) {
        best = c2;
        best_k = k;
        best_kind = kind;
      }
    }
    if (closer) continue;
    start = best_k;
    max_c2 = best;
    const double ang = edge_angle(P, y, best_k, best_kind);
    ang_v[idx] = ang;
    if (ang > best_ang) best_ang = ang;
  }
}

struct ArcTask {
  double A[3];
  double B[3];
  double da, db;
  double half_ang;
  double upper;
  std::vector<int> active;
  bool operator<(const ArcTask& o) const { return upper < o.upper; }
};

void directed_haus_1d_s2(const Unit& x, const Unit& y, bool is_polygon,
                         bool check_pip, const std::vector<double>& ang_v,
                         double tol_ang, double min_ang,
                         double& best_ang, double& unconv_upper) {
  const int ne_x = x.ne();
  const int ne_y = y.ne();
  if (ne_x == 0 || ne_y == 0) return;

  std::vector<int> all_y(ne_y);
  for (int k = 0; k < ne_y; k++) all_y[k] = k;

  std::priority_queue<ArcTask> pq;
  for (int k = 0; k < ne_x; k++) {
    if (x.b[k] < 0) continue;
    const int ia = x.a[k];
    const int ib = x.b[k];
    const double h = 0.5 * x.arc_len[k];
    const double da = ang_v[ia];
    const double db = ang_v[ib];
    const double upper = 0.5 * (da + db) + h;
    if (upper > best_ang + tol_ang && h > min_ang) {
      ArcTask task;
      for (int c = 0; c < 3; c++) {
        task.A[c] = x.v[3 * ia + c];
        task.B[c] = x.v[3 * ib + c];
      }
      task.da = da;
      task.db = db;
      task.half_ang = h;
      task.upper = upper;
      task.active = all_y;
      pq.push(std::move(task));
    }
  }

  int iter = 0;
  while (!pq.empty()) {
    if ((++iter & 1023) == 0) Rcpp::checkUserInterrupt();
    if (iter > kMaxBranchIter) {
      unconv_upper = std::max(unconv_upper, pq.top().upper);
      break;
    }
    ArcTask cur = pq.top();
    pq.pop();
    if (cur.upper <= best_ang + tol_ang) break;
    if (cur.half_ang <= min_ang) continue;

    double M[3] = {
      cur.A[0] + cur.B[0],
      cur.A[1] + cur.B[1],
      cur.A[2] + cur.B[2]
    };
    if (!normalize3(M)) continue;

    const double h = cur.half_ang;
    const bool check_early = (best_ang > h);
    const double early_c2 = check_early ? angle_to_chord2(best_ang - h) : -1.0;
    const double cover_lim = (best_ang + 0.5 * tol_ang < 0.49 * kPi && h < 0.45 * kPi)
      ? (std::sin(best_ang + 0.5 * tol_ang) * std::cos(h) - 1e-14)
      : -1.0;
    const bool in_y = is_polygon && check_pip && point_in_poly_s2(M, y);

    const int na = static_cast<int>(cur.active.size());
    std::vector<double> d_c2(na, INF);
    std::vector<int> d_kind(na, 0);
    double best_c2 = INF;
    int best_idx = -1;
    bool pruned = false;
    for (int i = 0; i < na; i++) {
      const int ky = cur.active[i];
      int kind;
      const double c2 = edge_chord2(M, y, ky, kind);
      d_c2[i] = c2;
      d_kind[i] = kind;
      if (c2 < best_c2) {
        best_c2 = c2;
        best_idx = i;
      }
      if (check_early && c2 <= early_c2) {
        pruned = true;
        break;
      }
      // Exact spherical arc-cover prune for coinciding or close arcs.
      if (y.b[ky] >= 0 && cover_lim > 0.0 &&
          dot3(&y.ta[3 * ky], cur.A) >= -1e-15 &&
          dot3(&y.tb[3 * ky], cur.A) >= -1e-15 &&
          dot3(&y.ta[3 * ky], cur.B) >= -1e-15 &&
          dot3(&y.tb[3 * ky], cur.B) >= -1e-15) {
        const double* N_ky = &y.n[3 * ky];
        if (std::max(std::fabs(dot3(N_ky, cur.A)),
                     std::fabs(dot3(N_ky, cur.B))) <= cover_lim) {
          pruned = true;
          break;
        }
      }
    }
    if (pruned || best_idx < 0) continue;

    const double d_edge_min = edge_angle(M, y, cur.active[best_idx], d_kind[best_idx]);
    if (in_y && d_edge_min >= h * (1.0 + 1e-12)) continue;

    const double dm = in_y ? 0.0 : d_edge_min;
    if (dm > best_ang) best_ang = dm;

    const double keep_lim = d_edge_min + 2.0 * h + 1e-12 * std::max(1.0, h);
    const double keep_c2 = angle_to_chord2(keep_lim);
    std::vector<int> child_active;
    child_active.reserve(na);
    for (int i = 0; i < na; i++) {
      if (d_c2[i] <= keep_c2) child_active.push_back(cur.active[i]);
    }

    const double child_h = 0.5 * h;
    if (child_h <= min_ang) continue;
    const double u_left = 0.5 * (cur.da + dm + h);
    if (u_left > best_ang + tol_ang) {
      ArcTask left;
      for (int c = 0; c < 3; c++) { left.A[c] = cur.A[c]; left.B[c] = M[c]; }
      left.da = cur.da; left.db = dm; left.half_ang = child_h;
      left.upper = u_left; left.active = child_active;
      pq.push(std::move(left));
    }
    const double u_right = 0.5 * (dm + cur.db + h);
    if (u_right > best_ang + tol_ang) {
      ArcTask right;
      for (int c = 0; c < 3; c++) { right.A[c] = M[c]; right.B[c] = cur.B[c]; }
      right.da = dm; right.db = cur.db; right.half_ang = child_h;
      right.upper = u_right; right.active = std::move(child_active);
      pq.push(std::move(right));
    }
  }
}

struct GnomonicBox {
  double u0, v0, u1, v1;
  double upper;
  std::vector<int> x_poly_edges;
  std::vector<int> y_edges;
  bool operator<(const GnomonicBox& o) const { return upper < o.upper; }
};

inline void gnomonic_to_s2(const Unit& x, double u, double v, double out[3]) {
  for (int c = 0; c < 3; c++) {
    out[c] = x.cap_c[c] + u * x.e1[c] + v * x.e2[c];
  }
  normalize3(out);
}

void directed_haus_poly_2d_s2(const Unit& x, const Unit& y, bool check_pip_y,
                              double tol_ang, double min_ang,
                              double& best_ang, double& unconv_upper) {
  if (!x.has_poly || !x.has_cap || y.ne() == 0) return;

  std::vector<int> x_poly_all;
  x_poly_all.reserve(x.ne());
  for (int k = 0; k < x.ne(); k++) {
    if (x.e_is_poly[k] && x.b[k] >= 0) x_poly_all.push_back(k);
  }
  if (x_poly_all.empty()) return;

  std::vector<int> y_all(y.ne());
  for (int k = 0; k < y.ne(); k++) y_all[k] = k;

  std::priority_queue<GnomonicBox> pq;
  GnomonicBox root;
  root.u0 = x.poly_min_u;
  root.v0 = x.poly_min_v;
  root.u1 = x.poly_max_u;
  root.v1 = x.poly_max_v;
  root.upper = kPi;
  root.x_poly_edges = std::move(x_poly_all);
  root.y_edges = std::move(y_all);
  pq.push(std::move(root));

  int iter = 0;
  while (!pq.empty()) {
    if ((++iter & 1023) == 0) Rcpp::checkUserInterrupt();
    if (iter > kMaxBranchIter) {
      unconv_upper = std::max(unconv_upper, pq.top().upper);
      break;
    }
    GnomonicBox cur = pq.top();
    pq.pop();
    if (cur.upper <= best_ang + tol_ang) break;

    const double um = 0.5 * (cur.u0 + cur.u1);
    const double vm = 0.5 * (cur.v0 + cur.v1);
    double C[3];
    gnomonic_to_s2(x, um, vm, C);

    double K[4][3];
    gnomonic_to_s2(x, cur.u0, cur.v0, K[0]);
    gnomonic_to_s2(x, cur.u1, cur.v0, K[1]);
    gnomonic_to_s2(x, cur.u1, cur.v1, K[2]);
    gnomonic_to_s2(x, cur.u0, cur.v1, K[3]);
    double r = 0.0;
    for (int j = 0; j < 4; j++) {
      const double a = gc_angle_pts(C, K[j]);
      if (a > r) r = a;
    }
    if (r <= min_ang) continue;

    // 1. Filter x polygon edges intersecting this gnomonic box.
    std::vector<int> x_cur;
    x_cur.reserve(cur.x_poly_edges.size());
    for (int k : cur.x_poly_edges) {
      if (seg_intersects_box(x.gu[x.a[k]], x.gv[x.a[k]],
                             x.gu[x.b[k]], x.gv[x.b[k]],
                             cur.u0, cur.v0, cur.u1, cur.v1)) {
        x_cur.push_back(k);
      }
    }
    const bool in_x = point_in_poly_gnomonic(um, vm, x);
    if (x_cur.empty() && !in_x) continue;

    // 2. Scan active y edges from spherical box center C.
    const bool valid_cap = (r < 0.49 * kPi);
    const bool in_y = check_pip_y && point_in_poly_s2(C, y);
    const bool check_early = valid_cap && (best_ang > r);
    const double early_c2 = check_early ? angle_to_chord2(best_ang - r) : -1.0;

    const int ny = static_cast<int>(cur.y_edges.size());
    std::vector<double> d_c2(ny, INF);
    std::vector<int> d_kind(ny, 0);
    double best_c2 = INF;
    int best_idx = -1;
    double best_poly_c2 = INF;
    int best_poly_idx = -1;
    bool pruned = false;

    for (int i = 0; i < ny; i++) {
      const int ky = cur.y_edges[i];
      int kind;
      const double c2 = edge_chord2(C, y, ky, kind);
      d_c2[i] = c2;
      d_kind[i] = kind;
      if (c2 < best_c2) {
        best_c2 = c2;
        best_idx = i;
      }
      if (y.e_is_poly[ky] && c2 < best_poly_c2) {
        best_poly_c2 = c2;
        best_poly_idx = i;
      }
      if (check_early && c2 <= early_c2) {
        pruned = true;
        break;
      }
    }
    if (pruned || best_idx < 0) continue;

    if (valid_cap && in_y && best_poly_idx >= 0) {
      const double poly_ang = edge_angle(C, y, cur.y_edges[best_poly_idx],
                                         d_kind[best_poly_idx]);
      if (poly_ang >= r * (1.0 + 1e-12)) continue;
    }

    const double d_edge_min = edge_angle(C, y, cur.y_edges[best_idx], d_kind[best_idx]);
    const double dc = in_y ? 0.0 : d_edge_min;
    if (in_x && dc > best_ang) best_ang = dc;

    const double upper = valid_cap ? std::min(kPi, dc + r) : kPi;
    if (upper <= best_ang + tol_ang) continue;

    std::vector<int> y_next;
    if (valid_cap) {
      const double keep_lim = d_edge_min + 2.0 * r + 1e-12 * std::max(1.0, r);
      const double keep_c2 = angle_to_chord2(keep_lim);
      y_next.reserve(ny);
      for (int i = 0; i < ny; i++) {
        if (d_c2[i] <= keep_c2) y_next.push_back(cur.y_edges[i]);
      }
    } else {
      y_next = cur.y_edges;
    }

    // 3. Subdivide gnomonic box along its longer side.
    GnomonicBox c1, c2;
    c1.u0 = cur.u0; c1.u1 = cur.u1; c1.v0 = cur.v0; c1.v1 = cur.v1;
    c2.u0 = cur.u0; c2.u1 = cur.u1; c2.v0 = cur.v0; c2.v1 = cur.v1;
    if (cur.u1 - cur.u0 >= cur.v1 - cur.v0) {
      c1.u1 = um;
      c2.u0 = um;
    } else {
      c1.v1 = vm;
      c2.v0 = vm;
    }
    c1.upper = upper;
    c2.upper = upper;
    c1.x_poly_edges = x_cur;
    c1.y_edges = y_next;
    c2.x_poly_edges = std::move(x_cur);
    c2.y_edges = std::move(y_next);
    pq.push(std::move(c1));
    pq.push(std::move(c2));
  }
}

double haus_s2_method_units(const Unit& x_in, const Unit& y_in,
                            int method_code, double tol, double earth_radius,
                            int& n_unconv, double& max_gap_m) {
  if (x_in.nv() == 0 || y_in.nv() == 0) return NA_REAL;
  if (method_code == 0) {
    return haus_units(x_in, y_in, earth_radius);
  }
  if (method_code == 2 && (x_in.poly_cap_exceeded || y_in.poly_cap_exceeded)) {
    return NA_REAL;
  }
  const bool is_polygon = (method_code == 2);
  if (same_undirected_edges_s2(x_in, y_in, is_polygon)) {
    return 0.0;
  }

  // Canonicalize unit order for bit-for-bit symmetry between (x, y) and (y, x).
  const bool swap_order =
    (y_in.nv() < x_in.nv()) ||
    (y_in.nv() == x_in.nv() && (y_in.ne() < x_in.ne() ||
      (y_in.ne() == x_in.ne() && std::make_pair(y_in.cap_rad, y_in.v[0]) <
                                 std::make_pair(x_in.cap_rad, x_in.v[0]))));
  const Unit& a = swap_order ? y_in : x_in;
  const Unit& b = swap_order ? x_in : y_in;

  const bool overlap = is_polygon && units_caps_overlap_s2(a, b);
  const bool check_pip_b = overlap && b.has_poly;
  const bool check_pip_a = overlap && a.has_poly;

  const double diam_a = a.has_cap ? 2.0 * a.cap_rad * earth_radius : kPi * earth_radius;
  const double diam_b = b.has_cap ? 2.0 * b.cap_rad * earth_radius : kPi * earth_radius;
  const double scale_m = std::max(1.0, std::max(diam_a, diam_b));
  const double tol_m = std::max(tol * scale_m, 64.0 * EPS * earth_radius);
  const double tol_ang = tol_m / earth_radius;
  const double min_ang = 8.0 * EPS;

  double best_ang = 0.0;
  double unconv_upper = 0.0;
  std::vector<double> ang_a, ang_b;
  vertex_pass_set_s2(a, b, is_polygon, check_pip_b, best_ang, ang_a);
  vertex_pass_set_s2(b, a, is_polygon, check_pip_a, best_ang, ang_b);

  directed_haus_1d_s2(a, b, is_polygon, check_pip_b, ang_a,
                      tol_ang, min_ang, best_ang, unconv_upper);
  directed_haus_1d_s2(b, a, is_polygon, check_pip_a, ang_b,
                      tol_ang, min_ang, best_ang, unconv_upper);

  if (is_polygon) {
    directed_haus_poly_2d_s2(a, b, check_pip_b, tol_ang, min_ang,
                             best_ang, unconv_upper);
    directed_haus_poly_2d_s2(b, a, check_pip_a, tol_ang, min_ang,
                             best_ang, unconv_upper);
  }

  if (unconv_upper > best_ang + tol_ang) {
    n_unconv++;
    max_gap_m = std::max(max_gap_m, (unconv_upper - best_ang) * earth_radius);
  }
  return best_ang * earth_radius;
}

void warn_poly_cap_exceeded(const std::vector<Unit>& u, const char* arg) {
  std::string idx_list;
  for (std::size_t i = 0; i < u.size(); i++) {
    if (u[i].poly_cap_exceeded) {
      if (!idx_list.empty()) idx_list += ", ";
      idx_list += std::to_string(i + 1);
    }
  }
  if (!idx_list.empty()) {
    Rcpp::warning(
      "These polygon units of %s exceed the spherical cap limit "
      "(radius >= 0.45 * pi * R): %s. Returning NA for pairs involving them; "
      "use method = \"discrete\" or \"boundary\" for global-scale polygons.",
      arg, idx_list.c_str()
    );
  }
}

void warn_unconverged(int n_unconv, double max_gap_m) {
  if (n_unconv > 0) {
    Rcpp::warning(
      "Continuous spherical Hausdorff distance did not reach tol for %d pair(s) "
      "within the iteration limit (largest remaining bound gap: %.3g m). "
      "Consider increasing tol or using method = \"discrete\".",
      n_unconv, max_gap_m
    );
  }
}

int parse_method_code(const std::string& method) {
  if (method == "discrete") return 0;
  if (method == "boundary") return 1;
  if (method == "polygon") return 2;
  Rcpp::stop("Unknown method: must be \"discrete\", \"boundary\", or \"polygon\".");
}

}  // namespace

//' @title Spherical Hausdorff distance between two units (internal)
//' @param x,y matrices with the columns longitude, latitude (degrees), and
//'   optionally a part id.
//' @param earth_radius radius of the sphere, in meters.
//' @noRd
// [[Rcpp::export]]
double haus_s2_aux_continuous(SEXP x, SEXP y,
                              double earth_radius = 6371010.0) {
  return haus_units(unit_from_r(x), unit_from_r(y), earth_radius);
}

//' @title Spherical Hausdorff distance matrix (internal)
//' @param poly_list list of matrices from .to_list_mat().
//' @param earth_radius radius of the sphere, in meters.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd dist_haus_s2(const Rcpp::List& poly_list,
                             double earth_radius = 6371010.0) {
  const int n = poly_list.size();
  const std::vector<Unit> u = make_units(poly_list);
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n);

  for (int i = 0; i < n; i++) {
    Rcpp::checkUserInterrupt();
    if (u[i].nv() == 0) out(i, i) = NA_REAL;
    for (int j = i + 1; j < n; j++) {
      out(i, j) = haus_units(u[i], u[j], earth_radius);
      out(j, i) = out(i, j);  // The matrix is symmetric
    }
  }
  return out;
}

//' @title Spherical Hausdorff distances from units to points (internal)
//' @param poly_list list of matrices from .to_list_mat().
//' @param pred_mat matrix with one point (longitude, latitude) for each row.
//' @param earth_radius radius of the sphere, in meters.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cdist_haus_s2_mp(const Rcpp::List& poly_list, SEXP pred_mat,
                                 double earth_radius = 6371010.0) {
  const Rcpp::NumericMatrix pm = as_unit_matrix(pred_mat);
  const Eigen::Map<const Eigen::MatrixXd> pred = as_eigen(pm);
  const int n1 = poly_list.size();
  const int n2 = static_cast<int>(pred.rows());
  const std::vector<Unit> u = make_units(poly_list);
  std::vector<Unit> p;
  p.reserve(n2);
  for (int j = 0; j < n2; j++) {
    p.push_back(make_unit(Eigen::MatrixXd(pred.row(j))));
  }
  Eigen::MatrixXd out(n1, n2);

  for (int i = 0; i < n1; i++) {
    Rcpp::checkUserInterrupt();
    for (int j = 0; j < n2; j++) {
      out(i, j) = haus_units(u[i], p[j], earth_radius);
    }
  }
  return out;
}

//' @title Spherical Hausdorff cross-distance matrix (internal)
//' @param poly_list,poly_pred lists of matrices from .to_list_mat().
//' @param earth_radius radius of the sphere, in meters.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cdist_haus_s2_lists(const Rcpp::List& poly_list,
                                    const Rcpp::List& poly_pred,
                                    double earth_radius = 6371010.0) {
  const int n1 = poly_list.size();
  const int n2 = poly_pred.size();
  const std::vector<Unit> u1 = make_units(poly_list);
  const std::vector<Unit> u2 = make_units(poly_pred);
  Eigen::MatrixXd out(n1, n2);

  for (int i = 0; i < n1; i++) {
    Rcpp::checkUserInterrupt();
    for (int j = 0; j < n2; j++) {
      out(i, j) = haus_units(u1[i], u2[j], earth_radius);
    }
  }
  return out;
}

//' @title Spherical Hausdorff distance matrix with method selection (internal)
//' @param poly_list list of matrices from .to_list_mat() or .to_polygon_list().
//' @param method one of "discrete", "boundary", or "polygon".
//' @param tol convergence tolerance for continuous methods.
//' @param earth_radius radius of the sphere, in meters.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd dist_haus_s2_method(const Rcpp::List& poly_list,
                                    const std::string& method = "discrete",
                                    double tol = 1e-8,
                                    double earth_radius = 6371010.0) {
  const int method_code = parse_method_code(method);
  if (method_code != 0 && (!std::isfinite(tol) || tol <= 0.0)) {
    Rcpp::stop("tol must be a single finite positive number.");
  }
  if (method_code == 0) return dist_haus_s2(poly_list, earth_radius);

  const int n = poly_list.size();
  const std::vector<Unit> u = make_units(poly_list);
  if (method_code == 2) warn_poly_cap_exceeded(u, "x");

  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n);
  int n_unconv = 0;
  double max_gap_m = 0.0;

  for (int i = 0; i < n; i++) {
    Rcpp::checkUserInterrupt();
    if (u[i].nv() == 0 || (method_code == 2 && u[i].poly_cap_exceeded)) {
      out(i, i) = NA_REAL;
    }
    for (int j = i + 1; j < n; j++) {
      const double d = haus_s2_method_units(u[i], u[j], method_code, tol,
                                            earth_radius, n_unconv, max_gap_m);
      out(i, j) = d;
      out(j, i) = d;
    }
  }
  warn_unconverged(n_unconv, max_gap_m);
  return out;
}

//' @title Spherical Hausdorff cross-distance matrix with method selection (internal)
//' @param poly_list,poly_pred lists of matrices from .to_list_mat() or
//'   .to_polygon_list().
//' @param method one of "discrete", "boundary", or "polygon".
//' @param tol convergence tolerance for continuous methods.
//' @param earth_radius radius of the sphere, in meters.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cdist_haus_s2_lists_method(const Rcpp::List& poly_list,
                                           const Rcpp::List& poly_pred,
                                           const std::string& method = "discrete",
                                           double tol = 1e-8,
                                           double earth_radius = 6371010.0) {
  const int method_code = parse_method_code(method);
  if (method_code != 0 && (!std::isfinite(tol) || tol <= 0.0)) {
    Rcpp::stop("tol must be a single finite positive number.");
  }
  if (method_code == 0) return cdist_haus_s2_lists(poly_list, poly_pred, earth_radius);

  const int n1 = poly_list.size();
  const int n2 = poly_pred.size();
  const std::vector<Unit> u1 = make_units(poly_list);
  const std::vector<Unit> u2 = make_units(poly_pred);
  if (method_code == 2) {
    warn_poly_cap_exceeded(u1, "x");
    warn_poly_cap_exceeded(u2, "y");
  }

  Eigen::MatrixXd out(n1, n2);
  int n_unconv = 0;
  double max_gap_m = 0.0;

  for (int i = 0; i < n1; i++) {
    Rcpp::checkUserInterrupt();
    for (int j = 0; j < n2; j++) {
      out(i, j) = haus_s2_method_units(u1[i], u2[j], method_code, tol,
                                       earth_radius, n_unconv, max_gap_m);
    }
  }
  warn_unconverged(n_unconv, max_gap_m);
  return out;
}
