#include "haus-common.h"
#include <array>
#include <cstdint>
#include <queue>

// [[Rcpp::depends(RcppEigen)]]

// Spherical (S^2) Hausdorff distance between spatial units.
//
// Supports three methods:
// - "discrete": discrete vertex-to-arc Hausdorff distance (matches s2).
// - "boundary": continuous Hausdorff distance between boundaries and linework.
// - "polygon":  continuous Hausdorff distance between closed spatial sets
//               (polygons include 2D interiors excluding holes, following s2
//               smaller-region orientation per ring).
//
// Edges are great-circle arcs, as in s2. Every distance in this file is in
// meters on a sphere of radius `earth_radius`.
//
// For method = "polygon", each polygon unit must fit inside a spherical cap of
// radius below kMaxPolyCapRad (0.45 * pi radians, about 81 degrees). Inside
// that cap, gnomonic projection maps great-circle arcs to exact straight line
// segments, so spherical point-in-polygon and 2D branch-and-bound use exact
// planar tangent coordinates while distances stay exact S^2 central angles.

namespace {

using smile_geom::Pt2;
using smile_geom::Vec3;
using smile_geom::angle_to_chord2;
using smile_geom::cross3;
using smile_geom::cross_dist_matrix;
using smile_geom::dist2_3d;
using smile_geom::dot3;
using smile_geom::gc_angle;
using smile_geom::golden_order;
using smile_geom::kEps;
using smile_geom::kInf;
using smile_geom::kPi;
using smile_geom::ll_to_xyz;
using smile_geom::norm3;
using smile_geom::normalize3;
using smile_geom::sym_dist_matrix;
using smile_geom::tangent_basis;
using smile_geom::vec3_less;

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

const double kMaxPolyCapRad = 0.45 * kPi;
const int kMaxBranchIter = 50000;

// 1. A spatial unit on S^2, prepared once for all its pairwise distances.
//
// Each edge is a great-circle arc from vertex a to vertex b. It stores the
// unit normal n of its great circle and the inward tangents ta = n x A and
// tb = B x n. For a point, or for an edge shorter than about 0.6 mm on Earth,
// b is -1, and the distance to the edge is the distance to vertex a.
struct UnitS2 {
  std::vector<Vec3> v;            // Unit vectors on S^2
  std::vector<int> order;         // Golden-ratio vertex order
  std::vector<int> a, b;          // Edge endpoints (b == -1 for points)
  std::vector<Vec3> n, ta, tb;    // Great-circle normal and arc tangents
  std::vector<double> arc_len;    // Edge central angle (continuous methods)
  std::vector<uint8_t> e_is_poly; // 1 if edge belongs to a polygon ring
  std::vector<PolyPart2D> parts;  // Per-part gnomonic metadata ("polygon")

  bool has_poly;
  bool has_cap;
  bool poly_cap_exceeded;
  Vec3 cap_c;
  double cap_rad;
  double cos_cap;

  Vec3 e1, e2;
  std::vector<Pt2> g;             // Gnomonic tangent-plane vertices
  double poly_min_u, poly_max_u, poly_min_v, poly_max_v;

  int nv() const { return static_cast<int>(order.size()); }
  int ne() const { return static_cast<int>(a.size()); }
  bool valid_for(HausMethod method) const {
    return nv() > 0 && !(method == kPolygon && poly_cap_exceeded);
  }
};

void add_edge_s2(UnitS2& u, int i, int j, bool is_poly, bool need_arc_len) {
  Vec3 nrm = {0.0, 0.0, 0.0};
  Vec3 ta = {0.0, 0.0, 0.0};
  Vec3 tb = {0.0, 0.0, 0.0};
  double alen = 0.0;
  if (j >= 0) {
    const Vec3& A = u.v[i];
    const Vec3& B = u.v[j];
    nrm = cross3(A, B);
    const double len = norm3(nrm);
    if (len < 1e-10) {
      j = -1;  // Endpoints are closer than about 0.6 mm on Earth
    } else {
      nrm.x /= len;
      nrm.y /= len;
      nrm.z /= len;
      ta = cross3(nrm, A);
      tb = cross3(B, nrm);
      if (need_arc_len) alen = std::atan2(len, dot3(A, B));
    }
  }
  u.a.push_back(i);
  u.b.push_back(j);
  u.n.push_back(nrm);
  u.ta.push_back(ta);
  u.tb.push_back(tb);
  if (need_arc_len) u.arc_len.push_back(alen);
  u.e_is_poly.push_back(is_poly ? 1 : 0);
}

bool has_three_distinct_vertices(const std::vector<Vec3>& v,
                                 int begin, int end) {
  int distinct = 0;
  int first_idx = -1;
  int last_idx = -1;
  for (int i = begin; i < end; i++) {
    if (last_idx < 0 || dist2_3d(v[last_idx], v[i]) > 1e-24) {
      if (first_idx < 0) first_idx = i;
      last_idx = i;
      distinct++;
    }
  }
  if (distinct >= 2 && dist2_3d(v[last_idx], v[first_idx]) <= 1e-24) {
    distinct--;
  }
  return distinct >= 3;
}

// Tests both the vertex mean and the 3D bounding-box midpoint as candidate cap
// centers and keeps whichever gives the smaller angular radius.
void compute_bounding_cap_s2(UnitS2& u) {
  const int nr = static_cast<int>(u.v.size());
  Vec3 mean_c = {0.0, 0.0, 0.0};
  Vec3 lo = {kInf, kInf, kInf};
  Vec3 hi = {-kInf, -kInf, -kInf};
  for (int i = 0; i < nr; i++) {
    const Vec3& p = u.v[i];
    mean_c.x += p.x; mean_c.y += p.y; mean_c.z += p.z;
    lo.x = std::min(lo.x, p.x); hi.x = std::max(hi.x, p.x);
    lo.y = std::min(lo.y, p.y); hi.y = std::max(hi.y, p.y);
    lo.z = std::min(lo.z, p.z); hi.z = std::max(hi.z, p.z);
  }
  Vec3 cands[2] = {
    mean_c,
    {0.5 * (lo.x + hi.x), 0.5 * (lo.y + hi.y), 0.5 * (lo.z + hi.z)}
  };

  Vec3 best_c = {0.0, 0.0, 1.0};
  double best_r = kPi;
  for (int t = 0; t < 2; t++) {
    if (!normalize3(cands[t])) continue;
    double r = 0.0;
    for (int i = 0; i < nr; i++) {
      r = std::max(r, gc_angle(cands[t], u.v[i]));
    }
    if (r < best_r) {
      best_r = r;
      best_c = cands[t];
    }
  }
  if (best_r < kMaxPolyCapRad) {
    u.has_cap = true;
    u.cap_c = best_c;
    u.cap_rad = best_r;
    u.cos_cap = std::cos(best_r);
  }
}

void project_gnomonic_s2(UnitS2& u) {
  tangent_basis(u.cap_c, u.e1, u.e2);
  const int nr = static_cast<int>(u.v.size());
  u.g.resize(nr);
  for (int i = 0; i < nr; i++) {
    const double denom = dot3(u.v[i], u.cap_c);
    u.g[i] = {dot3(u.v[i], u.e1) / denom, dot3(u.v[i], u.e2) / denom};
  }
  for (PolyPart2D& part : u.parts) {
    if (!part.is_poly) continue;
    for (int k = part.e_begin; k < part.e_end; k++) {
      const int ends[2] = {u.a[k], u.b[k]};
      for (int m = 0; m < 2; m++) {
        const int idx = ends[m];
        if (idx < 0) continue;
        part.xmin = std::min(part.xmin, u.g[idx].x);
        part.xmax = std::max(part.xmax, u.g[idx].x);
        part.ymin = std::min(part.ymin, u.g[idx].y);
        part.ymax = std::max(part.ymax, u.g[idx].y);
      }
    }
    u.poly_min_u = std::min(u.poly_min_u, part.xmin);
    u.poly_max_u = std::max(u.poly_max_u, part.xmax);
    u.poly_min_v = std::min(u.poly_min_v, part.ymin);
    u.poly_max_v = std::max(u.poly_max_v, part.ymax);
  }
}

UnitS2 make_unit_s2(const Eigen::Ref<const Eigen::MatrixXd>& poly,
                    HausMethod method) {
  UnitS2 u;
  u.has_poly = false;
  u.has_cap = false;
  u.poly_cap_exceeded = false;
  u.cap_c = {0.0, 0.0, 1.0};
  u.cap_rad = kPi;
  u.cos_cap = -1.0;
  u.e1 = {1.0, 0.0, 0.0};
  u.e2 = {0.0, 1.0, 0.0};
  u.poly_min_u = u.poly_min_v = kInf;
  u.poly_max_u = u.poly_max_v = -kInf;

  if (!has_valid_coords(poly)) return u;
  const int nr = static_cast<int>(poly.rows());
  u.v.resize(nr);
  for (int i = 0; i < nr; i++) {
    u.v[i] = ll_to_xyz(poly(i, 0), poly(i, 1));
  }

  const bool need_arc_len = (method != kDiscrete);
  const bool need_poly = (method == kPolygon);
  const std::vector<RawPart> raw_parts = extract_raw_parts(poly);

  for (const RawPart& rp : raw_parts) {
    bool part_poly =
      need_poly && rp.is_poly &&
      has_three_distinct_vertices(u.v, rp.begin, rp.end);

    const int e_begin = u.ne();
    if (rp.end - rp.begin == 1) {
      add_edge_s2(u, rp.begin, -1, false, need_arc_len);
    } else {
      for (int i = rp.begin; i + 1 < rp.end; i++) {
        add_edge_s2(u, i, i + 1, part_poly, need_arc_len);
      }
      if (rp.close_ring) {
        add_edge_s2(u, rp.end - 1, rp.begin, part_poly, need_arc_len);
      }
    }

    if (need_poly) {
      if (part_poly) u.has_poly = true;
      u.parts.push_back(
        {e_begin, u.ne(), part_poly, kInf, kInf, -kInf, -kInf}
      );
    }
  }

  if (method != kDiscrete) {
    compute_bounding_cap_s2(u);
    if (need_poly && u.has_poly) {
      if (u.has_cap) {
        project_gnomonic_s2(u);
      } else {
        u.poly_cap_exceeded = true;
      }
    }
  }

  u.order = golden_order(nr);
  return u;
}

// 2. Point-to-arc distance primitives on S^2.
//
// Squared chord length increases strictly with central angle on [0, pi], so
// inner loops rank and prune edges without inverse trigonometric calls. `kind`
// identifies the closest feature: 0 (vertex a), 1 (vertex b), or 2 (arc
// interior).
inline double edge_chord2(const Vec3& P, const UnitS2& y, int k, int& kind) {
  const double dA = dist2_3d(P, y.v[y.a[k]]);
  if (y.b[k] < 0) {
    kind = 0;
    return dA;
  }
  if (dot3(y.ta[k], P) >= 0.0 && dot3(y.tb[k], P) >= 0.0) {
    const Vec3& N = y.n[k];
    const double s = std::fabs(dot3(P, N));
    const Vec3 c = cross3(N, P);
    kind = 2;
    return 2.0 * s * s / (1.0 + norm3(c));
  }
  const double dB = dist2_3d(P, y.v[y.b[k]]);
  kind = dA <= dB ? 0 : 1;
  return std::min(dA, dB);
}

inline double edge_angle(const Vec3& P, const UnitS2& y, int k, int kind) {
  if (kind == 2) {
    const Vec3& N = y.n[k];
    const Vec3 c = cross3(N, P);
    return std::atan2(std::fabs(dot3(P, N)), norm3(c));
  }
  return gc_angle(P, y.v[kind == 0 ? y.a[k] : y.b[k]]);
}

// 3. Spherical point-in-polygon via gnomonic tangent coordinates.
inline bool point_in_poly_gnomonic(double pu, double pv, const UnitS2& u) {
  if (!u.has_poly || !u.has_cap) return false;
  return point_in_poly_2d(pu, pv, u.g, u.a, u.b, u.parts,
                          u.poly_min_u, u.poly_min_v,
                          u.poly_max_u, u.poly_max_v, 1e-12);
}

inline bool point_in_poly_s2(const Vec3& P, const UnitS2& u) {
  if (!u.has_poly || !u.has_cap) return false;
  const double denom = dot3(P, u.cap_c);
  if (denom < u.cos_cap - 1e-12 || denom <= 0.0) return false;
  return point_in_poly_gnomonic(dot3(P, u.e1) / denom,
                                dot3(P, u.e2) / denom, u);
}

inline bool units_caps_overlap_s2(const UnitS2& x, const UnitS2& y) {
  if (!x.has_cap || !y.has_cap) return true;
  return gc_angle(x.cap_c, y.cap_c) <= x.cap_rad + y.cap_rad + 1e-12;
}

// 4. Directed vertex-to-set pass with the early break of Taha and Hanbury
// (2015, IEEE TPAMI 37(11)). Serves both the exact discrete method (`ang_v ==
// nullptr`) and the initial lower-bound pass of the continuous methods.
void vertex_pass_s2(const UnitS2& x, const UnitS2& y, bool check_pip,
                    double& max_c2, double& best_ang,
                    std::vector<double>* ang_v = nullptr) {
  const int nv = x.nv();
  const int ne = y.ne();
  if (ang_v != nullptr) ang_v->assign(nv, 0.0);

  int start = 0;
  for (int t = 0; t < nv; t++) {
    const int idx = x.order[t];
    const Vec3& P = x.v[idx];
    if (check_pip && point_in_poly_s2(P, y)) {
      if (ang_v != nullptr) (*ang_v)[idx] = 0.0;
      continue;
    }
    double best = kInf;
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
        if (ang_v != nullptr) (*ang_v)[idx] = edge_angle(P, y, k, kind);
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
    if (ang_v != nullptr) (*ang_v)[idx] = ang;
    if (ang > best_ang) best_ang = ang;
  }
}

double haus_discrete_s2(const UnitS2& x, const UnitS2& y, double earth_radius) {
  double max_c2 = -1.0;
  double best_ang = 0.0;
  vertex_pass_s2(x, y, false, max_c2, best_ang);
  vertex_pass_s2(y, x, false, max_c2, best_ang);
  return best_ang * earth_radius;
}

// 5. Continuous 1D (arc) and 2D (gnomonic box) branch-and-bound solvers.

bool same_undirected_edges_s2(const UnitS2& x, const UnitS2& y,
                              bool check_poly) {
  const int ne = x.ne();
  if (ne != y.ne() || ne == 0) return false;
  if (check_poly && x.has_poly != y.has_poly) return false;
  typedef std::array<double, 7> EdgeKey;
  auto build_keys = [&](const UnitS2& u) {
    std::vector<EdgeKey> keys(ne);
    for (int k = 0; k < ne; k++) {
      const Vec3& A = u.v[u.a[k]];
      const Vec3& B = (u.b[k] >= 0) ? u.v[u.b[k]] : A;
      const bool swap_ab = vec3_less(B, A);
      const Vec3& P = swap_ab ? B : A;
      const Vec3& Q = swap_ab ? A : B;
      keys[k] = {{P.x, P.y, P.z, Q.x, Q.y, Q.z,
                  check_poly ? static_cast<double>(u.e_is_poly[k]) : 0.0}};
    }
    std::sort(keys.begin(), keys.end());
    return keys;
  };
  return build_keys(x) == build_keys(y);
}

// Shared priority-queue pop with user-interrupt check, iteration cap, and
// bound-gap stop condition.
template <typename TaskT>
bool pop_next_task_s2(std::priority_queue<TaskT>& pq, int& iter,
                      double stop_lim, double& unconv_upper, TaskT& out) {
  if (pq.empty()) return false;
  if ((++iter & 1023) == 0) Rcpp::checkUserInterrupt();
  if (iter > kMaxBranchIter) {
    unconv_upper = std::max(unconv_upper, pq.top().upper);
    return false;
  }
  out = pq.top();
  pq.pop();
  return out.upper > stop_lim;
}

struct EdgeScanS2 {
  std::vector<double> d_c2;
  double d_edge_min;
  double poly_ang;
  bool has_poly_edge;
  bool pruned;
};

// Scans active edges of `y` from query point `Q` (with enclosing radius
// `radius`), applying early distance pruning and optional 1D arc-cover pruning.
EdgeScanS2 scan_active_edges_s2(const Vec3& Q, const UnitS2& y,
                                const std::vector<int>& active,
                                double radius, double best_ang,
                                bool allow_early, bool track_poly,
                                const Vec3* arc_A = nullptr,
                                const Vec3* arc_B = nullptr,
                                double cover_lim = -1.0) {
  const int na = static_cast<int>(active.size());
  EdgeScanS2 res;
  res.d_c2.assign(na, kInf);
  res.d_edge_min = kInf;
  res.poly_ang = kInf;
  res.has_poly_edge = false;
  res.pruned = false;

  const double early_c2 =
    (allow_early && best_ang > radius) ? angle_to_chord2(best_ang - radius)
                                       : -1.0;
  double best_c2 = kInf;
  int best_idx = -1;
  int best_kind = 0;
  double best_poly_c2 = kInf;
  int best_poly_idx = -1;
  int best_poly_kind = 0;

  for (int i = 0; i < na; i++) {
    const int ky = active[i];
    int kind;
    const double c2 = edge_chord2(Q, y, ky, kind);
    res.d_c2[i] = c2;
    if (c2 < best_c2) {
      best_c2 = c2;
      best_idx = i;
      best_kind = kind;
    }
    if (track_poly && y.e_is_poly[ky] && c2 < best_poly_c2) {
      best_poly_c2 = c2;
      best_poly_idx = i;
      best_poly_kind = kind;
    }
    if (early_c2 >= 0.0 && c2 <= early_c2) {
      res.pruned = true;
      return res;
    }
    if (cover_lim > 0.0 && y.b[ky] >= 0 &&
        dot3(y.ta[ky], *arc_A) >= -1e-15 &&
        dot3(y.tb[ky], *arc_A) >= -1e-15 &&
        dot3(y.ta[ky], *arc_B) >= -1e-15 &&
        dot3(y.tb[ky], *arc_B) >= -1e-15) {
      const Vec3& N_ky = y.n[ky];
      if (std::max(std::fabs(dot3(N_ky, *arc_A)),
                   std::fabs(dot3(N_ky, *arc_B))) <= cover_lim) {
        res.pruned = true;
        return res;
      }
    }
  }

  if (best_idx < 0) {
    res.pruned = true;
    return res;
  }
  res.d_edge_min = edge_angle(Q, y, active[best_idx], best_kind);
  if (track_poly && best_poly_idx >= 0) {
    res.has_poly_edge = true;
    res.poly_ang = edge_angle(Q, y, active[best_poly_idx], best_poly_kind);
  }
  return res;
}

// Filters active edges of `y` for child tasks via the spherical triangle
// inequality: any edge farther than d_edge_min + 2 * radius from `Q` cannot be
// closest to any point inside the cap of radius `radius` around `Q`.
std::vector<int> filter_active_edges_s2(const std::vector<int>& active,
                                        const std::vector<double>& d_c2,
                                        double d_edge_min, double radius) {
  const double keep_lim =
    d_edge_min + 2.0 * radius + 1e-12 * std::max(1.0, radius);
  const double keep_c2 = angle_to_chord2(keep_lim);
  std::vector<int> next;
  next.reserve(active.size());
  for (std::size_t i = 0; i < active.size(); i++) {
    if (d_c2[i] <= keep_c2) next.push_back(active[i]);
  }
  return next;
}

struct ArcTask {
  Vec3 A, B;
  double da, db;
  double half_ang;
  double upper;
  std::vector<int> active;
  bool operator<(const ArcTask& o) const { return upper < o.upper; }
};

void directed_haus_1d_s2(const UnitS2& x, const UnitS2& y, bool check_pip,
                         const std::vector<double>& ang_v,
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
      pq.push({x.v[ia], x.v[ib], da, db, h, upper, all_y});
    }
  }

  int iter = 0;
  ArcTask cur;
  while (pop_next_task_s2(pq, iter, best_ang + tol_ang, unconv_upper, cur)) {
    const double h = cur.half_ang;
    if (h <= min_ang) continue;

    Vec3 M = {cur.A.x + cur.B.x, cur.A.y + cur.B.y, cur.A.z + cur.B.z};
    if (!normalize3(M)) continue;

    const double cover_lim =
      (best_ang + 0.5 * tol_ang < 0.49 * kPi && h < 0.45 * kPi)
        ? (std::sin(best_ang + 0.5 * tol_ang) * std::cos(h) - 1e-14)
        : -1.0;
    const bool in_y = check_pip && point_in_poly_s2(M, y);
    const EdgeScanS2 scan = scan_active_edges_s2(
      M, y, cur.active, h, best_ang, true, false, &cur.A, &cur.B, cover_lim
    );
    if (scan.pruned) continue;
    if (in_y && scan.d_edge_min >= h * (1.0 + 1e-12)) continue;

    const double dm = in_y ? 0.0 : scan.d_edge_min;
    if (dm > best_ang) best_ang = dm;

    const double child_h = 0.5 * h;
    if (child_h <= min_ang) continue;

    std::vector<int> child_active =
      filter_active_edges_s2(cur.active, scan.d_c2, scan.d_edge_min, h);
    const double u_left = 0.5 * (cur.da + dm + h);
    const double u_right = 0.5 * (dm + cur.db + h);
    const bool push_left = (u_left > best_ang + tol_ang);
    const bool push_right = (u_right > best_ang + tol_ang);
    if (push_left) {
      pq.push({cur.A, M, cur.da, dm, child_h, u_left,
               push_right ? child_active : std::move(child_active)});
    }
    if (push_right) {
      pq.push({M, cur.B, dm, cur.db, child_h, u_right,
               std::move(child_active)});
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

inline Vec3 gnomonic_to_s2(const UnitS2& x, double u, double v) {
  Vec3 out = {
    x.cap_c.x + u * x.e1.x + v * x.e2.x,
    x.cap_c.y + u * x.e1.y + v * x.e2.y,
    x.cap_c.z + u * x.e1.z + v * x.e2.z
  };
  normalize3(out);
  return out;
}

void directed_haus_poly_2d_s2(const UnitS2& x, const UnitS2& y, bool check_pip_y,
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
  pq.push({x.poly_min_u, x.poly_min_v, x.poly_max_u, x.poly_max_v,
           kPi, std::move(x_poly_all), std::move(y_all)});

  int iter = 0;
  GnomonicBox cur;
  while (pop_next_task_s2(pq, iter, best_ang + tol_ang, unconv_upper, cur)) {
    const double um = 0.5 * (cur.u0 + cur.u1);
    const double vm = 0.5 * (cur.v0 + cur.v1);
    const Vec3 C = gnomonic_to_s2(x, um, vm);

    const Vec3 K[4] = {
      gnomonic_to_s2(x, cur.u0, cur.v0),
      gnomonic_to_s2(x, cur.u1, cur.v0),
      gnomonic_to_s2(x, cur.u1, cur.v1),
      gnomonic_to_s2(x, cur.u0, cur.v1)
    };
    double r = 0.0;
    for (int j = 0; j < 4; j++) r = std::max(r, gc_angle(C, K[j]));
    if (r <= min_ang) continue;

    std::vector<int> x_cur;
    x_cur.reserve(cur.x_poly_edges.size());
    for (int k : cur.x_poly_edges) {
      const Pt2& A = x.g[x.a[k]];
      const Pt2& B = x.g[x.b[k]];
      if (seg_intersects_box(A.x, A.y, B.x, B.y,
                             cur.u0, cur.v0, cur.u1, cur.v1)) {
        x_cur.push_back(k);
      }
    }
    const bool in_x = point_in_poly_gnomonic(um, vm, x);
    if (x_cur.empty() && !in_x) continue;

    const bool valid_cap = (r < 0.49 * kPi);
    const bool in_y = check_pip_y && point_in_poly_s2(C, y);
    const EdgeScanS2 scan = scan_active_edges_s2(
      C, y, cur.y_edges, r, best_ang, valid_cap, valid_cap && in_y
    );
    if (scan.pruned) continue;
    if (valid_cap && in_y && scan.has_poly_edge &&
        scan.poly_ang >= r * (1.0 + 1e-12)) {
      continue;
    }

    const double dc = in_y ? 0.0 : scan.d_edge_min;
    if (in_x && dc > best_ang) best_ang = dc;

    const double upper = valid_cap ? std::min(kPi, dc + r) : kPi;
    if (upper <= best_ang + tol_ang) continue;

    std::vector<int> y_next = valid_cap
      ? filter_active_edges_s2(cur.y_edges, scan.d_c2, scan.d_edge_min, r)
      : cur.y_edges;

    GnomonicBox c1 = {cur.u0, cur.v0, cur.u1, cur.v1, upper, x_cur, y_next};
    GnomonicBox c2 = {cur.u0, cur.v0, cur.u1, cur.v1, upper,
                      std::move(x_cur), std::move(y_next)};
    if (cur.u1 - cur.u0 >= cur.v1 - cur.v0) {
      c1.u1 = um;
      c2.u0 = um;
    } else {
      c1.v1 = vm;
      c2.v0 = vm;
    }
    pq.push(std::move(c1));
    pq.push(std::move(c2));
  }
}

double haus_s2_pair(const UnitS2& x_in, const UnitS2& y_in,
                    HausMethod method, double tol, double earth_radius,
                    int& n_unconv, double& max_gap_m) {
  if (method == kDiscrete) {
    return haus_discrete_s2(x_in, y_in, earth_radius);
  }
  const bool is_polygon = (method == kPolygon);
  if (same_undirected_edges_s2(x_in, y_in, is_polygon)) {
    return 0.0;
  }

  const bool swap_order = should_swap_units(
    x_in.nv(), x_in.ne(), std::make_pair(x_in.cap_rad, x_in.v[0].x),
    y_in.nv(), y_in.ne(), std::make_pair(y_in.cap_rad, y_in.v[0].x)
  );
  const UnitS2& a = swap_order ? y_in : x_in;
  const UnitS2& b = swap_order ? x_in : y_in;

  const bool overlap = is_polygon && units_caps_overlap_s2(a, b);
  const bool check_pip_b = overlap && b.has_poly;
  const bool check_pip_a = overlap && a.has_poly;

  const double diam_a = a.has_cap ? 2.0 * a.cap_rad * earth_radius
                                  : kPi * earth_radius;
  const double diam_b = b.has_cap ? 2.0 * b.cap_rad * earth_radius
                                  : kPi * earth_radius;
  const double scale_m = std::max(1.0, std::max(diam_a, diam_b));
  const double tol_ang =
    std::max(tol * scale_m / earth_radius, 64.0 * kEps);
  const double min_ang = 8.0 * kEps;

  double max_c2 = -1.0;
  double best_ang = 0.0;
  double unconv_upper = 0.0;
  std::vector<double> ang_a, ang_b;
  vertex_pass_s2(a, b, check_pip_b, max_c2, best_ang, &ang_a);
  vertex_pass_s2(b, a, check_pip_a, max_c2, best_ang, &ang_b);

  directed_haus_1d_s2(a, b, check_pip_b, ang_a,
                      tol_ang, min_ang, best_ang, unconv_upper);
  directed_haus_1d_s2(b, a, check_pip_a, ang_b,
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

void warn_poly_cap_exceeded(const std::vector<UnitS2>& u, const char* arg) {
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

}  // namespace

//' @title Spherical Hausdorff distance matrix (internal)
//' @param poly_list list of matrices from .to_list_mat() or .to_polygon_list().
//' @param method one of "discrete", "boundary", or "polygon".
//' @param tol convergence tolerance for continuous methods.
//' @param earth_radius radius of the sphere, in meters.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd dist_haus_s2(const Rcpp::List& poly_list,
                             const std::string& method = "discrete",
                             double tol = 1e-8,
                             double earth_radius = 6371010.0) {
  const HausMethod code = parse_haus_method(method, tol);
  const std::vector<UnitS2> u = make_unit_vec<UnitS2>(
    poly_list, [&](const Eigen::Ref<const Eigen::MatrixXd>& m) {
      return make_unit_s2(m, code);
    }
  );
  if (code == kPolygon) warn_poly_cap_exceeded(u, "x");

  int n_unconv = 0;
  double max_gap_m = 0.0;
  Eigen::MatrixXd out = sym_dist_matrix(
    static_cast<int>(u.size()),
    [&](int i) { return u[i].valid_for(code); },
    [&](int i, int j) {
      return haus_s2_pair(u[i], u[j], code, tol, earth_radius,
                          n_unconv, max_gap_m);
    }
  );
  warn_unconverged(n_unconv, max_gap_m);
  return out;
}

//' @title Spherical Hausdorff cross-distance matrix (internal)
//' @param poly_list,poly_pred lists of matrices from .to_list_mat() or
//'   .to_polygon_list().
//' @param method one of "discrete", "boundary", or "polygon".
//' @param tol convergence tolerance for continuous methods.
//' @param earth_radius radius of the sphere, in meters.
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cdist_haus_s2(const Rcpp::List& poly_list,
                              const Rcpp::List& poly_pred,
                              const std::string& method = "discrete",
                              double tol = 1e-8,
                              double earth_radius = 6371010.0) {
  const HausMethod code = parse_haus_method(method, tol);
  const auto builder = [&](const Eigen::Ref<const Eigen::MatrixXd>& m) {
    return make_unit_s2(m, code);
  };
  const std::vector<UnitS2> u1 = make_unit_vec<UnitS2>(poly_list, builder);
  const std::vector<UnitS2> u2 = make_unit_vec<UnitS2>(poly_pred, builder);
  if (code == kPolygon) {
    warn_poly_cap_exceeded(u1, "x");
    warn_poly_cap_exceeded(u2, "y");
  }

  int n_unconv = 0;
  double max_gap_m = 0.0;
  Eigen::MatrixXd out = cross_dist_matrix(
    static_cast<int>(u1.size()),
    static_cast<int>(u2.size()),
    [&](int i) { return u1[i].valid_for(code); },
    [&](int j) { return u2[j].valid_for(code); },
    [&](int i, int j) {
      return haus_s2_pair(u1[i], u2[j], code, tol, earth_radius,
                          n_unconv, max_gap_m);
    }
  );
  warn_unconverged(n_unconv, max_gap_m);
  return out;
}
