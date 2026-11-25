#include "geom-common.h"
#include <cstdint>
#include <string>

// [[Rcpp::depends(RcppEigen)]]

// Center-radii (CR) distance between spatial units.
//
// The CR distance between units i and j is d(c_i, c_j) + |r_i - r_j|, where
// c and r are the center and radius of the minimum enclosing ball (the
// Chebyshev center and radius) of the vertices of a unit. Welzl's algorithm
// computes the exact ball. In the plane, d is the Euclidean distance in the
// units of the coordinates. On the sphere, d is the great-circle distance in
// meters on a sphere of radius `earth_radius`.
//
// Input: the coordinates of all units (x and y, or longitude and latitude in
// degrees), with `feature_id` (1-based) identifying the unit of each row, as
// returned by wk::wk_coords().

namespace {

using smile_geom::Pt2;
using smile_geom::Vec3;
using smile_geom::cross3;
using smile_geom::dist2_2d;
using smile_geom::dist2_3d;
using smile_geom::dot3;
using smile_geom::gc_angle;
using smile_geom::kGoldenRatio;
using smile_geom::kPi;
using smile_geom::ll_to_xyz;
using smile_geom::norm3;
using smile_geom::tangent_basis;
using smile_geom::xyz_to_ll;

// 1. Generic Welzl minimum enclosing ball solver and feature grouping.

template <typename PtT, typename BallT,
          typename Ball1Fn, typename Ball2Fn, typename Ball3Fn,
          typename OutsideFn>
BallT welzl_ball(const std::vector<PtT>& pts,
                 Ball1Fn ball1, Ball2Fn ball2, Ball3Fn ball3,
                 OutsideFn is_outside) {
  const int n = static_cast<int>(pts.size());
  auto with_2 = [&](int m, const PtT& p, const PtT& q) -> BallT {
    BallT b = ball2(p, q);
    for (int k = 0; k < m; k++) {
      if (is_outside(pts[k], b)) b = ball3(p, q, pts[k]);
    }
    return b;
  };
  auto with_1 = [&](int m, const PtT& p) -> BallT {
    BallT b = ball1(p);
    for (int j = 0; j < m; j++) {
      if (is_outside(pts[j], b)) b = with_2(j, p, pts[j]);
    }
    return b;
  };
  BallT b = ball1(pts[0]);
  for (int i = 1; i < n; i++) {
    if (is_outside(pts[i], b)) b = with_1(i, pts[i]);
  }
  return b;
}

inline int gcd_int(int a, int b) {
  while (b != 0) {
    const int t = a % b;
    a = b;
    b = t;
  }
  return a;
}

// Step coprime to `n` near `n * kGoldenRatio`. Stepping through vertices by
// this stride visits spread boundary points early in O(n) time without sorting.
inline int golden_step(int n) {
  if (n <= 1) return 1;
  int step = static_cast<int>(n * kGoldenRatio);
  if (step < 1) step = 1;
  if (step >= n) step = n - 1;
  while (gcd_int(step, n) != 1) {
    step++;
    if (step >= n) step = 1;
  }
  return step;
}

template <typename PtT>
struct GroupedPoints {
  std::vector<PtT> pts;
  std::vector<R_xlen_t> offsets;
};

// Groups coordinate vectors into one flat buffer indexed by `offsets` for
// 1-based `feature_id`. Any feature with zero rows or at least one non-finite
// coordinate gets an empty slice (`offsets[i] == offsets[i + 1]`).
template <typename PtT, typename MakePtFn>
GroupedPoints<PtT> group_feature_points(
    const Rcpp::NumericVector& c1,
    const Rcpp::NumericVector& c2,
    const Rcpp::IntegerVector& feature_id,
    int n_features,
    MakePtFn make_pt) {
  GroupedPoints<PtT> out;
  out.offsets.assign(n_features + 1, 0);
  if (n_features == 0) return out;

  const R_xlen_t m = feature_id.size();
  if (c1.size() != m || c2.size() != m) {
    Rcpp::stop("Coordinate vectors and feature_id must have the same length.");
  }
  const double* p1 = REAL(c1);
  const double* p2 = REAL(c2);
  const int* pfid = INTEGER(feature_id);

  std::vector<R_xlen_t> counts(n_features, 0);
  std::vector<uint8_t> has_na(n_features, 0);
  for (R_xlen_t k = 0; k < m; k++) {
    const int fid = pfid[k];
    if (fid == NA_INTEGER || fid < 1 || fid > n_features) {
      Rcpp::stop("feature_id contains out-of-range indices.");
    }
    const int idx = fid - 1;
    if (!std::isfinite(p1[k]) || !std::isfinite(p2[k])) {
      has_na[idx] = 1;
    } else {
      counts[idx]++;
    }
  }

  for (int i = 0; i < n_features; i++) {
    const R_xlen_t valid_count = has_na[i] ? 0 : counts[i];
    out.offsets[i + 1] = out.offsets[i] + valid_count;
  }
  out.pts.resize(out.offsets[n_features]);
  std::vector<R_xlen_t> cursor = out.offsets;
  for (R_xlen_t k = 0; k < m; k++) {
    const int idx = pfid[k] - 1;
    if (!has_na[idx]) out.pts[cursor[idx]++] = make_pt(p1[k], p2[k]);
  }
  return out;
}

// 2. Planar (R^2) minimum enclosing disk.

struct Disk2 {
  double cx;
  double cy;
  double r2;
  double lim2;
};

inline Disk2 make_disk(double cx, double cy, double r2, double tol) {
  const double r = std::sqrt(r2);
  const double lim = r + tol * (1.0 + r);
  return {cx, cy, r2, lim * lim};
}

inline Disk2 disk2(const Pt2& a, const Pt2& b, double tol) {
  const double cx = 0.5 * (a.x + b.x);
  const double cy = 0.5 * (a.y + b.y);
  return make_disk(cx, cy, dist2_2d(a.x, a.y, cx, cy), tol);
}

inline Disk2 max_diametral_disk(const Pt2& a, const Pt2& b, const Pt2& c,
                                double tol) {
  Disk2 d = disk2(a, b, tol);
  const Disk2 d_bc = disk2(b, c, tol);
  if (d_bc.r2 > d.r2) d = d_bc;
  const Disk2 d_ac = disk2(a, c, tol);
  if (d_ac.r2 > d.r2) d = d_ac;
  return d;
}

inline Disk2 disk3(const Pt2& a, const Pt2& b, const Pt2& c, double tol) {
  const double ux = b.x - a.x;
  const double uy = b.y - a.y;
  const double vx = c.x - a.x;
  const double vy = c.y - a.y;
  const double u2 = ux * ux + uy * uy;
  const double v2 = vx * vx + vy * vy;
  const double det = 2.0 * (ux * vy - uy * vx);
  const double scale2 = u2 + v2;
  if (scale2 == 0.0 || std::fabs(det) <= 1e-14 * scale2) {
    return max_diametral_disk(a, b, c, tol);
  }
  const double dx = (vy * u2 - uy * v2) / det;
  const double dy = (ux * v2 - vx * u2) / det;
  return make_disk(a.x + dx, a.y + dy, dx * dx + dy * dy, tol);
}

inline bool outside_disk(const Pt2& p, const Disk2& d) {
  return dist2_2d(p.x, p.y, d.cx, d.cy) > d.lim2;
}

Disk2 welzl_disk(const std::vector<Pt2>& pts, double tol = 1e-12) {
  const double tol2 = tol * tol;
  return welzl_ball<Pt2, Disk2>(
    pts,
    [tol2](const Pt2& p) -> Disk2 { return {p.x, p.y, 0.0, tol2}; },
    [tol](const Pt2& a, const Pt2& b) { return disk2(a, b, tol); },
    [tol](const Pt2& a, const Pt2& b, const Pt2& c) {
      return disk3(a, b, c, tol);
    },
    outside_disk
  );
}

struct Meb2Result {
  Eigen::ArrayXd cx;
  Eigen::ArrayXd cy;
  Eigen::ArrayXd radius;
  std::vector<int> invalid;
};

Meb2Result compute_meb_planar(const Rcpp::NumericVector& x,
                              const Rcpp::NumericVector& y,
                              const Rcpp::IntegerVector& feature_id,
                              int n_features) {
  Meb2Result res;
  res.cx = Eigen::ArrayXd::Zero(n_features);
  res.cy = Eigen::ArrayXd::Zero(n_features);
  res.radius = Eigen::ArrayXd::Zero(n_features);

  const GroupedPoints<Pt2> grouped = group_feature_points<Pt2>(
    x, y, feature_id, n_features,
    [](double px, double py) -> Pt2 { return {px, py}; }
  );

  std::vector<Pt2> shifted;
  for (int i = 0; i < n_features; i++) {
    const R_xlen_t beg = grouped.offsets[i];
    const int cnt = static_cast<int>(grouped.offsets[i + 1] - beg);
    if (cnt == 0) {
      res.invalid.push_back(i);
      continue;
    }
    if (cnt == 1) {
      res.cx[i] = grouped.pts[beg].x;
      res.cy[i] = grouped.pts[beg].y;
      res.radius[i] = 0.0;
      continue;
    }

    // Shift by the first vertex to preserve precision on large projected
    // coordinates such as UTM or Web Mercator, while scrambling in one pass.
    const double ox = grouped.pts[beg].x;
    const double oy = grouped.pts[beg].y;
    shifted.resize(cnt);
    const int step = golden_step(cnt);
    int idx = 0;
    for (int k = 0; k < cnt; k++) {
      const Pt2& p = grouped.pts[beg + idx];
      shifted[k] = {p.x - ox, p.y - oy};
      idx += step;
      if (idx >= cnt) idx -= cnt;
    }

    const Disk2 d = welzl_disk(shifted);
    const Pt2 center_s = {d.cx, d.cy};
    double max_r2 = 0.0;
    for (const Pt2& p : shifted) {
      max_r2 = std::max(max_r2, dist2_2d(p, center_s));
    }
    res.cx[i] = d.cx + ox;
    res.cy[i] = d.cy + oy;
    res.radius[i] = std::sqrt(max_r2);
  }
  return res;
}

Rcpp::NumericMatrix planar_sym_matrix(const Meb2Result& m) {
  const int n = static_cast<int>(m.cx.size());
  if (n == 0) return Rcpp::NumericMatrix(0, 0);
  Rcpp::NumericMatrix rmat(Rcpp::no_init_matrix(n, n));
  Eigen::Map<Eigen::MatrixXd> out(REAL(rmat), n, n);

  for (int j = 0; j < n; j++) {
    if ((j & 255) == 0) Rcpp::checkUserInterrupt();
    out(j, j) = 0.0;
    const int rem = n - j - 1;
    if (rem > 0) {
      out.col(j).tail(rem) =
        ((m.cx.tail(rem) - m.cx[j]).square() +
         (m.cy.tail(rem) - m.cy[j]).square()).sqrt() +
        (m.radius.tail(rem) - m.radius[j]).abs();
    }
  }
  out.triangularView<Eigen::StrictlyUpper>() =
    out.triangularView<Eigen::StrictlyLower>().transpose();
  for (int idx : m.invalid) {
    out.col(idx).setConstant(NA_REAL);
    out.row(idx).setConstant(NA_REAL);
  }
  return rmat;
}

Rcpp::NumericMatrix planar_cross_matrix(const Meb2Result& r1,
                                        const Meb2Result& r2) {
  const int n1 = static_cast<int>(r1.cx.size());
  const int n2 = static_cast<int>(r2.cx.size());
  if (n1 == 0 || n2 == 0) return Rcpp::NumericMatrix(n1, n2);
  Rcpp::NumericMatrix rmat(Rcpp::no_init_matrix(n1, n2));
  Eigen::Map<Eigen::MatrixXd> out(REAL(rmat), n1, n2);

  for (int j = 0; j < n2; j++) {
    if ((j & 255) == 0) Rcpp::checkUserInterrupt();
    out.col(j) =
      ((r1.cx - r2.cx[j]).square() +
       (r1.cy - r2.cy[j]).square()).sqrt() +
      (r1.radius - r2.radius[j]).abs();
  }
  for (int j : r2.invalid) out.col(j).setConstant(NA_REAL);
  for (int i : r1.invalid) out.row(i).setConstant(NA_REAL);
  return rmat;
}

// 3. Spherical (S^2) minimum enclosing cap.

struct Cap3 {
  Vec3 c;
  double r2;    // Squared chord radius
  double lim2;  // Squared chord tolerance bound
  bool ok;
};

inline Cap3 invalid_cap() {
  return {{0.0, 0.0, 0.0}, 4.0, 4.0, false};
}

inline Cap3 make_cap(const Vec3& c, double r2, double tol) {
  const double lim = std::sqrt(r2) + tol;
  return {c, r2, lim * lim, true};
}

inline Cap3 cap2(const Vec3& a, const Vec3& b, double tol) {
  const Vec3 s = {a.x + b.x, a.y + b.y, a.z + b.z};
  const double n = norm3(s);
  if (n < 1e-15) return invalid_cap();
  const Vec3 u = {s.x / n, s.y / n, s.z / n};
  return make_cap(u, std::max(dist2_3d(a, u), dist2_3d(b, u)), tol);
}

inline Cap3 max_chord_cap2(const Vec3& a, const Vec3& b, const Vec3& c,
                           double tol) {
  Cap3 out = cap2(a, b, tol);
  const Cap3 cb = cap2(b, c, tol);
  if (!out.ok || (cb.ok && cb.r2 > out.r2)) out = cb;
  const Cap3 ca = cap2(a, c, tol);
  if (!out.ok || (ca.ok && ca.r2 > out.r2)) out = ca;
  return out;
}

inline Cap3 cap3(const Vec3& a, const Vec3& b, const Vec3& c, double tol) {
  const Vec3 sum = {a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z};
  const double m_norm = norm3(sum);
  if (m_norm < 1e-15) return max_chord_cap2(a, b, c, tol);
  const Vec3 m = {sum.x / m_norm, sum.y / m_norm, sum.z / m_norm};

  Vec3 e1, e2;
  if (!tangent_basis(m, e1, e2)) return max_chord_cap2(a, b, c, tol);

  const Vec3 ba = {b.x - a.x, b.y - a.y, b.z - a.z};
  const Vec3 ca = {c.x - a.x, c.y - a.y, c.z - a.z};
  const double a11 = dot3(ba, e1);
  const double a12 = dot3(ba, e2);
  const double a21 = dot3(ca, e1);
  const double a22 = dot3(ca, e2);

  const Vec3 am = {a.x - m.x, a.y - m.y, a.z - m.z};
  const Vec3 bm = {b.x - m.x, b.y - m.y, b.z - m.z};
  const Vec3 cm = {c.x - m.x, c.y - m.y, c.z - m.z};
  const double da = dot3(am, am);
  const double r1 = 0.5 * (dot3(bm, bm) - da);
  const double r2 = 0.5 * (dot3(cm, cm) - da);
  const double det = a11 * a22 - a12 * a21;
  const double scale = dot3(ba, ba) + dot3(ca, ca);

  Vec3 u;
  if (scale > 0.0 && std::fabs(det) > 1e-14 * scale) {
    const double c1 = (r1 * a22 - a12 * r2) / det;
    const double c2 = (a11 * r2 - r1 * a21) / det;
    u = {m.x + c1 * e1.x + c2 * e2.x,
         m.y + c1 * e1.y + c2 * e2.y,
         m.z + c1 * e1.z + c2 * e2.z};
  } else {
    u = cross3(ba, ca);
    if (dot3(u, m) < 0.0) u = {-u.x, -u.y, -u.z};
  }
  const double u_norm = norm3(u);
  if (u_norm < 1e-15) return max_chord_cap2(a, b, c, tol);
  u = {u.x / u_norm, u.y / u_norm, u.z / u_norm};
  const double max_d2 =
    std::max(dist2_3d(a, u), std::max(dist2_3d(b, u), dist2_3d(c, u)));
  return make_cap(u, max_d2, tol);
}

inline bool outside_cap(const Vec3& p, const Cap3& cap) {
  return cap.ok && dist2_3d(p, cap.c) > cap.lim2;
}

Cap3 welzl_cap(const std::vector<Vec3>& pts, double tol = 1e-13) {
  const double tol2 = tol * tol;
  return welzl_ball<Vec3, Cap3>(
    pts,
    [tol2](const Vec3& p) -> Cap3 { return {p, 0.0, tol2, true}; },
    [tol](const Vec3& a, const Vec3& b) { return cap2(a, b, tol); },
    [tol](const Vec3& a, const Vec3& b, const Vec3& c) {
      return cap3(a, b, c, tol);
    },
    outside_cap
  );
}

struct Meb3Result {
  std::vector<double> lon;
  std::vector<double> lat;
  Eigen::ArrayXd ux;
  Eigen::ArrayXd uy;
  Eigen::ArrayXd uz;
  Eigen::ArrayXd radius;
  std::vector<int> invalid;
  std::vector<int> non_hemi;
};

Meb3Result compute_meb_sphere(const Rcpp::NumericVector& lon,
                              const Rcpp::NumericVector& lat,
                              const Rcpp::IntegerVector& feature_id,
                              int n_features,
                              double earth_radius,
                              bool include_ll) {
  Meb3Result res;
  if (include_ll) {
    res.lon.assign(n_features, NA_REAL);
    res.lat.assign(n_features, NA_REAL);
  }
  res.ux = Eigen::ArrayXd::Ones(n_features);
  res.uy = Eigen::ArrayXd::Zero(n_features);
  res.uz = Eigen::ArrayXd::Zero(n_features);
  res.radius = Eigen::ArrayXd::Zero(n_features);

  const GroupedPoints<Pt2> grouped = group_feature_points<Pt2>(
    lon, lat, feature_id, n_features,
    [](double plon, double plat) -> Pt2 { return {plon, plat}; }
  );

  std::vector<Vec3> pts;
  for (int i = 0; i < n_features; i++) {
    const R_xlen_t beg = grouped.offsets[i];
    const int cnt = static_cast<int>(grouped.offsets[i + 1] - beg);
    if (cnt == 0) {
      res.invalid.push_back(i);
      continue;
    }
    if (cnt == 1) {
      const Pt2& p = grouped.pts[beg];
      const Vec3 u = ll_to_xyz(p.x, p.y);
      if (include_ll) {
        if (p.x >= -180.0 && p.x <= 180.0 && p.y > -90.0 && p.y < 90.0) {
          res.lon[i] = p.x;
          res.lat[i] = p.y;
        } else {
          xyz_to_ll(u, res.lon[i], res.lat[i]);
        }
      }
      res.ux[i] = u.x;
      res.uy[i] = u.y;
      res.uz[i] = u.z;
      res.radius[i] = 0.0;
      continue;
    }

    pts.resize(cnt);
    const int step = golden_step(cnt);
    int idx = 0;
    for (int k = 0; k < cnt; k++) {
      const Pt2& p = grouped.pts[beg + idx];
      pts[k] = ll_to_xyz(p.x, p.y);
      idx += step;
      if (idx >= cnt) idx -= cnt;
    }

    const Cap3 cap = welzl_cap(pts);
    if (!cap.ok) {
      res.invalid.push_back(i);
      res.non_hemi.push_back(i + 1);
      continue;
    }

    // Chord distance increases strictly with central angle on S^2, so find the
    // farthest vertex with squared distances and compute one angle at the end.
    double max_d2 = -1.0;
    int far_k = 0;
    for (int k = 0; k < cnt; k++) {
      const double d2 = dist2_3d(pts[k], cap.c);
      if (d2 > max_d2) {
        max_d2 = d2;
        far_k = k;
      }
    }
    const double max_angle = gc_angle(pts[far_k], cap.c);
    if (!std::isfinite(max_angle) || max_angle >= 0.5 * kPi) {
      res.invalid.push_back(i);
      res.non_hemi.push_back(i + 1);
      continue;
    }

    if (include_ll) xyz_to_ll(cap.c, res.lon[i], res.lat[i]);
    res.ux[i] = cap.c.x;
    res.uy[i] = cap.c.y;
    res.uz[i] = cap.c.z;
    res.radius[i] = max_angle * earth_radius;
  }
  return res;
}

void warn_non_hemisphere(const Meb3Result& m, const char* arg) {
  if (m.non_hemi.empty()) return;
  std::string idx_list;
  for (std::size_t k = 0; k < m.non_hemi.size(); k++) {
    if (k > 0) idx_list += ", ";
    idx_list += std::to_string(m.non_hemi[k]);
  }
  Rcpp::warning(
    "These units of %s are not inside an open hemisphere: %s. Their "
    "Chebyshev radius is at least pi * R / 2. The algorithm needs a "
    "radius below pi * R / 2, so the function returns NA for them.",
    arg, idx_list.c_str()
  );
}

Rcpp::NumericMatrix sphere_sym_matrix(const Meb3Result& m,
                                      double earth_radius) {
  const int n = static_cast<int>(m.ux.size());
  if (n == 0) return Rcpp::NumericMatrix(0, 0);
  Rcpp::NumericMatrix rmat(Rcpp::no_init_matrix(n, n));
  Eigen::Map<Eigen::MatrixXd> out(REAL(rmat), n, n);

  Eigen::ArrayXd dm(n);
  Eigen::ArrayXd dp(n);
  const double scale = 2.0 * earth_radius;
  auto atan2_op = [](double y, double x) { return std::atan2(y, x); };

  for (int j = 0; j < n; j++) {
    if ((j & 255) == 0) Rcpp::checkUserInterrupt();
    out(j, j) = 0.0;
    const int rem = n - j - 1;
    if (rem > 0) {
      dm.head(rem) =
        ((m.ux.tail(rem) - m.ux[j]).square() +
         (m.uy.tail(rem) - m.uy[j]).square() +
         (m.uz.tail(rem) - m.uz[j]).square()).sqrt();
      dp.head(rem) =
        ((m.ux.tail(rem) + m.ux[j]).square() +
         (m.uy.tail(rem) + m.uy[j]).square() +
         (m.uz.tail(rem) + m.uz[j]).square()).sqrt();
      out.col(j).tail(rem) =
        scale * dm.head(rem).binaryExpr(dp.head(rem), atan2_op) +
        (m.radius.tail(rem) - m.radius[j]).abs();
    }
  }
  out.triangularView<Eigen::StrictlyUpper>() =
    out.triangularView<Eigen::StrictlyLower>().transpose();
  for (int idx : m.invalid) {
    out.col(idx).setConstant(NA_REAL);
    out.row(idx).setConstant(NA_REAL);
  }
  return rmat;
}

Rcpp::NumericMatrix sphere_cross_matrix(const Meb3Result& r1,
                                        const Meb3Result& r2,
                                        double earth_radius) {
  const int n1 = static_cast<int>(r1.ux.size());
  const int n2 = static_cast<int>(r2.ux.size());
  if (n1 == 0 || n2 == 0) return Rcpp::NumericMatrix(n1, n2);
  Rcpp::NumericMatrix rmat(Rcpp::no_init_matrix(n1, n2));
  Eigen::Map<Eigen::MatrixXd> out(REAL(rmat), n1, n2);

  Eigen::ArrayXd dm(n1);
  Eigen::ArrayXd dp(n1);
  const double scale = 2.0 * earth_radius;
  auto atan2_op = [](double y, double x) { return std::atan2(y, x); };

  for (int j = 0; j < n2; j++) {
    if ((j & 255) == 0) Rcpp::checkUserInterrupt();
    dm =
      ((r1.ux - r2.ux[j]).square() +
       (r1.uy - r2.uy[j]).square() +
       (r1.uz - r2.uz[j]).square()).sqrt();
    dp =
      ((r1.ux + r2.ux[j]).square() +
       (r1.uy + r2.uy[j]).square() +
       (r1.uz + r2.uz[j]).square()).sqrt();
    out.col(j) =
      scale * dm.binaryExpr(dp, atan2_op) +
      (r1.radius - r2.radius[j]).abs();
  }
  for (int j : r2.invalid) out.col(j).setConstant(NA_REAL);
  for (int i : r1.invalid) out.row(i).setConstant(NA_REAL);
  return rmat;
}

}  // namespace

//' @title Planar Chebyshev centers and radii (internal)
//' @noRd
// [[Rcpp::export]]
Rcpp::List mbc_planar_cpp(const Rcpp::NumericVector& x,
                          const Rcpp::NumericVector& y,
                          const Rcpp::IntegerVector& feature_id,
                          int n_features) {
  const Meb2Result res = compute_meb_planar(x, y, feature_id, n_features);
  Rcpp::NumericVector cx(n_features);
  Rcpp::NumericVector cy(n_features);
  Rcpp::NumericVector radius(n_features);
  for (int i = 0; i < n_features; i++) {
    cx[i] = res.cx[i];
    cy[i] = res.cy[i];
    radius[i] = res.radius[i];
  }
  for (int idx : res.invalid) {
    cx[idx] = NA_REAL;
    cy[idx] = NA_REAL;
    radius[idx] = NA_REAL;
  }
  return Rcpp::List::create(
    Rcpp::Named("cx") = cx,
    Rcpp::Named("cy") = cy,
    Rcpp::Named("radius") = radius
  );
}

//' @title Spherical Chebyshev centers and radii (internal)
//' @noRd
// [[Rcpp::export]]
Rcpp::List mbc_sphere_cpp(const Rcpp::NumericVector& lon,
                          const Rcpp::NumericVector& lat,
                          const Rcpp::IntegerVector& feature_id,
                          int n_features,
                          double earth_radius = 6371010.0) {
  const Meb3Result res =
    compute_meb_sphere(lon, lat, feature_id, n_features, earth_radius, true);
  warn_non_hemisphere(res, "x");
  Rcpp::NumericVector radius(n_features);
  for (int i = 0; i < n_features; i++) radius[i] = res.radius[i];
  for (int idx : res.invalid) radius[idx] = NA_REAL;
  return Rcpp::List::create(
    Rcpp::Named("lon") = Rcpp::wrap(res.lon),
    Rcpp::Named("lat") = Rcpp::wrap(res.lat),
    Rcpp::Named("radius") = radius
  );
}

//' @title Planar CR distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix cr_planar_cpp(const Rcpp::NumericVector& x,
                                  const Rcpp::NumericVector& y,
                                  const Rcpp::IntegerVector& feature_id,
                                  int n_features) {
  const Meb2Result res = compute_meb_planar(x, y, feature_id, n_features);
  return planar_sym_matrix(res);
}

//' @title Planar CR cross-distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix cr_cross_planar_cpp(const Rcpp::NumericVector& x1,
                                        const Rcpp::NumericVector& y1,
                                        const Rcpp::IntegerVector& fid1,
                                        int n1,
                                        const Rcpp::NumericVector& x2,
                                        const Rcpp::NumericVector& y2,
                                        const Rcpp::IntegerVector& fid2,
                                        int n2) {
  const Meb2Result r1 = compute_meb_planar(x1, y1, fid1, n1);
  const Meb2Result r2 = compute_meb_planar(x2, y2, fid2, n2);
  return planar_cross_matrix(r1, r2);
}

//' @title Spherical CR distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix cr_sphere_cpp(const Rcpp::NumericVector& lon,
                                  const Rcpp::NumericVector& lat,
                                  const Rcpp::IntegerVector& feature_id,
                                  int n_features,
                                  double earth_radius = 6371010.0) {
  const Meb3Result res =
    compute_meb_sphere(lon, lat, feature_id, n_features, earth_radius, false);
  warn_non_hemisphere(res, "x");
  return sphere_sym_matrix(res, earth_radius);
}

//' @title Spherical CR cross-distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Rcpp::NumericMatrix cr_cross_sphere_cpp(const Rcpp::NumericVector& lon1,
                                        const Rcpp::NumericVector& lat1,
                                        const Rcpp::IntegerVector& fid1,
                                        int n1,
                                        const Rcpp::NumericVector& lon2,
                                        const Rcpp::NumericVector& lat2,
                                        const Rcpp::IntegerVector& fid2,
                                        int n2,
                                        double earth_radius = 6371010.0) {
  const Meb3Result r1 =
    compute_meb_sphere(lon1, lat1, fid1, n1, earth_radius, false);
  const Meb3Result r2 =
    compute_meb_sphere(lon2, lat2, fid2, n2, earth_radius, false);
  warn_non_hemisphere(r1, "x");
  warn_non_hemisphere(r2, "y");
  return sphere_cross_matrix(r1, r2, earth_radius);
}
