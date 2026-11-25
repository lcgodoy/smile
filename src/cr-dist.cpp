#include "geom-common.h"
#include <cstdint>

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
using smile_geom::chord_dist;
using smile_geom::cross3;
using smile_geom::cross_dist_matrix;
using smile_geom::dist2_2d;
using smile_geom::dot3;
using smile_geom::gc_angle;
using smile_geom::golden_scramble;
using smile_geom::is_invalid_scalar;
using smile_geom::kPi;
using smile_geom::ll_to_xyz;
using smile_geom::norm3;
using smile_geom::sym_dist_matrix;
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

// Groups coordinate vectors by 1-based `feature_id`. Any feature with zero rows
// or at least one non-finite coordinate produces an empty vector.
template <typename PtT, typename MakePtFn>
std::vector<std::vector<PtT> > group_feature_points(
    const Rcpp::NumericVector& c1,
    const Rcpp::NumericVector& c2,
    const Rcpp::IntegerVector& feature_id,
    int n_features,
    MakePtFn make_pt) {
  std::vector<std::vector<PtT> > grouped(n_features);
  if (n_features == 0) return grouped;

  const R_xlen_t m = feature_id.size();
  std::vector<int> counts(n_features, 0);
  std::vector<uint8_t> has_na(n_features, 0);
  for (R_xlen_t k = 0; k < m; k++) {
    const int fid = feature_id[k];
    if (fid == NA_INTEGER || fid < 1 || fid > n_features) {
      Rcpp::stop("feature_id contains out-of-range indices.");
    }
    const int idx = fid - 1;
    if (is_invalid_scalar(c1[k]) || is_invalid_scalar(c2[k])) {
      has_na[idx] = 1;
    } else {
      counts[idx]++;
    }
  }

  for (int i = 0; i < n_features; i++) {
    if (!has_na[i] && counts[i] > 0) grouped[i].reserve(counts[i]);
  }
  for (R_xlen_t k = 0; k < m; k++) {
    const int idx = feature_id[k] - 1;
    if (!has_na[idx]) grouped[idx].push_back(make_pt(c1[k], c2[k]));
  }
  return grouped;
}

// 2. Planar (R^2) minimum enclosing disk.

struct Disk2 {
  double cx;
  double cy;
  double r2;
};

inline Disk2 disk2(const Pt2& a, const Pt2& b) {
  const double cx = 0.5 * (a.x + b.x);
  const double cy = 0.5 * (a.y + b.y);
  return {cx, cy, dist2_2d(a.x, a.y, cx, cy)};
}

inline Disk2 max_diametral_disk(const Pt2& a, const Pt2& b, const Pt2& c) {
  Disk2 d = disk2(a, b);
  const Disk2 d_bc = disk2(b, c);
  if (d_bc.r2 > d.r2) d = d_bc;
  const Disk2 d_ac = disk2(a, c);
  if (d_ac.r2 > d.r2) d = d_ac;
  return d;
}

inline Disk2 disk3(const Pt2& a, const Pt2& b, const Pt2& c) {
  const double ux = b.x - a.x;
  const double uy = b.y - a.y;
  const double vx = c.x - a.x;
  const double vy = c.y - a.y;
  const double u2 = ux * ux + uy * uy;
  const double v2 = vx * vx + vy * vy;
  const double det = 2.0 * (ux * vy - uy * vx);
  const double scale2 = u2 + v2;
  if (scale2 == 0.0 || std::fabs(det) <= 1e-14 * scale2) {
    return max_diametral_disk(a, b, c);
  }
  const double dx = (vy * u2 - uy * v2) / det;
  const double dy = (ux * v2 - vx * u2) / det;
  return {a.x + dx, a.y + dy, dx * dx + dy * dy};
}

inline bool outside_disk(const Pt2& p, const Disk2& d, double tol) {
  const double r = std::sqrt(d.r2);
  const double lim = r + tol * (1.0 + r);
  return dist2_2d(p.x, p.y, d.cx, d.cy) > (lim * lim);
}

Disk2 welzl_disk(const std::vector<Pt2>& pts, double tol = 1e-12) {
  return welzl_ball<Pt2, Disk2>(
    pts,
    [](const Pt2& p) -> Disk2 { return {p.x, p.y, 0.0}; },
    disk2,
    disk3,
    [tol](const Pt2& p, const Disk2& d) { return outside_disk(p, d, tol); }
  );
}

struct Meb2Result {
  std::vector<double> cx;
  std::vector<double> cy;
  std::vector<double> radius;
  std::vector<uint8_t> valid;
};

Meb2Result compute_meb_planar(const Rcpp::NumericVector& x,
                              const Rcpp::NumericVector& y,
                              const Rcpp::IntegerVector& feature_id,
                              int n_features) {
  Meb2Result res;
  res.cx.assign(n_features, NA_REAL);
  res.cy.assign(n_features, NA_REAL);
  res.radius.assign(n_features, NA_REAL);
  res.valid.assign(n_features, 0);

  const std::vector<std::vector<Pt2> > grouped = group_feature_points<Pt2>(
    x, y, feature_id, n_features,
    [](double px, double py) -> Pt2 { return {px, py}; }
  );

  for (int i = 0; i < n_features; i++) {
    const std::vector<Pt2>& raw = grouped[i];
    if (raw.empty()) continue;
    if (raw.size() == 1) {
      res.cx[i] = raw[0].x;
      res.cy[i] = raw[0].y;
      res.radius[i] = 0.0;
      res.valid[i] = 1;
      continue;
    }

    // Shift by the first vertex to preserve precision on large projected
    // coordinates such as UTM or Web Mercator.
    const double ox = raw[0].x;
    const double oy = raw[0].y;
    std::vector<Pt2> shifted;
    shifted.reserve(raw.size());
    for (const Pt2& p : raw) shifted.push_back({p.x - ox, p.y - oy});
    golden_scramble(shifted);

    const Disk2 d = welzl_disk(shifted);
    const Pt2 center_s = {d.cx, d.cy};
    double max_r2 = 0.0;
    for (const Pt2& p : shifted) {
      max_r2 = std::max(max_r2, dist2_2d(p, center_s));
    }
    res.cx[i] = d.cx + ox;
    res.cy[i] = d.cy + oy;
    res.radius[i] = std::sqrt(max_r2);
    res.valid[i] = 1;
  }
  return res;
}

inline double planar_cr_pair(const Meb2Result& a, int i,
                             const Meb2Result& b, int j) {
  return std::sqrt(dist2_2d(a.cx[i], a.cy[i], b.cx[j], b.cy[j])) +
         std::fabs(a.radius[i] - b.radius[j]);
}

// 3. Spherical (S^2) minimum enclosing cap.

struct Cap3 {
  Vec3 c;
  double r;  // Chord radius
  bool ok;
};

inline Cap3 cap2(const Vec3& a, const Vec3& b) {
  const Vec3 s = {a.x + b.x, a.y + b.y, a.z + b.z};
  const double n = norm3(s);
  if (n < 1e-15) {
    return {{0.0, 0.0, 0.0}, 2.0, false};
  }
  const Vec3 u = {s.x / n, s.y / n, s.z / n};
  return {u, std::max(chord_dist(a, u), chord_dist(b, u)), true};
}

inline Cap3 max_chord_cap2(const Vec3& a, const Vec3& b, const Vec3& c) {
  Cap3 out = cap2(a, b);
  const Cap3 cb = cap2(b, c);
  if (!out.ok || (cb.ok && cb.r > out.r)) out = cb;
  const Cap3 ca = cap2(a, c);
  if (!out.ok || (ca.ok && ca.r > out.r)) out = ca;
  return out;
}

inline Cap3 cap3(const Vec3& a, const Vec3& b, const Vec3& c) {
  const Vec3 sum = {a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z};
  const double m_norm = norm3(sum);
  if (m_norm < 1e-15) return max_chord_cap2(a, b, c);
  const Vec3 m = {sum.x / m_norm, sum.y / m_norm, sum.z / m_norm};

  Vec3 e1, e2;
  if (!tangent_basis(m, e1, e2)) return max_chord_cap2(a, b, c);

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
  if (u_norm < 1e-15) return max_chord_cap2(a, b, c);
  u = {u.x / u_norm, u.y / u_norm, u.z / u_norm};
  const double r =
    std::max(chord_dist(a, u), std::max(chord_dist(b, u), chord_dist(c, u)));
  return {u, r, true};
}

inline bool outside_cap(const Vec3& p, const Cap3& cap, double tol) {
  return cap.ok && chord_dist(p, cap.c) > (cap.r + tol);
}

Cap3 welzl_cap(const std::vector<Vec3>& pts, double tol = 1e-13) {
  return welzl_ball<Vec3, Cap3>(
    pts,
    [](const Vec3& p) -> Cap3 { return {p, 0.0, true}; },
    cap2,
    cap3,
    [tol](const Vec3& p, const Cap3& cap) { return outside_cap(p, cap, tol); }
  );
}

struct Meb3Result {
  std::vector<double> lon;
  std::vector<double> lat;
  std::vector<Vec3> u;
  std::vector<double> radius;
  std::vector<uint8_t> valid;
  std::vector<uint8_t> non_empty;
};

Meb3Result compute_meb_sphere(const Rcpp::NumericVector& lon,
                              const Rcpp::NumericVector& lat,
                              const Rcpp::IntegerVector& feature_id,
                              int n_features,
                              double earth_radius) {
  Meb3Result res;
  res.lon.assign(n_features, NA_REAL);
  res.lat.assign(n_features, NA_REAL);
  res.u.assign(n_features, {NA_REAL, NA_REAL, NA_REAL});
  res.radius.assign(n_features, NA_REAL);
  res.valid.assign(n_features, 0);
  res.non_empty.assign(n_features, 0);

  std::vector<std::vector<Vec3> > grouped = group_feature_points<Vec3>(
    lon, lat, feature_id, n_features, ll_to_xyz
  );

  for (int i = 0; i < n_features; i++) {
    std::vector<Vec3>& pts = grouped[i];
    if (pts.empty()) continue;
    res.non_empty[i] = 1;
    golden_scramble(pts);

    const Cap3 cap = welzl_cap(pts);
    if (!cap.ok) continue;

    double max_angle = 0.0;
    for (const Vec3& p : pts) {
      max_angle = std::max(max_angle, gc_angle(p, cap.c));
    }
    if (!std::isfinite(max_angle) || max_angle >= 0.5 * kPi) continue;

    xyz_to_ll(cap.c, res.lon[i], res.lat[i]);
    res.u[i] = cap.c;
    res.radius[i] = max_angle * earth_radius;
    res.valid[i] = 1;
  }
  return res;
}

inline double sphere_cr_pair(const Meb3Result& a, int i,
                             const Meb3Result& b, int j,
                             double earth_radius) {
  return gc_angle(a.u[i], b.u[j]) * earth_radius +
         std::fabs(a.radius[i] - b.radius[j]);
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
  return Rcpp::List::create(
    Rcpp::Named("cx") = Rcpp::wrap(res.cx),
    Rcpp::Named("cy") = Rcpp::wrap(res.cy),
    Rcpp::Named("radius") = Rcpp::wrap(res.radius)
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
    compute_meb_sphere(lon, lat, feature_id, n_features, earth_radius);
  Rcpp::LogicalVector non_empty(n_features);
  for (int i = 0; i < n_features; i++) non_empty[i] = (res.non_empty[i] != 0);
  return Rcpp::List::create(
    Rcpp::Named("lon") = Rcpp::wrap(res.lon),
    Rcpp::Named("lat") = Rcpp::wrap(res.lat),
    Rcpp::Named("radius") = Rcpp::wrap(res.radius),
    Rcpp::Named("non_empty") = non_empty
  );
}

//' @title Planar CR distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cr_planar_cpp(const Rcpp::NumericVector& x,
                              const Rcpp::NumericVector& y,
                              const Rcpp::IntegerVector& feature_id,
                              int n_features) {
  const Meb2Result res = compute_meb_planar(x, y, feature_id, n_features);
  return sym_dist_matrix(
    n_features,
    [&](int i) { return res.valid[i] != 0; },
    [&](int i, int j) { return planar_cr_pair(res, i, res, j); }
  );
}

//' @title Planar CR cross-distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cr_cross_planar_cpp(const Rcpp::NumericVector& x1,
                                    const Rcpp::NumericVector& y1,
                                    const Rcpp::IntegerVector& fid1,
                                    int n1,
                                    const Rcpp::NumericVector& x2,
                                    const Rcpp::NumericVector& y2,
                                    const Rcpp::IntegerVector& fid2,
                                    int n2) {
  const Meb2Result r1 = compute_meb_planar(x1, y1, fid1, n1);
  const Meb2Result r2 = compute_meb_planar(x2, y2, fid2, n2);
  return cross_dist_matrix(
    n1, n2,
    [&](int i) { return r1.valid[i] != 0; },
    [&](int j) { return r2.valid[j] != 0; },
    [&](int i, int j) { return planar_cr_pair(r1, i, r2, j); }
  );
}

//' @title Spherical CR distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cr_sphere_cpp(const Rcpp::NumericVector& lon,
                              const Rcpp::NumericVector& lat,
                              const Rcpp::IntegerVector& feature_id,
                              int n_features,
                              double earth_radius = 6371010.0) {
  const Meb3Result res =
    compute_meb_sphere(lon, lat, feature_id, n_features, earth_radius);
  return sym_dist_matrix(
    n_features,
    [&](int i) { return res.valid[i] != 0; },
    [&](int i, int j) { return sphere_cr_pair(res, i, res, j, earth_radius); }
  );
}

//' @title Spherical CR cross-distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cr_cross_sphere_cpp(const Rcpp::NumericVector& lon1,
                                    const Rcpp::NumericVector& lat1,
                                    const Rcpp::IntegerVector& fid1,
                                    int n1,
                                    const Rcpp::NumericVector& lon2,
                                    const Rcpp::NumericVector& lat2,
                                    const Rcpp::IntegerVector& fid2,
                                    int n2,
                                    double earth_radius = 6371010.0) {
  const Meb3Result r1 = compute_meb_sphere(lon1, lat1, fid1, n1, earth_radius);
  const Meb3Result r2 = compute_meb_sphere(lon2, lat2, fid2, n2, earth_radius);
  return cross_dist_matrix(
    n1, n2,
    [&](int i) { return r1.valid[i] != 0; },
    [&](int j) { return r2.valid[j] != 0; },
    [&](int i, int j) { return sphere_cr_pair(r1, i, r2, j, earth_radius); }
  );
}
