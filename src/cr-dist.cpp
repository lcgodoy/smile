#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

// Center-radii (CR) distance between spatial units.
//
// The CR distance between units i and j is d(c_i, c_j) + |r_i - r_j|, where
// c and r are the center and radius of the minimum enclosing ball (the
// Chebyshev center and radius) of the vertices of a unit. Welzl's algorithm
// computes the exact ball. In the plane, d is the Euclidean distance, in the
// units of the coordinates. On the sphere, d is the great-circle distance, in
// meters, on a sphere of radius `earth_radius`.
//
// Input: the coordinates of all units (x and y, or longitude and latitude in
// degrees), with `feature_id` (1-based) as the unit of each row, as from
// wk::wk_coords().

namespace {

const double kPi = 3.14159265358979323846;
const double GOLDEN_RATIO = 0.6180339887498949;

inline bool is_invalid_scalar(double v) {
  return R_IsNA(v) || std::isnan(v) || !std::isfinite(v);
}

inline double deg2rad(double degrees) {
  return degrees * kPi / 180.0;
}

inline double rad2deg(double radians) {
  return radians * 180.0 / kPi;
}

// Deterministic golden-ratio permutation used by Welzl routines.
template <typename T>
void golden_scramble(std::vector<T>& pts) {
  const size_t n = pts.size();
  if (n <= 1) return;
  std::vector<std::pair<double, size_t> > key(n);
  for (size_t i = 0; i < n; i++) {
    key[i] = std::make_pair(std::fmod(static_cast<double>(i + 1) * GOLDEN_RATIO, 1.0), i);
  }
  std::sort(key.begin(), key.end());
  std::vector<T> out;
  out.reserve(n);
  for (size_t i = 0; i < n; i++) out.push_back(pts[key[i].second]);
  pts.swap(out);
}

// ---- Planar (R^2) geometry & Welzl MEB -------------------------------------

struct Pt2 {
  double x;
  double y;
};

struct Disk2 {
  double cx;
  double cy;
  double r2;
};

inline double dist2_2d(const Pt2& a, const Pt2& b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  return dx * dx + dy * dy;
}

inline Disk2 disk2(const Pt2& a, const Pt2& b) {
  const double cx = 0.5 * (a.x + b.x);
  const double cy = 0.5 * (a.y + b.y);
  const double dx = a.x - cx;
  const double dy = a.y - cy;
  return {cx, cy, dx * dx + dy * dy};
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
  const double dx = p.x - d.cx;
  const double dy = p.y - d.cy;
  const double r = std::sqrt(d.r2);
  const double lim = r + tol * (1.0 + r);
  return (dx * dx + dy * dy) > (lim * lim);
}

Disk2 welzl_disk(const std::vector<Pt2>& pts, double tol = 1e-12) {
  const int n = static_cast<int>(pts.size());
  auto with_2 = [&](int m, const Pt2& p, const Pt2& q) -> Disk2 {
    Disk2 d = disk2(p, q);
    for (int k = 0; k < m; k++) {
      if (outside_disk(pts[k], d, tol)) {
        d = disk3(p, q, pts[k]);
      }
    }
    return d;
  };
  auto with_1 = [&](int m, const Pt2& p) -> Disk2 {
    Disk2 d = {p.x, p.y, 0.0};
    for (int j = 0; j < m; j++) {
      if (outside_disk(pts[j], d, tol)) {
        d = with_2(j, p, pts[j]);
      }
    }
    return d;
  };
  Disk2 d = {pts[0].x, pts[0].y, 0.0};
  for (int i = 1; i < n; i++) {
    if (outside_disk(pts[i], d, tol)) {
      d = with_1(i, pts[i]);
    }
  }
  return d;
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
  if (n_features == 0) return res;

  const R_xlen_t m = feature_id.size();
  std::vector<int> counts(n_features, 0);
  std::vector<uint8_t> has_na(n_features, 0);
  for (R_xlen_t k = 0; k < m; k++) {
    const int fid = feature_id[k];
    if (fid == NA_INTEGER || fid < 1 || fid > n_features) {
      Rcpp::stop("feature_id contains out-of-range indices.");
    }
    const int idx = fid - 1;
    if (is_invalid_scalar(x[k]) || is_invalid_scalar(y[k])) {
      has_na[idx] = 1;
    } else {
      counts[idx]++;
    }
  }

  std::vector<std::vector<Pt2> > grouped(n_features);
  for (int i = 0; i < n_features; i++) {
    if (!has_na[i] && counts[i] > 0) grouped[i].reserve(counts[i]);
  }
  for (R_xlen_t k = 0; k < m; k++) {
    const int idx = feature_id[k] - 1;
    if (!has_na[idx]) {
      grouped[idx].push_back({x[k], y[k]});
    }
  }

  for (int i = 0; i < n_features; i++) {
    if (has_na[i] || grouped[i].empty()) continue;
    const auto& raw = grouped[i];
    if (raw.size() == 1) {
      res.cx[i] = raw[0].x;
      res.cy[i] = raw[0].y;
      res.radius[i] = 0.0;
      res.valid[i] = 1;
      continue;
    }

    // Shift by the first vertex to preserve floating-point precision on large
    // projected coordinates (e.g., UTM / Web Mercator).
    const double ox = raw[0].x;
    const double oy = raw[0].y;
    std::vector<Pt2> shifted;
    shifted.reserve(raw.size());
    for (const auto& p : raw) {
      shifted.push_back({p.x - ox, p.y - oy});
    }
    golden_scramble(shifted);
    const Disk2 d = welzl_disk(shifted);
    double max_r2 = 0.0;
    const Pt2 center_s = {d.cx, d.cy};
    for (const auto& p : shifted) {
      const double d2 = dist2_2d(p, center_s);
      if (d2 > max_r2) max_r2 = d2;
    }
    res.cx[i] = d.cx + ox;
    res.cy[i] = d.cy + oy;
    res.radius[i] = std::sqrt(max_r2);
    res.valid[i] = 1;
  }
  return res;
}

// ---- Spherical (S^2) geometry & Welzl MEC ----------------------------------

struct Vec3 {
  double x;
  double y;
  double z;
};

struct Cap3 {
  Vec3 c;
  double r;  // chord radius
  bool ok;
};

inline Vec3 ll_to_xyz(double lon_deg, double lat_deg) {
  const double lon = deg2rad(lon_deg);
  const double lat = deg2rad(lat_deg);
  const double clat = std::cos(lat);
  return {clat * std::cos(lon), clat * std::sin(lon), std::sin(lat)};
}

inline void xyz_to_ll(const Vec3& u, double& lon_deg, double& lat_deg) {
  lon_deg = rad2deg(std::atan2(u.y, u.x));
  lat_deg = rad2deg(std::atan2(u.z, std::sqrt(u.x * u.x + u.y * u.y)));
}

inline double dot3(const Vec3& a, const Vec3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vec3 cross3(const Vec3& a, const Vec3& b) {
  return {
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  };
}

inline double chord_dist(const Vec3& p, const Vec3& u) {
  const double dx = p.x - u.x;
  const double dy = p.y - u.y;
  const double dz = p.z - u.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline double gc_angle(const Vec3& u, const Vec3& v) {
  const Vec3 c = cross3(u, v);
  return std::atan2(std::sqrt(dot3(c, c)), dot3(u, v));
}

inline Cap3 cap2(const Vec3& a, const Vec3& b) {
  const Vec3 s = {a.x + b.x, a.y + b.y, a.z + b.z};
  const double norm = std::sqrt(dot3(s, s));
  if (norm < 1e-15) {
    return {{0.0, 0.0, 0.0}, 2.0, false};
  }
  const Vec3 u = {s.x / norm, s.y / norm, s.z / norm};
  const double r = std::max(chord_dist(a, u), chord_dist(b, u));
  return {u, r, true};
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
  const double m_norm = std::sqrt(dot3(sum, sum));
  if (m_norm < 1e-15) {
    return max_chord_cap2(a, b, c);
  }
  const Vec3 m = {sum.x / m_norm, sum.y / m_norm, sum.z / m_norm};
  const Vec3 ref = (std::fabs(m.z) < 0.9) ? Vec3{0.0, 0.0, 1.0} : Vec3{1.0, 0.0, 0.0};
  Vec3 e1 = cross3(ref, m);
  const double e1_norm = std::sqrt(dot3(e1, e1));
  if (e1_norm < 1e-15) {
    return max_chord_cap2(a, b, c);
  }
  e1 = {e1.x / e1_norm, e1.y / e1_norm, e1.z / e1_norm};
  const Vec3 e2 = cross3(m, e1);

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
  const double scale = (dot3(ba, ba) + dot3(ca, ca));

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
  const double u_norm = std::sqrt(dot3(u, u));
  if (u_norm < 1e-15) {
    return max_chord_cap2(a, b, c);
  }
  u = {u.x / u_norm, u.y / u_norm, u.z / u_norm};
  const double r = std::max(chord_dist(a, u), std::max(chord_dist(b, u), chord_dist(c, u)));
  return {u, r, true};
}

inline bool outside_cap(const Vec3& p, const Cap3& cap, double tol) {
  if (!cap.ok) return false;
  return chord_dist(p, cap.c) > (cap.r + tol);
}

Cap3 welzl_cap(const std::vector<Vec3>& pts, double tol = 1e-13) {
  const int n = static_cast<int>(pts.size());
  auto with_2 = [&](int m, const Vec3& p, const Vec3& q) -> Cap3 {
    Cap3 cap = cap2(p, q);
    for (int k = 0; k < m; k++) {
      if (outside_cap(pts[k], cap, tol)) {
        cap = cap3(p, q, pts[k]);
      }
    }
    return cap;
  };
  auto with_1 = [&](int m, const Vec3& p) -> Cap3 {
    Cap3 cap = {p, 0.0, true};
    for (int j = 0; j < m; j++) {
      if (outside_cap(pts[j], cap, tol)) {
        cap = with_2(j, p, pts[j]);
      }
    }
    return cap;
  };
  Cap3 cap = {pts[0], 0.0, true};
  for (int i = 1; i < n; i++) {
    if (outside_cap(pts[i], cap, tol)) {
      cap = with_1(i, pts[i]);
    }
  }
  return cap;
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
  if (n_features == 0) return res;

  const R_xlen_t m = feature_id.size();
  std::vector<int> counts(n_features, 0);
  std::vector<uint8_t> has_na(n_features, 0);
  for (R_xlen_t k = 0; k < m; k++) {
    const int fid = feature_id[k];
    if (fid == NA_INTEGER || fid < 1 || fid > n_features) {
      Rcpp::stop("feature_id contains out-of-range indices.");
    }
    const int idx = fid - 1;
    if (is_invalid_scalar(lon[k]) || is_invalid_scalar(lat[k])) {
      has_na[idx] = 1;
    } else {
      counts[idx]++;
    }
  }

  std::vector<std::vector<Vec3> > grouped(n_features);
  for (int i = 0; i < n_features; i++) {
    if (!has_na[i] && counts[i] > 0) {
      grouped[i].reserve(counts[i]);
      res.non_empty[i] = 1;
    }
  }
  for (R_xlen_t k = 0; k < m; k++) {
    const int idx = feature_id[k] - 1;
    if (!has_na[idx]) {
      grouped[idx].push_back(ll_to_xyz(lon[k], lat[k]));
    }
  }

  for (int i = 0; i < n_features; i++) {
    if (has_na[i] || grouped[i].empty()) continue;
    std::vector<Vec3> pts = grouped[i];
    golden_scramble(pts);
    const Cap3 cap = welzl_cap(pts);
    if (!cap.ok) continue;

    double max_angle = 0.0;
    for (const auto& p : pts) {
      const double ang = gc_angle(p, cap.c);
      if (ang > max_angle) max_angle = ang;
    }
    if (!std::isfinite(max_angle) || max_angle >= 0.5 * kPi) {
      continue;
    }

    double c_lon = NA_REAL, c_lat = NA_REAL;
    xyz_to_ll(cap.c, c_lon, c_lat);
    res.lon[i] = c_lon;
    res.lat[i] = c_lat;
    res.u[i] = cap.c;
    res.radius[i] = max_angle * earth_radius;
    res.valid[i] = 1;
  }
  return res;
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
  const Meb3Result res = compute_meb_sphere(lon, lat, feature_id, n_features, earth_radius);
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
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n_features, n_features);
  for (int i = 0; i < n_features; i++) {
    for (int j = i + 1; j < n_features; j++) {
      if (!res.valid[i] || !res.valid[j]) {
        out(i, j) = NA_REAL;
        out(j, i) = NA_REAL;
      } else {
        const double dx = res.cx[i] - res.cx[j];
        const double dy = res.cy[i] - res.cy[j];
        const double d = std::sqrt(dx * dx + dy * dy) +
          std::fabs(res.radius[i] - res.radius[j]);
        out(i, j) = d;
        out(j, i) = d;
      }
    }
  }
  return out;
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
  Eigen::MatrixXd out(n1, n2);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      if (!r1.valid[i] || !r2.valid[j]) {
        out(i, j) = NA_REAL;
      } else {
        const double dx = r1.cx[i] - r2.cx[j];
        const double dy = r1.cy[i] - r2.cy[j];
        out(i, j) = std::sqrt(dx * dx + dy * dy) +
          std::fabs(r1.radius[i] - r2.radius[j]);
      }
    }
  }
  return out;
}

//' @title Spherical CR distance matrix (internal)
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cr_sphere_cpp(const Rcpp::NumericVector& lon,
                        const Rcpp::NumericVector& lat,
                        const Rcpp::IntegerVector& feature_id,
                        int n_features,
                        double earth_radius = 6371010.0) {
  const Meb3Result res = compute_meb_sphere(lon, lat, feature_id, n_features, earth_radius);
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n_features, n_features);
  for (int i = 0; i < n_features; i++) {
    for (int j = i + 1; j < n_features; j++) {
      if (!res.valid[i] || !res.valid[j]) {
        out(i, j) = NA_REAL;
        out(j, i) = NA_REAL;
      } else {
        const double d = gc_angle(res.u[i], res.u[j]) * earth_radius +
          std::fabs(res.radius[i] - res.radius[j]);
        out(i, j) = d;
        out(j, i) = d;
      }
    }
  }
  return out;
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
  Eigen::MatrixXd out(n1, n2);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      if (!r1.valid[i] || !r2.valid[j]) {
        out(i, j) = NA_REAL;
      } else {
        out(i, j) = gc_angle(r1.u[i], r2.u[j]) * earth_radius +
          std::fabs(r1.radius[i] - r2.radius[j]);
      }
    }
  }
  return out;
}

//' @title Spherical CR distance with the exceedance bound (internal)
//' @noRd
// [[Rcpp::export]]
Eigen::MatrixXd cr_sphere_meb_cpp(const Rcpp::NumericVector& lon,
                            const Rcpp::NumericVector& lat,
                            const Rcpp::IntegerVector& feature_id,
                            int n_features,
                            bool exceedance = false,
                            double earth_radius = 6371010.0) {
  const Meb3Result res = compute_meb_sphere(lon, lat, feature_id, n_features, earth_radius);
  const double half_circ = kPi * earth_radius;
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n_features, n_features);

  for (int i = 0; i < n_features; i++) {
    if (exceedance) {
      out(i, i) = res.valid[i] ? (half_circ - res.radius[i]) : NA_REAL;
    } else if (!res.valid[i]) {
      out(i, i) = NA_REAL;
    }
    for (int j = i + 1; j < n_features; j++) {
      if (!res.valid[i] || !res.valid[j]) {
        out(i, j) = NA_REAL;
        out(j, i) = NA_REAL;
      } else {
        const double exc = half_circ - std::min(res.radius[i], res.radius[j]);
        double val = exc;
        if (!exceedance) {
          const double cr_val = gc_angle(res.u[i], res.u[j]) * earth_radius +
            std::fabs(res.radius[i] - res.radius[j]);
          val = std::min(exc, cr_val);
        }
        out(i, j) = val;
        out(j, i) = val;
      }
    }
  }
  return out;
}

//' @title Spherical CR upper bound, exceedance bound, and CR (internal)
//' @noRd
// [[Rcpp::export]]
Rcpp::List cr_all_cpp(const Rcpp::NumericVector& lon,
                      const Rcpp::NumericVector& lat,
                      const Rcpp::IntegerVector& feature_id,
                      int n_features,
                      double earth_radius = 6371010.0) {
  const Meb3Result res = compute_meb_sphere(lon, lat, feature_id, n_features, earth_radius);
  const double half_circ = kPi * earth_radius;
  Eigen::MatrixXd upper = Eigen::MatrixXd::Zero(n_features, n_features);
  Eigen::MatrixXd exc = Eigen::MatrixXd::Zero(n_features, n_features);
  Eigen::MatrixXd cr_mat = Eigen::MatrixXd::Zero(n_features, n_features);

  for (int i = 0; i < n_features; i++) {
    exc(i, i) = res.valid[i] ? (half_circ - res.radius[i]) : NA_REAL;
    if (!res.valid[i]) {
      cr_mat(i, i) = NA_REAL;
    }
    for (int j = i + 1; j < n_features; j++) {
      if (!res.valid[i] || !res.valid[j]) {
        upper(i, j) = NA_REAL;
        upper(j, i) = NA_REAL;
        exc(i, j) = NA_REAL;
        exc(j, i) = NA_REAL;
        cr_mat(i, j) = NA_REAL;
        cr_mat(j, i) = NA_REAL;
      } else {
        const double u_val = gc_angle(res.u[i], res.u[j]) * earth_radius +
          std::fabs(res.radius[i] - res.radius[j]);
        const double e_val = half_circ - std::min(res.radius[i], res.radius[j]);
        const double c_val = std::min(e_val, u_val);
        upper(i, j) = u_val;
        upper(j, i) = u_val;
        exc(i, j) = e_val;
        exc(j, i) = e_val;
        cr_mat(i, j) = c_val;
        cr_mat(j, i) = c_val;
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("upper") = upper,
    Rcpp::Named("exceedance") = exc,
    Rcpp::Named("cr") = cr_mat
  );
}
