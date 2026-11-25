#ifndef SMILE_GEOM_COMMON_H_
#define SMILE_GEOM_COMMON_H_

#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

// Geometry primitives, golden-ratio ordering, and distance-matrix builders
// shared by src/cr-dist.cpp, src/haus-r2.cpp, and src/haus-s2.cpp.

namespace smile_geom {

const double kPi = 3.14159265358979323846;
const double kInf = std::numeric_limits<double>::infinity();
const double kEps = std::numeric_limits<double>::epsilon();
const double kGoldenRatio = 0.6180339887498949;

inline bool is_invalid_scalar(double v) {
  return R_IsNA(v) || std::isnan(v) || !std::isfinite(v);
}

inline double deg2rad(double degrees) {
  return degrees * kPi / 180.0;
}

inline double rad2deg(double radians) {
  return radians * 180.0 / kPi;
}

// Converts a central angle (radians) to its squared 3D chord length.
inline double angle_to_chord2(double ang) {
  if (ang <= 0.0) return 0.0;
  if (ang >= kPi) return 4.0;
  const double s = std::sin(0.5 * ang);
  return 4.0 * s * s;
}

struct Pt2 {
  double x;
  double y;
};

inline double dist2_2d(double x1, double y1, double x2, double y2) {
  const double dx = x1 - x2;
  const double dy = y1 - y2;
  return dx * dx + dy * dy;
}

inline double dist2_2d(const Pt2& a, const Pt2& b) {
  return dist2_2d(a.x, a.y, b.x, b.y);
}

struct Vec3 {
  double x;
  double y;
  double z;
};

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

inline double norm3(const Vec3& u) {
  return std::sqrt(dot3(u, u));
}

inline bool normalize3(Vec3& u) {
  const double n = norm3(u);
  if (n <= 0.0 || !std::isfinite(n)) return false;
  u.x /= n;
  u.y /= n;
  u.z /= n;
  return true;
}

inline double dist2_3d(const Vec3& a, const Vec3& b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

inline double chord_dist(const Vec3& a, const Vec3& b) {
  return std::sqrt(dist2_3d(a, b));
}

inline double gc_angle(const Vec3& u, const Vec3& v) {
  const Vec3 c = cross3(u, v);
  return std::atan2(norm3(c), dot3(u, v));
}

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

// Builds an orthonormal basis (e1, e2) for the plane tangent to S^2 at unit
// vector `c`.
inline bool tangent_basis(const Vec3& c, Vec3& e1, Vec3& e2) {
  const Vec3 ref = (std::fabs(c.z) < 0.9) ? Vec3{0.0, 0.0, 1.0}
                                          : Vec3{1.0, 0.0, 0.0};
  e1 = cross3(ref, c);
  if (!normalize3(e1)) return false;
  e2 = cross3(c, e1);
  return true;
}

inline bool vec3_less(const Vec3& a, const Vec3& b) {
  if (a.x != b.x) return a.x < b.x;
  if (a.y != b.y) return a.y < b.y;
  return a.z < b.z;
}

// Deterministic golden-ratio index order over 0, ..., n - 1. Visiting spread
// vertices early accelerates both Welzl's algorithm and Taha-Hanbury pruning.
inline std::vector<int> golden_order(int n) {
  std::vector<int> order(n);
  if (n <= 0) return order;
  std::vector<std::pair<double, int> > key(n);
  for (int i = 0; i < n; i++) {
    key[i] = std::make_pair(
      std::fmod(static_cast<double>(i + 1) * kGoldenRatio, 1.0), i
    );
  }
  std::sort(key.begin(), key.end());
  for (int i = 0; i < n; i++) order[i] = key[i].second;
  return order;
}

template <typename T>
void golden_scramble(std::vector<T>& pts) {
  const int n = static_cast<int>(pts.size());
  if (n <= 1) return;
  const std::vector<int> ord = golden_order(n);
  std::vector<T> out;
  out.reserve(n);
  for (int i = 0; i < n; i++) out.push_back(pts[ord[i]]);
  pts.swap(out);
}

// Builds a symmetric n x n distance matrix, setting rows and columns of
// invalid units (including the diagonal) to NA.
template <typename ValidFn, typename PairFn>
Eigen::MatrixXd sym_dist_matrix(int n, ValidFn is_valid, PairFn pair_dist) {
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n);
  for (int i = 0; i < n; i++) {
    Rcpp::checkUserInterrupt();
    const bool ok_i = is_valid(i);
    if (!ok_i) out(i, i) = NA_REAL;
    for (int j = i + 1; j < n; j++) {
      const double d = (ok_i && is_valid(j)) ? pair_dist(i, j) : NA_REAL;
      out(i, j) = d;
      out(j, i) = d;
    }
  }
  return out;
}

// Builds an n1 x n2 cross-distance matrix, setting entries with an invalid row
// or column unit to NA.
template <typename Valid1Fn, typename Valid2Fn, typename PairFn>
Eigen::MatrixXd cross_dist_matrix(int n1, int n2,
                                  Valid1Fn valid1, Valid2Fn valid2,
                                  PairFn pair_dist) {
  Eigen::MatrixXd out(n1, n2);
  for (int i = 0; i < n1; i++) {
    Rcpp::checkUserInterrupt();
    const bool ok_i = valid1(i);
    for (int j = 0; j < n2; j++) {
      out(i, j) = (ok_i && valid2(j)) ? pair_dist(i, j) : NA_REAL;
    }
  }
  return out;
}

}  // namespace smile_geom

#endif  // SMILE_GEOM_COMMON_H_
