#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

// Spherical (S^2) Hausdorff distance between spatial units.
//
// Edges are great-circle arcs, as in s2. Every distance in this file is in
// meters, on a sphere of radius `earth_radius` (the default is
// s2::s2_earth_radius_meters()).

namespace {

const double kPi = 3.14159265358979323846;
const double INF = std::numeric_limits<double>::infinity();
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

// Squared distance between two 3D vectors
inline double dist2(const double* u, const double* v) {
  const double d0 = u[0] - v[0];
  const double d1 = u[1] - v[1];
  const double d2 = u[2] - v[2];
  return d0 * d0 + d1 * d1 + d2 * d2;
}

// 1. A spatial unit, prepared once for all its distances.
//
// The input matrix has the columns longitude, latitude, and optionally a part
// id (see .to_list_mat() in R/dist-areal.R). A part is a ring, a line, or a
// point. An edge joins two consecutive rows only if they are in the same
// part, and a part with one row is a point. Without a part id, all rows are
// one ring, and the ring is closed if its first and last rows are different.
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
  int nv() const { return static_cast<int>(order.size()); }
  int ne() const { return static_cast<int>(a.size()); }
};

void add_edge(Unit& u, int i, int j) {
  double nrm[3] = {0.0, 0.0, 0.0};
  double ta[3] = {0.0, 0.0, 0.0};
  double tb[3] = {0.0, 0.0, 0.0};
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
    }
  }
  u.a.push_back(i);
  u.b.push_back(j);
  u.n.insert(u.n.end(), nrm, nrm + 3);
  u.ta.insert(u.ta.end(), ta, ta + 3);
  u.tb.insert(u.tb.end(), tb, tb + 3);
}

// A unit with a non-finite coordinate (for example, the NaN row that
// wk::wk_coords() gives for an EMPTY member of a MULTIPOINT) is empty, so its
// distances are NA, as in src/cr-dist.cpp.
Unit make_unit(const Eigen::Ref<const Eigen::MatrixXd>& poly) {
  Unit u;
  const int nr = static_cast<int>(poly.rows());
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
  for (int i = 0; i < nr; i++) {
    const bool edge_next = i + 1 < nr && (!has_part || poly(i, 2) == poly(i + 1, 2));
    const bool edge_prev = i > 0 && (!has_part || poly(i, 2) == poly(i - 1, 2));
    if (edge_next) {
      add_edge(u, i, i + 1);
    } else if (!edge_prev) {
      add_edge(u, i, -1);  // A part with one row is a point
    }
  }
  // Closing edge if needed (only for input without a part id)
  if (!has_part && nr > 1 &&
      (poly(0, 0) != poly(nr - 1, 0) || poly(0, 1) != poly(nr - 1, 1))) {
    add_edge(u, nr - 1, 0);
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
  // P projects on the inside of the arc: the closest point is on the arc, at
  // the cross-track angle t, with sin(t) = |P . n| and cos(t) = |n x P|. The
  // cross product keeps cos(t) precise near 90 degrees, where
  // sqrt(1 - sin(t)^2) and asin() lose about 1e-8 radians.
  if (dot3(&y.ta[3 * k], P) >= 0.0 && dot3(&y.tb[3 * k], P) >= 0.0) {
    const double* N = &y.n[3 * k];
    const double s = std::fabs(dot3(P, N));
    double c[3];
    cross3(N, P, c);
    kind = 2;
    return 2.0 * s * s / (1.0 + std::sqrt(dot3(c, c)));
  }
  // Otherwise, the closest point is one of the endpoints.
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
  // atan2(norm(cross), dot) is highly accurate for all angles
  return std::atan2(std::sqrt(dot3(c, c)), dot3(P, Q));
}

// 4. Directed Hausdorff distance from the vertices of x to the edges of y,
// with the early break of Taha and Hanbury (2015, IEEE TPAMI 37(11)).
// `max_c2` and `max_angle` hold the largest distance so far, as a squared
// chord and as an angle. If a vertex has an edge closer than that, the vertex
// cannot change the maximum, so the scan of its edges stops there. Each scan
// starts at the edge that ended the previous scan, which is often close. The
// result is exact. Both units must have at least one vertex (see
// haus_units()).
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
    // P is farther from y than every vertex before it: a new maximum
    start = best_k;
    max_c2 = best;
    max_angle = std::max(max_angle, edge_angle(P, y, best_k, best_kind));
  }
}

// 5. Hausdorff distance between two units, in meters. The second direction
// starts from the maximum of the first, so it also stops early. An empty unit
// gives NA, as in the planar code and in sf.
double haus_units(const Unit& x, const Unit& y, double earth_radius) {
  if (x.nv() == 0 || y.nv() == 0) return NA_REAL;
  double max_c2 = -1.0;
  double max_angle = 0.0;
  directed_haus(x, y, max_c2, max_angle);
  directed_haus(y, x, max_c2, max_angle);
  return max_angle * earth_radius;
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
  // Each row of pred_mat is a unit with one point
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
