#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

// Planar (R^2) Hausdorff distance between spatial units.
//
// Every distance in this file is in the units of the input coordinates.

namespace {

const double INF = std::numeric_limits<double>::infinity();
const double GOLDEN_RATIO = 0.6180339887498949;

inline double dist2_2d(double x1, double y1, double x2, double y2) {
  const double dx = x1 - x2;
  const double dy = y1 - y2;
  return dx * dx + dy * dy;
}

// A spatial unit, prepared once for all its distances.
//
// The input matrix has the columns x, y, and optionally a part id (see
// .to_list_mat() in R/dist-areal.R). A part is a ring, a line, or a point. An
// edge joins two consecutive rows only if they are in the same part, and a
// part with one row is a point. Without a part id, all rows are one ring, and
// the ring is closed if its first and last rows are different.
struct UnitR2 {
  std::vector<double> v;         // 2 * nr coordinates: (x, y)
  std::vector<int> order;        // Golden-ratio vertex traversal order
  std::vector<int> a, b;         // Edge endpoint indices (b == -1 for points)
  std::vector<double> dx, dy;    // Edge vector B - A
  std::vector<double> inv_len2;  // 1 / ||B - A||^2
  int nv() const { return static_cast<int>(order.size()); }
  int ne() const { return static_cast<int>(a.size()); }
};

void add_edge_r2(UnitR2& u, int i, int j) {
  double dx = 0.0;
  double dy = 0.0;
  double inv_len2 = 0.0;
  if (j >= 0) {
    dx = u.v[2 * j] - u.v[2 * i];
    dy = u.v[2 * j + 1] - u.v[2 * i + 1];
    const double len2 = dx * dx + dy * dy;
    if (len2 == 0.0) {
      j = -1;
    } else {
      inv_len2 = 1.0 / len2;
    }
  }
  u.a.push_back(i);
  u.b.push_back(j);
  u.dx.push_back(dx);
  u.dy.push_back(dy);
  u.inv_len2.push_back(inv_len2);
}

// A unit with a non-finite coordinate (for example, the NaN row that
// wk::wk_coords() gives for an EMPTY member of a MULTIPOINT) is empty, so its
// distances are NA, as in src/cr-dist.cpp.
UnitR2 make_unit_r2(const Eigen::Ref<const Eigen::MatrixXd>& poly) {
  UnitR2 u;
  const int nr = static_cast<int>(poly.rows());
  if (nr == 0) return u;
  for (int i = 0; i < nr; i++) {
    if (!std::isfinite(poly(i, 0)) || !std::isfinite(poly(i, 1))) return u;
  }

  u.v.resize(2 * nr);
  for (int i = 0; i < nr; i++) {
    u.v[2 * i] = poly(i, 0);
    u.v[2 * i + 1] = poly(i, 1);
  }
  const bool has_part = poly.cols() > 2;
  for (int i = 0; i < nr; i++) {
    const bool edge_next = i + 1 < nr && (!has_part || poly(i, 2) == poly(i + 1, 2));
    const bool edge_prev = i > 0 && (!has_part || poly(i, 2) == poly(i - 1, 2));
    if (edge_next) {
      add_edge_r2(u, i, i + 1);
    } else if (!edge_prev) {
      add_edge_r2(u, i, -1);  // A part with one row is a point
    }
  }
  // Closing edge if needed (only for input without a part id)
  if (!has_part && nr > 1 &&
      (poly(0, 0) != poly(nr - 1, 0) || poly(0, 1) != poly(nr - 1, 1))) {
    add_edge_r2(u, nr - 1, 0);
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
