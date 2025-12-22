#include <RcppEigen.h>
#include <functional>

// [[Rcpp::depends(RcppEigen)]]

// Define a type for our covariance kernel function for cleaner code
using CovKernel = std::function<double(double)>;

//--- cov functions for scalar distances ----

//--- * Matern Family ----
//' @rdname single-matern
//' @keywords internal
// [[Rcpp::export]]
double single_exp(double d, double sigsq, double phi) {
  if (d <= 0) return sigsq;
  return sigsq * exp(-d / phi);
}

//' @rdname single-matern
//' @keywords internal
// [[Rcpp::export]]
double single_matern3(double d, double sigsq, double phi) {
  if (d <= 0) return sigsq;
  const double d_phi = d / phi;
  return sigsq * (1.0 + d_phi) * exp(-d_phi);
}

//' @rdname single-matern
//' @keywords internal
// [[Rcpp::export]]
double single_matern5(double d, double sigsq, double phi) {
  if (d <= 0) return sigsq;
  const double d_phi = d / phi;
  return sigsq * (1.0 + d_phi + (d_phi * d_phi) / 3.0) * exp(-d_phi);
}

//' @title Matern covariance function (scalar - generic)
//'
//' @description Computing the Matern covariance function for a scalar distance,
//'   adapted from \code{geoR}.
//'
//' @details \code{single_matern3} and \code{single_matern5} are optimized for
//'   when \eqn{\nu} is 1.5 or 2.5, respectively. Similarly, \code{single_exp}
//'   and \code{single_gauss} represent the cases where \eqn{\nu = 0.5} or
//'   \eqn{\nu \to \infty}. In other words, they are the exponential and
//'   Gaussian covariance functions.
//' 
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the \eqn{\nu} parameter from the Matern covariance function,
//'   controls the differentiability of the process.
//' @name single-matern
//' 
//' @return a scalar representing the (matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_matern3}}, \code{\link{single_matern5}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_matern(double d, double sigsq, double phi, double nu) {
  if (d <= 0) return sigsq;
  const double d_phi = d / phi;
  return sigsq * (pow(2.0, 1.0 - nu) / ::tgamma(nu)) * pow(d_phi, nu) *
    Rf_bessel_k(d_phi, nu, 1.0);
}

//--- * Powered Exponential ----
//' @title Powered Exponential covariance function (scalar)
//'
//' @description Computing the Powered Exponential covariance function for a
//'   scalar distance.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Exponential covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Exponential covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the \eqn{\nu \in (0, 2]} parameter representing the "power"
//'
//' @return a scalar representing the (exponential) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{single_matern3}}, \code{\link{single_matern5}},
//'   \code{\link{mat_cov}}
//'
//' @keywords internal
// [[Rcpp::export]]
double single_pexp(double d, double sigsq, double phi, double nu) {
  if (d <= 0) return sigsq;
  return sigsq * exp(-pow(d / phi, nu));
}

//' @title Matern covariance function (scalar - generic)
//'
//' @description Computing the Matern covariance function for a scalar distance,
//'   adapted from \code{geoR}.
//'
//' @details \code{single_matern3} and \code{single_matern5} are optimized for
//'   when \eqn{\nu} is 1.5 or 2.5, respectively. Similarly, \code{single_exp}
//'   and \code{single_gauss} represent the cases where \eqn{\nu = 0.5} or
//'   \eqn{\nu \to \infty}. In other words, they are the exponential and
//'   Gaussian covariance functions.
//' 
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the \eqn{\nu} parameter from the Matern covariance function,
//'   controls the differentiability of the process.
//' @name single-matern
//' 
//' @return a scalar representing the (matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_matern3}}, \code{\link{single_matern5}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_dagum(double d, double sigsq, double phi, double nu) {
  if (d <= 0) return sigsq;
  const double d_phi = d / phi;
  return sigsq * (pow(2.0, 1.0 - nu) / ::tgamma(nu)) * pow(d_phi, nu) *
    Rf_bessel_k(d_phi, nu, 1.0);
}


//--- Internal generic functions ----

// Generic function to apply any covariance kernel to a distance matrix
Eigen::MatrixXd generic_cov_matrix(const Eigen::MatrixXd& dists,
                                   const CovKernel& kernel) {
  return dists.unaryExpr(kernel);
}

// Generic auxiliary function to compute the mean of a covariance matrix.
double generic_aux_mean(const Eigen::MatrixXd& dist,
                        const CovKernel& kernel) {
  if (dist.size() == 0) return 0.0;
  return generic_cov_matrix(dist, kernel).mean();
}

// Generalized function for computing polygon covariances.
Eigen::MatrixXd comp_cov_generic(const Rcpp::List& cross_dists,
                                 int n,
                                 int n2,
                                 const CovKernel& kernel) {
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n2);
  if (n == n2) {
    int k = 0;
    for (int j = 0; j < n; ++j) {
      for (int i = j; i < n; ++i) {
        Eigen::MatrixXd d_mat = Rcpp::as<Eigen::MatrixXd>(cross_dists[k]);
        double cov_val = generic_aux_mean(d_mat, kernel);
        out(i, j) = cov_val;
        out(j, i) = cov_val;
        k++;
      }
    }
  } else {
    for(int i = 0; i < cross_dists.size(); ++i) {
      Eigen::MatrixXd d_mat =
        Rcpp::as<Eigen::MatrixXd>(cross_dists[i]);
      *(out.data() + i) = generic_aux_mean(d_mat, kernel);
    }
  }
  return out;
}


//--- Exported fuctions ----

//' @title Matern covariance function for a given distance matrix.
//'
//' @description Computing the Matern covariance function for a matrix of
//'   distances.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the \eqn{\nu} parameter from the Matern covariance function,
//'   controls the differentiability of the process.
//' 
//' @return The matern covariance function (for a stationary and isotropic
//'   process) associated with the provided distances (\code{dists}) and the
//'   given set of parameters.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}}
//'
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd mat_cov(const Eigen::MatrixXd& dists,
                        double sigsq,
                        double phi,
                        double nu) {
  auto kernel = [=](double d) {
    if (nu == 0.5) return single_exp(d, sigsq, phi);
    if (nu == 1.5) return single_matern3(d, sigsq, phi);
    if (nu == 2.5) return single_matern5(d, sigsq, phi);
    return single_matern(d, sigsq, phi, nu);
  };
  return generic_cov_matrix(dists, kernel);
}

//' @title Matern covariance function for polygons.
//'
//' @description Computing the Matern covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @rdname mat_cov
//' @return The matern covariance matrix associated with a set of polygons.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd comp_mat_cov(const Rcpp::List& cross_dists,
                             int n,
                             int n2,
                             double sigsq,
                             double phi,
                             double nu) {
  auto kernel = [=](double d) {
    if (nu == 0.5) return single_exp(d, sigsq, phi);
    if (nu == 1.5) return single_matern3(d, sigsq, phi);
    if (nu == 2.5) return single_matern5(d, sigsq, phi);
    return single_matern(d, sigsq, phi, nu);
  };
  return comp_cov_generic(cross_dists, n, n2, kernel);
}

//' @title Powered Exponential Covariance
//' @rdname mat_cov
// [[Rcpp::export]]
Eigen::MatrixXd pexp_cov(const Eigen::MatrixXd& dists,
                         double sigsq,
                         double phi,
                         double nu) {
  auto kernel = [=](double d) { return single_pexp(d, sigsq, phi, nu); };
  return generic_cov_matrix(dists, kernel);
}

//' @title Powered Exponential Covariance for polygons
//' @rdname mat_cov
// [[Rcpp::export]]
Eigen::MatrixXd comp_pexp_cov(const Rcpp::List& cross_dists,
                              int n,
                              int n2,
                              double sigsq,
                              double phi,
                              double nu) {
  auto kernel = [=](double d) { return single_pexp(d, sigsq, phi, nu); };
  return comp_cov_generic(cross_dists, n, n2, kernel);
}
