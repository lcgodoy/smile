#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <numeric>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double mean_mat(arma::mat mat_aux) {
  return std::accumulate(mat_aux.begin(), mat_aux.end(), 0.0) /
    (mat_aux.n_rows * mat_aux.n_cols);
}

//' @title Matern covariance function (scalar - generic)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, adapted from \code{geoR}.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the \eqn{\nu} parameter from the Matern covariance function,
//'   controls the differentiability of the process.
//'
//' @return a scalar representing the (matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_matern}}, \code{\link{single_matern5}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_matern(double d, double sigsq,
		     double phi, double nu) {
  double out = sigsq;
  if(d > 0) {
    out *=
      (pow(2, 1 - nu) / ::tgamma(nu)) *
      pow(d / phi, nu) *
      Rf_bessel_k(d / phi, nu, 1);
  } 
  return out;
}

//' @title Matern covariance function (scalar - nu = 3/2)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, with \eqn{\nu = 3/2}.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//'
//' @return a scalar representing the (matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_matern}}, \code{\link{single_matern5}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_matern3(double d, double sigsq,
		      double phi) {
  double out = sigsq;
  if(d > 0) {
    out *= (1 + d/phi) * exp( - d / phi);
  } 
  return out;
}

//' @title Matern covariance function (scalar - nu = 5/2)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, with \eqn{\nu = 5/2}.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//'
//' @return a scalar representing the (matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' 
//' @seealso \code{\link{single_matern}}, \code{\link{single_matern3}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_matern5(double d, double sigsq,
		      double phi) {
  double out = sigsq;
  if(d > 0) {
    out *= (1 + (d/phi) + (pow(d / phi, 2)/3)) * exp( - d / phi);
  } 
  return out;
}

//' @title Exponential covariance function (scalar)
//'
//' @description Computing the Exponential covariance function for a single
//'   distance measure.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Exponential covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Exponential covariance function,
//'   controls the range of the spatial dependence.
//'
//' @return a scalar representing the (exponential) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' 
//' @seealso \code{\link{single_pexp}}, \code{\link{single_matern}},
//'   \code{\link{single_matern3}}, \code{\link{single_matern5}},
//'   \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_exp(double d, double sigsq, double phi) {
  double out = sigsq;
  if(d > 0) {
    out *= exp(- d / phi);
  }
  return out;
}

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
arma::mat mat_cov(const arma::mat& dists, double sigsq,
		  double phi, double nu) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  if(nu == 0.5) {
    std::transform(dists.begin(), dists.end(),
		   out.begin(), std::bind(&single_exp, std::placeholders::_1,
					  sigsq, phi));
  } else if(nu == 1.5) {
    std::transform(dists.begin(), dists.end(),
		   out.begin(), std::bind(&single_matern3, std::placeholders::_1,
					  sigsq, phi));
  } else if(nu == 2.5) {
    std::transform(dists.begin(), dists.end(),
		   out.begin(), std::bind(&single_matern5, std::placeholders::_1,
					  sigsq, phi));
  } else {
    std::transform(dists.begin(), dists.end(),
		   out.begin(), std::bind(&single_matern, std::placeholders::_1,
					  sigsq, phi, nu));
  }
  return out;
}

// This one is not as fast as `cov_mat_cpp`
// arma::mat cov_mat_cpp2(arma::mat dists,
// 		      double sigsq,
// 		      double phi,
// 		      double kappa) {
//   arma::mat out = dists;
//   if(kappa == 0.5) {
//     out.transform([=](double val) {return(single_exp(val, sigsq, phi)); });
//   } else {
//     out.transform([=](double val) {return(single_matern(val, sigsq, phi, kappa)); });
//   }
//   return out;
// }

//' @title Mean of a (Matern) covariance function (Internal use)
//'
//' @description This is an auxiliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polygon and speed-up the computations.
//'
//' @param dist a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the \eqn{\nu} parameter from the Matern covariance function,
//'   controls the differentiability of the process.
//' 
//' @return The mean of \code{mat_cov(dist, sigsq, phi, nu)}.
//' @keywords internal
// [[Rcpp::export]]
double aux_matern(arma::mat dist, double sigsq,
		  double phi, double nu) {
  return mean_mat(mat_cov(dist, sigsq, phi, nu));
}

//' @title Matern covariance function for a polygons.
//'
//' @description Computing the Matern covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an integer representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the \eqn{\nu} parameter from the Matern covariance function,
//'   controls the differentiability of the process. Note that, if we set
//'   \code{nu = .5}, then the calculations are based on the exponential
//'   covariance function.
//' 
//' @return The matern covariance matrix associated with a set of polygons.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_mat_cov(const List& cross_dists, int n,
		       int n2, double sigsq, double phi,
		       double nu) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(), cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_matern, std::placeholders::_1,
			       sigsq, phi, nu));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(), cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_matern, std::placeholders::_1,
			     sigsq, phi, nu));
  }
  return out;
}

//' @title Powered Exponential covariance function (scalar)
//'
//' @description Computing the Powered Exponential covariance function for a
//'   single distance measure.
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
  double out = sigsq;
  if(d > 0) {
    out *= exp( - pow( (d / phi) , nu));
  }
  return out;
}

//' @title Powered Exponential covariance function for a given distance matrix.
//'
//' @description Computing the Matern covariance function for a matrix of
//'   distances.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Exponential covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Exponential covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the \eqn{\nu \in (0, 2]} parameter representing the "power"
//' 
//' @return The powered exponential covariance function (for a stationary and
//'   isotropic process) associated with the provided distances (\code{dists})
//'   and the given set of parameters.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}}
//' @keywords internal
// [[Rcpp::export]]
arma::mat pexp_cov(const arma::mat& dists, double sigsq,
		   double phi, double nu) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  std::transform(dists.begin(), dists.end(),
		 out.begin(), std::bind(&single_pexp, std::placeholders::_1,
					sigsq, phi, nu));
  return out;
}

//' @title Mean of a (Powered Exponential) covariance function (Internal use)
//'
//' @description This is an auxiliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polygon and speed-up the computations.
//'
//' @param dist a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Powered Exponential
//'   covariance function.
//' @param phi the \eqn{\phi} parameter from the Powered Exponential covariance
//'   function, controls the range of the spatial dependence.
//' @param nu the \eqn{\nu} parameter from the Powered Exponential
//'   covariance function,
//'   controls the differentiability of the process.
//' 
//' @return The mean of \code{pexp_cov(dist, sigsq, phi, nu)}.
//' @keywords internal
// [[Rcpp::export]]
double aux_pexp(arma::mat dist, double sigsq,
	       double phi, double nu) {
  return mean_mat(pexp_cov(dist, sigsq, phi, nu));
}

//' @title Powered Exponential covariance function for a polygons.
//'
//' @description Computing the Powered Exponential covariance function between
//'   polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an integer representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Powered Exponential
//'   covariance function.
//' @param phi the \eqn{\phi} parameter from the Powered Exponential covariance
//'   function, controls the range of the spatial dependence.
//' @param nu the \eqn{\nu \in (0, 2]} parameter representing the "power"
//' 
//' @return The powered exponential covariance matrix associated with a set of
//'   polygons.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{mat_cov}}
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_pexp_cov(const List& cross_dists, int n,
			int n2, double sigsq, double phi, 
			double nu) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(), cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_pexp, std::placeholders::_1,
			       sigsq, phi, nu));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(), cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_pexp, std::placeholders::_1,
			     sigsq, phi, nu));
  }
  return out;
}

//' @title Gaussian covariance function (scalar)
//'
//' @description Computing the Gaussian covariance function for a single
//'   distance measure.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Gaussian covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Gaussian covariance function,
//'   controls the range of the spatial dependence.
//'
//' @return a scalar representing the (gaussian) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{single_matern3}}, \code{\link{single_matern5}},
//'   \code{\link{mat_cov}}
//'
//' @keywords internal
// [[Rcpp::export]]
double single_gauss(double d, double sigsq, double phi) {
  double out = sigsq;
  if(d > 0) {
    out *= exp( - .5 * pow( (d / phi) , 2));
  }
  return out;
}

//' @title Computing the Gaussian covariance function for a single distance
//'   measure.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Gaussian covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Gaussian covariance function,
//'   controls the range of the spatial dependence.
//'
//' @return The Gaussian covariance function (for a stationary and
//'   isotropic process) associated with the provided distances (\code{dists})
//'   and the given set of parameters.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}}
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat gauss_cov(const arma::mat& dists, double sigsq, double phi) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  std::transform(dists.begin(), dists.end(),
		 out.begin(), std::bind(&single_gauss, std::placeholders::_1,
					sigsq, phi));
  return out;
}


//' @title Mean of a (Gaussian) covariance function (Internal use)
//'
//' @description This is an auxiliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polygon and speed-up the computations.
//'
//' @param dist a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Gaussian covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Gaussian covariance function,
//'   controls the range of the spatial dependence.
//' 
//' @return The mean of \code{gauss_cov(dist, sigsq, phi)}.
//'
//' @keywords internal
// [[Rcpp::export]]
double aux_gauss(arma::mat dist, double sigsq, double phi) {
  return mean_mat(gauss_cov(dist, sigsq, phi));
}

//' @title Gaussian covariance function for a polygons.
//'
//' @description Computing the Gaussian covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an integer representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Gaussian covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Gaussian covariance function,
//'   controls the range of the spatial dependence.
//' 
//' @return The gaussian covariance matrix associated with a set of
//'   polygons.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_gauss_cov(const List& cross_dists, int n,
			 int n2, double sigsq, double phi) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(), cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_gauss, std::placeholders::_1,
			       sigsq, phi));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(), cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_gauss, std::placeholders::_1,
			     sigsq, phi));
  }
  return out;
}

//' @title Spherical covariance function (scalar)
//'
//' @description Computing the Spherical covariance function for a single
//'   distance measure.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Spherical covariance.
//'   function.
//' @param phi the \eqn{\phi} parameter from the Spherical covariance function,
//'   controls the range of the spatial dependence.
//'
//' @return a scalar representing the (gaussian) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{single_matern3}}, \code{\link{single_matern5}},
//'   \code{\link{mat_cov}}
//'
//' @keywords internal
// [[Rcpp::export]]
double single_spher(double d, double sigsq, double phi) {
  double out = sigsq;
  if(d > phi) {
    out = 0;
  } else if(d > 0) {
    out *= (1 - 1.5 * (d / phi) + .5 * pow(d / phi, 3));
  }
  return out;
}

//' @title Computing the Spherical covariance function for a single distance
//'   measure.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Spherical covariance.
//'   function.
//' @param phi the \eqn{\phi} parameter from the Spherical covariance function,
//'   controls the range of the spatial dependence.
//''
//' @return The Spherical covariance function (for a stationary and
//'   isotropic process) associated with the provided distances (\code{dists})
//'   and the given set of parameters.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}}
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat spher_cov(const arma::mat& dists, double sigsq, double phi) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  std::transform(dists.begin(), dists.end(),
		 out.begin(), std::bind(&single_spher, std::placeholders::_1,
					sigsq, phi));
  return out;
}

//' @title Mean of a (Spherical) covariance function (Internal use)
//'
//' @description This is an auxiliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polygon and speed-up the computations.
//'
//' @param dist a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Spherical covariance.
//'   function.
//' @param phi the \eqn{\phi} parameter from the Spherical covariance function,
//'   controls the range of the spatial dependence.
//' 
//' @return The mean of \code{spher_cov(dist, sigsq, phi)}.
//'
//' @keywords internal
// [[Rcpp::export]]
double aux_spher(arma::mat dist, double sigsq, double phi) {
  return mean_mat(spher_cov(dist, sigsq, phi));
}

//' @title Spherical covariance function for a polygons.
//'
//' @description Computing the Spherical covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an integer representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Spherical covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Spherical covariance function,
//'   controls the range of the spatial dependence.
//' 
//' @return The spherical covariance matrix associated with a set of
//'   polygons.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{mat_cov}}
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_spher_cov(const List& cross_dists, int n,
			 int n2, double sigsq, double phi) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(), cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_spher, std::placeholders::_1,
			       sigsq, phi));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(), cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_spher, std::placeholders::_1,
			     sigsq, phi));
  }
  return out;
}

//' @title Cubic spline covariance function (scalar)
//'
//' @description Computing the Spherical covariance function for a single
//'   distance measure.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Spherical covariance.
//'   function.
//' @param phi the \eqn{\phi} parameter from the Spherical covariance function,
//'   controls the range of the spatial dependence.
//'
//' @return a scalar representing the (gaussian) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{single_matern3}}, \code{\link{single_matern5}},
//'   \code{\link{mat_cov}}
//'
//' @keywords internal
// [[Rcpp::export]]
double single_cs(double d, double sigsq, double phi) {
  double out = sigsq;
  if(d < .5 * phi) {
    out *= (1 - 6 * pow(d / phi, 2) + 6 * pow(d / phi, 3));
  } else if(.5 * phi <= d && d < phi) {
    out *= (2 * pow(1 - (d / phi), 3));
  } else {
    out = 0;
  }
  return out;
}

//' @title Computing the Cubic spline covariance function for a single distance
//'   measure.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Spherical covariance.
//'   function.
//' @param phi the \eqn{\phi} parameter from the Spherical covariance function,
//'   controls the range of the spatial dependence.
//''
//' @return The Spherical covariance function (for a stationary and
//'   isotropic process) associated with the provided distances (\code{dists})
//'   and the given set of parameters.
//'
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat cs_cov(const arma::mat& dists, double sigsq, double phi) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  std::transform(dists.begin(), dists.end(),
		 out.begin(), std::bind(&single_cs, std::placeholders::_1,
					sigsq, phi));
  return out;
}

//' @title Mean of a (Cubic spline) covariance function (Internal use)
//'
//' @description This is an auxiliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polygon and speed-up the computations.
//'
//' @param dist a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Spherical covariance.
//'   function.
//' @param phi the \eqn{\phi} parameter from the Spherical covariance function,
//'   controls the range of the spatial dependence.
//' 
//' @return The mean of \code{spher_cov(dist, sigsq, phi)}.
//'
//' @keywords internal
// [[Rcpp::export]]
double aux_cs(arma::mat dist, double sigsq, double phi) {
  return mean_mat(cs_cov(dist, sigsq, phi));
}

//' @title Cubic spline covariance function for a polygons.
//'
//' @description Computing the Spherical covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an integer representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Spherical covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Spherical covariance function,
//'   controls the range of the spatial dependence.
//' 
//' @return The spherical covariance matrix associated with a set of
//'   polygons.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}},
//'   \code{\link{mat_cov}}
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_cs_cov(const List& cross_dists, int n,
		      int n2, double sigsq, double phi) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(), cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_cs, std::placeholders::_1,
			       sigsq, phi));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(), cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_cs, std::placeholders::_1,
			     sigsq, phi));
  }
  return out;
}

//' @title Matern Generalized Wendland (GW) covariance function with kappa = 0
//'   (scalar - generic)
//'
//' @description adapted from Bevilacqua et al. 2019.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param mu a parameter that controls the smoothness of the covariance
//'   function. Note that, \eqn{\mu \geq 1}.
//' 
//' @return a scalar representing the GW covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_gw0(double d, double sigsq, double phi,
		  double mu) {
  double out = sigsq, aux = d / phi,
    beta = mu + .5;
  if(aux < 1) {
    out *= pow(1 - aux, beta);
  } else {
    out *= 0;
  }
  return out;
}

//' @title Matern Generalized Wendland (GW) covariance function with kappa = 1
//'   (scalar - generic)
//'
//' @description adapted from Bevilacqua et al. 2019.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param mu a parameter that controls the smoothness of the covariance
//'   function. Note that, \eqn{\mu \geq 1}.
//' 
//' @return a scalar representing the GW covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_gw1(double d, double sigsq, double phi,
		  double mu) {
  double out = sigsq, aux = d / phi,
    beta = mu + 2.5;
  if(aux < 1) {
    out *= (1 + beta * aux) * pow(1 - aux, beta);
  } else {
    out *= 0;
  }
  return out;
}

//' @title Matern Generalized Wendland (GW) covariance function with kappa = 2
//'   (scalar - generic)
//'
//' @description adapted from Bevilacqua et al. 2019.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param mu a parameter that controls the smoothness of the covariance
//'   function. Note that, \eqn{\mu \geq 1}.
//' 
//' @return a scalar representing the GW covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_gw2(double d, double sigsq, double phi,
		  double mu) {
  double out = sigsq, aux = d / phi,
    beta = mu + 4.5;
  if(aux < 1) {
    out *= pow(1 - aux, beta) *
      (1 + beta * aux +
       ((beta * beta - 1) *
	aux * aux / 3));
  } else {
    out *= 0;
  }
  return out;
}

//' @title Matern Generalized Wendland (GW) covariance function with kappa = 3
//'   (scalar - generic)
//'
//' @description adapted from Bevilacqua et al. 2019.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param mu a parameter that controls the smoothness of the covariance
//'   function. Note that, \eqn{\mu \geq 1}.
//' 
//' @return a scalar representing the GW covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_gw3(double d, double sigsq, double phi,
		  double mu) {
  double out = sigsq, aux = d / phi,
    beta = mu + 6.5;
  if(aux < 1) {
    out *= pow(1 - aux, beta) *
      (1 + beta * aux + ((2 * beta * beta - 3) * aux * aux * .2 ) +
       ((beta * beta - 4) * beta * aux * aux * aux / 15));
  } else {
    out *= 0;
  }
  return out;
}

//' @title Matern Generalized Wendland (GW) covariance function with kappa = 0
//'   (scalar - generic)
//'
//' @description adapted from Bevilacqua et al. 2019.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param kappa \eqn{\kappa \in \{0, \ldots, 3 \}}.
//' @param mu a parameter that controls the smoothness of the covariance
//'   function. Note that, \eqn{\mu \geq 1}.
//' 
//' @return a scalar representing the GW covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_gw(double d, double sigsq, double phi,
		 int kappa, double mu) {
  double out = 0;
  if(kappa == 0) {
    out = single_gw0(d, sigsq, phi, mu);
  } else if(kappa == 1) {
    out = single_gw1(d, sigsq, phi, mu);
  } else 
  if(kappa == 2) {
    out = single_gw2(d, sigsq, phi, mu);
  } else {
    out = single_gw3(d, sigsq, phi, mu);
  }
  return out;
}

//' @title Generalized Wendland covariance function for a given distance matrix.
//'
//' @description Computing the Matern covariance function for a matrix of
//'   distances.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi \eqn{\phi} is the range of the covariance function.
//' @param kappa \eqn{\kappa \in \{0, \ldots, 3 \}}.
//' @param mu \eqn{\mu} controls the smoothness of the covariance function
//' @return The GW (isotropic) covariance function associated with the provided
//'   distances (\code{dists}) and the given set of parameters.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat gw_cov(const arma::mat& dists, double sigsq,
		 double phi, int kappa, double mu) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  if(kappa == 0) {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_gw0,
			     std::placeholders::_1,
			     sigsq, phi, mu));
  } else if(kappa == 1) {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_gw1,
			     std::placeholders::_1,
			     sigsq, phi, mu));
  } else if(kappa == 2) {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_gw2,
			     std::placeholders::_1,
			     sigsq, phi, mu));
  } else {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_gw3,
			     std::placeholders::_1,
			     sigsq, phi, mu));
  }  
  return out;
}

//' @title Mean of a Generalized Wendland covariance function (Internal use)
//'
//' @description This is an auxiliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polygon and speed-up the computations.
//'
//' @param dist a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} variance of the covariance function
//'   function.
//' @param phi the \eqn{\phi} is the range parameter of the covariance function.
//' @param kappa \eqn{\kappa \in \{0, \ldots, 3 \}}.
//' @param mu \eqn{\mu} controls the smoothness of the covariance function
//' @return The mean of \code{mat_cov(dist, sigsq, phi, kappa)}.
//' @keywords internal
// [[Rcpp::export]]
double aux_gw(arma::mat dist, double sigsq, double phi, int kappa, double mu) {
  return mean_mat(gw_cov(dist, sigsq, phi, kappa, mu));
}

//' @title Generalized Wendland covariance function for a polygons.
//'
//' @description Computing the Matern covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an integer representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} variance of the covariance function
//'   function.
//' @param phi the \eqn{\phi} is the range parameter of the covariance function.
//' @param kappa \eqn{\kappa \in \{0, \ldots, 3 \}}.
//' @param mu \eqn{\mu} controls the smoothness of the covariance function
//' @return The wendland-1 covariance matrix associated with a set of polygons.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_gw_cov(const List& cross_dists, int n,
		      int n2, double sigsq, double phi,
		      int kappa, double mu) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(),
		     cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_gw,
			       std::placeholders::_1,
			       sigsq, phi, kappa, mu));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(),
		   cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_gw, std::placeholders::_1,
			     sigsq, phi, kappa, mu));
  }
  return out;
}

//' @title Matern (tapered) covariance function (scalar - generic)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, adapted from \code{geoR} using Wendland-1 as a tapper.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param nu the smoothness parameter \eqn{\nu} from the Matern covariance
//'   function, controls the differentiability of the process.
//'
//' @param theta the \eqn{\theta} tapper range.
//'
//' @return a scalar representing the (tapered matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_matern}}, \code{\link{single_matern5}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_tapmat(double d, double sigsq,
		     double phi, double nu,
		     double theta) {
  double out = single_matern(d, sigsq, phi, nu) *
    single_gw1(d, 1.0, theta, 1.5);
  return out;
}

//' @title tapered Matern covariance function (scalar - nu = 1/2)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, with \eqn{\nu = 1/2}.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param theta \eqn{\theta} taper range.
//'
//' @return a scalar representing the (tapered matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_matern}}, \code{\link{single_matern5}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_tapmat1(double d, double sigsq,
		      double phi, double theta) {
  double out = single_exp(d, sigsq, phi) *
    single_gw1(d, 1.0, theta, 1.5);
  return out;
}

//' @title Tapered Matern covariance function (scalar - nu = 3/2)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, with \eqn{\nu = 3/2}.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param theta \eqn{\theta} taper range.
//'
//' @return a scalar representing the (tapered matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_matern}}, \code{\link{single_matern5}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_tapmat3(double d, double sigsq,
		      double phi, double theta) {
  double out = single_matern3(d, sigsq, phi) *
    single_gw1(d, 1.0, theta, 1.5);
  return out;
}

//' @title Tapered Matern covariance function for a given distance matrix.
//'
//' @description Computing the tapered Matern covariance function for a matrix
//'   of
//'   distances.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,'
//' @param nu smoothness parameter
//' @param theta \eqn{\theta} taper range.
//' 
//' @return The tapered matern covariance function (for a stationary and isotropic
//'   process) associated with the provided distances (\code{dists}) and the
//'   given set of parameters.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat tapmat_cov(const arma::mat& dists, double sigsq,
		     double phi, double nu, double theta) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  if(nu == .5) {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_tapmat1,
			     std::placeholders::_1,
			     sigsq, phi, theta));
  } else if(nu == 1.5) {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_tapmat3,
			     std::placeholders::_1,
			     sigsq, phi, theta));
  } else {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_tapmat,
			     std::placeholders::_1,
			     sigsq, phi, nu, theta));
  }
  return out;
}

//' @title Mean of a (Matern - Wendland-1 tapper) covariance function (Internal use)
//'
//' @description This is an auxiliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polygon and speed-up the computations.
//'
//' @param dist a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//' @param nu smoothness parameter
//' @param theta the \eqn{\theta} tapper range.
//' @return The mean of \code{mat_cov(dist, sigsq, phi, nu)}.
//' @keywords internal
// [[Rcpp::export]]
double aux_tapmat(arma::mat dist, double sigsq,
		  double phi, double nu,
		  double theta) {
  return mean_mat(tapmat_cov(dist, sigsq, phi, nu, theta));
}

//' @title Wendland-1 covariance function for a polygons.
//'
//' @description Computing the Matern covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an integer representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function
//' @param nu the smoothness parameter \eqn{\nu} parameter from the Matern
//'   covariance function
//' @param theta the taper distance.
//' @return The wendland-1 covariance matrix associated with a set of polygons.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_tapmat_cov(const List& cross_dists, int n,
			  int n2, double sigsq, double phi,
			  double nu, double theta) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(),
		     cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_tapmat,
			       std::placeholders::_1,
			       sigsq, phi,
			       nu, theta));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(),
		   cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_tapmat, std::placeholders::_1,
			     sigsq, phi,
			     nu, theta));
  }
  return out;
}
