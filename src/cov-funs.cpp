#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <numeric>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//' @title Mean of a matrix (Internal use)
//'
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon.
//'
//' @param mat_aux a numeric matrix.
//' @return The mean of \code{mat_aux}.
//' @keywords internal
double mean_mat(arma::mat mat_aux) {
  return std::accumulate(mat_aux.begin(), mat_aux.end(), 0.0) /
    (mat_aux.n_rows * mat_aux.n_cols);
}

//' @title Matern covariance function (scalar - generic)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, addapted from \code{geoR}.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param kappa the \eqn{\kappa} parameter from the Matern covariance function,
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
		     double phi, double kappa) {
  double out = sigsq;
  if(d > 0) {
    out *=
      (pow(2, 1 - kappa) / ::tgamma(kappa)) *
      pow(d / phi, kappa) *
      Rf_bessel_k(d / phi, kappa, 1);
  } 
  return out;
}

//' @title Matern covariance function (scalar - kappa = 3/2)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, with \eqn{\kappa = 3/2}.
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

//' @title Matern covariance function (scalar - kappa = 5/2)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, with \eqn{\kappa = 5/2}.
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
//' @param kappa the \eqn{\kappa} parameter from the Matern covariance function,
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
		  double phi, double kappa) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  if(kappa == 0.5) {
    std::transform(dists.begin(), dists.end(),
		   out.begin(), std::bind(&single_exp, std::placeholders::_1,
					  sigsq, phi));
  } else if(kappa == 1.5) {
    std::transform(dists.begin(), dists.end(),
		   out.begin(), std::bind(&single_matern3, std::placeholders::_1,
					  sigsq, phi));
  } else if(kappa == 2.5) {
    std::transform(dists.begin(), dists.end(),
		   out.begin(), std::bind(&single_matern5, std::placeholders::_1,
					  sigsq, phi));
  } else {
    std::transform(dists.begin(), dists.end(),
		   out.begin(), std::bind(&single_matern, std::placeholders::_1,
					  sigsq, phi, kappa));
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
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon and speed-up the computations.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param kappa the \eqn{\kappa} parameter from the Matern covariance function,
//'   controls the differentiability of the process.
//' 
//' @return The mean of \code{mat_cov(dist, sigsq, phi, kappa)}.
//' @keywords internal
// [[Rcpp::export]]
double aux_matern(arma::mat dist, double sigsq,
		  double phi, double kappa) {
  return mean_mat(mat_cov(dist, sigsq, phi, kappa));
}

//' @title Matern covariance function for a polygons.
//'
//' @description Computing the Matern covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an ingeger representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param kappa the \eqn{\kappa} parameter from the Matern covariance function,
//'   controls the differentiability of the process. Note that, if we set
//'   \code{kappa = .5}, then the calculations are based on the exponential
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
		       double kappa) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(), cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_matern, std::placeholders::_1,
			       sigsq, phi, kappa));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(), cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_matern, std::placeholders::_1,
			     sigsq, phi, kappa));
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
//' @param kappa the \eqn{\kappa \in (0, 2]} parameter representing the "power"
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
double single_pexp(double d, double sigsq, double phi, double kappa) {
  double out = sigsq;
  if(d > 0) {
    out *= exp( - pow( (d / phi) , kappa));
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
//' @param kappa the \eqn{\kappa \in (0, 2]} parameter representing the "power"
//' 
//' @return The powered exponential covariance function (for a stationary and
//'   isotropic process) associated with the provided distances (\code{dists})
//'   and the given set of parameters.
//'
//' @seealso \code{\link{single_exp}}, \code{\link{single_matern}}
//' @keywords internal
// [[Rcpp::export]]
arma::mat pexp_cov(const arma::mat& dists, double sigsq,
		   double phi, double kappa) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  std::transform(dists.begin(), dists.end(),
		 out.begin(), std::bind(&single_pexp, std::placeholders::_1,
					sigsq, phi, kappa));
  return out;
}

//' @title Mean of a (Powered Exponential) covariance function (Internal use)
//'
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon and speed-up the computations.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Powered Exponential
//'   covariance function.
//' @param phi the \eqn{\phi} parameter from the Powered Exponential covariance
//'   function, controls the range of the spatial dependence.
//' @param kappa the \eqn{\kappa} parameter from the Powered Exponential
//'   covariance function,
//'   controls the differentiability of the process.
//' 
//' @return The mean of \code{pexp_cov(dist, sigsq, phi, kappa)}.
//' @keywords internal
// [[Rcpp::export]]
double aux_pexp(arma::mat dist, double sigsq,
	       double phi, double kappa) {
  return mean_mat(pexp_cov(dist, sigsq, phi, kappa));
}

//' @title Powered Exponential covariance function for a polygons.
//'
//' @description Computing the Powered Exponential covariance function between
//'   polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an ingeger representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Powered Exponential
//'   covariance function.
//' @param phi the \eqn{\phi} parameter from the Powered Exponential covariance
//'   function, controls the range of the spatial dependence.
//' @param kappa the \eqn{\kappa \in (0, 2]} parameter representing the "power"
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
			double kappa) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(), cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_pexp, std::placeholders::_1,
			       sigsq, phi, kappa));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(), cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_pexp, std::placeholders::_1,
			     sigsq, phi, kappa));
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
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon and speed-up the computations.
//'
//' @param dists a numeric matrix representing the distance between spatial
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
//' @param n an ingeger representing number of polygons (note that, this is
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
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon and speed-up the computations.
//'
//' @param dists a numeric matrix representing the distance between spatial
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
//' @param n an ingeger representing number of polygons (note that, this is
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
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon and speed-up the computations.
//'
//' @param dists a numeric matrix representing the distance between spatial
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
//' @param n an ingeger representing number of polygons (note that, this is
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

//' @title Matern Wendland-1 covariance unction(scalar - generic)
//'
//' @description adapted from Furrer et al. 2006.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' 
//' @return a scalar representing the (wendland) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_w1(double d, double sigsq, double phi) {
  double out = sigsq, aux = d / phi;
  if(aux < 1) {
    out *= (1 - aux) * (1 - aux) * (1 + .5 * aux);
  } else {
    out *= 0;
  }
  return out;
}

//' @title Wendland-1 covariance function for a given distance matrix.
//'
//' @description Computing the Matern covariance function for a matrix of
//'   distances.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,' 
//' @return The matern covariance function (for a stationary and isotropic
//'   process) associated with the provided distances (\code{dists}) and the
//'   given set of parameters.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat w1_cov(const arma::mat& dists, double sigsq,
		 double phi) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  std::transform(dists.begin(), dists.end(), out.begin(),
		 std::bind(&single_w1,
			   std::placeholders::_1,
			   sigsq, phi));
  return out;
}

//' @title Mean of a (Wendland-1) covariance function (Internal use)
//'
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon and speed-up the computations.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//' @return The mean of \code{mat_cov(dist, sigsq, phi, kappa)}.
//' @keywords internal
// [[Rcpp::export]]
double aux_w1(arma::mat dist, double sigsq, double phi) {
  return mean_mat(w1_cov(dist, sigsq, phi));
}

//' @title Wendland-1 covariance function for a polygons.
//'
//' @description Computing the Matern covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an ingeger representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,' 
//' @return The wendland-1 covariance matrix associated with a set of polygons.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_w1_cov(const List& cross_dists, int n,
		      int n2, double sigsq, double phi) {
  arma::mat out(n, n2, arma::fill::zeros);
  if(n == n2) {
      arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
      arma::vec aux(cross_dists.size(), arma::fill::zeros);
      
      std::transform(cross_dists.begin(),
		     cross_dists.end(),
		     aux.begin(),
		     std::bind(&aux_w1,
			       std::placeholders::_1,
			       sigsq, phi));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(),
		   cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_w1, std::placeholders::_1,
			     sigsq, phi));
  }
  return out;
}

//' @title Matern (tappered) covariance function (scalar - generic)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, addapted from \code{geoR} using Wendland-1 as a tapper.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param kappa the \eqn{\kappa} parameter from the Matern covariance function,
//'   controls the differentiability of the process.
//'
//' @param theta the \eqn{\theta} tapper range.
//'
//' @return a scalar representing the (tappered matern) covariance between two
//'   observations \code{d} apart of each other.
//' 
//' @seealso \code{\link{single_matern}}, \code{\link{single_matern5}}
//'   \code{\link{single_exp}}, \code{\link{mat_cov}}
//' 
//' @keywords internal
// [[Rcpp::export]]
double single_tapmat(double d, double sigsq,
		     double phi, double kappa,
		     double theta) {
  double out = single_matern(d, sigsq, phi, kappa) *
    single_w1(d, 1.0, theta);
  return out;
}

//' @title Tappered Matern covariance function (scalar - kappa = 1/2)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, with \eqn{\kappa = 3/2}.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param theta \eqn{\theta} taper range.
//'
//' @return a scalar representing the (tappered matern) covariance between two
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
    single_w1(d, 1.0, theta);
  return out;
}

//' @title Tappered Matern covariance function (scalar - kappa = 3/2)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, with \eqn{\kappa = 3/2}.
//'
//' @param d a scalar representing the distance on which it is desired to
//'   evaluate the covariance function.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//'   controls the range of the spatial dependence.
//' @param theta \eqn{\theta} taper range.
//'
//' @return a scalar representing the (tappered matern) covariance between two
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
    single_w1(d, 1.0, theta);
  return out;
}


//' @title Tappered Matern  covariance function for a given distance matrix.
//'
//' @description Computing the tappered Matern covariance function for a matrix
//'   of
//'   distances.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,'
//' @param kappa smoothness parameter
//' @param theta \eqn{\theta} taper range.
//' 
//' @return The tappered matern covariance function (for a stationary and isotropic
//'   process) associated with the provided distances (\code{dists}) and the
//'   given set of parameters.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat tapmat_cov(const arma::mat& dists, double sigsq,
		     double phi, double kappa, double theta) {
  int nr = dists.n_rows, nc = dists.n_cols;
  arma::mat out(nr, nc, arma::fill::zeros);
  if(kappa == .5) {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_tapmat1,
			     std::placeholders::_1,
			     sigsq, phi, theta));
  } else if(kappa == 1.5) {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_tapmat3,
			     std::placeholders::_1,
			     sigsq, phi, theta));
  } else {
    std::transform(dists.begin(), dists.end(), out.begin(),
		   std::bind(&single_tapmat,
			     std::placeholders::_1,
			     sigsq, phi, kappa, theta));
  }
  return out;
}

//' @title Mean of a (Matern - Wendland-1 tapper) covariance function (Internal use)
//'
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon and speed-up the computations.
//'
//' @param dists a numeric matrix representing the distance between spatial
//'   entities.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,
//' @param kappa
//' @param theta
//' @return The mean of \code{mat_cov(dist, sigsq, phi, kappa)}.
//' @keywords internal
// [[Rcpp::export]]
double aux_tapmat(arma::mat dist, double sigsq,
		  double phi, double kappa,
		  double theta) {
  return mean_mat(tapmat_cov(dist, sigsq, phi, kappa, theta));
}

//' @title Wendland-1 covariance function for a polygons.
//'
//' @description Computing the Matern covariance function between polygons.
//'
//' @param cross_dists a \code{list} such that each position contains the cross
//'   distances between points within different polygons.
//' @param n an ingeger representing number of polygons (note that, this is
//'   different than the size of the list \code{cross_dists}
//' @param n2 usually, equal to \code{n}, except when the function is being used
//'   to calculate the "cross" covariance between two different partitions of
//'   the same space.
//' @param sigsq the \eqn{\sigma^2} parameter from the Matern covariance
//'   function.
//' @param phi the \eqn{\phi} parameter from the Matern covariance function,' 
//' @return The wendland-1 covariance matrix associated with a set of polygons.
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat comp_tapmat_cov(const List& cross_dists, int n,
			  int n2, double sigsq, double phi,
			  double kappa, double theta) {
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
			       kappa, theta));
      
      out.elem(lw_idx) = aux;
      out = arma::symmatl(out);
  } else {
    std::transform(cross_dists.begin(),
		   cross_dists.end(),
		   out.begin(),
		   std::bind(&aux_tapmat, std::placeholders::_1,
			     sigsq, phi,
			     kappa, theta));
  }
  return out;
}
