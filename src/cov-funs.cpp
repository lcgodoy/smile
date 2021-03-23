#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <numeric>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//' @title Pairwise distances between matrices
//' 
//' @description Internal use.
//'
//' @param m1 a matrix representing a grid of points within a polygon.
//' @param m2 a matrix representing a grid of points within a polygon.
//' 
// [[Rcpp::export]]
arma::mat crossdist(const arma::mat& m1, const arma::mat& m2) {

  int nrow1 = m1.n_rows, nrow2 = m2.n_rows, ncol = m1.n_cols;
 
  if (m1.n_cols != m2.n_cols) {
    throw std::runtime_error("matrices with different number of columns");
  }
 
  arma::mat out(nrow1, nrow2, arma::fill::zeros);

  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r2, c12), 2);
      }
      out(r1, r2) = sqrt(total);
    }
  }
  return out;
}

//' @title Mean of a matrix (Internal use)
//'
//' @description This is an auxilliary function for internal use. It helps to
//'   numerically integrate a covariance function evaluated at a grid of points
//'   within a polyigon.
//'
//' @param mat_aux a numeric matrix.
//' @return The mean of \code{mat_aux}.
double mean_mat(arma::mat mat_aux) {
  return std::accumulate(mat_aux.begin(), mat_aux.end(), 0.0) /
    (mat_aux.n_rows * mat_aux.n_cols);
}

//' @title Matern covariance function (scalar - generic)
//'
//' @description Computing the Matern covariance function for a single distance
//'   measure, addapted from \code{\link[geoR]{matern}}.
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
// [[Rcpp::export]]
double single_matern(double d, double sigsq,
		     double phi, double kappa) {
  double out = sigsq;
  if(d > 0) {
    out =
      sigsq * (pow(2, 1 - kappa) / ::tgamma(kappa)) *
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
// [[Rcpp::export]]
double single_matern3(double d, double sigsq,
		      double phi) {
  double out = sigsq;
  if(d > 0) {
    out =
      sigsq * (1 + (sqrt(3) * (d/phi))) * exp( - sqrt(3) * (d/phi));
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
// [[Rcpp::export]]
double single_matern5(double d, double sigsq,
		      double phi) {
  double out = sigsq;
  if(d > 0) {
    out =
      sigsq * (1 + ( sqrt(5) * ( d/phi ) ) + ( (5 / 3) * pow( ( d / phi ), 2 )) ) *
      exp( - sqrt(5) * (d/phi));
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
// [[Rcpp::export]]
double single_exp(double d, double sigsq, double phi) {
  double out = sigsq;
  if(d > 0) {
    out = sigsq * exp(- d / phi);
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
// [[Rcpp::export]]
double single_pexp(double d, double sigsq, double phi, double kappa) {
  double out = sigsq;
  if(d > 0) {
    out = sigsq * exp( - pow( (d / phi) , kappa));
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
//' 
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
//' @export
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
// [[Rcpp::export]]
double single_gauss(double d, double sigsq, double phi) {
  double out = sigsq;
  if(d > 0) {
    out = sigsq * exp( - .5 * pow( (d / phi) , 2));
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
// [[Rcpp::export]]
double single_spher(double d, double sigsq, double phi) {
  double out = sigsq;
  if(d < phi) {
    out = 0;
  } else if(d > 0) {
    out = sigsq *
      (1 - ( (2 / M_PI) * sqrt(1 - pow( (d / phi) , 2) ) ) + asin( (d / phi) )  );
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

