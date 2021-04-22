#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//' @name aux_mat
// [[Rcpp::export]]
double eucl_aux(double x, double y) {
  return sqrt( x*x + y*y );
}

//' @title Creatin a distance matrix
//' 
//' @description Internal use. For now it only supports euclidean distance. (May
//'   import parallelDist in the future).
//'
//' @param my_mat a matrix representing a grid of points.
//' 
// [[Rcpp::export]]
arma::mat distmat(const arma::mat& my_mat) {
  int n = my_mat.n_rows, k = 0;
  arma::vec aux( n * (n + 1) * .5,
		 arma::fill::zeros);

  arma::mat out(n, n, arma::fill::zeros);
  
  arma::uvec lw_idx = arma::trimatl_ind( arma::size(out) );
  
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      if(i == j) {
	aux(k) = 0;
      } else {
        aux(k) = eucl_aux( my_mat(i, 0) - my_mat(j, 0),
			   my_mat(i, 1) - my_mat(j, 1)  );
      }
      k += 1;
    }
  }
 
  out.elem(lw_idx) = aux;
  out = arma::symmatl(out);
  
  return out;
}

//' @title Pairwise distances between matrices
//' 
//' @description Internal use.
//'
//' @param m1 a matrix representing a grid of points within a polygon.
//' @param m2 a matrix representing a grid of points within a polygon.
//' 
// [[Rcpp::export]]
arma::mat crossdist(const arma::mat& m1, const arma::mat& m2) {

  int nrow1 = m1.n_rows, nrow2 = m2.n_rows;
 
  arma::mat out(nrow1, nrow2, arma::fill::zeros);

  for (int i = 0; i < nrow1; i++) {
    for (int j = 0; j < nrow2; j++) {
      out(i, j) = eucl_aux( m1(i, 0) - m2(j, 0),
			    m1(i, 1) - m2(j, 1)  );
    }
  }
  
  return out;
}

//' @title Internal use only
//' @param mat_list internal use
//' @param mat_list1 internal use
//' @param mat_list2 internal use
//' @param x internal use
//' @param y internal use
//' @param return_single internal use
//' @param pred_mat internal use
//' @param x_to_list internal use
//' @param by internal use
//' @param y_grid internal use
//' @param x_grid internal use
//' @name aux_mat
// [[Rcpp::export]]
List single_dists(const List& mat_list) {
  int n_out = mat_list.size(), k = 0;
  List out( n_out * (n_out + 1) / 2 );
  for(int i = 0; i < n_out; i++) {
    for(int j = i; j < n_out; j++) {
      out[k] = crossdist(mat_list[i], mat_list[j]);
      k += 1;
    }
  }
  
  return out;
}

//' @name aux_mat
// [[Rcpp::export]]
List mult_dists(const List& mat_list1, const List& mat_list2,
		const bool& return_single) {
  int n1 = mat_list1.size(), n2 = mat_list2.size(), k = 0;
  List out_cross( n1 * n2 );
  for(int i = 0; i < n1; i++) {
    for(int j = 0; j < n2; j++) {
      out_cross[k] = crossdist(mat_list1[i], mat_list2[j]);
      k += 1;
    }
  }
  if(return_single) {
    out_cross =  List::create(Rcpp::Named("dists_1") = single_dists(mat_list2),
			      Rcpp::Named("dists_2") = single_dists(mat_list2),
			      Rcpp::Named("cross")   = out_cross);
  }

  return out_cross;
}

//' @name aux_mat
// [[Rcpp::export]]
List pred_cdist(const List& mat_list, const arma::mat& pred_mat) {

  int n1 = mat_list.size(), n2 = pred_mat.n_rows, k = 0;
  List out_cross( n1 * n2 );

  for(int i = 0; i < n2; i++) {
    for(int j = 0; j < n1; j++) {
      out_cross[k] = crossdist(mat_list[j], pred_mat.row(i));
      k += 1;
    }
  }
  
  return out_cross;
}

