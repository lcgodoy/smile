#include <RcppArmadillo.h>
#include <Rcpp.h>

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

//' @title Internal use only
//' @param mat_list internal use
//' @param mat_list1 internal use
//' @param mat_list2 internal use
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

  for(int i = 0; i < n1; i++) {
    for(int j = 0; j < n2; j++) {
      out_cross[k] = crossdist(mat_list[i], pred_mat.row(j));
      k += 1;
    }
  }
  return out_cross;
}

