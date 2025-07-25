#include <RcppEigen.h>
#include <cmath> // For std::hypot

// [[Rcpp::depends(RcppEigen)]]

//' @title Creating a symmetric distance matrix (Eigen version)
//' @description Computes the Euclidean distance between all pairs of rows in a
//'   matrix.
//' @param my_mat A matrix where each row is a 2D coordinate.
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd distmat(const Eigen::MatrixXd& my_mat) {
  const int n = my_mat.rows();
  Eigen::MatrixXd out = Eigen::MatrixXd::Zero(n, n);
  for (int j = 0; j < n; ++j) {
    // Start from i = j + 1 to only compute the upper triangle
    for (int i = j + 1; i < n; ++i) {
      const double dist = std::hypot(my_mat(i, 0) - my_mat(j, 0),
                                     my_mat(i, 1) - my_mat(j, 1));
      out(i, j) = dist;
      out(j, i) = dist; // Enforce symmetry
    }
  }
  return out;
}

//' @title Pairwise distances between two matrices (Eigen version)
//' @description Computes Euclidean distance between rows of m1 and rows of m2.
//' @param m1 A matrix where each row is a 2D coordinate.
//' @param m2 A matrix where each row is a 2D coordinate.
//' @keywords internal
// [[Rcpp::export]]
Eigen::MatrixXd crossdist(const Eigen::MatrixXd& m1,
                          const Eigen::MatrixXd& m2) {
  const int nrow1 = m1.rows();
  const int nrow2 = m2.rows();
  Eigen::MatrixXd out(nrow1, nrow2);
  for (int j = 0; j < nrow2; ++j) {
    for (int i = 0; i < nrow1; ++i) {
      out(i, j) = std::hypot(m1(i, 0) - m2(j, 0),
                             m1(i, 1) - m2(j, 1));
    }
  }
  return out;
}

//' @title Pairwise distances for a list of matrices (Internal use)
//' @param mat_list internal use
//' @param mat_list1 internal use
//' @param mat_list2 internal use
//' @param return_single internal use
//' @param pred_mat internal use
//' @param x_to_list internal use
//' @param by internal use
//' @param y_grid internal use
//' @param x_grid internal use
//' @param tr_vec index for distance truncation
//' @param tr_inp truncation input
//' @name aux_mat
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List single_dists(const Rcpp::List& mat_list) {
  const int n_out = mat_list.size();
  Rcpp::List out(n_out * (n_out + 1) / 2);
  int k = 0;
  for (int j = 0; j < n_out; ++j) {
    for (int i = j; i < n_out; ++i) {
      out[k++] = crossdist(Rcpp::as<Eigen::MatrixXd>(mat_list[i]),
                           Rcpp::as<Eigen::MatrixXd>(mat_list[j]));
    }
  }
  return out;
}

//' @title Cross-distances for two lists of matrices (Eigen version)
//' @rdname aux_mat
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List mult_dists(const Rcpp::List& mat_list1,
                      const Rcpp::List& mat_list2,
                      const bool return_single) {
  const int n1 = mat_list1.size();
  const int n2 = mat_list2.size();
  Rcpp::List out_cross(n1 * n2);
  int k = 0;
  for (int j = 0; j < n2; ++j) {
    for (int i = 0; i < n1; ++i) {
      out_cross[k++] = crossdist(Rcpp::as<Eigen::MatrixXd>(mat_list1[i]),
                                 Rcpp::as<Eigen::MatrixXd>(mat_list2[j]));
    }
  }
  if (return_single) {
    // Avoid re-computing single_dists twice
    Rcpp::List dists2 = single_dists(mat_list2);
    return Rcpp::List::create(
                              Rcpp::Named("dists_1") = single_dists(mat_list1),
                              Rcpp::Named("dists_2") = dists2,
                              Rcpp::Named("cross")   = out_cross
                              );
  }
    
  return out_cross;
}

//' @title Prediction cross-distances (Eigen version)
//' @rdname aux_mat
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List pred_cdist(const Rcpp::List& mat_list,
                      const Eigen::MatrixXd& pred_mat) {
  const int n1 = mat_list.size();
  const int n2 = pred_mat.rows();
  Rcpp::List out_cross(n1 * n2);
  int k = 0;
  for (int j = 0; j < n2; ++j) {
    for (int i = 0; i < n1; ++i) {
      // .row(j) is efficient in Eigen, it returns a view not a copy
      out_cross[k++] = crossdist(Rcpp::as<Eigen::MatrixXd>(mat_list[i]),
                                 pred_mat.row(j));
    }
  }
  return out_cross;
}
