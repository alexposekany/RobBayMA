#ifndef MATRIX_RGAMMA_H
#define MATRIX_RGAMMA_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
arma::mat matrix_rgamma(int n_genes, int n_rows, double shape, arma::mat scaling_mat);

#endif


