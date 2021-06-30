#ifndef GROUP_DATA_H
#define GROUP_DATA_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
Rcpp::List group_data1(arma::mat& data, arma::vec& grouping_vec);
Rcpp::List group_data2(arma::mat& data, arma::vec& grouping_vec);
Rcpp::List substr(arma::mat& data, arma::vec& grouping_vec);
Rcpp::NumericMatrix get_design(Rcpp::NumericVector grouping_vec, int n_groups);

#endif


