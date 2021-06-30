#ifndef UPDATE_MHNUPHI_CPP_H
#define UPDATE_MHNUPHI_CPP_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
Rcpp::S4 update_mhnuphi_cpp(Rcpp::S4 &model_input, Rcpp::S4 &model_output);

#endif



