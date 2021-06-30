#ifndef UPDATE_LAMBDA_H
#define UPDATE_LAMBDA_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
double update_lambda_cpp(Rcpp::S4 model_input, Rcpp::S4 model_output);

#endif