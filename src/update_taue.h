#ifndef UPDATE_TAUE_CPP_H
#define UPDATE_TAUE_CPP_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
double update_taue_cpp(Rcpp::S4 &model_input, Rcpp::S4 &model_output);

#endif



