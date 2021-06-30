#ifndef UPDATE_BETAI_CPP_H
#define UPDATE_BETAI_CPP_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
Rcpp::S4 update_betai_cpp(Rcpp::S4 model_input, Rcpp::S4 model_output);
// void update_betai_cpp2(const Rcpp::S4 &model_input, Rcpp::S4 &model_output);

#endif


