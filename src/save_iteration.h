#ifndef SAVE_ITERATION_CPP_H
#define SAVE_ITERATION_CPP_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
Rcpp::S4 save_iteration_cpp(Rcpp::S4 &mcmc_output, Rcpp::S4 &model_output, arma::mat su);
    

#endif



