#ifndef UPDATE_BETAWD_H
#define UPDATE_BETAWD_H

#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
Rcpp::List update_wdbetai_cpp(arma::mat data, arma::mat design, double total_mu, int n_groups, arma::mat Betas, arma::mat phi_scaling, arma::vec Ig_vec, double taue, int current_gene, double lambda);

#endif
