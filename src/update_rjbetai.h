#ifndef UPDATE_BETARJ_H
#define UPDATE_BETARJ_H

#include <RcppArmadillo.h>
#include <cmath> 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
Rcpp::List update_rjbetai_cpp(arma::mat data, arma::mat design, double total_mu, int n_groups, arma::mat Betas, arma::mat phi_scaling, arma::vec Ig_vec, double taue, double p, arma::mat gene_mu_grouped, arma::mat m_tau, int current_gene, double lambda);

#endif