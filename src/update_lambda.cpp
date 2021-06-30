// --------------------------- Description ---------------------------------- //
// update lambda 
// Input: ModelInput, ModelOutput 
// output: updated lambda 
// -------------------------------------------------------------------------- //


#define ARMA_64BIT_WORD 1


// #include "gperftools/profiler.h"
#include <RcppArmadillo.h>
#include "update_lambda.h"
#include "group_data.h"

using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
double update_lambda_cpp(Rcpp::S4 model_input, Rcpp::S4 model_output) {
    
    // ------------------------ start wrapper -------------------------------- //
    // ModelInput
    arma::mat data = Rcpp::as<arma::mat>(model_input.slot("m_data"));
    double hyper_c = Rcpp::as<double>(model_input.slot("m_hyper_c"));
    double hyper_d = Rcpp::as<double>(model_input.slot("m_hyper_d"));
    double hyper_a = Rcpp::as<double>(model_input.slot("m_hyper_a"));
    double hyper_b = Rcpp::as<double>(model_input.slot("m_hyper_b"));
    double total_mu = Rcpp::as<double>(model_input.slot("m_total_mu"));
    double n_groups = Rcpp::as<double>(model_input.slot("m_n_groups"));

    // ModelOutput
    double total_Ig1 = Rcpp::as<double>(model_output.slot("m_total_Ig1"));
    arma::mat Betas = Rcpp::as<arma::mat>(model_output.slot("m_Betas"));
    arma::vec Ig_vec = Rcpp::as<arma::vec>(model_output.slot("m_Ig_vec"));
    
    // helper variables
    double n_genes = data.n_rows;
    double n_exp = data.n_cols;
    
    // ------------------------ start update ----------------------------------- //
    
    // c*
    double c = hyper_c + 0.5 * (n_genes - total_Ig1 + total_Ig1 * n_groups);
    
    // d*
    arma::vec Ig_vec_group_data = Ig_vec + 1;
    Rcpp::List sbi = group_data2(Betas, Ig_vec_group_data); // new version
    // Rcpp::List sbi = substr(Betas, Ig_vec_group_data); // old version
    double sum_g0;
    double sum_g1;
    
    // all Betas of 1st row of m_Betas Matrix split up in two lists where list 1 is Ig == 0 & list 2 is Ig == 1 
    arma::mat y = sbi[0];
    
    if(y.size() == 0){ // if all n_genes are differentially expressed i.g. Ig=1
        sum_g0 = 0;
    
    } else {
        arma::mat beta_g0;
        beta_g0 = trans(y.row(0)); // can be declared as vector
        arma::rowvec mu_vec_r(1);
        mu_vec_r.fill(total_mu);
        beta_g0.each_row() -= mu_vec_r;
        sum_g0 = arma::as_scalar(trans(beta_g0) * beta_g0);
    }
    
    
    if(y.size() == (n_genes *2)){
        sum_g1 = 0;
    } else {
        arma::mat beta_g1 = sbi[1];
        arma::rowvec mu_vec_r(beta_g1.n_cols);
        mu_vec_r.fill(total_mu);
        beta_g1.each_row() -= mu_vec_r;
        arma::vec e = vectorise(beta_g1);
        sum_g1 = arma::as_scalar(trans(e)*e);
    }
    
    double d = 1.0 / (hyper_d + 0.5* (sum_g0 + sum_g1));
    double lambda_updated = R::rgamma(c, d);
    return lambda_updated;
}

/*** R

*/
