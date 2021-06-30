// ---------------------------- description --------------------------------- //
// perform within dimension (wd) update
// Input: required variables are passed by refernce 

// Output: List of variables updated for the current gene (I, beta, mu, tau) 
// -------------------------------------------------------------------------- //

#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include "update_wdbetai.h"
#include <cmath>  

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List update_wdbetai_cpp(arma::mat data, arma::mat design, double total_mu, int n_groups, arma::mat Betas, arma::mat phi_scaling, arma::vec Ig_vec, double taue, int current_gene, double lambda){
    
    int n_exp = data.n_cols;

    arma::mat beta = arma::zeros(n_groups, 1);
    arma::mat mu = arma::zeros(n_groups, 1);
    arma::mat tau = arma::zeros(n_groups, 1);
    
    if( Ig_vec(current_gene) == 0.0 ){
        
        double sum_y = arma::as_scalar(data.row(current_gene) * trans(phi_scaling.row(current_gene)));  
        double sum_phi = arma::as_scalar(phi_scaling.row(current_gene) * arma::ones(n_exp,1)); // matlab: sphi 
        double uni_tau = taue * sum_phi + lambda; // update of lambda (lambda*) 
        double uni_mu = (taue * sum_y + lambda * total_mu) / uni_tau ; // uni_mu (update of mu to update wd betag Ig==0)
        // draw b0 from N(mu, tau^-1)
        double z = R::rnorm(0,1);
        double b0 = uni_mu + z / std::sqrt(uni_tau); 
        
        beta.fill(b0);
        mu.fill(uni_mu);
        tau.fill(uni_tau);

    } else {
        
        arma::mat xdphixt = design * diagmat(phi_scaling.row(current_gene)) * trans(design);
        arma::mat Lambda = (taue * xdphixt + lambda * arma::eye(n_groups, n_groups));
        arma::vec Lambda_extr;
        Lambda_extr = Lambda.elem( find(Lambda != 0));
        tau = arma::conv_to<arma::mat>::from(Lambda_extr); // store Lambda_extr as matrix, # Lambda in mathematischer Formel heißt im code tau 
        arma::mat sumxphiy = design * diagmat(phi_scaling.row(current_gene)) * trans(data.row(current_gene));
        mu = diagmat(1/tau) * (lambda * total_mu * arma::ones(n_groups, 1) + taue * sumxphiy );
        arma::mat beta_rnorm(n_groups, 1); // ToDo: use arma function
        beta_rnorm = Rcpp::rnorm(n_groups); // ToDo: use arma function
        beta = arma::diagmat(1/arma::sqrt(tau)) * beta_rnorm + mu;
    }
    
    double I = Ig_vec(current_gene);
    Rcpp::List wd_output = Rcpp::List::create(Rcpp::Named("I")=I  , Rcpp::_["beta"] = beta , Rcpp::_["mu"] = mu, Rcpp::_["tau"] = tau );
    
    return wd_output;
}


/*** R
*/
