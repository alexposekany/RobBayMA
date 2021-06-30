// ---------------------------- description --------------------------------- //
// Updates degrees of freedom (nu) of student t distribution in a reversible
// jump step

// Input: required variables are passed by refernce 
// Output: List of variables updated for the current gene (I, beta, mu, tau) 
// -------------------------------------------------------------------------- //
    


#define ARMA_64BIT_WORD 1
    
#include <RcppArmadillo.h>
#include "update_rjbetai.h"
#include <cmath>  
    
    
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::List update_rjbetai_cpp(arma::mat data, arma::mat design, double total_mu, int n_groups, arma::mat Betas, arma::mat phi_scaling, arma::vec Ig_vec, double taue, double p, arma::mat gene_mu_grouped, arma::mat m_tau, int current_gene, double lambda){
        
        int n_exp = data.n_cols;
        
        // ------------- start update ---------------------------------------- //
        arma::mat beta = arma::zeros(n_groups, 1);
        arma::mat mu = arma::zeros(n_groups, 1);
        arma::mat tau = arma::zeros(n_groups, 1);
        double I;
        
        // Multivariate proposal variables
        arma::mat xdphixt = design * diagmat(phi_scaling.row(current_gene)) * trans(design); 
        arma::mat Lambda = (taue * xdphixt) + (lambda *arma::eye(n_groups, n_groups));
        arma::vec Lambda_extr;
        double inf = std::numeric_limits<double>::infinity(); 
        if(taue == inf){
            Lambda_extr = Lambda.elem( find(Lambda == inf));
        } else {
            Lambda_extr = Lambda.elem( find(Lambda != 0)); // extract non zero-elements from Lambda 
        }
        arma::mat multi_tau = arma::conv_to<arma::mat>::from(Lambda_extr); // store Lambda_extr as matrix, # Lambda in mathematischer Formel heißt im code tau 
        arma::mat sumxphiy = design * diagmat(phi_scaling.row(current_gene)) * trans(data.row(current_gene));
        arma::mat multi_mu = diagmat(1/multi_tau) * (lambda * total_mu * arma::ones(n_groups, 1) + taue * sumxphiy);
        // Univariate proposal variables
        double sum_y = arma::as_scalar( data.row(current_gene)*trans(phi_scaling.row(current_gene)) ); //matlab: sy,  var1 ist ein vector, var2 eine matrix. Beim Kreuzprodukt von Vektor und Matrix macht R transponieren selber. Ergebniss ist eine Matrix 
        double sum_phi = arma::as_scalar( phi_scaling.row(current_gene)*arma::ones(n_exp,1) ); // matlab: sphi 
        double uni_tau = taue * sum_phi + lambda; // update of lambda (lambda*) 
        double uni_mu = (taue * sum_y + lambda * total_mu) / uni_tau ; // uni_mu (update of mu to update wd betag Ig==0)
        
        // Acceptance probability
        double diffp = std::log(p) - std::log(1 - p); // call to log is ambiguous
        double sumlogtau = arma::as_scalar(arma::ones(1, n_groups) * arma::log(multi_tau));
        double diffdet = std::log(uni_tau) - sumlogtau;
        double diffpr = (n_groups - 1) * std::log(lambda); 
        double sumsqpr = (n_groups - 1)* lambda * std::pow(total_mu, 2);
        double sumsqprop = arma::as_scalar( (uni_tau * std::pow(uni_mu, 2)) - (trans(multi_mu) * Lambda * multi_mu));
        double logaccp = diffp + 0.5 * (diffpr + diffdet - (sumsqpr + sumsqprop));
       
        if(Ig_vec(current_gene) == 0.0){
            
        // Proposal: univariate beta --> multivariate beta
        double u = R::runif(0,1);
        if(u<std::exp(logaccp)){
            I = 1;
            arma::mat beta_rnorm(n_groups, 1); // To Do: use arma function
            beta_rnorm = Rcpp::rnorm(n_groups); // To Do: use arma function
            beta = diagmat(1/arma::sqrt(multi_tau)) * beta_rnorm + multi_mu;
            mu = multi_mu;
            tau = multi_tau;
            
            } else { // reject proposal and keep univariate Betas
                I = 0;
                beta = Betas.col(current_gene);
                mu = gene_mu_grouped.col(current_gene);
                tau = m_tau.col(current_gene);
                }
    } else { // Proposal: univariate -> multivariate 
        double u = R::runif(0,1);
        if(u < std::exp(-logaccp)){ 
            I = 0;
            // draw new value 
            // draw b0 from N(mu, tau^-1)
            double z = R::rnorm(0,1);
            double b0 = uni_mu + z / std::sqrt(uni_tau); // To Do: Call to sqrt is ambiguous
            beta.fill(b0);
            mu.fill(uni_mu);
            tau.fill(uni_tau);
        } else { // reject proposal and keep old, multivariate Betas
            I = 1;
            beta = Betas.col(current_gene);
            mu = gene_mu_grouped.col(current_gene);
            tau = m_tau.col(current_gene);
        }
    }
        
    
    Rcpp::List rj_output = Rcpp::List::create(Rcpp::Named("I")=I  , Rcpp::_["beta"] = beta , Rcpp::_["mu"] = mu, Rcpp::_["tau"] = tau );
    return rj_output;
}


/*** R

*/
