// --------------------------- Description ---------------------------------- //
// Updates degrees of freedom (nu) of student t distribution in a reversible
// jump step
// 
// Input: ModelInput, ModelOutput 
// Output: ModelOutput with updated variables nu and phi 
// -------------------------------------------------------------------------- //


#define ARMA_64BIT_WORD 1

// #include <RcppArmadilloExtensions/sample.h>
#include <cmath> 
#include "matrix_rgamma.h"
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::S4 update_mhnuphi_cpp(Rcpp::S4 &model_input, Rcpp::S4 &model_output) {
    
    // ------------------ extract data and helper variables ----------------- //
    // ModelInput
    arma::mat data = Rcpp::as<arma::mat>(model_input.slot("m_data"));
    arma::mat design = Rcpp::as<arma::mat>(model_input.slot("m_design_mat"));
    double numax = Rcpp::as<double>(model_input.slot("m_numax")); 
    int n_groups = Rcpp::as<int>(model_input.slot("m_n_groups"));
    
    // ModelOutput
    arma::mat Betas = Rcpp::as<arma::mat>(model_output.slot("m_Betas"));
    arma::vec Ig_vec = Rcpp::as<arma::vec>(model_output.slot("m_Ig_vec"));
    double nu_old = Rcpp::as<double>(model_output.slot("m_nu"));
    double taue = Rcpp::as<double>(model_output.slot("m_taue"));
    
    // helper varaibles
    int n_genes = data.n_rows;
    int n_exp = data.n_cols;
    
    // ------------------------ start update ----------------------------------- //
    
    // nuf: factor for prop(nu_o|nu_n)/prop(nu_n|nu_o), mostly 1, if other ->switch
    double nu_factor = 1.0;
    double nu_new;
    // new value for nu
    if (nu_old == 1.0){
        nu_new = 2.0;
        nu_factor = 0.5;
    } else if (nu_old == 2.0){
        Rcpp::NumericVector r = Rcpp::sample(Rcpp::NumericVector::create(1.0,2.0), 2, false);
        if(r(0) == 1.0){
            nu_new = 3.0;
        } else {
            nu_new = 1.0;
            nu_factor = 2.0;
        }
    } else if(nu_old == (numax - 1.0)) {
        Rcpp::NumericVector r = Rcpp::sample(Rcpp::NumericVector::create(1.0,2.0), 2, false);

        if(r(0) == 1.0){
            nu_new = numax - 2.0 ;
        } else { // t -> normal
            nu_new = numax;
            nu_factor = 2.0;
        }
    } else if (nu_old == numax) {
        nu_new = numax - 1.0;
        nu_factor = 0.5;

    } else {
        Rcpp::NumericVector r = Rcpp::sample(Rcpp::NumericVector::create(1.0,2.0), 2, false);
        if(r(0) == 1.0){
            nu_new = nu_old + 1.0;
        } else {
            nu_new = nu_old - 1.0;
        }
    }

    // calculate parameters for full conditional of phi #Achtung: was bedeutet das
    model_output.slot("m_total_Ig1") = arma::as_scalar( arma::ones( 1, n_genes) * arma::conv_to<arma::mat>::from(Ig_vec) );
    arma::mat yxtb = data - trans(trans(design) * Betas);
    arma::mat var_F = pow(yxtb, 2); // .^ matlab equivalent zu ^2 R , benötigt zur berechnung von h*
    // shape parameter g*, see update phi
    double g_old = (nu_old + 1) * 0.5;
    double g_new = (nu_new + 1) * 0.5;
    // rate parameter h* -> GxN matrix , see update phi

    arma::mat h_old = 0.5*(nu_old* arma::ones(var_F.n_rows, var_F.n_cols) + taue*var_F);
    arma::mat h_new = 0.5*(nu_new* arma::ones(var_F.n_rows, var_F.n_cols) + taue*var_F);

    // Acceptance probability
    double logaccp;
    
    // Two version to calculate logaccp are given which differ in how 
    // logaccp is calculated when nu-old = numax | nu-new = numax
    // Where nu lies between 1 and nu-max where nu-max is the cutoff between a t-distributed or a normal distributed noise model. \textbf{ToDo: Why is nu discrete, was ist grid}. 
    // The prior for the compared nu values which is non-informative except for nu-old = 1 or nu-old = nu-max where nu-new is weighted higher resulting in a higher probability to be accepted but remaining the possibility to be rejected. 
    // 
    // During re-implementation, the opposite behavior was detected. 
    // code within original programming start and stop represents this original behavior 
    // The algorithm firstly, always rejects nu-new after accepting nu-max and secondly always accepts nu-new when it is nu-max.
    // 
    // Behavior originates in these cases being calculated differently. 
    // This is necessary because changing between different distributions requires substantial parameter changes.
    // Attempts to resolve this, were partially sucessfull . 
    // The current algorithm version does always accept nu-new when nu-old = nu-max or nu-new = nu-max. 
    // See current programming 16.6.21
    // This can be seen as sufficient work in progress behavior. Implementation of the intended behavior is planned for the future.  
    
    
    // current programming: start (16.6.21)
    if( nu_old == numax | nu_new == numax ){
        logaccp = arma::as_scalar( arma::ones(1, n_genes) * (g_old * arma::log(h_old)) * arma::ones(n_exp, 1) - std::lgamma(g_old) + std::lgamma(0.5*nu_old) - n_exp*n_genes*0.5*nu_old*std::log(0.5*nu_old) ); // lgamma: Rcpp sugar function
    } else {
        double gammdiff = std::lgamma(g_new) - std::lgamma(g_old) + std::lgamma(0.5*nu_old) - std::lgamma(0.5*nu_new);
        double nudiff = 0.5*(nu_new * std::log(0.5*nu_new) - nu_old * std::log(0.5*nu_old));
        arma::mat Hdiff = (g_old*log(h_old) - g_new*log(h_new));
        double hdiff = arma::as_scalar(arma::ones(1, n_genes) * Hdiff * arma::ones(n_exp, 1));
        logaccp = std::log(nu_factor) + n_exp*n_genes*(gammdiff+nudiff)+hdiff;
    }
    // current programming: stop (16.6.21)
    
    // original programming start
    // if(nu_new == numax){
    //     logaccp = arma::as_scalar( arma::ones(1, n_genes) * (g_old * log(h_old)) * arma::ones(n_exp, 1) - std::lgamma(g_old) + std::lgamma(0.5*nu_old) - n_exp*n_genes*0.5*nu_old*std::log(0.5*nu_old) ); // lgamma: Rcpp sugar function
    // 
    // } else if(nu_old == numax){
    //     logaccp = arma::as_scalar( arma::ones(1, n_genes) * (-g_new * log(h_new)) * arma::ones(n_exp, 1) + std::lgamma(g_new) - std::lgamma(0.5*nu_new) + n_exp*n_genes*0.5*nu_new*std::log(0.5*nu_new) ); // std::lgamma: Rcpp sugar function
    // } else {
    //     double gammdiff = std::lgamma(g_new) - std::lgamma(g_old) + std::lgamma(0.5*nu_old) - std::lgamma(0.5*nu_new);
    //     double nudiff = 0.5*(nu_new * std::log(0.5*nu_new) - nu_old * std::log(0.5*nu_old));
    //     arma::mat Hdiff = (g_old*log(h_old) - g_new*log(h_new));
    //     double hdiff = arma::as_scalar(arma::ones(1, n_genes) * Hdiff * arma::ones(n_exp, 1));
    //     logaccp = std::log(nu_factor) + n_exp*n_genes*(gammdiff+nudiff)+hdiff;
    // }
    // original programming stop
   

    ////////////////////////////////////////////////////////////////////////
    double u = R::runif(0,1);
    
    if (u < std::exp(logaccp)){
        model_output.slot("m_nu") = nu_new;

        if(nu_new == numax){
            model_output.slot("m_phi_scaling") = arma::ones(n_genes, n_exp);
        } else {
            // Hnew ist eine matrix --> jeder wert soll für berechnung von phi_new genommen werden mit funktion matrix_rgamma()
            h_new = 1/h_new;
            arma::mat phinew = matrix_rgamma(n_genes, n_exp, g_new , h_new); // Achtung: dirty fix!!!
            model_output.slot("m_phi_scaling") = phinew;
        }
    }
    return model_output;
    
}

/*** R

    
*/
