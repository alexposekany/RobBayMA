#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
// #include <gperftools/profiler.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// -------------------------- Description ----------------------------------- //
// Helper function saving output of current iteration in McmcOutput class 
// 
// Input:       mcmc_output, storage class collecting result of every iteration
//              model_output, storage class collecting reult of current itaeration
//              su, storage variable counting how often indicator variable Ig was 1 
//              over all runs
// 
// Output:      mcmc_output, storage class collecting result of every iteration
// -------------------------------------------------------------------------- //

// [[Rcpp::export]]
Rcpp::S4 save_iteration_cpp(Rcpp::S4 &mcmc_output, Rcpp::S4 &model_output, arma::mat su) {
        // nus
        Rcpp::NumericVector nu_vec = Rcpp::as<Rcpp::NumericVector>( mcmc_output.slot("nus"));
        nu_vec.push_back(Rcpp::as<double>(model_output.slot("m_nu")));
        mcmc_output.slot("nus") = nu_vec;

        // taues
        Rcpp::NumericVector taue_vec = Rcpp::as<Rcpp::NumericVector>( mcmc_output.slot("taues"));
        taue_vec.push_back(Rcpp::as<double>(model_output.slot("m_taue")));
        mcmc_output.slot("taues") = taue_vec;

        // ps
        Rcpp::NumericVector p_vec = Rcpp::as<Rcpp::NumericVector>( mcmc_output.slot("ps"));
        p_vec.push_back(Rcpp::as<double>(model_output.slot("m_p")));
        mcmc_output.slot("ps") = p_vec;

        // Igs
        Rcpp::List list_Ig_vec = Rcpp::as<Rcpp::List>( mcmc_output.slot("Igs"));
        list_Ig_vec.push_back(Rcpp::as<Rcpp::NumericVector>(model_output.slot("m_Ig_vec")));
        mcmc_output.slot("Igs") = list_Ig_vec;

        // lambdas
        Rcpp::NumericVector lambda_vec = Rcpp::as<Rcpp::NumericVector>( mcmc_output.slot("lambdas"));
        lambda_vec.push_back(Rcpp::as<double>(model_output.slot("m_lambda")));
        mcmc_output.slot("lambdas") = lambda_vec;

        // sus
         mcmc_output.slot("sus") = su;

    return mcmc_output;

}





/*** R

*/

  