// ---------------------------- description --------------------------------- //
// decide if update is within dimension (wd) or cross-dimension (rj) and
// perform the respective update
// Input: S4 classes structurizing data in/output (ModelInput, ModelOutput)  

// Output: ModelOutput with updated variables
// m_Ig_vec             G-dim vector of indicators I_g
// m_Betas              SxG-dim matrix with column vectors beta_g
// m_gene_mu_grouped    SxG-dim matrix  with column vectors mu_g
// m_tau                SxG-dim matrix  with column vectors tau_g
// m_lambda             double 
// -------------------------------------------------------------------------- //

#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include "update_lambda.h"
#include "group_data.h"
#include "update_wdbetai.h"
#include "update_rjbetai.h"
// #include "gperftools/profiler.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::S4 update_betai_cpp(Rcpp::S4 model_input, Rcpp::S4 model_output){
    // ------------------------ start wrapper & helper variables ------------ //
    // ModelInput
    arma::mat data = Rcpp::as<arma::mat>(model_input.slot("m_data"));
    arma::mat design = Rcpp::as<arma::mat>(model_input.slot("m_design_mat"));
    arma::vec hyper_g = Rcpp::as<arma::vec>(model_input.slot("m_hyper_g"));
    arma::vec hyper_h = Rcpp::as<arma::vec>(model_input.slot("m_hyper_h"));
    double total_mu = Rcpp::as<double>(model_input.slot("m_total_mu"));
    int n_groups = Rcpp::as<int>(model_input.slot("m_n_groups"));
    
    // ModelOutput
    arma::mat Betas = Rcpp::as<arma::mat>(model_output.slot("m_Betas"));
    arma::mat phi_scaling = Rcpp::as<arma::mat>(model_output.slot("m_phi_scaling"));
    arma::vec Ig_vec = Rcpp::as<arma::vec>(model_output.slot("m_Ig_vec"));
    double taue = Rcpp::as<double>(model_output.slot("m_taue"));
    double p = Rcpp::as<double>(model_output.slot("m_p"));
    arma::mat gene_mu_grouped = Rcpp::as<arma::mat>(model_output.slot("m_gene_mu_grouped"));
    arma::mat m_tau = Rcpp::as<arma::mat>(model_output.slot("m_tau"));
    
    // things added for lambda  
    double hyper_c = Rcpp::as<double>(model_input.slot("m_hyper_c"));
    double hyper_d = Rcpp::as<double>(model_input.slot("m_hyper_d"));
    double n_i1 = Rcpp::as<double>(model_output.slot("m_total_Ig1"));
    
    int n_genes = data.n_rows;
    
    // --------------------------- lambda update ---------------------------- //
    double lambda_updated = update_lambda_cpp(model_input, model_output);

    // -------------- update betas, Is and mus for every gene --------------- //
    // initialize local storage variables
    Rcpp::NumericVector updated_Igs(n_genes);
    arma::mat updated_Betas = arma::zeros(n_groups, n_genes);
    arma::mat updated_mus = arma::zeros(n_groups, n_genes);
    arma::mat updated_taus = arma::zeros(n_groups, n_genes);

    Rcpp::List current_gene_update = Rcpp::List::create();

    // choose randomly between reversible jump (rj) and within dimension (wd) update
    for(int current_gene = 0; current_gene < n_genes; ++current_gene){
        double r = R::runif(0,1);
        if(r>0.5){
            current_gene_update = update_wdbetai_cpp(data, design, total_mu, n_groups, Betas, phi_scaling, Ig_vec, taue, current_gene, lambda_updated);
        } else {
            current_gene_update = update_rjbetai_cpp(data, design, total_mu, n_groups, Betas, phi_scaling, Ig_vec, taue, p, gene_mu_grouped, m_tau, current_gene, lambda_updated);
        }

        // save results in local storage variables
        updated_Igs[current_gene] = current_gene_update["I"];
        updated_Betas.col(current_gene) = Rcpp::as<arma::vec>(current_gene_update["beta"]);
        updated_mus.col(current_gene) = Rcpp::as<arma::vec>(current_gene_update["mu"]);
        updated_taus.col(current_gene) = Rcpp::as<arma::vec>(current_gene_update["tau"]);

    }

    // save results in storage class
    model_output.slot("m_Ig_vec") = updated_Igs;
    model_output.slot("m_Betas") = updated_Betas;
    model_output.slot("m_gene_mu_grouped") = updated_mus;
    model_output.slot("m_tau") = updated_taus;
    model_output.slot("m_lambda") = lambda_updated;
    
    return model_output;

    }




/*** R
*/