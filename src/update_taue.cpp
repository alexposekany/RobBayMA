// --------------------------- Description ---------------------------------- //
// Update error precision parameter taue
// Input: ModelInput, ModelOutput 
// Output: double taue
// -------------------------------------------------------------------------- //


#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double update_taue_cpp(Rcpp::S4 &model_input, Rcpp::S4 &model_output){
    // ------------------- extract data from input -------------------------- //
    // Rcpp::Rcout << "Address of a: %model_output: " <<  &model_output << std::endl;
    arma::mat data = Rcpp::as<arma::mat>(model_input.slot("m_data"));
    arma::mat design = Rcpp::as<arma::mat>(model_input.slot("m_design_mat"));
    double hyper_g = Rcpp::as<double>(model_input.slot("m_hyper_g"));
    double hyper_h = Rcpp::as<double>(model_input.slot("m_hyper_h"));

    arma::mat Betas = Rcpp::as<arma::mat>(model_output.slot("m_Betas"));
    arma::mat phi_scaling = Rcpp::as<arma::mat>(model_output.slot("m_phi_scaling"));
    
    int n_genes = data.n_rows;
    int n_exp = data.n_cols;
    
    // ------------------------ start update ----------------------------------- //
    // g*
    double g = hyper_g + 0.5 * n_exp * n_genes ;

    // h*
    arma::mat yxtb = data -  trans(trans(design) * Betas); // yxtb: Term in der Klammer zum Update von hyperparameter h
    
    arma::mat e = square(yxtb) % phi_scaling; // e: elementwise multiplication von yxtb2 und phi
    
    double f = arma::as_scalar(arma::ones(1, n_genes) * e * arma::ones(n_exp, 1)); // f: calculate sum

    // finalize calculation
    double h = 1 / (hyper_h + 0.5 * f);

    // update parameter taue, taue follows Ga(g,1/h) distribution
    double taue = R::rgamma(g, h);
    
    return taue;
    }




/*** R

*/
