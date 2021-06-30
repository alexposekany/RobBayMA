#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include "matrix_rgamma.h"
// #include <gperftools/profiler.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
arma::mat matrix_rgamma(int n_genes, int n_exp, double shape, arma::mat scaling_mat) {
    // -------------------------- Description ------------------------------- //
    // Helper function creating matrix containing random numbers following gamma
    // distribution which are parameterized with constant shape parameter and 
    // scale parameter stored in matrix 
    // 
    // Input:       matrix dimensions (n_genes, n_exp), constant shape parameter (shape), 
    //              matrix containing all scale parameters (scaling_mat)
    // Output:      n_genes x n_exp matrix containing gamma distributed random variables   
    // ---------------------------------------------------------------------- //
    
    int vec_length = n_genes*n_exp;
    
    arma::mat output = arma::mat(n_genes, n_exp, arma::fill::zeros);
    arma::vec temp(output.memptr(), output.n_elem, false, false);
    
    for (int i = 0; i<vec_length; i++){
        temp(i) = R::rgamma(shape, scaling_mat(i) );
    }
    return output;
}




/*** R


*/
