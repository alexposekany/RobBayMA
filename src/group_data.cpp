#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <map>
#include "group_data.h"
// #include <gperftools/profiler.h>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]



// [[Rcpp::export]]
Rcpp::List group_data1(arma::mat& data, arma::vec& grouping_vec) {
    // -------------------------- Description ------------------------------- //
    // Helper function extracting list of submatrices of input array according to 
    // grouping defined in input grouping vector
    // Input:       Data matrix, grouping vector
    // Output:      Nested list containing matrices of grouped data, 
    //              1st list element --> data of group 1 etc.    
    // ---------------------------------------------------------------------- //
    
    int n_groups = grouping_vec.max();
    int n_exp = data.n_cols;
    
    // create empty list with sublist equal to number of groups 
    Rcpp::List data_grouped = Rcpp::List::create() ;
    for (int sublist_i = 0; sublist_i<n_groups; ++sublist_i){
        data_grouped.push_back( Rcpp::List::create() );
    }
    
    // fill sublists empty matrices
    for (int sublist_i = 0; sublist_i<n_groups; ++sublist_i){
        arma::mat m;
        data_grouped[sublist_i] = m;
    }
    
    
    for (int current_exp = 0; current_exp<n_exp; ++current_exp ){ 
        for (int current_group = 0; current_group<n_groups; current_group++){
            // Ig vector is passed to function as grouping_vec +1
            int current_Ig = grouping_vec(current_exp) - 1; 
            
            if ( current_Ig == current_group){
                arma::mat temp_mat = data_grouped[current_group];
                data_grouped[current_group] = arma::join_rows(temp_mat, data.col(current_exp));
            }
        }
    }
    
    return data_grouped;
}

// [[Rcpp::export]]
Rcpp::List group_data2(arma::mat& data, arma::vec& grouping_vec) {
    // -------------------------- Description ------------------------------- //
    // Helper function extracting list of submatrices of input array according to 
    // grouping defined in input grouping vector
    // Input:       Data matrix, grouping vector
    // Output:      Nested list containing matrices of grouped data, 
    //              1st list element --> data of group 1 etc.    
    // ---------------------------------------------------------------------- //
    
    int n_groups = grouping_vec.max();
    if (n_groups == 1){
        n_groups = 2; 
    }
    
    int n_exp = data.n_cols;
    
    // create empty list with sublist equal to number of groups 
    Rcpp::List data_grouped = Rcpp::List::create() ;
    for (int sublist_i = 0; sublist_i<n_groups; ++sublist_i){
        data_grouped.push_back( Rcpp::List::create() );
    }
    
    // fill sublists empty matrices
    for (int sublist_i = 0; sublist_i<n_groups; ++sublist_i){
        arma::mat m;
        data_grouped[sublist_i] = m;
    }
    
    
    for (int current_exp = 0; current_exp<n_exp; ++current_exp ){ 
        for (int current_group = 0; current_group<n_groups; current_group++){
            // Rcpp::Rcout << "current_group: " << current_group << std::endl;
            // Rcpp::Rcout << "n_groups: " << n_groups << std::endl;
            // Ig vector is passed to function as grouping_vec +1
            int current_Ig = grouping_vec(current_exp) - 1; 
            
            if ( current_Ig == current_group){
                arma::mat temp_mat = data_grouped[current_group];
                data_grouped[current_group] = arma::join_rows(temp_mat, data.col(current_exp));
            }
        }
    }
    
    return data_grouped;
}



// [[Rcpp::export]]
Rcpp::List substr(arma::mat& Betas, arma::vec& Ig_vec) { 
    // extracts subarrays of array Y according to grouping defined in grouping
    // vector s
    int S_groups= Ig_vec.max();
    int N_exp = Betas.n_cols;
    
    // create empty list with sublist equal to number of conditions S_groups
    Rcpp::List YS = Rcpp::List::create() ;
    for (int sublist_i = 0; sublist_i<S_groups; ++sublist_i){
        YS.push_back( Rcpp::List::create() );
    }
    
    // fill sublists empty matrices
    for (int sublist_i = 0; sublist_i<S_groups; ++sublist_i){
        arma::mat m;
        YS[sublist_i] = m;
    }
    
    
    // int co = 0;
    for (int current_exp = 0; current_exp<N_exp; ++current_exp ){ 
        
        for (int current_group = 0; current_group<S_groups; current_group++){
            // Ig vector is passed to function as Ig_vec +1
            int current_Ig = Ig_vec(current_exp) - 1; 
            
            if ( current_Ig == current_group){
                arma::mat temp_mat = YS[current_group];
                YS[current_group] = arma::join_rows(temp_mat, Betas.col(current_exp));
            }
        }
        // co = co + 1;
    }
    
    
    
    return YS;
}





// [[Rcpp::export]]
Rcpp::NumericMatrix get_design(Rcpp::NumericVector grouping_vec, int n_groups) {
    // -------------------------- Description ------------------------------- //
    // Helper function creating matrix encoding experimental design based on grouping vector
    // Input:       grouping vector, number of groups 
    // Output:      groups x experiments matrix encoding design with 0/1
    // ---------------------------------------------------------------------- //
    
    int n_exp = grouping_vec.length();
    Rcpp::NumericMatrix design_mat(n_groups, n_exp); 
    
    for(int current_group = 1; current_group <= n_groups; current_group++){
        for(int current_exp = 0; current_exp < n_exp; current_exp++){
            if(grouping_vec[current_exp]  == current_group){
                design_mat(current_group - 1, current_exp) = 1;
            }
        }
    }
    
    return design_mat;
    
}


/*** R

*/
