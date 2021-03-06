// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// approximatelyEqualRel
bool approximatelyEqualRel(double a, double b, double relEpsilon);
RcppExport SEXP _robbayma_approximatelyEqualRel(SEXP aSEXP, SEXP bSEXP, SEXP relEpsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type relEpsilon(relEpsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(approximatelyEqualRel(a, b, relEpsilon));
    return rcpp_result_gen;
END_RCPP
}
// approximatelyEqualAbsRel
bool approximatelyEqualAbsRel(double a, double b, double absEpsilon, double relEpsilon);
RcppExport SEXP _robbayma_approximatelyEqualAbsRel(SEXP aSEXP, SEXP bSEXP, SEXP absEpsilonSEXP, SEXP relEpsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type absEpsilon(absEpsilonSEXP);
    Rcpp::traits::input_parameter< double >::type relEpsilon(relEpsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(approximatelyEqualAbsRel(a, b, absEpsilon, relEpsilon));
    return rcpp_result_gen;
END_RCPP
}
// group_data1
Rcpp::List group_data1(arma::mat& data, arma::vec& grouping_vec);
RcppExport SEXP _robbayma_group_data1(SEXP dataSEXP, SEXP grouping_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grouping_vec(grouping_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(group_data1(data, grouping_vec));
    return rcpp_result_gen;
END_RCPP
}
// group_data2
Rcpp::List group_data2(arma::mat& data, arma::vec& grouping_vec);
RcppExport SEXP _robbayma_group_data2(SEXP dataSEXP, SEXP grouping_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type grouping_vec(grouping_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(group_data2(data, grouping_vec));
    return rcpp_result_gen;
END_RCPP
}
// substr
Rcpp::List substr(arma::mat& Betas, arma::vec& Ig_vec);
RcppExport SEXP _robbayma_substr(SEXP BetasSEXP, SEXP Ig_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type Betas(BetasSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Ig_vec(Ig_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(substr(Betas, Ig_vec));
    return rcpp_result_gen;
END_RCPP
}
// get_design
Rcpp::NumericMatrix get_design(Rcpp::NumericVector grouping_vec, int n_groups);
RcppExport SEXP _robbayma_get_design(SEXP grouping_vecSEXP, SEXP n_groupsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type grouping_vec(grouping_vecSEXP);
    Rcpp::traits::input_parameter< int >::type n_groups(n_groupsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_design(grouping_vec, n_groups));
    return rcpp_result_gen;
END_RCPP
}
// matrix_rgamma
arma::mat matrix_rgamma(int n_genes, int n_exp, double shape, arma::mat scaling_mat);
RcppExport SEXP _robbayma_matrix_rgamma(SEXP n_genesSEXP, SEXP n_expSEXP, SEXP shapeSEXP, SEXP scaling_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_genes(n_genesSEXP);
    Rcpp::traits::input_parameter< int >::type n_exp(n_expSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type scaling_mat(scaling_matSEXP);
    rcpp_result_gen = Rcpp::wrap(matrix_rgamma(n_genes, n_exp, shape, scaling_mat));
    return rcpp_result_gen;
END_RCPP
}
// save_iteration_cpp
Rcpp::S4 save_iteration_cpp(Rcpp::S4& mcmc_output, Rcpp::S4& model_output, arma::mat su);
RcppExport SEXP _robbayma_save_iteration_cpp(SEXP mcmc_outputSEXP, SEXP model_outputSEXP, SEXP suSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type mcmc_output(mcmc_outputSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4& >::type model_output(model_outputSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type su(suSEXP);
    rcpp_result_gen = Rcpp::wrap(save_iteration_cpp(mcmc_output, model_output, su));
    return rcpp_result_gen;
END_RCPP
}
// update_betai_cpp
Rcpp::S4 update_betai_cpp(Rcpp::S4 model_input, Rcpp::S4 model_output);
RcppExport SEXP _robbayma_update_betai_cpp(SEXP model_inputSEXP, SEXP model_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type model_input(model_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type model_output(model_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(update_betai_cpp(model_input, model_output));
    return rcpp_result_gen;
END_RCPP
}
// update_lambda_cpp
double update_lambda_cpp(Rcpp::S4 model_input, Rcpp::S4 model_output);
RcppExport SEXP _robbayma_update_lambda_cpp(SEXP model_inputSEXP, SEXP model_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type model_input(model_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type model_output(model_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(update_lambda_cpp(model_input, model_output));
    return rcpp_result_gen;
END_RCPP
}
// update_mhnuphi_cpp
Rcpp::S4 update_mhnuphi_cpp(Rcpp::S4& model_input, Rcpp::S4& model_output);
RcppExport SEXP _robbayma_update_mhnuphi_cpp(SEXP model_inputSEXP, SEXP model_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type model_input(model_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4& >::type model_output(model_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(update_mhnuphi_cpp(model_input, model_output));
    return rcpp_result_gen;
END_RCPP
}
// update_rjbetai_cpp
Rcpp::List update_rjbetai_cpp(arma::mat data, arma::mat design, double total_mu, int n_groups, arma::mat Betas, arma::mat phi_scaling, arma::vec Ig_vec, double taue, double p, arma::mat gene_mu_grouped, arma::mat m_tau, int current_gene, double lambda);
RcppExport SEXP _robbayma_update_rjbetai_cpp(SEXP dataSEXP, SEXP designSEXP, SEXP total_muSEXP, SEXP n_groupsSEXP, SEXP BetasSEXP, SEXP phi_scalingSEXP, SEXP Ig_vecSEXP, SEXP taueSEXP, SEXP pSEXP, SEXP gene_mu_groupedSEXP, SEXP m_tauSEXP, SEXP current_geneSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type design(designSEXP);
    Rcpp::traits::input_parameter< double >::type total_mu(total_muSEXP);
    Rcpp::traits::input_parameter< int >::type n_groups(n_groupsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Betas(BetasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phi_scaling(phi_scalingSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ig_vec(Ig_vecSEXP);
    Rcpp::traits::input_parameter< double >::type taue(taueSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gene_mu_grouped(gene_mu_groupedSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type m_tau(m_tauSEXP);
    Rcpp::traits::input_parameter< int >::type current_gene(current_geneSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(update_rjbetai_cpp(data, design, total_mu, n_groups, Betas, phi_scaling, Ig_vec, taue, p, gene_mu_grouped, m_tau, current_gene, lambda));
    return rcpp_result_gen;
END_RCPP
}
// update_rjnufine_cpp
Rcpp::S4 update_rjnufine_cpp(Rcpp::S4& model_input, Rcpp::S4& model_output);
RcppExport SEXP _robbayma_update_rjnufine_cpp(SEXP model_inputSEXP, SEXP model_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type model_input(model_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4& >::type model_output(model_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(update_rjnufine_cpp(model_input, model_output));
    return rcpp_result_gen;
END_RCPP
}
// update_taue_cpp
double update_taue_cpp(Rcpp::S4& model_input, Rcpp::S4& model_output);
RcppExport SEXP _robbayma_update_taue_cpp(SEXP model_inputSEXP, SEXP model_outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4& >::type model_input(model_inputSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4& >::type model_output(model_outputSEXP);
    rcpp_result_gen = Rcpp::wrap(update_taue_cpp(model_input, model_output));
    return rcpp_result_gen;
END_RCPP
}
// update_wdbetai_cpp
Rcpp::List update_wdbetai_cpp(arma::mat data, arma::mat design, double total_mu, int n_groups, arma::mat Betas, arma::mat phi_scaling, arma::vec Ig_vec, double taue, int current_gene, double lambda);
RcppExport SEXP _robbayma_update_wdbetai_cpp(SEXP dataSEXP, SEXP designSEXP, SEXP total_muSEXP, SEXP n_groupsSEXP, SEXP BetasSEXP, SEXP phi_scalingSEXP, SEXP Ig_vecSEXP, SEXP taueSEXP, SEXP current_geneSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type design(designSEXP);
    Rcpp::traits::input_parameter< double >::type total_mu(total_muSEXP);
    Rcpp::traits::input_parameter< int >::type n_groups(n_groupsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Betas(BetasSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phi_scaling(phi_scalingSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Ig_vec(Ig_vecSEXP);
    Rcpp::traits::input_parameter< double >::type taue(taueSEXP);
    Rcpp::traits::input_parameter< int >::type current_gene(current_geneSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(update_wdbetai_cpp(data, design, total_mu, n_groups, Betas, phi_scaling, Ig_vec, taue, current_gene, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_robbayma_approximatelyEqualRel", (DL_FUNC) &_robbayma_approximatelyEqualRel, 3},
    {"_robbayma_approximatelyEqualAbsRel", (DL_FUNC) &_robbayma_approximatelyEqualAbsRel, 4},
    {"_robbayma_group_data1", (DL_FUNC) &_robbayma_group_data1, 2},
    {"_robbayma_group_data2", (DL_FUNC) &_robbayma_group_data2, 2},
    {"_robbayma_substr", (DL_FUNC) &_robbayma_substr, 2},
    {"_robbayma_get_design", (DL_FUNC) &_robbayma_get_design, 2},
    {"_robbayma_matrix_rgamma", (DL_FUNC) &_robbayma_matrix_rgamma, 4},
    {"_robbayma_save_iteration_cpp", (DL_FUNC) &_robbayma_save_iteration_cpp, 3},
    {"_robbayma_update_betai_cpp", (DL_FUNC) &_robbayma_update_betai_cpp, 2},
    {"_robbayma_update_lambda_cpp", (DL_FUNC) &_robbayma_update_lambda_cpp, 2},
    {"_robbayma_update_mhnuphi_cpp", (DL_FUNC) &_robbayma_update_mhnuphi_cpp, 2},
    {"_robbayma_update_rjbetai_cpp", (DL_FUNC) &_robbayma_update_rjbetai_cpp, 13},
    {"_robbayma_update_rjnufine_cpp", (DL_FUNC) &_robbayma_update_rjnufine_cpp, 2},
    {"_robbayma_update_taue_cpp", (DL_FUNC) &_robbayma_update_taue_cpp, 2},
    {"_robbayma_update_wdbetai_cpp", (DL_FUNC) &_robbayma_update_wdbetai_cpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_robbayma(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
