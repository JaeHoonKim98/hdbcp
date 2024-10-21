// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_cov2_mxPBF_single
arma::mat cpp_cov2_mxPBF_single(const arma::mat& X, const arma::mat& Y, double a0, double b0, double log_gamma);
RcppExport SEXP _hdbcp_cpp_cov2_mxPBF_single(SEXP XSEXP, SEXP YSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP log_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< double >::type log_gamma(log_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_cov2_mxPBF_single(X, Y, a0, b0, log_gamma));
    return rcpp_result_gen;
END_RCPP
}
// cpd_cov_mxPBF
arma::vec cpd_cov_mxPBF(const arma::mat& X, double a0, double b0, int nw, double alp, int n_threads);
RcppExport SEXP _hdbcp_cpd_cov_mxPBF(SEXP XSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP nwSEXP, SEXP alpSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< int >::type nw(nwSEXP);
    Rcpp::traits::input_parameter< double >::type alp(alpSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpd_cov_mxPBF(X, a0, b0, nw, alp, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// simulate_mxPBF_cov
arma::mat simulate_mxPBF_cov(const arma::mat& data, double a0, double b0, int nw, const arma::vec& alps, int n_samples, int n_threads);
RcppExport SEXP _hdbcp_simulate_mxPBF_cov(SEXP dataSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP nwSEXP, SEXP alpsSEXP, SEXP n_samplesSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< int >::type nw(nwSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alps(alpsSEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_mxPBF_cov(data, a0, b0, nw, alps, n_samples, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// cpp_mean2_mxPBF_single
arma::vec cpp_mean2_mxPBF_single(const arma::mat& X, const arma::mat& Y, double log_gamma);
RcppExport SEXP _hdbcp_cpp_mean2_mxPBF_single(SEXP XSEXP, SEXP YSEXP, SEXP log_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type log_gamma(log_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mean2_mxPBF_single(X, Y, log_gamma));
    return rcpp_result_gen;
END_RCPP
}
// cpd_mean_mxPBF
arma::vec cpd_mean_mxPBF(const arma::mat& X, int nw, double alp, int n_threads);
RcppExport SEXP _hdbcp_cpd_mean_mxPBF(SEXP XSEXP, SEXP nwSEXP, SEXP alpSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nw(nwSEXP);
    Rcpp::traits::input_parameter< double >::type alp(alpSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpd_mean_mxPBF(X, nw, alp, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// simulate_mxPBF_mean
arma::mat simulate_mxPBF_mean(const arma::mat& data, int nw, const arma::vec& alps, int n_samples, int n_threads);
RcppExport SEXP _hdbcp_simulate_mxPBF_mean(SEXP dataSEXP, SEXP nwSEXP, SEXP alpsSEXP, SEXP n_samplesSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type nw(nwSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type alps(alpsSEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_mxPBF_mean(data, nw, alps, n_samples, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// cpp_mvrnorm
arma::mat cpp_mvrnorm(int n, const arma::vec& mu, const arma::mat& sigma);
RcppExport SEXP _hdbcp_cpp_mvrnorm(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_mvrnorm(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_hdbcp_cpp_cov2_mxPBF_single", (DL_FUNC) &_hdbcp_cpp_cov2_mxPBF_single, 5},
    {"_hdbcp_cpd_cov_mxPBF", (DL_FUNC) &_hdbcp_cpd_cov_mxPBF, 6},
    {"_hdbcp_simulate_mxPBF_cov", (DL_FUNC) &_hdbcp_simulate_mxPBF_cov, 7},
    {"_hdbcp_cpp_mean2_mxPBF_single", (DL_FUNC) &_hdbcp_cpp_mean2_mxPBF_single, 3},
    {"_hdbcp_cpd_mean_mxPBF", (DL_FUNC) &_hdbcp_cpd_mean_mxPBF, 4},
    {"_hdbcp_simulate_mxPBF_mean", (DL_FUNC) &_hdbcp_simulate_mxPBF_mean, 5},
    {"_hdbcp_cpp_mvrnorm", (DL_FUNC) &_hdbcp_cpp_mvrnorm, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_hdbcp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
