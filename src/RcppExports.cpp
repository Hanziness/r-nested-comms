// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// is_2k2_neigh
bool is_2k2_neigh(Rcpp::NumericVector n1, Rcpp::NumericVector n2, int v1, int v2);
RcppExport SEXP _nested_comms_is_2k2_neigh(SEXP n1SEXP, SEXP n2SEXP, SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< int >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(is_2k2_neigh(n1, n2, v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// nested_direction
int nested_direction(Rcpp::NumericVector n1, Rcpp::NumericVector n2, int v1, int v2);
RcppExport SEXP _nested_comms_nested_direction(SEXP n1SEXP, SEXP n2SEXP, SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< int >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(nested_direction(n1, n2, v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// nestedness_value
double nestedness_value(Rcpp::NumericVector n1, Rcpp::NumericVector n2, int v1, int v2);
RcppExport SEXP _nested_comms_nestedness_value(SEXP n1SEXP, SEXP n2SEXP, SEXP v1SEXP, SEXP v2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type v1(v1SEXP);
    Rcpp::traits::input_parameter< int >::type v2(v2SEXP);
    rcpp_result_gen = Rcpp::wrap(nestedness_value(n1, n2, v1, v2));
    return rcpp_result_gen;
END_RCPP
}
// get_lk_all_topsort
List get_lk_all_topsort(List neighs);
RcppExport SEXP _nested_comms_get_lk_all_topsort(SEXP neighsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type neighs(neighsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_lk_all_topsort(neighs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nested_comms_is_2k2_neigh", (DL_FUNC) &_nested_comms_is_2k2_neigh, 4},
    {"_nested_comms_nested_direction", (DL_FUNC) &_nested_comms_nested_direction, 4},
    {"_nested_comms_nestedness_value", (DL_FUNC) &_nested_comms_nestedness_value, 4},
    {"_nested_comms_get_lk_all_topsort", (DL_FUNC) &_nested_comms_get_lk_all_topsort, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_nested_comms(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
