// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif


RcppExport SEXP _rcpp_module_boot_stan_fit4Hierarchical_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4STAN_VB_pandr_raw_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4STAN_VB_pxr_raw_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4STAN_nested_1_random_effect_pandr_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4STAN_nested_1_random_effect_pxr_raw_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4raw_source_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4Hierarchical_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4Hierarchical_mod, 0},
    {"_rcpp_module_boot_stan_fit4STAN_VB_pandr_raw_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4STAN_VB_pandr_raw_mod, 0},
    {"_rcpp_module_boot_stan_fit4STAN_VB_pxr_raw_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4STAN_VB_pxr_raw_mod, 0},
    {"_rcpp_module_boot_stan_fit4STAN_nested_1_random_effect_pandr_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4STAN_nested_1_random_effect_pandr_mod, 0},
    {"_rcpp_module_boot_stan_fit4STAN_nested_1_random_effect_pxr_raw_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4STAN_nested_1_random_effect_pxr_raw_mod, 0},
    {"_rcpp_module_boot_stan_fit4raw_source_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4raw_source_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_cosimmrSTAN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
