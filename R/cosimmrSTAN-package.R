#' cosimmrSTAN: An R package for Stable Isotope Mixing Models
#'
#' cosimmrSTAN is a package that has been developed to allow for running of Stable
#' Isotope Mixing Models in R. It allows for the inclusion of covariates and
#' has been designed to be easy to use for non-expert users. cosimmrSTAN uses
#' Variational Bayes to run SIMMs, instead of MCMC. This allows for faster
#' running of models without any issues with convergence
#'
#' @name cosimmrSTAN
#' @author Emma Govan <emmagovan@@gmail.com>, Andrew Parnell
#' @import Rcpp
#' @import methods
#' @importFrom Rcpp evalCpp
#' @useDynLib cosimmr
#' @aliases cosimmr-package
## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib cosimmrSTAN, .registration = TRUE
## usethis namespace: end
NULL
NULL

