#' Prints convergence check for the output created with \code{\link{cosimmr_stan}}
#'
#' Produces table of convergence-checking values for an object created
#' with  \code{\link{cosimmr_stan}}. Depending on what algorithm has been
#' run this will either produce a table of pareto-k-diagnostic values
#' (if the algorithm run wasVariational Bayes) or rhat values (if MCMC
#' was used)
#'
#'
#'
#' @param object An object of class \code{cosimmrSTAN_output} produced by the
#' function \code{\link{cosimmr_stan}}
#' @param ...  Not used
#' @return A print out showing either pareto-k-diagnostic values or rhat values summarised
#' @author Emma Govan <emmagovan@@gmail.com> Andrew Parnell
#' @seealso See \code{\link{cosimmr_stan}}for creating objects suitable for
#' this function, and many more examples.
#' See also \code{\link{cosimmrSTAN_load}} for creating cosimmrSTAN objects,
#' \code{\link{plot.cosimmrSTAN_input}} for creating isospace plots,
#' \code{\link{plot.cosimmrSTAN_output}} for plotting output.
#'
#' @importFrom stats sd cor
#'
#' @examples
#' \donttest{
#' # A simple example with 10 observations, 2 tracers and 4 sources
#'
#' # The data
#' data(geese_data_day1)
#' cosimmrSTAN_1 <- with(
#'   geese_data_day1,
#'   cosimmrSTAN_load(
#'     formula = mixtures ~ 1,
#'     source_names = source_names,
#'     source_means = source_means,
#'     source_sds = source_sds,
#'     correction_means = correction_means,
#'     correction_sds = correction_sds,
#'     concentration_means = concentration_means
#'   )
#' )
#'
#' # Plot
#' plot(cosimmrSTAN_1)
#'
#'
#' # FFVB run
#' cosimmrSTAN_1_out <- cosimmr_stan(cosimmrSTAN_1)
#'
#' # Check convergence
#' convergence_check(cosimmrSTAN1_out)
#' }
#' @export
convergence_check <-
  function(object, ...) {
    if(inherits(object, "cosimmrSTAN_output") == TRUE) {
if(is.null(object$output$rhat_summary)){
  #Then we want pareto_k_summary
  print(object$output$pareto_k_summary)
}else if(is.null(object$output$pareto_k_summary)){
  print(object$output$rhat_summary)
}

    } else {
      (return(message("incorrect object passed to function")))
    }
  }


