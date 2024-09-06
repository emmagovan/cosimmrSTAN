#' Print simmr input object
#'
#' @param x An object of class \code{cosimmr_input}
#' @param ... Other arguments (not supported)
#'
#'#' @author Emma Govan <emmagovan@@gmail.com> Andrew Parnell
#'
#' @seealso \code{\link{cosimmrSTAN_load}} for creating objects suitable for this
#' function
#'
#' @return A neat presentation of your simmr object.
#' @export print.cosimmrSTAN_input
print.cosimmrSTAN_input <-
  function(x, ...) {
    if(inherits(x, "cosimmrSTAN_input") == TRUE){
    message("This is a valid cosimmrSTAN input object with ")
    message(paste(x$n_obs, " observations, "))
    message(paste(ncol(x$x_scaled), " covariates, "))
    message(paste(x$n_tracers, "tracers, and "))
    message(paste(x$n_sources, "sources.\n"))
    message(" The formula is ")
    print((x$formula))
    message("\nThe source names are: ")
    print(x$source_names, sep = ", ")
    message(".\n")
    message("The tracer names are: ")
    print(colnames(x$mixtures), sep = ", ")
    message("\n\n")
  }
}
