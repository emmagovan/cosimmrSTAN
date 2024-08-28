#' Print a cosimmrSTAN output object
#'
#' @param x An object of class \code{cosimmrSTAN_output}
#' @param ... Other arguments (not supported)
#'
#' @return Returns a neat summary of the object
#'
#' @seealso  \code{\link{cosimmr_stan}} for creating
#' \code{cosimmrSTAN_output} objects
#' @export
print.cosimmr_output <-
  function(x, ...) {
    if (inherits(x, "cosimmrSTAN_output") == TRUE) {
        print(x$input)
        message("The input data has been run via cosimmr_stan and has produced ")
        message(nrow(x$output$BUGSoutput$sims.list$sigma), " samples.")
      }
  }
