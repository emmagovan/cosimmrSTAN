#' Plot different features of an object created from  \code{\link{predict.cosimmrSTAN_output}}.
#'
#' This function allows for different types of plots of the simmr output
#' created from \code{\link{cosimmr_stan}}.
#'
#' The matrix plot should form a necessary part of any SIMM analysis since it
#' allows the user to judge which sources are identifiable by the model.
#' Further detail about these plots is provided in the vignette.
#'
#' @param x An object of class \code{cosimmrSTAN_pred_out} created via
#'  \code{\link{predict}}.
#' @param type The type of plot required. Can be one or more of 'isospace',
#' 'prob_histogram', 'prob_density'
#' @param binwidth The width of the bins for the histogram. Defaults to 0.05
#' @param alpha The degree of transparency of the plots. Not relevant for
#' matrix plots
#' @param title The title of the plot.
#' @param obs The observation you wish to plot
#' @param ...  Currently not used
#'
#' @return one or more of 'isospace', 'prop_histogram', 'prop_density'
#'
#' @import ggplot2
#' @import graphics
#' @import viridis
#' @importFrom reshape2 "melt"
#' @importFrom stats "cor"
#'
#' @author Emma Govan <emmagovan@@gmail.com>>, Andrew Parnell
#' @seealso See  \code{\link{cosimmr_stan}} for
#' creating objects suitable for this function, and many more examples. See
#' also \code{\link{cosimmrSTAN_load}} for creating simmr objects,
#' \code{\link{plot.cosimmrSTAN_input}} for creating isospace plots.
#'
#' @examples
#'
#' \donttest{
#' # A simple example with 10 observations, 2 tracers and 4 sources
#'
#' # The data
#' data(geese_data_day1)
#'
#' # Load into simmr
#' cosimmr_1 <- with(
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
#' # Plot
#' plot(cosimmr_1)
#'
#'
#' # FFVB run
#' cosimmr_1_out <- cosimmr_stan(cosimmr_1)
#'
#' }
#' @export
plot.cosimmrSTAN_pred_out <-
  function(x,
           type = c(
             "isospace",
             "prop_histogram",
             "prop_density"
           ),
           obs = 1,
           binwidth = 0.05,
           alpha = 0.5,
           title = NULL,
           ...) {
    if(inherits(x, "cosimmrSTAN_pred_out") == TRUE){
      title_input = title
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)

      # Iso-space plot is special as all groups go on one plot
      # Add in extra dots here as they can be sent to this plot function
      if ("isospace" %in% type) {
        if(is.null(title_input) == TRUE){
          title = "isospace plot"
        } else{title = title_input}
        graphics::plot(x$input, title = title, ...)

      }


      for(i in 1:length(obs)){


        #Need to have a separate matrix for each ind value
        #So do all this in loop and repeat I think is easiest

        if(is.null(title_input) == TRUE){
          title = c(rep(NA, length(obs)))
          for(j in 1:length(obs)){
            title[j] = paste("Proportions: Prediction", obs[j])
          }
        } else{title = rep(title_input, length(obs))}

        curr_ind = obs[i]
        out_all_p = x$p[curr_ind,,]


        colnames(out_all_p) = x$input$source_names

        df <- reshape2::melt(out_all_p)


        colnames(df) = c("Num", "Source", "Proportion")


        if("prop_histogram" %in% type){
          g <- ggplot(df, aes(
            x = Proportion,
            fill = Source
          )) +
            scale_fill_viridis(discrete = TRUE) +
            geom_histogram(binwidth = binwidth, alpha = alpha) +
            theme_bw() +
            ggtitle(title[i]) +
            facet_wrap("~ Source") +
            theme(legend.position = "none")
          print(g)


        }

        if("prop_density" %in% type){
          g <- ggplot(df, aes(
            x = Proportion,
            fill = Source
          )) +
            scale_fill_viridis(discrete = TRUE) +
            geom_density(aes(y = after_stat(density)), alpha = alpha, linetype = 0) +
            theme_bw() +
            theme(legend.position = "none") +
            ggtitle(title[i]) +
            ylab("Density") +
            facet_wrap("~ Source")
          print(g)
        }

      }




    if (exists("g")) invisible(g)
    }
  }

