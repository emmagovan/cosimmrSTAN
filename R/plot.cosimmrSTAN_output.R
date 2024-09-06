#' Plot different features of an object created from  \code{\link{cosimmr_stan}}.
#'
#' This function allows for 4 different types of plots of the simmr output
#' created from \code{\link{cosimmr_stan}}. The
#' types are: plot of beta values
#'
#' The matrix plot should form a necessary part of any SIMM analysis since it
#' allows the user to judge which sources are identifiable by the model.
#' Further detail about these plots is provided in the vignette.
#'
#' @param x An object of class \code{cosimmrSTAN_output} created via
#'  \code{\link{cosimmr_stan}}.
#' @param type The type of plot required. Can be one or more of 'isospace',
#' 'prop_histogram', 'prop_density', 'beta_fixed_histogram', 'beta_fixed_boxplot',
#' 'beta_random_histogram', 'beta_random_boxplot'
#' @param obs The observation number you wish to plot
#' @param binwidth The width of the bins for the histogram. Defaults to 0.05
#' @param alpha The degree of transparency of the plots. Not relevant for
#' matrix plots
#' @param title The title of the plot.
#' @param ...  Currently not used
#'
#' @return one or more of 'isospace', 'beta_fixed_histogram',
#' 'beta_fixed_boxplot', 'prop_histogram', 'prop_density', 'beta_random_boxplot'
#' or 'beta_random_histogram'
#'
#' @import ggplot2
#' @import graphics
#' @import viridis
#' @importFrom reshape2 "melt"
#' @importFrom stats "cor"
#'
#' @author Emma Govan <emmagovan@@gmail.com>, Andrew Parnell
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
#' simmr_1 <- with(
#'   geese_data_day1,
#'   cosimmr_load(
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
#' plot(simmr_1)
#'
#'
#' # FFVB run
#' simmr_1_out <- cosimmr_ffvb(simmr_1)
#'
#'plot(simmr_1_out, type = c("isospace", "beta_hist"))
#' }
#' @export
plot.cosimmrSTAN_output <-
  function(x,
           type = c(
             "isospace",
             "prop_histogram",
             "prop_density",
             "beta_fixed_histogram",
             "beta_fixed_boxplot",
             "beta_random_histogram",
             "beta_random_boxplot"
           ),
           obs = 1,
           binwidth = 0.05,
           alpha = 0.5,
           title = NULL,
           ...) {
    if(inherits(x, "cosimmrSTAN_output") == TRUE){
      title_input = title
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)

#want to extract the right covariate number using the name?

      if(is.null(x$output$beta_random)  & (("beta_random_histogram" %in% type)|("beta_random_boxplot" %in% type))){
        message("You have selected a plot type that requires random effects but
the model does not contain random effects. Please reselect
model or plots to create and rerun.")
      } else{



        # Iso-space plot is special as all groups go on one plot
        # Add in extra dots here as they can be sent to this plot function
        if ("isospace" %in% type) {
          if(is.null(title_input) == TRUE){
            title = "isospace plot"
          } else{title = title_input}
          graphics::plot(x$input, title = title, ...)

        }





        #Need to have a separate matrix for each ind value
        #So do all this in loop and repeat I think is easiest
        for(i in 1:length(obs)){
          if(is.null(title_input) == TRUE){
            title = c(rep(NA, length(obs)))
            for(j in 1:length(obs)){
              title[j] = paste("Proportions: Observation", obs[j])
            }
          } else{title = rep(title_input, length(obs))}
          curr_ind = obs[i]
          out_all_p = x$output$p[,curr_ind,]


          colnames(out_all_p) = x$input$source_names

          df <- reshape2::melt(out_all_p)


          colnames(df) = c("Num", "Source", "Proportion")

          #add other plot types here maybe
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



        #Prep data
        #Data needs to be edited I think to make life easier

      ## EAsiest thing for betas might just be to plot all of them? And then they can choose what they care about
    #  x$output$beta_fixed # n samples x K x L fixed covariates
     # x$output$beta_random # n samples x K x L random covariates - these will only ever be groups


        total_cov = ncol(x$input$original_x)
        if(x$input$intercept == TRUE){
        n_fixed_cov = ncol(x$input$X_fixed) -1
        }else{
          n_fixed_cov = ncol(x$input$X_fixed)
        }
        n_random_cov = total_cov - n_fixed_cov
      levels_random_cov = x$input$re_levels



total_levels = sum(levels_random_cov)
cov_name_random = c(rep(NA, total_levels))

random_names = c(rep(NA, n_random_cov))
random_names = x$input$re_names_order#colnames(x$input$original_x)[(n_fixed_cov+1):(n_fixed_cov + n_random_cov)]
cov_name_random = rep(random_names, levels_random_cov)


      #Could just do a loop and plot all betas
      for(l in 1:n_fixed_cov){


          if(x$input$intercept == TRUE){
            c_name = colnames(x$input$X_fixed)[1+l]
          }else{
            c_name = colnames(x$input$X_fixed)[l]
          }

          if(x$input$intercept == TRUE){

           out_all_beta = x$output$beta_fixed[,(l+1),] #n samples x K x L

          } else{
            out_all_beta = x$output$beta_fixed[,(l),] #n samples x K x L

}

          # out_all_beta = beta[cov_ind,,]
          colnames(out_all_beta) = x$input$source_names
          #I don't actually understand what this is doing
          df_beta <- reshape2::melt(out_all_beta)
          colnames(df_beta) = c("Num", "Source", "Beta")

          if("beta_fixed_histogram" %in% type){

          if(is.null(title_input) == TRUE){

              title = paste("beta histogram plot for", c_name, "covariate")

          } else{title = title_input}


          #Histograms
          g <- ggplot(df_beta, aes(x = Beta)) +
            scale_fill_viridis(discrete = TRUE) +
            geom_histogram(binwidth = binwidth, alpha = alpha) +
            theme_bw() +
            ggtitle(title) +
            facet_wrap("~ Source") +
            theme(legend.position = "none",
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) +
            geom_vline(xintercept = 0, colour = "red")


          #Boxplot

          print(g)
        }




          if("beta_fixed_boxplot" %in% type){

          if(is.null(title_input) == TRUE){

            title = paste("beta boxplot for", c_name, "covariate")

          } else{title = title_input}

          g <- ggplot(df_beta, aes(x = Beta)) +
            scale_fill_viridis(discrete = TRUE) +
            geom_boxplot() +
            theme_bw() +
            ggtitle(title) +
            facet_wrap("~ Source")+
            geom_vline(xintercept = 0, colour = "red") +
            theme(axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())

          print(g)

        }
      }


      for(l in 1:total_levels){


          # if(x$input$intercept == TRUE){
          #   c_name = colnames(x$input$original_x)[n_fixed_cov+l+1]
          # }else{
          #   c_name = colnames(x$input$original_x)[n_fixed_cov+l]
          # }

          c_name = cov_name_random[l]
        level_name = colnames(x$input$X_random)[l]
          out_all_beta = x$output$beta_random[,l,] #n samples x K x L



          # out_all_beta = beta[cov_ind,,]
          colnames(out_all_beta) = x$input$source_names
          #I don't actually understand what this is doing
          df_beta <- reshape2::melt(out_all_beta)
          colnames(df_beta) = c("Num", "Source", "Beta")

          if("beta_random_histogram" %in% type){

          if(is.null(title_input) == TRUE){

            title = paste("beta histogram plot for", c_name, "covariate, level", level_name)

          } else{title = title_input}


          #Histograms
          g <- ggplot(df_beta, aes(x = Beta)) +
            scale_fill_viridis(discrete = TRUE) +
            geom_histogram(binwidth = binwidth, alpha = alpha) +
            theme_bw() +
            ggtitle(title) +
            facet_wrap("~ Source") +
            theme(legend.position = "none",
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) +
            geom_vline(xintercept = 0, colour = "red")


          #Boxplot

          print(g)
        }

        if("beta_random_boxplot" %in% type){


          if(is.null(title_input) == TRUE){

            title = paste("beta boxplot for", c_name, "covariate, level", level_name)

          } else{title = title_input}

          g <- ggplot(df_beta, aes(x = Beta)) +
            scale_fill_viridis(discrete = TRUE) +
            geom_boxplot() +
            theme_bw() +
            ggtitle(title) +
            facet_wrap("~ Source")+
            geom_vline(xintercept = 0, colour = "red") +
            theme(axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())

          print(g)

        }
      }



        #cov loop bracket

        if (exists("g")) invisible(g)
      }
    }

  }
