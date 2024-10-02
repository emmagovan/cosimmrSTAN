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
             "beta_random_boxplot",
             "covariates_plot"
           ),
           obs = 1,
           binwidth = 0.05,
           alpha = 0.5,
           title = NULL,
           cov_name = NULL,
           one_plot = TRUE,
           n_pred_samples = 1000,
           ...) {
    if(inherits(x, "cosimmrSTAN_output") == TRUE){
      title_input = title
      # Get the specified type
      type <- match.arg(type, several.ok = TRUE)

      # if(x$input$intercept == TRUE){
      #   n_fixed_cov =  dim(x$output$beta_fixed)[2]#ncol(x$input$x_scaled) -1
      # }else{
      #   n_fixed_cov = ncol(x$input$x_scaled)
      # }
      n_fixed_cov =  dim(x$output$beta_fixed)[2]

      if(is.null(x$output$beta_random)  & (("beta_random_histogram" %in% type)|("beta_random_boxplot" %in% type))){
        message("You have selected a plot type that requires random effects but
the model does not contain random effects. Please reselect
model or plots to create and rerun.")
      } else if(n_fixed_cov == 0 & (("beta_fixed_histogram" %in% type)|("beta_fixed_boxplot" %in% type))){
      message("You have selected a plot type that requires fixed effects but
the model does not contain fixed effects. Please reselect
model or plots to create and rerun.")}else{



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


      #  total_cov = ncol(x$input$x_scaled)


   if(is.null(x$output$beta_random)){
     total_levels = NULL
   }else{

      n_random_cov = length(x$input$re_names_order)
      total_cov = n_random_cov + n_fixed_cov
      levels_random_cov = x$input$re_levels



total_levels = sum(levels_random_cov)
cov_name_random = c(rep(NA, total_levels))

random_names = c(rep(NA, n_random_cov))
random_names = x$input$re_names_order#colnames(x$input$original_x)[(n_fixed_cov+1):(n_fixed_cov + n_random_cov)]
cov_name_random = rep(random_names, levels_random_cov)
}


        if("beta_fixed_histogram" %in% type){
      #Could just do a loop and plot all betas - want to loop over number of 2nd dimension in x$output$beta_fixed

          for(l in 1:n_fixed_cov){


          # if(x$input$intercept == TRUE){
          #   c_name = colnames(x$input$x_scaled)[1+l]
          # }else{
          #
          # }
            c_name = colnames(x$input$x_scaled)[l]

          # if(x$input$intercept == TRUE){
          #
          #  out_all_beta = x$output$beta_fixed[,(l+1),] #n samples x K x L
          #
          # } else{
            out_all_beta = x$output$beta_fixed[,(l),] #n samples x K x L

#}

          # out_all_beta = beta[cov_ind,,]
          colnames(out_all_beta) = x$input$source_names
          #I don't actually understand what this is doing
          df_beta <- reshape2::melt(out_all_beta)
          colnames(df_beta) = c("Num", "Source", "Beta")



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

        }



          if("beta_fixed_boxplot" %in% type){
            if(x$input$intercept == TRUE){
              c_name = colnames(x$input$x_scaled)[1+l]
            }else{
              c_name = colnames(x$input$x_scaled)[l]
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



            for(l in 1:n_fixed_cov){

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

        if("beta_random_histogram" %in% type){
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
        }

        if("beta_random_boxplot" %in% type){
          for(l in 1:total_levels){
            c_name = cov_name_random[l]
            level_name = colnames(x$input$X_random)[l]
            out_all_beta = x$output$beta_random[,l,] #n samples x K x L



            # out_all_beta = beta[cov_ind,,]
            colnames(out_all_beta) = x$input$source_names
            #I don't actually understand what this is doing
            df_beta <- reshape2::melt(out_all_beta)
            colnames(df_beta) = c("Num", "Source", "Beta")


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


  if("covariates_plot" %in% type){
    ### The matrix we want to basically make is x$input$original_x
    for(i in 1:length(cov_name)){
      cov_name_in_use = cov_name[i]
      if (!(cov_name_in_use %in% colnames(x$input$covariates_df))) {
        stop("Covariate name not found in the original matrix.")
      }

      # Copy the original matrix and select the column of interest
      target_column <- x$input$covariates_df[[cov_name_in_use]]


      if (is.numeric(target_column)) {
        # For numeric columns: create a sequence of 5000 values from min to max
        new_column <- seq(min(target_column, na.rm = TRUE), max(target_column, na.rm = TRUE), length.out = n_pred_samples)

      } else if (is.factor(target_column)) {
        # For factor columns: repeat each level 1000 times
        new_column <- as.factor(rep(levels(target_column), each = 1)) #This is 1 because p/beta already have 10000 samples built into them, increasing this does nothing bc its just the same groups.

      } else {
        stop("Unsupported column type.")
      }

      n_rows <- length(new_column)

      # Initialize a placeholder for the new prediction matrix
      pred_matrix <- as.data.frame(matrix(ncol = ncol(x$input$covariates_df), nrow = n_rows))
      colnames(pred_matrix) <- colnames(x$input$covariates_df)

    for (i in seq_len(ncol(x$input$covariates_df))) {
      if (colnames(x$input$covariates_df)[i] == cov_name_in_use) {
        pred_matrix[,i] <- new_column
      } else {
        if (is.numeric(x$input$covariates_df[[i]])) {
          # Set 1s for numeric columns - set to use the median ???
          pred_matrix[,i] <- rep(median(x$input$covariates_df[[i]]), n_rows)
        } else if (is.factor(x$input$covariates_df[[i]])) {
          # Set first level for factor columns
          pred_matrix[,i] <- as.factor(rep(levels(x$input$covariates_df[[i]])[1], n_rows))
        }
      }
    }

    ## Now we want to use the predict function
      p_predict = array(NA, dim = c(n_rows, dim(x$output$p)[1], x$input$n_sources))

      for(i in 1:n_rows){
        p_predict[i,,] = (predict(x, as.data.frame(pred_matrix[i,])))$p[1,,] #The 1 here is bc its doing 1 row at a time so p only has 1 entry in that dimension
      }


      # if(is.factor(x$input$original_x[[cov_name_in_use]])){
      #   #This is if it is a factor so we do boxplots
      #
      # } else{
        #This is if its not a factor so we do line plots
        source_chosen = c(x$input$source_names) #Just do all the sources

          source_loop = x$input$n_sources


        if(one_plot == TRUE){

          # #so now we have the source name (s_name) and the number of the source (source_n)
          # n_samples = length(x$output$BUGSoutput$sims.list$p[1,,1])
          # n_ind = x$input$n_obs

          #check which plot we want to make - boxplot or lineplot
          box_or_line = NULL


          cov_selected_col = matrix(x$input$covariates_df[[cov_name_in_use]], ncol = 1)#matrix(c((x$input$covariates_df)[,l]), ncol = 1)

          colnames(cov_selected_col) = cov_name_in_use#colnames(x$input$covariates_df)[l]

          #Just chooses which plot we're making
          if(is.numeric(cov_selected_col)){
            box_or_line = "LINE"
          }else {
            box_or_line = "BOX"
          }

          if(box_or_line == "LINE"){
            #Want to use predict function and a regular grid of the covariate that we want




            n_sources = x$input$n_sources
            #line_array = array(NA, dim = c(n_rows, x$output$vb_control$n_samples, n_sources))
            mean_line_mat = matrix(NA, nrow = n_rows, ncol = n_sources)
            save_sd =  matrix(NA, nrow = n_rows, ncol = n_sources)

            for(s in 1:n_sources){
             # line_array[,,s] = matrix(p_predict[,,s], nrow = n_rows, ncol = x$output$vb_control$n_samples)
              mean_line_mat[,s] = rowMeans(p_predict[,,s])
              for(i in 1:n_rows){
                save_sd[i,s] = sd(p_predict[i,,s])
              }
            }

            out_all_mean <- mean_line_mat
            colnames(out_all_mean) <- x$input$source_names
            df_mean <- reshape2::melt(out_all_mean)

            colnames(df_mean) <- c("Num", "Source", "Mean")

            out_all_sd <- save_sd
            colnames(out_all_sd) <- x$input$source_names
            df_sd <- reshape2::melt(out_all_sd)

            colnames(df_sd) <- c("Num", "Source", "SD")

            df_plot = data.frame(mean = df_mean$Mean,
                                 sd = df_sd$SD,
                                 cov = new_column,
                                 psd = df_mean$Mean + 2*df_sd$SD,
                                 nsd = df_mean$Mean - 2*df_sd$SD,
                                 num = df_mean$Num,
                                 Source = df_mean$Source)



            g <- ggplot(data = df_plot, aes(x = cov, y = mean, colour = Source)) +
              geom_ribbon(data = df_plot, aes(ymin = nsd, ymax = psd, fill = Source), alpha = alpha) +
              geom_line() + xlab(colnames(cov_selected_col)) +
              ylab(paste("Proportion (\u00B1 2sd)")) + ggtitle(paste0("Change in consumption over ", colnames(cov_selected_col)))

            print(g)


          }else if(box_or_line == "BOX"){
            message("Cannot facet wrap boxplot")

          }
        }
        else if(one_plot == FALSE){
          for(s in 1:source_loop){
              s_name = x$input$source_names[s]


            #so now we have the source name (s_name) and the number of the source (source_n)
            n_samples = dim(x$output$p)[1]
            n_ind = x$input$n_obs

            #check which plot we want to make - boxplot or lineplot
            box_or_line = NULL


            cov_selected_col = matrix(c((x$input$covariates_df)[[cov_name_in_use]]), ncol = 1)

            colnames(cov_selected_col) = cov_name_in_use

            #Just chooses which plot we're making
            if(is.numeric(cov_selected_col)){
              box_or_line = "LINE"
            }else {
              box_or_line = "BOX"
            }

            if(box_or_line == "LINE"){

              line_mat= matrix(p_predict[,,s], nrow = n_pred_samples, ncol = n_samples)
              mean_line_mat = c(rep(NA, n_pred_samples))
              mean_line_mat = rowMeans(line_mat)
              save_sd = c(rep(NA, n_pred_samples))
              for(i in 1:n_pred_samples){
                save_sd[i] = sd(line_mat[i,])
              }


              df_plot = data.frame(mean = mean_line_mat,
                                   sd = save_sd,
                                   cov = new_column,
                                   psd = (mean_line_mat + 2*save_sd),
                                   nsd = (mean_line_mat - 2*save_sd))




              g <- ggplot(data = df_plot, aes(x = cov, y = mean)) +
                geom_ribbon(data = df_plot, aes(ymin = nsd, ymax = psd), alpha = alpha) +
                geom_line() +
                ggtitle(paste0("Proportion changing for ", s_name, " consumption over ", colnames(cov_selected_col))) +
                xlab(paste0(colnames(cov_selected_col))) + ylab(paste("Proportion (\u00B1 2sd)"))

              print(g)


            }else if(box_or_line == "BOX"){
              #First make this the length of each group
              n_groups = length(unique(cov_selected_col))
              g_names = c(unique(cov_selected_col))
              ind_mat = matrix(nrow = n_samples, ncol = n_groups)

          for(i in 1:n_groups){
            ind_mat[,i] = p_predict[i,,s]
          }

              ind_vec =c(ind_mat)

              g_names_mat = matrix(nrow = n_samples, ncol = n_groups)
              for(i in 1:n_groups){
                g_names_mat[,i] = c(rep(paste0(g_names[i]), n_samples))
              }




              g_names_rep_vec = c(g_names_mat)

              df_plot = data.frame(samples = ind_vec, Group = g_names_rep_vec)


              g = ggplot(data = df_plot, aes(x = Group, y = samples, colour = Group))  +
                geom_boxplot() +
                ggtitle(paste0(s_name, " consumption over ",colnames(cov_selected_col), " covariate")) +
                xlab(paste0(colnames(cov_selected_col))) + ylab("Proportion")

              print(g)

            }

          }




        }

      #}



    }






  }


        #cov loop bracket

        if (exists("g")) invisible(g)
      }
    }

  }
