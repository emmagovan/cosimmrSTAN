#'Predicts proportion of each source in a mixture, based on values provided for covariates
#'
#'
#'
#' @param object An object of class \code{cosimmrSTAN_output} created via the
#' function \code{\link{cosimmr_stan}}
#' @param x_pred A data.frame of covariate values that the user wishes
#' to predict source proportions for, provided in the same order that the
#' original covariance matrix was. Important for this to be a data.frame otherwise
#' numeric values can be set as characters and this causes incorrect calculations.
#' @param n_output the number of posterior samples to generate. Defaults to 3600.
#' @param ... Other arguments (not used)
#'
#' @return object of class 'cosimmrSTAN_pred_out'
#'
#' @author Emma Govan <emmagovan@@gmail.com> Andrew Parnell
#'
#' @seealso \code{\link{cosimmr_load}} for creating objects suitable for this
#' function,  and
#' \code{\link{plot.cosimmr_output}} for plotting output.
#'
#' @references Andrew C. Parnell, Donald L. Phillips, Stuart Bearhop, Brice X.
#' Semmens, Eric J. Ward, Jonathan W. Moore, Andrew L. Jackson, Jonathan Grey,
#' David J. Kelly, and Richard Inger. Bayesian stable isotope mixing models.
#' Environmetrics, 24(6):387–399, 2013.
#'
#' @importFrom R2jags jags
#'
#' @examples
#' \donttest{
#' ## See the package vignette for a detailed run through of these 4 examples
#'
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' data(geese_data_day1)
#' cov_1 = c(1,2,3,2,3,1,1,1,2)
#' simmr_1 <- with(
#'   geese_data_day1,
#'   cosimmrSTAN_load(
#'     formula = mixtures ~ cov_1,
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
#' plot(simmr_1)
#'
#' # Print
#' simmr_1
#'
#' # FFVB run
#' simmr_1_out <- cosimmr_stan(simmr_1)
#'
#' # Print it
#' print(simmr_1_out)
#'
#'
#' # Plot
#' plot(simmr_1_out, type = "isospace")
#' plot(simmr_1_out, type = "beta_fixed_histogram")
#'
#' x_pred = data.frame(cov_1 = c(1,5))
#'
#'pred_array<-predict(simmr_1_out, x_pred)
#'
#' }
#' @export
predict.cosimmrSTAN_output <- function(object,
                                   x_pred, ...) {

  #Makes sure the object is the correct class
  if(inherits(object, "cosimmrSTAN_output") == TRUE){

    #It has to be a data.frame if theres numeric and categorical data otherwise matrices
    #just turn everything into categorical so it wont work
    if(inherits(x_pred, "data.frame") == FALSE) stop("x_pred must be of type `data.frame` to make predictions with")

    if(object$input$random_effects == TRUE){

      #Deal with the fixed effects first
      scale_x = object$input$scale_x
      #Creating a max and min vector so we can check if the new x_pred falls outside the range of the original data
      #Create vectors here but do comparison after all data has been scaled
      max_vec = c(rep(NA, ncol(object$input$x_scaled)))
      min_vec = c(rep(NA, ncol(object$input$x_scaled)))

      for(i in 1:(ncol(object$input$x_scaled))){
        max_vec[i] = max(object$input$x_scaled[,i])
        min_vec[i] = min(object$input$x_scaled[,i])
      }
      K = object$input$n_sources
      n_tracers = object$input$n_tracers
      n_covariates = ncol(object$input$x_scaled)
      mixtures = object$input$mixtures

      #So now what we want to do is to add x_pred onto the bottom of the original x matrix,
      #scale all, then remove original x mat

      original_x = data.frame(object$input$covariates_df)

      #Add check that col names match
      #Otherwise next line wont work

      # if(colnames(x_pred) != colnames(original_x))stop("Column names of original data and data you wish to predict with do not match. Please fix and rerun.")
      colnames(original_x) = colnames(x_pred)

      new_x = rbind(original_x, x_pred)


      # This adds a column of ones for predicting with if the original had an
      #intercept - want to do this after
      # if(simmr_out$input$intercept == TRUE){
      #   new_x = cbind(c(rep(1,nrow(new_x))), new_x)
      # } else if(simmr_out$input$intercept == FALSE){
      #   new_x = new_x
      # }

      #Now scale


      #So now what we want to do is to add x_pred onto the bottom of the original x matrix,
      #scale all, then remove original x mat


      if(object$input$intercept == TRUE){
        if(ncol(object$input$X_fixed) == 1){
          #Then there is nothing to do for this one - we only have an intercept
          new_x = object$input$X_fixed
        }else{
        new_x = rbind(original_x, x_pred)[1:(ncol(object$input$X_fixed) -1)]
        }
      } else {
        new_x = rbind(original_x, x_pred)[1:(ncol(object$input$X_fixed))]
      }

      # This adds a column of ones for predicting with if the original had an
      #intercept - want to do this after
      # if(simmr_out$input$intercept == TRUE){
      #   new_x = cbind(c(rep(1,nrow(new_x))), new_x)
      # } else if(simmr_out$input$intercept == FALSE){
      #   new_x = new_x
      # }

      #Now scale

      if(nrow(mixtures) == 1){
        #This is if its just 1 entry
        new_x == new_x
      } else{
        if(scale_x == TRUE){
          if(object$input$intercept == TRUE){
            if(ncol(object$input$X_fixed) == 1){
              #This is if there is only an intercept
              x_pred_mat = as.matrix(object$input$X_fixed[1:nrow(x_pred),])
            } else{
            # Original code
            ncol_scaled =  (ncol(stats::model.matrix(~ ., data=new_x))) - 1
            scaled_full_mat = matrix(scale(stats::model.matrix(~ ., data=new_x),
                                           center = c(1,object$input$scaled_center),
                                           scale = c(1, object$input$scaled_scale))[,-c(1)], ncol = ncol_scaled)
            scaled_full_mat = cbind(c(rep(1,nrow(scaled_full_mat))), scaled_full_mat)

            x_pred_mat = matrix(scaled_full_mat[-c(1:nrow(original_x)),], ncol = ncol(scaled_full_mat))
}
          }else if(object$input$intercept == FALSE){
            scaled_full_mat = scale(stats::model.matrix(~ . -1, data=new_x),
                                    center = object$input$scaled_center,
                                    scale = object$input$scaled_scale)


            x_pred_mat = matrix(scaled_full_mat[-c(1:nrow(original_x)),], ncol = ncol(scaled_full_mat))

          }

        }else if(scale_x == FALSE){
          if(object$input$intercept == TRUE){
            if(ncol(new_x) ==1){ #This is if the fixed part is just an intercept, we want to keep this as is??
              x_pred_mat = matrix(1, ncol = 1)

              }else{
            scaled_full_mat = (stats::model.matrix(~ ., data=new_x))

            x_pred_mat = scaled_full_mat[-c(1:nrow(original_x)),]
              }
          }else if(object$input$intercept == FALSE){
            scaled_full_mat = stats::model.matrix(~ .-1, data=new_x)

            x_pred_mat = scaled_full_mat[-c(1:nrow(original_x)),]
          }

        }
      }



      #Checks that all the values are above or equal to the min and below or equal to the max
      for(j in 1:(nrow(x_pred_mat))){
        for(i in 1:(ncol(object$input$x_scaled))){
          if(x_pred_mat[j,i] >= min_vec[i] & x_pred_mat[j,i] <= max_vec[i]){
            #message("Data falls within range of data used in original model, okay to predict with")
            print_err = FALSE
          } else(print_err = TRUE)
        }
      }

      #This is separate because otherwise its inside the loop and it prints a bunch of times
      if(print_err){message("Please note: The data you wish to predict with falls outside the range of data used in the original model")}



      #### NOW DO RANDOM EFFECTS
      #I think I can basically do the same thing, add the correct columns to the bottom of the matrix? And then just pull out the correct ones??

      re_names = object$input$re_names_order
      n_re = length(re_names)
      levels = object$input$re_levels
      nlevels = sum(levels)


      if(object$input$intercept == TRUE){
        n_fixed = ncol(x_pred_mat) - 1
      } else {
        n_fixed = ncol(x_pred_mat)
      }

      re_x_pred = as.matrix(x_pred[,(n_fixed+1):ncol(x_pred)]) #This is the random effects part that we want to estimate from x_pred


      col_numbers = matrix(NA, ncol = n_re, nrow = nrow(x_pred))
      re_mat = as.matrix(object$input$covariates_df[,(n_fixed+1):ncol(x_pred)])
      colnames(col_numbers) = colnames(object$input$covariates_df)[(n_fixed+1):ncol(x_pred)]



      for(i in 1:n_re){
        for(j in 1:nrow(x_pred)){

          col_numbers[j,i] = which((re_mat[,i]) == (re_x_pred[j,i]))[1]
        }
      }

      reordered_matrix <- col_numbers[, re_names, drop = FALSE]

      #So now we want to make the actual x matrix by extracting the right rows
      extract_and_bind_rows <- function(effect_matrix, x_random, rows_range) {
        # Initialize a list to store the extracted rows
        extracted_rows <- list()

        # Iterate over each entry in the effect matrix
        for (i in 1:nrow(effect_matrix)) {
          for (j in 1:ncol(effect_matrix)) {
            # Get the indices for the row to extract and the range
            row_index <- effect_matrix[i, j]
            range <- rows_range[[j]]  # Get the range for the current column

            # Extract the specified rows from x_random
            extracted <- x_random[row_index, range, drop = FALSE]

            # Add the extracted rows to the list
            extracted_rows[[length(extracted_rows) + 1]] <- extracted
          }
        }

        # Combine all extracted rows into a single matrix


        # Ensure the result matrix has the correct number of columns
        # if (ncol(result_matrix) != length(unlist(rows_range))) {
        #   stop("Error: The number of columns in the result matrix does not match the expected number.")
        # }

        return(extracted_rows)
      }

      # Function to create cumulative column ranges from a vector
      create_cumulative_ranges <- function(range_vector) {
        start <- 1
        rows_range <- list()

        for (i in seq_along(range_vector)) {
          end <- start + range_vector[i] - 1
          rows_range[[i]] <- start:end
          start <- end + 1
        }

        return(rows_range)
      }

      rows_range = create_cumulative_ranges(levels)

      x_rows = extract_and_bind_rows(reordered_matrix, object$input$X_random, rows_range)

      Z_pred_mat = matrix(data = NA, nrow = nrow(x_pred), ncol = nlevels)
      combined_rows <- list()
      for (i in seq(1, length(x_rows), by = n_re)) {
        # Bind the rows in the current group
        group <- x_rows[i:(i + n_re - 1)]
        combined <- do.call(cbind, group)
        combined_rows[[length(combined_rows) + 1]] <- combined
      }

      for(i in 1:nrow(x_pred)){
        Z_pred_mat[i,] = combined_rows[[i]]
      }

      sigma = object$output$sigma #n_samples x n_isotopes
      n_output = nrow(sigma)
      p_sample = array(NA, dim =  c(nrow(x_pred_mat), n_output, K))

      beta_fixed = object$output$beta_fixed #n_samples * n_cov x K
      beta_random = object$output$beta_random

      f <- array(NA, dim = c(nrow(x_pred_mat), K, n_output))

      for(i in 1:nrow(x_pred_mat)){
        for(k in 1:K){
          for(s in 1:n_output){
            f[i,k,s] = x_pred_mat[i,] %*% beta_fixed[s,,k] + Z_pred_mat[i,] %*% beta_random[s,,k]
          }
        }
      }

      for(j in 1:n_output){
        for (n_obs in 1:nrow(x_pred_mat)) {
          p_sample[n_obs,j, ] <- exp(f[n_obs,1:K, j]) / (sum((exp(f[n_obs,1:K, j]))))
        }
      }

      out<-list(
        p = p_sample,
        sigma = sigma,
        input = object$input
      )

      class(out) = "cosimmrSTAN_pred_out"
      return(out)


    }else{

    scale_x = object$input$scale_x


    #Creating a max and min vector so we can check if the new x_pred falls outside the range of the original data
    #Create vectors here but do comparison after all data has been scaled
    max_vec = c(rep(NA, ncol(object$input$x_scaled)))
    min_vec = c(rep(NA, ncol(object$input$x_scaled)))

    for(i in 1:(ncol(object$input$x_scaled))){
      max_vec[i] = max(object$input$x_scaled[,i])
      min_vec[i] = min(object$input$x_scaled[,i])
    }

    K = object$input$n_sources
    n_tracers = object$input$n_tracers
    n_covariates = ncol(object$input$x_scaled)
    mixtures = object$input$mixtures

    #So now what we want to do is to add x_pred onto the bottom of the original x matrix,
    #scale all, then remove original x mat

    original_x = data.frame(object$input$covariates_df)

    #Add check that col names match
    #Otherwise next line wont work

    # if(colnames(x_pred) != colnames(original_x))stop("Column names of original data and data you wish to predict with do not match. Please fix and rerun.")
    colnames(original_x) = colnames(x_pred)

    new_x = rbind(original_x, x_pred)


    #Now scale

    if(nrow(mixtures) == 1){
      #This is if its just 1 entry
      new_x == new_x
    } else{
      if(scale_x == TRUE){
        if(object$input$intercept == TRUE){
          # Original code
          ncol_scaled =  (ncol(stats::model.matrix(~ ., data=new_x))) - 1
          scaled_full_mat = matrix(scale(stats::model.matrix(~ ., data=new_x),
                                         center = c(1,object$input$scaled_center),
                                         scale = c(1, object$input$scaled_scale))[,-c(1)], ncol = ncol_scaled)
          scaled_full_mat = cbind(c(rep(1,nrow(scaled_full_mat))), scaled_full_mat)

          x_pred_mat = matrix(scaled_full_mat[-c(1:nrow(original_x)),], ncol = ncol(scaled_full_mat))

        }else if(object$input$intercept == FALSE){
          scaled_full_mat = scale(stats::model.matrix(~ . -1, data=new_x),
                                  center = object$input$scaled_center,
                                  scale = object$input$scaled_scale)


          x_pred_mat = matrix(scaled_full_mat[-c(1:nrow(original_x)),], ncol = ncol(scaled_full_mat))

        }

      }else if(scale_x == FALSE){
        if(object$input$intercept == TRUE){
          scaled_full_mat = stats::model.matrix(~ ., data=new_x)

          x_pred_mat = matrix(scaled_full_mat[-c(1:nrow(original_x)),], nrow = nrow(x_pred))
        }else if(object$input$intercept == FALSE){
          scaled_full_mat = stats::model.matrix(~ .-1, data=new_x)

          x_pred_mat = matrix(scaled_full_mat[-c(1:nrow(original_x)),], nrow = nrow(x_pred))
        }

      }
    }



    #Checks that all the values are above or equal to the min and below or equal to the max
    for(j in 1:(nrow(x_pred_mat))){
      for(i in 1:(ncol(object$input$x_scaled))){
        if(x_pred_mat[j,i] >= min_vec[i] & x_pred_mat[j,i] <= max_vec[i]){
          #message("Data falls within range of data used in original model, okay to predict with")
          print_err = FALSE
        } else(print_err = TRUE)
      }
    }

    #This is separate because otherwise its inside the loop and it prints a bunch of times
    if(print_err){message("Please note: The data you wish to predict with falls outside the range of data used in the original model")}


    sigma = object$output$sigma #n_samples x n_isotopes
    n_output = nrow(sigma)
    p_sample = array(NA, dim =  c(nrow(x_pred_mat), n_output, K))

    beta_fixed = object$output$beta_fixed #n_samples * n_cov x K


    f <- array(NA, dim = c(nrow(x_pred_mat), K, n_output))

    for(i in 1:nrow(x_pred_mat)){
    for(k in 1:K){
    for(s in 1:n_output){
      f[i,k,s] = x_pred_mat[i,] %*% beta_fixed[s,,k]
    }
    }
    }

    for(j in 1:n_output){
      for (n_obs in 1:nrow(x_pred_mat)) {
        p_sample[n_obs,j, ] <- exp(f[n_obs,1:K, j]) / (sum((exp(f[n_obs,1:K, j]))))
      }
    }

    out<-list(
      p = p_sample,
      sigma = sigma,
      input = object$input
    )

    class(out) = "cosimmrSTAN_pred_out"
    return(out)
  }
  }else (message("Can only predict using cosimmrSTAN_output object generated from
cosimmr_stan function"))
}
