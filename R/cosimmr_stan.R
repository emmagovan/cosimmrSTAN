#' Run a \code{cosimmrSTAN_input} object through STAN
#'
#' This is the main function of cosimmrSTAN. It takes a \code{cosimmrSTAN_input} object
#' created via \code{\link{cosimmrSTAN_load}}, runs it in fixed form
#' Variational Bayes  via STAN to determine the dietary proportions, and then
#' outputs a \code{cosimmrSTAN_output} object for further analysis and plotting
#' via \code{\link{plot.cosimmrSTAN_output}}.
#'
#'@param cosimmrSTAN_in An object created via the function \code{\link{cosimmrSTAN_load}}
#'@param type What type of model to run using STAN. Options are 'STAN_VB
#'@param error_type Whether to use 'processxresidual' error term or
#''process+residual' term. Defaults to 'processxresidual'
#'@param prior_control A list of values including arguments named \code{sigma_shape}
#'(prior values for sigma shape), \code{sigma_rate} (prior values for sigma rate)
#'@param n_samples Number of samples to output. Defaults to 3600.
#'nested within each other. Defaults to FALSE.
#'
#'@return an object of class \code{cosimmrSTAN_output} with two named top-level
#'components: \item{input}{The \code{cosimmrSTAN_input} object given to the
#' \code{cosimmrSTAN_ffvb} function} \item{output}{A set of outputs produced by
#' the FFVB function. These can be analysed using the
#' \code{\link{summary.cosimmrSTAN_output}} and \code{\link{plot.cosimmrSTAN_output}}
#' functions.}
#'
#' @author Emma Govan <emmagovan@@gmail.com>, Andrew Parnell
#'
#' @seealso \code{\link{cosimmrSTAN_load}} for creating objects suitable for this
#' function, \code{\link{plot.cosimmrSTAN_input}} for creating isospace plots,
#' \code{\link{summary.cosimmrSTAN_output}} for summarising output, and
#' \code{\link{plot.cosimmrSTAN_output}} for plotting output.
#'
#' @references Andrew C. Parnell, Donald L. Phillips, Stuart Bearhop, Brice X.
#' Semmens, Eric J. Ward, Jonathan W. Moore, Andrew L. Jackson, Jonathan Grey,
#' David J. Kelly, and Richard Inger. Bayesian stable isotope mixing models.
#' Environmetrics, 24(6):387–399, 2013.
#'
#' Andrew C Parnell, Richard Inger, Stuart Bearhop, and Andrew L Jackson.
#' Source partitioning using stable isotopes: coping with too much variation.
#' PLoS ONE, 5(3):5, 2010.
#'
#'
#' @examples
#' \donttest{
#' ## See the package vignette for a detailed run through of these examples
#'
#' # Data set 1: 10 obs on 2 isos, 4 sources, with tefs and concdep
#' data(geese_data_day1)
#' x =  c(1,2,3,2,1,3,2,1,2)
#' cosimmrSTAN_1 <- with(
#'   geese_data_day1,
#'   cosimmrSTAN_load(
#'     formula = mixtures ~ x,
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
#' # Print
#' cosimmrSTAN_1
#'
#' # FFVB run
#' cosimmrSTAN_1_out <- cosimmr_stan(cosimmrSTAN_1)
#'
#' # Print it
#' print(cosimmrSTAN_1_out)
#'
#' # Summary
#' summary(cosimmrSTAN_1_out, type = "correlations")
#' summary(cosimmrSTAN_1_out, type = "statistics")
#' ans <- summary(cosimmrSTAN_1_out, type = c("quantiles", "statistics"))
#'
#' # Plot
#' plot(cosimmrSTAN_1_out, type = "beta_boxplot")
#' plot(cosimmrSTAN_1_out, type = "beta_histogram")
#'
#'}
#' @export cosimmr_stan
cosimmr_stan <- function(cosimmrSTAN_in,
                         type = "STAN_VB",
                         error_type = "processxresidual",
                         prior_control = list(
                          sigma_shape = c(rep(1, cosimmrSTAN_in$n_tracers)),
                          sigma_rate = c(rep(1, cosimmrSTAN_in$n_tracers)),
                          omicron_shape = c(rep(1, cosimmrSTAN_in$n_tracers)),
                          omicron_rate = c(rep(1, cosimmrSTAN_in$n_tracers))
                         ),
                         n_samples = 3600
                         ){

  #Core detection - potentially have this as an option people can turn on and off?
  # options(mc.cores = parallel::detectCores())

  #Need loops for error type


  #Need to add solo check here
  if (nrow(cosimmrSTAN_in$mixtures) == 1) {
    message("Only 1 mixture value, performing a simmr solo run...\n")
    solo <- 0
  } else {
    solo <- 1
  }


  if(type == "STAN_VB"){


    if(cosimmrSTAN_in$random_effects == TRUE){
        if(error_type == "processxresidual"){
          #This is not nested random effects pxr WORKING
          stan_dat = list(
            J = cosimmrSTAN_in$n_tracers,
            N = cosimmrSTAN_in$n_obs,
            K = cosimmrSTAN_in$n_sources,
            L1 = ncol(cosimmrSTAN_in$X_inner),
            y = cosimmrSTAN_in$mixtures,
            q = cosimmrSTAN_in$concentration_means,
            s_mean = cosimmrSTAN_in$source_means,
            c_mean = cosimmrSTAN_in$correction_means,
            s_sd = cosimmrSTAN_in$source_sds,
            c_sd = cosimmrSTAN_in$correction_sds,
            sigma_shape = prior_control$sigma_shape,
            sigma_rate = prior_control$sigma_rate,
            not_solo = solo,
            X_int = cosimmrSTAN_in$X_intercept,
            X_inner = cosimmrSTAN_in$X_inner,
            omicron_shape = prior_control$omicron_shape,
            omicron_rate = prior_control$omicron_rate
          )


          model = stanmodels$STAN_nested_1_random_effect_pxr_raw


          fit_opt <- rstan::optimizing(model,
                                       data = stan_dat)

          which_beta0 <- grep('beta0', names(fit_opt$par))
          which_beta1 <- grep('beta1', names(fit_opt$par))


          which_sigma <- grep('sigma', names(fit_opt$par))
          which_omicron <- grep('omicron', names(fit_opt$par))

          beta0_start_opt <- structure(fit_opt$par[which_beta0],
                                       dim  = c(stan_dat$K, 1))

          beta1_start_opt <- structure(fit_opt$par[which_beta1],
                                       dim  = c(stan_dat$K, stan_dat$L1))



          sigma_raw_start_opt <- fit_opt$par[which_sigma][1:stan_dat$J]
          omicron_start_opt <- fit_opt$par[which_omicron][1:stan_dat$J]



          fit_vb <- rstan::vb(
            model, data = stan_dat, output_samples = n_samples,
            adapt_iter = 30000,
            init = list('beta0' = beta0_start_opt,
                        'beta1' = beta1_start_opt,
                        'sigma_raw' = sigma_raw_start_opt,
                        'omicron' = omicron_start_opt),
            refresh = FALSE
          )

          extracted_samples = rstan::extract(fit_vb)

          #Want to extract all the betas in a sensible way first I think
          #alpha_ans = extracted_samples$alpha
          beta1_ans = extracted_samples$beta1 # This is n_samples * K * n_covariates
          sigma_ans = extracted_samples$sigma
          omicron_ans = extracted_samples$omicron
          p_sample = extracted_samples$p

          #CONVERT TO F
          #f is x_inner * beta1 + x_outer * beta_2
          #p should be n_ind * n_samples * K
          # f = array(NA, dim = c(cosimmrSTAN_in$n_obs, cosimmrSTAN_in$n_sources, n_samples))
          #
          # for(i in 1:cosimmrSTAN_in$n_obs){
          #   for(k in 1: cosimmrSTAN_in$n_sources){
          #     for(s in 1:n_samples){
          #       f[i,k,s] =  X_inner[i,] %*% beta1_ans[s,k,]
          #     }
          #   }
          # }

          #Then each p is sum f / sum exp f

          # p_sample = array(NA, dim = c(cosimmrSTAN_in$n_obs, n_samples,  cosimmrSTAN_in$n_sources))
          #
          # for(j in 1:n_samples){
          #   for (n_obs in 1:cosimmrSTAN_in$n_obs) {
          #     p_sample[n_obs,j, ] <- exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]) / (sum((exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]))))
          #   }
          # }
          #




        }else if(error_type == "process+residual"){
          #not nested random effects p+r WORKING
          stan_dat = list(
            J = cosimmrSTAN_in$n_tracers,
            N = cosimmrSTAN_in$n_obs,
            K = cosimmrSTAN_in$n_sources,
            L1 = ncol(cosimmrSTAN_in$X_inner),
            y = cosimmrSTAN_in$mixtures,
            q = cosimmrSTAN_in$concentration_means,
            s_mean = cosimmrSTAN_in$source_means,
            c_mean = cosimmrSTAN_in$correction_means,
            s_sd = cosimmrSTAN_in$source_sds,
            c_sd = cosimmrSTAN_in$correction_sds,
            sigma_shape = prior_control$sigma_shape,
            sigma_rate = prior_control$sigma_rate,
            not_solo = solo,
            X_int = cosimmrSTAN_in$X_intercept,
            X_inner = cosimmrSTAN_in$X_inner
          )
          model = stanmodels$STAN_nested_1_random_effect_pandr

          fit_opt <- rstan::optimizing(model,
                                       data = stan_dat)

          which_beta0 <- grep('beta0', names(fit_opt$par))
          which_beta1 <- grep('beta1', names(fit_opt$par))


          which_sigma <- grep('sigma', names(fit_opt$par))
          which_omicron <- grep('omicron', names(fit_opt$par))

          beta0_start_opt <- structure(fit_opt$par[which_beta0],
                                       dim  = c(stan_dat$K, 1))

          beta1_start_opt <- structure(fit_opt$par[which_beta1],
                                       dim  = c(stan_dat$K, stan_dat$L1))



          sigma_raw_start_opt <- fit_opt$par[which_sigma][1:stan_dat$J]
          omicron_start_opt <- fit_opt$par[which_omicron][1:stan_dat$J]



          fit_vb <- rstan::vb(
            model, data = stan_dat, output_samples = n_samples,
            adapt_iter = 30000,
            init = list('beta0' = beta0_start_opt,
                        'beta1' = beta1_start_opt,
                        'sigma_raw' = sigma_raw_start_opt,
                        'omicron' = omicron_start_opt),
            refresh = FALSE
          )

          extracted_samples = rstan::extract(fit_vb)

          #Want to extract all the betas in a sensible way first I think
          #alpha_ans = extracted_samples$alpha
          beta0_ans = extracted_samples$beta0
          beta1_ans = extracted_samples$beta1 # This is n_samples * K * n_covariates
          sigma_ans = extracted_samples$sigma
          omicron_ans = NULL
          p_sample = extracted_samples$p

          #CONVERT TO F
          #f is x_inner * beta1 + x_outer * beta_2
          #p should be n_ind * n_samples * K
          # f = array(NA, dim = c(cosimmrSTAN_in$n_obs, cosimmrSTAN_in$n_sources, n_samples))
          #
          # for(i in 1:cosimmrSTAN_in$n_obs){
          #   for(k in 1: cosimmrSTAN_in$n_sources){
          #     for(s in 1:n_samples){
          #       f[i,k,s] =  cosimmrSTAN_in$X_inner[i,] %*% beta1_ans[s,k,]
          #     }
          #   }
          # }
          #
          # #Then each p is sum f / sum exp f
          #
          # p_sample = array(NA, dim = c(cosimmrSTAN_in$n_obs, n_samples,  cosimmrSTAN_in$n_sources))
          #
          # for(j in 1:n_samples){
          #   for (n_obs in 1:cosimmrSTAN_in$n_obs) {
          #     p_sample[n_obs,j, ] <- exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]) / (sum((exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]))))
          #   }
         # }
        }


    }else if(cosimmrSTAN_in$random_effects == FALSE){
      if(error_type == "processxresidual"){
      #This is not nested fixed effects pxr WORKING
      stan_dat = list(
        J = cosimmrSTAN_in$n_tracers,
        N = cosimmrSTAN_in$n_obs,
        K = cosimmrSTAN_in$n_sources,
        L = cosimmrSTAN_in$n_covariates,
        y = cosimmrSTAN_in$mixtures,
        q = cosimmrSTAN_in$concentration_means,
        s_mean = cosimmrSTAN_in$source_means,
        c_mean = cosimmrSTAN_in$correction_means,
        s_sd = cosimmrSTAN_in$source_sds,
        c_sd = cosimmrSTAN_in$correction_sds,
        sigma_shape = prior_control$sigma_shape,
        sigma_rate = prior_control$sigma_rate,
        not_solo = solo,
        X = cosimmrSTAN_in$x_scaled,
        omicron_shape = prior_control$omicron_shape,
        omicron_rate = prior_control$omicron_rate
      )

      model = stanmodels$STAN_VB_pxr_raw

      fit_opt <- rstan::optimizing(model,
                                   data = stan_dat)


      which_beta <- grep('beta', names(fit_opt$par))


      which_sigma <- grep('sigma', names(fit_opt$par))
      which_omicron <- grep('omicron', names(fit_opt$par))


      beta1_start_opt <- structure(fit_opt$par[which_beta],
                                   dim  = c(stan_dat$K, stan_dat$L))



      sigma_raw_start_opt <- fit_opt$par[which_sigma][1:stan_dat$J]
      omicron_start_opt <- fit_opt$par[which_omicron][1:stan_dat$J]


      # fit_vb <- rstan::vb(
      #   model, data = stan_dat,
      #   algorithm = 'fullrank',
      #   pars = c('beta', 'sigma', 'omicron'),
      #   init = list('beta' = beta1_start_opt,
      #               'sigma' = sigma_raw_start_opt,
      #               'omicron' = omicron_start_opt),
      #   tol_rel_obj = 0.0000001, #convergence tolerance on the relative norm of the objective
      #   output_samples = n_samples
      # )

      fit_vb <- rstan::vb(
        model, data = stan_dat, output_samples = n_samples,
        adapt_iter = 30000,
        init = list('beta' = beta1_start_opt,
                    'sigma_raw' = sigma_raw_start_opt,
                    'omicron' = omicron_start_opt),
        refresh = FALSE
      )



      extracted_samples = rstan::extract(fit_vb)

      #Want to extract all the betas in a sensible way first I think
      #alpha_ans = extracted_samples$alpha
      beta1_ans = extracted_samples$beta # This is n_samples * K * n_covariates
      sigma_ans = extracted_samples$sigma
      omicron_ans = extracted_samples$omicron
      beta0_ans = NULL
      p_sample = extracted_samples$p

      #CONVERT TO F
      #f is x_inner * beta1 + x_outer * beta_2
      #p should be n_ind * n_samples * K
      # f = array(NA, dim = c(cosimmrSTAN_in$n_obs, cosimmrSTAN_in$n_sources, n_samples))
      #
      # for(i in 1:cosimmrSTAN_in$n_obs){
      #   for(k in 1: cosimmrSTAN_in$n_sources){
      #     for(s in 1:n_samples){
      #       f[i,k,s] =  cosimmrSTAN_in$x_scaled[i,] %*% beta1_ans[s,k,]
      #     }
      #   }
      # }
      #
      # #Then each p is sum f / sum exp f
      #
      # p_sample = array(NA, dim = c(cosimmrSTAN_in$n_obs, n_samples,  cosimmrSTAN_in$n_sources))
      #
      # for(j in 1:n_samples){
      #   for (n_obs in 1:cosimmrSTAN_in$n_obs) {
      #     p_sample[n_obs,j, ] <- exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]) / (sum((exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]))))
      #   }
      # }



    }else if(error_type == "process+residual"){
      #not nested fixed effects p+r WORKING

      stan_dat = list(
        J = cosimmrSTAN_in$n_tracers,
        N = cosimmrSTAN_in$n_obs,
        K = cosimmrSTAN_in$n_sources,
        L = cosimmrSTAN_in$n_covariates,
        y = cosimmrSTAN_in$mixtures,
        q = cosimmrSTAN_in$concentration_means,
        s_mean = cosimmrSTAN_in$source_means,
        c_mean = cosimmrSTAN_in$correction_means,
        s_sd = cosimmrSTAN_in$source_sds,
        c_sd = cosimmrSTAN_in$correction_sds,
        sigma_shape = prior_control$sigma_shape,
        sigma_rate = prior_control$sigma_rate,
        not_solo = solo,
        X = cosimmrSTAN_in$x_scaled
      )

      model = stanmodels$STAN_VB_pandr_raw


      fit_opt <- rstan::optimizing(model,
                                   data = stan_dat)
      #
      # which_alpha <- grep('alpha', names(fit_opt$par))
       which_beta1 <- grep('beta', names(fit_opt$par))
      #
      #
       which_sigma <- grep('sigma', names(fit_opt$par))
      #
      # alpha_start_opt <- structure(fit_opt$par[which_alpha],
      #                              dim  = c(stan_dat$K, 1))
      #
       beta1_start_opt <- structure(fit_opt$par[which_beta1],
                                    dim  = c(stan_dat$K, stan_dat$L))

       sigma_raw_start_opt <- fit_opt$par[which_sigma][1:2]


      # fit_vb <- rstan::vb(
      #   model, data = stan_dat,
      #   algorithm = 'meanfield',
      #   pars = c('beta', 'sigma'),
      #   init = list('beta' = beta1_start_opt,
      #               'sigma' = sigma_raw_start_opt),
      #   tol_rel_obj = 0.00001, #convergence tolerance on the relative norm of the objective
      #   output_samples = n_samples,
      #   elbo_samples = 100,
      #   iter = 10000
      # )

       # fit_vb = rstan::vb(model, data = stan_dat,
       #    pars = c('beta', 'sigma'),
       #    init = list('beta' = beta1_start_opt,
       #               'sigma' = sigma_raw_start_opt),
       #    output_samples = n_samples,
       #    adapt_iter = 30000,
       #    refresh = FALSE)

       fit_vb <- rstan::vb(
         model, data = stan_dat, output_samples = n_samples,
         adapt_iter = 30000,
         init = list('beta' = beta1_start_opt,
                     'sigma_raw' = sigma_raw_start_opt),
         refresh = FALSE
       )


      extracted_samples = rstan::extract(fit_vb)

      #Want to extract all the betas in a sensible way first I think
      #alpha_ans = extracted_samples$alpha
      beta1_ans = extracted_samples$beta # This is n_samples * K * n_covariates
      sigma_ans = extracted_samples$sigma
      omicron_ans = NULL
      beta0_ans = NULL
      p_sample = extracted_samples$p

      #CONVERT TO F
      #f is x_inner * beta1 + x_outer * beta_2
      #p should be n_ind * n_samples * K
      # f = array(NA, dim = c(cosimmrSTAN_in$n_obs, cosimmrSTAN_in$n_sources, n_samples))
      #
      # for(i in 1:cosimmrSTAN_in$n_obs){
      #   for(k in 1: cosimmrSTAN_in$n_sources){
      #     for(s in 1:n_samples){
      #       f[i,k,s] =  cosimmrSTAN_in$x_scaled[i,] %*% beta1_ans[s,k,]
      #     }
      #   }
      # }
      #
      # #Then each p is sum f / sum exp f
      #
      # p_sample = array(NA, dim = c(cosimmrSTAN_in$n_obs, n_samples,  cosimmrSTAN_in$n_sources))
      #
      # for(j in 1:n_samples){
      #   for (n_obs in 1:cosimmrSTAN_in$n_obs) {
      #     p_sample[n_obs,j, ] <- exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]) / (sum((exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]))))
      #   }
      # }




    }




  }
    } else if(type == "STAN_MCMC"){


      if(cosimmrSTAN_in$random_effects == TRUE){
        if(error_type == "processxresidual"){
          #This is not nested random effects pxr WORKING
          stan_dat = list(
            J = cosimmrSTAN_in$n_tracers,
            N = cosimmrSTAN_in$n_obs,
            K = cosimmrSTAN_in$n_sources,
            L1 = ncol(cosimmrSTAN_in$X_inner),
            y = cosimmrSTAN_in$mixtures,
            q = cosimmrSTAN_in$concentration_means,
            s_mean = cosimmrSTAN_in$source_means,
            c_mean = cosimmrSTAN_in$correction_means,
            s_sd = cosimmrSTAN_in$source_sds,
            c_sd = cosimmrSTAN_in$correction_sds,
            sigma_shape = prior_control$sigma_shape,
            sigma_rate = prior_control$sigma_rate,
            not_solo = solo,
            X_int = cosimmrSTAN_in$X_intercept,
            X_inner = cosimmrSTAN_in$X_inner,
            omicron_shape = prior_control$omicron_shape,
            omicron_rate = prior_control$omicron_rate
          )


          model = stanmodels$STAN_nested_1_random_effect_pxr_raw

          fit_mcmc <- sampling(
            model,
            data = stan_dat,
            seed = 1,
            iter = 5000,
            cores = 1
          )

          extracted_samples = rstan::extract(fit_mcmc)


          #Want to extract all the betas in a sensible way first I think
          #alpha_ans = extracted_samples$alpha
          beta0_ans = extracted_samples$beta0
          beta1_ans = extracted_samples$beta1 # This is n_samples * K * n_covariates
          sigma_ans = extracted_samples$sigma
          omicron_ans = extracted_samples$omicron
          p_sample = extracted_samples$p






        }else if(error_type == "process+residual"){
          #not nested random effects p+r WORKING
          stan_dat = list(
            J = cosimmrSTAN_in$n_tracers,
            N = cosimmrSTAN_in$n_obs,
            K = cosimmrSTAN_in$n_sources,
            L1 = ncol(cosimmrSTAN_in$X_inner),
            y = cosimmrSTAN_in$mixtures,
            q = cosimmrSTAN_in$concentration_means,
            s_mean = cosimmrSTAN_in$source_means,
            c_mean = cosimmrSTAN_in$correction_means,
            s_sd = cosimmrSTAN_in$source_sds,
            c_sd = cosimmrSTAN_in$correction_sds,
            sigma_shape = prior_control$sigma_shape,
            sigma_rate = prior_control$sigma_rate,
            not_solo = solo,
            X_int = cosimmrSTAN_in$X_intercept,
            X_inner = cosimmrSTAN_in$X_inner
          )
          model = stanmodels$STAN_nested_1_random_effect_pandr

          fit_mcmc <- sampling(
            model,
            data = stan_dat,
            seed = 1,
            iter = 5000,
            cores = 1
          )

          extracted_samples = rstan::extract(fit_mcmc)

          #Want to extract all the betas in a sensible way first I think
          #alpha_ans = extracted_samples$alpha
          beta1_ans = extracted_samples$beta1 # This is n_samples * K * n_covariates
          beta0_ans = extracted_samples$beta0
          sigma_ans = extracted_samples$sigma
          omicron_ans = NULL
          p_sample = extracted_samples$p

          #CONVERT TO F
          #f is x_inner * beta1 + x_outer * beta_2
          #p should be n_ind * n_samples * K
        #   f = array(NA, dim = c(cosimmrSTAN_in$n_obs, cosimmrSTAN_in$n_sources, n_samples))
        #
        #   for(i in 1:cosimmrSTAN_in$n_obs){
        #     for(k in 1: cosimmrSTAN_in$n_sources){
        #       for(s in 1:n_samples){
        #         f[i,k,s] =  cosimmrSTAN_in$X_inner[i,] %*% beta1_ans[s,k,]
        #       }
        #     }
        #   }
        #
        #   #Then each p is sum f / sum exp f
        #
        #   p_sample = array(NA, dim = c(cosimmrSTAN_in$n_obs, n_samples,  cosimmrSTAN_in$n_sources))
        #
        #   for(j in 1:n_samples){
        #     for (n_obs in 1:cosimmrSTAN_in$n_obs) {
        #       p_sample[n_obs,j, ] <- exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]) / (sum((exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]))))
        #     }
        #   }
        }


      }else if(cosimmrSTAN_in$random_effects == FALSE){
        if(error_type == "processxresidual"){
          #This is not nested fixed effects pxr WORKING
          stan_dat = list(
            J = cosimmrSTAN_in$n_tracers,
            N = cosimmrSTAN_in$n_obs,
            K = cosimmrSTAN_in$n_sources,
            L = cosimmrSTAN_in$n_covariates,
            y = cosimmrSTAN_in$mixtures,
            q = cosimmrSTAN_in$concentration_means,
            s_mean = cosimmrSTAN_in$source_means,
            c_mean = cosimmrSTAN_in$correction_means,
            s_sd = cosimmrSTAN_in$source_sds,
            c_sd = cosimmrSTAN_in$correction_sds,
            sigma_shape = prior_control$sigma_shape,
            sigma_rate = prior_control$sigma_rate,
            not_solo = solo,
            X = cosimmrSTAN_in$x_scaled,
            omicron_shape = prior_control$omicron_shape,
            omicron_rate = prior_control$omicron_rate
          )

          model = stanmodels$STAN_VB_pxr_raw
          fit_mcmc <- sampling(
            model,
            data = stan_dat,
            seed = 1,
            iter = 5000,
            cores = 1
          )

          extracted_samples = rstan::extract(fit_mcmc)

          #Want to extract all the betas in a sensible way first I think
          #alpha_ans = extracted_samples$alpha
          beta1_ans = extracted_samples$beta # This is n_samples * K * n_covariates
          beta0_ans = NULL
          sigma_ans = extracted_samples$sigma
          omicron_ans = extracted_samples$omicron
          p_sample = extracted_samples$p

          #CONVERT TO F
          #f is x_inner * beta1 + x_outer * beta_2
          #p should be n_ind * n_samples * K
          # f = array(NA, dim = c(cosimmrSTAN_in$n_obs, cosimmrSTAN_in$n_sources, n_samples))
          #
          # for(i in 1:cosimmrSTAN_in$n_obs){
          #   for(k in 1: cosimmrSTAN_in$n_sources){
          #     for(s in 1:n_samples){
          #       f[i,k,s] =  cosimmrSTAN_in$x_scaled[i,] %*% beta1_ans[s,k,]
          #     }
          #   }
          # }
          #
          # #Then each p is sum f / sum exp f
          #
          # p_sample = array(NA, dim = c(cosimmrSTAN_in$n_obs, n_samples,  cosimmrSTAN_in$n_sources))
          #
          # for(j in 1:n_samples){
          #   for (n_obs in 1:cosimmrSTAN_in$n_obs) {
          #     p_sample[n_obs,j, ] <- exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]) / (sum((exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]))))
          #   }
          # }



        }else if(error_type == "process+residual"){
          #not nested fixed effects p+r WORKING

          stan_dat = list(
            J = cosimmrSTAN_in$n_tracers,
            N = cosimmrSTAN_in$n_obs,
            K = cosimmrSTAN_in$n_sources,
            L = cosimmrSTAN_in$n_covariates,
            y = cosimmrSTAN_in$mixtures,
            q = cosimmrSTAN_in$concentration_means,
            s_mean = cosimmrSTAN_in$source_means,
            c_mean = cosimmrSTAN_in$correction_means,
            s_sd = cosimmrSTAN_in$source_sds,
            c_sd = cosimmrSTAN_in$correction_sds,
            sigma_shape = prior_control$sigma_shape,
            sigma_rate = prior_control$sigma_rate,
            not_solo = solo,
            X = cosimmrSTAN_in$x_scaled
          )

          model = stanmodels$STAN_VB_pandr_raw


          fit_mcmc <- sampling(
            model,
            data = stan_dat,
            seed = 1,
            iter = 5000,
            cores = 1
          )

          extracted_samples = rstan::extract(fit_mcmc)
          #Want to extract all the betas in a sensible way first I think
          #alpha_ans = extracted_samples$alpha
          beta1_ans = extracted_samples$beta # This is n_samples * K * n_covariates
          beta0_ans = NULL
          sigma_ans = extracted_samples$sigma
          omicron_ans = NULL
          p_sample = extracted_samples$p

          #CONVERT TO F
          #f is x_inner * beta1 + x_outer * beta_2
          #p should be n_ind * n_samples * K
          # f = array(NA, dim = c(cosimmrSTAN_in$n_obs, cosimmrSTAN_in$n_sources, n_samples))
          #
          # for(i in 1:cosimmrSTAN_in$n_obs){
          #   for(k in 1: cosimmrSTAN_in$n_sources){
          #     for(s in 1:n_samples){
          #       f[i,k,s] =  cosimmrSTAN_in$x_scaled[i,] %*% beta1_ans[s,k,]
          #     }
          #   }
          # }
          #
          # #Then each p is sum f / sum exp f
          #
          # p_sample = array(NA, dim = c(cosimmrSTAN_in$n_obs, n_samples,  cosimmrSTAN_in$n_sources))
          #
          # for(j in 1:n_samples){
          #   for (n_obs in 1:cosimmrSTAN_in$n_obs) {
          #     p_sample[n_obs,j, ] <- exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]) / (sum((exp(f[n_obs,1: cosimmrSTAN_in$n_sources, j]))))
          #   }
          # }
          #



        }




   ### COPY VB WHEN WORKING AND JUST CHANGE TO MCMC


  }
}
  mylist <- list(
    source_names = cosimmrSTAN_in$source_names,
    p = p_sample,
    beta1_ans = beta1_ans,
    sigma_ans = sigma_ans,
    omicron_ans = omicron_ans
  )


  output_all <- list(input = cosimmrSTAN_in, output = mylist)

  class(output_all) <- c("cosimmrSTAN_output")

  return(output_all)



}

