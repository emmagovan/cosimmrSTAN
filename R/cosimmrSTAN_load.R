#' Function to load in cosimmrSTAN data and check for errors
#'
#' This function takes in the mixture data, food source means and standard
#' deviations, and (optionally) correction factor means and standard
#' deviations, and concentration proportions. It performs some (non-exhaustive)
#' checking of the data to make sure it will run through cosimmrSTAN. It outputs an
#' object of class \code{cosimmrSTAN_input}.
#'
#' For standard stable isotope mixture modelling, the mixture matrix will
#' contain a row for each individual and a column for each isotopic value.
#' \code{cosimmrSTAN} will allow for any number of isotopes and any number of
#' observations, within computational limits. The source means/sds should be
#' provided for each food source on each isotope. The correction means (usually
#' trophic enrichment factors) can be set as zero if required, and should be of
#' the same shape as the source values. The concentration dependence means
#' should be estimated values of the proportion of each element in the food
#' source in question and should be given in proportion format between 0 and 1.
#' At present there is no means to include concentration standard deviations.
#'
#' @param formula Formula giving in form y ~ x where y is a vector or matrix
#' of mixture values and x is a vector or matrix of covariates
#' @param source_names The names of the sources given as a character string
#' @param source_means The means of the source values, given as a matrix where
#' the number of rows is the number of sources and the number of columns is the
#' number of tracers
#' @param source_sds The standard deviations of the source values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers
#' @param source The raw source data if you wish to use raw source data.
#' @param n_each_source The number of sources collected. Needed for hierarchical
#' fitting
#' @param correction_means The means of the correction values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers. If not provided these are set to 0.
#' @param correction_sds The standard deviations of the correction values,
#' given as a matrix where the number of rows is the number of sources and the
#' number of columns is the number of tracers. If not provided these are set to
#' 0.
#' @param concentration_means The means of the concentration values, given as a
#' matrix where the number of rows is the number of sources and the number of
#' columns is the number of tracers. These should be between 0 and 1. If not
#' provided these are all set to 1.
#' @param scale_x Whether or not you wish to scale the x values provided, or run
#' the model using the original x values. Defaults to TRUE.
#' @param raw_source If you are using raw source data or means. Defaults to FALSE.
#' @param hierarchical_fitting Whether to hierarchically fit the source data. Defaults
#' to FALSE. Requires n_sources to be provided.
#' @param shape_sig Shape parameter for fitting raw source data. Defaults to 1.
#'
#' @import checkmate
#'
#'
#' @return An object of class \code{cosimmrSTAN_input} with the following elements:
#' \item{mixtures }{The mixture data} \item{source_names }{Source means}
#' \item{sources_sds }{Source standard deviations} \item{correction_means
#' }{Correction means} \item{correction_sds }{Correction standard deviations}
#' \item{concentration_means }{Concentration dependence means} \item{n_obs
#' }{The number of observations} \item{n_tracers }{The number of
#' tracers/isotopes} \item{n_sources }{The number of sources} \item{n_groups
#' }{The number of groups}
#' @author Emma Govan <emmagovan@@gmail.com>, Andrew Parnell
#' @seealso See \code{\link{cosimmr_stan}} for complete examples.
#' @examples
#' \donttest{
#'
#' # A simple example with 10 observations, 2 tracers and 4 sources
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
#'
#'
#' print(cosimmrSTAN_1)
#' }
#' @export cosimmrSTAN_load
cosimmrSTAN_load <- function(formula,
source_names,
source_means = NULL,
source_sds = NULL,
source = NULL,
n_each_source = NULL,
correction_means = NULL,
correction_sds = NULL,
concentration_means = NULL,
scale_x = TRUE,
raw_source = FALSE,
hierarchical_fitting = FALSE,
shape_sig = 1) {
  # Function to load in data for cosimmrSTAN
formula = formula
  #This just checks whether theres a | in the formula or not
  r_e <- lme4::findbars(formula)

  if (length(r_e) == 0) {
    random_effects = FALSE
  } else {
    random_effects = TRUE
  }

  if(random_effects == TRUE){

    a = lme4::glFormula(formula)
    X_fixed = a$X #This extracts the fixed parts

    X_random = t(as.matrix(a$reTrms$Zt))
    re_names_order = names(a$reTrms$flist)
    re_levels <- unlist(lapply(a$reTrms$flist, nlevels))

    if(scale_x == TRUE){
      if(stats::sd(X_fixed[,1]) == 0){
        if(ncol(X_fixed) == 1){
          intercept = TRUE
          scaled_mat = X_fixed
          x_scaled = scaled_mat
          colnames(x_scaled) = c(colnames(X_fixed))

          scaled_center = attr(scaled_mat, "scaled:center")

          scaled_scale = attr(scaled_mat, "scaled:scale")
          message("Cannot scale intercept")

        } else if(ncol(X_fixed) != 1){
          # Original code
          intercept = TRUE
          scaled_mat = scale(X_fixed[,(2:ncol(X_fixed))])
          x_scaled = cbind(X_fixed[,1],
                           scaled_mat)
          colnames(x_scaled) = c(colnames(X_fixed))

          scaled_center = attr(scaled_mat, "scaled:center")

          scaled_scale = attr(scaled_mat, "scaled:scale")

        }

      }else if(stats::sd(X_fixed[,1]) != 0){
        intercept = FALSE
        scaled_mat = scale(X_fixed)
        x_scaled = scaled_mat
        colnames(x_scaled) = c(colnames(X_fixed))

        scaled_center = attr(scaled_mat, "scaled:center")

        scaled_scale = attr(scaled_mat, "scaled:scale")

      }
    }else if(scale_x == FALSE){

      if(stats::sd(X_fixed[,1]) == 0){
        intercept = TRUE}else(intercept = FALSE)
      x_scaled = X_fixed
      scaled_center = NULL

      scaled_scale = NULL


    }
    # Write a function that generically tests for any 2D numeric data shape such as matrix, data frame or tibble
    assert_2D_numeric <- function(x,
                                  nrows = NULL,
                                  ncols = NULL,
                                  null.ok = FALSE) {
      assert(
        test_data_frame(x,
                        types = c("double", "numeric"),
                        nrows = nrows,
                        ncols = ncols,
                        null.ok = null.ok
        ),
        test_matrix(x,
                    mode = "numeric",
                    nrows = nrows,
                    ncols = ncols,
                    null.ok = null.ok
        ),
        test_tibble(x,
                    types = c("double", "numeric"),
                    nrows = nrows,
                    ncols = ncols,
                    null.ok = null.ok
        )
      )
    }

    original_x = a$fr[,(2):ncol(a$fr)]

    cnames = colnames(a$fr)[-1]
    covariates = data.frame(a$fr[-1])
    colnames(covariates) = cnames

    mixtures = as.matrix(a$fr[1])
    n_tracers = ncol(mixtures)

    # Mixtures must be a matrix - the number of rows is the number of observations and the number of columns is the number of tracers
    # assert_matrix(mixtures)
    assert_2D_numeric(mixtures)
    n_obs <- nrow(mixtures)
    n_tracers <- ncol(mixtures)

    # Add column names if they're not there
    if (is.null(colnames(mixtures))) {
      colnames(mixtures) <- paste0("tracer", 1:n_tracers)
    }

    if(is.null(x_scaled) == TRUE){

    }else{
      if (is.null(colnames(x_scaled))) {
        colnames(x_scaled) <- paste0("covariate", 1:n_tracers)
      }

    }



  }else if(random_effects == FALSE){
    re_names_order = NULL
    re_levels = FALSE
    X_intercept = NULL
    X_random = NULL
    X_fixed = NULL
    #Possibly need covariate data frame
    #c_df = data.frame(colour_cat = colour_cat, letter_cat = letter_cat, numeric_cov = numeric_cov)
    # Go through each object and check that it matches the requirements
    mixtures = as.matrix(stats::model.frame(formula)[,1])
    original_x = stats::model.matrix(formula)

    #Need to add some method here to keep the column names
    cnames = colnames(stats::model.frame(formula))[-1]
    covariates = data.frame((stats::model.frame(formula)[,-1]))
    colnames(covariates) = cnames




    if(nrow(mixtures) == 1){
      #This is if its just 1 entry
      x_scaled = stats::model.matrix(formula)
    }
    if(scale_x == TRUE){
      if(stats::sd(stats::model.matrix(formula)[,1]) == 0){
        if(ncol(stats::model.matrix(formula)) == 1){
          intercept = TRUE
          scaled_mat = (stats::model.matrix(formula))
          x_scaled = scaled_mat
          colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))

          scaled_center = attr(scaled_mat, "scaled:center")

          scaled_scale = attr(scaled_mat, "scaled:scale")
          message("Cannot scale when using mixtures ~1")
        } else if(ncol(stats::model.matrix(formula)) != 1){
          # Original code
          intercept = TRUE
          scaled_mat = scale(stats::model.matrix(formula)[,(2:ncol(stats::model.matrix(formula)))])
          x_scaled = cbind(stats::model.matrix(formula)[,1],
                           scaled_mat)
          colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))

          scaled_center = attr(scaled_mat, "scaled:center")

          scaled_scale = attr(scaled_mat, "scaled:scale")

        }

      }else if(stats::sd(stats::model.matrix(formula)[,1]) != 0){
        intercept = FALSE
        scaled_mat = scale(stats::model.matrix(formula))
        x_scaled = scaled_mat
        colnames(x_scaled) = c(colnames(stats::model.matrix(formula)))

        scaled_center = attr(scaled_mat, "scaled:center")

        scaled_scale = attr(scaled_mat, "scaled:scale")

      }
    } else if(scale_x == FALSE){

      if(stats::sd(stats::model.matrix(formula)[,1]) == 0){
        intercept = TRUE}else(intercept = FALSE)
      x_scaled = stats::model.matrix(formula)
      scaled_center = NULL

      scaled_scale = NULL


    }

    # Write a function that generically tests for any 2D numeric data shape such as matrix, data frame or tibble
    assert_2D_numeric <- function(x,
                                  nrows = NULL,
                                  ncols = NULL,
                                  null.ok = FALSE) {
      assert(
        test_data_frame(x,
                        types = c("double", "numeric"),
                        nrows = nrows,
                        ncols = ncols,
                        null.ok = null.ok
        ),
        test_matrix(x,
                    mode = "numeric",
                    nrows = nrows,
                    ncols = ncols,
                    null.ok = null.ok
        ),
        test_tibble(x,
                    types = c("double", "numeric"),
                    nrows = nrows,
                    ncols = ncols,
                    null.ok = null.ok
        )
      )
    }

    # Mixtures must be a matrix - the number of rows is the number of observations and the number of columns is the number of tracers
    # assert_matrix(mixtures)
    assert_2D_numeric(mixtures)
    n_obs <- nrow(mixtures)
    n_tracers <- ncol(mixtures)

    # Add column names if they're not there
    if (is.null(colnames(mixtures))) {
      colnames(mixtures) <- paste0("tracer", 1:n_tracers)
    }

    if(is.null(x_scaled) == TRUE){

    }else{
    if (is.null(colnames(x_scaled))) {
      colnames(x_scaled) <- paste0("covariate", 1:n_tracers)
    }

    }
  }
  assert_2D_numeric <- function(x,
                                nrows = NULL,
                                ncols = NULL,
                                null.ok = FALSE) {
    assert(
      test_data_frame(x,
                      types = c("double", "numeric"),
                      nrows = nrows,
                      ncols = ncols,
                      null.ok = null.ok
      ),
      test_matrix(x,
                  mode = "numeric",
                  nrows = nrows,
                  ncols = ncols,
                  null.ok = null.ok
      ),
      test_tibble(x,
                  types = c("double", "numeric"),
                  nrows = nrows,
                  ncols = ncols,
                  null.ok = null.ok
      )
    )
  }

  # source_names must be a character vector - the length of it is the number of sources
  assert_character(source_names)
  n_sources <- length(source_names)
  n_obs <- nrow(mixtures)
  n_tracers <- ncol(mixtures)



  # assert_matrix(source_sds, nrows = n_sources, ncols = n_tracers)
  assert_2D_numeric(correction_means,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = ifelse(is.null(correction_sds),
                                     TRUE, FALSE
                    )
  )

  # assert_matrix(correction_means,
  #   nrows = n_sources,
  #   ncols = n_tracers,
  #   null.ok = ifelse(is.null(correction_sds),
  #     TRUE, FALSE
  #   )
  # )
  assert_2D_numeric(correction_sds,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = ifelse(is.null(correction_sds),
                                     TRUE, FALSE
                    )
  )
  # assert_matrix(correction_sds,
  #   nrows = n_sources,
  #   ncols = n_tracers,
  #   null.ok = ifelse(is.null(correction_means),
  #     TRUE, FALSE
  #   )
  # )
  assert_2D_numeric(concentration_means,
                    nrows = n_sources,
                    ncols = n_tracers,
                    null.ok = TRUE
  )
  # assert_matrix(concentration_means,
  #   nrows = n_sources,
  #   ncols = n_tracers, null.ok = TRUE
  # )

  # Fill in correction means
  if (is.null(correction_means)) {
    correction_means <- matrix(0, ncol = n_tracers, nrow = n_sources)
    correction_sds <- matrix(0, ncol = n_tracers, nrow = n_sources)
  }

  # concentration_means must be a matrix where all elements are less than 1
  if (is.null(concentration_means)) {
    concentration_means <- matrix(1, ncol = n_tracers, nrow = n_sources)
  } else {
    assert_true(all(concentration_means < 1) & all(concentration_means > 0))
  }

  # Check the groups are the right length and structure if given


  ####### USE RAW DATA OR HIERARCHICAL HERE???????

  if(raw_source == TRUE){
    stan_dat = list(
      J = (ncol(source)-1),
      K = length(unique(source$Source)),
      N = length(source$Source),
      Y = source[,2:ncol(source)],
      source = as.numeric(factor(source$Source, levels = unique(source$Source))),
      shape_sig = shape_sig
    )

    fit_mcmc<- sampling(
      stanmodels$raw_source,
      data = stan_dat,
      seed = 1,
      iter = 1000,
      cores = 1)

    #extracted_mcmc = extract(fit_mcmc)

    extracted_samples = rstan::extract(fit_mcmc)

    #Want to extract all the betas in a sensible way first I think
    colnames = colnames(source)[-1]
    mean_names <- paste0("mean_", colnames)
    sd_names <- paste0("sd_", colnames)
    mu_out = extracted_samples$mu_jk # This is n_samples * K * n_covariates


    Sigma_out = extracted_samples$Sigma_k



    source_means_out = t(apply(mu_out, c(2,3), mean))
    colnames(source_means_out) = mean_names
    source_sds_out = apply(Sigma_out, c(2,3), mean)
    colnames(source_sds_out) = sd_names



  } else if(raw_source == FALSE){
    if(hierarchical_fitting == TRUE){
    # source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
    assert_2D_numeric(source_means,
                      nrows = n_sources,
                      ncols = n_tracers)

    # assert_matrix(source_means, nrows = n_sources, ncols = n_tracers)
    assert_2D_numeric(source_sds,
                      nrows = n_sources,
                      ncols = n_tracers)

    stan_dat = list(
      J = ncol(source_means),
      K = length(source_names),
      n = n_each_source,
      source_mean = source_means,
      source_sd = source_sds
    )


model = stanmodels$Hierarchical
    #rstan::rstan_options(auto_write = TRUE)
    fit_mcmc<- rstan::sampling(
      model,
      data = stan_dat,
      seed = 1,
      iter = 1000,
      cores = 1
    )

    extracted_mcmc = rstan::extract(fit_mcmc)

    extracted_samples = rstan::extract(fit_mcmc)

    #Want to extract all the betas in a sensible way first I think
    mu_out = extracted_samples$mu # This is n_samples * K * n_covariates
    Sigma_out = extracted_samples$sigma


  source_means_out = apply(mu_out, c(2,3), mean)
  colnames(source_means_out) = colnames(source_means)

  source_sds_out = apply(Sigma_out, c(2,3), mean)
  colnames(source_sds_out) = colnames(source_sds)
    ##Need to check if these need to be edited before they can be used??

  #  stanmodels$Hierarchical
    } else{
      # source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
      assert_2D_numeric(source_means,
                        nrows = n_sources,
                        ncols = n_tracers)

      # assert_matrix(source_means, nrows = n_sources, ncols = n_tracers)
      assert_2D_numeric(source_sds,
                        nrows = n_sources,
                        ncols = n_tracers)

      source_means_out = source_means
      source_sds_out = source_sds
}
  }









  # Prepare output and give class
  out <- list(
    mixtures = mixtures,
    x_scaled = x_scaled,
    X_fixed = X_fixed,
    X_random = X_random,
    source_names = source_names,
    source_means = source_means_out,
    source_sds = source_sds_out,
    correction_means = correction_means,
    correction_sds = correction_sds,
    concentration_means = concentration_means,
    n_obs = n_obs,
    n_tracers = n_tracers,
    n_sources = n_sources,
    scale_x = scale_x,
    scaled_center = scaled_center,
    scaled_scale = scaled_scale,
    intercept = intercept,
    covariates_df = covariates,
    n_covariates =  ncol(x_scaled),
    original_x = original_x,
    random_effects = random_effects,
    hierarchical_fitting = hierarchical_fitting,
    formula = deparse(formula),
    re_names_order = re_names_order,
    re_levels = re_levels

  )

  # Look through to see whether there are any missing values in anything
  if (any(unlist(lapply(out, "is.na")))) {
    warning("Missing values provided for some values. Check your inputs")
  }

  class(out) <- "cosimmrSTAN_input"

  return(out)
}

