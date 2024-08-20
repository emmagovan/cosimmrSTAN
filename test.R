#I want to test the cosimmrSTAN_load function


#Setup 1 - do the same as cosimmr basically
formula = alligator_data$mixtures ~ alligator_data$long

out<-cosimmrSTAN_load(formula,
                 source_names = alligator_data$source_names,
                 source_means = alligator_data$source_means,
                 source_sds = alligator_data$source_sds,
                 source = NULL,
                 n_each_source = c(5,10),
                 correction_means = alligator_data$TEF_means,
                 correction_sds = alligator_data$TEF_sds,
                 concentration_means = NULL,
                 scale_x = FALSE,
                 raw_source = FALSE,
                 random_effects = FALSE,
                 hierarchical_fitting = TRUE,
                 shape_sig = 1)


##Setup 2 random effect
formula = alligator_data$mixtures ~ 1|as.factor(alligator_data$sex)

out<-cosimmrSTAN_load(formula,
                      source_names = alligator_data$source_names,
                      source_means = alligator_data$source_means,
                      source_sds = alligator_data$source_sds,
                      source = NULL,
                      n_each_source = c(5,10),
                      correction_means = alligator_data$TEF_means,
                      correction_sds = alligator_data$TEF_sds,
                      concentration_means = NULL,
                      scale_x = FALSE,
                      raw_source = FALSE,
                      random_effects = TRUE,
                      hierarchical_fitting = TRUE,
                      shape_sig = 1)


## Nested
pack = as.factor(c(rep(1,21), rep(2,50), rep(3, 70), rep(4,40)))
region = as.factor(c(rep(1, 71), rep(2,110)))
formula2 = alligator_data$mixtures ~ 1|region/pack

out_nest<-cosimmrSTAN_load(formula2,
                      source_names = alligator_data$source_names,
                      source_means = alligator_data$source_means,
                      source_sds = alligator_data$source_sds,
                      source = NULL,
                      n_each_source = c(5,10),
                      correction_means = alligator_data$TEF_means,
                      correction_sds = alligator_data$TEF_sds,
                      concentration_means = NULL,
                      scale_x = FALSE,
                      raw_source = FALSE,
                      random_effects = TRUE,
                      hierarchical_fitting = TRUE,
                      shape_sig = 1)



##SETUP wolves
load("~/Documents/GitHub/testing_STAN_cosimmr/data/wolves_data.rda")
wolves_consumer = wolves[[1]]
wolves_discrimination = wolves[[2]]
wolves_sources = wolves[[3]]
y = as.matrix(wolves_consumer[,c(1:2)])
pack = as.factor(wolves_consumer$Pack)
region = as.factor(wolves_consumer$Region)
formula = y ~ (1|region/pack) -1
wolves_sources = wolves[[3]][c(1,4,7),] #Small for now - need to accept different sources at some point
q = matrix(rep(1, 6), ncol = 2)
s_mean = as.matrix(wolves_sources[,c(3,5)])
s_sd = as.matrix(wolves_sources[,c(4,6)])
c_mean = as.matrix(wolves_discrimination[,c(2,4)])
c_sd = as.matrix(wolves_discrimination[,c(3,5)])

out_wolves = cosimmrSTAN_load(formula,
                              source_names = wolves_sources$...1,
                              source_means = s_mean,
                              source_sds = s_sd,
                              n_each_source = wolves_sources$n,
                              correction_means = c_mean,
                              correction_sds = c_sd,
                              concentration_means = NULL,
                              scale_x = FALSE,
                              raw_source = FALSE,
                              random_effects = TRUE,
                              hierarchical_fitting = TRUE

)

a = lme4::glFormula(formula)

X_inner = t(as.matrix(a$reTrms$Zt))[,1:a$reTrms$nl[1]]
X_region = t(as.matrix(a$reTrms$Zt))[,(a$reTrms$nl[1] +1):(a$reTrms$nl[1] + a$reTrms$nl[2])]

stan_dat = list(
  J = 2,
  N = nrow(y),
  K = nrow(s_mean),
  L1 = ncol(X_inner),
  L2 = ncol(X_region),
  y = y,
  q = q,
  s_mean = s_mean,
  c_mean = c_mean,
  s_sd = s_sd,
  c_sd = c_sd,
  sigma_shape = c(1,1),
  sigma_rate = c(1,1),
  not_solo = 1,
  X_intercept = a$X,
  X_inner = X_inner,
  X_outer = X_region,
  region_for_pack = region_for_pack_output,
  omicron_shape = c(1,1),
  omicron_rate = c(1,1)
)


fit <- stan(
  file = "inst/stan/STAN_nested_2_random_effect_pxr_raw.stan",
  data = stan_dat,
  init = function() list(
    sigma_raw = rep(1, J),
    sigma_pack = matrix(1, K, L1),
    sigma_region = rep(1, K),
    mu_region = rep(1, K),
    omicron = rep(1, J),
    beta1 = matrix(1, K, L1),
    beta2 = matrix(1, K, L2)
  ),
  chains = 1,
  iter = 2000,
  warmup = 1000,
  verbose = TRUE
)


fit <- sampling(
  object = stanmodels$STAN_nested_2_random_effect_pxr_raw,     # Use the precompiled model
  data = stan_dat,
  init = function() list(
    sigma_raw = rep(1, J),
    sigma_pack = matrix(1, K, L1),
    sigma_region = rep(1, K),
    mu_region = rep(1, K),
    omicron = rep(1, J),
    beta1 = matrix(1, K, L1),
    beta2 = matrix(1, K, L2)
  ),
  chains = 1,
  iter = 2000,
  warmup = 1000,
  verbose = TRUE
)











