#I want to test the cosimmrSTAN_load function


#Setup 1 - do the same as cosimmr basically--------
Length = alligator_data$length
formula = alligator_data$mixtures ~ Length
source_names = alligator_data$source_names
source_means = alligator_data$source_means
source_sds = alligator_data$source_sds
source = NULL
n_each_source = c(5,10)
correction_means = alligator_data$TEF_means
correction_sds = alligator_data$TEF_sds
concentration_means = NULL
scale_x = FALSE
raw_source = FALSE
random_effects = FALSE
hierarchical_fitting = TRUE
shape_sig = 1

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
                 hierarchical_fitting = FALSE,
                 shape_sig = 1)

out_1 = cosimmr_stan(out, error_type = "processxresidual")


out2 = cosimmr_stan(out, type = "STAN_MCMC", error_type = "process+residual")


##Setup random effect----------------
formula = alligator_data$mixtures ~ 1|as.factor(alligator_data$sex)

in_random<-cosimmrSTAN_load(formula,
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


out_random = cosimmr_stan(out, error_type = "process+residual")


out_mcmc = cosimmr_stan(out, type = "STAN_MCMC", error_type = "process+residual")
out_1 = cosimmr_stan(out)


out2 = cosimmr_stan(out, type = "STAN_MCMC")





##SETUP wolves--------
load("~/Documents/GitHub/testing_STAN_cosimmr/data/wolves_data.rda")
wolves_consumer = wolves[[1]]
wolves_discrimination = wolves[[2]]
wolves_sources = wolves[[3]]
y = as.matrix(wolves_consumer[,c(1:2)])
pack = as.factor(wolves_consumer$Pack)
region = as.factor(wolves_consumer$Region)
formula = y ~ (1|pack)
wolves_sources = wolves[[3]][c(1,4,7),] #Small for now - need to accept different sources at some point
q = matrix(rep(1, 6), ncol = 2)
s_mean = as.matrix(wolves_sources[,c(3,5)])
s_sd = as.matrix(wolves_sources[,c(4,6)])
c_mean = as.matrix(wolves_discrimination[,c(2,4)])
c_sd = as.matrix(wolves_discrimination[,c(3,5)])

in_wolves = cosimmrSTAN_load(formula,
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

out_wolves = cosimmr_stan(in_wolves)

### Geese example------------------------


# Get the data
library(readxl)
path <- system.file("extdata", "geese_data.xls", package = "simmr")
geese_data <- lapply(excel_sheets(path), read_excel, path = path)
targets <- geese_data[[1]] #[1:20,]
sources <- geese_data[[2]]
TEFs <- geese_data[[3]]
concdep <- geese_data[[4]]
N <- nrow(targets)
K <- nrow(sources)
J <- 2
L <- length(unique(targets$Time))
Time = as.factor(targets$Time)
y = as.matrix(targets[,1:2]) # Data matri
q = concdep[,2:3] # cond dep matri
s_mean = sources[,2:3] # source mean
c_mean = TEFs[,2:3] # corr mean
s_sd = sources[,4:5] # s_sd matri
c_sd = TEFs[,4:5] # c_sd matri
sigma_shape = rep(1, 2)
sigma_rate = rep(1, 2)

formula = y ~ (1|Time)

source_means = s_mean
source_sds = s_sd
source = NULL
correction_means = c_mean
correction_sds = c_sd
concentration_means = q
scale_x = FALSE
raw_source = FALSE
random_effects = TRUE
hierarchical_fitting = FALSE
shape_sig = 1

in_geese<-cosimmrSTAN_load(formula,
                            source_names =sources$Sources,
                            source_means = s_mean,
                            source_sds = s_sd,
                            source = NULL,
                            correction_means = c_mean,
                            correction_sds = c_sd,
                            concentration_means = q,
                            scale_x = FALSE,
                            raw_source = FALSE,
                            random_effects = TRUE,
                            hierarchical_fitting = FALSE,
                            shape_sig = 1)

out_geese = cosimmr_stan(in_geese, error_type = "process+residual")

Weight = targets$`Net Wt`
formula2 = y ~ 1|Weight
in_geese2<-cosimmrSTAN_load(formula2,
                           source_names =sources$Sources,
                           source_means = s_mean,
                           source_sds = s_sd,
                           source = NULL,
                           correction_means = c_mean,
                           correction_sds = c_sd,
                           concentration_means = q,
                           scale_x = FALSE,
                           raw_source = FALSE,
                           random_effects = FALSE,
                           hierarchical_fitting = FALSE,
                           shape_sig = 1)

out_geese3 = cosimmr_stan(in_geese2, error_type = "process+residual")




