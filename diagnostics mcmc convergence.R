library(bhmbasket)
library(coda)

############################################################
# 1. Create one continuous trial
############################################################

trial_one <- createTrial(
  n_subjects = c(30, 30, 30),
  means      = c(1.0, 1.5, 2.0),
  sds        = c(0.5, 0.5, 0.5),
  endpoint   = "normal"
)

scenario <- trial_one[[1]]

############################################################
# 2. Prior parameters for YOUR current normal model
############################################################

prior_parameters_list <- getPriorParameters(
  endpoint     = "normal",
  method_names = "normal",
  target_means = c(1.0, 1.5, 2.0),
  n_worth      = 1,
  tau_scale    = 1,
  w_j          = 0.5,
  sigma_shape  = 1,
  sigma_rate   = 1
)

prep <- prepareAnalysis(
  method_name      = "normal",
  prior_parameters = prior_parameters_list[["normal"]]
)

############################################################
# 3. Build JAGS data for one trial realization
############################################################

j_data <- prep$j_data
j_data$y <- scenario$y[1, ]
j_data$n <- scenario$n_subjects[1, ]

############################################################
# 4. Run posterior sampling
############################################################

samples <- getPosteriors(
  j_parameters      = prep$j_parameters,
  j_model_file      = prep$j_model_file,
  j_data            = j_data,
  n_mcmc_iterations = 2000
)

samples <- as.matrix(samples)

colnames(samples)


############################################################
# 5. Split into two chains
############################################################

chain1 <- samples[seq(1, nrow(samples), 2), , drop = FALSE]
chain2 <- samples[seq(2, nrow(samples), 2), , drop = FALSE]

mcmc_list <- mcmc.list(
  mcmc(chain1),
  mcmc(chain2)
)

############################################################
# 6A. Traceplot
############################################################

traceplot(mcmc_list[, "theta_1"])

############################################################
# 6B. Cumulative mean
############################################################

theta <- samples[, "theta_1"]
cma <- cumsum(theta) / seq_along(theta)

plot(cma, type = "l", main = "Cumulative Mean of theta_1")

############################################################
# 6C. Cumulative 97.5% quantile
############################################################

cum_q975 <- sapply(seq_along(theta), function(k) {
  quantile(theta[1:k], 0.975)
})

plot(cum_q975, type = "l", main = "Cumulative 97.5% Quantile of theta_1")

############################################################
# 6D. Gelman-Rubin Rhat
############################################################

gelman.diag(mcmc_list[, c("theta_1", "theta_2", "theta_3")])
traceplot(mcmc_list[, "mu"])
traceplot(mcmc_list[, "tau"])
