library(bhmbasket)
library(coda)

############################################################
# 1. Create trial
############################################################

trial_one <- createTrial(
  n_subjects = c(30,30,30),
  means      = c(1.0, 1.5, 2.0),
  sds        = c(0.5, 0.5, 0.5),
  endpoint   = "normal"
)
scenario <- trial_one[[1]]

############################################################
# 2. Correct prior parameters for your JAGS model
############################################################

priors <- list(
  mu_pop_mean = 0,
  prec_mu_pop = 1/10^2,
  
  tau_shape   = 0.001,
  tau_rate    = 0.001,
  
  sigma_shape = 0.001,
  sigma_rate  = 0.001
)

prep <- prepareAnalysis(
  method_name      = "normal",
  prior_parameters = priors
)

############################################################
# 3. Build correct JAGS data for THIS model
############################################################

j_data <- prep$j_data
j_data$y <- scenario$y[1, ]
j_data$n <- scenario$n_subjects[1, ]
j_data$J <- length(j_data$y)

############################################################
# 4. Get posterior samples
############################################################

samples <- getPosteriors(
  j_parameters      = prep$j_parameters,
  j_model_file      = prep$j_model_file,
  j_data            = j_data,
  n_mcmc_iterations = 200
)


samples <- as.matrix(samples)

# Two chains: getPosteriors includes both chains interleaved
chain1 <- samples[seq(1, nrow(samples), 2), ]
chain2 <- samples[seq(2, nrow(samples), 2), ]

mcmc_list <- mcmc.list(mcmc(chain1), mcmc(chain2))

############################################################
# 4. Supervisor-requested diagnostics
############################################################

# A) Traceplot
plot(mcmc_list[, "theta_1"], main = "Traceplot for theta_1 (two chains)")

# B) Cumulative mean
theta <- samples[, "theta_1"]
cma <- cumsum(theta) / seq_along(theta)
plot(cma, type="l", main="Cumulative Mean of theta_1")

# C) Cumulative 97.5% quantile
cum_q975 <- sapply(seq_along(theta), function(k)
  quantile(theta[1:k], 0.975)
)
plot(cum_q975, type="l", main="Cumulative 97.5% Quantile of theta_1")

# D) Gelman–Rubin Rhat
gelman.diag(mcmc_list)
