library(bhmbasket)
library(future)

################################################################################
# PART 2 — MONTE CARLO STABILITY (GO PROBABILITY)
################################################################################

set.seed(123)
Nsim <- 1000

################################################################################
# 1. Simulate scenarios
################################################################################

scenario_list <- simulateScenarios(
  n_subjects_list = list(rep(30, 3)),
  means_list      = list(c(0.4, 0.5, 0.6)),
  sds_list        = list(rep(2.5, 3)),
  n_trials        = Nsim,
  endpoint        = "normal"
)

################################################################################
# 2. Prior for current normal model
################################################################################

prior_parameters_list <- getPriorParameters(
  endpoint     = "normal",
  method_names = "normal",
  target_means = c(0, 0, 0),   # use 0 if you want theta_j to equal actual means
  n_worth      = 1,
  tau_scale    = 1,
  w_j          = 0.5,
  sigma_shape  = 1,
  sigma_rate   = 1
)

################################################################################
# 3. Run analysis
################################################################################

plan(multisession)

ana <- performAnalyses(
  scenario_list         = scenario_list,
  method_names          = "normal",
  prior_parameters_list = prior_parameters_list,
  n_mcmc_iterations     = 3000,
  verbose               = FALSE
)

################################################################################
# 4. Apply go-decision rule
################################################################################
# With evidence_levels = 0.95, x[j] corresponds to the lower 5% posterior quantile.
# So x[j] > c means P(theta_j > c) > 0.95.

decision_rules <- quote(c(
  x[1] > 0.1,
  x[2] > 0.2,
  x[3] > 0.3
))

go_list <- getGoDecisions(
  analyses_list   = ana,
  cohort_names    = c("theta_1", "theta_2", "theta_3"),
  evidence_levels = c(0.95, 0.95, 0.95),
  boundary_rules  = decision_rules,
  overall_min_gos = 1
)

go_matrix <- go_list$scenario_1$decisions_list$normal
go_vec    <- go_matrix[, "overall"]

################################################################################
# 5. Moving average of GO probability
################################################################################

mov <- cumsum(go_vec) / seq_along(go_vec)

plot(
  mov, type = "l", lwd = 2,
  main = "Stabilization of GO Probability",
  xlab = "Number of Simulated Trials",
  ylab = "Estimated GO Probability"
)

abline(h = mov[length(mov)], col = "red", lwd = 2)

