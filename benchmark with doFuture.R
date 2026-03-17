library(bhmbasket)
library(future)
library(doFuture)

# --- Create 100 scenarios ---
set.seed(999)
scenario_list <- list()

for (i in 1:500) {
  trial <- createTrial(
    n_subjects   = matrix(c(30,30,30), nrow = 1),
    n_responders = matrix(rbinom(3, size = 30, prob = c(0.2, 0.25, 0.3)), nrow = 1)
  )
  
  # createTrial returns a scenario_list of length 1
  scenario_list[[i]] <- trial[[1]]
}

class(scenario_list) <- "scenario_list"

# -------------------------------------------------------------
# Benchmark 1: sequential (single worker)
# -------------------------------------------------------------
plan(sequential)
t_seq <- system.time(
  performAnalyses(
    scenario_list       = scenario_list,
    target_rates        = rep(0.3, 3),
    method_names        = "berry",
    n_mcmc_iterations   = 10000,
    verbose             = FALSE
  )
)

# -------------------------------------------------------------
# Benchmark 2: parallel (multisession)
# -------------------------------------------------------------
plan(multisession)  # automatically uses all cores
t_par <- system.time(
  performAnalyses(
    scenario_list       = scenario_list,
    target_rates        = rep(0.3, 3),
    method_names        = "berry",
    n_mcmc_iterations   = 10000,
    verbose             = FALSE
  )
)

rbind(sequential = t_seq, multisession = t_par)


# ------------------------------------------------------------------------------

plan(sequential)     # run sequential
set.seed(111)

res_seq <- performAnalyses(
  scenario_list     = scenario_list,
  target_rates      = rep(0.3, 3),
  method_names      = "berry",
  n_mcmc_iterations = 5000,
  verbose           = FALSE
)

plan(multisession)   # run parallel
set.seed(111)

res_par <- performAnalyses(
  scenario_list     = scenario_list,
  target_rates      = rep(0.3, 3),
  method_names      = "berry",
  n_mcmc_iterations = 5000,
  verbose           = FALSE
)

# Direct comparison
identical(res_seq, res_par)

# ------------------------------------------------------------------------------



