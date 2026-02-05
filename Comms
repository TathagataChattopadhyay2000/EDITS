# 02_oc_calculations.R
# Functions to compute OC grids (coarse + fine) for given inputs
# ---------------------------------------------------------------------------

run_oc_calculations <- function(p0,
                                p1,
                                p_beta_vec,
                                type_1_error,
                                power_min,
                                n_trials_oc,
                                method_names,
                                n_grid_coarse,
                                gamma_grid_coarse,
                                n_fine_margin_n     = 2,
                                n_fine_margin_gamma = 0.05) {
  n_coh <- length(p0)

  if (!identical(length(p1),       n_coh) ||
      !identical(length(p_beta_vec), n_coh)) {
    stop("p0, p1, and p_beta_vec must all have the same length (number of cohorts).")
  }

  cohort_names <- paste0("p_", seq_len(n_coh))

  # Helper: run OC for one (n, gamma) pair across all methods
  run_oc_for_pair <- function(n_i, gamma_i) {
    # Type I error (H0)
    scen0 <- simulateScenarios(
      n_subjects_list     = list(rep(n_i, n_coh)),
      response_rates_list = list(p0),
      n_trials            = n_trials_oc
    )

    analyses0 <- performAnalyses(
      scenario_list      = scen0,
      evidence_levels    = seq(0.5, 0.95, by = 0.01),
      target_rates       = p_beta_vec,
      method_names       = method_names,
      n_mcmc_iterations  = 100,
      verbose            = FALSE
    )

    decisions0 <- getGoDecisions(
      analyses_list   = analyses0,
      cohort_names    = cohort_names,
      evidence_levels = rep(gamma_i, n_coh),
      boundary_rules  = quote(x > p_beta_vec),
      overall_min_gos = 1
    )

    go_probs0 <- getGoProbabilities(decisions0)

    # Power (H1)
    scen1 <- simulateScenarios(
      n_subjects_list     = list(rep(n_i, n_coh)),
      response_rates_list = list(p1),
      n_trials            = n_trials_oc
    )

    analyses1 <- performAnalyses(
      scenario_list      = scen1,
      evidence_levels    = seq(0.5, 0.95, by = 0.01),
      target_rates       = p_beta_vec,
      method_names       = method_names,
      n_mcmc_iterations  = 100,
      verbose            = FALSE
    )

    decisions1 <- getGoDecisions(
      analyses_list   = analyses1,
      cohort_names    = cohort_names,
      evidence_levels = rep(gamma_i, n_coh),
      boundary_rules  = quote(x > p_beta_vec),
      overall_min_gos = 1
    )

    go_probs1 <- getGoProbabilities(decisions1)

    out <- data.frame(
      method           = method_names,
      type_1_error_hat = NA_real_,
      power_hat        = NA_real_,
      stringsAsFactors = FALSE
    )

    for (m in method_names) {
      go_mat0 <- go_probs0[[m]][["scenario_1"]]
      go_mat1 <- go_probs1[[m]][["scenario_1"]]
      out$type_1_error_hat[out$method == m] <- go_mat0["Go", "overall"]
      out$power_hat[out$method == m]        <- go_mat1["Go", "overall"]
    }

    out
  }

  # Coarse grid
  oc_results_coarse <- expand.grid(
    n      = n_grid_coarse,
    gamma  = gamma_grid_coarse,
    method = method_names,
    stringsAsFactors = FALSE
  )

  oc_results_coarse$type_1_error_hat <- NA_real_
  oc_results_coarse$power_hat        <- NA_real_

  pairs_coarse <- unique(oc_results_coarse[c("n", "gamma")])

  for (i in seq_len(nrow(pairs_coarse))) {
    n_i     <- pairs_coarse$n[i]
    gamma_i <- pairs_coarse$gamma[i]

    oc_pair <- run_oc_for_pair(n_i, gamma_i)

    for (m in method_names) {
      rid <- which(
        oc_results_coarse$n == n_i &
          oc_results_coarse$gamma == gamma_i &
          oc_results_coarse$method == m
      )
      oc_results_coarse$type_1_error_hat[rid] <- oc_pair$type_1_error_hat[oc_pair$method == m]
      oc_results_coarse$power_hat[rid]        <- oc_pair$power_hat[oc_pair$method == m]
    }
  }

  oc_results_coarse$feasible <- with(
    oc_results_coarse,
    type_1_error_hat <= type_1_error & power_hat >= power_min
  )
  oc_results_coarse$resolution <- "Coarse"

  # Fine grid: zoom into feasible coarse region
  feas_coarse <- oc_results_coarse[oc_results_coarse$feasible, , drop = FALSE]

  if (nrow(feas_coarse) > 0L) {
    n_range     <- range(feas_coarse$n)
    gamma_range <- range(feas_coarse$gamma)

    n_grid_fine <- seq(
      from = max(min(n_grid_coarse),  n_range[1] - n_fine_margin_n),
      to   = min(max(n_grid_coarse),  n_range[2] + n_fine_margin_n),
      by   = 1
    )

    gamma_grid_fine <- seq(
      from = max(min(gamma_grid_coarse),  gamma_range[1] - n_fine_margin_gamma),
      to   = min(max(gamma_grid_coarse),  gamma_range[2] + n_fine_margin_gamma),
      by   = 0.01
    )

    oc_results_fine <- expand.grid(
      n      = n_grid_fine,
      gamma  = gamma_grid_fine,
      method = method_names,
      stringsAsFactors = FALSE
    )

    oc_results_fine$type_1_error_hat <- NA_real_
    oc_results_fine$power_hat        <- NA_real_

    pairs_fine <- unique(oc_results_fine[c("n", "gamma")])

    for (i in seq_len(nrow(pairs_fine))) {
      n_i     <- pairs_fine$n[i]
      gamma_i <- pairs_fine$gamma[i]

      oc_pair <- run_oc_for_pair(n_i, gamma_i)

      for (m in method_names) {
        rid <- which(
          oc_results_fine$n == n_i &
            oc_results_fine$gamma == gamma_i &
            oc_results_fine$method == m
        )
        oc_results_fine$type_1_error_hat[rid] <- oc_pair$type_1_error_hat[oc_pair$method == m]
        oc_results_fine$power_hat[rid]        <- oc_pair$power_hat[oc_pair$method == m]
      }
    }

    oc_results_fine$feasible <- with(
      oc_results_fine,
      type_1_error_hat <= type_1_error & power_hat >= power_min
    )
    oc_results_fine$resolution <- "Fine"

  } else {
    oc_results_fine <- oc_results_coarse[FALSE, , drop = FALSE]
    oc_results_fine$resolution <- "Fine"
  }

  # Combine coarse + fine
  oc_all <- dplyr::bind_rows(oc_results_coarse, oc_results_fine)

  oc_all$feasible_factor <- ifelse(oc_all$feasible, "Feasible", "Not feasible")
  oc_all$feasible_factor <- factor(
    oc_all$feasible_factor,
    levels = c("Not feasible", "Feasible")
  )

  list(
    oc_results_coarse = oc_results_coarse,
    oc_results_fine   = oc_results_fine,
    oc_all            = oc_all
  )
}
