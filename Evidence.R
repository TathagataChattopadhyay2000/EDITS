## ------------------------------------------------------------
## Shrinkage plot for selected designs (per method, no ESS)
## ------------------------------------------------------------

# 1) Choose one representative design per method from the fine grid:
#    the feasible design with the smallest n for each method.
ref_designs <- oc_results_fine %>%
  dplyr::filter(feasible) %>%
  dplyr::group_by(method) %>%
  dplyr::slice_min(n, with_ties = FALSE) %>%
  dplyr::ungroup()

ref_designs
# columns of interest: method, n, gamma

# 2) For each chosen design, simulate one trial under p1 and extract posterior summaries
shrinkage_list <- list()

for (i in seq_len(nrow(ref_designs))) {
  m         <- as.character(ref_designs$method[i])
  n_ref     <- ref_designs$n[i]
  gamma_ref <- ref_designs$gamma[i]   # for labelling only
  
  cohort_names <- c("p_1", "p_2")
  n_coh        <- length(cohort_names)
  n_per_cohort <- n_ref
  
  # One scenario under p1 for this design
  scen_borrow <- simulateScenarios(
    n_subjects_list     = list(rep(n_per_cohort, n_coh)),
    response_rates_list = list(p1),   # alternative scenario
    n_trials            = 1
  )
  
  analyses_borrow <- performAnalyses(
    scenario_list      = scen_borrow,
    evidence_levels    = seq(0.5, 0.95, by = 0.01),
    target_rates       = p_beta_vec,
    method_names       = m,           # single method here
    n_mcmc_iterations  = 2000,
    verbose            = FALSE
  )
  
  # Extract posterior quantiles for this method
  quantiles_list_all <- analyses_borrow$scenario_1$quantiles_list
  q_list_m           <- quantiles_list_all[[m]]
  if (length(q_list_m) == 0L) next
  
  # qmat is a quantile matrix, rows include "Mean", "2.5%", "97.5%", etc.
  qmat <- q_list_m[[1L]]
  
  # Find the cohort columns (p_1, p_2)
  cohort_cols <- intersect(cohort_names, colnames(qmat))
  if (length(cohort_cols) == 0L) next
  
  mu    <- qmat["Mean",  cohort_cols]
  lower <- qmat["2.5%",  cohort_cols]
  upper <- qmat["97.5%", cohort_cols]
  
  df_m <- data.frame(
    method = m,
    n      = n_ref,
    gamma  = gamma_ref,
    cohort = factor(cohort_cols, levels = cohort_names),
    mean   = as.numeric(mu),
    lower  = as.numeric(lower),
    upper  = as.numeric(upper)
  )
  
  shrinkage_list[[length(shrinkage_list) + 1L]] <- df_m
}

shrinkage_df <- dplyr::bind_rows(shrinkage_list)

shrinkage_df
# method | n | gamma | cohort | mean | lower | upper

# 3) Shrinkage plot: posterior means + 95% CIs per cohort and method
library(ggplot2)

shrinkage_plot <- ggplot(
  shrinkage_df,
  aes(x = cohort,
      y = mean,
      ymin = lower,
      ymax = upper,
      colour = method,
      group = interaction(method, cohort))
) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  geom_errorbar(
    width = 0.15,
    position = position_dodge(width = 0.4)
  ) +
  facet_grid(
    . ~ paste0(
      "method = ", method,
      ", n = ", n,
      ", gamma = ", round(gamma, 2)
    )
  ) +
  labs(
    x = "Cohort",
    y = "Posterior mean response (95% CI)",
    title = "Shrinkage of cohort estimates for selected feasible designs"
  ) +
  theme_minimal()

shrinkage_plot
