# Borrowing / shrinkage integrated with OC grid ----------------

# We assume you already have:
#   oc_results_fine  with columns n, gamma, method, feasible, ...
#   p1, p_beta_vec, method_names, etc. defined above.

# 1) Choose one representative design per method from the fine grid
#    Here: the feasible design with the smallest n for each method.
ref_designs <- oc_results_fine %>%
  dplyr::filter(feasible) %>%
  dplyr::group_by(method) %>%
  dplyr::slice_min(n, with_ties = FALSE) %>%
  dplyr::ungroup()

ref_designs
# Columns of interest: method, n, gamma (and type_1_error_hat, power_hat if you want)

# 2) For each chosen design, simulate one trial and get posterior summaries
shrinkage_list <- list()

for (i in seq_len(nrow(ref_designs))) {
  m       <- as.character(ref_designs$method[i])
  n_ref   <- ref_designs$n[i]
  gamma_ref <- ref_designs$gamma[i]   # for labelling only; does not affect posterior
  
  cohort_names   <- c("p_1", "p_2")
  n_coh_borrow   <- length(cohort_names)
  n_borrow       <- n_ref
  
  # One scenario under p1 for this design
  scen_borrow <- simulateScenarios(
    n_subjects_list     = list(rep(n_borrow, n_coh_borrow)),
    response_rates_list = list(p1),        # alternative scenario (can change if you want)
    n_trials            = 1
  )
  
  analyses_borrow <- performAnalyses(
    scenario_list      = scen_borrow,
    evidence_levels    = seq(0.5, 0.95, by = 0.01),
    target_rates       = p_beta_vec,
    method_names       = m,               # single method here
    n_mcmc_iterations  = 2000,
    verbose            = FALSE
  )
  
  # Extract posterior quantiles for this method
  quantiles_list_all <- analyses_borrow$scenario_1$quantiles_list
  
  # If the method name is used as list name:
  q_list_m <- quantiles_list_all[[m]]
  if (length(q_list_m) == 0L) next
  
  qmat <- q_list_m[[1L]]   # matrix [51 x #params], rows include "Mean", "SD", "2.5%", "97.5%"
  
  # Find p_1, p_2 columns (cohorts)
  idx <- match(cohort_names, colnames(qmat))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0L) next
  
  mu    <- qmat["Mean",   idx]
  sd    <- qmat["SD",     idx]
  lower <- qmat["2.5%",   idx]
  upper <- qmat["97.5%",  idx]
  
  # Moment-match to Beta(alpha,beta) to get ESS = alpha + beta
  var <- sd^2
  S   <- (mu * (1 - mu) / var) - 1
  ESS <- pmax(S, 0)                   # guard against numerical issues
  borrowed <- pmax(ESS - n_borrow, 0) # "extra" pseudo-patients
  
  df_m <- data.frame(
    method   = m,
    n        = n_ref,
    gamma    = gamma_ref,
    cohort   = seq_along(idx),
    mean     = as.numeric(mu),
    lower    = as.numeric(lower),
    upper    = as.numeric(upper),
    sd       = as.numeric(sd),
    ESS      = as.numeric(ESS),
    borrowed = as.numeric(borrowed),
    n_obs    = n_borrow
  )
  
  shrinkage_list[[length(shrinkage_list) + 1L]] <- df_m
}

shrinkage_df <- dplyr::bind_rows(shrinkage_list)

shrinkage_df
# You should now see rows like:
# method | n | gamma | cohort | mean | lower | upper | ESS | borrowed | n_obs

# 3) Plot shrinkage of posterior means (per method / cohort / design)
shrinkage_plot <- ggplot(
  shrinkage_df,
  aes(x = factor(cohort),
      y = mean,
      ymin = lower,
      ymax = upper,
      colour = method)
) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  geom_errorbar(width = 0.2,
                position = position_dodge(width = 0.4)) +
  facet_grid(. ~ paste0("n = ", n, ", gamma = ", round(gamma, 2))) +
  labs(
    x = "Cohort",
    y = "Posterior mean response (95% CI)",
    title = "Shrinkage of cohort estimates for selected designs"
  ) +
  theme_minimal()

# 4) Plot implied borrowed ESS (ESS − n_obs) per cohort and method
borrowed_plot <- ggplot(
  shrinkage_df,
  aes(x = factor(cohort),
      y = borrowed,
      fill = method)
) +
  geom_col(position = "dodge") +
  facet_grid(. ~ paste0("n = ", n, ", gamma = ", round(gamma, 2))) +
  labs(
    x = "Cohort",
    y = "Borrowed effective sample size (ESS − n_obs)",
    title = "Implied borrowing for selected designs"
  ) +
  theme_minimal()

shrinkage_plot
borrowed_plot
