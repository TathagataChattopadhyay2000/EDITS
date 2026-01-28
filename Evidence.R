doFuture::registerDoFuture()
future::plan(future::multisession())

library(dplyr)
library(plotly)

# Generic Go-probability function (single final analysis, stratified backend)
go_prob_overall <- function(n,
                            p_true_vec,
                            p_beta_vec,
                            gamma,
                            overall_min_gos = 1,
                            n_trials,
                            method_name = "stratified") {
  
  scen <- simulateScenarios(
    n_subjects_list     = list(c(n, n)),
    response_rates_list = list(p_true_vec),
    n_trials            = n_trials
  )
  
  analyses <- performAnalyses(
    scenario_list      = scen,
    evidence_levels    = seq(0.5, 0.95, by = 0.05),
    target_rates       = p_beta_vec,
    method_names       = method_name,
    n_mcmc_iterations  = 100,
    verbose            = FALSE
  )
  
  decisions <- getGoDecisions(
    analyses_list   = analyses,
    cohort_names    = c("p_1", "p_2"),
    evidence_levels = c(gamma, gamma),
    boundary_rules  = quote(c(x[1] > p_beta_vec[1],
                              x[2] > p_beta_vec[2])),
    overall_min_gos = overall_min_gos
  )
  
  go_probs_list <- getGoProbabilities(decisions)
  go_mat        <- go_probs_list[[method_name]][["scenario_1"]]
  
  go_mat["Go", "overall"]
}

# Design settings --------------------------------------------------------------
p0           <- c(0.30, 0.30)
p1           <- c(0.45, 0.45)
p_beta_vec   <- c(0.40, 0.40)
type_1_error <- 0.20
power_min    <- 0.80
n_trials_oc  <- 100

# Coarse and fine grids (n from 10 to 40) --------------------------------------
n_grid_coarse     <- seq(10, 40, by = 5)
gamma_grid_coarse <- seq(0.5, 0.95, by = 0.05)

# Helper: compute an OC grid for given n/gamma sets and a Go-probability function
compute_oc_grid <- function(n_grid,
                            gamma_grid,
                            p0,
                            p1,
                            p_beta_vec,
                            type_1_error,
                            power_min,
                            n_trials_oc,
                            go_fun) {
  
  oc_results <- expand.grid(
    n     = n_grid,
    gamma = gamma_grid
  )
  
  oc_results$type_1_error_hat <- NA_real_
  oc_results$power_hat        <- NA_real_
  
  for (i in seq_len(nrow(oc_results))) {
    n_i     <- oc_results$n[i]
    gamma_i <- oc_results$gamma[i]
    
    oc_results$type_1_error_hat[i] <- go_fun(
      n            = n_i,
      p_true_vec   = p0,
      p_beta_vec   = p_beta_vec,
      gamma        = gamma_i,
      n_trials     = n_trials_oc
    )
    
    oc_results$power_hat[i] <- go_fun(
      n            = n_i,
      p_true_vec   = p1,
      p_beta_vec   = p_beta_vec,
      gamma        = gamma_i,
      n_trials     = n_trials_oc
    )
  }
  
  oc_results$feasible <- with(
    oc_results,
    type_1_error_hat <= type_1_error & power_hat >= power_min
  )
  
  oc_results
}

set.seed(2026)

# 1) Coarse grid with stratified model -----------------------------------------
oc_results_coarse <- compute_oc_grid(
  n_grid        = n_grid_coarse,
  gamma_grid    = gamma_grid_coarse,
  p0            = p0,
  p1            = p1,
  p_beta_vec    = p_beta_vec,
  type_1_error  = type_1_error,
  power_min     = power_min,
  n_trials_oc   = n_trials_oc,
  go_fun        = go_prob_overall
)

oc_results_coarse$resolution <- "Coarse"

# 2) Define region-of-interest from coarse run and run a fine grid -------------
feas_coarse <- oc_results_coarse[oc_results_coarse$feasible, , drop = FALSE]

if (nrow(feas_coarse) > 0L) {
  n_range     <- range(feas_coarse$n)
  gamma_range <- range(feas_coarse$gamma)
  
  n_grid_fine <- seq(
    from = max(10, n_range[1] - 2),
    to   = min(40, n_range[2] + 2),
    by   = 1
  )
  
  gamma_grid_fine <- seq(
    from = max(0.5,  gamma_range[1] - 0.05),
    to   = min(0.95, gamma_range[2] + 0.05),
    by   = 0.01
  )
  
  oc_results_fine <- compute_oc_grid(
    n_grid        = n_grid_fine,
    gamma_grid    = gamma_grid_fine,
    p0            = p0,
    p1            = p1,
    p_beta_vec    = p_beta_vec,
    type_1_error  = type_1_error,
    power_min     = power_min,
    n_trials_oc   = n_trials_oc,
    go_fun        = go_prob_overall
  )
  
  oc_results_fine$resolution <- "Fine"
  
} else {
  oc_results_fine <- oc_results_coarse[FALSE, , drop = FALSE]
  oc_results_fine$resolution <- "Fine"
}

# 3) Combine results and compute boundaries for each resolution ----------------
oc_all <- bind_rows(oc_results_coarse, oc_results_fine)

oc_all$feasible_factor <- ifelse(oc_all$feasible, "Feasible", "Not feasible")
oc_all$feasible_factor <- factor(
  oc_all$feasible_factor,
  levels = c("Not feasible", "Feasible")
)

calc_boundaries <- function(df) {
  feas_only <- df[df$feasible, , drop = FALSE]
  if (nrow(feas_only) == 0L) {
    return(df[FALSE, , drop = FALSE])
  }
  
  tmp_gamma <- aggregate(n ~ gamma, data = feas_only, FUN = min)
  min_n_per_gamma <- merge(
    tmp_gamma,
    feas_only,
    by = c("gamma", "n"),
    all.x = TRUE,
    sort = TRUE
  )
  
  tmp_n <- aggregate(gamma ~ n, data = feas_only, FUN = min)
  min_gamma_per_n <- merge(
    tmp_n,
    feas_only,
    by = c("n", "gamma"),
    all.x = TRUE,
    sort = TRUE
  )
  
  list(
    min_n_per_gamma   = min_n_per_gamma,
    min_gamma_per_n   = min_gamma_per_n
  )
}

bnd_coarse <- calc_boundaries(oc_results_coarse)
bnd_fine   <- calc_boundaries(oc_results_fine)

md_coarse <- oc_all[oc_all$resolution == "Coarse", , drop = FALSE]
md_fine   <- oc_all[oc_all$resolution == "Fine",   , drop = FALSE]

# 4) Plotly with toggle between coarse / fine / both ---------------------------

fig <- plot_ly() %>%
  add_trace(
    data  = md_coarse,
    x     = ~n,
    y     = ~gamma,
    type  = "scatter",
    mode  = "markers",
    color = ~feasible_factor,
    colors = c("Not feasible" = "red", "Feasible" = "green"),
    text  = ~paste0(
      "Resolution = Coarse",
      "<br>n = ", n,
      "<br>gamma = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    hoverinfo = "text",
    name  = "Coarse designs",
    visible = TRUE
  ) %>%
  add_trace(
    data   = bnd_coarse$min_n_per_gamma,
    x      = ~n,
    y      = ~gamma,
    type   = "scatter",
    mode   = "lines+markers",
    inherit = FALSE,
    line   = list(width = 2),
    marker = list(size = 8, symbol = "x"),
    name   = "Coarse boundary",
    hoverinfo = "text",
    text = ~paste0(
      "Resolution = Coarse",
      "<br>gamma = ", round(gamma, 3),
      "<br>min n = ", n,
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3)
    ),
    visible = TRUE
  ) %>%
  add_trace(
    data  = md_fine,
    x     = ~n,
    y     = ~gamma,
    type  = "scatter",
    mode  = "markers",
    color = ~feasible_factor,
    colors = c("Not feasible" = "red", "Feasible" = "green"),
    text  = ~paste0(
      "Resolution = Fine",
      "<br>n = ", n,
      "<br>gamma = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    hoverinfo = "text",
    name  = "Fine designs",
    visible = FALSE,
    showlegend = FALSE
  ) %>%
  add_trace(
    data   = bnd_fine$min_n_per_gamma,
    x      = ~n,
    y      = ~gamma,
    type   = "scatter",
    mode   = "lines+markers",
    inherit = FALSE,
    line   = list(width = 2),
    marker = list(size = 8, symbol = "x"),
    name   = "Fine boundary",
    hoverinfo = "text",
    text = ~paste0(
      "Resolution = Fine",
      "<br>gamma = ", round(gamma, 3),
      "<br>min n = ", n,
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3)
    ),
    visible = FALSE,
    showlegend = FALSE
  ) %>%
  layout(
    title = "Feasible vs Non-feasible Designs (stratified, 2 cohorts)",
    xaxis = list(title = "Sample size per cohort (n)", range = c(10, 40)),
    yaxis = list(title = "Evidence level (gamma)", range = c(0.5, 0.95)),
    updatemenus = list(
      list(
        type = "dropdown",
        x    = 1.05,
        y    = 1,
        buttons = list(
          list(
            label = "Coarse only",
            method = "restyle",
            args = list("visible", list(TRUE, TRUE, FALSE, FALSE))
          ),
          list(
            label = "Fine only",
            method = "restyle",
            args = list("visible", list(FALSE, FALSE, TRUE, TRUE))
          ),
          list(
            label = "Both",
            method = "restyle",
            args = list("visible", list(TRUE, TRUE, TRUE, TRUE))
          )
        )
      )
    )
  )

fig
