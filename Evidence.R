doFuture::registerDoFuture()
future::plan(future::sequential())

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
    evidence_levels    = seq(0.5, 0.95, by = 0.01),
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

md_coarse_feas <- md_coarse %>% filter(feasible_factor == "Feasible")
md_coarse_nfeas <- md_coarse %>% filter(feasible_factor == "Not feasible")
md_fine_feas <- md_fine %>% filter(feasible_factor == "Feasible")
md_fine_nfeas <- md_fine %>% filter(feasible_factor == "Not feasible")


# 4) Plotly with toggle between coarse / fine / both ---------------------------


bnd_coarse <- calc_boundaries(oc_results_coarse)
bnd_fine   <- calc_boundaries(oc_results_fine)

md_coarse <- oc_all[oc_all$resolution == "Coarse", , drop = FALSE]
md_fine   <- oc_all[oc_all$resolution == "Fine",   , drop = FALSE]

# Split data explicitly to control legend and toggles
md_coarse_feas  <- md_coarse %>% filter(feasible_factor == "Feasible")
md_coarse_nfeas <- md_coarse %>% filter(feasible_factor == "Not feasible")
md_fine_feas    <- md_fine   %>% filter(feasible_factor == "Feasible")
md_fine_nfeas   <- md_fine   %>% filter(feasible_factor == "Not feasible")

fig <- plot_ly() %>%
  # 1) Coarse — Feasible (trace 1)
  add_trace(
    data  = md_coarse_feas,
    x     = ~n, y = ~gamma,
    type  = "scatter", mode = "markers",
    marker = list(color = "green"),
    name  = "Coarse — Feasible",
    legendgroup = "coarse",
    hoverinfo = "text",
    text  = ~paste0(
      "Resolution = Coarse",
      "<br>n = ", n,
      "<br>gamma = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    visible = TRUE
  ) %>%
  # 2) Coarse — Not feasible (trace 2)
  add_trace(
    data  = md_coarse_nfeas,
    x     = ~n, y = ~gamma,
    type  = "scatter", mode = "markers",
    marker = list(color = "red"),
    name  = "Coarse — Not feasible",
    legendgroup = "coarse",
    hoverinfo = "text",
    text  = ~paste0(
      "Resolution = Coarse",
      "<br>n = ", n,
      "<br>gamma = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    visible = TRUE
  ) %>%
  # 3) Coarse boundary (trace 3)
  add_trace(
    data   = bnd_coarse$min_n_per_gamma,
    x      = ~n, y = ~gamma,
    type   = "scatter", mode = "lines+markers",
    inherit = FALSE,
    line   = list(width = 2),
    marker = list(size = 8, symbol = "x", color = "black"),
    name   = "Coarse boundary",
    legendgroup = "coarse",
    hoverinfo = "text",
    text = ~paste0(
      "Resolution = Coarse",
      "<br>gamma = ", round(gamma, 3),
      "<br>min n = ", n,
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3)
    ),
    showlegend = FALSE,
    visible = TRUE
  ) %>%
  # 4) Fine — Feasible (trace 4)
  add_trace(
    data  = md_fine_feas,
    x     = ~n, y = ~gamma,
    type  = "scatter", mode = "markers",
    marker = list(color = "green"),
    name  = "Fine — Feasible",
    legendgroup = "fine",
    hoverinfo = "text",
    text  = ~paste0(
      "Resolution = Fine",
      "<br>n = ", n,
      "<br>gamma = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    visible = FALSE
  ) %>%
  # 5) Fine — Not feasible (trace 5)
  add_trace(
    data  = md_fine_nfeas,
    x     = ~n, y = ~gamma,
    type  = "scatter", mode = "markers",
    marker = list(color = "red"),
    name  = "Fine — Not feasible",
    legendgroup = "fine",
    hoverinfo = "text",
    text  = ~paste0(
      "Resolution = Fine",
      "<br>n = ", n,
      "<br>gamma = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    visible = FALSE
  ) %>%
  # 6) Fine boundary (trace 6)
  add_trace(
    data   = bnd_fine$min_n_per_gamma,
    x      = ~n, y = ~gamma,
    type   = "scatter", mode = "lines+markers",
    inherit = FALSE,
    line   = list(width = 2),
    marker = list(size = 8, symbol = "x", color = "black"),
    name   = "Fine boundary",
    legendgroup = "fine",
    hoverinfo = "text",
    text = ~paste0(
      "Resolution = Fine",
      "<br>gamma = ", round(gamma, 3),
      "<br>min n = ", n,
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3)
    ),
    showlegend = FALSE,
    visible = FALSE
  ) %>%
  layout(
    title = "Feasible vs Non-feasible Designs (stratified, 2 cohorts)",
    xaxis = list(title = "Sample size per cohort (n)", range = c(10, 40), fixedrange = TRUE),
    yaxis = list(title = "Evidence level (gamma)", range = c(0.5, 0.95), fixedrange = TRUE),
    legend = list(title = list(text = "Designs")),
    updatemenus = list(
      # (1) Resolution selector: controls VISIBILITY
      list(
        type = "dropdown", x = 1.12, y = 1,
        buttons = list(
          list(
            label = "Coarse only",
            method = "update",
            args = list(
              # [1 C Feas, 2 C Not, 3 C Bound, 4 F Feas, 5 F Not, 6 F Bound]
              list(visible = c(TRUE, TRUE, TRUE,  FALSE, FALSE, FALSE)),
              list(title = "Feasible vs Non-feasible — Coarse")
            )
          ),
          list(
            label = "Fine only",
            method = "update",
            args = list(
              list(visible = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)),
              list(title = "Feasible vs Non-feasible — Fine")
            )
          ),
          list(
            label = "Both",
            method = "update",
            args = list(
              list(visible = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)),
              list(title = "Feasible vs Non-feasible — Coarse & Fine")
            )
          )
        )
      ),
      # (2) Feasibility selector: controls MARKER OPACITY ONLY (keeps resolution choice intact)
      list(
        type = "dropdown", x = 1.12, y = 0.92,
        buttons = list(
          list(
            label = "Feasible only",
            method = "restyle",
            args = list("marker.opacity", list(1, 0, 1, 1, 0, 1))
          ),
          list(
            label = "Not feasible only",
            method = "restyle",
            args = list("marker.opacity", list(0, 1, 1, 0, 1, 1))
          ),
          list(
            label = "Both feasibilities",
            method = "restyle",
            args = list("marker.opacity", list(1, 1, 1, 1, 1, 1))
          )
        )
      )
    )
  )

fig

################################################################################

# ---- Robust extractor for per-trial 'overall' decisions ----
extract_overall_vec <- function(decisions, method_name = NULL) {
  # If it's already an atomic vector (logical/numeric), return it as logical
  if (!is.list(decisions)) {
    return(as.logical(decisions))
  }
  
  # Helper to try a node and get overall from it
  try_overall <- function(node) {
    if (is.null(node)) return(NULL)
    if (!is.list(node)) return(as.logical(node))
    if (!is.null(node$overall_vec)) return(as.logical(node$overall_vec))
    if (!is.null(node$overall))     return(as.logical(node$overall))
    NULL
  }
  
  # 1) decisions$scenario_1$decisions_list[[method_name]]
  if (!is.null(decisions$scenario_1) && !is.null(decisions$scenario_1$decisions_list)) {
    if (!is.null(method_name) && !is.null(decisions$scenario_1$decisions_list[[method_name]])) {
      out <- try_overall(decisions$scenario_1$decisions_list[[method_name]])
      if (!is.null(out)) return(out)
    }
    # If method_name is not inside decisions_list, try any element
    for (nm in names(decisions$scenario_1$decisions_list)) {
      out <- try_overall(decisions$scenario_1$decisions_list[[nm]])
      if (!is.null(out)) return(out)
    }
  }
  
  # 2) decisions$decisions_list[[method_name]]$scenario_1
  if (!is.null(decisions$decisions_list)) {
    if (!is.null(method_name) && !is.null(decisions$decisions_list[[method_name]])) {
      out <- try_overall(decisions$decisions_list[[method_name]]$scenario_1)
      if (!is.null(out)) return(out)
    }
    for (nm in names(decisions$decisions_list)) {
      out <- try_overall(decisions$decisions_list[[nm]]$scenario_1)
      if (!is.null(out)) return(out)
    }
  }
  
  # 3) decisions[[method_name]]$scenario_1
  if (!is.null(method_name) && !is.null(decisions[[method_name]])) {
    out <- try_overall(decisions[[method_name]]$scenario_1)
    if (!is.null(out)) return(out)
  }
  
  # 4) Top-level fields (sometimes returned directly)
  out <- try_overall(decisions$overall_vec); if (!is.null(out)) return(out)
  out <- try_overall(decisions$overall);     if (!is.null(out)) return(out)
  
  stop("Could not locate an 'overall' logical vector in decisions. Inspect structure with str(decisions, max.level = 3).")
}

# ---- Safer two-stage function using extractor ----
go_prob_overall_interim <- function(n,
                                    p_true_vec,
                                    p_beta_vec,
                                    gamma_interim,
                                    gamma_final,
                                    overall_min_gos = 1L,
                                    n_trials,
                                    method_name = "stratified") {
  stopifnot(length(p_true_vec) == 2L, length(p_beta_vec) == 2L)
  
  n_int  <- floor(n / 2)
  n_add  <- n - n_int
  cohort_names <- c("p_1", "p_2")
  
  # ---------- Stage 1 ----------
  scen_stage1 <- simulateScenarios(
    n_subjects_list     = list(rep(n_int, 2L)),
    response_rates_list = list(p_true_vec),
    n_trials            = n_trials
  )
  
  analysis_stage1 <- performAnalyses(
    scenario_list      = scen_stage1,
    evidence_levels    = c(gamma_interim, gamma_interim),
    target_rates       = p_beta_vec,
    method_names       = method_name,
    n_mcmc_iterations  = 100,
    verbose            = FALSE
  )
  
  decisions_stage1 <- getGoDecisions(
    analyses_list   = analysis_stage1,
    cohort_names    = cohort_names,
    evidence_levels = c(gamma_interim, gamma_interim),
    boundary_rules  = quote(c(x[1] > p_beta_vec[1],
                              x[2] > p_beta_vec[2])),
    overall_min_gos = overall_min_gos
  )
  
  overall_stage1 <- extract_overall_vec(decisions_stage1, method_name)
  
  # ---------- Stage 2 ----------
  scen_stage2 <- continueRecruitment(
    n_subjects_add_list = list(rep(n_add, 2L)),
    decisions_list      = decisions_stage1,
    method_name         = method_name
  )
  
  analysis_stage2 <- performAnalyses(
    scenario_list      = scen_stage2,
    evidence_levels    = c(gamma_final, gamma_final),
    target_rates       = p_beta_vec,
    method_names       = method_name,
    n_mcmc_iterations  = 100,
    verbose            = FALSE
  )
  
  decisions_stage2 <- getGoDecisions(
    analyses_list   = analysis_stage2,
    cohort_names    = cohort_names,
    evidence_levels = c(gamma_final, gamma_final),
    boundary_rules  = quote(c(x[1] > p_beta_vec[1],
                              x[2] > p_beta_vec[2])),
    overall_min_gos = overall_min_gos
  )
  
  overall_final <- extract_overall_vec(decisions_stage2, method_name)
  
  # Early futility: require interim overall TRUE + final overall TRUE
  # (elementwise AND over n_trials replicates)
  overall_final_adj <- as.logical(overall_stage1) & as.logical(overall_final)
  
  mean(overall_final_adj)
}

gamma_interim_fix <- 0.60

go_fun_interim <- function(n, p_true_vec, p_beta_vec, gamma, n_trials) {
  go_prob_overall_interim(
    n              = n,
    p_true_vec     = p_true_vec,
    p_beta_vec     = p_beta_vec,
    gamma_interim  = gamma_interim_fix,
    gamma_final    = gamma,
    overall_min_gos = 1L,
    n_trials       = n_trials,
    method_name    = "stratified"
  )
}

set.seed(3030)

oc_results_coarse_int <- compute_oc_grid(
  n_grid        = n_grid_coarse,
  gamma_grid    = gamma_grid_coarse,
  p0            = p0,
  p1            = p1,
  p_beta_vec    = p_beta_vec,
  type_1_error  = type_1_error,
  power_min     = power_min,
  n_trials_oc   = n_trials_oc,
  go_fun        = go_fun_interim
)

oc_results_coarse_int$resolution  <- "Coarse"
oc_results_coarse_int$design_type <- "Interim"

feas_coarse_int <- oc_results_coarse_int[oc_results_coarse_int$feasible, , drop = FALSE]

if (nrow(feas_coarse_int) > 0L) {
  n_range_int     <- range(feas_coarse_int$n)
  gamma_range_int <- range(feas_coarse_int$gamma)
  
  n_grid_fine_int <- seq(
    from = max(10, n_range_int[1] - 2),
    to   = min(40, n_range_int[2] + 2),
    by   = 1
  )
  
  gamma_grid_fine_int <- seq(
    from = max(0.5,  gamma_range_int[1] - 0.05),
    to   = min(0.95, gamma_range_int[2] + 0.05),
    by   = 0.01
  )
  
  oc_results_fine_int <- compute_oc_grid(
    n_grid        = n_grid_fine_int,
    gamma_grid    = gamma_grid_fine_int,
    p0            = p0,
    p1            = p1,
    p_beta_vec    = p_beta_vec,
    type_1_error  = type_1_error,
    power_min     = power_min,
    n_trials_oc   = n_trials_oc,
    go_fun        = go_fun_interim
  )
  
  oc_results_fine_int$resolution  <- "Fine"
  oc_results_fine_int$design_type <- "Interim"
  
} else {
  oc_results_fine_int <- oc_results_coarse_int[FALSE, , drop = FALSE]
  oc_results_fine_int$resolution  <- character(0)    # length 0, OK for 0 rows
  oc_results_fine_int$design_type <- character(0)    # length 0, OK for 0 rows
}

oc_all_int <- dplyr::bind_rows(oc_results_coarse_int, oc_results_fine_int)

oc_all_int$feasible_factor <- ifelse(oc_all_int$feasible, "Feasible", "Not feasible")
oc_all_int$feasible_factor <- factor(
  oc_all_int$feasible_factor,
  levels = c("Not feasible", "Feasible")
)

bnd_coarse_int <- calc_boundaries(oc_results_coarse_int)
bnd_fine_int   <- calc_boundaries(oc_results_fine_int)

md_coarse_int <- oc_all_int[oc_all_int$resolution == "Coarse", , drop = FALSE]
md_fine_int   <- oc_all_int[oc_all_int$resolution == "Fine",   , drop = FALSE]

md_coarse_feas_int  <- md_coarse_int %>% dplyr::filter(feasible_factor == "Feasible")
md_coarse_nfeas_int <- md_coarse_int %>% dplyr::filter(feasible_factor == "Not feasible")
md_fine_feas_int    <- md_fine_int   %>% dplyr::filter(feasible_factor == "Feasible")
md_fine_nfeas_int   <- md_fine_int   %>% dplyr::filter(feasible_factor == "Not feasible")

fig_interim <- plot_ly() %>%
  add_trace(
    data  = md_coarse_feas_int,
    x     = ~n, y = ~gamma,
    type  = "scatter", mode = "markers",
    marker = list(color = "green"),
    name  = "Coarse — Feasible (interim)",
    legendgroup = "coarse_int",
    hoverinfo = "text",
    text  = ~paste0(
      "Design = Interim",
      "<br>Resolution = Coarse",
      "<br>n = ", n,
      "<br>gamma_final = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    visible = TRUE
  ) %>%
  add_trace(
    data  = md_coarse_nfeas_int,
    x     = ~n, y = ~gamma,
    type  = "scatter", mode = "markers",
    marker = list(color = "red"),
    name  = "Coarse — Not feasible (interim)",
    legendgroup = "coarse_int",
    hoverinfo = "text",
    text  = ~paste0(
      "Design = Interim",
      "<br>Resolution = Coarse",
      "<br>n = ", n,
      "<br>gamma_final = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    visible = TRUE
  ) %>%
  add_trace(
    data   = bnd_coarse_int$min_n_per_gamma,
    x      = ~n, y = ~gamma,
    type   = "scatter", mode = "lines+markers",
    inherit = FALSE,
    line   = list(width = 2),
    marker = list(size = 8, symbol = "x", color = "black"),
    name   = "Coarse boundary (interim)",
    legendgroup = "coarse_int",
    hoverinfo = "text",
    text = ~paste0(
      "Design = Interim",
      "<br>Resolution = Coarse",
      "<br>gamma_final = ", round(gamma, 3),
      "<br>min n = ", n,
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3)
    ),
    showlegend = FALSE,
    visible = TRUE
  ) %>%
  add_trace(
    data  = md_fine_feas_int,
    x     = ~n, y = ~gamma,
    type  = "scatter", mode = "markers",
    marker = list(color = "green"),
    name  = "Fine — Feasible (interim)",
    legendgroup = "fine_int",
    hoverinfo = "text",
    text  = ~paste0(
      "Design = Interim",
      "<br>Resolution = Fine",
      "<br>n = ", n,
      "<br>gamma_final = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    visible = FALSE
  ) %>%
  add_trace(
    data  = md_fine_nfeas_int,
    x     = ~n, y = ~gamma,
    type  = "scatter", mode = "markers",
    marker = list(color = "red"),
    name  = "Fine — Not feasible (interim)",
    legendgroup = "fine_int",
    hoverinfo = "text",
    text  = ~paste0(
      "Design = Interim",
      "<br>Resolution = Fine",
      "<br>n = ", n,
      "<br>gamma_final = ", round(gamma, 3),
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3),
      "<br>Feasible = ", feasible
    ),
    visible = FALSE
  ) %>%
  add_trace(
    data   = bnd_fine_int$min_n_per_gamma,
    x      = ~n, y = ~gamma,
    type   = "scatter", mode = "lines+markers",
    inherit = FALSE,
    line   = list(width = 2),
    marker = list(size = 8, symbol = "x", color = "black"),
    name   = "Fine boundary (interim)",
    legendgroup = "fine_int",
    hoverinfo = "text",
    text = ~paste0(
      "Design = Interim",
      "<br>Resolution = Fine",
      "<br>gamma_final = ", round(gamma, 3),
      "<br>min n = ", n,
      "<br>Type I error = ", round(type_1_error_hat, 3),
      "<br>Power = ", round(power_hat, 3)
    ),
    showlegend = FALSE,
    visible = FALSE
  ) %>%
  layout(
    title = paste0(
      "Feasible vs Non-feasible Designs (stratified, 2 cohorts, interim; ",
      "gamma_interim = ", gamma_interim_fix, ")"
    ),
    xaxis = list(title = "Sample size per cohort (n)", range = c(10, 40), fixedrange = TRUE),
    yaxis = list(title = "Final evidence level (gamma)", range = c(0.5, 0.95), fixedrange = TRUE),
    legend = list(title = list(text = "Designs")),
    updatemenus = list(
      list(
        type = "dropdown", x = 1.12, y = 1,
        buttons = list(
          list(
            label = "Coarse only",
            method = "update",
            args = list(
              list(visible = c(TRUE, TRUE, TRUE,  FALSE, FALSE, FALSE)),
              list(title = paste0(
                "Feasible vs Non-feasible — Coarse (interim; gamma_interim = ",
                gamma_interim_fix, ")"
              ))
            )
          ),
          list(
            label = "Fine only",
            method = "update",
            args = list(
              list(visible = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)),
              list(title = paste0(
                "Feasible vs Non-feasible — Fine (interim; gamma_interim = ",
                gamma_interim_fix, ")"
              ))
            )
          ),
          list(
            label = "Both",
            method = "update",
            args = list(
              list(visible = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)),
              list(title = paste0(
                "Feasible vs Non-feasible — Coarse & Fine (interim; gamma_interim = ",
                gamma_interim_fix, ")"
              ))
            )
          )
        )
      ),
      list(
        type = "dropdown", x = 1.12, y = 0.92,
        buttons = list(
          list(
            label = "Feasible only",
            method = "restyle",
            args = list("marker.opacity", list(1, 0, 1, 1, 0, 1))
          ),
          list(
            label = "Not feasible only",
            method = "restyle",
            args = list("marker.opacity", list(0, 1, 1, 0, 1, 1))
          ),
          list(
            label = "Both feasibilities",
            method = "restyle",
            args = list("marker.opacity", list(1, 1, 1, 1, 1, 1))
          )
        )
      )
    )
  )

fig_interim
