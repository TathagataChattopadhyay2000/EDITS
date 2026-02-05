# 03_oc_plots.R
# Plot functions: main OC plot and boundary-only plot
# ---------------------------------------------------------------------------

plot_oc_main <- function(oc_all, method_names) {

  # Decide resolution order
  oc_all$resolution <- factor(
    oc_all$resolution,
    levels = c("Coarse", "Fine")
  )

  # Compute boundaries (min n per gamma within each method + resolution)
  calc_boundaries <- function(df) {
    feas_only <- df[df$feasible, , drop = FALSE]
    if (nrow(feas_only) == 0L) {
      empty <- df[FALSE, , drop = FALSE]
      return(empty)
    }

    tmp_gamma <- aggregate(n ~ gamma + method + resolution, data = feas_only, FUN = min)
    min_n_per_gamma <- merge(
      tmp_gamma,
      feas_only,
      by = c("gamma", "method", "resolution", "n"),
      all.x = TRUE,
      sort = TRUE
    )
    min_n_per_gamma
  }

  bnd_all <- calc_boundaries(oc_all)

  oc_plot_gg <- ggplot(oc_all, aes(x = n, y = gamma)) +
    geom_point(aes(
      color = feasible_factor,
      shape = resolution
    ),
    size = 2,
    alpha = 0.8) +
    geom_line(
      data = bnd_all,
      aes(
        x = n,
        y = gamma,
        linetype = resolution,
        group = interaction(method, resolution)
      ),
      inherit.aes = FALSE
    ) +
    geom_point(
      data = bnd_all,
      aes(
        x = n,
        y = gamma,
        shape = resolution,
        group = interaction(method, resolution)
      ),
      inherit.aes = FALSE,
      size = 2
    ) +
    facet_wrap(~ method) +
    scale_color_manual(
      values = c("Not feasible" = "red", "Feasible" = "darkgreen"),
      name   = "Design feasibility"
    ) +
    scale_shape_manual(
      values = c("Coarse" = 16, "Fine" = 17),
      name   = "Grid resolution"
    ) +
    labs(
      title = "Operating characteristics for bhmbasket basket designs",
      subtitle = paste(
        "Methods:", paste(method_names, collapse = ", "),
        "| cohorts:", length(unique(gsub("p_", "", grep("^p_", names(oc_all), value = TRUE)))),
        "| points:", nrow(oc_all)
      ),
      x = "Sample size per cohort (n)",
      y = "Evidence level (gamma)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom"
    )

  ggplotly(oc_plot_gg)
}

plot_boundary_only <- function(oc_results_coarse,
                               oc_results_fine,
                               method_names) {

  calc_boundaries_basic <- function(df, res_label) {
    df$resolution <- res_label
    feas_only <- df[df$feasible, , drop = FALSE]
    if (nrow(feas_only) == 0L) {
      return(df[FALSE, , drop = FALSE])
    }

    tmp_gamma <- aggregate(n ~ gamma + method, data = feas_only, FUN = min)
    min_n_per_gamma <- merge(
      tmp_gamma,
      feas_only,
      by = c("gamma", "method", "n"),
      all.x = TRUE,
      sort = TRUE
    )
    min_n_per_gamma$resolution <- res_label
    min_n_per_gamma
  }

  bnd_coarse <- calc_boundaries_basic(oc_results_coarse, "Coarse")
  bnd_fine   <- calc_boundaries_basic(oc_results_fine,   "Fine")

  bnd_all <- dplyr::bind_rows(bnd_coarse, bnd_fine)

  bnd_all$resolution <- factor(
    bnd_all$resolution,
    levels = c("Coarse", "Fine")
  )

  bnd_plot_gg <- ggplot(bnd_all, aes(x = gamma, y = n)) +
    geom_line(aes(
      linetype = resolution,
      color    = resolution,
      group    = interaction(method, resolution)
    )) +
    geom_point(aes(
      shape = resolution,
      color = resolution
    ), size = 2) +
    facet_wrap(~ method) +
    scale_color_manual(
      values = c("Coarse" = "black", "Fine" = "blue"),
      name   = "Grid resolution"
    ) +
    scale_shape_manual(
      values = c("Coarse" = 16, "Fine" = 17),
      name   = "Grid resolution"
    ) +
    labs(
      title = "Decision boundaries: minimum n per evidence level",
      subtitle = paste(
        "Methods:", paste(method_names, collapse = ", ")
      ),
      x = "Evidence level (gamma)",
      y = "Minimum sample size per cohort (n)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom"
    )

  ggplotly(bnd_plot_gg)
}
