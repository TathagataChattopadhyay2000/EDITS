# 01_user_input.R
# User-facing script: set design inputs and run everything
# ---------------------------------------------------------------------------

## Packages
library(bhmbasket)
library(dplyr)
library(ggplot2)
library(plotly)
library(future)
library(doFuture)
library(foreach)

doFuture::registerDoFuture()
future::plan(future::sequential())  # user can also choose multisession()

## Source helper files (functions only, do not run anything yet)
source("02_oc_calculations.R")
source("03_oc_plots.R")

## -------------------- USER INPUT BLOCK ---------------------------------------
## Only edit this block to change the design.

# True response rates under H0 and H1
# Number of cohorts = length of these vectors.
# To add a cohort, just append another value to p0, p1 and p_beta_vec.
p0         <- c(0.40, 0.40)   # e.g. 2 cohorts
p1         <- c(0.55, 0.55)
p_beta_vec <- c(0.50, 0.50)   # clinical boundary per cohort

# Design targets
type_1_error <- 0.20          # maximum allowed type I error (overall Go)
power_min    <- 0.80          # required power (overall Go)

# Number of simulated trials per (n, gamma) combo
n_trials_oc  <- 100

# Methods to compare (must be available in bhmbasket)
# E.g. "stratified", "berry", "pooled", "exnex", ...
method_names <- c("stratified", "berry")

# Coarse grid for n and gamma (broad scan)
n_grid_coarse     <- seq(10, 40, by = 5)        # candidate sample sizes per cohort
gamma_grid_coarse <- seq(0.5, 0.95, by = 0.05)  # candidate evidence levels

# Fine-grid margins (for automatic zoom around feasible region)
n_fine_margin_n     <- 2     # extend feasible n-range by ± this many patients
n_fine_margin_gamma <- 0.05  # extend feasible gamma-range by ± this

# Random seed for reproducibility
set.seed(2026)

## -------------------- RUN CALCULATIONS ---------------------------------------

oc_list <- run_oc_calculations(
  p0                 = p0,
  p1                 = p1,
  p_beta_vec         = p_beta_vec,
  type_1_error       = type_1_error,
  power_min          = power_min,
  n_trials_oc        = n_trials_oc,
  method_names       = method_names,
  n_grid_coarse      = n_grid_coarse,
  gamma_grid_coarse  = gamma_grid_coarse,
  n_fine_margin_n    = n_fine_margin_n,
  n_fine_margin_gamma = n_fine_margin_gamma
)

oc_results_coarse <- oc_list$oc_results_coarse
oc_results_fine   <- oc_list$oc_results_fine
oc_all            <- oc_list$oc_all

## -------------------- CREATE PLOTS -------------------------------------------

fig_main <- plot_oc_main(
  oc_all       = oc_all,
  method_names = method_names
)

fig_boundary <- plot_boundary_only(
  oc_results_coarse = oc_results_coarse,
  oc_results_fine   = oc_results_fine,
  method_names      = method_names
)

## Show plots in RStudio Viewer / browser
fig_main
fig_boundary
