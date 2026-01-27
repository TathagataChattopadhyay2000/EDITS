# Tests for getEstimates -------------------------------------------------------

scenarios_list <- simulateScenarios(
  n_subjects_list     = list(c(10, 20, 30)),
  response_rates_list = list(c(0.1, 0.2, 3)),  # 3 triggers historic-rate branch
  n_trials            = 10
)

analyses_list <- performAnalyses(
  scenario_list       = scenarios_list,
  target_rates        = c(0.5, 0.5, 0.5),
  calc_differences    = matrix(c(3, 2, 2, 1), ncol = 2),  # ensures diff cohorts exist
  n_mcmc_iterations   = 100,
  method_names = "berry"
)

outcome <- createTrial(
  n_subjects   = c(10, 20, 30),
  n_responders = c( 1,  2,  3)
)

outcome_analysis <- performAnalyses(
  scenario_list     = outcome,
  target_rates      = c(0.5, 0.5, 0.5),
  n_mcmc_iterations = 100
)


## 1. Baseline functional behaviour & structure  -------------------------------

test_that("getEstimates works for simulated scenarios and returns sensible structure", {
  
  res <- getEstimates(analyses_list)
  
  # function must run and return a list-like structure
  expect_type(res, "list")
  
  expect_true(length(res) > 0)
  
  # Handle both shapes: listPerMethod can give list-of-lists or list-of-matrices
  # TODO seems like it gives list all the time because the name of the method is 
  # also printed so we can skip it
  
  first_obj <- res[[1]]
  
  if (is.list(first_obj)) {
    
    first_mat <- first_obj[[1]]
    
  } else {
    
    first_mat <- first_obj
    
  }
  
  # Inner object must be a matrix of estimates
  expect_true(is.matrix(first_mat))
  
  # For simulated trials (n_trials > 1): expect posterior summaries + Bias + MSE
  expect_true(all(c("Mean", "SD", "2.5%", "50%", "97.5%", "Bias", "MSE") %in% colnames(first_mat)))
  
  # Row names should include:
  #  - p_* cohorts (response rate parameters)
  #  - 'diff' rows from calc_differences matrix
  expect_true(any(grepl("^p_",   rownames(first_mat))))
  
  expect_true(any(grepl("diff", rownames(first_mat))))
})


# 2. Additional parameters + historic-rate handling ----------------------------

test_that("additional parameters are added and have NA bias/MSE", {
  
  # Baseline without additional parameters: used to compare row/column structure
  res_base <- getEstimates(analyses_list)
  
  base_obj <- res_base[[1]]
  
  if (is.list(base_obj)) base_mat <- base_obj[[1]] else base_mat <- base_obj
  
  # Same analyses_list, but now request additional hierarchical parameters
  res_add <- getEstimates(
    
    analyses_list   = analyses_list,
    add_parameters  = c("mu", "tau", "w_1", "w_2", "w_3")
    
  )
  add_obj <- res_add[[1]]
  
  if (is.list(add_obj)) add_mat <- add_obj[[1]] else add_mat <- add_obj
  
  # ---- Historic rate branch check for p_3 ----
  # Row for p_3 must exist (this cohort had response_rate = 3 in simulateScenarios)
  p3_row <- grep("^p_3$", rownames(add_mat))
  
  expect_length(p3_row, 1)
  
  bias_p3   <- add_mat[p3_row, "Bias"]
  
  median_p3 <- add_mat[p3_row, "50%"]
  
  # In scenario_data: response_rates = c(0.1, 0.2, 3) and n_subjects[1,3] = 30
  # getEstimates transforms true_rr[3] to 3 / 30 = 0.1 in the 'historic' branch
  true_rr_raw       <- analyses_list[[1]]$scenario_data$response_rates
  
  n_subj            <- analyses_list[[1]]$scenario_data$n_subjects[1, ]
  
  expected_true_rr3 <- true_rr_raw[3] / n_subj[[3]]
  
  # Bias is defined as: Bias = (point_estimate - true_rr_used)
  # -> point_estimate = Bias + true_rr_used
  point_estimate3 <- bias_p3 + expected_true_rr3
  
  # For median point_estimator (default), the "50%" column should match that estimate
  expect_equal(point_estimate3, median_p3)
  
  # Column structure must be identical with and without add_parameters
  expect_identical(colnames(base_mat), colnames(add_mat))
  
  # Extra rows should correspond exactly to mu/tau/w_* (hierarchical parameters)
  extra_rows <- setdiff(rownames(add_mat), rownames(base_mat))
  
  expect_true(length(extra_rows) > 0)
  
  expect_true(all(grepl("mu|tau|w_", extra_rows)))
  
  # By design: Bias and MSE are *not* computed for non-p_* parameters -> should be NA
  expect_true(all(is.na(add_mat[extra_rows, "Bias"])))
  
  expect_true(all(is.na(add_mat[extra_rows, "MSE"])))
  
})

# 3. Single-trial: only posterior summaries, no Bias/MSE -----------------------

test_that("single-trial outcome returns only posterior summaries (no bias/MSE)", {
  
  # outcome_analysis has a single trial generated with createTrial()
  res_single <- getEstimates(outcome_analysis)
  
  single_obj <- res_single[[1]]
  
  if (is.list(single_obj)) single_mat <- single_obj[[1]] else single_mat <- single_obj
  
  # For a single trial, the function should skip bias/MSE and return only posterior summaries
  expect_true(is.matrix(single_mat))
  
  expect_identical(
    
    colnames(single_mat),
    c("Mean", "SD", "2.5%", "50%", "97.5%")
    
  )
})

# 4. Validation of alpha_level and add_parameters ------------------------------

test_that("alpha_level and add_parameters are validated correctly", {
  
  # (a) alpha_level must be in (0,1) AND match stored quantiles
  #     This hits: numeric constraints + 'stored quantiles' check.
  expect_error(
    
    getEstimates(analyses_list, alpha_level = 0.07)
    
  )
  
  # (b) add_parameters that never appear in any method must trigger the specific error
  expect_error(
    
    getEstimates(analyses_list, add_parameters = c("totally_unknown_param"))
    
  )
})

# 5. point_estimator behaviour and basic input type assertions -----------------

test_that("point_estimator argument is respected and input type is validated", {
  
  # (a) Both 'median' and 'mean' must work and produce matrices of identical shape
  res_median <- getEstimates(analyses_list, point_estimator = "median")
  
  res_mean   <- getEstimates(analyses_list, point_estimator = "mean")
  
  med_obj  <- res_median[[1]]
  
  mean_obj <- res_mean[[1]]
  
  if (is.list(med_obj))  med_mat  <- med_obj[[1]]  else med_mat  <- med_obj
  
  if (is.list(mean_obj)) mean_mat <- mean_obj[[1]] else mean_mat <- mean_obj
  
  expect_identical(dim(med_mat),   dim(mean_mat))
  
  expect_identical(colnames(med_mat), colnames(mean_mat))
  
  # (b) Invalid point_estimator choice -> checkmate::assertChoice should fail
  expect_error(
    
    getEstimates(analyses_list, point_estimator = "mode")

  )
  
  # (c) Non-analysis_list input must fail the class assertion on analyses_list
  expect_error(
    
    getEstimates(list(a = 1))
  
  )
})

# 6. Additional explicit validation tests --------------------------------------
#     (Some overlap with test 4 & 5)

test_that("throws error for invalid alpha_level (outside 0,1)", {
  # This hits the numeric constraint (0 < alpha_level < 1)
  expect_error(getEstimates(analyses_list, alpha_level = 1.5), "alpha_level")
})

test_that("throws error if alpha_level quantiles not available", {
  # This specifically hits: 'alpha_level must be among the stored quantiles'
  expect_error(
    getEstimates(analyses_list, alpha_level = 0.123),
    "must be among the stored quantiles"
  )
})

test_that("throws error if add_parameters not found in any method", {
  # This isolates the 'all(!occurences)' error branch for add_parameters
  expect_error(
    getEstimates(analyses_list, add_parameters = c("nonexistent")),
    "do not occur"
  )
})

test_that("works for single trial outcome (basic check)", {
  # Simpler version of the single-trial test: only structure is checked here
  result <- getEstimates(outcome_analysis)
  expect_type(result, "list")
})

test_that("invalid point_estimator throws error explicitly", {
  # Explicit repeat of the assertChoice branch for point_estimator
  expect_error(
    getEstimates(analyses_list, point_estimator = "mode"),
    "point_estimator"
  )
})

test_that("invalid analyses_list class throws error explicitly", {
  # Explicit repeat of the assertClass(analyses_list, 'analysis_list')
  expect_error(
    getEstimates(list()),
    "analyses_list"
  )
})

# Tests for getGoDecisions -----------------------------------------------------

set.seed(123)

# Use a standard 3-cohort simulation as the main test fixture.

scenarios_list <- simulateScenarios(
  n_subjects_list     = list(c(10, 20, 30)),
  response_rates_list = list(c(0.3, 0.5, 0.7)),
  n_trials            = 10
)

analyses_list <- performAnalyses(
  scenario_list     = scenarios_list,
  target_rates      = c(0.3, 0.3, 0.3),
  n_mcmc_iterations = 100
)

default_cohorts <- c("p_1", "p_2", "p_3")


## 1. Basic argument validation ------------------------------------------------

test_that("getGoDecisions: errors if analyses_list is not of class 'analysis_list'", {
  # exercise checkmate::assertClass(analyses_list, "analysis_list").
  # expect a clear error mentioning 'analysis_list'.
  expect_error(
    getGoDecisions(
      analyses_list   = list(),          # wrong class
      cohort_names    = "p_1",
      evidence_levels = 0.5,
      boundary_rules  = quote(c(TRUE))
    ),
    "analysis_list"
  )
})

test_that("getGoDecisions: errors when evidence_levels or boundary_rules are missing", {
  # triggers the custom simpleError objects for missing evidence_levels / boundary_rules.
  # we expect messages referring to evidence_levels / boundary_rules.
  expect_error(
    getGoDecisions(
      analyses_list = analyses_list,
      cohort_names  = "p_1",
      boundary_rules = quote(c(TRUE))
    ),
    "evidence_levels"
  )
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = "p_1",
      evidence_levels = 0.5
    ),
    "boundary_rules"
  )
})

test_that("getGoDecisions: errors if cohort_names are not valid posterior parameters", {
  # exercises checkmate::assertSubset(cohort_names, choices = ...).
  # an invalid cohort name should cause an immediate error.
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = "invalid",
      evidence_levels = 0.5,
      boundary_rules  = quote(c(TRUE))
    ),
    "cohort_names"
  )
})

test_that("getGoDecisions: errors if overall_min_gos is not a positive integer", {
  # exercises checkmate::assertInt(overall_min_gos, lower = 1).
  # zero or negative overall_min_gos is not allowed.
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = "p_1",
      evidence_levels = 0.5,
      boundary_rules  = quote(c(TRUE)),
      overall_min_gos = 0L
    ),
    "overall_min_gos"
  )
})


## 2. evidence_levels as numeric and as list -----------------------------------

test_that("getGoDecisions: numeric evidence_levels outside (0,1) cause an error", {
  # still in the non-list branch, but with invalid values.
  # we expect check.evidence.levels() to reject them.
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(-0.1, 0.5, 1.3),
      boundary_rules  = quote(c(TRUE, TRUE, TRUE))
    )
  )
})

test_that("getGoDecisions: list evidence_levels all elements valid", {
  # exercises the list(evidence_levels) branch and the for-loop over elements.
  # all elements are valid vectors in (0,1), so we expect no error.
  ev_list_ok <- list(
    rep(0.5, length(default_cohorts)), # here length is 3
    rep(0.5, length(default_cohorts))
  )
  
  expect_silent(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = ev_list_ok,
      boundary_rules  = quote(c(TRUE, TRUE, TRUE))
    )
  )
})

## 3. boundary_rules validation (non-list branch) ------------------------------

test_that("getGoDecisions: non-list boundary_rules must be a language object", {
  # exercises the non-list boundary_rules branch in check_boundary_rules.
  # a plain number is not a language object and should fail is.language().
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = 1  # not language
    ),
    "boundary_rules"
  )
})

test_that("getGoDecisions: non-list boundary_rules must start with c()", {
  # boundary_rules must have head c() (e.g., quote(c(...))).
  # here the call is list(...), so the [1] check against quote(c()) fails.
  expect_error(
  getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = c(0.5, 0.5, 0.8),
    boundary_rules  = quote(list(TRUE, TRUE, TRUE))  # language, but not c()
  ),
  "boundary_rules"
  )
})

test_that("getGoDecisions: non-list boundary_rules must have one entry per cohort", {
  # exercises the non-list length check in check_boundary_rules.
  # with 3 cohorts, c(...) must have length 3, length 2 should trigger an error.
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,           # 3 cohorts
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = quote(c(TRUE, TRUE))       # only 2 entries
    ),
    "boundary_rules"
  )
})

## 4. boundary_rules validation (list branch) ----------------------------------

test_that("getGoDecisions: list boundary_rules each element must be a language object", {
  # Exercises the list(boundary_rules) branch.
  # Non-language element (123) should fail is.language(boundary_rules[[i]]).
  br <- list(123)
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = br
    ),
    "boundary_rules"
  )
})

test_that("getGoDecisions: list boundary_rules each element must start with c()", {
  # In the list case, each element still must be c(...).
  # list(TRUE, TRUE, TRUE) fails the head == quote(c()) check.
  br <- list(quote(list(TRUE, TRUE, TRUE)))
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = br
    ),
    "boundary_rules"
  )
})

test_that("getGoDecisions: list boundary_rules must match number of cohorts", {
  # Exercises list + length check: for 3 cohorts, c(...) must have 3 entries.
  # Here we only provide 2, which should trigger the error.
  br <- list(quote(c(TRUE, TRUE)))
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = br
    ),
    "boundary_rules"
  )
})

## 5. Length consistency: method_names vs boundary_rules / gamma_levels --------

test_that("errors if boundary_rules list is longer than method_names", {
  
  m    <- analyses_list[[1]]$analysis_parameters$method_names
  coh  <- c("p_1", "p_2", "p_3")
  ev   <- c(0.5, 0.5, 0.5)
  
  # list longer than number of methods
  br   <- rep(list(quote(c(TRUE, TRUE, TRUE))), length(m) + 1L)
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = coh,
      evidence_levels = ev,
      boundary_rules  = br
    ),
    "boundary_rules"
  )
})


test_that("errors if evidence_levels list is longer than method_names", {
  
  m    <- analyses_list[[1]]$analysis_parameters$method_names
  coh  <- c("p_1", "p_2", "p_3")
  
  # 0.5 is safe (1 - 0.5 = 0.5 in stored quantiles)
  ev_vec <- c(0.5, 0.5, 0.5)
  ev     <- rep(list(ev_vec), length(m) + 1L)
  
  expect_error(
    getGoDecisions(
      analyses_list   = analyses_list,
      cohort_names    = coh,
      evidence_levels = ev,
      boundary_rules  = quote(c(TRUE, TRUE, TRUE))
    ),
    "evidence_levels"
  )
})


test_that("single boundary_rules expression is recycled to all methods", {
  
  m    <- analyses_list[[1]]$analysis_parameters$method_names
  coh  <- c("p_1", "p_2", "p_3")
  ev   <- c(0.5, 0.5, 0.5)
  rule <- quote(c(TRUE, TRUE, TRUE))
  
  dec <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = coh,
    evidence_levels = ev,
    boundary_rules  = rule
  )
  
  br <- dec[[1]]$decision_rules$boundary_rules
  
  expect_identical(length(br), length(m))
  for (i in seq_along(br)) {
    expect_true(identical(br[[i]], rule))
  }
})


test_that("single evidence_levels vector is recycled to all methods", {
  
  m    <- analyses_list[[1]]$analysis_parameters$method_names
  coh  <- c("p_1", "p_2", "p_3")
  ev   <- c(0.5, 0.5, 0.5)
  
  dec <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = coh,
    evidence_levels = ev,                # not a list -> will be wrapped + recycled
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  gamma <- dec[[1]]$decision_rules$gamma_levels
  
  expect_identical(length(gamma), length(m))
  for (i in seq_along(gamma)) {
    expect_equal(gamma[[i]], ev)
  }
})


## 6. Method names consistency across scenarios --------------------------------

test_that("getGoDecisions: succeeds when all scenarios use identical method_names", {
  # If all scenarios use the same method_names, no error is thrown.
  
  res <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = c(0.5, 0.5, 0.8),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  expect_s3_class(res, "decision_list")
})

test_that("getGoDecisions: errors when scenarios were analysed with different methods", {
  # We clone scenario_1 as scenario_2 and then reverse its method_names.
  
  bad <- analyses_list
  bad$scenario_2 <- bad$scenario_1
  bad$scenario_2$analysis_parameters$method_names <-
    rev(bad$scenario_2$analysis_parameters$method_names)
  
  expect_error(
    getGoDecisions(
      analyses_list   = bad,
      cohort_names    = default_cohorts,
      evidence_levels = c(0.5, 0.5, 0.8),
      boundary_rules  = quote(c(TRUE, TRUE, TRUE))
    ),
    "analysed with different methods"
  )
})


## 7. Core structure: decision_list and 'overall' column -----------------------

test_that("getGoDecisions: returns decision_list with overall and cohort decisions", {
  # Exercises the full decision-generation code path for a valid case.
  # We inspect the structure of one scenario/method for basic sanity.
  decisions <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = c(0.5, 0.5, 0.8),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  scen1  <- decisions[[1]]
  m_dec  <- as.matrix(scen1$decisions_list[[1]])
  
  # Expect one 'overall' column and at least one cohort-level decision column.
  expect_true("overall" %in% colnames(m_dec))
  cohort_cols <- setdiff(colnames(m_dec), "overall")
  expect_true(length(cohort_cols) >= 1)
  
  # decision_rules should remember the cohort_names we passed in.
  expect_identical(scen1$decision_rules$cohort_names, default_cohorts)
  
  # gamma_levels may be a vector or list; flatten and check values are present.
  stored_gamma_flat <- unlist(scen1$decision_rules$gamma_levels, use.names = FALSE)
  expect_true(all(c(0.5, 0.5, 0.8) %in% stored_gamma_flat))
})


## 8. Semantics: overall_min_gos behaviour -------------------------------------

test_that("overall_min_gos = 1 means overall Go if at least one cohort is Go", {
  
  dec <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = rep(0.5, length(default_cohorts)),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE)),  # simple always-true rule
    overall_min_gos = 1L
  )
  
  m <- as.matrix(dec[[1]]$decisions_list[[1]])
  
  # All columns except 'overall' are cohort-level decisions.
  coh_cols <- setdiff(colnames(m), "overall")
  expect_true(length(coh_cols) >= 1)
  
  coh <- m[, coh_cols, drop = FALSE] > 0   # logical matrix of cohort decisions
  
  # Definition: overall is TRUE if at least one cohort is TRUE.
  overall_calc <- apply(coh, 1, function(x) sum(x) >= 1L)
  
  overall <- as.logical(m[, "overall"])
  expect_identical(overall, overall_calc)
})


test_that("overall_min_gos = 2 means overall Go if at least two cohorts are Go", {
  
  dec <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = rep(0.5, length(default_cohorts)),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE)),
    overall_min_gos = 2L
  )
  
  m <- as.matrix(dec[[1]]$decisions_list[[1]])
  
  coh_cols <- setdiff(colnames(m), "overall")
  expect_true(length(coh_cols) >= 1)
  
  coh <- m[, coh_cols, drop = FALSE] > 0
  
  # Definition: overall is TRUE if at least two cohorts are TRUE.
  overall_calc <- apply(coh, 1, function(x) sum(x) >= 2L)
  
  overall <- as.logical(m[, "overall"])
  expect_identical(overall, overall_calc)
})

## 9. There should not be a case where a cohort was NoGo before and go now -----

test_that("getGoDecisions: previous go_decisions prevent resurrection of stopped cohorts", {
  # This test exercises 
  #   go_decisions <- go_decisions * previous_gos > 0
  # It ensures that once a cohort was no-go previously, it cannot become go now.
  if (is.null(analyses_list[[1]]$scenario_data$previous_analyses)) {
    skip("previous_analyses not available in scenario_data")
  }
  
  decisions <- getGoDecisions(
    analyses_list   = analyses_list,
    cohort_names    = default_cohorts,
    evidence_levels = c(0.5, 0.5, 0.8),
    boundary_rules  = quote(c(TRUE, TRUE, TRUE))
  )
  
  new_mat     <- as.matrix(decisions[[1]]$decisions_list[[1]])
  new_cohcols <- setdiff(colnames(new_mat), "overall")
  new_gos     <- new_mat[, new_cohcols, drop = FALSE] > 0
  
  prev_mat <- analyses_list[[1]]$scenario_data$previous_analyses$go_decisions
  # Drop first column (overall) to get per-cohort previous decisions.
  prev_gos <- as.matrix(prev_mat[, -1, drop = FALSE]) > 0
  
  # Ensure matrices are conformable for comparison.
  expect_identical(dim(prev_gos), dim(new_gos))
  
  # A 'resurrected' decision would be TRUE now but FALSE previously.
  resurrected <- new_gos & !prev_gos
  expect_false(any(resurrected))
})

# Tests for getGoProbabilities -------------------------------------------------

set.seed(456)

# We simulate a simple 2-cohort scenario and run getGoDecisions() twice:
#   - once for Go rules
#   - once for NoGo rules
# These two decision_list objects are then reused across all tests.

scenarios_list_prob <- simulateScenarios(
  n_subjects_list     = list(c(10, 20)),
  response_rates_list = list(rep(0.9, 2)),
  n_trials            = 10
)

analyses_list_prob <- performAnalyses(
  scenario_list     = scenarios_list_prob,
  target_rates      = rep(0.5, 2),
  n_mcmc_iterations = 100
)

prob_cohorts <- c("p_1", "p_2")

go_decisions_list <- getGoDecisions(
  analyses_list   = analyses_list_prob,
  cohort_names    = prob_cohorts,
  evidence_levels = c(0.5, 0.5),
  boundary_rules  = quote(c(x[1] > 0.8, x[2] > 0.6))
)

nogo_decisions_list <- getGoDecisions(
  analyses_list   = analyses_list_prob,
  cohort_names    = prob_cohorts,
  evidence_levels = c(0.5, 0.5),
  boundary_rules  = quote(c(x[1] < 0.5, x[2] < 0.3))
)


## 1. Basic argument validation ------------------------------------------------

test_that("getGoProbabilities: errors if go_decisions_list is not a decision_list", {
  # This hits checkmate::assertClass(go_decisions_list, "decision_list").
  # Passing a plain list is invalid and should trigger an informative error.
  
  expect_error(
    
    getGoProbabilities(go_decisions_list = list()),
    "go_decisions_list"
    
  )
})

test_that("getGoProbabilities: errors if nogo_decisions_list has wrong class", {
  # This hits checkmate::assertClass(nogo_decisions_list, ..., null.ok = TRUE).
  # go_decisions_list is OK, but nogo_decisions_list is a raw list -> error.
  
  expect_error(
    
    getGoProbabilities(
      go_decisions_list   = go_decisions_list,
      nogo_decisions_list = list()
    ),
    "decision_list"
    
  )
})

test_that("getGoProbabilities: errors if Go and NoGo matrices have different dimensions", {
  # This targets the explicit dimension check:
  #   identical(dim(go_decisions_list[[1]]$decisions_list[[1]]),
  #             dim(nogo_decisions_list[[1]]$decisions_list[[1]]))
  # We deliberately break the dimensions of the first NoGo matrix.
  
  bad_nogo <- nogo_decisions_list
  bad_nogo[[1]]$decisions_list[[1]] <-
    bad_nogo[[1]]$decisions_list[[1]][, -1, drop = FALSE]  # drop one column
  
  expect_error(
    
    getGoProbabilities(
      go_decisions_list   = go_decisions_list,
      nogo_decisions_list = bad_nogo
    ),
  )
})


## 2. Go-only case: basic structure and semantics ------------------------------

test_that("getGoProbabilities: Go-only call returns list-of-lists of matrices", {
  # This exercises the main loop without nogo_decisions_list:
  #   - only "Go" row is created
  #   - listPerMethod is applied at the end
  
  probs <- getGoProbabilities(go_decisions_list)
  
  # Top-level: one element per method (after listPerMethod)
  expect_type(probs, "list")
  expect_true(all(sapply(probs, is.list)))
  
  # Inspect first method, first scenario
  first_method   <- names(probs)[1]
  first_scenario <- names(probs[[first_method]])[1]
  mat            <- probs[[first_method]][[first_scenario]]
  
  # In Go-only mode, we expect a single "Go" row and ≥1 column (overall + cohorts)
  expect_true(is.matrix(mat))
  expect_identical(rownames(mat), "Go")
  expect_true(ncol(mat) >= 1)
})

test_that("getGoProbabilities: Go row equals column means of go_decisions", {
  # This test directly checks the core numeric logic:
  # decisions_matrix <- t(as.matrix(colMeans(go_decisions)))
  
  probs <- getGoProbabilities(go_decisions_list)
  
  # Method and scenario names derived from the go_decisions_list structure
  method_names   <- names(go_decisions_list[[1]]$decisions_list)
  scenario_names <- names(go_decisions_list)
  
  for (m in method_names) {
    for (s in scenario_names) {
      mat_prob <- probs[[m]][[s]]
      go_mat   <- go_decisions_list[[s]]$decisions_list[[m]]
      
      # Row "Go" should match the column means of the original go_decisions matrix
      expected <- colMeans(go_mat)
      expect_equal(
        mat_prob["Go", ],
        expected,
        info = paste("Mismatch in method", m, "scenario", s)
      )
    }
  }
})


## 3. Go + NoGo case: structure and probability --------------------------------

test_that("getGoProbabilities: Go and NoGo rows match colMeans of input decisions", {
  # This checks that:
  #   - "Go" row = colMeans(go_decisions)
  #   - "NoGo" row = colMeans(nogo_decisions)
  
  probs <- getGoProbabilities(
    go_decisions_list   = go_decisions_list,
    nogo_decisions_list = nogo_decisions_list
  )
  
  method_names   <- names(go_decisions_list[[1]]$decisions_list)
  scenario_names <- names(go_decisions_list)
  
  for (m in method_names) {
    for (s in scenario_names) {
      mat_prob <- probs[[m]][[s]]
      
      go_mat   <- go_decisions_list[[s]]$decisions_list[[m]]
      nogo_mat <- nogo_decisions_list[[s]]$decisions_list[[m]]
      
      expected_go   <- colMeans(go_mat)
      expected_nogo <- colMeans(nogo_mat)
      
      expect_equal(
        mat_prob["Go", ],
        expected_go,
        info = paste("Go row mismatch in method", m, "scenario", s)
      )
      
      expect_equal(
        mat_prob["NoGo", ],
        expected_nogo,
        info = paste("NoGo row mismatch in method", m, "scenario", s)
      )
    }
  }
})

test_that("getGoProbabilities: columns sum to 1 when NoGo decisions are provided", {
  # This checks the probability normalization:
  #   Consider = round(1 - p(Go) - p(NoGo), 9)
  # so row sums per column should be ~1.
  
  probs <- getGoProbabilities(
    go_decisions_list   = go_decisions_list,
    nogo_decisions_list = nogo_decisions_list
  )
  
  first_method <- names(probs)[1]
  for (scen_name in names(probs[[first_method]])) {
    mat <- probs[[first_method]][[scen_name]]
    
    col_sums <- colSums(mat)
    expect_true(
      all(abs(col_sums - 1) == 0),
      info = paste("Column sums not 1 for scenario", scen_name)
    )
  }
})


## 4. Overlap check: Go and NoGo cannot both be TRUE for the same cell ---------

test_that("getGoProbabilities: errors if any decision is both Go and NoGo", {
  # This directly hits:
  #   if (!isTRUE(all.equal(sum(go_decisions * nogo_decisions), 0))) stop(...)
  #
  # We artificially modify one scenario/method so that at least one entry
  # is TRUE in both go_decisions and nogo_decisions, forcing the error.
  overlap_go   <- go_decisions_list
  overlap_nogo <- nogo_decisions_list
  
  # Pick scenario_1, first method, first column (arbitrary but consistent).
  scen1_name <- names(overlap_go)[1]
  meth1_name <- names(overlap_go[[scen1_name]]$decisions_list)[1]
  
  overlap_go[[scen1_name]]$decisions_list[[meth1_name]][1, 1]   <- TRUE
  overlap_nogo[[scen1_name]]$decisions_list[[meth1_name]][1, 1] <- TRUE
  
  expect_error(
    getGoProbabilities(
      go_decisions_list   = overlap_go,
      nogo_decisions_list = overlap_nogo
    ),
    "both go and nogo decisions"
  )
})

# Tests for print.decision.list ------------------------------------------------

set.seed(789)

# We build a decision_list with:
# - 2 scenarios (so we can see pluralisation in the header)
# - potentially multiple methods (whatever performAnalyses uses)
scenarios_print <- simulateScenarios(
  n_subjects_list     = list(c(10, 20), c(10, 20)),
  response_rates_list = list(c(0.3, 0.5), c(0.4, 0.6)),
  n_trials            = 5
)

analyses_print <- performAnalyses(
  scenario_list     = scenarios_print,
  target_rates      = c(0.3, 0.3),
  n_mcmc_iterations = 80
)

coh_print <- c("p_1", "p_2")

decisions_dl <- getGoDecisions(
  analyses_list   = analyses_print,
  cohort_names    = coh_print,
  evidence_levels = c(0.5, 0.5),
  boundary_rules  = quote(c(TRUE, TRUE))
)

# Convenience values for expected header
n_scenarios  <- length(decisions_dl)
n_methods    <- length(decisions_dl[[1]]$decisions_list)


## 1. Header line: scenario/method counts, singular/plural ---------------------

test_that("print.decision_list: header shows correct number of scenarios and methods", {
  # This hits the initial cat() that prints:
  #   "decision_list of <N> scenario(s) with <M> method(s)"
  # We build the exact expected phrase based on our fixture.
  header_pattern <- sprintf(
    "decision_list of %d scenario%s with %d method%s",
    n_scenarios,
    ifelse(n_scenarios == 1, "", "s"),
    n_methods,
    ifelse(n_methods == 1, "", "s")
  )
  
  expect_output(
    print(decisions_dl),
    header_pattern
  )
})

## 2. Scenario sections: each scenario name is printed -------------------------

test_that("print.decision_list: prints a section for each scenario", {
  # After the header, the function prints:
  #   "  - <scenario_name>"
  # for each scenario in x.
  scen_names <- names(decisions_dl)
  expect_true(length(scen_names) >= 1)
  
  output <- capture.output(print(decisions_dl))
  
  # Every scenario name should appear at least once in the output.
  for (nm in scen_names) {
    expect_true(
      any(grepl(paste0("\\b", nm, "\\b"), output)),
      info = paste("Scenario name", nm, "not found in printed output")
    )
  }
})

## 3. Method labels: each method appears in the row labels ---------------------

test_that("print.decision_list: each method name appears in printed matrix rows", {
  # print.decision_list builds rownames from method_names using firstUpper() and some padding.
  # We don't check the exact spacing, but we do require that the method names are visible.
  method_names <- names(decisions_dl[[1]]$decisions_list)
  expect_true(length(method_names) >= 1)
  
  output <- capture.output(print(decisions_dl))
  
  # Each method name (case-insensitive) should appear at least once.
  for (m in method_names) {
    # firstUpper(m) is what's actually used, so we mimic that in the search
    m_upper <- paste0(toupper(substr(m, 1, 1)), substr(m, 2, nchar(m)))
    expect_true(
      any(grepl(m_upper, output, fixed = TRUE)),
      info = paste("Method label", m_upper, "not found in printed output")
    )
  }
})

## 4. Digits argument: different digits should still print without error -------

test_that("print.decision_list: digits argument is accepted and does not error", {
  # This hits the round(mat_out, digits = digits) call.
  # We don't try to parse exact numeric formatting; we just assert that
  # different digits values are accepted and produce some usable output.
  expect_output(
    print(decisions_dl, digits = 1),
    "decision_list of"
  )
  
  expect_output(
    print(decisions_dl, digits = 4),
    "decision_list of"
  )
})

## 5. decision_rules branch: print works when decision_rules are present -------

test_that("print.decision_list: handles non-NULL decision_rules without error", {
  # This test specifically exercises the big 'if (!is.null(x[[1]]$decision_rules))' block:
  # - It rewrites boundary_rules expressions
  # - Uses decision_rules$cohort_names and decision_rules$gamma_levels[[n]]
  #
  # Any bug in that transformation would show up as an error here.
  expect_output(
    print(decisions_dl),
    "decision_list of"
  )
})

## 6. decision_rules branch: print also works when decision_rules are NULL -----

test_that("print.decision_list: works when decision_rules are NULL", {
  # Here we remove decision_rules from the decision_list to ensure
  # the function safely skips the formatting block and still prints the summary.
  dl_no_rules <- decisions_dl
  for (i in seq_along(dl_no_rules)) {
    dl_no_rules[[i]]$decision_rules <- NULL
  }
  
  # Expect the same kind of header, and no error.
  header_pattern <- sprintf(
    "decision_list of %d scenario%s with %d method%s",
    n_scenarios,
    ifelse(n_scenarios == 1, "", "s"),
    n_methods,
    ifelse(n_methods == 1, "", "s")
  )
  
  expect_output(
    print(dl_no_rules),
    header_pattern
  )
})

## 7. Multi-method aggregation: rows for each method from getGoProbabilities ---

test_that("print.decision_list: for multiple methods, one row per method is printed per scenario", {
  # This test uses:
  #   go_probs <- getGoProbabilities(x)
  #   mat_out <- do.call(rbind, lapply(go_probs, function(y) y[[n]]))
  #
  # For each scenario, mat_out has one row per method.
  # print.decision_list loops over *all scenarios* and prints that block each time,
  # so across the full output we expect:
  #   total method label lines = n_scenarios * n_methods.
  go_probs <- getGoProbabilities(decisions_dl)
  
  # For scenario 1, reconstruct mat_out as in print.decision_list
  scenario_index <- 1L
  mat_ref <- do.call(
    rbind,
    lapply(go_probs, function(y) y[[scenario_index]])
  )
  
  output <- capture.output(print(decisions_dl))
  
  # Lines starting with "    - " are the method label rows we created as rownames.
  method_label_lines <- grep("^    - ", output, value = TRUE)
  
  expected_total_rows <- n_scenarios * nrow(mat_ref)
  
  expect_equal(
    length(method_label_lines),
    expected_total_rows,
    info = paste(
      "Number of printed method label lines (", length(method_label_lines),
      ") does not match expected total rows (", expected_total_rows, ")"
    )
  )
})

# Tests for getAverageNSubjects ------------------------------------------------

scenarios_list <- simulateScenarios(
  n_subjects_list     = list(c(10, 20, 30), c(15, 25, 35)),
  response_rates_list = list(c(0.1, 0.2, 0.3), c(0.15, 0.25, 0.35)),
  n_trials            = 2
)

test_that("returns correct structure for multiple scenarios", {
  result <- getAverageNSubjects(scenarios_list)
  expect_type(result, "list")
  expect_equal(length(result), length(scenarios_list))
  expect_true(all(sapply(result, is.numeric)))
})

test_that("computes correct column means", {
  # Manually compute expected means for first scenario
  expected <- colMeans(scenarios_list[[1]]$n_subjects)
  result <- getAverageNSubjects(scenarios_list)
  expect_equal(result[[1]], expected)
})

test_that("works for single scenario", {
  single_scenario <- simulateScenarios(
    n_subjects_list     = list(c(10, 20, 30)),
    response_rates_list = list(c(0.1, 0.2, 0.3)),
    n_trials            = 3
  )
  result <- getAverageNSubjects(single_scenario)
  expect_equal(names(result), names(single_scenario))
  expect_equal(result[[1]], colMeans(single_scenario[[1]]$n_subjects))
})

test_that("throws error for invalid input class", {
  expect_error(getAverageNSubjects(list()), "scenario_list")
})

test_that("returns empty list for empty scenario_list", {
  empty_list <- structure(list(), class = "scenario_list")
  result <- getAverageNSubjects(empty_list)
  expect_equal(result, list())
})

# Test for negateGoDecisions ---------------------------------------------------

set.seed(101)

#  - 2 cohorts
#  - several trials

scenarios_neg <- simulateScenarios(
  n_subjects_list     = list(c(10, 15)),
  response_rates_list = list(c(0.4, 0.7)),
  n_trials            = 6
)

analyses_neg <- performAnalyses(
  scenario_list     = scenarios_neg,
  target_rates      = c(0.4, 0.4),
  n_mcmc_iterations = 80
)

coh_neg <- c("p_1", "p_2")

go_decisions_list_neg <- getGoDecisions(
  analyses_list   = analyses_neg,
  cohort_names    = coh_neg,
  evidence_levels = c(0.5, 0.5),            
  boundary_rules  = quote(c(TRUE, TRUE))    # always-go rule for simplicity
)

n_scen_neg <- length(go_decisions_list_neg)

n_meth_neg <- length(go_decisions_list_neg[[1]]$decisions_list)

## 1. Basic argument validation ------------------------------------------------

test_that("negateGoDecisions: go_decisions_list must be a decision_list", {
  # Hits assertClass(go_decisions_list, 'decision_list').
  expect_error(
    negateGoDecisions(go_decisions_list = list()),
    "go_decisions_list"
  )
})

test_that("negateGoDecisions: overall_min_nogos must be 'all' or non-negative integer", {
  # Hits the combined assertion:
  #   checkChoice('all') OR checkInt(lower = 0).
  expect_error(
    negateGoDecisions(
      go_decisions_list = go_decisions_list_neg,
      overall_min_nogos = -1L
    ),
    "overall_min_nogos"
  )
  
  expect_error(
    negateGoDecisions(
      go_decisions_list = go_decisions_list_neg,
      overall_min_nogos = "foo"
    ),
    "overall_min_nogos"
  )
})

## 2. negation of cohort-level decisions ---------------------------------------


test_that("negateGoDecisions: cohort-level entries are logically negated", {
  # Checks:
  #   nogo_decisions_list[[s]]$decisions_list[[m]] <- !go_decisions_list[[s]]$decisions_list[[m]]
  # then ignores the overwritten overall column and only checks cohort columns.
  
  nogo_list <- negateGoDecisions(
    go_decisions_list = go_decisions_list_neg,
    overall_min_nogos = "all"
  )
  
  for (s in seq_len(n_scen_neg)) {
    for (m in seq_len(n_meth_neg)) {
      go_mat   <- go_decisions_list_neg[[s]]$decisions_list[[m]]
      nogo_mat <- nogo_list[[s]]$decisions_list[[m]]
      
      if (ncol(go_mat) > 1L) {
        go_coh   <- go_mat[, -1, drop = FALSE]
        nogo_coh <- nogo_mat[, -1, drop = FALSE]
        
        expect_identical(
          nogo_coh,
          !go_coh,
          info = paste("Cohort decisions not negated in scenario", s, "method", m)
        )
      }
    }
  }
})


## 3. overall_min_nogos = 'all' ------------------------------------------------

test_that("negateGoDecisions: overall_min_nogos = 'all' means 'all cohorts NoGo'", {
  # Hits:
  #   overall_min_nogos_org <- "all" => overall_min_nogos <- n_decisions - 1L
  # Overall column is TRUE iff all cohort columns are TRUE.
  
  nogo_list <- negateGoDecisions(
    go_decisions_list = go_decisions_list_neg,
    overall_min_nogos = "all"
  )
  
  for (s in seq_len(n_scen_neg)) {
    for (m in seq_len(n_meth_neg)) {
      nogo_mat <- nogo_list[[s]]$decisions_list[[m]]
      
      if (ncol(nogo_mat) <= 1L) next
      
      n_decisions   <- ncol(nogo_mat)
      nogo_coh      <- nogo_mat[, -1, drop = FALSE]
      expected_over <- apply(
        nogo_coh, 1,
        function(x) sum(x) >= (n_decisions - 1L)
      )
      actual_over   <- as.logical(nogo_mat[, 1])
      
      expect_identical(
        actual_over,
        expected_over
      )
    }
  }
})


## 4. Numeric overall_min_nogos ------------------------------------------------

test_that("negateGoDecisions: overall_min_nogos = 0 makes overall always TRUE", {
  # With threshold 0, apply(sum(x) >= 0) is always TRUE.
  nogo_list <- negateGoDecisions(
    go_decisions_list = go_decisions_list_neg,
    overall_min_nogos = 0L
  )
  
  for (s in seq_len(n_scen_neg)) {
    for (m in seq_len(n_meth_neg)) {
      nogo_mat <- nogo_list[[s]]$decisions_list[[m]]
      
      if (ncol(nogo_mat) <= 1L) next
      
      overall_col <- as.logical(nogo_mat[, 1])
      expect_true(
        all(overall_col)
      )
    }
  }
})

test_that("negateGoDecisions: numeric overall_min_nogos = k means 'at least k cohorts NoGo'", {
  # For k = 1, overall is TRUE iff at least one cohort column is TRUE.
  k <- 1L
  nogo_list <- negateGoDecisions(
    go_decisions_list = go_decisions_list_neg,
    overall_min_nogos = k
  )
  
  for (s in seq_len(n_scen_neg)) {
    for (m in seq_len(n_meth_neg)) {
      nogo_mat <- nogo_list[[s]]$decisions_list[[m]]
      if (ncol(nogo_mat) <= 1L) next
      
      nogo_coh      <- nogo_mat[, -1, drop = FALSE]
      expected_over <- apply(nogo_coh, 1, function(x) sum(x) >= k)
      actual_over   <- as.logical(nogo_mat[, 1])
      
      expect_identical(
        actual_over,
        expected_over,
        info = paste("Overall NoGo != '>= k cohorts' for k =", k,
                     "in scenario", s, "method", m)
      )
    }
  }
})

## 5. Edge case: matrices with only one column ---------------------------------

test_that("negateGoDecisions: single-column decisions are just negated, no overall logic", {
  # Hits the n_decisions <= 1 branch: overall column is not recomputed.
  simple_dl <- list(
    scenario_1 = list(
      decisions_list = list(
        method_1 = matrix(c(TRUE, FALSE, TRUE), ncol = 1)
      )
    )
  )
  class(simple_dl) <- "decision_list"
  
  nogo_simple <- negateGoDecisions(
    go_decisions_list = simple_dl,
    overall_min_nogos = "all"
  )
  
  orig_mat <- simple_dl$scenario_1$decisions_list$method_1
  new_mat  <- nogo_simple$scenario_1$decisions_list$method_1
  
  expect_identical(
    new_mat,
    !orig_mat,
    info = "Single-column decision matrix not fully negated"
  )
})
