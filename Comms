# Tests for getMuVar -----------------------------------------------------------

test_that("mu_var calculation is correct for simple case", {
  
  # Manual calculation: (n_worth * p * (1 - p))^-1 - tau_scale^2
  expected <- (2 * 0.5 * (1 - 0.5))^-1 - 1^2
  
  result <- getMuVar(response_rate = 0.5, tau_scale = 1, n_worth = 2)
  
  expect_equal(result, expected)
  
  expect_type(result, "double")
  
  expect_true(is.finite(result))
  
})

test_that("error if response_rate out of bounds", {
  
  expect_error(getMuVar(response_rate = -0.1, tau_scale = 1))
  
  expect_error(getMuVar(response_rate = 1.5, tau_scale = 1))
  
})

test_that("error if tau_scale negative", {
  
  expect_error(getMuVar(response_rate = 0.5, tau_scale = -1))
  
})

test_that("error if n_worth not positive integer", {
  
  expect_error(getMuVar(response_rate = 0.5, tau_scale = 1, n_worth = 0))
  
  expect_error(getMuVar(response_rate = 0.5, tau_scale = 1, n_worth = 1.5))
  
})

# Tests for getPriorParametersBerry --------------------------------------------

test_that("valid input returns prior_parameters_list with correct structure", {
  result <- getPriorParametersBerry(target_rates = c(0.2, 0.8), tau_scale = 1, n_worth = 2)
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "berry")
  
  # Check inner structure
  expect_true(is.list(result$berry))
  expect_named(result$berry, c("mu_mean", "mu_sd", "tau_scale"))
})

test_that("mu_sd is computed correctly for simple case", {
  target_rates <- c(0.2, 0.8)
  tau_scale <- 1
  n_worth <- 2
  
  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]
  expected_mu_var <- (n_worth * target_rate_max_var * (1 - target_rate_max_var))^-1 - tau_scale^2
  expected_mu_sd <- sqrt(expected_mu_var)
  
  result <- getPriorParametersBerry(target_rates, tau_scale, n_worth)
  expect_equal(result$berry$mu_sd, expected_mu_sd)
})

test_that("error if mu_var <= 0", {
  # Large tau_scale will make mu_var negative
  expect_error(
    getPriorParametersBerry(target_rates = c(0.5), tau_scale = 100, n_worth = 1)
  )
})

test_that("single target_rate works", {
  result <- getPriorParametersBerry(target_rates = c(0.5), tau_scale = 1, n_worth = 1)
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "berry")
})


# Tests for setPriorParametersBerry --------------------------------------------

test_that("valid input returns prior_parameters_list with correct structure", {
  result <- setPriorParametersBerry(
    mu_mean = 0.1,
    mu_sd = 0.5,
    tau_scale = 1
  )
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "berry")
  expect_true(is.list(result$berry))
  expect_named(result$berry, c("mu_mean", "mu_sd", "tau_scale"))
})

test_that("error if mu_sd is non-positive", {
  expect_error(
    setPriorParametersBerry(mu_mean = 0.1, mu_sd = 0, tau_scale = 1),
    "mu_sd"
  )
})

test_that("error if tau_scale is non-positive", {
  expect_error(
    setPriorParametersBerry(mu_mean = 0.1, mu_sd = 0.5, tau_scale = 0),
    "tau_scale"
  )
})

test_that("error if mu_mean is not numeric", {
  expect_error(
    setPriorParametersBerry(mu_mean = "a", mu_sd = 0.5, tau_scale = 1),
    "mu_mean"
  )
})

# Tests for getPriorParametersExNeX --------------------------------------------

test_that("valid input returns prior_parameters_list with correct structure", {
  result <- getPriorParametersExNex(target_rates = c(0.3, 0.9), tau_scale = 1, n_worth = 2, w_j = 0.5)
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "exnex")
  
  # Check inner structure
  expect_true(is.list(result$exnex))
  expect_named(result$exnex, c("mu_mean", "mu_sd", "tau_scale", "mu_j", "tau_j", "w_j"))
})

test_that("mu_mean and mu_j use logit transform", {
  target_rates <- 0.8
  tau_scale <- 1
  n_worth <- 2
  w_j <- 0.5
  
  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]
  expected_mu_mean <- logit(target_rate_max_var)
  expected_mu_j <- logit(target_rates)
  result <- getPriorParametersExNex(target_rates, tau_scale, n_worth, w_j)
  
  expect_equal(result$exnex$mu_mean, expected_mu_mean)
  expect_equal(result$exnex$mu_j, expected_mu_j)
})

test_that("error if mu_var <= 0", {
  # Large tau_scale will make mu_var negative
  expect_error(
    getPriorParametersExNex(target_rates = c(0.5), tau_scale = 100, n_worth = 1)
  )
})

test_that("mu_sd and tau_j are computed correctly", {
  target_rates <- 0.8
  tau_scale <- 1
  n_worth <- 2
  
  target_rate_max_var <- target_rates[abs(target_rates - 0.5) == max(abs(target_rates - 0.5))][1]
  expected_mu_var <- (n_worth * target_rate_max_var * (1 - target_rate_max_var))^-1 - tau_scale^2
  expected_mu_sd <- sqrt(expected_mu_var)
  
  expected_tau_j <- sqrt((n_worth * target_rates * (1 - target_rates))^-1)
  
  result <- getPriorParametersExNex(target_rates, tau_scale, n_worth, w_j = 0.5)
  
  expect_equal(result$exnex$mu_sd, expected_mu_sd)
  expect_equal(result$exnex$tau_j, expected_tau_j)
})

# Tests for setPriorParametersExNeX --------------------------------------------

test_that("valid input returns prior_parameters_list with correct structure", {
  result <- setPriorParametersExNex(
    mu_mean = c(0.1, 0.2),
    mu_sd = c(0.5, 0.6),
    tau_scale = 1,
    mu_j = c(0.1, 0.2),
    tau_j = c(0.3, 0.4),
    w_j = c(0.3, 0.3, 0.4)
  )

  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "exnex")
  expect_true(is.list(result$exnex))
  expect_named(result$exnex, c("mu_mean", "mu_sd", "tau_scale", "mu_j", "tau_j", "w_j"))
})


test_that("error if w_j length is invalid when mu_mean length = 1", {
  expect_error(
    setPriorParametersExNex(
      mu_mean = c(0.1),
      mu_sd = c(0.5),
      tau_scale = 1,
      mu_j = c(0.1),
      tau_j = c(0.3),
      w_j = c(0.3, 0.3, 0.4)  # length 3 is invalid for mu_mean length 1
    )
  )
})

# Tests for getPriorParametersExNeXAdj -----------------------------------------

test_that("valid input returns prior_parameters_list with correct structure", {
  result <- getPriorParametersExNexAdj(
    target_rates = c(0.2, 0.3),
    tau_scale = 1,
    n_worth = 2,
    w_j = 0.5
  )
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "exnex_adj")
  expect_true(is.list(result$exnex_adj))
  expect_named(result$exnex_adj, c("mu_mean", "mu_sd", "tau_scale", "mu_j", "tau_j", "w_j"))
  
  # Check adjusted values
  expect_equal(result$exnex_adj$mu_mean, 0)
  expect_equal(result$exnex_adj$mu_j, rep(0, length(c(0.2, 0.3))))
})

# Tests for setPriorParametersExNeXAdj -----------------------------------------

test_that("valid input returns prior_parameters_list with correct structure", {
  result <- setPriorParametersExNexAdj(
    mu_mean = c(0.1, 0.2),
    mu_sd = c(0.5, 0.6),
    tau_scale = 1,
    mu_j = c(0.1, 0.2),
    tau_j = c(0.3, 0.4),
    w_j = c(0.3, 0.3, 0.4)
  )
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "exnex_adj")
  expect_true(is.list(result$exnex_adj))
  expect_named(result$exnex_adj, c("mu_mean", "mu_sd", "tau_scale", "mu_j", "tau_j", "w_j"))
})


test_that("'w_j' length rule works when mu_mean length = 1", {
  # Valid cases: w_j length 1 or 2
  expect_silent(
    setPriorParametersExNexAdj(
      mu_mean = 0.1,
      mu_sd = 0.5,
      tau_scale = 1,
      mu_j = c(0.1),
      tau_j = c(0.3),
      w_j = c(0.5)  # length 1
    )
  )
})
# Tests for getPriorParametersPooled --------------------------------------------

test_that("valid input returns prior_parameters_list with correct structure", {
  result <- getPriorParametersPooled(
    target_rates = c(0.2, 0.4, 0.6),
    n_worth = 2
  )
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "pooled")
  expect_true(is.list(result$pooled))
  expect_named(result$pooled, c("a", "b"))
  expect_true(result$pooled$a > 0)
  expect_true(result$pooled$b > 0)
})

test_that("error if target_rates contains invalid values", {
  expect_error(getPriorParametersPooled(target_rates = c(0, 0.5), n_worth = 1), "target_rates")
  expect_error(getPriorParametersPooled(target_rates = c(1, 0.5), n_worth = 1), "target_rates")
})

test_that("error if n_worth is non-positive", {
  expect_error(getPriorParametersPooled(target_rates = c(0.3, 0.4), n_worth = 0), "n_worth")
})

test_that("selects target_rate closest to 0.5", {
  result <- getPriorParametersPooled(target_rates = c(0.1, 0.4, 0.8), n_worth = 3)
  expect_equal(result$pooled$a, 0.4 * 3)
  expect_equal(result$pooled$b, (1 - 0.4) * 3)
})

# Tests for setPriorParametersPooled -------------------------------------------


test_that("valid input returns prior_parameters_list with correct structure", {
  result <- setPriorParametersPooled(a = 2, b = 3)
  
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "pooled")
  expect_true(is.list(result$pooled))
  expect_named(result$pooled, c("a", "b"))
  expect_equal(result$pooled$a, 2)
  expect_equal(result$pooled$b, 3)
})

# Tests for getPriorParametersStratified ---------------------------------------

test_that("computes a_j and b_j correctly for multiple target rates", {
  target_rates <- c(0.2, 0.4, 0.6)
  n_worth <- 2
  result <- getPriorParametersStratified(target_rates, n_worth)
  
  expect_equal(result$stratified$a_j, target_rates * n_worth)
  expect_equal(result$stratified$b_j, (1 - target_rates) * n_worth)
})

test_that("handles single target rate correctly", {
  result <- getPriorParametersStratified(target_rates = 0.3, n_worth = 5)
  expect_equal(result$stratified$a_j, 0.3 * 5)
  expect_equal(result$stratified$b_j, (1 - 0.3) * 5)
})

test_that("result structure is consistent", {
  result <- getPriorParametersStratified(c(0.25, 0.75), 3)
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, "stratified")
  expect_named(result$stratified, c("a_j", "b_j"))
})

# Tests for setPriorParametersStratified ---------------------------------------

test_that("returns correct structure and preserves input values", {
  a_j <- c(1, 2, 3)
  b_j <- c(4, 5, 6)
  result <- setPriorParametersStratified(a_j, b_j)
  
  expect_equal(result$stratified$a_j, a_j)
  expect_equal(result$stratified$b_j, b_j)
  expect_s3_class(result, "prior_parameters_list")
})

test_that("works with single-element vectors", {
  result <- setPriorParametersStratified(a_j = 10, b_j = 20)
  expect_equal(result$stratified$a_j, 10)
  expect_equal(result$stratified$b_j, 20)
})

test_that("length consistency check works", {
  expect_error(setPriorParametersStratified(a_j = c(1, 2), b_j = c(3)), "a_j and b_j")
})

# Tests for getPriorParameters -------------------------------------------------

test_that("valid input returns a prior_parameters_list with correct names and class", {
  result <- getPriorParameters(
    method_names = c("berry", "pooled", "stratified", "exnex", "exnex_adj"),
    target_rates = c(0.2, 0.1, 0.3, 0.2, 0.4, 0.6),
    n_worth = 2,
    tau_scale = 1,
    w_j = 0.5
  )
  expect_s3_class(result, "prior_parameters_list")
  expect_named(result, c("berry", "exnex", "exnex_adj", "pooled", "stratified"))
  expect_type(result, "list")
  expect_length(result, 5)

})

test_that("invalid method_names throws error", {
  expect_error(getPriorParameters(
    method_names = "invalid",
    target_rates = c(0.2, 0.3)
  ))
})

test_that("invalid target_rates throws error", {
  expect_error(getPriorParameters(
    method_names = "berry",
    target_rates = c(-0.1, 1.2)
  ))
})

test_that("invalid tau_scale throws error", {
  expect_error(getPriorParameters(
    method_names = "berry",
    target_rates = c(0.2, 0.3),
    tau_scale = -1
  ))
})

test_that("invalid n_worth throws error", {
  expect_error(getPriorParameters(
    method_names = "berry",
    target_rates = c(0.2, 0.3),
    n_worth = 0
  ))
})

test_that("invalid w_j throws error", {
  expect_error(getPriorParameters(
    method_names = "berry",
    target_rates = c(0.2, 0.3),
    w_j = 2
  ))
})

# Tests for combinePriorParameters ---------------------------------------------

test_that("combinePriorParameters returns correct structure with real objects", {
  prior_parameters_stratified <- setPriorParametersStratified(c(1, 2), c(3, 4))
  prior_parameters_berry      <- setPriorParametersBerry(1, 1, 2)
  
  result <- combinePriorParameters(list(prior_parameters_berry, prior_parameters_stratified))
  
  # Check class and type
  expect_s3_class(result, "prior_parameters_list")
  expect_type(result, "list")
  
  # Names should match method names from input
  expect_named(result, sort(c(names(prior_parameters_berry), names(prior_parameters_stratified))))
  
  # Length should equal number of input elements
  expect_length(result, 2)
  
  # Each element should be a list (inner structure)
  expect_true(all(vapply(result, is.list, logical(1))))
})

test_that("combinePriorParameters sorts names alphabetically", {
  prior_parameters_stratified <- setPriorParametersStratified(c(1, 2), c(3, 4))
  prior_parameters_berry      <- setPriorParametersBerry(0, 1, 2)
  
  result <- combinePriorParameters(list(prior_parameters_stratified, prior_parameters_berry))
  
  expect_equal(names(result), sort(c(names(prior_parameters_stratified), names(prior_parameters_berry))))
})

test_that("combinePriorParameters errors on duplicate method names", {
  prior_parameters_berry1 <- setPriorParametersBerry(0, 1, 2)
  prior_parameters_berry2 <- setPriorParametersBerry(0, 1, 2)
  
  expect_error(
    combinePriorParameters(list(prior_parameters_berry1, prior_parameters_berry2))
  )
})

test_that("combinePriorParameters errors if input is not a list", {
  expect_error(combinePriorParameters("not_a_list"))
})

test_that("combinePriorParameters errors if any element is not prior_parameters_list", {
  prior_parameters_berry <- setPriorParametersBerry(0, 1, 2)
  bad_input <- list(prior_parameters_berry, list(dummy = TRUE))
  
  expect_error(combinePriorParameters(bad_input))
})
