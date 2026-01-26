library("testthat")# HTML report for that subset
library("covr")
library("foreach")

# From package root (bhmbasket/)
src <- "R/AnalysisFunctions.R" # just this R file
tst <- "tests/testthat/test-AnalysisFunctions.R"             # just this test file
cov1 <- covr::file_coverage(source_files = src, test_files = tst)
covr::report(cov1)

# cov <- covr::package_coverage(path = "/fsx/home/chattota/bhmbasket", type = "tests")
getwd()
devtools::load_all()

cov <- covr::code_coverage(code = testthat::test_local(stop_on_failure = FALSE, stop_on_warning= FALSE, reporter = "summary"), type = "tests")
packageVersion("covr")
?covr

testthat::test_local()

source("R/DataFunctions.R")
source("R/AnalysisFunctions.R")
source("R/CheckingFunctions.R")
source("R/OCFunctions.R")
source("R/PriorFunctions.R")
source("R/Misc.R")

library(doFuture)
library(future)

doFuture::registerDoFuture()
future::plan(multisession)

source("https://raw.githubusercontent.com/MangoTheCat/remotes/master/install-github.R")$value("mangothecat/functionMap")
library(functionMap)

# Create a function map from an R script
map <- functionMap::map_cran_package("bhmbasket")

# github actions - we can run tests automatically, push a function one by one and then check that the output is the same



















## 6. Integration: AND/OR example from documentation ---------------------------

test_that("negateGoDecisions: integrates correctly with AND/OR design from docs", {
  set.seed(202)

  scenarios_list <- simulateScenarios(
    n_subjects_list     = list(c(10, 20, 30)),
    response_rates_list = list(rep(0.9, 3)),
    n_trials            = 10
  )

  analysis_list <- performAnalyses(
    scenario_list     = scenarios_list,
    target_rates      = rep(0.5, 3),
    n_mcmc_iterations = 100
  )

  cohort_names <- c("p_1", "p_2", "p_3",
                    "p_1", "p_2", "p_3")

  evidence_levels <- c(0.5,  0.5,  0.5,
                       0.95, 0.95, 0.95)

  # Go: Significance AND Relevance for each cohort
  go_decisions_list <- getGoDecisions(
    analyses_list   = analysis_list,
    cohort_names    = cohort_names,
    evidence_levels = evidence_levels,
    boundary_rules  = quote(c(
      x[1] > 0.8 & x[4] > 0.4,
      x[2] > 0.8 & x[5] > 0.4,
      x[3] > 0.8 & x[6] > 0.4
    ))
  )

  # NoGo: negation of (Significance OR Relevance) for each cohort
  nogo_decisions <- negateGoDecisions(
    getGoDecisions(
      analyses_list   = analysis_list,
      cohort_names    = cohort_names,
      evidence_levels = evidence_levels,
      boundary_rules  = quote(c(
        x[1] > 0.8 | x[4] > 0.4,
        x[2] > 0.8 | x[5] > 0.4,
        x[3] > 0.8 | x[6] > 0.4
      ))
    )
  )

  expect_s3_class(nogo_decisions, "decision_list")
  expect_identical(names(nogo_decisions), names(go_decisions_list))

  for (s in seq_along(go_decisions_list)) {
    meth_names <- names(go_decisions_list[[s]]$decisions_list)

    for (m in meth_names) {
      go_mat   <- go_decisions_list[[s]]$decisions_list[[m]]
      nogo_mat <- nogo_decisions[[s]]$decisions_list[[m]]

      # Same shape is required for later use in getGoProbabilities()
      expect_identical(dim(go_mat), dim(nogo_mat))

      # No cell may be both Go and NoGo at the same time.
      overlap <- go_mat & nogo_mat
      expect_false(
        any(overlap),
        info = paste("Overlap of Go and NoGo decisions in scenario", s, "method", m)
      )
    }
  }
})

