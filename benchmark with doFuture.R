library(bhmbasket)
library(future)
library(doFuture)

# --- Create 100 scenarios ---
set.seed(999)
scenario_list <- list()

for (i in 1:500) {
  trial <- createTrial(
    n_subjects   = matrix(c(30,30,30), nrow = 1),
    n_responders = matrix(rbinom(3, size = 30, prob = c(0.2, 0.25, 0.3)), nrow = 1)
  )
  
  # createTrial returns a scenario_list of length 1
  scenario_list[[i]] <- trial[[1]]
}

class(scenario_list) <- "scenario_list"

# -------------------------------------------------------------
# Benchmark 1: sequential (single worker)
# -------------------------------------------------------------
plan(sequential)
t_seq <- system.time(
  performAnalyses(
    scenario_list       = scenario_list,
    target_rates        = rep(0.3, 3),
    method_names        = "berry",
    n_mcmc_iterations   = 10000,
    verbose             = FALSE
  )
)

# -------------------------------------------------------------
# Benchmark 2: parallel (multisession)
# -------------------------------------------------------------
plan(multisession)  # automatically uses all cores
t_par <- system.time(
  performAnalyses(
    scenario_list       = scenario_list,
    target_rates        = rep(0.3, 3),
    method_names        = "berry",
    n_mcmc_iterations   = 10000,
    verbose             = FALSE
  )
)

rbind(sequential = t_seq, multisession = t_par)


# ------------------------------------------------------------------------------

plan(sequential)     # run sequential
set.seed(111)

res_seq <- performAnalyses(
  scenario_list     = scenario_list,
  target_rates      = rep(0.3, 3),
  method_names      = "berry",
  n_mcmc_iterations = 5000,
  verbose           = FALSE
)

plan(multisession)   # run parallel
set.seed(111)

res_par <- performAnalyses(
  scenario_list     = scenario_list,
  target_rates      = rep(0.3, 3),
  method_names      = "berry",
  n_mcmc_iterations = 5000,
  verbose           = FALSE
)

# Direct comparison
identical(res_seq, res_par)

# ------------------------------------------------------------------------------

library(bhmbasket)
library(future)
library(doFuture)
library(ggplot2)

set.seed(999)
scenario_list <- list()

for (i in 1:500) {
  trial <- createTrial(
    n_subjects   = matrix(c(30, 30, 30), nrow = 1),
    n_responders = matrix(rbinom(3, size = 30, prob = c(0.2, 0.25, 0.3)), nrow = 1)
  )
  scenario_list[[i]] <- trial[[1]]
}

class(scenario_list) <- "scenario_list"

run_benchmark <- function(workers,
                          scheduling = NULL,
                          chunk_size = NULL,
                          n_mcmc_iterations = 5000) {
  
  plan(multisession, workers = workers)
  
  opts <- list(seed = TRUE)
  
  if (!is.null(chunk_size)) {
    opts$chunk.size <- chunk_size
  } else if (!is.null(scheduling)) {
    opts$scheduling <- scheduling
  }
  
  t <- system.time({
    res <- foreach(
      i = seq_along(scenario_list),
      .combine = c,
      .options.future = opts,
      .packages = "bhmbasket"
    ) %dofuture% {
      performAnalyses(
        scenario_list     = structure(list(scenario_list[[i]]), class = "scenario_list"),
        target_rates      = rep(0.3, 3),
        method_names      = "berry",
        n_mcmc_iterations = n_mcmc_iterations,
        verbose           = FALSE
      )
    }
  })
  
  plan(sequential)
  
  data.frame(
    workers = workers,
    scheduling = if (is.null(scheduling)) NA else scheduling,
    chunk_size = if (is.null(chunk_size)) NA else chunk_size,
    elapsed = unname(t["elapsed"]),
    user.self = unname(t["user.self"]),
    sys.self = unname(t["sys.self"])
  )
}

worker_grid <- c(2, 4, 8, 16, 32)
scheduling_grid <- c(1, 2, 4)
chunk_grid <- c(5, 10, 20)

results_sched <- do.call(
  rbind,
  lapply(worker_grid, function(w) {
    do.call(
      rbind,
      lapply(scheduling_grid, function(s) {
        run_benchmark(workers = w, scheduling = s, chunk_size = NULL)
      })
    )
  })
)

results_chunk <- do.call(
  rbind,
  lapply(worker_grid, function(w) {
    do.call(
      rbind,
      lapply(chunk_grid, function(cs) {
        run_benchmark(workers = w, scheduling = NULL, chunk_size = cs)
      })
    )
  })
)

results_all <- rbind(results_sched, results_chunk)
print(results_all)

best_run <- results_all[which.min(results_all$elapsed), ]
print(best_run)

p1 <- ggplot(results_sched, aes(x = workers, y = elapsed, color = factor(scheduling), group = factor(scheduling))) +
  geom_line() +
  geom_point(size = 2) +
  labs(color = "scheduling")

p2 <- ggplot(results_chunk, aes(x = workers, y = elapsed, color = factor(chunk_size), group = factor(chunk_size))) +
  geom_line() +
  geom_point(size = 2) +
  labs(color = "chunk_size")

print(p1)
print(p2)

#------------


library(future)
library(parallelly)

availableCores()

plan(multisession, workers = 4L)
nbrOfWorkers()

plan(multisession, workers = 8L)
nbrOfWorkers()

plan(sequential)

