applicablePreviousTrials <- function(
    
  scenario_list,
  method_names,
  quantiles,
  n_cohorts,
  calc_differences
  
) {
  
  ## analyze only unique trials that have not been previously analyzed,
  ## i.e. trials that were updated with continueRecruitment(),
  ## i.e. trials that have an overall go decision from a previous decision rule
  ## i.e. trials that have quantiles stored for each cohort that is to be analysed
  applicable_previous_trials <-
    ## check that in each scenario the same analysis methods were analyzed previously
    all(sapply(seq_along(scenario_list), function (i) {
      isTRUE(all.equal(names(scenario_list[[i]]$previous_analyses$post_quantiles),
                       names(scenario_list[[1]]$previous_analyses$post_quantiles)))
    })) &
    ## check that the current analysis method names match the method names of the previous analyses
    isTRUE(all.equal(names(scenario_list[[1]]$previous_analyses$post_quantiles), method_names)) &
    ## check that the stored quantiles are the same across all scenarios
    all(sapply(seq_along(scenario_list), function (i) {
      isTRUE(all.equal(rownames(scenario_list[[i]]$previous_analyses$post_quantiles[[1]][[1]]),
                       rownames(scenario_list[[1]]$previous_analyses$post_quantiles[[1]][[1]])))
    })) &
    ## check that the new quantiles are within the stored quantiles
    all(paste0(as.character(quantiles * 100), "%") %in%
          rownames(scenario_list[[1]]$previous_analyses$post_quantiles[[1]][[1]])) &
    ## check that there are stored quantiles for each cohort that is to be analysed
    all(paste0("p_", seq_len(n_cohorts)) %in%
          colnames(scenario_list[[1]]$previous_analyses$post_quantiles[[1]][[1]]))
  
  ## check that all differences have been previously calculated
  if (!is.null(calc_differences)) {
    
    applicable_previous_trials <- applicable_previous_trials &
      all(apply(calc_differences, 1, function (x) {
        paste0("p_diff_", paste0(as.character(x), collapse = ""))
      }) %in% colnames(scenario_list[[1]]$previous_analyses$post_quantiles[[1]][[1]]))
    
  }
  
  return (applicable_previous_trials)
  
}

calcDiffsMCMCNormal <- function (
    
  posterior_samples,
  calc_differences
  
) {
  
  org_names <- colnames(posterior_samples)
  
  diffs <- apply(calc_differences, 1, function (x) {
    
    matrix(posterior_samples[, grepl(x[1], org_names) & grepl("p", org_names)] -
             posterior_samples[, grepl(x[2], org_names) & grepl("p", org_names)],
           ncol = 1)
    
  })
  
  diff_names <- apply(calc_differences, 1, function (x) {
    
    paste0("p_diff_", paste0(as.character(x), collapse = ""))
    
  })
  
  colnames(diffs) <- diff_names
  
  posterior_samples <- cbind(posterior_samples, diffs)
  
  return (posterior_samples)
  
}

getPosteriorsNormal <- function (
    
  j_parameters,
  j_model_file,
  j_data,
  
  n_mcmc_iterations
  
) {
  
  ## Adaption and burn-in not included
  posterior_samples <- performJags(
    data               = j_data,
    parameters_to_save = j_parameters,
    model_file         = j_model_file, 
    n_iter             = n_mcmc_iterations)
  
  ## replace squarebrackets provided by rjags with workable characters
  colnames(posterior_samples) <- gsub("\\[", "_", colnames(posterior_samples))
  colnames(posterior_samples) <- gsub("\\]", "", colnames(posterior_samples))
  
  return (posterior_samples)
  
}


getPostQuantilesNormal <- function (
    
  ## The method to be applied to the likelihood and the quantiles of the posterior
  method_name,
  quantiles,
  
  ## Scenario data
  scenario_data,
  
  ## Differences between cohorts
  calc_differences  = NULL,
  
  ## JAGS parameters
  j_parameters,
  j_model_file,
  j_data,
  
  ## MCMC Parameters
  n_mcmc_iterations = 1e4,
  
  ## Where to save one of the posterior response rates approximations provided by JAGS
  save_path         = NULL,
  save_trial        = NULL
  
) {

  n_analyses <- nrow(scenario_data$n_subjects)
  
  ## Create random index for saving one of the posterior response rates
  if (is.null(save_trial) && !is.null(save_path)) {
    # set.seed(seed)
    save_trial <- sample(seq_len(n_analyses), size = 1)
  }
  
  ## Run parallel loops
  ## prepare foreach loop over 
  
  exported_stuff <- c(
    "posteriors2QuantilesNormal", "getPosteriorsNormal", "getPostQuantilesOfTrialNormal",
    "qbetaDiffNormal", "chunkVector")
  
  ## prepare chunking
  chunks_outer <- chunkVector(seq_len(n_analyses), foreach::getDoParWorkers())
  
  "%dorng%" <- doRNG::"%dorng%"
  "%dopar%" <- foreach::"%dopar%"
  posterior_quantiles_list <- suppressMessages(
    foreach::foreach(
      k = chunks_outer,
      .combine  = c,
      .verbose  = FALSE,
      .packages = c("rjags"),
      .export   = exported_stuff) %dorng% {
        
        chunks_inner <- chunkVector(k, foreach::getDoParWorkers())
        
        foreach::foreach(i = chunks_inner, .combine = c) %dorng% {
          
          lapply(i, function (j) {
            
            ## Calculate the posterior quantiles for the kth unique trial outcome
            getPostQuantilesOfTrialNormal(
              y_list            = scenario_data$trials[[j]]$y_list,
              j_data            = j_data,
              j_parameters      = j_parameters,
              j_model_file      = j_model_file,
              method_name       = method_name,
              quantiles         = quantiles,
              calc_differences  = calc_differences,
              n_mcmc_iterations = n_mcmc_iterations,
              save_path         = save_path,
              save_trial        = save_trial)
            
          })
          
        }
        
      })
  
  return (posterior_quantiles_list)
  
}

getPostQuantilesOfTrialNormal <- function(
    y_list,
    j_data,
    j_parameters,
    j_model_file,
    method_name,
    quantiles,
    calc_differences,
    n_mcmc_iterations,
    save_path,
    save_trial
) {
  # number of cohorts
  K <- length(y_list)
  if (K < 1L) stop("y_list must be a non-empty list of numeric vectors.")
  
  # subjects per cohort
  n_vec <- vapply(y_list, length, integer(1L))
  if (any(n_vec <= 0L)) {
    stop("All cohorts in y_list must have at least one observation.")
  }
  
  # pad to rectangular y[K, max_n]
  max_n <- max(n_vec)
  y_mat <- matrix(NA_real_, nrow = K, ncol = max_n)
  for (k in seq_len(K)) {
    y_mat[k, seq_len(n_vec[k])] <- y_list[[k]]
  }
  
  # 🔧 IMPORTANT: JAGS expects J, not K
  j_data_full <- c(
    j_data,
    list(
      J = K,      # <- number of cohorts for bhm_normal.txt
      n = n_vec,  # length-J integer vector
      y = y_mat   # J x max_n matrix
    )
  )
  
  # run JAGS
  posterior_samples <- getPosteriorsNormal(
    j_parameters      = j_parameters,
    j_model_file      = j_model_file,
    j_data            = j_data_full,
    n_mcmc_iterations = n_mcmc_iterations
  )
  
  # extract mu_k columns (cohort means)
  mu_cols <- grepl("^mu_[0-9]+$", colnames(posterior_samples))
  mu_post <- posterior_samples[, mu_cols, drop = FALSE]
  
  if (ncol(mu_post) != K) {
    stop("Could not identify exactly K mu columns in posterior_samples.")
  }
  
  # for compatibility with decision code, call them p_1, p_2, ...
  colnames(mu_post) <- paste0("p_", seq_len(K))
  
  # optional differences (mu_i - mu_j)
  if (!is.null(calc_differences)) {
    calc_differences <- convertVector2Matrix(calc_differences)
    mu_post <- calcDiffsMCMCNormal(
      posterior_samples = mu_post,
      calc_differences  = calc_differences
    )
  }
  
  # optional save of raw MCMC sample
  if (!is.null(save_path) && !is.null(save_trial)) {
    saveRDS(
      posterior_samples,
      file = file.path(
        save_path,
        paste0("posterior_samples_", save_trial, "_", method_name, "_rds")
      )
    )
  }
  
  # turn posterior sample into quantile table + mean + sd
  posterior2quantilesNormal(
    quantiles  = quantiles,
    posteriors = mu_post
  )
}

#' @title getUniqueTrialsNormal
#' @description Create unique trial groups for normal endpoints by binning
#'   trial-level sample means per cohort.
#' @param scenario_list Object of class `scenario_list_normal`.
#' @param nbins Positive integer: number of bins per cohort (if `bin_breaks` is NULL).
#' @param bin_breaks Optional list of numeric vectors of breaks per cohort
#'   (names must match the cohort names in `y_list`), used instead of `nbins`.
#' @return A list with
#'   \item{groups}{data.frame with one row per group:
#'         `group_id`, `signature`, `n_members`.}
#'   \item{map}{data.frame mapping each (scenario, trial) to its `group_id`.}
#'   \item{breaks}{list of breaks per cohort actually used.}
getUniqueTrialsNormal <- function(
    scenario_list,
    nbins      = 5,
    bin_breaks = NULL
) {
  if (missing(scenario_list)) {
    stop("Please provide 'scenario_list' for getUniqueTrialsNormal().")
  }
  if (!is.scenario_list_normal(scenario_list)) {
    stop("scenario_list must be of class 'scenario_list_normal'.")
  }
  if (!is.null(bin_breaks) && !is.list(bin_breaks)) {
    stop("'bin_breaks' must be a list or NULL.")
  }
  if (is.null(bin_breaks) && (!is.numeric(nbins) || length(nbins) != 1L || nbins < 1)) {
    stop("'nbins' must be a positive integer if 'bin_breaks' is NULL.")
  }
  
  ybar_list <- list()
  scen_idx  <- integer(0)
  trial_idx <- integer(0)
  
  # 1) Collect per-cohort summary means (already computed in simulateScenariosNormal)
  for (s in seq_along(scenario_list)) {
    trials <- scenario_list[[s]]$trials
    if (is.null(trials) || !length(trials)) {
      stop("scenario_list[[", s, "]] has no trials.")
    }
    
    for (t in seq_along(trials)) {
      y_bar <- trials[[t]]$y_bar
      
      if (is.null(y_bar) || !is.numeric(y_bar) || length(y_bar) == 0L) {
        stop(
          "scenario_list[[", s, "]]$trials[[", t, "]]$y_bar must be a non-empty numeric vector.\n",
          "Tip: compute y_bar during simulation and store it per trial."
        )
      }
      
      # ensure names exist
      if (is.null(names(y_bar))) {
        names(y_bar) <- paste0("c_", seq_along(y_bar))
      }
      
      ybar_list[[length(ybar_list) + 1L]] <- y_bar
      scen_idx  <- c(scen_idx, s)
      trial_idx <- c(trial_idx, t)
    }
  }
  
  # 2) Align cohorts across all trials (safe even if names differ)
  all_cohorts <- unique(unlist(lapply(ybar_list, names), use.names = FALSE))
  if (is.null(all_cohorts) || !length(all_cohorts)) {
    stop("Could not determine cohort names from y_bar.")
  }
  
  ybar_mat <- t(vapply(ybar_list, function(v) {
    out <- rep(NA_real_, length(all_cohorts))
    names(out) <- all_cohorts
    out[names(v)] <- v
    out
  }, numeric(length(all_cohorts))))
  colnames(ybar_mat) <- all_cohorts
  
  index_df <- data.frame(
    scenario_index   = scen_idx,
    trial_index      = trial_idx,
    stringsAsFactors = FALSE
  )
  
  # 3) Build breaks per cohort
  if (is.null(bin_breaks)) {
    breaks <- lapply(seq_len(ncol(ybar_mat)), function(j) {
      r <- range(ybar_mat[, j], finite = TRUE)
      if (!all(is.finite(r))) {
        stop("Non-finite y_bar values found; cannot create bins.")
      }
      if (r[1L] == r[2L]) r <- r + c(-0.5, 0.5)
      seq(from = r[1L], to = r[2L], length.out = nbins + 1L)
    })
    names(breaks) <- all_cohorts
  } else {
    if (!all(all_cohorts %in% names(bin_breaks))) {
      stop(
        "All cohorts must appear as names in 'bin_breaks'. Missing: ",
        paste(setdiff(all_cohorts, names(bin_breaks)), collapse = ", ")
      )
    }
    breaks <- bin_breaks[all_cohorts]
  }
  
  # 4) Bin each cohort mean and create a signature per trial
  bin_codes <- lapply(seq_len(ncol(ybar_mat)), function(j) {
    cut(
      ybar_mat[, j],
      breaks         = breaks[[j]],
      include.lowest = TRUE,
      labels         = FALSE
    )
  })
  bin_codes <- as.data.frame(bin_codes, check.names = FALSE)
  colnames(bin_codes) <- all_cohorts
  
  signature <- apply(bin_codes, 1L, function(z) paste(z, collapse = "-"))
  
  sig_levels <- unique(signature)
  group_id   <- match(signature, sig_levels)
  
  groups_df <- data.frame(
    group_id  = seq_along(sig_levels),
    signature = sig_levels,
    n_members = as.integer(tabulate(group_id, nbins = length(sig_levels))),
    stringsAsFactors = FALSE
  )
  
  map_df <- data.frame(
    scenario_index   = index_df$scenario_index,
    trial_index      = index_df$trial_index,
    group_id         = group_id,
    signature        = signature,
    stringsAsFactors = FALSE
  )
  
  list(
    groups = groups_df,
    map    = map_df,
    breaks = breaks
  )
}

#' @title mapUniqueTrialsNormal
#' @description Map group-level posterior summaries back to all trials
#'   for normal endpoints.
#' @param scenario_list Object of class `scenario_list_normal`.
#' @param method_quantiles_by_group List of length = number of groups,
#'   each element is the posterior summary for that group's pooled data.
#' @param unique_info Output from getUniqueTrialsNormal() (`groups`, `map`).
#' @return A list of length(scenario_list); each element is a list of length
#'   `n_trials` with the corresponding posterior quantiles.
mapUniqueTrialsNormal <- function(
    scenario_list,
    method_quantiles_by_group,
    unique_info
) {
  if (missing(scenario_list)) {
    stop("Please provide 'scenario_list' for mapUniqueTrialsNormal().")
  }
  if (!is.scenario_list_normal(scenario_list)) {
    stop("scenario_list must be of class 'scenario_list_normal'.")
  }
  if (missing(method_quantiles_by_group)) {
    stop("Please provide 'method_quantiles_by_group' (one result per group).")
  }
  if (missing(unique_info) || !is.list(unique_info) ||
      is.null(unique_info$map) || is.null(unique_info$groups)) {
    stop("Please provide 'unique_info' as returned by getUniqueTrialsNormal().")
  }
  
  mapping  <- unique_info$map
  groups   <- unique_info$groups
  n_groups <- nrow(groups)
  
  if (length(method_quantiles_by_group) != n_groups) {
    stop("Length of 'method_quantiles_by_group' must equal nrow(unique_info$groups).")
  }
  
  scenario_numbers <- sapply(scenario_list, function(x) x$scenario_number)
  out <- vector("list", length(scenario_list))
  names(out) <- paste0("scenario_", scenario_numbers)
  
  for (s in seq_along(scenario_list)) {
    n_trials_s <- length(scenario_list[[s]]$trials)
    qlist_s    <- vector("list", n_trials_s)
    
    map_s <- mapping[mapping$scenario_index == s, , drop = FALSE]
    
    for (i in seq_len(nrow(map_s))) {
      t_idx <- map_s$trial_index[i]
      g     <- map_s$group_id[i]
      qlist_s[[t_idx]] <- method_quantiles_by_group[[g]]
    }
    
    out[[s]] <- qlist_s
  }
  
  out
}


loadAnalysesNormal <- function(
    scenario_numbers,
    analysis_numbers = rep(1, length(scenario_numbers)),
    load_path        = tempdir()
) {
  loadAnalyses(
    scenario_numbers = scenario_numbers,
    analysis_numbers = analysis_numbers,
    load_path        = load_path
  )
}


performAnalysesNormal <- function(
    scenario_list,
    evidence_levels         = c(0.025, 0.05, 0.5, 0.8, 0.9, 0.95, 0.975),
    calc_differences        = NULL,
    n_mcmc_iterations       = 1e4,
    prior_parameters_normal = NULL,
    nbins                   = 5,
    bin_breaks              = NULL,
    verbose                 = TRUE
) {
  error_scenario_list <- simpleError(
    "Please provide an object of class 'scenario_list_normal' for 'scenario_list'")
  error_evidence_levels <- simpleError(
    "Please provide a vector of numerics in (0, 1) for the argument 'evidence_levels'")
  error_n_mcmc_iterations <- simpleError(
    "Please provide a positive integer for the argument 'n_mcmc_iterations'")
  error_verbose <- simpleError(
    "Please provide a logical for the argument 'verbose'")
  
  if (missing(scenario_list)) stop(error_scenario_list)
  if (!is.scenario_list_normal(scenario_list)) stop(error_scenario_list)
  
  if (!is.numeric.in.zero.one(evidence_levels)) stop(error_evidence_levels)
  if (!is.single.positive.wholenumber(n_mcmc_iterations)) stop(error_n_mcmc_iterations)
  if (!is.logical(verbose) || length(verbose) != 1L) stop(error_verbose)
  
  if (!is.null(calc_differences)) {
    calc_differences <- convertVector2Matrix(calc_differences)
  }
  
  quantiles <- sort(unique(round(
    1 - c(0.025, 0.05, 0.5, 0.8, 0.9, 0.95, 0.975, evidence_levels),
    9
  )))
  
  if (is.null(prior_parameters_normal)) {
    prior_parameters_normal <- getPriorParametersNormal()
  }
  
  prep         <- prepareAnalysisNormal(prior_parameters = prior_parameters_normal)
  j_parameters <- prep$j_parameters
  j_model_file <- prep$j_model_file
  j_data_fixed <- prep$j_data
  
  scenario_numbers <- sapply(scenario_list, function(x) x$scenario_number)
  
  # --- 1) Build unique groups based on per-cohort means & bins ---
  unique_info <- getUniqueTrialsNormal(
    scenario_list = scenario_list,
    nbins         = nbins,
    bin_breaks    = bin_breaks
  )
  
  groups    <- unique_info$groups
  mapping   <- unique_info$map
  n_groups  <- nrow(groups)
  
  # Get cohort structure from first trial
  y_example    <- scenario_list[[1]]$trials[[1]]$y_list
  n_cohorts    <- length(y_example)
  cohort_names <- names(y_example)
  if (is.null(cohort_names)) {
    cohort_names <- paste0("c_", seq_len(n_cohorts))
  }
  
  if (verbose) {
    message(
      format(Sys.time(), "%d-%b-%Y"),
      " Performing Analyses (normal endpoint)"
    )
    message(
      "   Analyzing ",
      length(scenario_numbers), " scenario",
      ifelse(length(scenario_numbers) == 1, "", "s"),
      " (", n_groups, " unique trial group",
      ifelse(n_groups == 1, "", "s"), ")"
    )
  }
  
  # --- 2) Build pooled y_list per group and run JAGS once per group ---
  if (verbose) {
    start_time <- Sys.time()
    message("   Running bhm_normal model on pooled groups ...")
  }
  
  method_quantiles_by_group <- vector("list", n_groups)
  
  for (g in seq_len(n_groups)) {
    members <- mapping$group_id == g
    scen_g  <- mapping$scenario_index[members]
    trial_g <- mapping$trial_index[members]
    
    # pool values within each cohort over all member trials
    y_group <- vector("list", n_cohorts)
    for (k in seq_len(n_cohorts)) {
      vals <- numeric(0)
      for (idx in seq_along(scen_g)) {
        y_list_k <- scenario_list[[scen_g[idx]]]$trials[[trial_g[idx]]]$y_list[[k]]
        vals     <- c(vals, y_list_k)
      }
      y_group[[k]] <- vals
    }
    names(y_group) <- cohort_names
    
    # run JAGS once for this group's pooled data
    method_quantiles_by_group[[g]] <- getPostQuantilesOfTrialNormal(
      y_list            = y_group,
      j_data            = j_data_fixed,
      j_parameters      = j_parameters,
      j_model_file      = j_model_file,
      method_name       = "normal",
      quantiles         = quantiles,
      calc_differences  = calc_differences,
      n_mcmc_iterations = n_mcmc_iterations,
      save_path         = NULL,
      save_trial        = FALSE
    )
  }
  
  if (verbose) {
    message(
      "       finished after ",
      round(Sys.time() - start_time, 1), " ",
      units(Sys.time() - start_time), "."
    )
  }
  
  # --- 3) Map group-level posteriors back to all trials ---
  if (verbose) {
    start_time <- Sys.time()
    message("   Processing scenarios (mapping group results) ...")
  }
  
  scenario_method_quantiles_list <- mapUniqueTrialsNormal(
    scenario_list            = scenario_list,
    method_quantiles_by_group = method_quantiles_by_group,
    unique_info              = unique_info
  )
  
  if (verbose) {
    message(
      "       finished after ",
      round(Sys.time() - start_time, 1), " ",
      units(Sys.time() - start_time), "."
    )
  }
  
  # --- 4) Wrap up in analysis_list_normal (unchanged structure) ---
  analyses_list        <- vector("list", length(scenario_numbers))
  names(analyses_list) <- paste0("scenario_", scenario_numbers)
  
  for (s in seq_along(scenario_numbers)) {
    analyses_list[[s]] <- list(
      quantiles_list      = list(normal = scenario_method_quantiles_list[[s]]),
      scenario_data       = scenario_list[[s]],
      analysis_parameters = list(
        quantiles               = quantiles,
        evidence_levels         = evidence_levels,
        method_names            = "normal",
        prior_parameters_normal = prior_parameters_normal,
        n_mcmc_iterations       = n_mcmc_iterations,
        nbins                   = nbins,
        bin_breaks              = unique_info$breaks
      )
    )
  }
  
  class(analyses_list) <- "analysis_list_normal"
  analyses_list
}

## based on R2jags::jags
## stripped down to improve performance
performJags <- function (
    
  data,
  parameters_to_save,
  model_file, 
  n_chains = 2,
  n_iter   = 1e4,
  n_burnin = floor(n_iter/3)
  
) {
  
  n_adapt <- ifelse(n_burnin > 0, n_burnin, 100)
  
  inits <- vector("list", n_chains)
  for (i in 1:n_chains) {
    inits[[i]]$.RNG.name <- "base::Wichmann-Hill"
    inits[[i]]$.RNG.seed <- stats::runif(1, 0, 2^31)
  }
  
  j_model <- rjags::jags.model(file     = model_file,
                               data     = data,
                               inits    = inits, 
                               n.chains = n_chains,
                               n.adapt  = 0,
                               quiet    = TRUE)
  
  rjags::adapt(object         = j_model,
               n.iter         = n_adapt,
               progress.bar   = "none",
               end.adaptation = TRUE)
  
  samples <- rjags::coda.samples(model          = j_model,
                                 variable.names = parameters_to_save, 
                                 n.iter         = n_iter - n_burnin,
                                 thin           = 1, 
                                 progress.bar   = "none")
  
  return(do.call(rbind, samples))
  
}

posterior2quantilesNormal <- function (
    
  quantiles,
  posteriors
  
) {
  
  posterior_quantiles <- apply(posteriors, 2, function (x) stats::quantile(x, probs = quantiles))
  
  posterior_mean      <- apply(posteriors, 2, mean)
  posterior_sd        <- apply(posteriors, 2, stats::sd)
  posterior_quantiles <- rbind(posterior_quantiles,
                               Mean = posterior_mean,
                               SD   = posterior_sd)
  
  return (posterior_quantiles)
  
}
prepareAnalysisNormal <- function(prior_parameters = NULL) {
  
  if (is.null(prior_parameters)) {
    prior_parameters <- getPriorParametersNormal()
  } else {
    if (!inherits(prior_parameters, "prior_parameters_normal")) {
      stop("prior_parameters must have class 'prior_parameters_normal'")
    }
  }
  
  j_data <- list(
    mu_pop_mean = prior_parameters$mu_pop_mean,
    prec_mu_pop = 1 / prior_parameters$mu_pop_sd^2,
    tau_shape   = prior_parameters$tau_shape,
    tau_rate    = prior_parameters$tau_rate,
    sigma_shape = prior_parameters$sigma_shape,
    sigma_rate  = prior_parameters$sigma_rate
  )
  
  j_model_file <- system.file(
    package   = "bhmbasket",
    "jags_models",
    "normal.txt",
    mustWork  = TRUE
  )
  
  j_parameters <- c("mu", "mu_pop", "prec_tau", "prec_sigma")
  
  list(
    j_parameters = j_parameters,
    j_model_file = j_model_file,
    j_data       = j_data
  )
}

print.analysis_listnormal <- function(x, digits = 2, ...) {
  print.analysis_list(x, digits = digits, ...)
}

qbetaDiffNormal <- function(
    quantiles,
    x_1_shape1,
    x_1_shape2,
    x_2_shape1,
    x_2_shape2,
    n_mcmc = 1e6
) {
  qbetaDiff(
    quantiles  = quantiles,
    x_1_shape1 = x_1_shape1,
    x_1_shape2 = x_1_shape2,
    x_2_shape1 = x_2_shape1,
    x_2_shape2 = x_2_shape2,
    n_mcmc     = n_mcmc
  )
}

saveAnalysesNormal <- function(
    analyses_list,
    save_path        = tempdir(),
    analysis_numbers = NULL
) {
  saveAnalyses(
    analyses_list    = analyses_list,
    save_path        = save_path,
    analysis_numbers = analysis_numbers
  )
}

