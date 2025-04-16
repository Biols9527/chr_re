#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Model Selection Module
# Author: Bioinformatics Team
# Date: 2023-08-10
# Description: Methods for fitting, comparing, and selecting evolutionary models
#              for chromosome number evolution, including continuous and discrete
#              models, with criteria for model selection and visualization tools
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(phytools)
  library(ggplot2)
  library(parallel)
  library(viridis)
  library(dplyr)
})

#===============================================================================
# Model Fitting Functions
#===============================================================================

#' Fit a suite of continuous trait evolution models to chromosome count data
#' 
#' Fits multiple models of chromosome evolution including Brownian Motion,
#' Ornstein-Uhlenbeck, Early Burst, and others
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param models Vector of model names to fit (default fits all available models)
#' @param criterion Model selection criterion: "AIC", "BIC", "AICc"
#' @param discrete Whether to treat chromosome counts as discrete or continuous
#' @param n_cores Number of cores for parallel model fitting (NULL = auto-detect)
#' @return List with model fitting results and comparisons
#' @export
fit_chromosome_models <- function(tree, 
                               chr_counts, 
                               models = c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta"),
                               criterion = "AICc",
                               discrete = TRUE,
                               n_cores = NULL) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Match taxa between tree and data
  shared_taxa <- intersect(tree$tip.label, names(chr_counts))
  
  if(length(shared_taxa) < 4) {
    stop("At least 4 shared taxa required for model fitting")
  }
  
  # Prune tree and data to shared taxa
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, shared_taxa))
  pruned_counts <- chr_counts[pruned_tree$tip.label]
  
  # Check if geiger is available
  if(!requireNamespace("geiger", quietly = TRUE)) {
    stop("The geiger package is required for model fitting")
  }
  
  # Set up parallel computing if requested
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  if(n_cores > 1) {
    # Set up cluster
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("pruned_tree", "pruned_counts", "discrete"), 
                       envir = environment())
    
    # Load required packages on all nodes
    parallel::clusterEvalQ(cl, {
      library(geiger)
    })
    
    # Fit models in parallel
    model_fits <- parallel::parLapply(cl, models, function(model) {
      tryCatch({
        if(discrete) {
          # Round data to integers for discrete models
          discrete_counts <- round(pruned_counts)
          fit <- geiger::fitContinuous(pruned_tree, discrete_counts, model = model)
        } else {
          fit <- geiger::fitContinuous(pruned_tree, pruned_counts, model = model)
        }
        return(fit)
      }, error = function(e) {
        message(paste("Error fitting", model, "model:", e$message))
        return(NULL)
      })
    })
    names(model_fits) <- models
  } else {
    # Fit models sequentially
    model_fits <- list()
    for(model in models) {
      model_fits[[model]] <- tryCatch({
        if(discrete) {
          # Round data to integers for discrete models
          discrete_counts <- round(pruned_counts)
          geiger::fitContinuous(pruned_tree, discrete_counts, model = model)
        } else {
          geiger::fitContinuous(pruned_tree, pruned_counts, model = model)
        }
      }, error = function(e) {
        message(paste("Error fitting", model, "model:", e$message))
        return(NULL)
      })
    }
  }
  
  # Remove failed fits
  model_fits <- model_fits[!sapply(model_fits, is.null)]
  
  if(length(model_fits) == 0) {
    stop("All models failed to fit")
  }
  
  # Extract model comparison metrics
  model_comparison <- data.frame(
    Model = names(model_fits),
    LogLik = sapply(model_fits, function(x) x$opt$lnL),
    AIC = sapply(model_fits, function(x) x$opt$aic),
    AICc = sapply(model_fits, function(x) x$opt$aicc),
    Parameters = sapply(model_fits, function(x) length(x$opt$npars)),
    stringsAsFactors = FALSE
  )
  
  # Calculate BIC
  model_comparison$BIC <- model_comparison$AIC + 
                         (log(length(pruned_counts)) - 2) * model_comparison$Parameters
  
  # Calculate delta values and weights
  if(criterion == "AIC") {
    model_comparison$Delta <- model_comparison$AIC - min(model_comparison$AIC)
    model_comparison$Weight <- exp(-0.5 * model_comparison$Delta) / sum(exp(-0.5 * model_comparison$Delta))
    model_comparison <- model_comparison[order(model_comparison$AIC), ]
  } else if(criterion == "AICc") {
    model_comparison$Delta <- model_comparison$AICc - min(model_comparison$AICc)
    model_comparison$Weight <- exp(-0.5 * model_comparison$Delta) / sum(exp(-0.5 * model_comparison$Delta))
    model_comparison <- model_comparison[order(model_comparison$AICc), ]
  } else if(criterion == "BIC") {
    model_comparison$Delta <- model_comparison$BIC - min(model_comparison$BIC)
    model_comparison$Weight <- exp(-0.5 * model_comparison$Delta) / sum(exp(-0.5 * model_comparison$Delta))
    model_comparison <- model_comparison[order(model_comparison$BIC), ]
  } else {
    stop("criterion must be one of: 'AIC', 'AICc', 'BIC'")
  }
  
  # Record best model
  best_model <- model_comparison$Model[1]
  
  # Generate model parameter summary
  param_summary <- extract_model_parameters(model_fits)
  
  # Create results object
  results <- list(
    model_fits = model_fits,
    model_comparison = model_comparison,
    best_model = best_model,
    param_summary = param_summary,
    tree = pruned_tree,
    chr_counts = pruned_counts,
    discrete = discrete,
    criterion = criterion
  )
  
  class(results) <- c("chr_model_fits", class(results))
  
  return(results)
}

#' Extract parameter estimates from fitted models
#' 
#' @param model_fits List of fitted models from fitContinuous
#' @return Data frame with parameter estimates for each model
#' @keywords internal
extract_model_parameters <- function(model_fits) {
  # Initialize parameter list
  all_params <- list()
  
  # Extract parameters from each model
  for(model_name in names(model_fits)) {
    fit <- model_fits[[model_name]]
    
    # All models have sigsq (Brownian rate)
    params <- data.frame(
      Model = model_name,
      Parameter = "sigsq",
      Value = fit$opt$sigsq,
      stringsAsFactors = FALSE
    )
    
    # Add additional parameters specific to each model
    if(model_name == "OU") {
      # Ornstein-Uhlenbeck model has alpha parameter
      params <- rbind(params, data.frame(
        Model = model_name,
        Parameter = "alpha",
        Value = fit$opt$alpha,
        stringsAsFactors = FALSE
      ))
    } else if(model_name == "EB") {
      # Early Burst model has a parameter
      params <- rbind(params, data.frame(
        Model = model_name,
        Parameter = "a",
        Value = fit$opt$a,
        stringsAsFactors = FALSE
      ))
    } else if(model_name == "trend") {
      # Trend model has mu parameter
      params <- rbind(params, data.frame(
        Model = model_name,
        Parameter = "mu",
        Value = fit$opt$mu,
        stringsAsFactors = FALSE
      ))
    } else if(model_name == "lambda") {
      # Lambda model (phylogenetic signal)
      params <- rbind(params, data.frame(
        Model = model_name,
        Parameter = "lambda",
        Value = fit$opt$lambda,
        stringsAsFactors = FALSE
      ))
    } else if(model_name == "kappa") {
      # Kappa model (punctuational vs. gradual evolution)
      params <- rbind(params, data.frame(
        Model = model_name,
        Parameter = "kappa",
        Value = fit$opt$kappa,
        stringsAsFactors = FALSE
      ))
    } else if(model_name == "delta") {
      # Delta model (time-dependent evolution)
      params <- rbind(params, data.frame(
        Model = model_name,
        Parameter = "delta",
        Value = fit$opt$delta,
        stringsAsFactors = FALSE
      ))
    }
    
    # Add to all parameters
    all_params[[model_name]] <- params
  }
  
  # Combine into one data frame
  param_summary <- do.call(rbind, all_params)
  rownames(param_summary) <- NULL
  
  return(param_summary)
}

#' Fit discrete state transition models for chromosome evolution
#' 
#' Fits discrete Markov models for chromosome number evolution using
#' different rate matrices and transition constraints
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param model_types Vector of model types to fit: "ER", "SYM", "ARD", "custom"
#' @param custom_matrix Custom rate matrix for transitions (for "custom" model)
#' @param max_difference Maximum difference in chromosome number to allow transitions
#' @param allow_jumps Whether to allow transitions between non-adjacent states
#' @param criterion Model selection criterion: "AIC", "BIC", "AICc"
#' @param n_cores Number of cores for parallel model fitting
#' @return List with discrete model fitting results
#' @export
fit_discrete_models <- function(tree, 
                              chr_counts,
                              model_types = c("ER", "SYM", "ARD"),
                              custom_matrix = NULL,
                              max_difference = 2,
                              allow_jumps = FALSE,
                              criterion = "AICc",
                              n_cores = NULL) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Match taxa between tree and data
  shared_taxa <- intersect(tree$tip.label, names(chr_counts))
  
  if(length(shared_taxa) < 4) {
    stop("At least 4 shared taxa required for model fitting")
  }
  
  # Prune tree and data to shared taxa
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, shared_taxa))
  pruned_counts <- chr_counts[pruned_tree$tip.label]
  
  # Round to integers if not already
  pruned_counts <- round(pruned_counts)
  
  # Convert to factor for discrete analysis
  count_levels <- sort(unique(pruned_counts))
  count_factor <- factor(pruned_counts, levels = count_levels)
  names(count_factor) <- names(pruned_counts)
  
  # Check if phytools is available
  if(!requireNamespace("phytools", quietly = TRUE)) {
    stop("phytools package is required for discrete model fitting")
  }
  
  # Set up parallel computing if requested
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Create model specifications including custom models
  model_specs <- list()
  
  for(model_type in model_types) {
    if(model_type == "custom" && is.null(custom_matrix)) {
      warning("custom_matrix must be provided for 'custom' model type. Skipping.")
      next
    }
    
    if(model_type %in% c("ER", "SYM", "ARD")) {
      # Standard models
      model_specs[[model_type]] <- list(
        type = model_type,
        matrix = NULL  # Use default matrix from fitDiscrete
      )
    } else if(model_type == "custom") {
      # Custom rate matrix
      model_specs[[model_type]] <- list(
        type = "custom",
        matrix = custom_matrix
      )
    } else if(model_type == "adjacent") {
      # Only allow transitions between adjacent states
      # Create adjacent-only matrix
      adj_matrix <- create_adjacent_matrix(count_levels, max_difference)
      model_specs[["adjacent"]] <- list(
        type = "custom",
        matrix = adj_matrix
      )
    } else if(model_type == "asymmetric") {
      # Asymmetric rates for increases vs. decreases
      asym_matrix <- create_asymmetric_matrix(count_levels, max_difference)
      model_specs[["asymmetric"]] <- list(
        type = "custom",
        matrix = asym_matrix
      )
    }
  }
  
  # If no valid models, return error
  if(length(model_specs) == 0) {
    stop("No valid models to fit")
  }
  
  # Fit models using fitDiscrete from geiger
  if(n_cores > 1) {
    # Set up cluster
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("pruned_tree", "count_factor"), envir = environment())
    
    # Load required packages on all nodes
    parallel::clusterEvalQ(cl, {
      library(geiger)
    })
    
    # Fit models in parallel
    model_fits <- parallel::parLapply(cl, model_specs, function(spec) {
      tryCatch({
        fit <- geiger::fitDiscrete(pruned_tree, count_factor, model = spec$type, 
                                rate.mat = spec$matrix, ncores = 1)
        return(fit)
      }, error = function(e) {
        message(paste("Error fitting", spec$type, "model:", e$message))
        return(NULL)
      })
    })
  } else {
    # Fit models sequentially
    model_fits <- list()
    for(model_name in names(model_specs)) {
      spec <- model_specs[[model_name]]
      model_fits[[model_name]] <- tryCatch({
        geiger::fitDiscrete(pruned_tree, count_factor, model = spec$type,
                         rate.mat = spec$matrix)
      }, error = function(e) {
        message(paste("Error fitting", spec$type, "model:", e$message))
        return(NULL)
      })
    }
  }
  
  # Remove failed fits
  model_fits <- model_fits[!sapply(model_fits, is.null)]
  
  if(length(model_fits) == 0) {
    stop("All discrete models failed to fit")
  }
  
  # Extract model comparison metrics
  model_comparison <- data.frame(
    Model = names(model_fits),
    LogLik = sapply(model_fits, function(x) x$opt$lnL),
    AIC = sapply(model_fits, function(x) x$opt$aic),
    Parameters = sapply(model_fits, function(x) length(x$opt$index.matrix[!is.na(x$opt$index.matrix)])),
    stringsAsFactors = FALSE
  )
  
  # Calculate AICc and BIC
  n_samples <- length(pruned_counts)
  model_comparison$AICc <- model_comparison$AIC + 
                         (2 * model_comparison$Parameters * (model_comparison$Parameters + 1)) / 
                         (n_samples - model_comparison$Parameters - 1)
  model_comparison$BIC <- model_comparison$AIC + model_comparison$Parameters * log(n_samples) - 
                         2 * model_comparison$Parameters
  
  # Calculate delta values and weights
  if(criterion == "AIC") {
    model_comparison$Delta <- model_comparison$AIC - min(model_comparison$AIC)
    model_comparison$Weight <- exp(-0.5 * model_comparison$Delta) / sum(exp(-0.5 * model_comparison$Delta))
    model_comparison <- model_comparison[order(model_comparison$AIC), ]
  } else if(criterion == "AICc") {
    model_comparison$Delta <- model_comparison$AICc - min(model_comparison$AICc)
    model_comparison$Weight <- exp(-0.5 * model_comparison$Delta) / sum(exp(-0.5 * model_comparison$Delta))
    model_comparison <- model_comparison[order(model_comparison$AICc), ]
  } else if(criterion == "BIC") {
    model_comparison$Delta <- model_comparison$BIC - min(model_comparison$BIC)
    model_comparison$Weight <- exp(-0.5 * model_comparison$Delta) / sum(exp(-0.5 * model_comparison$Delta))
    model_comparison <- model_comparison[order(model_comparison$BIC), ]
  } else {
    stop("criterion must be one of: 'AIC', 'AICc', 'BIC'")
  }
  
  # Record best model
  best_model <- model_comparison$Model[1]
  
  # Create results object
  results <- list(
    model_fits = model_fits,
    model_comparison = model_comparison,
    best_model = best_model,
    tree = pruned_tree,
    chr_counts = pruned_counts,
    count_levels = count_levels,
    criterion = criterion
  )
  
  class(results) <- c("chr_discrete_fits", class(results))
  
  return(results)
}

#' Create an adjacent-only transition matrix
#' 
#' @param states Vector of states (chromosome numbers)
#' @param max_diff Maximum allowed difference for transitions
#' @return Matrix with allowed transitions
#' @keywords internal
create_adjacent_matrix <- function(states, max_diff = 1) {
  n_states <- length(states)
  mat <- matrix(0, nrow = n_states, ncol = n_states)
  
  for(i in 1:n_states) {
    for(j in 1:n_states) {
      # Only allow transitions within max_diff
      if(i != j && abs(states[i] - states[j]) <= max_diff) {
        mat[i, j] <- 1
      }
    }
  }
  
  return(mat)
}

#' Create an asymmetric transition matrix
#' 
#' @param states Vector of states (chromosome numbers)
#' @param max_diff Maximum allowed difference for transitions
#' @return Matrix with asymmetric transition costs
#' @keywords internal
create_asymmetric_matrix <- function(states, max_diff = 1) {
  n_states <- length(states)
  mat <- matrix(0, nrow = n_states, ncol = n_states)
  
  for(i in 1:n_states) {
    for(j in 1:n_states) {
      if(i != j) {
        diff <- states[j] - states[i]
        if(abs(diff) <= max_diff) {
          if(diff > 0) {
            # Increase in chromosome number (fission)
            mat[i, j] <- 2
          } else {
            # Decrease in chromosome number (fusion)
            mat[i, j] <- 1
          }
        }
      }
    }
  }
  
  # Find unique non-zero values in the matrix
  unique_vals <- sort(unique(as.vector(mat)[as.vector(mat) > 0]))
  
  # Renumber values from 1 to length(unique_vals)
  for(i in 1:n_states) {
    for(j in 1:n_states) {
      if(mat[i, j] > 0) {
        mat[i, j] <- which(unique_vals == mat[i, j])
      }
    }
  }
  
  return(mat)
}

#===============================================================================
# Model Comparison and Evaluation Functions
#===============================================================================

#' Compare chromosome evolution models using information criteria
#' 
#' Compares different evolutionary models fitted to chromosome count data,
#' calculates model weights, and provides model-averaged parameter estimates
#' 
#' @param model_fits Result from fit_chromosome_models or fit_discrete_models
#' @param criterion Information criterion for model comparison: "AIC", "AICc", "BIC"
#' @param delta_threshold Threshold for models to be included in confidence set
#' @param weight_threshold Cumulative weight threshold for confidence set
#' @return List with model comparison results and confidence set
#' @export
compare_chromosome_models <- function(model_fits, 
                                    criterion = "AICc", 
                                    delta_threshold = 2,
                                    weight_threshold = 0.95) {
  # Validate input
  if(!inherits(model_fits, "chr_model_fits") && !inherits(model_fits, "chr_discrete_fits")) {
    stop("model_fits must be from fit_chromosome_models or fit_discrete_models")
  }
  
  # Use model comparison table if available
  if(!is.null(model_fits$model_comparison)) {
    comparison <- model_fits$model_comparison
  } else {
    stop("model_fits does not contain model comparison information")
  }
  
  # Check criterion
  if(!criterion %in% c("AIC", "AICc", "BIC")) {
    stop("criterion must be one of: 'AIC', 'AICc', 'BIC'")
  }
  
  # Calculate delta values and weights if not already present
  if(!all(c("Delta", "Weight") %in% colnames(comparison))) {
    if(criterion == "AIC") {
      comparison$Delta <- comparison$AIC - min(comparison$AIC)
      comparison$Weight <- exp(-0.5 * comparison$Delta) / sum(exp(-0.5 * comparison$Delta))
    } else if(criterion == "AICc") {
      comparison$Delta <- comparison$AICc - min(comparison$AICc)
      comparison$Weight <- exp(-0.5 * comparison$Delta) / sum(exp(-0.5 * comparison$Delta))
    } else if(criterion == "BIC") {
      comparison$Delta <- comparison$BIC - min(comparison$BIC)
      comparison$Weight <- exp(-0.5 * comparison$Delta) / sum(exp(-0.5 * comparison$Delta))
    }
    
    # Sort by criterion
    comparison <- comparison[order(comparison[[criterion]]), ]
  }
  
  # Determine confidence set based on delta values
  confidence_set <- comparison$Model[comparison$Delta <= delta_threshold]
  
  # Also determine a confidence set based on cumulative weight
  cum_weight <- cumsum(comparison$Weight)
  weight_set <- comparison$Model[cum_weight <= weight_threshold]
  
  # If weight_set is empty (first model already exceeds threshold), take first model
  if(length(weight_set) == 0) {
    weight_set <- comparison$Model[1]
  }
  
  # Combine confidence sets (union)
  combined_set <- unique(c(confidence_set, weight_set))
  
  # Extract evidence ratios
  best_weight <- comparison$Weight[1]
  comparison$EvidenceRatio <- best_weight / comparison$Weight
  
  # Calculate relative likelihoods
  comparison$RelativeLikelihood <- exp(-0.5 * comparison$Delta)
  
  # Record results
  results <- list(
    comparison_table = comparison,
    criterion = criterion,
    best_model = comparison$Model[1],
    delta_confidence_set = confidence_set,
    weight_confidence_set = weight_set,
    combined_confidence_set = combined_set,
    cumulative_weight = cum_weight
  )
  
  # Calculate model-averaged parameters if available (only for continuous models)
  if(inherits(model_fits, "chr_model_fits") && !is.null(model_fits$param_summary)) {
    # Get parameter summary
    param_summary <- model_fits$param_summary
    
    # Add weights
    param_df <- merge(param_summary, 
                    data.frame(Model = comparison$Model, Weight = comparison$Weight),
                    by = "Model")
    
    # Group by parameter and calculate weighted average
    if(requireNamespace("dplyr", quietly = TRUE)) {
      avg_params <- param_df %>%
        group_by(Parameter) %>%
        summarize(ModelAveraged = sum(Value * Weight),
                WeightedSD = sqrt(sum(Weight * (Value - sum(Value * Weight))^2)),
                .groups = "drop")
      
      results$model_averaged_params <- avg_params
    }
  }
  
  # Add class
  class(results) <- c("chr_model_comparison", class(results))
  
  return(results)
}

#' Perform model adequacy tests for chromosome evolution models
#' 
#' Tests the adequacy of fitted models by comparing observed data to
#' simulations under the fitted model parameters
#' 
#' @param model_fits Result from fit_chromosome_models or fit_discrete_models
#' @param model_name Name of model to test (if NULL, uses best model)
#' @param test_statistics Vector of test statistics to compute
#' @param n_simulations Number of simulations for parametric bootstrapping
#' @param alpha Significance level for tests
#' @param n_cores Number of cores for parallel processing
#' @return List with model adequacy test results
#' @export
test_model_adequacy <- function(model_fits,
                              model_name = NULL,
                              test_statistics = c("mean", "variance", "skewness", "range"),
                              n_simulations = 1000,
                              alpha = 0.05,
                              n_cores = NULL) {
  # Validate input
  if(!inherits(model_fits, "chr_model_fits") && !inherits(model_fits, "chr_discrete_fits")) {
    stop("model_fits must be from fit_chromosome_models or fit_discrete_models")
  }
  
  # Set model to test
  if(is.null(model_name)) {
    model_name <- model_fits$best_model
  } else if(!model_name %in% names(model_fits$model_fits)) {
    stop(paste("Model", model_name, "not found in model_fits"))
  }
  
  # Extract tree and data
  tree <- model_fits$tree
  chr_counts <- model_fits$chr_counts
  
  # Extract parameters from the model
  model_fit <- model_fits$model_fits[[model_name]]
  
  # Set up parallel computing if requested
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Calculate observed test statistics
  obs_stats <- calculate_test_statistics(chr_counts, test_statistics)
  
  # Simulate data under the model
  if(inherits(model_fits, "chr_model_fits")) {
    # Continuous model
    sim_results <- simulate_from_continuous_model(
      tree = tree,
      model_fit = model_fit,
      model_name = model_name,
      n_simulations = n_simulations,
      n_cores = n_cores
    )
  } else if(inherits(model_fits, "chr_discrete_fits")) {
    # Discrete model
    sim_results <- simulate_from_discrete_model(
      tree = tree,
      model_fit = model_fit,
      model_name = model_name,
      count_levels = model_fits$count_levels,
      n_simulations = n_simulations,
      n_cores = n_cores
    )
  } else {
    stop("Unrecognized model type")
  }
  
  # Calculate test statistics for simulated data
  sim_stats <- t(sapply(sim_results, function(x) {
    calculate_test_statistics(x, test_statistics)
  }))
  
  # Calculate p-values (two-tailed)
  p_values <- sapply(test_statistics, function(stat) {
    # Calculate proportion of simulations with more extreme values
    sim_values <- sim_stats[, stat]
    obs_value <- obs_stats[stat]
    
    # Calculate two-tailed p-value
    if(obs_value >= mean(sim_values)) {
      p_value <- 2 * sum(sim_values >= obs_value) / length(sim_values)
    } else {
      p_value <- 2 * sum(sim_values <= obs_value) / length(sim_values)
    }
    
    # Ensure p-value is not greater than 1
    min(p_value, 1)
  })
  
  # Determine if model is adequate
  is_adequate <- all(p_values >= alpha)
  
  # Calculate quantiles for each statistic
  quantiles <- t(apply(sim_stats, 2, function(x) {
    quantile(x, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
  }))
  
  # Create summary table
  summary_table <- data.frame(
    Statistic = test_statistics,
    Observed = sapply(test_statistics, function(x) obs_stats[x]),
    PValue = p_values,
    Q2.5 = quantiles[, "2.5%"],
    Q5 = quantiles[, "5%"],
    Median = quantiles[, "50%"],
    Q95 = quantiles[, "95%"],
    Q97.5 = quantiles[, "97.5%"],
    stringsAsFactors = FALSE
  )
  
  # Add pass/fail column
  summary_table$Pass <- summary_table$PValue >= alpha
  
  # Create results object
  results <- list(
    model_name = model_name,
    summary_table = summary_table,
    observed_stats = obs_stats,
    simulated_stats = sim_stats,
    p_values = p_values,
    alpha = alpha,
    is_adequate = is_adequate,
    n_simulations = n_simulations,
    test_statistics = test_statistics
  )
  
  class(results) <- c("chr_model_adequacy", class(results))
  
  return(results)
}

#' Calculate test statistics for chromosome count data
#' 
#' @param chr_counts Vector of chromosome counts
#' @param statistics Vector of statistic names to calculate
#' @return Named vector of test statistics
#' @keywords internal
calculate_test_statistics <- function(chr_counts, statistics) {
  # Initialize results
  result <- numeric(length(statistics))
  names(result) <- statistics
  
  # Calculate requested statistics
  for(stat in statistics) {
    if(stat == "mean") {
      result[stat] <- mean(chr_counts)
    } else if(stat == "variance") {
      result[stat] <- var(chr_counts)
    } else if(stat == "range") {
      result[stat] <- max(chr_counts) - min(chr_counts)
    } else if(stat == "skewness") {
      # Calculate skewness
      m3 <- mean((chr_counts - mean(chr_counts))^3)
      s3 <- sd(chr_counts)^3
      result[stat] <- m3 / s3
    } else if(stat == "kurtosis") {
      # Calculate excess kurtosis
      m4 <- mean((chr_counts - mean(chr_counts))^4)
      s4 <- sd(chr_counts)^4
      result[stat] <- m4 / s4 - 3
    } else if(stat == "median") {
      result[stat] <- median(chr_counts)
    } else if(stat == "CV") {
      # Coefficient of variation
      result[stat] <- sd(chr_counts) / mean(chr_counts)
    } else {
      warning(paste("Unknown statistic:", stat))
      result[stat] <- NA
    }
  }
  
  return(result)
}

#' Simulate data from a continuous evolutionary model
#' 
#' @param tree Phylogenetic tree
#' @param model_fit Fitted model from fitContinuous
#' @param model_name Name of the model
#' @param n_simulations Number of simulations
#' @param n_cores Number of cores for parallel processing
#' @return List of simulated datasets
#' @keywords internal
simulate_from_continuous_model <- function(tree, model_fit, model_name, n_simulations, n_cores) {
  # Get model parameters
  params <- list(sig2 = model_fit$opt$sigsq)
  
  # Add model-specific parameters
  if(model_name == "OU") {
    params$alpha <- model_fit$opt$alpha
  } else if(model_name == "EB") {
    params$a <- model_fit$opt$a
  } else if(model_name == "lambda") {
    params$lambda <- model_fit$opt$lambda
  } else if(model_name == "kappa") {
    params$kappa <- model_fit$opt$kappa
  } else if(model_name == "delta") {
    params$delta <- model_fit$opt$delta
  }
  
  # Simulate data using geiger::simulate.xxx functions
  if(n_cores > 1) {
    # Set up cluster
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("tree", "model_name", "params"), envir = environment())
    
    # Load required packages on all nodes
    parallel::clusterEvalQ(cl, {
      library(geiger)
    })
    
    # Simulate in parallel
    sim_results <- parallel::parLapply(cl, 1:n_simulations, function(i) {
      # Use appropriate simulation function based on model
      if(model_name == "BM") {
        sim <- geiger::sim.char(tree, params$sig2, model = "BM")[1, ]
      } else if(model_name == "OU") {
        sim <- geiger::sim.char(tree, params$sig2, model = "OU", alpha = params$alpha)[1, ]
      } else if(model_name == "EB") {
        sim <- geiger::sim.char(tree, params$sig2, model = "EB", a = params$a)[1, ]
      } else if(model_name == "lambda") {
        # Modify tree with lambda
        tree_lambda <- tree
        tree_lambda$edge.length <- tree_lambda$edge.length * params$lambda
        sim <- geiger::sim.char(tree_lambda, params$sig2, model = "BM")[1, ]
      } else if(model_name == "kappa") {
        # Modify tree with kappa
        tree_kappa <- tree
        tree_kappa$edge.length <- tree_kappa$edge.length ^ params$kappa
        sim <- geiger::sim.char(tree_kappa, params$sig2, model = "BM")[1, ]
      } else if(model_name == "delta") {
        # Modify tree with delta
        tree_delta <- tree
        tree_delta$edge.length <- tree_delta$edge.length ^ params$delta
        sim <- geiger::sim.char(tree_delta, params$sig2, model = "BM")[1, ]
      } else {
        # Default to BM
        sim <- geiger::sim.char(tree, params$sig2, model = "BM")[1, ]
      }
      
      # Round to integers if the original model was discrete
      return(sim)
    })
  } else {
    # Simulate sequentially
    sim_results <- vector("list", n_simulations)
    for(i in 1:n_simulations) {
      # Use appropriate simulation function based on model
      if(model_name == "BM") {
        sim_results[[i]] <- geiger::sim.char(tree, params$sig2, model = "BM")[1, ]
      } else if(model_name == "OU") {
        sim_results[[i]] <- geiger::sim.char(tree, params$sig2, model = "OU", alpha = params$alpha)[1, ]
      } else if(model_name == "EB") {
        sim_results[[i]] <- geiger::sim.char(tree, params$sig2, model = "EB", a = params$a)[1, ]
      } else if(model_name == "lambda") {
        # Modify tree with lambda
        tree_lambda <- tree
        tree_lambda$edge.length <- tree_lambda$edge.length * params$lambda
        sim_results[[i]] <- geiger::sim.char(tree_lambda, params$sig2, model = "BM")[1, ]
      } else if(model_name == "kappa") {
        # Modify tree with kappa
        tree_kappa <- tree
        tree_kappa$edge.length <- tree_kappa$edge.length ^ params$kappa
        sim_results[[i]] <- geiger::sim.char(tree_kappa, params$sig2, model = "BM")[1, ]
      } else if(model_name == "delta") {
        # Modify tree with delta
        tree_delta <- tree
        tree_delta$edge.length <- tree_delta$edge.length ^ params$delta
        sim_results[[i]] <- geiger::sim.char(tree_delta, params$sig2, model = "BM")[1, ]
      } else {
        # Default to BM
        sim_results[[i]] <- geiger::sim.char(tree, params$sig2, model = "BM")[1, ]
      }
    }
  }
  
  return(sim_results)
}

#' Simulate data from a discrete evolutionary model
#' 
#' @param tree Phylogenetic tree
#' @param model_fit Fitted model from fitDiscrete
#' @param model_name Name of the model
#' @param count_levels Vector of possible chromosome count values
#' @param n_simulations Number of simulations
#' @param n_cores Number of cores for parallel processing
#' @return List of simulated datasets
#' @keywords internal
simulate_from_discrete_model <- function(tree, model_fit, model_name, count_levels, n_simulations, n_cores) {
  # Get transition matrix from model fit
  Q <- model_fit$opt$Q
  
  # Check if phytools is available
  if(!requireNamespace("phytools", quietly = TRUE)) {
    stop("phytools package is required for simulating from discrete models")
  }
  
  # Simulate data using sim.history function in phytools
  if(n_cores > 1) {
    # Set up cluster
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("tree", "Q", "count_levels"), envir = environment())
    
    # Load required packages on all nodes
    parallel::clusterEvalQ(cl, {
      library(phytools)
    })
    
    # Simulate in parallel
    sim_results <- parallel::parLapply(cl, 1:n_simulations, function(i) {
      # Simulate character history
      sim <- phytools::sim.history(tree, Q, anc = 1)
      
      # Extract tip states
      tip_states <- sim$states[1:length(tree$tip.label)]
      
      # Convert state indices to actual chromosome counts
      chr_counts <- count_levels[tip_states]
      names(chr_counts) <- tree$tip.label
      
      return(chr_counts)
    })
  } else {
    # Simulate sequentially
    sim_results <- vector("list", n_simulations)
    for(i in 1:n_simulations) {
      # Simulate character history
      sim <- phytools::sim.history(tree, Q, anc = 1)
      
      # Extract tip states
      tip_states <- sim$states[1:length(tree$tip.label)]
      
      # Convert state indices to actual chromosome counts
      chr_counts <- count_levels[tip_states]
      names(chr_counts) <- tree$tip.label
      
      sim_results[[i]] <- chr_counts
    }
  }
  
  return(sim_results)
}

#===============================================================================
# Visualization Functions
#===============================================================================

#' Plot model comparison results
#' 
#' Creates a visualization of model comparison metrics for chromosome
#' evolution models, highlighting the best-fitting model
#' 
#' @param model_comparison Result from compare_chromosome_models function
#' @param plot_type Type of plot: "weights", "delta", "both"
#' @param highlight_confidence Whether to highlight the confidence set
#' @param models Vector of model names to include (NULL for all)
#' @param color_palette Color palette for bars
#' @return ggplot object with model comparison visualization
#' @export
plot_model_comparison <- function(model_comparison,
                               plot_type = "both",
                               highlight_confidence = TRUE,
                               models = NULL,
                               color_palette = NULL) {
  # Validate input
  if(!inherits(model_comparison, "chr_model_comparison")) {
    stop("model_comparison must be from compare_chromosome_models function")
  }
  
  # Get comparison table
  comparison <- model_comparison$comparison_table
  
  # Filter to specified models if provided
  if(!is.null(models)) {
    comparison <- comparison[comparison$Model %in% models, ]
    if(nrow(comparison) == 0) {
      stop("No specified models found in comparison table")
    }
  }
  
  # Define confidence set for highlighting
  confidence_set <- model_comparison$combined_confidence_set
  comparison$InConfidenceSet <- comparison$Model %in% confidence_set
  
  # Define color palette if not provided
  if(is.null(color_palette)) {
    if(highlight_confidence) {
      color_palette <- c("steelblue", "gray70")
    } else {
      color_palette <- viridis(nrow(comparison), option = "D")
    }
  }
  
  # Create plots based on type
  if(plot_type == "weights" || plot_type == "both") {
    # Create weight plot
    weight_plot <- ggplot(comparison, aes(x = reorder(Model, -Weight), y = Weight, 
                                       fill = InConfidenceSet)) +
      geom_bar(stat = "identity") +
      labs(title = "Model Weights",
           x = "Model",
           y = paste(model_comparison$criterion, "Weights")) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Add color scale based on whether highlighting confidence set
    if(highlight_confidence) {
      weight_plot <- weight_plot +
        scale_fill_manual(values = color_palette, 
                        labels = c("Outside Confidence Set", "In Confidence Set"),
                        name = "Confidence Set")
    } else {
      weight_plot <- weight_plot +
        scale_fill_viridis_d(option = "D", guide = "none")
    }
  }
  
  if(plot_type == "delta" || plot_type == "both") {
    # Create delta plot
    delta_plot <- ggplot(comparison, aes(x = reorder(Model, Delta), y = Delta, 
                                      fill = InConfidenceSet)) +
      geom_bar(stat = "identity") +
      labs(title = paste("Model", model_comparison$criterion, "Differences"),
           x = "Model",
           y = paste("Delta", model_comparison$criterion)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Add horizontal line at delta threshold
    delta_plot <- delta_plot +
      geom_hline(yintercept = 2, linetype = "dashed", color = "red")
    
    # Add color scale based on whether highlighting confidence set
    if(highlight_confidence) {
      delta_plot <- delta_plot +
        scale_fill_manual(values = color_palette, 
                        labels = c("Outside Confidence Set", "In Confidence Set"),
                        name = "Confidence Set")
    } else {
      delta_plot <- delta_plot +
        scale_fill_viridis_d(option = "D", guide = "none")
    }
  }
  
  # Return plots based on type
  if(plot_type == "weights") {
    return(weight_plot)
  } else if(plot_type == "delta") {
    return(delta_plot)
  } else {
    # Combine plots if patchwork is available
    if(requireNamespace("patchwork", quietly = TRUE)) {
      combined_plot <- weight_plot + delta_plot + patchwork::plot_layout(guides = "collect")
      return(combined_plot)
    } else {
      return(list(weight_plot = weight_plot, delta_plot = delta_plot))
    }
  }
}

#' Plot model adequacy test results
#' 
#' Creates a visualization of model adequacy test results showing
#' observed vs. simulated distributions of test statistics
#' 
#' @param model_adequacy Result from test_model_adequacy function
#' @param statistics Vector of statistics to plot (NULL for all)
#' @param plot_type Plot type: "histogram", "density", "boxplot"
#' @param add_observed Whether to show observed value
#' @param ncol Number of columns for grid arrangement
#' @return ggplot object with model adequacy visualization
#' @export
plot_model_adequacy <- function(model_adequacy,
                              statistics = NULL,
                              plot_type = "histogram",
                              add_observed = TRUE,
                              ncol = 2) {
  # Validate input
  if(!inherits(model_adequacy, "chr_model_adequacy")) {
    stop("model_adequacy must be from test_model_adequacy function")
  }
  
  # Get test statistics
  if(is.null(statistics)) {
    statistics <- model_adequacy$test_statistics
  } else if(!all(statistics %in% model_adequacy$test_statistics)) {
    missing_stats <- setdiff(statistics, model_adequacy$test_statistics)
    stop(paste("Statistics not found in model_adequacy:", paste(missing_stats, collapse = ", ")))
  }
  
  # Get data
  sim_stats <- model_adequacy$simulated_stats
  obs_stats <- model_adequacy$observed_stats
  
  # Create plots for each statistic
  plots <- list()
  
  for(stat in statistics) {
    # Create data frame for plotting
    plot_data <- data.frame(
      Statistic = stat,
      Value = sim_stats[, stat],
      Type = "Simulated"
    )
    
    # Create title with p-value
    p_value <- model_adequacy$p_values[stat]
    plot_title <- paste0(stat, " (p = ", round(p_value, 3), ")")
    
    # Create plot
    if(plot_type == "histogram") {
      p <- ggplot(plot_data, aes(x = Value)) +
        geom_histogram(fill = "steelblue", color = "black", alpha = 0.7, bins = 30) +
        labs(title = plot_title,
             x = stat,
             y = "Frequency")
    } else if(plot_type == "density") {
      p <- ggplot(plot_data, aes(x = Value)) +
        geom_density(fill = "steelblue", color = "black", alpha = 0.7) +
        labs(title = plot_title,
             x = stat,
             y = "Density")
    } else if(plot_type == "boxplot") {
      p <- ggplot(plot_data, aes(y = Value, x = 1)) +
        geom_boxplot(fill = "steelblue", color = "black") +
        labs(title = plot_title,
             y = stat) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    } else {
      stop("Unsupported plot_type. Use 'histogram', 'density', or 'boxplot'")
    }
    
    # Add observed value
    if(add_observed) {
      obs_value <- obs_stats[stat]
      
      if(plot_type %in% c("histogram", "density")) {
        p <- p + geom_vline(xintercept = obs_value, color = "red", size = 1, linetype = "dashed") +
          annotate("text", x = obs_value, y = 0, label = "Observed", 
                 color = "red", angle = 90, hjust = -0.2, vjust = -0.5)
      } else if(plot_type == "boxplot") {
        p <- p + geom_point(aes(x = 1, y = obs_value), color = "red", size = 3) +
          annotate("text", x = 1.2, y = obs_value, label = "Observed", 
                 color = "red", hjust = 0)
      }
    }
    
    # Add theme
    p <- p + theme_minimal() +
      theme(plot.title = element_text(size = 10))
    
    # Store plot
    plots[[stat]] <- p
  }
  
  # Combine plots
  if(requireNamespace("patchwork", quietly = TRUE) && length(plots) > 1) {
    # Combine using patchwork
    combined_plot <- patchwork::wrap_plots(plots, ncol = ncol) +
      patchwork::plot_annotation(
        title = paste("Model Adequacy Tests for", model_adequacy$model_name, "Model"),
        subtitle = paste(model_adequacy$n_simulations, "simulations")
      )
    return(combined_plot)
  } else if(requireNamespace("gridExtra", quietly = TRUE) && length(plots) > 1) {
    # Combine using gridExtra as fallback
    combined_plot <- gridExtra::grid.arrange(
      grobs = plots, 
      ncol = ncol,
      top = paste("Model Adequacy Tests for", model_adequacy$model_name, "Model")
    )
    return(combined_plot)
  } else if(length(plots) == 1) {
    # Return single plot
    return(plots[[1]])
  } else {
    # Return list of plots
    return(plots)
  }
}

#' Plot ancestral state reconstructions under different models
#' 
#' Creates a visualization comparing ancestral state reconstructions
#' under different evolutionary models
#' 
#' @param tree Phylogenetic tree
#' @param model_fits List of model fits from different methods
#' @param model_names Vector of model names to include
#' @param layout Tree layout: "rectangular", "circular", "fan"
#' @param show_labels Whether to show node labels
#' @param highlight_nodes Vector of nodes to highlight
#' @param color_scheme Color scheme for mapping values
#' @return ggplot object with ancestral state comparison
#' @export
plot_model_reconstructions <- function(tree,
                                     model_fits,
                                     model_names = NULL,
                                     layout = "rectangular",
                                     show_labels = TRUE,
                                     highlight_nodes = NULL,
                                     color_scheme = "viridis") {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggtree package is required for this function")
  }
  
  # Check if model_fits is a list or single object
  if(is.list(model_fits) && !inherits(model_fits, "chr_model_fits")) {
    if(is.null(model_names)) {
      model_names <- names(model_fits)
      if(is.null(model_names)) {
        model_names <- paste("Model", 1:length(model_fits))
      }
    } else if(length(model_names) != length(model_fits)) {
      warning("Length of model_names does not match length of model_fits. Using default names.")
      model_names <- paste("Model", 1:length(model_fits))
    }
  } else {
    # Convert single model fit to list
    model_fits <- list(model_fits)
    if(is.null(model_names)) {
      if(inherits(model_fits[[1]], "chr_model_fits")) {
        model_names <- model_fits[[1]]$best_model
      } else {
        model_names <- "Model 1"
      }
    }
  }
  
  # Reconstruct ancestral states for each model
  reconstructions <- list()
  
  for(i in 1:length(model_fits)) {
    fit <- model_fits[[i]]
    name <- model_names[i]
    
    # Extract best model if fit is a chr_model_fits object
    if(inherits(fit, "chr_model_fits")) {
      model <- fit$best_model
      best_fit <- fit$model_fits[[model]]
    } else {
      best_fit <- fit
    }
    
    # Get chromosome counts
    chr_counts <- fit$chr_counts
    
    # Reconstruct ancestral states
    if(inherits(fit, "chr_model_fits")) {
      # Continuous model
      ace_result <- ace(chr_counts, tree, type = "continuous", method = "ML", model = best_fit$model)
      
      # Extract ancestral states
      n_tips <- length(tree$tip.label)
      internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
      
      # Combine tip and node states
      all_states <- c(chr_counts, ace_result$ace)
      names(all_states) <- c(1:n_tips, internal_nodes)
    } else if(inherits(fit, "chr_discrete_fits")) {
      # Discrete model
      # Convert to factor for reconstruction
      count_levels <- fit$count_levels
      count_factor <- factor(chr_counts, levels = count_levels)
      
      # Reconstruct using fitDiscrete
      anc_states <- geiger::ancThresh(tree = tree, 
                                    x = count_factor,
                                    pars = best_fit$opt$pars,
                                    Q = best_fit$opt$Q)
      
      # Convert back to numeric
      n_tips <- length(tree$tip.label)
      internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
      
      # Extract ML estimates for internal nodes
      node_states <- sapply(anc_states, function(x) count_levels[which.max(x)])
      
      # Combine tip and node states
      all_states <- c(chr_counts, node_states)
      names(all_states) <- c(1:n_tips, internal_nodes)
    } else {
      # Use ace as fallback
      ace_result <- ace(chr_counts, tree, type = "continuous", method = "ML")
      
      # Extract ancestral states
      n_tips <- length(tree$tip.label)
      internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
      
      # Combine tip and node states
      all_states <- c(chr_counts, ace_result$ace)
      names(all_states) <- c(1:n_tips, internal_nodes)
    }
    
    reconstructions[[name]] <- all_states
  }
  
  # Create data frame for visualization
  plot_data <- data.frame(
    node = as.numeric(row.names(data.frame(reconstructions))),
    stringsAsFactors = FALSE
  )
  
  # Add reconstructed values for each model
  for(name in model_names) {
    plot_data[[name]] <- reconstructions[[name]]
  }
  
  # Create tree visualization using ggtree
  if(layout == "circular") {
    p <- ggtree::ggtree(tree, layout = "circular")
  } else if(layout == "fan") {
    p <- ggtree::ggtree(tree, layout = "fan")
  } else {
    p <- ggtree::ggtree(tree)
  }
  
  # Add node data
  p <- p %<+% plot_data
  
  # Create list to store individual tree plots
  tree_plots <- list()
  
  # Plot each reconstruction separately
  for(name in model_names) {
    # Create tree with mapped values
    tree_plot <- p + 
      ggtree::geom_nodepoint(aes_string(color = name, size = name), alpha = 0.8) +
      ggtree::geom_tippoint(aes_string(color = name, size = name), alpha = 0.8) +
      scale_color_viridis_c(name = "Chr Number", option = color_scheme) +
      scale_size_continuous(range = c(1, 4), guide = "none") +
      ggtree::theme_tree2() +
      labs(title = name)
    
    # Add tip labels if requested
    if(show_labels) {
      tree_plot <- tree_plot + ggtree::geom_tiplab(size = 3, hjust = -0.1)
    }
    
    # Highlight specified nodes if any
    if(!is.null(highlight_nodes)) {
      # Find the rows corresponding to highlight nodes
      highlight_idx <- which(plot_data$node %in% highlight_nodes)
      
      if(length(highlight_idx) > 0) {
        highlight_data <- plot_data[highlight_idx, ]
        
        # Add highlight points
        tree_plot <- tree_plot + 
          geom_point2(data = highlight_data, 
                    aes(x = x, y = y), 
                    shape = 21, size = 5, color = "black", fill = NA)
      }
    }
    
    tree_plots[[name]] <- tree_plot
  }
  
  # Combine plots if patchwork is available
  if(requireNamespace("patchwork", quietly = TRUE) && length(tree_plots) > 1) {
    combined_plot <- patchwork::wrap_plots(tree_plots) +
      patchwork::plot_annotation(
        title = "Ancestral State Reconstruction Comparison",
        theme = theme(plot.title = element_text(hjust = 0.5))
      )
    return(combined_plot)
  } else if(length(tree_plots) == 1) {
    # Return single plot
    return(tree_plots[[1]])
  } else {
    # Return list of plots
    return(tree_plots)
  }
}
