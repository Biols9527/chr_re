#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Maximum Likelihood Module
# Author: Bioinformatics Team
# Date: 2025-03-18
# Description: Implements maximum likelihood methods for ancestral chromosome
#              number reconstruction and evolutionary model comparison
#===============================================================================

suppressPackageStartupMessages({
  library(ape)       # For phylogenetic operations
  library(phangorn)  # For ML reconstruction
  library(geiger)    # For fitContinuous
  library(stats4)    # For mle function
  library(numDeriv)  # For gradient and hessian calculations
})

#===============================================================================
# Maximum Likelihood Reconstruction Functions
#===============================================================================

#' Initialize ML Run (Internal Helper)
#'
#' Performs initial input validation, tree pruning, and data ordering.
#'
#' @param tree Phylogenetic tree object.
#' @param chr_counts Chromosome counts, named vector.
#' @param model Evolutionary model string.
#' @return A list: `list(pruned_tree = pruned_tree, ordered_counts = ordered_counts)`.
#' @keywords internal
initialize_ml_run_internal <- function(tree, chr_counts, model) {
  # Check input
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  if (is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector with species names")
  }
  supported_models <- c("BM", "OU", "EB", "ACDC", "lambda")
  if (!model %in% supported_models) {
    stop(paste("Unsupported model. Choose from:", paste(supported_models, collapse = ", ")))
  }

  message(sprintf("Initializing ML reconstruction with %s model...", model)) # Message moved here

  # Ensure tree and data contain same species
  common_species <- intersect(names(chr_counts), tree$tip.label)
  if (length(common_species) < 3) {
    stop("Fewer than 3 species in common, cannot perform ancestral state reconstruction")
  }
  
  # Prune tree to match data
  pruned_tree <- ape::keep.tip(tree, common_species)
  
  # Order data to match tree
  ordered_counts <- chr_counts[pruned_tree$tip.label]
  
  return(list(pruned_tree = pruned_tree, ordered_counts = ordered_counts))
}

#' Estimate Ancestral States using ML (Internal Helper)
#'
#' Performs ancestral state estimation based on the chosen model and fitted parameters.
#'
#' @param pruned_tree Pruned phylogenetic tree object.
#' @param ordered_counts Ordered chromosome counts for tips.
#' @param model Evolutionary model string.
#' @param model_fit The result from `fit_evolutionary_model`.
#' @return The `ace_result` object from `ape::ace` or `reconstruct_ou_model`.
#' @keywords internal
estimate_ancestral_states_ml_internal <- function(pruned_tree, ordered_counts, model, model_fit) {
  ace_result <- NULL
  
  # For different models, use appropriate method
  if (model == "BM") {
    message("Estimating ancestral states using Brownian Motion model (ape::ace)...")
    ace_result <- ape::ace(ordered_counts, pruned_tree, type = "continuous", method = "ML")
  } else if (model == "OU") {
    message("Estimating ancestral states using Ornstein-Uhlenbeck model (custom reconstruction)...")
    # reconstruct_ou_model is an existing helper in this file
    ace_result <- reconstruct_ou_model(pruned_tree, ordered_counts, model_fit$parameters$alpha)
  } else if (model %in% c("EB", "ACDC", "lambda")) {
    message(sprintf("Estimating ancestral states using %s model with tree transformation...", model))
    # Use transformations of the tree according to models
    rate_param <- if (model == "EB") model_fit$parameters$r else
                  if (model == "ACDC") model_fit$parameters$beta else
                  model_fit$parameters$lambda
                  
    # Call the moved transform_tree_by_model function
    # Assuming 'transform_tree_by_model' from core/simulation_tools.R is loaded into the environment.
    # If it were in a package or specific environment, it would be e.g., simulation_tools::transform_tree_by_model
    if (!exists("transform_tree_by_model", mode = "function")) {
        stop("transform_tree_by_model function not found. Ensure core/simulation_tools.R is sourced/loaded.")
    }
    transformed_tree <- transform_tree_by_model(pruned_tree, model, rate_param)
    
    # Use standard ACE on transformed tree
    ace_result <- ape::ace(ordered_counts, transformed_tree, type = "continuous", method = "ML")
  } else {
    # This case should ideally be caught by initial validation, but as a safeguard:
    stop(sprintf("Ancestral state estimation for model '%s' is not implemented here.", model))
  }
  
  return(ace_result)
}

#' Format ML Reconstruction Results (Internal Helper)
#'
#' Formats the results from `ace` into the standard output structure.
#'
#' @param ace_result The result from `estimate_ancestral_states_ml_internal`.
#' @param pruned_tree Pruned phylogenetic tree object.
#' @param ordered_counts Ordered chromosome counts for tips.
#' @param model Evolutionary model string.
#' @param model_fit The result from `fit_evolutionary_model`.
#' @param original_chr_counts Original chromosome counts (before pruning/ordering).
#' @param confidence_level Confidence level for CI calculation (currently noted but not used for dynamic CI).
#' @return The final formatted result list.
#' @keywords internal
format_ml_results_internal <- function(ace_result, pruned_tree, ordered_counts, model, 
                                     model_fit, original_chr_counts, confidence_level) {
  
  # Extract ancestral states
  node_ids <- (Ntip(pruned_tree) + 1):(Ntip(pruned_tree) + Nnode(pruned_tree))
  
  # Create ancestral states data frame
  # Note: confidence_level parameter is available. Currently, 1.96 (95% CI) is hardcoded.
  # To use confidence_level dynamically: z_score <- qnorm(1 - (1 - confidence_level) / 2)
  # Then use z_score instead of 1.96.
  ancestors <- data.frame(
    node_id = node_ids,
    state = ace_result$ace,
    ci_lower = ace_result$ace - 1.96 * sqrt(diag(ace_result$var.anc)), # 1.96 for 95% CI
    ci_upper = ace_result$ace + 1.96 * sqrt(diag(ace_result$var.anc)), # 1.96 for 95% CI
    stringsAsFactors = FALSE
  )
  
  # Ensure non-negative lower bound
  ancestors$ci_lower <- pmax(0, ancestors$ci_lower)
  
  message(paste("Completed ancestral chromosome reconstruction formatting for", nrow(ancestors), "nodes."))
  
  # Compile final result
  final_result_structure <- list(
    ancestral_states = ancestors,
    tree = pruned_tree,
    method = "ML",
    model = model,
    model_fit = model_fit,
    tip_states = ordered_counts,
    ace_object = ace_result, # Storing the raw ACE result can be useful
    original_data = list(
      chr_counts = original_chr_counts # Store the original, full chr_counts
    )
  )
  
  # Calculate reconstruction quality metrics using the existing helper
  final_result_structure$quality <- calculate_ml_quality(final_result_structure) # Existing helper
  
  return(final_result_structure)
}

#' Reconstruct ancestral chromosome numbers using maximum likelihood
#' 
#' Applies maximum likelihood inference to reconstruct ancestral chromosome numbers
#' 
#' @param tree Phylogenetic tree object
#' @param chr_counts Chromosome counts, named vector with species names
#' @param model Evolutionary model: "BM" (Brownian Motion), "OU" (Ornstein-Uhlenbeck),
#'              "EB" (Early Burst), "ACDC" (Accelerating/Decelerating), or "lambda"
#' @param estimation_method Parameter estimation method: "reml" or "ml"
#' @param bounded Whether to use bounded optimization for parameters
#' @param confidence_level Confidence level for parameter CIs
#' @return Ancestral chromosome reconstruction results
#' @export
reconstruct_chromosomes_ml <- function(tree, chr_counts, 
                                    model = "BM", 
                                    estimation_method = "reml",
                                    bounded = TRUE,
                                    confidence_level = 0.95) {
  
  # Step 1: Initialization (validation, pruning, ordering)
  # The message "Initializing ML reconstruction with %s model..." is in initialize_ml_run_internal
  init_data <- initialize_ml_run_internal(tree, chr_counts, model)
  # init_data contains $pruned_tree and $ordered_counts
  
  # Step 2: Fit evolutionary model (existing helper)
  # The message "Fitting %s model to chromosome count data..." is in fit_evolutionary_model
  model_fit <- fit_evolutionary_model(
    tree = init_data$pruned_tree, 
    chr_counts = init_data$ordered_counts, 
    model = model, 
    method = estimation_method, 
    bounded = bounded
  )
  # The message "Model fitting completed..." is in fit_evolutionary_model
  
  # Step 3: Estimate ancestral states using the appropriate ML method
  # Messages related to specific estimation methods are in estimate_ancestral_states_ml_internal
  ace_result <- estimate_ancestral_states_ml_internal(
    pruned_tree = init_data$pruned_tree,
    ordered_counts = init_data$ordered_counts,
    model = model,
    model_fit = model_fit
  )
  
  # Step 4: Format results
  # The message "Completed ancestral chromosome reconstruction formatting..." is in format_ml_results_internal
  final_result <- format_ml_results_internal(
    ace_result = ace_result,
    pruned_tree = init_data$pruned_tree,
    ordered_counts = init_data$ordered_counts,
    model = model,
    model_fit = model_fit,
    original_chr_counts = chr_counts, # Pass the original chr_counts
    confidence_level = confidence_level
  )
  
  return(final_result)
}

#' Fit evolutionary model to chromosome count data
#' 
#' Fit continuous trait evolutionary model to chromosome count data
#' 
#' @param tree Phylogenetic tree object
#' @param chr_counts Chromosome counts, named vector
#' @param model Evolutionary model: "BM", "OU", "EB", "ACDC", or "lambda"
#' @param method Parameter estimation method: "reml" or "ml"
#' @param bounded Whether to use bounded optimization
#' @return Fitted model object
#' @keywords internal
fit_evolutionary_model <- function(tree, chr_counts, model = "BM", 
                                method = "reml", bounded = TRUE) {
  # Use geiger's fitContinuous to fit models
  model_name <- switch(model,
                     "BM" = "BM",
                     "OU" = "OU",
                     "EB" = "EB",
                     "ACDC" = "ACDC",
                     "lambda" = "lambda",
                     stop("Unsupported model"))
  
  message(sprintf("Fitting %s model to chromosome count data...", model_name))
  
  # Prepare data for fitContinuous
  if(!is.null(names(chr_counts))) {
    names_match <- all(names(chr_counts) %in% tree$tip.label)
    if(!names_match) {
      warning("Not all names in chr_counts match tree tip labels")
    }
  }
  
  # Fit model using geiger
  tryCatch({
    fit_result <- geiger::fitContinuous(
      phy = tree,
      dat = chr_counts,
      model = model_name,
      method = method,
      bounds = if(bounded) "bounded" else "unbounded"
    )
    
    # Extract parameters
    parameters <- fit_result$opt
    
    message(sprintf("Model fitting completed. Log-likelihood: %.4f, AIC: %.4f", 
                   parameters$lnL, parameters$aic))
    
    # Return complete result
    return(list(
      parameters = parameters,
      full_result = fit_result,
      model = model_name
    ))
    
  }, error = function(e) {
    warning(sprintf("Model fitting failed: %s. Falling back to Brownian Motion model.", e$message))
    
    # Fallback to simple Brownian Motion using ace
    ace_bm <- ape::ace(chr_counts, tree, type = "continuous", method = "ML")
    
    # Construct minimal result
    return(list(
      parameters = list(
        lnL = ace_bm$loglik,
        aic = -2 * ace_bm$loglik + 2 * 2,  # Simple AIC calculation (k=2 for BM)
        sigsq = ace_bm$sigma2,
        z0 = ace_bm$ace[1]  # Root state
      ),
      full_result = NULL,
      model = "BM"
    ))
  })
}

#' Transform phylogenetic tree according to model parameters
#' 
#' Reconstruct ancestral states under the OU model
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param alpha OU model alpha parameter
#' @return Reconstructed ancestral states
#' @keywords internal
reconstruct_ou_model <- function(tree, chr_counts, alpha) {
  # Implement ML reconstruction under OU model
  # This is a simplified version - a full implementation would use a more
  # complete implementation of the OU process
  
  # Fall back to standard ace, which works reasonably well for small alpha
  message("Using approximation for OU model reconstruction")
  
  # Standard ace implementation
  ace_result <- ape::ace(chr_counts, tree, type = "continuous", method = "ML")
  
  # For more accurate results, we would need to properly implement
  # reconstruction under the OU model considering the alpha parameter
  
  return(ace_result)
}

#' Calculate ML reconstruction quality metrics
#' 
#' @param ml_result Maximum likelihood reconstruction result
#' @return Quality assessment metrics
#' @keywords internal
calculate_ml_quality <- function(ml_result) {
  # Extract parameters
  tree <- ml_result$tree
  model_fit <- ml_result$model_fit
  ancestors <- ml_result$ancestral_states
  
  # Calculate confidence interval widths
  ci_widths <- ancestors$ci_upper - ancestors$ci_lower
  mean_ci_width <- mean(ci_widths)
  
  # Calculate normalized AIC (accounting for sample size)
  n_samples <- length(ml_result$tip_states)
  n_params <- length(model_fit$parameters) - 2  # Subtract lnL and AIC from count
  aic <- model_fit$parameters$aic
  aicc <- aic + (2 * n_params * (n_params + 1)) / (n_samples - n_params - 1)
  aic_per_sample <- aic / n_samples
  
  # Root state uncertainty
  root_node <- Ntip(tree) + 1
  root_state_idx <- which(ancestors$node_id == root_node)
  if(length(root_state_idx) > 0) {
    root_uncertainty <- ancestors$ci_upper[root_state_idx] - ancestors$ci_lower[root_state_idx]
  } else {
    root_uncertainty <- NA
  }
  
  # Return quality metrics
  quality <- list(
    mean_ci_width = mean_ci_width,
    max_ci_width = max(ci_widths),
    min_ci_width = min(ci_widths),
    lnL = model_fit$parameters$lnL,
    aic = aic,
    aicc = aicc,
    aic_per_sample = aic_per_sample,
    root_uncertainty = root_uncertainty
  )
  
  return(quality)
}

#===============================================================================
# Model Selection Functions
#===============================================================================

#' Compare multiple evolutionary models for chromosome evolution
#' 
#' Fits and compares different evolutionary models to chromosome count data
#' 
#' @param tree Phylogenetic tree object
#' @param chr_counts Chromosome counts, named vector
#' @param models Vector of models to compare
#' @param criterion Model selection criterion: "aic", "aicc", or "bic"
#' @param method Parameter estimation method: "reml" or "ml"
#' @return Model comparison results
#' @export
compare_chromosome_models <- function(tree, chr_counts, 
                                   models = c("BM", "OU", "EB", "ACDC", "lambda"),
                                   criterion = "aicc",
                                   method = "reml") {
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Ensure valid models
  valid_models <- c("BM", "OU", "EB", "ACDC", "lambda")
  models <- models[models %in% valid_models]
  
  if(length(models) == 0) {
    stop("No valid models specified")
  }
  
  message(sprintf("Comparing %d evolutionary models: %s", 
                 length(models), paste(models, collapse = ", ")))
  
  # Ensure tree and data contain same species
  common_species <- intersect(names(chr_counts), tree$tip.label)
  if(length(common_species) < 3) {
    stop("Fewer than 3 species in common, cannot perform model comparison")
  }
  
  # Prune tree to match data
  pruned_tree <- ape::keep.tip(tree, common_species)
  
  # Order data to match tree
  ordered_counts <- chr_counts[pruned_tree$tip.label]
  
  # Fit all models
  model_fits <- lapply(models, function(model) {
    fit_evolutionary_model(pruned_tree, ordered_counts, model, method = method)
  })
  names(model_fits) <- models
  
  # Extract model selection criteria
  n_samples <- length(ordered_counts)
  
  model_selection <- data.frame(
    model = models,
    lnL = sapply(model_fits, function(x) x$parameters$lnL),
    aic = sapply(model_fits, function(x) x$parameters$aic),
    aicc = sapply(model_fits, function(x) {
      aic <- x$parameters$aic
      n_params <- length(x$parameters) - 2  # Subtract lnL and AIC from count
      aic + (2 * n_params * (n_params + 1)) / (n_samples - n_params - 1)
    }),
    bic = sapply(model_fits, function(x) {
      n_params <- length(x$parameters) - 2  # Subtract lnL and AIC from count
      -2 * x$parameters$lnL + n_params * log(n_samples)
    }),
    stringsAsFactors = FALSE
  )
  
  # Calculate AIC weights
  if(criterion %in% c("aic", "aicc")) {
    criterion_values <- if(criterion == "aic") model_selection$aic else model_selection$aicc
    min_criterion <- min(criterion_values)
    delta_criterion <- criterion_values - min_criterion
    rel_likelihood <- exp(-0.5 * delta_criterion)
    model_selection$weight <- rel_likelihood / sum(rel_likelihood)
  } else if(criterion == "bic") {
    # BIC weights
    min_criterion <- min(model_selection$bic)
    delta_criterion <- model_selection$bic - min_criterion
    rel_likelihood <- exp(-0.5 * delta_criterion)
    model_selection$weight <- rel_likelihood / sum(rel_likelihood)
  }
  
  # Sort by selected criterion
  criterion_col <- match(criterion, c("aic", "aicc", "bic"))
  if(is.na(criterion_col)) criterion_col <- 3  # Default to aicc
  
  model_selection <- model_selection[order(model_selection[, criterion_col + 1]), ]
  
  # Determine best model
  best_model <- model_selection$model[1]
  
  message(sprintf("Model comparison completed. Best model according to %s: %s", 
                 toupper(criterion), best_model))
  
  # Put together comprehensive result
  result <- list(
    model_selection = model_selection,
    best_model = best_model,
    best_model_fit = model_fits[[best_model]],
    model_fits = model_fits,
    criterion = criterion,
    parameter_estimates = lapply(models, function(model) {
      # Extract key parameters for each model
      params <- model_fits[[model]]$parameters
      # Remove standard entries
      params$lnL <- params$aic <- NULL
      return(params)
    }),
    n_samples = n_samples,
    tree = pruned_tree,
    data = ordered_counts
  )
  
  # Generate model comparison plots
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Create model selection plot
    criterion_names <- c("AIC", "AICc", "BIC")
    criterion_name <- criterion_names[match(criterion, c("aic", "aicc", "bic"))]
    
    plot_data <- data.frame(
      Model = model_selection$model,
      Criterion = model_selection[, criterion],
      Weight = model_selection$weight
    )
    
    selection_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(Model, -Criterion), y = Criterion)) +
      ggplot2::geom_bar(stat = "identity", fill = "skyblue") +
      ggplot2::labs(title = paste("Model Selection by", criterion_name),
                   x = "Model", y = criterion_name) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    # Create weights plot
    weights_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(Model, -Weight), y = Weight)) +
      ggplot2::geom_bar(stat = "identity", fill = "salmon") +
      ggplot2::labs(title = paste("Model Weights by", criterion_name),
                   x = "Model", y = "Weight") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    result$plots <- list(
      selection = selection_plot,
      weights = weights_plot
    )
  }
  
  # Generate interpretation
  model_descriptions <- list(
    BM = "Brownian Motion (random walk with constant rate)",
    OU = "Ornstein-Uhlenbeck (random walk with central tendency)",
    EB = "Early Burst (exponential slowing through time)",
    ACDC = "Accelerating/Decelerating (exponential change in rate)",
    lambda = "Pagel's Lambda (phylogenetic signal scaling)"
  )
  
  result$explanation <- sprintf(
    "The best-supported model for chromosome number evolution is %s (%s), with %s = %.2f and weight = %.2f.",
    best_model, model_descriptions[[best_model]], criterion_name, 
    model_selection[1, criterion], model_selection$weight[1]
  )
  
  return(result)
}

#' Simulate chromosome evolution under different models
#' 
#' Simulate chromosome counts under fitted evolutionary models
#' 
#' @param tree Phylogenetic tree object
#' @param model_fits Model comparison results from compare_chromosome_models
#' @param n_simulations Number of simulations per model
#' @param true_counts Optional known true counts for comparison
#' @return Simulation results
#' @export
simulate_chromosome_evolution <- function(tree, model_fits, n_simulations = 100, true_counts = NULL) {
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.list(model_fits) || is.null(model_fits$model_fits)) {
    stop("model_fits must be the result from compare_chromosome_models")
  }
  
  message(sprintf("Simulating chromosome evolution under %d different models...", 
                 length(model_fits$model_fits)))
  
  # Extract models
  models <- names(model_fits$model_fits)
  
  # Initialize results
  sim_results <- list()
  
  # Simulate under each model
  for(model in models) {
    message(sprintf("  Simulating under %s model...", model))
    
    # Extract model parameters
    model_fit <- model_fits$model_fits[[model]]
    model_params <- model_fit$parameters
    
    # Setup simulation parameters
    root_state <- model_params$z0  # Ancestral state at root
    sigma_sq <- model_params$sigsq  # Rate parameter
    
    # Additional model-specific parameters
    extra_params <- list()
    if(model == "OU") {
      extra_params$alpha <- model_params$alpha
      extra_params$theta <- model_params$theta
    } else if(model == "EB") {
      extra_params$r <- model_params$r
    } else if(model == "ACDC") {
      extra_params$beta <- model_params$beta
    } else if(model == "lambda") {
      extra_params$lambda <- model_params$lambda
    }
    
    # Run simulations
    simulations <- replicate(n_simulations, {
      # Prepare parameters for run_chromosome_simulation
      simulation_params <- list(
        sigma_sq = sigma_sq,
        extra_params = extra_params # Contains alpha, theta, r, beta, lambda as needed
      )
      
      # Call the unified simulation function from core/simulation_tools.R
      # Assuming run_chromosome_simulation is available in the environment (e.g. sourced)
      sim_result_obj <- run_chromosome_simulation(
        tree = tree,
        simulation_type = "continuous_fitted",
        model_name = model,
        root_state = root_state,
        params = simulation_params
      )
      
      sim_counts_tips <- sim_result_obj$tip_states
      
      # Ensure non-negative counts (impose floor at 1)
      sim_counts_tips <- pmax(1, round(sim_counts_tips))
      
      return(sim_counts_tips)
    }, simplify = FALSE)
    
    # Store results
    sim_results[[model]] <- list(
      simulations = simulations,
      parameters = model_params
    )
    
    # Compare with true counts if provided
    if(!is.null(true_counts)) {
      common_species <- intersect(names(true_counts), tree$tip.label)
      
      if(length(common_species) > 0) {
        true_subset <- true_counts[common_species]
        
        # Calculate errors for each simulation
        errors <- sapply(simulations, function(sim) {
          sim_subset <- sim[common_species]
          
          # Calculate metrics
          mae <- mean(abs(sim_subset - true_subset))
          rmse <- sqrt(mean((sim_subset - true_subset)^2))
          cor_val <- cor(sim_subset, true_subset, method = "spearman")
          
          return(c(mae = mae, rmse = rmse, correlation = cor_val))
        })
        
        # Summarize errors
        error_summary <- data.frame(
          metric = c("MAE", "RMSE", "Correlation"),
          mean = c(mean(errors["mae", ]), mean(errors["rmse", ]), mean(errors["correlation", ])),
          sd = c(sd(errors["mae", ]), sd(errors["rmse", ]), sd(errors["correlation", ])),
          min = c(min(errors["mae", ]), min(errors["rmse", ]), min(errors["correlation", ])),
          max = c(max(errors["mae", ]), max(errors["rmse", ]), max(errors["correlation", ]))
        )
        
        sim_results[[model]]$error_metrics <- error_summary
      }
    }
  }
  
  # Create overall summary
  if(!is.null(true_counts)) {
    model_performance <- data.frame(
      model = models,
      mae = sapply(models, function(m) sim_results[[m]]$error_metrics$mean[1]),
      rmse = sapply(models, function(m) sim_results[[m]]$error_metrics$mean[2]),
      correlation = sapply(models, function(m) sim_results[[m]]$error_metrics$mean[3]),
      stringsAsFactors = FALSE
    )
    
    # Determine best performing model
    best_model_mae <- model_performance$model[which.min(model_performance$mae)]
    best_model_rmse <- model_performance$model[which.min(model_performance$rmse)]
    best_model_cor <- model_performance$model[which.max(model_performance$correlation)]
    
    message(sprintf("Simulation results: Best model by MAE: %s, by RMSE: %s, by correlation: %s", 
                   best_model_mae, best_model_rmse, best_model_cor))
    
    # Add to results
    sim_results$performance <- list(
      model_performance = model_performance,
      best_model_mae = best_model_mae,
      best_model_rmse = best_model_rmse,
      best_model_cor = best_model_cor
    )
  }
  
  return(sim_results)
}

#===============================================================================
# Model-Based Event Detection Functions
#===============================================================================

#' Detect chromosome events based on model-based ancestral reconstruction
#' 
#' Detects significant chromosome number changes based on ML reconstruction and
#' statistical evidence
#' 
#' @param tree Phylogenetic tree object
#' @param ml_result Maximum likelihood reconstruction result
#' @param fusion_threshold Fusion event threshold
#' @param fission_threshold Fission event threshold
#' @param wgd_threshold Whole genome duplication threshold
#' @param confidence_level Statistical confidence level
#' @return Detected chromosome events
#' @export
detect_events_ml <- function(tree, ml_result,
                          fusion_threshold = 0.67,
                          fission_threshold = 1.5,
                          wgd_threshold = 2.0,
                          confidence_level = 0.95) {
  # Check if reconstruction is ML result
  if(!is.list(ml_result) || is.null(ml_result$ancestral_states) || ml_result$method != "ML") {
    stop("ml_result must be the result from reconstruct_chromosomes_ml")
  }
  
  message("Detecting chromosome events using ML-based ancestral reconstruction...")
  
  # Extract reconstruction and tree
  ancestral_states <- ml_result$ancestral_states
  tip_states <- ml_result$tip_states
  tree <- ml_result$tree
  
  # Create all nodes data frame
  all_nodes <- data.frame(
    node_id = c(1:length(tree$tip.label), ancestral_states$node_id),
    state = c(tip_states, ancestral_states$state),
    ci_lower = c(rep(NA, length(tip_states)), ancestral_states$ci_lower),
    ci_upper = c(rep(NA, length(tip_states)), ancestral_states$ci_upper),
    is_tip = c(rep(TRUE, length(tip_states)), rep(FALSE, nrow(ancestral_states))),
    stringsAsFactors = FALSE
  )
  
  # Prepare edges data
  edges <- tree$edge
  n_edges <- nrow(edges)
  
  events <- data.frame(
    edge_id = 1:n_edges,
    parent_node = edges[,1],
    child_node = edges[,2],
    parent_count = numeric(n_edges),
    child_count = numeric(n_edges),
    change_ratio = numeric(n_edges),
    abs_change = numeric(n_edges),
    rel_change = numeric(n_edges),
    branch_length = tree$edge.length,
    event_type = character(n_edges),
    event_confidence = numeric(n_edges),
    event_magnitude = numeric(n_edges),
    p_value = numeric(n_edges),
    significant = logical(n_edges),
    stringsAsFactors = FALSE
  )
  
  # Calculate alpha critical value
  alpha <- 1 - confidence_level
  
  # Fill event data
  for(i in 1:n_edges) {
    parent_id <- edges[i, 1]
    child_id <- edges[i, 2]
    
    # Get node states
    parent_count <- all_nodes$state[all_nodes$node_id == parent_id]
    child_count <- all_nodes$state[all_nodes$node_id == child_id]
    
    events$parent_count[i] <- parent_count
    events$child_count[i] <- child_count
    
    # Calculate ratio and changes
    events$abs_change[i] <- child_count - parent_count
    
    if(parent_count > 0) {
      events$change_ratio[i] <- child_count / parent_count
      events$rel_change[i] <- (child_count - parent_count) / parent_count
    } else {
      events$change_ratio[i] <- NA
      events$rel_change[i] <- NA
    }
    
    # Determine if change is statistically significant
    # For internal nodes, use confidence intervals
    # For tips, assume observed values are fixed
    if(child_id <= length(tip_states)) {
      # Tip node - compare to parent CI
      parent_ci_lower <- all_nodes$ci_lower[all_nodes$node_id == parent_id]
      parent_ci_upper <- all_nodes$ci_upper[all_nodes$node_id == parent_id]
      
      # Check if tip value falls outside parent CI
      events$significant[i] <- child_count < parent_ci_lower || child_count > parent_ci_upper
      
      # Approximate p-value
      if(!is.na(parent_ci_lower) && !is.na(parent_ci_upper)) {
        parent_ci_width <- parent_ci_upper - parent_ci_lower
        parent_sd <- parent_ci_width / (2 * qnorm(confidence_level + (1-confidence_level)/2))
        
        # Two-tailed p-value
        events$p_value[i] <- 2 * pnorm(abs(child_count - parent_count) / parent_sd, lower.tail = FALSE)
      } else {
        events$p_value[i] <- NA
      }
    } else {
      # Internal node - compare CIs
      child_ci_lower <- all_nodes$ci_lower[all_nodes$node_id == child_id]
      child_ci_upper <- all_nodes$ci_upper[all_nodes$node_id == child_id]
      parent_ci_lower <- all_nodes$ci_lower[all_nodes$node_id == parent_id]
      parent_ci_upper <- all_nodes$ci_upper[all_nodes$node_id == parent_id]
      
      # Check if CIs don't overlap
      events$significant[i] <- (child_ci_lower > parent_ci_upper) || (child_ci_upper < parent_ci_lower)
      
      # Approximate p-value from overlap of CIs
      if(!is.na(child_ci_lower) && !is.na(child_ci_upper) && 
         !is.na(parent_ci_lower) && !is.na(parent_ci_upper)) {
        
        child_ci_width <- child_ci_upper - child_ci_lower
        parent_ci_width <- parent_ci_upper - parent_ci_lower
        
        child_sd <- child_ci_width / (2 * qnorm(confidence_level + (1-confidence_level)/2))
        parent_sd <- parent_ci_width / (2 * qnorm(confidence_level + (1-confidence_level)/2))
        
        # Combined standard error
        combined_se <- sqrt(child_sd^2 + parent_sd^2)
        
        # Two-tailed p-value
        events$p_value[i] <- 2 * pnorm(abs(child_count - parent_count) / combined_se, lower.tail = FALSE)
      } else {
        events$p_value[i] <- NA
      }
    }
    
    # Determine event type and confidence
    if(is.na(events$change_ratio[i])) {
      events$event_type[i] <- "unknown"
      events$event_confidence[i] <- 0
      events$event_magnitude[i] <- 0
    } else if(events$change_ratio[i] <= fusion_threshold && events$significant[i]) {
      # Fusion event
      events$event_type[i] <- "fusion"
      events$event_magnitude[i] <- parent_count - child_count
      
      # Calculate confidence based on both ratio and statistical significance
      ratio_conf <- min(1.0, (1 - events$change_ratio[i]) / (1 - fusion_threshold))
      stat_conf <- 1 - events$p_value[i]
      events$event_confidence[i] <- sqrt(ratio_conf * stat_conf)  # Geometric mean
      
    } else if(events$change_ratio[i] >= wgd_threshold && events$significant[i]) {
      # WGD event
      events$event_type[i] <- "wgd"
      events$event_magnitude[i] <- child_count - parent_count
      
      # Calculate confidence
      nearest_mult <- round(events$change_ratio[i])
      mult_conf <- 1 - abs(events$change_ratio[i] - nearest_mult) / 0.5
      stat_conf <- 1 - events$p_value[i]
      events$event_confidence[i] <- sqrt(mult_conf * stat_conf)  # Geometric mean
      
    } else if(events$change_ratio[i] >= fission_threshold && events$significant[i]) {
      # Fission event
      events$event_type[i] <- "fission"
      events$event_magnitude[i] <- child_count - parent_count
      
      # Calculate confidence
      ratio_conf <- min(1.0, (events$change_ratio[i] - fission_threshold) / fission_threshold)
      stat_conf <- 1 - events$p_value[i]
      events$event_confidence[i] <- sqrt(ratio_conf * stat_conf)  # Geometric mean
      
    } else {
      # No significant event
      events$event_type[i] <- "none"
      events$event_confidence[i] <- 1 - events$p_value[i]  # Higher for non-significant changes
      events$event_magnitude[i] <- abs(child_count - parent_count)
    }
  }
  
  # Add species names (for tip nodes)
  events$species <- NA
  tip_nodes <- events$child_node[events$child_node <= length(tree$tip.label)]
  events$species[events$child_node %in% tip_nodes] <- 
    tree$tip.label[events$child_node[events$child_node %in% tip_nodes]]
  
  # Summarize events
  event_counts <- table(events$event_type)
  message("Detected chromosome events:")
  for(event_type in names(event_counts)) {
    if(event_type != "none") {
      message(sprintf("  %s: %d", event_type, event_counts[event_type]))
    }
  }
  
  # Create result
  result <- list(
    events = events,
    tree = tree,
    ml_result = ml_result,
    parameters = list(
      fusion_threshold = fusion_threshold,
      fission_threshold = fission_threshold,
      wgd_threshold = wgd_threshold,
      confidence_level = confidence_level
    ),
    summary = list(
      event_counts = event_counts,
      avg_confidence = tapply(events$event_confidence, events$event_type, mean, na.rm = TRUE),
      avg_magnitude = tapply(events$event_magnitude, events$event_type, mean, na.rm = TRUE),
      significant_events = sum(events$significant & events$event_type != "none", na.rm = TRUE)
    )
  )
  
  return(result)
}
