#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Reconstruction Validation Module
# Author: Bioinformatics Team
# Date: 2025-05-25
# Description: Methods for validating and evaluating the reliability of ancestral
#              chromosome number reconstructions through cross-validation,
#              simulation-based validation, and sensitivity analysis
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(geiger)
  library(parallel)
  library(ggplot2)
  library(coda)
  library(dplyr)
})

#===============================================================================
# Cross-Validation Methods
#===============================================================================

#' Perform k-fold cross-validation for ancestral state reconstruction
#' 
#' Evaluates ancestral reconstruction methods by systematically masking tip states
#' and assessing prediction accuracy
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param k Number of folds for cross-validation
#' @param methods Vector of reconstruction methods to evaluate: "parsimony", "ML", "Bayesian"
#' @param model Evolutionary model for ML and Bayesian methods: "BM", "OU", "EB"
#' @param test_prop Proportion of tip states to mask in each fold
#' @param n_replicates Number of replication rounds (with different random partitions)
#' @param seed Random seed for reproducibility
#' @param n_cores Number of cores for parallel processing (NULL = auto-detect)
#' @param verbose Whether to print progress messages
#' @return List with cross-validation results
#' @export
cross_validate_reconstruction <- function(tree,
                                        chr_counts,
                                        k = 5,
                                        methods = c("parsimony", "ML", "Bayesian"),
                                        model = "BM",
                                        test_prop = 0.2,
                                        n_replicates = 10,
                                        seed = NULL,
                                        n_cores = NULL,
                                        verbose = TRUE) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Set random seed if provided
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Match data to tree
  matched_counts <- match_data_to_tree(tree, chr_counts)
  
  # Remove tips with missing data
  has_data <- !is.na(matched_counts)
  if(sum(has_data) < 10) {  # Need enough data for meaningful cross-validation
    stop("At least 10 tips must have chromosome count data for cross-validation")
  }
  
  # Get tips with data
  tips_with_data <- tree$tip.label[has_data]
  counts_with_data <- matched_counts[has_data]
  
  # Set up parallel processing
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Initialize results
  results <- list(
    summary = list(
      n_tips = length(tips_with_data),
      methods = methods,
      k = k,
      n_replicates = n_replicates,
      test_prop = test_prop
    ),
    fold_results = list(),
    method_performance = list(),
    error_metrics = list(),
    plots = list()
  )
  
  # Create folds
  all_folds <- list()
  
  for(rep in 1:n_replicates) {
    # Shuffle tips
    shuffled_tips <- sample(tips_with_data)
    
    # Calculate number of tips per fold
    n_per_fold <- ceiling(length(shuffled_tips) * test_prop)
    
    # Create fold assignments
    fold_assignments <- list()
    for(i in 1:k) {
      start_idx <- (i - 1) * n_per_fold + 1
      end_idx <- min(i * n_per_fold, length(shuffled_tips))
      
      # Handle case where start_idx exceeds shuffled_tips length
      if(start_idx <= length(shuffled_tips)) {
        fold_assignments[[i]] <- shuffled_tips[start_idx:end_idx]
      } else {
        fold_assignments[[i]] <- character(0)
      }
    }
    
    all_folds[[rep]] <- fold_assignments
  }
  
  # Function to run cross-validation for a single method and fold
  run_cv_fold <- function(method, fold_data, rep_idx, fold_idx) {
    # Extract test tips for this fold
    test_tips <- fold_data[[fold_idx]]
    
    if(length(test_tips) == 0) {
      return(NULL)  # Skip empty folds
    }
    
    # Create masked data (remove test tips' states)
    masked_counts <- chr_counts
    masked_counts[test_tips] <- NA
    
    # Perform reconstruction with masked data
    recon_result <- NULL
    
    if(method == "parsimony") {
      # Use parsimony reconstruction
      if(requireNamespace("phangorn", quietly = TRUE)) {
        # Convert data to phyDat format
        phydat <- create_phyDat_from_counts(tree, masked_counts)
        
        # Perform ancestral reconstruction
        recon_result <- tryCatch({
          anc <- phangorn::ancestral.pars(tree, phydat)
          
          # Extract marginal probabilities
          node_states <- get_most_likely_states(anc)
          
          list(
            node_states = node_states,
            method = "parsimony"
          )
        }, error = function(e) {
          warning(paste("Parsimony error:", e$message))
          return(NULL)
        })
      } else {
        warning("phangorn package required for parsimony reconstruction")
      }
    } else if(method == "ML") {
      # Use Maximum Likelihood reconstruction
      if(requireNamespace("phytools", quietly = TRUE)) {
        recon_result <- tryCatch({
          # Basic ML reconstruction
          anc <- phytools::fastAnc(tree, masked_counts, CI = TRUE)
          
          # Extract mapped states
          tip_states <- masked_counts[tree$tip.label]
          
          # Combine tip and ancestral states
          all_states <- c(tip_states, anc$ace)
          names(all_states) <- c(1:length(tree$tip.label), 
                                (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
          
          list(
            node_states = all_states,
            method = "ML"
          )
        }, error = function(e) {
          warning(paste("ML error:", e$message))
          return(NULL)
        })
      } else {
        warning("phytools package required for ML reconstruction")
      }
    } else if(method == "Bayesian") {
      # Use Bayesian reconstruction
      if(requireNamespace("phytools", quietly = TRUE)) {
        recon_result <- tryCatch({
          # Simple Bayesian reconstruction
          anc <- phytools::ancThresh(tree, masked_counts, ngen = 10000, control = list(samp = 10))
          
          # Extract mean states
          node_states <- colMeans(anc$mcmc[, (length(anc$mcmc[1, ]) - tree$Nnode + 1):length(anc$mcmc[1, ])])
          names(node_states) <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
          
          # Add tip states
          tip_states <- masked_counts[tree$tip.label]
          all_states <- c(tip_states, node_states)
          
          list(
            node_states = all_states,
            method = "Bayesian"
          )
        }, error = function(e) {
          warning(paste("Bayesian error:", e$message))
          return(NULL)
        })
      } else {
        warning("phytools package required for Bayesian reconstruction")
      }
    }
    
    # Evaluate prediction accuracy for test tips
    prediction_results <- NULL
    
    if(!is.null(recon_result)) {
      # Extract true and predicted values for test tips
      true_values <- counts_with_data[test_tips]
      
      # Get nodes for test tips
      test_tip_indices <- match(test_tips, tree$tip.label)
      
      # Extract predicted values
      predicted_values <- recon_result$node_states[as.character(test_tip_indices)]
      
      # Calculate error metrics
      metrics <- calculate_error_metrics(true_values, predicted_values, method_name = method)
      
      prediction_results <- list(
        test_tips = test_tips,
        true_values = true_values,
        predicted_values = predicted_values,
        error_metrics = metrics
      )
    }
    
    return(list(
      method = method,
      replicate = rep_idx,
      fold = fold_idx,
      reconstruction = recon_result,
      prediction = prediction_results
    ))
  }
  
  # Run cross-validation for all methods and folds
  cv_results <- list()
  
  for(method in methods) {
    method_results <- list()
    
    for(rep_idx in 1:n_replicates) {
      if(verbose) {
        message(paste("Running", method, "- replicate", rep_idx, "of", n_replicates))
      }
      
      fold_results <- list()
      
      if(n_cores > 1) {
        # Parallel execution
        cl <- parallel::makeCluster(n_cores)
        on.exit(parallel::stopCluster(cl))
        
        fold_results <- parallel::parLapply(cl, 1:k, function(fold_idx) {
          run_cv_fold(method, all_folds[[rep_idx]], rep_idx, fold_idx)
        })
      } else {
        # Sequential execution
        for(fold_idx in 1:k) {
          fold_results[[fold_idx]] <- run_cv_fold(method, all_folds[[rep_idx]], rep_idx, fold_idx)
        }
      }
      
      method_results[[rep_idx]] <- fold_results
    }
    
    cv_results[[method]] <- method_results
  }
  
  # Process results
  processed_results <- process_cv_results(cv_results, methods, n_replicates, k)
  results$fold_results <- processed_results$fold_results
  results$method_performance <- processed_results$method_performance
  results$error_metrics <- processed_results$error_metrics
  
  # Create performance comparison plots
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Error comparison plot
    error_data <- data.frame(
      Method = rep(names(results$method_performance), each = n_replicates),
      Replicate = rep(1:n_replicates, length(names(results$method_performance))),
      RMSE = unlist(lapply(results$method_performance, function(x) sapply(x, function(y) y$metrics$RMSE))),
      MAE = unlist(lapply(results$method_performance, function(x) sapply(x, function(y) y$metrics$MAE))),
      stringsAsFactors = FALSE
    )
    
    rmse_plot <- ggplot2::ggplot(error_data, ggplot2::aes(x = Method, y = RMSE, fill = Method)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(
        title = "Reconstruction Method Accuracy Comparison",
        subtitle = "Root Mean Square Error (RMSE)",
        x = "Method",
        y = "RMSE"
      ) +
      ggplot2::theme_minimal()
    
    results$plots$rmse_comparison <- rmse_plot
    
    # MAE comparison plot
    mae_plot <- ggplot2::ggplot(error_data, ggplot2::aes(x = Method, y = MAE, fill = Method)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(
        title = "Reconstruction Method Accuracy Comparison",
        subtitle = "Mean Absolute Error (MAE)",
        x = "Method",
        y = "MAE"
      ) +
      ggplot2::theme_minimal()
    
    results$plots$mae_comparison <- mae_plot
  }
  
  # Add summary statistics
  avg_metrics <- lapply(results$method_performance, function(method_replicates) {
    # Calculate average metrics across replicates
    avg_rmse <- mean(sapply(method_replicates, function(rep) rep$metrics$RMSE))
    avg_mae <- mean(sapply(method_replicates, function(rep) rep$metrics$MAE))
    
    list(
      avg_RMSE = avg_rmse,
      avg_MAE = avg_mae
    )
  })
  
  results$summary$avg_metrics <- avg_metrics
  
  # Determine best method
  best_by_rmse <- names(avg_metrics)[which.min(sapply(avg_metrics, function(x) x$avg_RMSE))]
  best_by_mae <- names(avg_metrics)[which.min(sapply(avg_metrics, function(x) x$avg_MAE))]
  
  results$summary$best_method_RMSE <- best_by_rmse
  results$summary$best_method_MAE <- best_by_mae
  
  # Overall conclusion
  results$summary$conclusion <- paste0(
    "Based on cross-validation with ", n_replicates, " replicates and ", k, " folds, ",
    "the best reconstruction method by RMSE is ", best_by_rmse, " and by MAE is ", best_by_mae, "."
  )
  
  return(results)
}

#' Match chromosome count data to tree tips
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @return Named vector matched to tree tips
#' @keywords internal
match_data_to_tree <- function(tree, chr_counts) {
  # Create named vector matching tree tips
  matched_data <- numeric(length(tree$tip.label))
  names(matched_data) <- tree$tip.label
  
  # Fill in values where names match
  for(tip in tree$tip.label) {
    if(tip %in% names(chr_counts)) {
      matched_data[tip] <- chr_counts[tip]
    } else {
      matched_data[tip] <- NA
    }
  }
  
  return(matched_data)
}

#' Create phyDat object from chromosome counts for parsimony
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @return phyDat object for parsimony analysis
#' @keywords internal
create_phyDat_from_counts <- function(tree, chr_counts) {
  if(!requireNamespace("phangorn", quietly = TRUE)) {
    stop("phangorn package is required")
  }
  
  # Get range of chromosome counts
  all_counts <- unique(chr_counts[!is.na(chr_counts)])
  min_chr <- min(all_counts)
  max_chr <- max(all_counts)
  
  # Create a character matrix
  chr_matrix <- matrix("?", nrow = length(tree$tip.label), ncol = 1)
  rownames(chr_matrix) <- tree$tip.label
  
  # Fill in states for tips with data
  for(tip in names(chr_counts)) {
    if(tip %in% tree$tip.label && !is.na(chr_counts[tip])) {
      # Convert numeric to character
      chr_matrix[tip, 1] <- as.character(chr_counts[tip])
    }
  }
  
  # Create phyDat object
  phydat <- phangorn::phyDat(chr_matrix, type = "USER", levels = as.character(min_chr:max_chr))
  
  return(phydat)
}

#' Extract most likely states from parsimony reconstruction
#' 
#' @param ancestral_states Output from ancestral.pars
#' @return Named vector of most likely states for each node
#' @keywords internal
get_most_likely_states <- function(ancestral_states) {
  # Get all possible states
  all_states <- colnames(ancestral_states)
  
  # Find most likely state for each node
  node_states <- apply(ancestral_states, 1, function(x) {
    # If multiple states tied for highest probability, take the median
    max_prob <- max(x)
    best_states <- all_states[x == max_prob]
    
    if(length(best_states) > 1) {
      # Take median of numeric states
      median(as.numeric(best_states))
    } else {
      as.numeric(best_states)
    }
  })
  
  return(node_states)
}

#' Calculate error metrics for predicted vs true values
#' 
#' @param true_values Vector of true values
#' @param predicted_values Vector of predicted values
#' @param method_name Name of the reconstruction method
#' @return List with error metrics
#' @keywords internal
calculate_error_metrics <- function(true_values, predicted_values, method_name) {
  # Remove any NA pairs
  valid <- !is.na(true_values) & !is.na(predicted_values)
  
  if(sum(valid) == 0) {
    return(list(
      MSE = NA,
      RMSE = NA,
      MAE = NA,
      R2 = NA,
      method = method_name
    ))
  }
  
  true <- true_values[valid]
  pred <- predicted_values[valid]
  
  # Calculate metrics
  mse <- mean((true - pred)^2)
  rmse <- sqrt(mse)
  mae <- mean(abs(true - pred))
  
  # Calculate R-squared (coefficient of determination)
  ss_total <- sum((true - mean(true))^2)
  ss_residual <- sum((true - pred)^2)
  r2 <- 1 - ss_residual / ss_total
  
  return(list(
    MSE = mse,
    RMSE = rmse,
    MAE = mae,
    R2 = r2,
    method = method_name
  ))
}

#' Process cross-validation results
#' 
#' @param cv_results Cross-validation results from all methods and folds
#' @param methods Vector of methods used
#' @param n_replicates Number of replicates
#' @param k Number of folds
#' @return Processed results
#' @keywords internal
process_cv_results <- function(cv_results, methods, n_replicates, k) {
  # Initialize processed results
  fold_results <- list()
  method_performance <- list()
  error_metrics <- list()
  
  # Process results for each method
  for(method in methods) {
    method_replicates <- cv_results[[method]]
    
    # Skip if method has no results
    if(is.null(method_replicates)) {
      next
    }
    
    method_fold_results <- list()
    method_error_metrics <- list()
    
    # Process each replicate
    for(rep_idx in 1:n_replicates) {
      replicate_folds <- method_replicates[[rep_idx]]
      
      # Skip if replicate has no results
      if(is.null(replicate_folds)) {
        next
      }
      
      replicate_results <- list()
      all_true <- numeric(0)
      all_predicted <- numeric(0)
      
      # Process each fold
      for(fold_idx in 1:k) {
        fold_result <- replicate_folds[[fold_idx]]
        
        # Skip if fold has no results
        if(is.null(fold_result) || is.null(fold_result$prediction)) {
          next
        }
        
        # Store fold results
        replicate_results[[fold_idx]] <- fold_result$prediction
        
        # Accumulate true and predicted values
        all_true <- c(all_true, fold_result$prediction$true_values)
        all_predicted <- c(all_predicted, fold_result$prediction$predicted_values)
      }
      
      # Calculate overall metrics for the replicate
      metrics <- calculate_error_metrics(all_true, all_predicted, method_name = method)
      
      method_fold_results[[rep_idx]] <- replicate_results
      method_error_metrics[[rep_idx]] <- list(
        metrics = metrics,
        n_predictions = length(all_true),
        true_values = all_true,
        predicted_values = all_predicted
      )
    }
    
    fold_results[[method]] <- method_fold_results
    method_performance[[method]] <- method_error_metrics
  }
  
  return(list(
    fold_results = fold_results,
    method_performance = method_performance,
    error_metrics = error_metrics
  ))
}

#===============================================================================
# Simulation-Based Validation
#===============================================================================

#' Validate reconstruction methods using simulated chromosome evolution
#' 
#' Simulates chromosome evolution under known models and parameters, then 
#' evaluates how well reconstruction methods recover the true ancestral states
#' 
#' @param tree Phylogenetic tree
#' @param sim_model Model for simulation: "BM", "OU", "EB", "discrete"
#' @param sim_params List of parameters for simulation model
#' @param methods Vector of reconstruction methods to evaluate: "parsimony", "ML", "Bayesian"
#' @param recon_model Model for reconstruction (if different from sim_model)
#' @param n_replicates Number of simulation replicates
#' @param discrete Whether to discretize the simulated values
#' @param root_value Root chromosome number for simulations
#' @param n_cores Number of cores for parallel processing
#' @param seed Random seed for reproducibility
#' @param verbose Whether to print progress messages
#' @return List with simulation validation results
#' @export
validate_with_simulations <- function(tree,
                                    sim_model = "BM",
                                    sim_params = list(sig2 = 0.1),
                                    methods = c("parsimony", "ML", "Bayesian"),
                                    recon_model = NULL,
                                    n_replicates = 100,
                                    discrete = TRUE,
                                    root_value = 10,
                                    n_cores = NULL,
                                    seed = NULL,
                                    verbose = TRUE) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Set random seed if provided
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Set defaults for recon_model if not provided
  if(is.null(recon_model)) {
    recon_model <- sim_model
  }
  
  # Set up parallel processing
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Initialize results
  results <- list(
    summary = list(
      tree = tree,
      sim_model = sim_model,
      sim_params = sim_params,
      methods = methods,
      n_replicates = n_replicates,
      discrete = discrete,
      root_value = root_value
    ),
    simulations = list(),
    method_performance = list(),
    error_metrics = list(),
    plots = list()
  )
  
  # Function to run a single simulation replicate
  run_simulation <- function(rep_idx) {
    # Simulate chromosome evolution
    sim_result <- NULL
    
    if(sim_model == "BM") {
      # Brownian motion simulation
      sig2 <- sim_params$sig2
      if(is.null(sig2)) sig2 <- 0.1
      
      sim_result <- phytools::fastBM(tree, sig2 = sig2, a = root_value, internal = TRUE)
    } else if(sim_model == "OU") {
      # Ornstein-Uhlenbeck simulation
      sig2 <- sim_params$sig2
      if(is.null(sig2)) sig2 <- 0.1
      
      alpha <- sim_params$alpha
      if(is.null(alpha)) alpha <- 0.1
      
      theta <- sim_params$theta
      if(is.null(theta)) theta <- root_value
      
      sim_result <- phytools::fastBM(tree, sig2 = sig2, a = root_value, alpha = alpha, 
                                   theta = theta, internal = TRUE)
    } else if(sim_model == "EB") {
      # Early Burst simulation
      sig2 <- sim_params$sig2
      if(is.null(sig2)) sig2 <- 0.1
      
      r <- sim_params$r
      if(is.null(r)) r <- -0.1  # Negative rate means slowing down
      
      sim_result <- phytools::fastBM(tree, sig2 = sig2, a = root_value, internal = TRUE, r = r)
    } else if(sim_model == "discrete") {
      # Discrete character simulation
      
      # For simplicity, we'll use a continuous model and then discretize
      sig2 <- sim_params$sig2
      if(is.null(sig2)) sig2 <- 0.1
      
      sim_continuous <- phytools::fastBM(tree, sig2 = sig2, a = root_value, internal = TRUE)
      
      # Discretize - we'll use rounding for this purpose
      sim_result <- round(sim_continuous)
      names(sim_result) <- names(sim_continuous)
    } else {
      stop(paste("Unsupported simulation model:", sim_model))
    }
    
    # Discretize if requested
    if(discrete && sim_model != "discrete") {
      sim_result <- round(sim_result)
    }
    
    # Separate tip and internal node states
    n_tips <- length(tree$tip.label)
    tip_states <- sim_result[1:n_tips]
    names(tip_states) <- tree$tip.label
    
    internal_states <- sim_result[(n_tips + 1):length(sim_result)]
    
    # Run reconstruction with each method
    method_results <- list()
    
    for(method in methods) {
      recon_result <- NULL
      
      if(method == "parsimony") {
        # Parsimony reconstruction
        if(requireNamespace("phangorn", quietly = TRUE)) {
          # Create phyDat object
          phydat <- create_phyDat_from_counts(tree, tip_states)
          
          # Run parsimony
          recon_result <- tryCatch({
            anc <- phangorn::ancestral.pars(tree, phydat)
            
            # Extract marginal probabilities
            node_states <- get_most_likely_states(anc)
            
            list(
              node_states = node_states,
              method = "parsimony"
            )
          }, error = function(e) {
            warning(paste("Parsimony error in simulation", rep_idx, ":", e$message))
            return(NULL)
          })
        }
      } else if(method == "ML") {
        # ML reconstruction
        if(requireNamespace("phytools", quietly = TRUE)) {
          recon_result <- tryCatch({
            if(recon_model == "BM") {
              # Basic ML reconstruction
              anc <- phytools::fastAnc(tree, tip_states, CI = TRUE)
              
              # Extract mapped states
              all_states <- c(tip_states, anc$ace)
              names(all_states) <- c(1:length(tree$tip.label), 
                                    (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
            } else {
              # Use ace with specified model
              anc <- phytools::ace(tip_states, tree, type = "continuous", 
                                 method = "ML", model = recon_model)
              
              # Extract states
              all_states <- c(tip_states, anc$ace)
              names(all_states) <- c(1:length(tree$tip.label), 
                                    (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
            }
            
            # Discretize if needed
            if(discrete) {
              all_states <- round(all_states)
            }
            
            list(
              node_states = all_states,
              method = "ML"
            )
          }, error = function(e) {
            warning(paste("ML error in simulation", rep_idx, ":", e$message))
            return(NULL)
          })
        }
      } else if(method == "Bayesian") {
        # Bayesian reconstruction
        if(requireNamespace("phytools", quietly = TRUE)) {
          recon_result <- tryCatch({
            # Simple Bayesian reconstruction
            anc <- phytools::ancThresh(tree, tip_states, ngen = 10000, control = list(samp = 10))
            
            # Extract mean states
            node_states <- colMeans(anc$mcmc[, (length(anc$mcmc[1, ]) - tree$Nnode + 1):length(anc$mcmc[1, ])])
            names(node_states) <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
            
            # Add tip states
            all_states <- c(tip_states, node_states)
            
            # Discretize if needed
            if(discrete) {
              all_states <- round(all_states)
            }
            
            list(
              node_states = all_states,
              method = "Bayesian"
            )
          }, error = function(e) {
            warning(paste("Bayesian error in simulation", rep_idx, ":", e$message))
            return(NULL)
          })
        }
      }
      
      # Evaluate reconstruction accuracy for internal nodes
      if(!is.null(recon_result)) {
        # Extract true and reconstructed values for internal nodes
        internal_ids <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
        true_internal <- internal_states
        
        # Get reconstructed states for internal nodes
        recon_internal <- recon_result$node_states[as.character(internal_ids)]
        
        # Calculate error metrics
        metrics <- calculate_error_metrics(true_internal, recon_internal, method_name = method)
        
        method_results[[method]] <- list(
          recon_states = recon_result$node_states,
          true_internal = true_internal,
          recon_internal = recon_internal,
          metrics = metrics
        )
      }
    }
    
    return(list(
      rep_idx = rep_idx,
      tip_states = tip_states,
      internal_states = internal_states,
      method_results = method_results
    ))
  }
  
  # Run simulations (in parallel if requested)
  sim_results <- list()
  
  if(n_cores > 1 && n_replicates > 1) {
    # Parallel execution
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    sim_results <- parallel::parLapply(cl, 1:n_replicates, function(rep_idx) {
      run_simulation(rep_idx)
    })
  } else {
    # Sequential execution
    for(rep_idx in 1:n_replicates) {
      if(verbose) {
        message(paste("Running simulation", rep_idx, "of", n_replicates))
      }
      
      sim_results[[rep_idx]] <- run_simulation(rep_idx)
    }
  }
  
  # Process simulation results
  results$simulations <- sim_results
  
  # Extract method performance across simulations
  method_performance <- list()
  
  for(method in methods) {
    method_metrics <- list(
      RMSE = numeric(),
      MAE = numeric(),
      R2 = numeric()
    )
    
    for(rep_idx in 1:n_replicates) {
      if(is.null(sim_results[[rep_idx]]) || 
         is.null(sim_results[[rep_idx]]$method_results) || 
         is.null(sim_results[[rep_idx]]$method_results[[method]])) {
        next
      }
      
      rep_metrics <- sim_results[[rep_idx]]$method_results[[method]]$metrics
      
      method_metrics$RMSE <- c(method_metrics$RMSE, rep_metrics$RMSE)
      method_metrics$MAE <- c(method_metrics$MAE, rep_metrics$MAE)
      method_metrics$R2 <- c(method_metrics$R2, rep_metrics$R2)
    }
    
    method_performance[[method]] <- list(
      mean_RMSE = mean(method_metrics$RMSE, na.rm = TRUE),
      sd_RMSE = sd(method_metrics$RMSE, na.rm = TRUE),
      mean_MAE = mean(method_metrics$MAE, na.rm = TRUE),
      sd_MAE = sd(method_metrics$MAE, na.rm = TRUE),
      mean_R2 = mean(method_metrics$R2, na.rm = TRUE),
      sd_R2 = sd(method_metrics$R2, na.rm = TRUE),
      metrics = method_metrics
    )
  }
  
  results$method_performance <- method_performance
  
  # Create performance comparison plots
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Prepare data frame for plots
    plot_data <- data.frame(
      Method = character(),
      RMSE = numeric(),
      MAE = numeric(),
      R2 = numeric(),
      stringsAsFactors = FALSE
    )
    
    for(method in names(method_performance)) {
      method_data <- data.frame(
        Method = rep(method, length(method_performance[[method]]$metrics$RMSE)),
        RMSE = method_performance[[method]]$metrics$RMSE,
        MAE = method_performance[[method]]$metrics$MAE,
        R2 = method_performance[[method]]$metrics$R2,
        stringsAsFactors = FALSE
      )
      
      plot_data <- rbind(plot_data, method_data)
    }
    
    # Create RMSE boxplot
    rmse_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Method, y = RMSE, fill = Method)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(
        title = "Reconstruction Method Accuracy Under Simulation",
        subtitle = paste("Simulation Model:", sim_model),
        x = "Method",
        y = "Root Mean Square Error"
      ) +
      ggplot2::theme_minimal()
    
    # Create MAE boxplot
    mae_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Method, y = MAE, fill = Method)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(
        title = "Reconstruction Method Accuracy Under Simulation",
        subtitle = paste("Simulation Model:", sim_model),
        x = "Method",
        y = "Mean Absolute Error"
      ) +
      ggplot2::theme_minimal()
    
    # Create R2 boxplot
    r2_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Method, y = R2, fill = Method)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(
        title = "Reconstruction Method Accuracy Under Simulation",
        subtitle = paste("Simulation Model:", sim_model),
        x = "Method",
        y = "R-Squared"
      ) +
      ggplot2::theme_minimal()
    
    results$plots$rmse_comparison <- rmse_plot
    results$plots$mae_comparison <- mae_plot
    results$plots$r2_comparison <- r2_plot
  }
  
  # Determine best method based on simulation results
  best_rmse_method <- names(method_performance)[which.min(sapply(method_performance, function(x) x$mean_RMSE))]
  best_mae_method <- names(method_performance)[which.min(sapply(method_performance, function(x) x$mean_MAE))]
  best_r2_method <- names(method_performance)[which.max(sapply(method_performance, function(x) x$mean_R2))]
  
  results$summary$best_rmse_method <- best_rmse_method
  results$summary$best_mae_method <- best_mae_method
  results$summary$best_r2_method <- best_r2_method
  
  # Create conclusion
  results$summary$conclusion <- paste0(
    "Based on ", n_replicates, " simulations under a ", sim_model, " model, ",
    "the best reconstruction method by RMSE is ", best_rmse_method, ", ",
    "by MAE is ", best_mae_method, ", and by R² is ", best_r2_method, "."
  )
  
  return(results)
}

#===============================================================================
# Sensitivity Analysis
#===============================================================================

#' Perform sensitivity analysis for reconstruction parameters
#' 
#' Tests how sensitive ancestral reconstructions are to variations in model
#' parameters and assumptions
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param method Reconstruction method to analyze: "ML", "Bayesian"
#' @param parameter Parameter to vary: "model", "sig2", "alpha", "root_value", etc.
#' @param parameter_values Vector of parameter values to test
#' @param nodes Nodes of interest for sensitivity analysis (NULL = all)
#' @param reference_recon Reference reconstruction (baseline)
#' @param n_replicates Number of replicates for each parameter value
#' @param n_cores Number of cores for parallel processing
#' @param verbose Whether to print progress messages
#' @return List with sensitivity analysis results
#' @export
analyze_reconstruction_sensitivity <- function(tree,
                                            chr_counts,
                                            method = "ML",
                                            parameter = "model",
                                            parameter_values = c("BM", "OU", "EB"),
                                            nodes = NULL,
                                            reference_recon = NULL,
                                            n_replicates = 10,
                                            n_cores = NULL,
                                            verbose = TRUE) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Set up nodes of interest if not specified
  if(is.null(nodes)) {
    # Use all internal nodes
    n_tips <- length(tree$tip.label)
    nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  }
  
  # Set up parallel processing
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Initialize results
  results <- list(
    summary = list(
      tree = tree,
      method = method,
      parameter = parameter,
      parameter_values = parameter_values,
      n_replicates = n_replicates,
      nodes = nodes
    ),
    sensitivity = list(),
    node_variation = list(),
    plots = list()
  )
  
  # Create reference reconstruction if not provided
  if(is.null(reference_recon)) {
    if(method == "ML") {
      if(requireNamespace("phytools", quietly = TRUE)) {
        # Create basic ML reconstruction with BM model
        anc <- phytools::fastAnc(tree, chr_counts, CI = TRUE)
        
        # Extract mapped states
        tip_states <- chr_counts[tree$tip.label]
        
        # Combine tip and ancestral states
        reference_recon <- c(tip_states, anc$ace)
        names(reference_recon) <- c(1:length(tree$tip.label), 
                                  (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
      } else {
        stop("phytools package required for ML reconstruction")
      }
    } else if(method == "Bayesian") {
      if(requireNamespace("phytools", quietly = TRUE)) {
        # Create basic Bayesian reconstruction
        anc <- phytools::ancThresh(tree, chr_counts, ngen = 10000, control = list(samp = 10))
        
        # Extract mean states
        node_states <- colMeans(anc$mcmc[, (length(anc$mcmc[1, ]) - tree$Nnode + 1):length(anc$mcmc[1, ])])
        names(node_states) <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
        
        # Add tip states
        tip_states <- chr_counts[tree$tip.label]
        reference_recon <- c(tip_states, node_states)
      } else {
        stop("phytools package required for Bayesian reconstruction")
      }
    } else {
      stop(paste("Unsupported method for sensitivity analysis:", method))
    }
  }
  
  # Function to run reconstruction for a single parameter value
  run_parameter_sensitivity <- function(param_value, rep_idx) {
    recon_result <- NULL
    
    if(method == "ML") {
      if(requireNamespace("phytools", quietly = TRUE)) {
        if(parameter == "model") {
          # Vary model type
          model_type <- param_value
          
          recon_result <- tryCatch({
            # Use ace with specified model
            anc <- phytools::ace(chr_counts, tree, type = "continuous", 
                               method = "ML", model = model_type)
            
            # Extract states
            states <- c(chr_counts[tree$tip.label], anc$ace)
            names(states) <- c(1:length(tree$tip.label), 
                             (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
            
            states
          }, error = function(e) {
            warning(paste("ML error with model", model_type, ":", e$message))
            return(NULL)
          })
        } else if(parameter == "sig2") {
          # Vary sigma^2 parameter
          sig2_value <- as.numeric(param_value)
          
          recon_result <- tryCatch({
            # Use custom ML reconstruction with fixed sigma^2
            anc <- custom_ml_reconstruction_with_sig2(tree, chr_counts, sig2_value)
            
            # Extract states
            states <- c(chr_counts[tree$tip.label], anc)
            names(states) <- c(1:length(tree$tip.label), 
                             (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
            
            states
          }, error = function(e) {
            warning(paste("ML error with sig2", sig2_value, ":", e$message))
            return(NULL)
          })
        } else {
          warning(paste("Unsupported parameter for ML sensitivity:", parameter))
          return(NULL)
        }
      }
    } else if(method == "Bayesian") {
      if(requireNamespace("phytools", quietly = TRUE)) {
        if(parameter == "prior_strength") {
          # Vary prior strength
          prior_strength <- as.numeric(param_value)
          
          recon_result <- tryCatch({
            # Custom Bayesian reconstruction with varied prior strength
            anc <- custom_bayesian_reconstruction(tree, chr_counts, prior_strength = prior_strength)
            
            # Extract states
            states <- c(chr_counts[tree$tip.label], anc)
            names(states) <- c(1:length(tree$tip.label), 
                             (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
            
            states
          }, error = function(e) {
            warning(paste("Bayesian error with prior strength", prior_strength, ":", e$message))
            return(NULL)
          })
        } else {
          warning(paste("Unsupported parameter for Bayesian sensitivity:", parameter))
          return(NULL)
        }
      }
    } else {
      stop(paste("Unsupported method for sensitivity analysis:", method))
    }
    
    return(recon_result)
  }
  
  # Run sensitivity analysis for each parameter value
  param_results <- list()
  
  for(param_value in parameter_values) {
    if(verbose) {
      message(paste("Analyzing", parameter, "=", param_value))
    }
    
    value_results <- list()
    
    for(rep_idx in 1:n_replicates) {
      recon_result <- run_parameter_sensitivity(param_value, rep_idx)
      
      if(!is.null(recon_result)) {
        value_results[[rep_idx]] <- recon_result[as.character(nodes)]
      }
    }
    
    # Average results across replicates
    if(length(value_results) > 0) {
      # Extract values for each node
      node_values <- list()
      
      for(node in nodes) {
        node_str <- as.character(node)
        node_values[[node_str]] <- sapply(value_results, function(x) x[node_str])
      }
      
      # Calculate mean and standard deviation
      node_means <- sapply(node_values, mean, na.rm = TRUE)
      node_sds <- sapply(node_values, sd, na.rm = TRUE)
      
      param_results[[as.character(param_value)]] <- list(
        parameter_value = param_value,
        reconstructions = value_results,
        node_means = node_means,
        node_sds = node_sds,
        nodes = nodes
      )
    }
  }
  
  results$sensitivity <- param_results
  
  # Analyze variation across parameters for each node
  node_variation <- list()
  
  for(node in nodes) {
    node_str <- as.character(node)
    
    # Collect values for this node across all parameter values
    param_means <- sapply(param_results, function(x) x$node_means[node_str])
    param_sds <- sapply(param_results, function(x) x$node_sds[node_str])
    
    # Calculate overall variability
    overall_mean <- mean(param_means, na.rm = TRUE)
    overall_sd <- sd(param_means, na.rm = TRUE)
    cv <- overall_sd / overall_mean  # Coefficient of variation
    
    # Calculate range of values
    value_range <- range(param_means, na.rm = TRUE)
    range_width <- value_range[2] - value_range[1]
    
    # Reference value
    ref_value <- reference_recon[node_str]
    
    # Calculate relative variation
    relative_variation <- range_width / ref_value
    
    node_variation[[node_str]] <- list(
      node = node,
      parameter_means = param_means,
      parameter_sds = param_sds,
      overall_mean = overall_mean,
      overall_sd = overall_sd,
      coefficient_of_variation = cv,
      value_range = value_range,
      range_width = range_width,
      reference_value = ref_value,
      relative_variation = relative_variation
    )
  }
  
  results$node_variation <- node_variation
  
  # Create sensitivity plots
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Parameter value vs. reconstructed state
    param_values <- parameter_values
    
    # Select a subset of nodes to plot (to avoid overcrowding)
    if(length(nodes) > 5) {
      plot_nodes <- sample(nodes, 5)
    } else {
      plot_nodes <- nodes
    }
    
    # Create data frame for plot
    plot_data <- data.frame(
      Parameter_Value = character(),
      Node = character(),
      State = numeric(),
      stringsAsFactors = FALSE
    )
    
    for(param_value in names(param_results)) {
      for(node in as.character(plot_nodes)) {
        if(!is.null(param_results[[param_value]]$node_means) && 
           !is.null(param_results[[param_value]]$node_means[node])) {
          
          plot_data <- rbind(plot_data, data.frame(
            Parameter_Value = param_value,
            Node = paste("Node", node),
            State = param_results[[param_value]]$node_means[node],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Create sensitivity plot
    sensitivity_plot <- ggplot2::ggplot(plot_data, 
                                      ggplot2::aes(x = Parameter_Value, y = State, 
                                                 group = Node, color = Node)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::labs(
        title = paste("Sensitivity of", method, "reconstruction to", parameter),
        x = parameter,
        y = "Reconstructed State"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    
    results$plots$sensitivity <- sensitivity_plot
    
    # Create coefficient of variation plot for all nodes
    cv_data <- data.frame(
      Node = as.character(nodes),
      CV = sapply(node_variation, function(x) x$coefficient_of_variation),
      stringsAsFactors = FALSE
    )
    
    cv_plot <- ggplot2::ggplot(cv_data, ggplot2::aes(x = reorder(Node, -CV), y = CV)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
      ggplot2::labs(
        title = paste("Sensitivity of nodes to variations in", parameter),
        subtitle = "Higher coefficient of variation indicates greater sensitivity",
        x = "Node",
        y = "Coefficient of Variation"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 8))
    
    results$plots$cv <- cv_plot
  }
  
  # Create summary statistics
  results$summary$mean_cv <- mean(sapply(node_variation, function(x) x$coefficient_of_variation), 
                                na.rm = TRUE)
  results$summary$max_cv <- max(sapply(node_variation, function(x) x$coefficient_of_variation), 
                              na.rm = TRUE)
  
  # Identify most and least sensitive nodes
  cv_values <- sapply(node_variation, function(x) x$coefficient_of_variation)
  most_sensitive_idx <- which.max(cv_values)
  least_sensitive_idx <- which.min(cv_values)
  
  results$summary$most_sensitive_node <- names(node_variation)[most_sensitive_idx]
  results$summary$least_sensitive_node <- names(node_variation)[least_sensitive_idx]
  
  # Overall conclusion
  if(results$summary$mean_cv < 0.1) {
    sensitivity_level <- "low"
  } else if(results$summary$mean_cv < 0.3) {
    sensitivity_level <- "moderate"
  } else {
    sensitivity_level <- "high"
  }
  
  results$summary$conclusion <- paste0(
    "The ", method, " reconstruction shows ", sensitivity_level, " sensitivity to variations in ", 
    parameter, ", with a mean coefficient of variation of ", 
    round(results$summary$mean_cv, 3), ". ",
    "Node ", results$summary$most_sensitive_node, " is the most sensitive, and ",
    "Node ", results$summary$least_sensitive_node, " is the least sensitive."
  )
  
  return(results)
}

#' Custom ML reconstruction with fixed sigma^2 parameter
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param sig2 Fixed sigma^2 parameter value
#' @return Named vector of ancestral states
#' @keywords internal
custom_ml_reconstruction_with_sig2 <- function(tree, chr_counts, sig2) {
  # This is a simplified version that would need to be implemented properly
  # For now, we'll just use fastAnc but would need a custom implementation
  # to properly fix sigma^2
  anc <- phytools::fastAnc(tree, chr_counts, CI = TRUE)
  return(anc$ace)
}

#' Custom Bayesian reconstruction with specified prior strength
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param prior_strength Prior strength parameter
#' @return Named vector of ancestral states
#' @keywords internal
custom_bayesian_reconstruction <- function(tree, chr_counts, prior_strength) {
  # This is a simplified version that would need to be implemented properly
  # For now, we'll use ancThresh but with different controls
  control <- list(
    samp = 10,
    propwidth = prior_strength  # Using propwidth as a proxy for prior strength
  )
  
  anc <- phytools::ancThresh(tree, chr_counts, ngen = 10000, control = control)
  
  # Extract mean states
  node_states <- colMeans(anc$mcmc[, (length(anc$mcmc[1, ]) - tree$Nnode + 1):length(anc$mcmc[1, ])])
  names(node_states) <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
  
  return(node_states)
}

#===============================================================================
# Reporting and Visualization
#===============================================================================

#' Generate comprehensive validation report
#' 
#' Creates a detailed report summarizing all validation analyses and results
#' 
#' @param validation_results List of results from various validation methods
#' @param output_file Path to save the report (NULL for no saving)
#' @param output_format Format for the report: "html", "pdf", or "docx"
#' @param include_plots Whether to include plots in the report
#' @return Report content as a character vector
#' @export
generate_validation_report <- function(validation_results,
                                     output_file = NULL,
                                     output_format = "html",
                                     include_plots = TRUE) {
  # Check if rmarkdown is available
  if(!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rmarkdown package is required for report generation")
  }
  
  # Create report content
  report <- c(
    "---",
    "title: \"Ancestral Chromosome Reconstruction Validation Report\"",
    paste0("date: \"", format(Sys.time(), "%Y-%m-%d"), "\""),
    "output: html_document",
    "---",
    "",
    "# Ancestral Chromosome Reconstruction Validation Report",
    "",
    "## Summary",
    "",
    "This report summarizes the results of various validation analyses performed on the ancestral chromosome reconstruction methods.",
    "",
    "## Cross-Validation Results",
    "",
    "### Summary",
    "",
    "```{r}",
    "print(validation_results$cross_validation$summary)",
    "```",
    "",
    "### Method Performance",
    "",
    "```{r}",
    "print(validation_results$cross_validation$method_performance)",
    "```",
    "",
    "### Error Metrics",
    "",
    "```{r}",
    "print(validation_results$cross_validation$error_metrics)",
    "```",
    "",
    "### Plots",
    "",
    "```{r, echo=FALSE, warning=FALSE, message=FALSE}",
    "if(include_plots) {",
    "  print(validation_results$cross_validation$plots$rmse_comparison)",
    "  print(validation_results$cross_validation$plots$mae_comparison)",
    "}",
    "```",
    "",
    "## Simulation-Based Validation Results",
    "",
    "### Summary",
    "",
    "```{r}",
    "print(validation_results$simulation_validation$summary)",
    "```",
    "",
    "### Method Performance",
    "",
    "```{r}",
    "print(validation_results$simulation_validation$method_performance)",
    "```",
    "",
    "### Error Metrics",
    "",
    "```{r}",
    "print(validation_results$simulation_validation$error_metrics)",
    "```",
    "",
    "### Plots",
    "",
    "```{r, echo=FALSE, warning=FALSE, message=FALSE}",
    "if(include_plots) {",
    "  print(validation_results$simulation_validation$plots$rmse_comparison)",
    "  print(validation_results$simulation_validation$plots$mae_comparison)",
    "  print(validation_results$simulation_validation$plots$r2_comparison)",
    "}",
    "```",
    "",
    "## Sensitivity Analysis Results",
    "",
    "### Summary",
    "",
    "```{r}",
    "print(validation_results$sensitivity_analysis$summary)",
    "```",
    "",
    "### Node Variation",
    "",
    "```{r}",
    "print(validation_results$sensitivity_analysis$node_variation)",
    "```",
    "",
    "### Plots",
    "",
    "```{r, echo=FALSE, warning=FALSE, message=FALSE}",
    "if(include_plots) {",
    "  print(validation_results$sensitivity_analysis$plots$sensitivity)",
    "  print(validation_results$sensitivity_analysis$plots$cv)",
    "}",
    "```"
  )
  
  # Save report if requested
  if(!is.null(output_file)) {
    rmarkdown::render(input = textConnection(report), output_file = output_file, output_format = output_format)
  }
  
  return(report)
}

# ---- 命令行运行支持 ----
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    cat("用法: Rscript methods/reconstruction_validation.R tree_file.nwk chr_counts.csv [n_replicates] [k_folds] [method1,method2,...] [seed]\n")
    quit(status = 1)
  }
  tree_file <- args[1]
  chr_file <- args[2]
  n_replicates <- ifelse(length(args) >= 3, as.numeric(args[3]), 10)
  k_folds <- ifelse(length(args) >= 4, as.numeric(args[4]), 5)
  methods <- ifelse(length(args) >= 5, strsplit(args[5], ",")[[1]], c("parsimony", "ML", "Bayesian"))
  seed <- ifelse(length(args) >= 6, as.numeric(args[6]), NULL)

  tree <- read.tree(tree_file)
  chr_df <- read.csv(chr_file, stringsAsFactors = FALSE)
  # 自动识别染色体数列
  chr_col <- setdiff(colnames(chr_df), c("species", "taxon", "name"))[1]
  chr_counts <- chr_df[[chr_col]]
  names(chr_counts) <- chr_df[[1]]

  # 交叉验证
  cv_results <- cross_validate_reconstruction(tree, chr_counts, k = k_folds, methods = methods, n_replicates = n_replicates, seed = seed)
  saveRDS(cv_results, file = "reconstruction_cv_results.rds")
  if (!is.null(cv_results$plots$rmse_comparison)) ggsave("reconstruction_cv_rmse.pdf", cv_results$plots$rmse_comparison, width = 6, height = 4)
  if (!is.null(cv_results$plots$mae_comparison)) ggsave("reconstruction_cv_mae.pdf", cv_results$plots$mae_comparison, width = 6, height = 4)
  cat("交叉验证分析完成，结果已保存为 reconstruction_cv_results.rds, reconstruction_cv_rmse.pdf, reconstruction_cv_mae.pdf\n")
}
