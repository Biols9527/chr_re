#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Ensemble Methods Module
# Author: Bioinformatics Team
# Date: 2025-03-26
# Description: Implements ensemble approaches for combining results from multiple
#              ancestral chromosome reconstruction methods
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(phytools)
})

#===============================================================================
# Ensemble Reconstruction Functions
#===============================================================================

#' Create ensemble reconstruction from multiple methods
#' 
#' Combines ancestral chromosome number reconstructions from multiple methods
#' 
#' @param reconstructions List of reconstruction results from different methods
#' @param weights Optional named vector of weights for different methods
#' @param method Ensemble method: "weighted_mean", "median", "ml", "bayesian_model_averaging"
#' @param quality_weighted Whether to weight by reconstruction quality metrics
#' @param tree Phylogenetic tree (if not included in reconstructions)
#' @return Ensemble reconstruction result
#' @export
create_ensemble_reconstruction <- function(reconstructions, 
                                         weights = NULL, 
                                         method = "weighted_mean",
                                         quality_weighted = TRUE,
                                         tree = NULL) {
  # Check input
  if(!is.list(reconstructions) || length(reconstructions) < 2) {
    stop("At least two reconstruction results are required for ensemble analysis")
  }
  
  method_names <- names(reconstructions)
  if(is.null(method_names)) {
    method_names <- paste0("method", 1:length(reconstructions))
    names(reconstructions) <- method_names
  }
  
  message(sprintf("Creating ensemble reconstruction using %s method from %d reconstruction methods: %s", 
                 method, length(reconstructions), paste(method_names, collapse = ", ")))
  
  # Extract tree from first reconstruction if not provided
  if(is.null(tree)) {
    tree <- reconstructions[[1]]$tree
    if(is.null(tree)) {
      stop("Tree not found in reconstructions and not provided separately")
    }
  }
  
  # Verify all reconstructions use same tree topology
  for(i in 1:length(reconstructions)) {
    if(!is.null(reconstructions[[i]]$tree)) {
      # Compare tree topologies
      if(!ape::all.equal.phylo(tree, reconstructions[[i]]$tree, use.edge.length = FALSE)) {
        warning("Tree topology in method '", method_names[i], "' differs from the reference tree")
      }
    }
  }
  
  # Get node IDs for all internal nodes
  node_ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
  
  # Determine weights for each method
  if(is.null(weights)) {
    if(quality_weighted) {
      # Derive weights from quality metrics
      weights <- calculate_quality_weights(reconstructions)
    } else {
      # Equal weights
      weights <- rep(1/length(reconstructions), length(reconstructions))
      names(weights) <- method_names
    }
  } else {
    # Normalize user-provided weights
    weights <- weights / sum(weights)
  }
  
  message("Using method weights:")
  for(i in 1:length(weights)) {
    message(sprintf("  %s: %.3f", names(weights)[i], weights[i]))
  }
  
  # Apply appropriate ensemble method
  if(method == "weighted_mean") {
    ensemble_result <- ensemble_weighted_mean(reconstructions, weights, node_ids, tree)
  } else if(method == "median") {
    ensemble_result <- ensemble_median(reconstructions, node_ids, tree)
  } else if(method == "ml") {
    ensemble_result <- ensemble_maximum_likelihood(reconstructions, weights, node_ids, tree)
  } else if(method == "bayesian_model_averaging") {
    ensemble_result <- ensemble_bayesian_model_averaging(reconstructions, weights, node_ids, tree)
  } else {
    stop("Unsupported ensemble method: ", method)
  }
  
  # Add information about source reconstructions
  ensemble_result$source_methods <- method_names
  ensemble_result$method_weights <- weights
  ensemble_result$tree <- tree
  
  # Calculate ensemble quality metrics
  ensemble_result$quality <- calculate_ensemble_quality(ensemble_result, reconstructions)
  
  return(ensemble_result)
}

#' Calculate quality-based weights for methods
#' 
#' Derives weights for different methods based on their quality metrics
#' 
#' @param reconstructions List of reconstruction results
#' @return Named vector of weights
#' @keywords internal
calculate_quality_weights <- function(reconstructions) {
  method_names <- names(reconstructions)
  n_methods <- length(reconstructions)
  
  # Initialize with equal weights as fallback
  weights <- rep(1/n_methods, n_methods)
  names(weights) <- method_names
  
  # Try to extract quality metrics
  qualities <- list()
  
  for(i in 1:n_methods) {
    method <- method_names[i]
    recon <- reconstructions[[i]]
    
    if(!is.null(recon$quality)) {
      # Extract useful quality metrics
      q <- recon$quality
      
      # Store common quality metrics
      qualities[[method]] <- list()
      
      # Uncertainty measure (lower is better)
      if(!is.null(q$mean_ci_width)) {
        qualities[[method]]$uncertainty <- q$mean_ci_width
      } else if(!is.null(q$root_uncertainty)) {
        qualities[[method]]$uncertainty <- q$root_uncertainty
      }
      
      # Model fit measure (higher is better)
      if(!is.null(q$lnL)) {
        qualities[[method]]$fit <- q$lnL
      }
      
      # AIC measure (lower is better)
      if(!is.null(q$aic)) {
        qualities[[method]]$aic <- q$aic
      } else if(!is.null(q$aicc)) {
        qualities[[method]]$aic <- q$aicc
      }
    }
  }
  
  # If we have uncertainty measures for all methods
  if(all(sapply(qualities, function(q) !is.null(q$uncertainty)))) {
    # Convert uncertainty to weights (smaller uncertainty = larger weight)
    uncertainties <- sapply(qualities, function(q) q$uncertainty)
    max_uncertainty <- max(uncertainties)
    
    # Invert and normalize to get weights
    inv_uncertainty <- 1 - (uncertainties / max_uncertainty)
    uncertainty_weights <- inv_uncertainty / sum(inv_uncertainty)
    
    # Update weights with uncertainty-based weights
    weights <- uncertainty_weights
  }
  
  # If we have AIC measures for all methods
  if(all(sapply(qualities, function(q) !is.null(q$aic)))) {
    # Convert AIC to weights using Akaike weights
    aics <- sapply(qualities, function(q) q$aic)
    min_aic <- min(aics)
    delta_aic <- aics - min_aic
    
    # Calculate AIC weights
    aic_weights <- exp(-0.5 * delta_aic)
    aic_weights <- aic_weights / sum(aic_weights)
    
    # Update weights using both uncertainty and AIC
    if(all(sapply(qualities, function(q) !is.null(q$uncertainty)))) {
      # Combine uncertainty and AIC weights
      weights <- 0.5 * weights + 0.5 * aic_weights
    } else {
      # Use only AIC weights
      weights <- aic_weights
    }
  }
  
  # Ensure weights are positive and sum to 1
  weights <- pmax(weights, 0.01)  # Minimum weight of 1%
  weights <- weights / sum(weights)
  
  return(weights)
}

#' Weighted mean ensemble method
#' 
#' @param reconstructions List of reconstruction results
#' @param weights Method weights
#' @param node_ids Internal node IDs
#' @param tree Phylogenetic tree
#' @return Ensemble reconstruction result
#' @keywords internal
ensemble_weighted_mean <- function(reconstructions, weights, node_ids, tree) {
  message("Applying weighted mean ensemble method...")
  
  # Initialize data structures for results
  n_nodes <- length(node_ids)
  ensemble_states <- numeric(n_nodes)
  ensemble_ci_lower <- numeric(n_nodes)
  ensemble_ci_upper <- numeric(n_nodes)
  method_states <- matrix(NA, nrow = n_nodes, ncol = length(reconstructions))
  
  # For each node, calculate weighted mean of ancestral states
  for(i in 1:n_nodes) {
    node_id <- node_ids[i]
    node_states <- numeric(length(reconstructions))
    node_ci_lower <- numeric(length(reconstructions))
    node_ci_upper <- numeric(length(reconstructions))
    
    # Extract state and CIs for this node from each method
    for(j in 1:length(reconstructions)) {
      recon <- reconstructions[[j]]
      
      if(!is.null(recon$ancestral_states)) {
        # Find this node in ancestral states
        node_row <- which(recon$ancestral_states$node_id == node_id)
        
        if(length(node_row) > 0) {
          node_states[j] <- recon$ancestral_states$state[node_row]
          method_states[i, j] <- node_states[j]
          
          # Get confidence intervals if available
          if("ci_lower" %in% colnames(recon$ancestral_states) && 
             "ci_upper" %in% colnames(recon$ancestral_states)) {
            node_ci_lower[j] <- recon$ancestral_states$ci_lower[node_row]
            node_ci_upper[j] <- recon$ancestral_states$ci_upper[node_row]
          }
        } else {
          warning("Node ", node_id, " not found in method ", names(reconstructions)[j])
        }
      }
    }
    
    # Calculate weighted mean for node state
    ensemble_states[i] <- sum(node_states * weights, na.rm = TRUE)
    
    # Calculate confidence intervals
    # Method 1: weighted mean of individual CIs
    weighted_range <- sum((node_ci_upper - node_ci_lower) * weights, na.rm = TRUE)
    ensemble_ci_lower[i] <- max(0, ensemble_states[i] - weighted_range/2)
    ensemble_ci_upper[i] <- ensemble_states[i] + weighted_range/2
    
    # Method 2: observed variability between methods plus weighted CI
    state_sd <- sd(node_states, na.rm = TRUE)
    if(!is.na(state_sd) && state_sd > 0) {
      # Add some weight to between-method variance
      ensemble_ci_lower[i] <- max(0, ensemble_states[i] - sqrt(weighted_range^2/4 + state_sd^2))
      ensemble_ci_upper[i] <- ensemble_states[i] + sqrt(weighted_range^2/4 + state_sd^2)
    }
  }
  
  # Create ancestral states data frame
  ancestral_states <- data.frame(
    node_id = node_ids,
    state = ensemble_states,
    ci_lower = ensemble_ci_lower,
    ci_upper = ensemble_ci_upper,
    stringsAsFactors = FALSE
  )
  
  # Compile result
  result <- list(
    ancestral_states = ancestral_states,
    method = "ensemble_weighted_mean",
    method_states = method_states,
    method_names = names(reconstructions)
  )
  
  return(result)
}

#' Median ensemble method
#' 
#' @param reconstructions List of reconstruction results
#' @param node_ids Internal node IDs
#' @param tree Phylogenetic tree
#' @return Ensemble reconstruction result
#' @keywords internal
ensemble_median <- function(reconstructions, node_ids, tree) {
  message("Applying median ensemble method...")
  
  # Initialize data structures for results
  n_nodes <- length(node_ids)
  ensemble_states <- numeric(n_nodes)
  ensemble_ci_lower <- numeric(n_nodes)
  ensemble_ci_upper <- numeric(n_nodes)
  method_states <- matrix(NA, nrow = n_nodes, ncol = length(reconstructions))
  
  # For each node, calculate median of ancestral states
  for(i in 1:n_nodes) {
    node_id <- node_ids[i]
    node_states <- numeric(length(reconstructions))
    node_ci_lower <- numeric(length(reconstructions))
    node_ci_upper <- numeric(length(reconstructions))
    
    # Extract state and CIs for this node from each method
    for(j in 1:length(reconstructions)) {
      recon <- reconstructions[[j]]
      
      if(!is.null(recon$ancestral_states)) {
        # Find this node in ancestral states
        node_row <- which(recon$ancestral_states$node_id == node_id)
        
        if(length(node_row) > 0) {
          node_states[j] <- recon$ancestral_states$state[node_row]
          method_states[i, j] <- node_states[j]
          
          # Get confidence intervals if available
          if("ci_lower" %in% colnames(recon$ancestral_states) && 
             "ci_upper" %in% colnames(recon$ancestral_states)) {
            node_ci_lower[j] <- recon$ancestral_states$ci_lower[node_row]
            node_ci_upper[j] <- recon$ancestral_states$ci_upper[node_row]
          }
        }
      }
    }
    
    # Calculate median for node state
    ensemble_states[i] <- median(node_states, na.rm = TRUE)
    
    # Calculate confidence intervals based on the range of estimates
    valid_states <- node_states[!is.na(node_states)]
    if(length(valid_states) >= 2) {
      # Use quantiles of method results for CI
      ensemble_ci_lower[i] <- max(0, quantile(valid_states, 0.025, na.rm = TRUE))
      ensemble_ci_upper[i] <- quantile(valid_states, 0.975, na.rm = TRUE)
    } else if(length(valid_states) == 1) {
      # Only one valid state, use its CI if available, or make a reasonable estimate
      idx <- which(!is.na(node_states))[1]
      if(!is.na(node_ci_lower[idx]) && !is.na(node_ci_upper[idx])) {
        ensemble_ci_lower[i] <- node_ci_lower[idx]
        ensemble_ci_upper[i] <- node_ci_upper[idx]
      } else {
        # No CI available, use a default range based on the state value
        ensemble_ci_lower[i] <- max(0, ensemble_states[i] - ensemble_states[i] * 0.2)
        ensemble_ci_upper[i] <- ensemble_states[i] + ensemble_states[i] * 0.2
      }
    } else {
      # No valid states, use defaults
      ensemble_ci_lower[i] <- 0
      ensemble_ci_upper[i] <- 0
    }
  }
  
  # Create ancestral states data frame
  ancestral_states <- data.frame(
    node_id = node_ids,
    state = ensemble_states,
    ci_lower = ensemble_ci_lower,
    ci_upper = ensemble_ci_upper,
    stringsAsFactors = FALSE
  )
  
  # Compile result
  result <- list(
    ancestral_states = ancestral_states,
    method = "ensemble_median",
    method_states = method_states,
    method_names = names(reconstructions)
  )
  
  return(result)
}

#' Maximum likelihood ensemble method
#' 
#' @param reconstructions List of reconstruction results
#' @param weights Method weights
#' @param node_ids Internal node IDs
#' @param tree Phylogenetic tree
#' @return Ensemble reconstruction result
#' @keywords internal
ensemble_maximum_likelihood <- function(reconstructions, weights, node_ids, tree) {
  message("Applying maximum likelihood ensemble method...")
  
  # This method would optimally use a more sophisticated approach involving 
  # the likelihood functions from each method, but as a simpler implementation
  # we'll use a weighted approach based on model fit
  
  # Initialize data structures for results
  n_nodes <- length(node_ids)
  ensemble_states <- numeric(n_nodes)
  ensemble_ci_lower <- numeric(n_nodes)
  ensemble_ci_upper <- numeric(n_nodes)
  method_states <- matrix(NA, nrow = n_nodes, ncol = length(reconstructions))
  
  # Get likelihood-based weights if available
  likelihood_weights <- weights
  ml_methods <- c()
  
  # Identify methods with ML characteristics
  for(j in 1:length(reconstructions)) {
    recon <- reconstructions[[j]]
    
    if((!is.null(recon$method) && recon$method == "ML") || 
       (!is.null(recon$model_fit) && !is.null(recon$model_fit$parameters$lnL))) {
      ml_methods <- c(ml_methods, j)
    }
  }
  
  # If we have ML methods, adjust weights to favor them
  if(length(ml_methods) > 0) {
    # Boost weights for ML methods
    boost_factor <- 1.5
    for(j in ml_methods) {
      likelihood_weights[j] <- likelihood_weights[j] * boost_factor
    }
    # Normalize weights
    likelihood_weights <- likelihood_weights / sum(likelihood_weights)
  }
  
  # For each node, calculate weighted mean with emphasis on ML methods
  for(i in 1:n_nodes) {
    node_id <- node_ids[i]
    node_states <- numeric(length(reconstructions))
    node_ci_lower <- numeric(length(reconstructions))
    node_ci_upper <- numeric(length(reconstructions))
    
    # Extract state and CIs for this node from each method
    for(j in 1:length(reconstructions)) {
      recon <- reconstructions[[j]]
      
      if(!is.null(recon$ancestral_states)) {
        # Find this node in ancestral states
        node_row <- which(recon$ancestral_states$node_id == node_id)
        
        if(length(node_row) > 0) {
          node_states[j] <- recon$ancestral_states$state[node_row]
          method_states[i, j] <- node_states[j]
          
          # Get confidence intervals if available
          if("ci_lower" %in% colnames(recon$ancestral_states) && 
             "ci_upper" %in% colnames(recon$ancestral_states)) {
            node_ci_lower[j] <- recon$ancestral_states$ci_lower[node_row]
            node_ci_upper[j] <- recon$ancestral_states$ci_upper[node_row]
          }
        }
      }
    }
    
    # Calculate weighted mean with likelihood-informed weights
    ensemble_states[i] <- sum(node_states * likelihood_weights, na.rm = TRUE)
    
    # Calculate confidence intervals
    # ML-informed confidence intervals
    ml_ranges <- numeric(0)
    for(j in ml_methods) {
      if(!is.na(node_ci_lower[j]) && !is.na(node_ci_upper[j])) {
        ml_ranges <- c(ml_ranges, node_ci_upper[j] - node_ci_lower[j])
      }
    }
    
    if(length(ml_ranges) > 0) {
      # Use ML-based confidence intervals if available
      mean_ml_range <- mean(ml_ranges)
      ensemble_ci_lower[i] <- max(0, ensemble_states[i] - mean_ml_range/2)
      ensemble_ci_upper[i] <- ensemble_states[i] + mean_ml_range/2
    } else {
      # Fallback to weighted average of all CI ranges
      weighted_range <- sum((node_ci_upper - node_ci_lower) * likelihood_weights, na.rm = TRUE)
      if(is.na(weighted_range) || weighted_range == 0) {
        # No valid CI ranges, use a default based on state variability
        state_sd <- sd(node_states, na.rm = TRUE)
        if(!is.na(state_sd) && state_sd > 0) {
          ensemble_ci_lower[i] <- max(0, ensemble_states[i] - 1.96 * state_sd)
          ensemble_ci_upper[i] <- ensemble_states[i] + 1.96 * state_sd
        } else {
          # Use a percentage of the state value
          ensemble_ci_lower[i] <- max(0, ensemble_states[i] - ensemble_states[i] * 0.2)
          ensemble_ci_upper[i] <- ensemble_states[i] + ensemble_states[i] * 0.2
        }
      } else {
        ensemble_ci_lower[i] <- max(0, ensemble_states[i] - weighted_range/2)
        ensemble_ci_upper[i] <- ensemble_states[i] + weighted_range/2
      }
    }
  }
  
  # Create ancestral states data frame
  ancestral_states <- data.frame(
    node_id = node_ids,
    state = ensemble_states,
    ci_lower = ensemble_ci_lower,
    ci_upper = ensemble_ci_upper,
    stringsAsFactors = FALSE
  )
  
  # Compile result
  result <- list(
    ancestral_states = ancestral_states,
    method = "ensemble_maximum_likelihood",
    method_states = method_states,
    method_names = names(reconstructions),
    ml_methods = ml_methods,
    likelihood_weights = likelihood_weights
  )
  
  return(result)
}

#' Bayesian model averaging ensemble method
#' 
#' @param reconstructions List of reconstruction results
#' @param weights Method weights
#' @param node_ids Internal node IDs
#' @param tree Phylogenetic tree
#' @return Ensemble reconstruction result
#' @keywords internal
ensemble_bayesian_model_averaging <- function(reconstructions, weights, node_ids, tree) {
  message("Applying Bayesian model averaging ensemble method...")
  
  # Identify Bayesian methods
  bayesian_methods <- c()
  for(j in 1:length(reconstructions)) {
    recon <- reconstructions[[j]]
    if(!is.null(recon$method) && grepl("Bayesian", recon$method)) {
      bayesian_methods <- c(bayesian_methods, j)
    }
  }
  
  # If we have Bayesian methods, adjust weights
  bayes_weights <- weights
  if(length(bayesian_methods) > 0) {
    # Boost weights for Bayesian methods
    boost_factor <- 2.0
    for(j in bayesian_methods) {
      bayes_weights[j] <- bayes_weights[j] * boost_factor
    }
    # Normalize weights
    bayes_weights <- bayes_weights / sum(bayes_weights)
  } else {
    message("No Bayesian methods found in reconstructions, using standard weighted average")
  }
  
  # Initialize data structures for results
  n_nodes <- length(node_ids)
  ensemble_states <- numeric(n_nodes)
  ensemble_ci_lower <- numeric(n_nodes)
  ensemble_ci_upper <- numeric(n_nodes)
  method_states <- matrix(NA, nrow = n_nodes, ncol = length(reconstructions))
  
  # For each node, apply Bayesian model averaging
  for(i in 1:n_nodes) {
    node_id <- node_ids[i]
    node_states <- numeric(length(reconstructions))
    node_ci_lower <- numeric(length(reconstructions))
    node_ci_upper <- numeric(length(reconstructions))
    
    # Extract state and CIs for this node from each method
    for(j in 1:length(reconstructions)) {
      recon <- reconstructions[[j]]
      
      if(!is.null(recon$ancestral_states)) {
        # Find this node in ancestral states
        node_row <- which(recon$ancestral_states$node_id == node_id)
        
        if(length(node_row) > 0) {
          node_states[j] <- recon$ancestral_states$state[node_row]
          method_states[i, j] <- node_states[j]
          
          # Get confidence intervals if available
          if("ci_lower" %in% colnames(recon$ancestral_states) && 
             "ci_upper" %in% colnames(recon$ancestral_states)) {
            node_ci_lower[j] <- recon$ancestral_states$ci_lower[node_row]
            node_ci_upper[j] <- recon$ancestral_states$ci_upper[node_row]
          }
        }
      }
    }
    
    # Calculate weighted average using Bayesian-adjusted weights
    ensemble_states[i] <- sum(node_states * bayes_weights, na.rm = TRUE)
    
    # Calculate confidence intervals
    if(length(bayesian_methods) > 0) {
      # If we have Bayesian methods, use them for confidence intervals
      bayes_ci_lower <- numeric()
      bayes_ci_upper <- numeric()
      
      for(j in bayesian_methods) {
        if(!is.na(node_ci_lower[j]) && !is.na(node_ci_upper[j])) {
          bayes_ci_lower <- c(bayes_ci_lower, node_ci_lower[j])
          bayes_ci_upper <- c(bayes_ci_upper, node_ci_upper[j])
        }
      }
      
      if(length(bayes_ci_lower) > 0 && length(bayes_ci_upper) > 0) {
        # Use Bayesian posterior credible intervals
        ensemble_ci_lower[i] <- max(0, min(bayes_ci_lower))
        ensemble_ci_upper[i] <- max(bayes_ci_upper)
      } else {
        # Fallback to weighted ranges
        weighted_range <- sum((node_ci_upper - node_ci_lower) * bayes_weights, na.rm = TRUE)
        ensemble_ci_lower[i] <- max(0, ensemble_states[i] - weighted_range/2)
        ensemble_ci_upper[i] <- ensemble_states[i] + weighted_range/2
      }
    } else {
      # No Bayesian methods, use weighted ranges
      weighted_range <- sum((node_ci_upper - node_ci_lower) * bayes_weights, na.rm = TRUE)
      if(is.na(weighted_range) || weighted_range == 0) {
        # No valid ranges, use state distribution
        state_sd <- sd(node_states, na.rm = TRUE)
        if(!is.na(state_sd) && state_sd > 0) {
          ensemble_ci_lower[i] <- max(0, ensemble_states[i] - 1.96 * state_sd)
          ensemble_ci_upper[i] <- ensemble_states[i] + 1.96 * state_sd
        } else {
          # Use a percentage of state value
          ensemble_ci_lower[i] <- max(0, ensemble_states[i] - ensemble_states[i] * 0.2)
          ensemble_ci_upper[i] <- ensemble_states[i] + ensemble_states[i] * 0.2
        }
      } else {
        ensemble_ci_lower[i] <- max(0, ensemble_states[i] - weighted_range/2)
        ensemble_ci_upper[i] <- ensemble_states[i] + weighted_range/2
      }
    }
  }
  
  # Create ancestral states data frame
  ancestral_states <- data.frame(
    node_id = node_ids,
    state = ensemble_states,
    ci_lower = ensemble_ci_lower,
    ci_upper = ensemble_ci_upper,
    stringsAsFactors = FALSE
  )
  
  # Compile result
  result <- list(
    ancestral_states = ancestral_states,
    method = "ensemble_bayesian_model_averaging",
    method_states = method_states,
    method_names = names(reconstructions),
    bayesian_methods = bayesian_methods,
    bayes_weights = bayes_weights
  )
  
  return(result)
}

#' Calculate ensemble reconstruction quality metrics
#' 
#' @param ensemble_result Ensemble reconstruction result
#' @param reconstructions List of original reconstruction results
#' @return Quality assessment metrics
#' @keywords internal
calculate_ensemble_quality <- function(ensemble_result, reconstructions) {
  # Extract ensemble states
  ancestral_states <- ensemble_result$ancestral_states
  
  # Calculate basic quality metrics
  ci_widths <- ancestral_states$ci_upper - ancestral_states$ci_lower
  mean_ci_width <- mean(ci_widths)
  
  # Calculate consistency metrics
  method_states <- ensemble_result$method_states
  
  # Average pairwise differences between methods
  n_methods <- ncol(method_states)
  if(n_methods >= 2) {
    pairwise_diffs <- matrix(0, nrow = n_methods, ncol = n_methods)
    
    for(i in 1:(n_methods-1)) {
      for(j in (i+1):n_methods) {
        # Calculate mean absolute difference between methods
        diffs <- abs(method_states[, i] - method_states[, j])
        pairwise_diffs[i, j] <- mean(diffs, na.rm = TRUE)
        pairwise_diffs[j, i] <- pairwise_diffs[i, j]
      }
    }
    
    # Overall consistency measure (lower values = more consistent)
    consistency <- mean(pairwise_diffs[upper.tri(pairwise_diffs)])
  } else {
    consistency <- NA
  }
  
  # Calculate method agreement as correlation
  method_correlations <- matrix(NA, nrow = n_methods, ncol = n_methods)
  
  for(i in 1:n_methods) {
    for(j in 1:n_methods) {
      if(i != j) {
        non_na <- !is.na(method_states[, i]) & !is.na(method_states[, j])
        if(sum(non_na) >= 3) {
          method_correlations[i, j] <- cor(method_states[non_na, i], 
                                         method_states[non_na, j], 
                                         method = "spearman")
        }
      } else {
        method_correlations[i, j] <- 1.0
      }
    }
  }
  
  mean_correlation <- mean(method_correlations, na.rm = TRUE)
  
  # Calculate variance explained by the ensemble
  total_variance <- sum(apply(method_states, 1, var, na.rm = TRUE))
  ensemble_variance <- sum((sweep(method_states, 1, ancestral_states$state))^2, na.rm = TRUE)
  if(total_variance > 0) {
    explained_variance <- 1 - (ensemble_variance / total_variance)
  } else {
    explained_variance <- 1.0  # All methods give identical results
  }
  
  # Return quality metrics
  quality <- list(
    mean_ci_width = mean_ci_width,
    max_ci_width = max(ci_widths),
    min_ci_width = min(ci_widths),
    method_consistency = consistency,
    mean_correlation = mean_correlation,
    explained_variance = explained_variance,
    method_correlations = method_correlations
  )
  
  return(quality)
}

#===============================================================================
# Reconstruction Comparison and Evaluation Functions
#===============================================================================

#' Compare reconstruction results from different methods
#' 
#' Analyzes and compares ancestral reconstruction results from different methods
#' 
#' @param reconstructions List of reconstruction results
#' @param only_nodes Vector of node IDs to compare (NULL for all)
#' @param reference_method Name of reference method to compare against
#' @return Comparison analysis results
#' @export
compare_reconstructions <- function(reconstructions, only_nodes = NULL, reference_method = NULL) {
  # Check input
  if(!is.list(reconstructions) || length(reconstructions) < 2) {
    stop("At least two reconstruction results required for comparison")
  }
  
  method_names <- names(reconstructions)
  if(is.null(method_names)) {
    method_names <- paste0("method", 1:length(reconstructions))
    names(reconstructions) <- method_names
  }
  
  message(sprintf("Comparing reconstructions from %d methods: %s", 
                 length(reconstructions), paste(method_names, collapse = ", ")))
  
  # If reference method specified, verify it exists
  if(!is.null(reference_method)) {
    if(!reference_method %in% method_names) {
      stop("Reference method '", reference_method, "' not found in reconstructions")
    }
  } else {
    # Default to first method as reference
    reference_method <- method_names[1]
  }
  
  # Extract tree from first reconstruction
  tree <- reconstructions[[1]]$tree
  
  # Get all node IDs
  all_nodes <- sapply(reconstructions, function(r) {
    if(!is.null(r$ancestral_states)) {
      return(r$ancestral_states$node_id)
    } else {
      return(NULL)
    }
  })
  node_ids <- Reduce(intersect, all_nodes)
  
  # Filter to requested nodes if specified
  if(!is.null(only_nodes)) {
    node_ids <- intersect(node_ids, only_nodes)
  }
  
  if(length(node_ids) == 0) {
    stop("No common nodes found across all reconstruction methods")
  }
  
  # Create comparison matrix
  n_nodes <- length(node_ids)
  n_methods <- length(reconstructions)
  states_matrix <- matrix(NA, nrow = n_nodes, ncol = n_methods)
  ci_lower_matrix <- matrix(NA, nrow = n_nodes, ncol = n_methods)
  ci_upper_matrix <- matrix(NA, nrow = n_nodes, ncol = n_methods)
  
  # Fill matrices with state data
  for(i in 1:n_nodes) {
    node_id <- node_ids[i]
    
    for(j in 1:n_methods) {
      recon <- reconstructions[[j]]
      
      if(!is.null(recon$ancestral_states)) {
        # Find this node in ancestral states
        node_row <- which(recon$ancestral_states$node_id == node_id)
        
        if(length(node_row) > 0) {
          states_matrix[i, j] <- recon$ancestral_states$state[node_row]
          
          # Get confidence intervals if available
          if("ci_lower" %in% colnames(recon$ancestral_states) && 
             "ci_upper" %in% colnames(recon$ancestral_states)) {
            ci_lower_matrix[i, j] <- recon$ancestral_states$ci_lower[node_row]
            ci_upper_matrix[i, j] <- recon$ancestral_states$ci_upper[node_row]
          }
        }
      }
    }
  }
  
  # Calculate comparison metrics
  ref_method_idx <- which(method_names == reference_method)
  
  # Calculate difference to reference method
  diff_to_ref <- matrix(NA, nrow = n_nodes, ncol = n_methods)
  rel_diff_to_ref <- matrix(NA, nrow = n_nodes, ncol = n_methods)
  agreement_with_ref <- matrix(NA, nrow = n_nodes, ncol = n_methods)
  
  for(j in 1:n_methods) {
    # Absolute difference to reference
    diff_to_ref[, j] <- states_matrix[, j] - states_matrix[, ref_method_idx]
    
    # Relative difference to reference
    valid_rows <- states_matrix[, ref_method_idx] > 0
    rel_diff_to_ref[valid_rows, j] <- diff_to_ref[valid_rows, j] / states_matrix[valid_rows, ref_method_idx]
    
    # Check if one method's estimate falls within the CI of the other
    for(i in 1:n_nodes) {
      # Check if method j's estimate falls within reference CI
      if(!is.na(ci_lower_matrix[i, ref_method_idx]) && !is.na(ci_upper_matrix[i, ref_method_idx]) && 
         !is.na(states_matrix[i, j])) {
        if(states_matrix[i, j] >= ci_lower_matrix[i, ref_method_idx] && 
           states_matrix[i, j] <= ci_upper_matrix[i, ref_method_idx]) {
          agreement_with_ref[i, j] <- TRUE
        } else {
          agreement_with_ref[i, j] <- FALSE
        }
      }
    }
  }
  
  # Calculate pairwise correlations between methods
  method_correlations <- matrix(NA, nrow = n_methods, ncol = n_methods)
  rownames(method_correlations) <- method_names
  colnames(method_correlations) <- method_names
  
  for(i in 1:n_methods) {
    for(j in 1:n_methods) {
      valid_rows <- !is.na(states_matrix[, i]) & !is.na(states_matrix[, j])
      if(sum(valid_rows) >= 3) {
        method_correlations[i, j] <- cor(states_matrix[valid_rows, i], 
                                       states_matrix[valid_rows, j], 
                                       method = "spearman")
      }
    }
  }
  
  # Calculate method agreement metrics
  method_agreement <- matrix(0, nrow = n_methods, ncol = n_methods)
  rownames(method_agreement) <- method_names
  colnames(method_agreement) <- method_names
  
  for(i in 1:n_methods) {
    for(j in 1:n_methods) {
      if(i != j) {
        # Count nodes where one method's estimate falls within the other's CI
        agreement_count <- 0
        valid_count <- 0
        
        for(node_idx in 1:n_nodes) {
          # i's estimate within j's CI
          if(!is.na(ci_lower_matrix[node_idx, j]) && !is.na(ci_upper_matrix[node_idx, j]) && 
             !is.na(states_matrix[node_idx, i])) {
            valid_count <- valid_count + 1
            if(states_matrix[node_idx, i] >= ci_lower_matrix[node_idx, j] && 
               states_matrix[node_idx, i] <= ci_upper_matrix[node_idx, j]) {
              agreement_count <- agreement_count + 1
            }
          }
        }
        
        if(valid_count > 0) {
          method_agreement[i, j] <- agreement_count / valid_count
        }
      } else {
        method_agreement[i, j] <- 1.0
      }
    }
  }
  
  # Calculate summary statistics
  mean_diff <- colMeans(abs(diff_to_ref), na.rm = TRUE)
  mean_rel_diff <- colMeans(abs(rel_diff_to_ref), na.rm = TRUE)
  agreement_rate <- colMeans(agreement_with_ref, na.rm = TRUE)
  
  # Generate node-by-node comparison data frame
  node_comparison <- data.frame(
    node_id = node_ids,
    stringsAsFactors = FALSE
  )
  
  for(j in 1:n_methods) {
    node_comparison[[method_names[j]]] <- states_matrix[, j]
  }
  
  for(j in 1:n_methods) {
    if(j != ref_method_idx) {
      diff_col_name <- paste0("diff_", method_names[j], "_vs_", reference_method)
      node_comparison[[diff_col_name]] <- diff_to_ref[, j]
    }
  }
  
  # Compile results
  comparison_result <- list(
    node_comparison = node_comparison,
    states_matrix = states_matrix,
    ci_lower_matrix = ci_lower_matrix,
    ci_upper_matrix = ci_upper_matrix,
    method_names = method_names,
    reference_method = reference_method,
    summary = list(
      mean_abs_diff = mean_diff,
      mean_rel_diff = mean_rel_diff,
      agreement_rate = agreement_rate
    ),
    method_correlations = method_correlations,
    method_agreement = method_agreement
  )
  
  return(comparison_result)
}
