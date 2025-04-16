#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Rate Heterogeneity Module
# Author: Bioinformatics Team
# Date: 2025-04-20
# Description: Implements methods for detecting and analyzing heterogeneity in
#              chromosome evolution rates across lineages
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(geiger)
  library(l1ou)  # For OU process shift detection
  library(bayou) # For Bayesian analysis of OU shifts
  library(parallel)
  library(ggplot2)
})

#===============================================================================
# Rate Heterogeneity Detection Functions
#===============================================================================

#' Analyze rate heterogeneity in chromosome evolution
#' 
#' Identifies clades or lineages with different evolutionary rates
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts, named vector with species names
#' @param method Rate analysis method: "brownie" (MEDUSA-like), 
#'               "shift" (rate shift detection), "clade" (clade comparison),
#'               "ou_shift" (OU process shifts), or "bayesian"
#' @param test_models Vector of models to compare, e.g., "BM1", "BMS"
#' @param predefined_clades List of user-defined clades to compare (NULL for automatic)
#' @param n_shifts Maximum number of rate shifts to detect
#' @param min_clade_size Minimum size of clade for rate estimation
#' @param criterion Model selection criterion: "AIC", "AICc", or "BIC"
#' @param n_cores Number of cores for parallel processing (NULL = auto-detect)
#' @return List with results of rate heterogeneity analysis
#' @export
analyze_rate_heterogeneity <- function(tree, chr_counts, 
                                     method = "brownie",
                                     test_models = c("BM1", "BMS"),
                                     predefined_clades = NULL,
                                     n_shifts = 5,
                                     min_clade_size = 5,
                                     criterion = "AICc",
                                     n_cores = NULL) {
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Ensure chr_counts is a named vector
  if(is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector with species names")
  }
  
  # Prune tree and data for common species
  common_taxa <- intersect(tree$tip.label, names(chr_counts))
  if(length(common_taxa) < 5) {
    stop("At least 5 species required for rate heterogeneity analysis")
  }
  
  # Prune tree and data
  tree_pruned <- ape::keep.tip(tree, common_taxa)
  counts_pruned <- chr_counts[common_taxa]
  
  # Initialize results structure
  results <- list(
    tree = tree_pruned,
    data = counts_pruned,
    method = method,
    models = test_models,
    clades = predefined_clades,
    rate_estimates = NULL,
    shifts = NULL,
    model_selection = NULL,
    best_model = NULL
  )
  
  # Proceed with chosen method
  if(method == "brownie") {
    # MEDUSA-like approach using Brownie
    brownie_results <- detect_rate_shifts_brownie(
      tree_pruned, counts_pruned, n_shifts, criterion, n_cores
    )
    results$rate_estimates <- brownie_results$rate_estimates
    results$shifts <- brownie_results$shifts
    results$model_selection <- brownie_results$model_selection
    results$best_model <- brownie_results$best_model
    
  } else if(method == "shift") {
    # Rate shift detection
    shift_results <- detect_rate_shifts(
      tree_pruned, counts_pruned, n_shifts, min_clade_size, criterion
    )
    results$rate_estimates <- shift_results$rate_estimates
    results$shifts <- shift_results$shifts
    results$model_selection <- shift_results$model_selection
    results$best_model <- shift_results$best_model
    
  } else if(method == "clade") {
    # Clade comparison
    if(is.null(predefined_clades)) {
      # Automatically identify major clades
      clades <- identify_major_clades(tree_pruned, min_clade_size)
      results$clades <- clades
    } else {
      results$clades <- predefined_clades
    }
    
    clade_results <- compare_clade_rates(
      tree_pruned, counts_pruned, results$clades, test_models, criterion
    )
    results$rate_estimates <- clade_results$rate_estimates
    results$model_selection <- clade_results$model_selection
    results$best_model <- clade_results$best_model
    results$clade_comparison <- clade_results$clade_comparison
    
  } else if(method == "ou_shift") {
    # OU process shift detection using l1ou
    if(requireNamespace("l1ou", quietly = TRUE)) {
      ou_results <- detect_ou_shifts(
        tree_pruned, counts_pruned, n_shifts, criterion
      )
      results$shifts <- ou_results$shifts
      results$rate_estimates <- ou_results$rate_estimates
      results$alpha_estimate <- ou_results$alpha_estimate
      results$theta_estimates <- ou_results$theta_estimates
      results$model_selection <- ou_results$model_selection
      results$best_model <- "OU_shifts"
    } else {
      warning("l1ou package required for OU shift detection. Falling back to rate shift detection.")
      shift_results <- detect_rate_shifts(
        tree_pruned, counts_pruned, n_shifts, min_clade_size, criterion
      )
      results$rate_estimates <- shift_results$rate_estimates
      results$shifts <- shift_results$shifts
      results$model_selection <- shift_results$model_selection
      results$best_model <- shift_results$best_model
    }
    
  } else if(method == "bayesian") {
    # Bayesian analysis of rate shifts using bayou
    if(requireNamespace("bayou", quietly = TRUE)) {
      bayesian_results <- detect_bayesian_shifts(
        tree_pruned, counts_pruned, n_shifts, n_cores
      )
      results$shifts <- bayesian_results$shifts
      results$rate_estimates <- bayesian_results$rate_estimates
      results$posterior_samples <- bayesian_results$posterior_samples
      results$best_model <- "Bayesian_shifts"
    } else {
      warning("bayou package required for Bayesian shift detection. Falling back to rate shift detection.")
      shift_results <- detect_rate_shifts(
        tree_pruned, counts_pruned, n_shifts, min_clade_size, criterion
      )
      results$rate_estimates <- shift_results$rate_estimates
      results$shifts <- shift_results$shifts
      results$model_selection <- shift_results$model_selection
      results$best_model <- shift_results$best_model
    }
  } else {
    stop(paste("Unsupported method:", method))
  }
  
  # Generate visualizations
  results$plots <- create_rate_visualizations(results)
  
  # Create summary of results
  results$summary <- summarize_rate_analysis(results)
  
  return(results)
}

#' Detect rate shifts using a MEDUSA-like approach
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param n_shifts Maximum number of rate shifts
#' @param criterion Model selection criterion
#' @param n_cores Number of cores for parallel processing
#' @return List with rate shift results
#' @keywords internal
detect_rate_shifts_brownie <- function(tree, chr_counts, n_shifts = 5, 
                                     criterion = "AICc", n_cores = NULL) {
  # Initialize results
  results <- list(
    rate_estimates = NULL,
    shifts = NULL,
    model_selection = NULL,
    best_model = NULL
  )
  
  message("Detecting rate shifts using Brownie-like approach...")
  
  # First fit a single-rate BM model as baseline
  bm1_fit <- geiger::fitContinuous(tree, chr_counts, model = "BM")
  
  # Calculate model fit criteria
  bm1_aic <- bm1_fit$opt$aic
  bm1_aicc <- bm1_fit$opt$aicc
  n_params_bm1 <- 2  # sigma2, z0
  bm1_bic <- -2 * bm1_fit$opt$lnL + n_params_bm1 * log(length(chr_counts))
  
  # Create model selection table
  model_selection <- data.frame(
    Model = "BM1",
    LogLik = bm1_fit$opt$lnL,
    Parameters = n_params_bm1,
    AIC = bm1_aic,
    AICc = bm1_aicc,
    BIC = bm1_bic,
    Rate_Shifts = 0,
    stringsAsFactors = FALSE
  )
  
  # Store baseline rate estimate
  rate_estimates <- data.frame(
    Clade = "All",
    Shift_Node = NA,
    Sigma2 = bm1_fit$opt$sigsq,
    Rate_Relative = 1.0,
    Shift_Support = NA,
    stringsAsFactors = FALSE
  )
  
  # Initialize shift nodes table
  shift_nodes <- data.frame(
    Node = integer(0),
    Improvement = numeric(0),
    Criterion = numeric(0),
    Model = character(0),
    stringsAsFactors = FALSE
  )
  
  # Sequentially test adding shifts
  current_best <- criterion
  current_best_val <- model_selection[[current_best]][1]
  
  # Identify potential shift nodes (internal nodes only)
  internal_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
  
  # Skip root node as it's redundant with the global rate
  potential_shift_nodes <- internal_nodes[-1]
  
  # Setup shift regime list - initially all edges in regime 1
  edge_regimes <- rep(1, nrow(tree$edge))
  
  for(shift in 1:n_shifts) {
    message(sprintf("  Testing for shift %d of %d...", shift, n_shifts))
    
    shift_improvements <- data.frame(
      Node = integer(0),
      LogLik = numeric(0),
      AIC = numeric(0),
      AICc = numeric(0),
      BIC = numeric(0),
      stringsAsFactors = FALSE
    )
    
    # Test adding a shift at each potential node
    for(node in potential_shift_nodes) {
      # Skip nodes where a shift is already assigned
      if(node %in% shift_nodes$Node) {
        next
      }
      
      # Get all descendants of this node
      descendants <- c(node, phytools::getDescendants(tree, node))
      
      # Create a copy of current regimes
      test_regimes <- edge_regimes
      
      # Set new regime for all edges leading to descendants
      for(i in 1:nrow(tree$edge)) {
        if(tree$edge[i, 2] %in% descendants) {
          test_regimes[i] <- shift + 1
        }
      }
      
      # Create regime map
      regime_map <- list(
        edges = tree$edge,
        regimes = test_regimes
      )
      
      # Fit multi-rate model
      tryCatch({
        bms_fit <- geiger::fitContinuous(tree, chr_counts, model = "BMS", regimes = test_regimes)
        
        # Calculate model fit criteria
        bms_aic <- bms_fit$opt$aic
        bms_aicc <- bms_fit$opt$aicc
        n_params_bms <- shift + 2  # sigma2 for each regime, plus z0
        bms_bic <- -2 * bms_fit$opt$lnL + n_params_bms * log(length(chr_counts))
        
        # Record improvement
        shift_improvements <- rbind(shift_improvements, data.frame(
          Node = node,
          LogLik = bms_fit$opt$lnL,
          AIC = bms_aic,
          AICc = bms_aicc,
          BIC = bms_bic,
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        message(sprintf("    Error testing shift at node %d: %s", node, as.character(e)))
      })
    }
    
    # If no valid improvements, stop
    if(nrow(shift_improvements) == 0) {
      message("  No valid shift locations found. Stopping search.")
      break
    }
    
    # Find best shift according to criterion
    best_shift_idx <- which.min(shift_improvements[[criterion]])
    best_shift_node <- shift_improvements$Node[best_shift_idx]
    best_shift_val <- shift_improvements[[criterion]][best_shift_idx]
    
    # Check if this improves the model
    improvement <- current_best_val - best_shift_val
    
    if(improvement > 0) {
      message(sprintf("  Found significant shift at node %d (improvement = %.2f)", 
                     best_shift_node, improvement))
      
      # Update current best
      current_best_val <- best_shift_val
      
      # Record this shift
      shift_nodes <- rbind(shift_nodes, data.frame(
        Node = best_shift_node,
        Improvement = improvement,
        Criterion = best_shift_val,
        Model = paste0("BMS", shift),
        stringsAsFactors = FALSE
      ))
      
      # Update edge regimes for next iteration
      descendants <- c(best_shift_node, phytools::getDescendants(tree, best_shift_node))
      for(i in 1:nrow(tree$edge)) {
        if(tree$edge[i, 2] %in% descendants) {
          edge_regimes[i] <- shift + 1
        }
      }
      
      # Add this model to selection table
      model_selection <- rbind(model_selection, data.frame(
        Model = paste0("BMS", shift),
        LogLik = shift_improvements$LogLik[best_shift_idx],
        Parameters = shift + 2,
        AIC = shift_improvements$AIC[best_shift_idx],
        AICc = shift_improvements$AICc[best_shift_idx],
        BIC = shift_improvements$BIC[best_shift_idx],
        Rate_Shifts = shift,
        stringsAsFactors = FALSE
      ))
      
      # Get rate estimates for this model
      bms_fit <- geiger::fitContinuous(tree, chr_counts, model = "BMS", regimes = edge_regimes)
      
      # Update rate estimates
      for(regime in 1:(shift + 1)) {
        if(regime == 1) {
          clade_name <- "Background"
        } else {
          shift_node <- shift_nodes$Node[regime - 1]
          clade_name <- paste0("Clade_", regime - 1)
        }
        
        # Get edges in this regime
        regime_edges <- which(edge_regimes == regime)
        
        # Find rate for this regime
        sigma2_val <- bms_fit$opt$sigsq[regime]
        
        # Add to rate estimates
        rate_estimates <- rbind(rate_estimates, data.frame(
          Clade = clade_name,
          Shift_Node = ifelse(regime == 1, NA, shift_nodes$Node[regime - 1]),
          Sigma2 = sigma2_val,
          Rate_Relative = sigma2_val / rate_estimates$Sigma2[1],
          Shift_Support = ifelse(regime == 1, NA, shift_nodes$Improvement[regime - 1]),
          stringsAsFactors = FALSE
        ))
      }
    } else {
      message("  No further significant shifts found.")
      break
    }
  }
  
  # Remove the "All" row if we found better models
  if(nrow(model_selection) > 1) {
    rate_estimates <- rate_estimates[-1, ]
  }
  
  # Determine best model
  best_idx <- which.min(model_selection[[criterion]])
  best_model <- model_selection$Model[best_idx]
  
  # Add delta criterion and weights
  model_selection[[paste0("d", criterion)]] <- model_selection[[criterion]] - min(model_selection[[criterion]])
  model_selection[[paste0(criterion, "_weight")]] <- exp(-0.5 * model_selection[[paste0("d", criterion)]]) / 
    sum(exp(-0.5 * model_selection[[paste0("d", criterion)]]))
  
  # Store results
  results$rate_estimates <- rate_estimates
  results$shifts <- shift_nodes
  results$model_selection <- model_selection
  results$best_model <- best_model
  results$edge_regimes <- edge_regimes
  
  return(results)
}

#' Detect rate shifts using direct shift detection
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param n_shifts Maximum number of rate shifts
#' @param min_clade_size Minimum size of clade for rate estimation
#' @param criterion Model selection criterion
#' @return List with rate shift results
#' @keywords internal
detect_rate_shifts <- function(tree, chr_counts, n_shifts = 5, 
                             min_clade_size = 5, criterion = "AICc") {
  # Initialize results
  results <- list(
    rate_estimates = NULL,
    shifts = NULL,
    model_selection = NULL,
    best_model = NULL
  )
  
  message("Detecting rate shifts with direct approach...")
  
  # First fit a single-rate BM model as baseline
  bm1_fit <- geiger::fitContinuous(tree, chr_counts, model = "BM")
  
  # Calculate contrast values along each branch
  contrasts <- calculate_branch_contrasts(tree, chr_counts)
  
  # Identify potential shifts using contrasts
  shift_candidates <- identify_shift_candidates(tree, contrasts, min_clade_size)
  
  # Create model selection table
  model_selection <- data.frame(
    Model = "BM1",
    LogLik = bm1_fit$opt$lnL,
    Parameters = 2,  # sigma2, z0
    AIC = bm1_fit$opt$aic,
    AICc = bm1_fit$opt$aicc,
    BIC = -2 * bm1_fit$opt$lnL + 2 * log(length(chr_counts)),
    Rate_Shifts = 0,
    stringsAsFactors = FALSE
  )
  
  # Store baseline rate estimate
  rate_estimates <- data.frame(
    Clade = "All",
    Shift_Node = NA,
    Sigma2 = bm1_fit$opt$sigsq,
    Rate_Relative = 1.0,
    Shift_Support = NA,
    stringsAsFactors = FALSE
  )
  
  # Define shift nodes
  shift_nodes <- data.frame(
    Node = shift_candidates$Node[1:min(n_shifts, nrow(shift_candidates))],
    Contrast_Value = shift_candidates$Contrast[1:min(n_shifts, nrow(shift_candidates))],
    Criterion = NA,
    Model = NA,
    stringsAsFactors = FALSE
  )
  
  # Test models with increasing numbers of shifts
  for(shift in 1:nrow(shift_nodes)) {
    # Get shift nodes up to current shift
    current_shift_nodes <- shift_nodes$Node[1:shift]
    
    # Create regime map
    edge_regimes <- create_regime_map(tree, current_shift_nodes)
    
    # Fit multi-rate model
    tryCatch({
      bms_fit <- geiger::fitContinuous(tree, chr_counts, model = "BMS", regimes = edge_regimes)
      
      # Calculate model fit criteria
      bms_aic <- bms_fit$opt$aic
      bms_aicc <- bms_fit$opt$aicc
      n_params_bms <- shift + 2  # sigma2 for each regime, plus z0
      bms_bic <- -2 * bms_fit$opt$lnL + n_params_bms * log(length(chr_counts))
      
      # Update shift node info
      shift_nodes$Criterion[shift] <- bms_fit$opt$aicc
      shift_nodes$Model[shift] <- paste0("BMS", shift)
      
      # Add this model to selection table
      model_selection <- rbind(model_selection, data.frame(
        Model = paste0("BMS", shift),
        LogLik = bms_fit$opt$lnL,
        Parameters = n_params_bms,
        AIC = bms_aic,
        AICc = bms_aicc,
        BIC = bms_bic,
        Rate_Shifts = shift,
        stringsAsFactors = FALSE
      ))
      
      # Update rate estimates
      if(criterion == "AICc" && bms_aicc < min(model_selection$AICc[model_selection$Rate_Shifts < shift])) {
        # For each regime
        for(regime in 1:(shift + 1)) {
          # Generate clade name
          if(regime == 1) {
            clade_name <- "Background"
          } else {
            shift_node <- current_shift_nodes[regime - 1]
            clade_name <- paste0("Clade_", regime - 1)
          }
          
          # Find rate for this regime
          sigma2_val <- bms_fit$opt$sigsq[regime]
          
          # Add to rate estimates
          rate_estimates <- rbind(rate_estimates, data.frame(
            Clade = clade_name,
            Shift_Node = ifelse(regime == 1, NA, current_shift_nodes[regime - 1]),
            Sigma2 = sigma2_val,
            Rate_Relative = sigma2_val / bm1_fit$opt$sigsq,
            Shift_Support = shift_candidates$Contrast[shift],
            stringsAsFactors = FALSE
          ))
        }
      }
      
    }, error = function(e) {
      message(sprintf("  Error fitting model with %d shifts: %s", shift, as.character(e)))
    })
  }
  
  # Determine best model
  best_idx <- which.min(model_selection[[criterion]])
  best_model <- model_selection$Model[best_idx]
  
  # Add delta criterion and weights
  model_selection[[paste0("d", criterion)]] <- model_selection[[criterion]] - min(model_selection[[criterion]])
  model_selection[[paste0(criterion, "_weight")]] <- exp(-0.5 * model_selection[[paste0("d", criterion)]]) / 
    sum(exp(-0.5 * model_selection[[paste0("d", criterion)]]))
  
  # Store results
  results$rate_estimates <- rate_estimates
  results$shifts <- shift_nodes
  results$model_selection <- model_selection
  results$best_model <- best_model
  
  return(results)
}

#' Calculate phylogenetic contrasts for each branch
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @return Data frame with branch contrast values
#' @keywords internal
calculate_branch_contrasts <- function(tree, chr_counts) {
  # Calculate ancestral states using ML
  anc <- ape::ace(chr_counts, tree, type = "continuous", method = "ML")
  
  # Extract ancestral state estimates
  node_states <- c(chr_counts, anc$ace)
  
  # Initialize contrast data frame
  contrasts <- data.frame(
    Edge = 1:nrow(tree$edge),
    Parent = tree$edge[, 1],
    Child = tree$edge[, 2],
    Parent_State = node_states[tree$edge[, 1]],
    Child_State = node_states[tree$edge[, 2]],
    Branch_Length = tree$edge.length,
    Contrast = NA,
    Contrast_Abs = NA,
    stringsAsFactors = FALSE
  )
  
  # Calculate contrast for each branch
  contrasts$Contrast <- (contrasts$Child_State - contrasts$Parent_State) / sqrt(contrasts$Branch_Length)
  contrasts$Contrast_Abs <- abs(contrasts$Contrast)
  
  # Calculate Z-scores for contrasts
  contrasts$Z_Score <- (contrasts$Contrast_Abs - mean(contrasts$Contrast_Abs)) / sd(contrasts$Contrast_Abs)
  
  return(contrasts)
}

#' Identify candidate shift nodes from contrasts
#' 
#' @param tree Phylogenetic tree
#' @param contrasts Branch contrast data
#' @param min_clade_size Minimum size of clade for rate estimation
#' @return Data frame of candidate shift nodes
#' @keywords internal
identify_shift_candidates <- function(tree, contrasts, min_clade_size = 5) {
  # Calculate threshold for significant contrasts (e.g., 95th percentile)
  threshold <- quantile(contrasts$Contrast_Abs, 0.95)
  
  # Identify significant contrasts
  sig_contrasts <- contrasts[contrasts$Contrast_Abs > threshold, ]
  
  # Sort by absolute contrast value (descending)
  sig_contrasts <- sig_contrasts[order(-sig_contrasts$Contrast_Abs), ]
  
  # Initialize candidate nodes
  candidates <- data.frame(
    Node = integer(0),
    Contrast = numeric(0),
    Clade_Size = integer(0),
    stringsAsFactors = FALSE
  )
  
  # For each significant contrast, get the child node
  for(i in 1:nrow(sig_contrasts)) {
    node <- sig_contrasts$Child[i]
    
    # Skip tip nodes - we need internal nodes for shifts
    if(node <= length(tree$tip.label)) {
      next
    }
    
    # Get descendants
    descendants <- phytools::getDescendants(tree, node)
    
    # Get number of tips in this clade
    tip_descendants <- descendants[descendants <= length(tree$tip.label)]
    clade_size <- length(tip_descendants)
    
    # Only include if clade is large enough
    if(clade_size >= min_clade_size) {
      candidates <- rbind(candidates, data.frame(
        Node = node,
        Contrast = sig_contrasts$Contrast_Abs[i],
        Clade_Size = clade_size,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(candidates)
}

#' Create regime map for multi-rate model
#' 
#' @param tree Phylogenetic tree
#' @param shift_nodes Vector of nodes where shifts occur
#' @return Vector of regime assignments for each edge
#' @keywords internal
create_regime_map <- function(tree, shift_nodes) {
  # Start with all edges in regime 1
  edge_regimes <- rep(1, nrow(tree$edge))
  
  # For each shift node, assign a new regime to all descendants
  for(i in 1:length(shift_nodes)) {
    # Get shift node
    shift_node <- shift_nodes[i]
    
    # Get all descendants
    descendants <- c(shift_node, phytools::getDescendants(tree, shift_node))
    
    # Set regime for edges leading to descendants
    for(j in 1:nrow(tree$edge)) {
      if(tree$edge[j, 2] %in% descendants) {
        edge_regimes[j] <- i + 1
      }
    }
  }
  
  return(edge_regimes)
}

#' Identify major clades in a phylogenetic tree
#' 
#' @param tree Phylogenetic tree
#' @param min_size Minimum clade size to consider
#' @return List of clades (each is a vector of tip names)
#' @keywords internal
identify_major_clades <- function(tree, min_size = 5) {
  # Find nodes that might define major clades
  n_tips <- length(tree$tip.label)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  
  # Skip root node
  internal_nodes <- internal_nodes[-1]
  
  # Initialize clades list
  clades <- list()
  clade_id <- 1
  
  # For each internal node
  for(node in internal_nodes) {
    # Get tips in this clade
    descendants <- phytools::getDescendants(tree, node)
    tip_descendants <- descendants[descendants <= n_tips]
    
    # Check if clade is large enough
    if(length(tip_descendants) >= min_size) {
      # Get tip names
      clade_tips <- tree$tip.label[tip_descendants]
      
      # Add to clades list
      clades[[clade_id]] <- clade_tips
      clade_id <- clade_id + 1
    }
  }
  
  # Name the clades
  names(clades) <- paste0("Clade_", 1:length(clades))
  
  return(clades)
}

#' Compare evolutionary rates between predefined clades
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param clades List of clades (each is a vector of tip names)
#' @param test_models Vector of models to compare
#' @param criterion Model selection criterion
#' @return List with clade comparison results
#' @keywords internal
compare_clade_rates <- function(tree, chr_counts, clades, 
                               test_models = c("BM1", "BMS"), criterion = "AICc") {
  # Initialize results
  results <- list(
    rate_estimates = NULL,
    model_selection = NULL,
    best_model = NULL,
    clade_comparison = NULL
  )
  
  message("Comparing evolutionary rates between clades...")
  
  # First fit a single-rate BM model as baseline
  bm1_fit <- geiger::fitContinuous(tree, chr_counts, model = "BM")
  
  # Create model selection table
  model_selection <- data.frame(
    Model = "BM1",
    LogLik = bm1_fit$opt$lnL,
    Parameters = 2,  # sigma2, z0
    AIC = bm1_fit$opt$aic,
    AICc = bm1_fit$opt$aicc,
    BIC = -2 * bm1_fit$opt$lnL + 2 * log(length(chr_counts)),
    Clades = 1,
    stringsAsFactors = FALSE
  )
  
  # Store baseline rate estimate
  rate_estimates <- data.frame(
    Clade = "All",
    Sigma2 = bm1_fit$opt$sigsq,
    Rate_Relative = 1.0,
    Tips = length(chr_counts),
    stringsAsFactors = FALSE
  )
  
  # Create clade comparison table
  clade_comparison <- data.frame(
    Clade1 = character(0),
    Clade2 = character(0),
    Rate_Ratio = numeric(0),
    P_Value = numeric(0),
    Significant = logical(0),
    stringsAsFactors = FALSE
  )
  
  # Fit a model with different rates for each clade
  # Create regimes based on clades
  edge_regimes <- rep(1, nrow(tree$edge))
  
  for(i in 1:length(clades)) {
    # Get clade tips
    clade_tips <- clades[[i]]
    
    # Find edges leading to tips in this clade
    for(j in 1:nrow(tree$edge)) {
      child <- tree$edge[j, 2]
      # If child is a tip and in this clade
      if(child <= length(tree$tip.label) && tree$tip.label[child] %in% clade_tips) {
        edge_regimes[j] <- i + 1
      }
    }
  }
  
  # Fit multi-rate model
  tryCatch({
    bms_fit <- geiger::fitContinuous(tree, chr_counts, model = "BMS", regimes = edge_regimes)
    
    # Calculate model fit criteria
    bms_aic <- bms_fit$opt$aic
    bms_aicc <- bms_fit$opt$aicc
    n_params_bms <- length(clades) + 2  # sigma2 for each clade + background, plus z0
    bms_bic <- -2 * bms_fit$opt$lnL + n_params_bms * log(length(chr_counts))
    
    # Add this model to selection table
    model_selection <- rbind(model_selection, data.frame(
      Model = "BMS_clades",
      LogLik = bms_fit$opt$lnL,
      Parameters = n_params_bms,
      AIC = bms_aic,
      AICc = bms_aicc,
      BIC = bms_bic,
      Clades = length(clades) + 1,
      stringsAsFactors = FALSE
    ))
    
    # Add rate estimates for each clade
    for(i in 1:(length(clades) + 1)) {
      if(i == 1) {
        clade_name <- "Background"
        clade_tips <- setdiff(tree$tip.label, unlist(clades))
      } else {
        clade_name <- names(clades)[i - 1]
        clade_tips <- clades[[i - 1]]
      }
      
      sigma2_val <- bms_fit$opt$sigsq[i]
      
      rate_estimates <- rbind(rate_estimates, data.frame(
        Clade = clade_name,
        Sigma2 = sigma2_val,
        Rate_Relative = sigma2_val / bm1_fit$opt$sigsq,
        Tips = length(clade_tips),
        stringsAsFactors = FALSE
      ))
    }
    
    # Compare rates between all pairs of clades
    for(i in 1:(length(clades) + 1)) {
      for(j in (i + 1):(length(clades) + 1)) {
        if(j > (length(clades) + 1)) break
        
        # Get clade names
        if(i == 1) {
          clade1 <- "Background"
        } else {
          clade1 <- names(clades)[i - 1]
        }
        
        if(j == 1) {
          clade2 <- "Background"
        } else {
          clade2 <- names(clades)[j - 1]
        }
        
        # Calculate rate ratio
        rate_ratio <- bms_fit$opt$sigsq[i] / bms_fit$opt$sigsq[j]
        
        # Perform statistical test (simple chi-squared)
        k <- n_params_bms - 1  # Degrees of freedom = number of rates - 1
        chi_stat <- 2 * (bms_fit$opt$lnL - bm1_fit$opt$lnL)
        p_value <- 1 - pchisq(chi_stat, df = k)
        
        # Add to comparison table
        clade_comparison <- rbind(clade_comparison, data.frame(
          Clade1 = clade1,
          Clade2 = clade2,
          Rate_Ratio = rate_ratio,
          P_Value = p_value,
          Significant = p_value < 0.05,
          stringsAsFactors = FALSE
        ))
      }
    }
    
  }, error = function(e) {
    message(sprintf("  Error fitting multi-rate model: %s", as.character(e)))
  })
  
  # Determine best model
  best_idx <- which.min(model_selection[[criterion]])
  best_model <- model_selection$Model[best_idx]
  
  # Add delta criterion and weights
  model_selection[[paste0("d", criterion)]] <- model_selection[[criterion]] - min(model_selection[[criterion]])
  model_selection[[paste0(criterion, "_weight")]] <- exp(-0.5 * model_selection[[paste0("d", criterion)]]) / 
    sum(exp(-0.5 * model_selection[[paste0("d", criterion)]]))
  
  # Store results
  results$rate_estimates <- rate_estimates
  results$model_selection <- model_selection
  results$best_model <- best_model
  results$clade_comparison <- clade_comparison
  
  return(results)
}

#' Detect shifts in OU process using l1ou
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param n_shifts Maximum number of shifts
#' @param criterion Model selection criterion
#' @return List with OU shift results
#' @keywords internal
detect_ou_shifts <- function(tree, chr_counts, n_shifts = 5, criterion = "AICc") {
  # Initialize results
  results <- list(
    shifts = NULL,
    rate_estimates = NULL,
    alpha_estimate = NULL,
    theta_estimates = NULL,
    model_selection = NULL
  )
  
  message("Detecting OU process shifts using l1ou...")
  
  # Check if l1ou package is available
  if(!requireNamespace("l1ou", quietly = TRUE)) {
    stop("l1ou package is required for OU shift detection")
  }
  
  # First fit a single-optimum OU model as baseline
  ou1_fit <- geiger::fitContinuous(tree, chr_counts, model = "OU")
  
  # Convert tree and data to format required by l1ou
  l1ou_data <- l1ou::transf.ouch(tree, chr_counts, scale = FALSE)
  
  # Run l1ou shift detection
  tryCatch({
    l1ou_fit <- l1ou::estimate.shift(l1ou_data$tree, l1ou_data$Y, 
                                   variance = TRUE, criterion = criterion,
                                   nCores = 1, kmax = n_shifts)
    
    # Organize shifts detected by l1ou
    shifts <- l1ou::get.shift.configuration(l1ou_fit)
    
    # Create shift table
    shift_nodes <- data.frame(
      Node = as.integer(names(shifts)),
      Shift_Size = shifts,
      Optimum = NA,
      Regime = NA,
      stringsAsFactors = FALSE
    )
    
    # Get parameter estimates
    theta_estimates <- l1ou_fit$theta
    alpha_estimate <- l1ou_fit$alpha
    sigma2_estimate <- l1ou_fit$sigma2
    
    # Assign optima to shift nodes
    for(i in 1:nrow(shift_nodes)) {
      shift_nodes$Optimum[i] <- theta_estimates[i + 1]  # +1 because first theta is for root regime
      shift_nodes$Regime[i] <- i + 1  # +1 because first regime is root regime
    }
    
    # Store shifts
    results$shifts <- shift_nodes
    
    # Create rate estimates table
    rate_estimates <- data.frame(
      Regime = 1:(length(theta_estimates)),
      Theta = theta_estimates,
      Node = c(NA, shift_nodes$Node),  # NA for root regime
      Sigma2 = sigma2_estimate,
      stringsAsFactors = FALSE
    )
    
    results$rate_estimates <- rate_estimates
    results$alpha_estimate <- alpha_estimate
    results$theta_estimates <- theta_estimates
    
    # Create model selection table
    model_selection <- data.frame(
      Model = c("OU1", paste0("OU", 1:length(shifts), "shifts")),
      LogLik = c(ou1_fit$opt$lnL, l1ou_fit$loglik),
      Parameters = c(3, 3 + 1:length(shifts)),  # alpha, sigma2, theta for OU1, additional theta for each shift
      AIC = c(ou1_fit$opt$aic, l1ou_fit$aic),
      AICc = c(ou1_fit$opt$aicc, l1ou_fit$aicc),
      BIC = c(-2 * ou1_fit$opt$lnL + 3 * log(length(chr_counts)), 
             -2 * l1ou_fit$loglik + (3 + 1:length(shifts)) * log(length(chr_counts))),
      Shifts = c(0, 1:length(shifts)),
      stringsAsFactors = FALSE
    )
    
    results$model_selection <- model_selection
    
  }, error = function(e) {
    message(sprintf("  Error in l1ou analysis: %s", as.character(e)))
  })
  
  return(results)
}

#' Detect shifts using Bayesian MCMC with bayou
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param n_shifts Expected number of shifts (prior)
#' @param n_cores Number of cores for parallel processing
#' @return List with Bayesian shift results
#' @keywords internal
detect_bayesian_shifts <- function(tree, chr_counts, n_shifts = 5, n_cores = NULL) {
  # Initialize results
  results <- list(
    shifts = NULL,
    rate_estimates = NULL,
    posterior_samples = NULL
  )
  
  message("Detecting shifts using Bayesian MCMC with bayou...")
  
  # Check if bayou package is available
  if(!requireNamespace("bayou", quietly = TRUE)) {
    stop("bayou package is required for Bayesian shift detection")
  }
  
  # Set up prior
  prior <- list(
    alpha = list(
      type = "halfnorm",
      param = list(sigma = 0.1)
    ),
    sig2 = list(
      type = "halfnorm", 
      param = list(sigma = 0.1)
    ),
    k = list(
      type = "fixed",
      param = list(k = n_shifts)
    ),
    theta = list(
      type = "normal", 
      param = list(mean = mean(chr_counts), sigma = 2 * sd(chr_counts))
    ),
    slide = list(
      type = "fixed",
      param = list(slide = 0)
    )
  )
  
  # Set up MCMC parameters
  mcmc_params <- list(
    ngen = 20000,
    burnin = 5000,
    sample.freq = 10,
    print.freq = 1000
  )
  
  # Run Bayesian MCMC
  tryCatch({
    # Start with random chain
    set.seed(42)  # For reproducibility
    bayou_chain <- bayou::bayou.mcmc(tree, chr_counts, prior, 
                                   store.mcmc = TRUE, ngen = mcmc_params$ngen,
                                   samp = mcmc_params$sample.freq, 
                                   chunk = mcmc_params$print.freq)
    
    # Discard burn-in
    bayou_chain <- bayou::set.burnin(bayou_chain, mcmc_params$burnin)
    
    # Calculate posterior probabilities for shifts
    bayou_results <- bayou::summary.bayou(bayou_chain)
    
    # Extract shifts with high posterior probability
    shift_summary <- bayou_results$branch.posteriors
    
    # Keep only shifts with posterior probability > 0.3
    shift_nodes <- data.frame(
      Branch = as.integer(rownames(shift_summary)[shift_summary$pp > 0.3]),
      Posterior = shift_summary$pp[shift_summary$pp > 0.3],
      Theta = NA,
      stringsAsFactors = FALSE
    )
    
    # Calculate parameter estimates from posterior
    alpha_estimate <- mean(bayou_chain$alpha)
    sigma2_estimate <- mean(bayou_chain$sig2)
    
    # Extract parameter regimes
    regimes <- bayou::bayou.makeghosts(bayou_chain, tree, 0.3)
    
    # Store results
    results$shifts <- shift_nodes
    results$rate_estimates <- data.frame(
      Alpha = alpha_estimate,
      Sigma2 = sigma2_estimate,
      stringsAsFactors = FALSE
    )
    results$posterior_samples <- bayou_chain
    results$regimes <- regimes
    
  }, error = function(e) {
    message(sprintf("  Error in bayou analysis: %s", as.character(e)))
  })
  
  return(results)
}

#' Create visualizations of rate heterogeneity results
#' 
#' @param results Rate heterogeneity analysis results
#' @return List of plot objects
#' @keywords internal
create_rate_visualizations <- function(results) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package is required for visualizations")
    return(list())
  }
  
  plots <- list()
  
  # Create rate comparison plot
  if(!is.null(results$rate_estimates) && nrow(results$rate_estimates) > 1) {
    rate_data <- results$rate_estimates
    
    # Adjust plot based on method
    if(results$method %in% c("brownie", "shift")) {
      plots$rate_comparison <- ggplot2::ggplot(rate_data, ggplot2::aes(x = Clade, y = Rate_Relative, fill = Clade)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          title = "Relative Evolutionary Rates by Clade",
          subtitle = paste("Method:", results$method),
          x = "Clade",
          y = "Relative Rate (Background = 1.0)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    } else if(results$method == "clade") {
      # Add tip count to label
      rate_data$Label <- paste0(rate_data$Clade, " (n=", rate_data$Tips, ")")
      
      plots$rate_comparison <- ggplot2::ggplot(rate_data, ggplot2::aes(x = reorder(Label, -Rate_Relative), y = Rate_Relative, fill = Clade)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          title = "Relative Evolutionary Rates by Clade",
          subtitle = "Clade comparison method",
          x = "Clade",
          y = "Relative Rate (Background = 1.0)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  }
  
  # Create model selection plot
  if(!is.null(results$model_selection)) {
    model_sel <- results$model_selection
    
    # Determine criterion
    if("AICc_weight" %in% colnames(model_sel)) {
      weight_col <- "AICc_weight"
      criterion <- "AICc"
    } else if("AIC_weight" %in% colnames(model_sel)) {
      weight_col <- "AIC_weight"
      criterion <- "AIC"
    } else if("BIC_weight" %in% colnames(model_sel)) {
      weight_col <- "BIC_weight"
      criterion <- "BIC"
    } else {
      weight_col <- NULL
      criterion <- "AICc"
    }
    
    # Create model selection plot
    if(!is.null(weight_col)) {
      plots$model_selection <- ggplot2::ggplot(model_sel, ggplot2::aes(x = Model, y = get(weight_col), fill = Model)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(
          title = "Model Comparison",
          subtitle = paste("Based on", criterion, "weights"),
          x = "Model",
          y = paste(criterion, "Weight")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  }
  
  # Create rate visualization on tree
  if(requireNamespace("ggtree", quietly = TRUE) && 
     !is.null(results$tree) && !is.null(results$rate_estimates)) {
    
    # For brownie method, we have edge regimes
    if(results$method == "brownie" && !is.null(results$edge_regimes)) {
      # Create data frame with edge colors based on regimes
      edge_data <- data.frame(
        edge = 1:nrow(results$tree$edge),
        regime = results$edge_regimes,
        stringsAsFactors = FALSE
      )
      
      # Create tree with colored edges
      p <- ggtree::ggtree(results$tree)
      
      # Add edge colors
      p <- p + ggtree::aes(color = factor(edge_data$regime)) +
        ggplot2::scale_color_brewer(name = "Rate Regime", palette = "Set1") +
        ggplot2::labs(title = "Chromosome Evolution Rate Regimes") +
        ggplot2::theme(legend.position = "right")
      
      plots$rate_tree <- p
    }
    
    # For ou_shift method, show optimum values
    if(results$method == "ou_shift" && !is.null(results$shifts)) {
      # TO BE IMPLEMENTED: Visualization of OU shifts on tree
    }
  }
  
  return(plots)
}

#' Create summary of rate heterogeneity analysis
#' 
#' @param results Rate heterogeneity analysis results
#' @return Character string with analysis summary
#' @keywords internal
summarize_rate_analysis <- function(results) {
  # Initialize summary
  summary_text <- paste0("Rate Heterogeneity Analysis Summary (", results$method, " method)\n")
  summary_text <- paste0(summary_text, "==================================================\n\n")
  
  # Add model selection information
  if(!is.null(results$best_model)) {
    summary_text <- paste0(summary_text, "Best model: ", results$best_model, "\n\n")
  }
  
  # Add rate information
  if(!is.null(results$rate_estimates) && nrow(results$rate_estimates) > 0) {
    summary_text <- paste0(summary_text, "Rate estimates:\n")
    
    for(i in 1:nrow(results$rate_estimates)) {
      row <- results$rate_estimates[i, ]
      summary_text <- paste0(summary_text, "  ", row$Clade, ": sigma² = ", 
                            round(row$Sigma2, 4), ", relative rate = ", 
                            round(row$Rate_Relative, 2), "\n")
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  # Add shift information
  if(!is.null(results$shifts) && nrow(results$shifts) > 0) {
    summary_text <- paste0(summary_text, "Detected rate shifts:\n")
    
    for(i in 1:nrow(results$shifts)) {
      row <- results$shifts[i, ]
      summary_text <- paste0(summary_text, "  Node ", row$Node, 
                            ifelse("Improvement" %in% colnames(results$shifts), 
                                  paste0(" (improvement = ", round(row$Improvement, 2), ")"), 
                                  ""), "\n")
    }
    summary_text <- paste0(summary_text, "\n")
  }
  
  # Add clade comparison information
  if(!is.null(results$clade_comparison) && nrow(results$clade_comparison) > 0) {
    summary_text <- paste0(summary_text, "Clade rate comparisons:\n")
    
    for(i in 1:nrow(results$clade_comparison)) {
      row <- results$clade_comparison[i, ]
      summary_text <- paste0(summary_text, "  ", row$Clade1, " vs ", row$Clade2, 
                            ": rate ratio = ", round(row$Rate_Ratio, 2),
                            ", p-value = ", round(row$P_Value, 3),
                            ifelse(row$Significant, " *", ""), "\n")
    }
    summary_text <- paste0(summary_text, "\n* = significant difference\n\n")
  }
  
  # Add model fit information
  if(!is.null(results$model_selection)) {
    summary_text <- paste0(summary_text, "Model comparison:\n")
    
    # Get AICc or AIC column
    if("AICc" %in% colnames(results$model_selection)) {
      ic_col <- "AICc"
    } else {
      ic_col <- "AIC"
    }
    
    for(i in 1:nrow(results$model_selection)) {
      row <- results$model_selection[i, ]
      summary_text <- paste0(summary_text, "  ", row$Model, 
                            ": LogLik = ", round(row$LogLik, 2),
                            ", ", ic_col, " = ", round(row[[ic_col]], 2),
                            ifelse(row$Model == results$best_model, " (best model)", ""), "\n")
    }
  }
  
  return(summary_text)
}
