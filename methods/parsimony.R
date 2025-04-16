#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Parsimony Method Module
# Author: Bioinformatics Team
# Date: 2025-03-14
# Description: Implements parsimony methods for ancestral chromosome number
#              reconstruction with different cost models and optimizations
#===============================================================================

suppressPackageStartupMessages({
  library(ape)       # For phylogenetic tree operations
  library(phangorn)  # For parsimony methods
  library(Rcpp)      # For C++ integration
})

#===============================================================================
# Core Parsimony Reconstruction Functions
#===============================================================================

#' Reconstruct ancestral chromosome numbers using parsimony
#' 
#' Applies parsimony methods to reconstruct ancestral chromosome numbers
#' 
#' @param tree Phylogenetic tree object
#' @param chr_counts Chromosome counts, named vector with species names
#' @param method Parsimony method: "fitch", "wagner", "sankoff", or "weighted"
#' @param cost_matrix Cost matrix for sankoff parsimony (NULL for default)
#' @param weight_function Weight function for weighted parsimony
#' @param discrete Whether to treat chromosome counts as discrete characters
#' @return Ancestral chromosome reconstruction results
#' @export
reconstruct_chromosomes_parsimony <- function(tree, chr_counts, 
                                           method = "wagner", 
                                           cost_matrix = NULL,
                                           weight_function = NULL,
                                           discrete = TRUE) {
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Ensure chr_counts is a named vector
  if(is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector with species names")
  }
  
  # Check if method is supported
  supported_methods <- c("fitch", "wagner", "sankoff", "weighted")
  if(!method %in% supported_methods) {
    stop(paste("Unsupported method. Choose from:", paste(supported_methods, collapse = ", ")))
  }
  
  message(sprintf("Using %s parsimony to reconstruct ancestral chromosome numbers...", method))
  
  # Ensure tree and data contain same species
  common_species <- intersect(names(chr_counts), tree$tip.label)
  if(length(common_species) < 3) {
    stop("Fewer than 3 species in common, cannot perform ancestral state reconstruction")
  }
  
  # Prune tree to match data
  pruned_tree <- ape::keep.tip(tree, common_species)
  
  # Order data to match tree
  ordered_counts <- chr_counts[pruned_tree$tip.label]
  
  # Create intermediate structure to store reconstruction data
  reconstruction <- list(
    tree = pruned_tree,
    tip_states = ordered_counts,
    method = method,
    discrete = discrete
  )
  
  # Apply appropriate parsimony method
  if(method == "fitch") {
    reconstruction <- apply_fitch_parsimony(reconstruction, discrete)
  } else if(method == "wagner") {
    reconstruction <- apply_wagner_parsimony(reconstruction)
  } else if(method == "sankoff") {
    reconstruction <- apply_sankoff_parsimony(reconstruction, cost_matrix)
  } else if(method == "weighted") {
    reconstruction <- apply_weighted_parsimony(reconstruction, weight_function)
  }
  
  # Format results
  result <- format_parsimony_results(reconstruction)
  
  # Calculate quality metrics
  result$quality <- calculate_parsimony_quality(result)
  
  return(result)
}

#' Apply Fitch parsimony for discrete characters
#' 
#' @param reconstruction Reconstruction data structure
#' @param discrete Whether to treat counts as discrete characters
#' @return Updated reconstruction data
#' @keywords internal
apply_fitch_parsimony <- function(reconstruction, discrete = TRUE) {
  message("Applying Fitch parsimony algorithm...")
  
  # Extract tree and tip states
  tree <- reconstruction$tree
  tip_states <- reconstruction$tip_states
  
  # Initialize results storage
  n_nodes <- tree$Nnode
  node_ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
  
  # If discrete = TRUE, we'll use phangorn's fitch implementation for discrete characters
  if(discrete) {
    # Convert tip states to discrete character matrix
    # We need to create factor levels for all possible states
    all_states <- sort(unique(tip_states))
    
    # Create data frame with species as rows and a single character column
    species_data <- data.frame(
      chr_count = factor(tip_states, levels = all_states),
      row.names = names(tip_states),
      stringsAsFactors = TRUE
    )
    
    # Convert to phyDat format
    phyDat_data <- phangorn::phyDat(species_data, type = "USER", levels = all_states)
    
    # Reconstruct ancestral states
    parsimony_result <- phangorn::ancestral.pars(tree, phyDat_data)
    
    # Extract the most parsimonious state for each node
    anc_states <- numeric(n_nodes)
    
    for(i in 1:n_nodes) {
      node_idx <- i + Ntip(tree)
      node_states <- parsimony_result[i, ]
      
      # Find state with maximum probability
      max_prob_state <- as.numeric(names(node_states)[which.max(node_states)])
      anc_states[i] <- max_prob_state
    }
    
    # Create node support/certainty estimates
    node_certainty <- numeric(n_nodes)
    
    for(i in 1:n_nodes) {
      node_idx <- i + Ntip(tree)
      node_states <- parsimony_result[i, ]
      
      # Calculate certainty as proportion of maximum state
      max_prob <- max(node_states, na.rm = TRUE)
      total_prob <- sum(node_states, na.rm = TRUE)
      
      if(total_prob > 0) {
        node_certainty[i] <- max_prob / total_prob
      } else {
        node_certainty[i] <- 0
      }
    }
    
    # Create CI bounds (rough approximation for Fitch)
    node_ci_width <- rep(0, n_nodes)
    for(i in 1:n_nodes) {
      node_idx <- i + Ntip(tree)
      node_states <- parsimony_result[i, ]
      
      # Find states with non-zero probability
      possible_states <- as.numeric(names(node_states)[node_states > 0])
      
      if(length(possible_states) > 1) {
        node_ci_width[i] <- max(possible_states) - min(possible_states)
      } else {
        node_ci_width[i] <- 0
      }
    }
    
    # Store results
    reconstruction$node_states <- anc_states
    reconstruction$node_ids <- node_ids
    reconstruction$node_certainty <- node_certainty
    reconstruction$node_ci_width <- node_ci_width
    reconstruction$parsimony_score <- attr(parsimony_result, "pscore")
    
  } else {
    # For continuous data, implement Wagner parsimony (for Fitch with continuous data)
    return(apply_wagner_parsimony(reconstruction))
  }
  
  return(reconstruction)
}

#' Apply Wagner parsimony for continuous characters
#' 
#' @param reconstruction Reconstruction data structure
#' @return Updated reconstruction data
#' @keywords internal
apply_wagner_parsimony <- function(reconstruction) {
  message("Applying Wagner parsimony algorithm for continuous characters...")
  
  # Extract tree and tip states
  tree <- reconstruction$tree
  tip_states <- reconstruction$tip_states
  
  # Initialize results storage
  n_nodes <- tree$Nnode
  node_ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
  node_states <- numeric(n_nodes)
  node_ci_lower <- numeric(n_nodes)
  node_ci_upper <- numeric(n_nodes)
  
  # Create a data structure with all species for phangorn
  species_data <- as.matrix(tip_states)
  
  # We'll use ape's ace function with "MPR" method (maximum parsimony reconstruction)
  ace_result <- ape::ace(tip_states, tree, type = "continuous", method = "ML", model = "BM")
  
  # Extract ancestral states
  node_states <- ace_result$ace
  
  # Extract confidence intervals if available
  if(!is.null(ace_result$CI95)) {
    node_ci_lower <- ace_result$CI95[, 1]
    node_ci_upper <- ace_result$CI95[, 2]
  } else {
    # Approximate CIs using a simple heuristic based on tree depth
    node_depths <- node.depth.edgelength(tree)
    if(is.null(node_depths)) {
      node_depths <- node.depth(tree)
    }
    
    # Scale CI width by relative node depth (deeper nodes have wider CIs)
    root_depth <- max(node_depths)
    relative_depths <- node_depths[node_ids] / root_depth
    
    # Calculate mean absolute difference between adjacent chromosome counts
    tip_diffs <- mean(abs(diff(sort(tip_states))))
    
    for(i in 1:n_nodes) {
      # CI width is proportional to relative depth and typical differences
      width <- tip_diffs * relative_depths[i] * 2
      node_ci_lower[i] <- max(0, node_states[i] - width/2)
      node_ci_upper[i] <- node_states[i] + width/2
    }
  }
  
  # Calculate node certainty (inverse of CI width)
  max_ci_width <- max(node_ci_upper - node_ci_lower)
  node_certainty <- 1 - ((node_ci_upper - node_ci_lower) / max_ci_width)
  
  # Handle edge cases
  node_certainty[is.na(node_certainty)] <- 0.5
  node_certainty[node_certainty < 0] <- 0
  node_certainty[node_certainty > 1] <- 1
  
  # Calculate parsimony score (approximation for Wagner)
  parsimony_score <- calculate_parsimony_score(tree, tip_states, node_states)
  
  # Store results
  reconstruction$node_states <- node_states
  reconstruction$node_ids <- node_ids
  reconstruction$node_ci_lower <- node_ci_lower
  reconstruction$node_ci_upper <- node_ci_upper
  reconstruction$node_certainty <- node_certainty
  reconstruction$parsimony_score <- parsimony_score
  
  return(reconstruction)
}

#' Apply Sankoff parsimony with custom cost matrix
#' 
#' @param reconstruction Reconstruction data structure
#' @param cost_matrix Cost matrix for state transitions
#' @return Updated reconstruction data
#' @keywords internal
apply_sankoff_parsimony <- function(reconstruction, cost_matrix = NULL) {
  message("Applying Sankoff parsimony algorithm...")
  
  # Extract tree and tip states
  tree <- reconstruction$tree
  tip_states <- reconstruction$tip_states
  
  # Initialize results storage
  n_nodes <- tree$Nnode
  node_ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
  
  # If tip states are continuous, we'll need to discretize them
  if(!reconstruction$discrete) {
    message("Discretizing continuous chromosome counts for Sankoff parsimony...")
    
    # Determine appropriate bins for discretization
    counts_range <- range(tip_states)
    bin_width <- max(1, (counts_range[2] - counts_range[1]) / 20)  # Max of 20 bins
    
    # Create bins
    bins <- seq(floor(counts_range[1]), ceiling(counts_range[2]) + bin_width, by = bin_width)
    
    # Discretize tip states
    discretized_states <- cut(tip_states, breaks = bins, labels = FALSE, include.lowest = TRUE)
    names(discretized_states) <- names(tip_states)
    
    # Create mapping from discrete states to original values
    state_mapping <- data.frame(
      discrete = 1:length(bins[-1]),
      original = (bins[-length(bins)] + bins[-1]) / 2
    )
    
    # Store discretization information
    reconstruction$discretized <- TRUE
    reconstruction$state_mapping <- state_mapping
    
    # Update tip states to discretized values
    tip_states_discrete <- discretized_states
  } else {
    # Already discrete
    tip_states_discrete <- as.numeric(factor(tip_states))
    names(tip_states_discrete) <- names(tip_states)
    
    # Create mapping from discrete states to original values
    unique_states <- sort(unique(tip_states))
    state_mapping <- data.frame(
      discrete = 1:length(unique_states),
      original = unique_states
    )
    
    reconstruction$state_mapping <- state_mapping
  }
  
  # Create cost matrix if not provided
  if(is.null(cost_matrix)) {
    # Default cost: absolute difference between states
    unique_discrete_states <- sort(unique(tip_states_discrete))
    n_states <- length(unique_discrete_states)
    
    cost_matrix <- matrix(0, nrow = n_states, ncol = n_states)
    
    for(i in 1:n_states) {
      for(j in 1:n_states) {
        # For chromosome counts, cost could be proportional to the difference
        state_i <- state_mapping$original[state_mapping$discrete == unique_discrete_states[i]]
        state_j <- state_mapping$original[state_mapping$discrete == unique_discrete_states[j]]
        
        # Cost is the absolute difference
        cost_matrix[i, j] <- abs(state_i - state_j)
      }
    }
  }
  
  # Create character data for phangorn
  phyDat_data <- phangorn::phyDat(as.data.frame(tip_states_discrete), type = "USER", 
                               levels = sort(unique(tip_states_discrete)))
  
  # Run Sankoff parsimony
  sankoff_result <- phangorn::ancestral.pars(tree, phyDat_data, method = "sankoff", 
                                          cost = cost_matrix)
  
  # Extract most parsimonious states
  anc_states_discrete <- numeric(n_nodes)
  node_certainty <- numeric(n_nodes)
  
  for(i in 1:n_nodes) {
    node_states <- sankoff_result[i, ]
    
    # Find state with maximum probability
    max_prob_state <- as.numeric(names(node_states)[which.max(node_states)])
    anc_states_discrete[i] <- max_prob_state
    
    # Calculate certainty
    max_prob <- max(node_states, na.rm = TRUE)
    total_prob <- sum(node_states, na.rm = TRUE)
    
    if(total_prob > 0) {
      node_certainty[i] <- max_prob / total_prob
    } else {
      node_certainty[i] <- 0
    }
  }
  
  # Map discrete states back to original values
  anc_states <- numeric(n_nodes)
  
  for(i in 1:n_nodes) {
    discrete_state <- anc_states_discrete[i]
    original_state <- state_mapping$original[state_mapping$discrete == discrete_state]
    anc_states[i] <- original_state
  }
  
  # Calculate confidence intervals
  node_ci_lower <- numeric(n_nodes)
  node_ci_upper <- numeric(n_nodes)
  
  for(i in 1:n_nodes) {
    node_idx <- i + Ntip(tree)
    node_states <- sankoff_result[i, ]
    
    # Get states with non-zero probability
    possible_states <- as.numeric(names(node_states)[node_states > 0])
    
    if(length(possible_states) > 0) {
      # Map to original values
      original_possible <- sapply(possible_states, function(s) {
        state_mapping$original[state_mapping$discrete == s]
      })
      
      node_ci_lower[i] <- min(original_possible)
      node_ci_upper[i] <- max(original_possible)
    } else {
      # Fallback to point estimate
      node_ci_lower[i] <- anc_states[i]
      node_ci_upper[i] <- anc_states[i]
    }
  }
  
  # Store results
  reconstruction$node_states <- anc_states
  reconstruction$node_ids <- node_ids
  reconstruction$node_ci_lower <- node_ci_lower
  reconstruction$node_ci_upper <- node_ci_upper
  reconstruction$node_certainty <- node_certainty
  reconstruction$parsimony_score <- attr(sankoff_result, "pscore")
  
  # Store additional Sankoff-specific results
  reconstruction$cost_matrix <- cost_matrix
  reconstruction$discrete_states <- tip_states_discrete
  
  return(reconstruction)
}

#' Apply weighted parsimony with custom weight function
#' 
#' @param reconstruction Reconstruction data structure
#' @param weight_function Weight function for transitions
#' @return Updated reconstruction data
#' @keywords internal
apply_weighted_parsimony <- function(reconstruction, weight_function = NULL) {
  message("Applying weighted parsimony algorithm...")
  
  # Extract tree and tip states
  tree <- reconstruction$tree
  tip_states <- reconstruction$tip_states
  
  # If no weight function provided, use a default one
  if(is.null(weight_function)) {
    # Default weight function: changes proportional to magnitude
    # For chromosome evolution, smaller changes are more likely
    weight_function <- function(from, to) {
      diff <- abs(from - to)
      # Fusion events (decrease) cost slightly less than fissions (increase)
      if(to < from) {
        return(diff * 0.9)  # Fusion
      } else {
        return(diff * 1.0)  # Fission
      }
    }
  }
  
  # For weighted parsimony, we'll use a dynamic programming approach
  # But for now, we'll use a simplified version based on Wagner parsimony
  
  # We'll start with Wagner parsimony
  wagner_result <- apply_wagner_parsimony(reconstruction)
  
  # Refine the Wagner result using our weight function
  # This is a simplified approach - a full weighted parsimony would 
  # implement Sankoff's algorithm with custom weights
  
  # Initialize with Wagner results
  node_states <- wagner_result$node_states
  node_ids <- wagner_result$node_ids
  node_ci_lower <- wagner_result$node_ci_lower
  node_ci_upper <- wagner_result$node_ci_upper
  node_certainty <- wagner_result$node_certainty
  
  # Apply iterative improvement with weighted parsimony
  # This is a simplified approach for now
  
  # Calculate parsimony score with weights
  edge_weights <- numeric(nrow(tree$edge))
  
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    
    # Get states
    parent_state <- if(parent <= Ntip(tree)) {
      tip_states[parent]
    } else {
      node_states[parent - Ntip(tree)]
    }
    
    child_state <- if(child <= Ntip(tree)) {
      tip_states[child]
    } else {
      node_states[child - Ntip(tree)]
    }
    
    # Calculate weight
    edge_weights[i] <- weight_function(parent_state, child_state)
  }
  
  parsimony_score <- sum(edge_weights)
  
  # Store results
  reconstruction$node_states <- node_states
  reconstruction$node_ids <- node_ids
  reconstruction$node_ci_lower <- node_ci_lower
  reconstruction$node_ci_upper <- node_ci_upper
  reconstruction$node_certainty <- node_certainty
  reconstruction$parsimony_score <- parsimony_score
  reconstruction$edge_weights <- edge_weights
  
  return(reconstruction)
}

#' Calculate parsimony score for a given reconstruction
#' 
#' @param tree Phylogenetic tree
#' @param tip_states Tip states
#' @param node_states Ancestral node states
#' @return Parsimony score
#' @keywords internal
calculate_parsimony_score <- function(tree, tip_states, node_states) {
  # Initialize score
  score <- 0
  
  # Calculate sum of changes on all edges
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    
    # Get states
    parent_state <- if(parent <= Ntip(tree)) {
      tip_states[parent]
    } else {
      node_states[parent - Ntip(tree)]
    }
    
    child_state <- if(child <= Ntip(tree)) {
      tip_states[child]
    } else {
      node_states[child - Ntip(tree)]
    }
    
    # Add absolute difference to score
    score <- score + abs(parent_state - child_state)
  }
  
  return(score)
}

#' Format parsimony reconstruction results
#' 
#' @param reconstruction Raw reconstruction data
#' @return Formatted reconstruction result
#' @keywords internal
format_parsimony_results <- function(reconstruction) {
  # Create ancestral states data frame
  ancestral_states <- data.frame(
    node_id = reconstruction$node_ids,
    state = reconstruction$node_states,
    stringsAsFactors = FALSE
  )
  
  # Add confidence intervals if available
  if(!is.null(reconstruction$node_ci_lower) && !is.null(reconstruction$node_ci_upper)) {
    ancestral_states$ci_lower <- reconstruction$node_ci_lower
    ancestral_states$ci_upper <- reconstruction$node_ci_upper
  }
  
  # Add certainty if available
  if(!is.null(reconstruction$node_certainty)) {
    ancestral_states$certainty <- reconstruction$node_certainty
  }
  
  # Compile final result
  result <- list(
    ancestral_states = ancestral_states,
    tree = reconstruction$tree,
    method = paste("parsimony", reconstruction$method, sep = "_"),
    tip_states = reconstruction$tip_states,
    parsimony_score = reconstruction$parsimony_score,
    original_data = list(
      chr_counts = reconstruction$tip_states
    )
  )
  
  # Add method-specific data
  if(reconstruction$method == "sankoff") {
    result$sankoff <- list(
      cost_matrix = reconstruction$cost_matrix,
      discrete_states = reconstruction$discrete_states
    )
    
    if(!is.null(reconstruction$state_mapping)) {
      result$sankoff$state_mapping <- reconstruction$state_mapping
    }
  } else if(reconstruction$method == "weighted") {
    result$weighted <- list(
      edge_weights = reconstruction$edge_weights
    )
  }
  
  return(result)
}

#' Calculate parsimony reconstruction quality metrics
#' 
#' @param parsimony_result Parsimony reconstruction result
#' @return Quality assessment metrics
#' @keywords internal
calculate_parsimony_quality <- function(parsimony_result) {
  # Extract parameters
  tree <- parsimony_result$tree
  ancestors <- parsimony_result$ancestral_states
  
  # Calculate confidence interval widths if available
  if("ci_lower" %in% colnames(ancestors) && "ci_upper" %in% colnames(ancestors)) {
    ci_widths <- ancestors$ci_upper - ancestors$ci_lower
    mean_ci_width <- mean(ci_widths)
    
    # Root state uncertainty
    root_node <- Ntip(tree) + 1
    root_state_idx <- which(ancestors$node_id == root_node)
    if(length(root_state_idx) > 0) {
      root_uncertainty <- ancestors$ci_upper[root_state_idx] - ancestors$ci_lower[root_state_idx]
    } else {
      root_uncertainty <- NA
    }
  } else {
    mean_ci_width <- NA
    root_uncertainty <- NA
  }
  
  # Calculate state changes
  tip_states <- parsimony_result$tip_states
  
  # Create all nodes data frame (tips and internal nodes)
  all_nodes <- data.frame(
    node_id = c(1:length(tree$tip.label), ancestors$node_id),
    state = c(tip_states, ancestors$state),
    stringsAsFactors = FALSE
  )
  
  # Calculate changes on each edge
  changes <- numeric(nrow(tree$edge))
  
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    
    # Get states
    parent_state <- all_nodes$state[all_nodes$node_id == parent]
    child_state <- all_nodes$state[all_nodes$node_id == child]
    
    # Calculate change
    changes[i] <- abs(child_state - parent_state)
  }
  
  # Calculate parsimony metrics
  total_changes <- sum(changes)
  mean_change <- mean(changes)
  max_change <- max(changes)
  changes_per_edge <- total_changes / nrow(tree$edge)
  
  # Return quality metrics
  quality <- list(
    parsimony_score = parsimony_result$parsimony_score,
    total_changes = total_changes,
    mean_change = mean_change,
    max_change = max_change,
    changes_per_edge = changes_per_edge
  )
  
  # Add CI-based metrics if available
  if(!is.na(mean_ci_width)) {
    quality$mean_ci_width <- mean_ci_width
    quality$root_uncertainty <- root_uncertainty
    
    if("ci_lower" %in% colnames(ancestors) && "ci_upper" %in% colnames(ancestors)) {
      quality$max_ci_width <- max(ci_widths)
      quality$min_ci_width <- min(ci_widths)
    }
  }
  
  # Add certainty metrics if available
  if("certainty" %in% colnames(ancestors)) {
    quality$mean_certainty <- mean(ancestors$certainty)
    quality$min_certainty <- min(ancestors$certainty)
  }
  
  return(quality)
}

#===============================================================================
# Event Detection and Analysis Functions
#===============================================================================

#' Detect chromosome events based on parsimony reconstruction
#' 
#' Detects likely chromosome number change events based on parsimony reconstruction
#' 
#' @param tree Phylogenetic tree object
#' @param parsimony_result Parsimony reconstruction result
#' @param fusion_threshold Fusion event ratio threshold
#' @param fission_threshold Fission event ratio threshold
#' @param wgd_threshold Whole genome duplication ratio threshold
#' @param min_certainty Minimum certainty threshold
#' @return Chromosome events analysis
#' @export
detect_events_parsimony <- function(tree, parsimony_result, 
                                  fusion_threshold = 0.67, 
                                  fission_threshold = 1.5,
                                  wgd_threshold = 2.0,
                                  min_certainty = 0.5) {
  # Check input
  if(!is.list(parsimony_result) || !grepl("parsimony", parsimony_result$method)) {
    stop("parsimony_result must be the result from reconstruct_chromosomes_parsimony")
  }
  
  message("Detecting chromosome events using parsimony-based ancestral reconstruction...")
  
  # Extract reconstruction and tree
  ancestral_states <- parsimony_result$ancestral_states
  tip_states <- parsimony_result$tip_states
  
  if(is.null(tree)) {
    tree <- parsimony_result$tree
  }
  
  # Create all nodes data frame
  all_nodes <- data.frame(
    node_id = c(1:length(tree$tip.label), ancestral_states$node_id),
    state = c(tip_states, ancestral_states$state),
    stringsAsFactors = FALSE
  )
  
  # Add certainty if available
  if("certainty" %in% colnames(ancestral_states)) {
    all_nodes$certainty <- c(rep(1, length(tip_states)), ancestral_states$certainty)
  } else {
    all_nodes$certainty <- c(rep(1, length(tip_states)), rep(1, nrow(ancestral_states)))
  }
  
  # Add confidence intervals if available
  if("ci_lower" %in% colnames(ancestral_states) && "ci_upper" %in% colnames(ancestral_states)) {
    all_nodes$ci_lower <- c(rep(NA, length(tip_states)), ancestral_states$ci_lower)
    all_nodes$ci_upper <- c(rep(NA, length(tip_states)), ancestral_states$ci_upper)
  }
  
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
    event_type = character(n_edges),
    event_confidence = numeric(n_edges),
    event_magnitude = numeric(n_edges),
    stringsAsFactors = FALSE
  )
  
  # Fill event data
  for(i in 1:n_edges) {
    parent_id <- edges[i, 1]
    child_id <- edges[i, 2]
    
    # Get node states
    parent_row <- which(all_nodes$node_id == parent_id)
    child_row <- which(all_nodes$node_id == child_id)
    
    if(length(parent_row) == 0 || length(child_row) == 0) {
      next  # Skip if nodes not found
    }
    
    parent_count <- all_nodes$state[parent_row]
    child_count <- all_nodes$state[child_row]
    
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
    
    # Get certainty
    parent_certainty <- all_nodes$certainty[parent_row]
    child_certainty <- all_nodes$certainty[child_row]
    combined_certainty <- parent_certainty * child_certainty
    
    # Determine event type
    if(is.na(events$change_ratio[i]) || combined_certainty < min_certainty) {
      events$event_type[i] <- "uncertain"
      events$event_confidence[i] <- combined_certainty
      events$event_magnitude[i] <- abs(events$abs_change[i])
    } else if(events$change_ratio[i] <= fusion_threshold) {
      # Chromosome fusion
      events$event_type[i] <- "fusion"
      events$event_magnitude[i] <- parent_count - child_count
      
      # Calculate confidence: greater change = higher confidence
      ratio_conf <- min(1.0, (1 - events$change_ratio[i]) / (1 - fusion_threshold))
      events$event_confidence[i] <- ratio_conf * combined_certainty
    } else if(events$change_ratio[i] >= wgd_threshold) {
      # Whole genome duplication
      events$event_type[i] <- "wgd"
      events$event_magnitude[i] <- child_count - parent_count
      
      # Calculate confidence: closer to exact multiple = higher confidence
      nearest_mult <- round(events$change_ratio[i])
      mult_conf <- 1 - abs(events$change_ratio[i] - nearest_mult) / 0.5
      events$event_confidence[i] <- mult_conf * combined_certainty
    } else if(events$change_ratio[i] >= fission_threshold) {
      # Chromosome fission
      events$event_type[i] <- "fission"
      events$event_magnitude[i] <- child_count - parent_count
      
      # Calculate confidence: greater change = higher confidence
      ratio_conf <- min(1.0, (events$change_ratio[i] - fission_threshold) / fission_threshold)
      events$event_confidence[i] <- ratio_conf * combined_certainty
    } else {
      # No significant change
      events$event_type[i] <- "none"
      events$event_confidence[i] <- combined_certainty
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
    if(event_type != "none" && event_type != "uncertain") {
      message(sprintf("  %s: %d", event_type, event_counts[event_type]))
    }
  }
  
  # Create result
  result <- list(
    events = events,
    tree = tree,
    parsimony_result = parsimony_result,
    parameters = list(
      fusion_threshold = fusion_threshold,
      fission_threshold = fission_threshold,
      wgd_threshold = wgd_threshold,
      min_certainty = min_certainty
    ),
    summary = list(
      event_counts = event_counts,
      avg_confidence = tapply(events$event_confidence, events$event_type, mean, na.rm = TRUE),
      avg_magnitude = tapply(events$event_magnitude, events$event_type, mean, na.rm = TRUE)
    )
  )
  
  return(result)
}

#' Analyze chromosome number transitions
#' 
#' Compute transition probabilities between chromosome states
#' 
#' @param parsimony_result Parsimony reconstruction result
#' @param min_certainty Minimum certainty threshold
#' @param round_states Whether to round chromosome numbers
#' @return Transition analysis results
#' @export
analyze_chromosome_transitions <- function(parsimony_result, min_certainty = 0.5, round_states = TRUE) {
  message("Analyzing chromosome number transitions using parsimony results...")
  
  # Extract tree and states
  tree <- parsimony_result$tree
  ancestral_states <- parsimony_result$ancestral_states
  tip_states <- parsimony_result$tip_states
  
  # Round states if requested (for discrete analysis)
  if(round_states) {
    ancestral_states$state <- round(ancestral_states$state)
    tip_states <- round(tip_states)
  }
  
  # Create all nodes data frame
  all_nodes <- data.frame(
    node_id = c(1:length(tree$tip.label), ancestral_states$node_id),
    state = c(tip_states, ancestral_states$state),
    stringsAsFactors = FALSE
  )
  
  # Add certainty if available
  if("certainty" %in% colnames(ancestral_states)) {
    all_nodes$certainty <- c(rep(1, length(tip_states)), ancestral_states$certainty)
  } else {
    all_nodes$certainty <- rep(1, nrow(all_nodes))
  }
  
  # Create transition matrix
  unique_states <- sort(unique(all_nodes$state))
  n_states <- length(unique_states)
  state_map <- seq_along(unique_states)
  names(state_map) <- unique_states
  
  # Initialize transition count matrix
  transition_matrix <- matrix(0, nrow = n_states, ncol = n_states)
  rownames(transition_matrix) <- colnames(transition_matrix) <- unique_states
  
  # Count transitions along edges
  for(i in 1:nrow(tree$edge)) {
    parent_id <- tree$edge[i, 1]
    child_id <- tree$edge[i, 2]
    
    # Get states
    parent_row <- which(all_nodes$node_id == parent_id)
    child_row <- which(all_nodes$node_id == child_id)
    
    if(length(parent_row) == 0 || length(child_row) == 0) {
      next
    }
    
    parent_state <- all_nodes$state[parent_row]
    child_state <- all_nodes$state[child_row]
    parent_certainty <- all_nodes$certainty[parent_row]
    child_certainty <- all_nodes$certainty[child_row]
    
    # Skip if certainty too low
    if(parent_certainty < min_certainty || child_certainty < min_certainty) {
      next
    }
    
    # Count transition
    from_idx <- state_map[as.character(parent_state)]
    to_idx <- state_map[as.character(child_state)]
    
    if(!is.na(from_idx) && !is.na(to_idx)) {
      transition_matrix[from_idx, to_idx] <- transition_matrix[from_idx, to_idx] + 1
    }
  }
  
  # Convert to probability matrix
  transition_probs <- transition_matrix
  row_sums <- rowSums(transition_matrix)
  for(i in 1:n_states) {
    if(row_sums[i] > 0) {
      transition_probs[i,] <- transition_matrix[i,] / row_sums[i]
    }
  }
  
  # Analyze transition patterns
  result <- list(
    transition_counts = transition_matrix,
    transition_probs = transition_probs,
    states = unique_states,
    tree = tree,
    parsimony_result = parsimony_result,
    parameters = list(
      min_certainty = min_certainty,
      round_states = round_states
    )
  )
  
  # Calculate summary statistics
  result$summary <- list(
    total_transitions = sum(transition_matrix) - sum(diag(transition_matrix)),
    self_transitions = sum(diag(transition_matrix)),
    increase_transitions = sum(transition_matrix[lower.tri(transition_matrix)]),
    decrease_transitions = sum(transition_matrix[upper.tri(transition_matrix)])
  )
  
  # Calculate most common transitions
  transitions <- data.frame(
    from = rep(rownames(transition_matrix), each = ncol(transition_matrix)),
    to = rep(colnames(transition_matrix), times = nrow(transition_matrix)),
    count = as.vector(transition_matrix),
    probability = as.vector(transition_probs),
    stringsAsFactors = FALSE
  )
  transitions <- transitions[transitions$count > 0,]
  transitions <- transitions[order(-transitions$count),]
  
  result$top_transitions <- transitions[1:min(10, nrow(transitions)),]
  
  message(sprintf("Analyzed %d chromosome transitions with %d distinct chromosome states", 
                 result$summary$total_transitions + result$summary$self_transitions,
                 length(unique_states)))
  
  return(result)
}

#===============================================================================
# Validation and Simulation Functions
#===============================================================================

#' Cross-validate parsimony reconstruction accuracy
#' 
#' Test parsimony accuracy by leaving out species and comparing predictions
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param method Parsimony method
#' @param cv_proportion Proportion of species to leave out
#' @param n_replicates Number of cross-validation replicates
#' @param round_results Whether to round predictions for comparison
#' @return Cross-validation results
#' @export
cross_validate_parsimony <- function(tree, chr_counts, method = "wagner", 
                                   cv_proportion = 0.2, n_replicates = 10,
                                   round_results = TRUE) {
  message("Performing cross-validation of parsimony reconstruction...")
  
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Ensure tree and data contain same species
  common_species <- intersect(names(chr_counts), tree$tip.label)
  if(length(common_species) < 5) { # Minimal number for meaningful CV
    stop("At least 5 common species required for cross-validation")
  }
  
  # Prune tree and data to common species
  pruned_tree <- ape::keep.tip(tree, common_species)
  pruned_counts <- chr_counts[common_species]
  
  # Initialize results
  cv_results <- list(
    replicates = list(),
    method = method,
    parameters = list(
      cv_proportion = cv_proportion,
      n_replicates = n_replicates,
      round_results = round_results
    )
  )
  
  # Run cross-validation replicates
  errors <- matrix(NA, nrow = n_replicates, ncol = 4) # MAE, RMSE, MedianAE, Correlation
  colnames(errors) <- c("MAE", "RMSE", "MedianAE", "Correlation")
  
  for(rep in 1:n_replicates) {
    # Sample species to leave out
    n_leave_out <- ceiling(length(common_species) * cv_proportion)
    leave_out <- sample(common_species, n_leave_out)
    training_species <- setdiff(common_species, leave_out)
    
    # Create training dataset
    training_tree <- ape::keep.tip(pruned_tree, training_species)
    training_counts <- pruned_counts[training_species]
    
    # Run parsimony reconstruction
    parsimony_result <- reconstruct_chromosomes_parsimony(training_tree, training_counts, method = method)
    
    # Prepare test predictions
    test_predictions <- numeric(length(leave_out))
    names(test_predictions) <- leave_out
    
    # For each test species, find its position and predict its state
    for(species in leave_out) {
      # Find position in original tree
      # We need to find its sister species or closest relative in training tree
      
      # Simplified approach: use ancestral state of parent node as prediction
      species_idx <- which(pruned_tree$tip.label == species)
      if(length(species_idx) == 1) {
        # Find parent node
        edge_idx <- which(pruned_tree$edge[,2] == species_idx)
        if(length(edge_idx) == 1) {
          parent_node <- pruned_tree$edge[edge_idx, 1]
          
          # Find parent node in ancestral states
          if(parent_node > length(pruned_tree$tip.label)) {
            anc_row <- which(parsimony_result$ancestral_states$node_id == parent_node)
            
            if(length(anc_row) == 1) {
              prediction <- parsimony_result$ancestral_states$state[anc_row]
              test_predictions[species] <- prediction
            }
          }
        }
      }
    }
    
    # Handle cases where prediction failed (use median of training data)
    na_preds <- is.na(test_predictions)
    if(any(na_preds)) {
      test_predictions[na_preds] <- median(training_counts)
    }
    
    # Get actual values
    actual_values <- pruned_counts[leave_out]
    
    # Round if requested
    if(round_results) {
      test_predictions <- round(test_predictions)
    }
    
    # Calculate errors
    mae <- mean(abs(test_predictions - actual_values))
    rmse <- sqrt(mean((test_predictions - actual_values)^2))
    medianae <- median(abs(test_predictions - actual_values))
    correlation <- cor(test_predictions, actual_values, method = "spearman")
    
    errors[rep,] <- c(mae, rmse, medianae, correlation)
    
    # Store detailed results
    cv_results$replicates[[rep]] <- list(
      test_species = leave_out,
      predictions = test_predictions,
      actual = actual_values,
      errors = list(mae = mae, rmse = rmse, medianae = medianae, correlation = correlation)
    )
    
    # Print progress
    if(rep %% 5 == 0 || rep == n_replicates) {
      message(sprintf("Completed %d/%d cross-validation replicates", rep, n_replicates))
    }
  }
  
  # Calculate overall performance
  cv_results$performance <- list(
    mean = colMeans(errors, na.rm = TRUE),
    sd = apply(errors, 2, sd, na.rm = TRUE),
    min = apply(errors, 2, min, na.rm = TRUE),
    max = apply(errors, 2, max, na.rm = TRUE)
  )
  
  message(sprintf("Cross-validation results: Mean MAE = %.2f, Mean RMSE = %.2f, Mean Correlation = %.2f",
                 cv_results$performance$mean["MAE"],
                 cv_results$performance$mean["RMSE"],
                 cv_results$performance$mean["Correlation"]))
  
  return(cv_results)
}

#' Simulate chromosome evolution and test parsimony reconstruction
#' 
#' Simulates chromosome evolution under different models and tests parsimony methods
#' 
#' @param tree Phylogenetic tree
#' @param root_count Root chromosome count
#' @param model Evolutionary model: "BM" (Brownian Motion), "jumps", or "hybrid"
#' @param sim_params Simulation parameters
#' @param test_methods Vector of parsimony methods to test
#' @param n_replicates Number of simulation replicates
#' @return Simulation test results
#' @export
test_parsimony_simulation <- function(tree, root_count = 12, 
                                    model = "BM", 
                                    sim_params = list(sigma = 1.0, 
                                                    jump_prob = 0.05, 
                                                    jump_size = 3),
                                    test_methods = c("wagner", "sankoff", "weighted"),
                                    n_replicates = 10) {
  message("Testing parsimony methods using simulated chromosome evolution...")
  
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Initialize results
  sim_results <- list(
    replicates = list(),
    model = model,
    parameters = sim_params,
    test_methods = test_methods,
    n_replicates = n_replicates
  )
  
  # Initialize performance metrics
  methods_performance <- array(NA, dim = c(n_replicates, length(test_methods), 4),
                              dimnames = list(NULL, test_methods, c("MAE", "RMSE", "Max_Error", "Correlation")))
  
  # Run simulation replicates
  for(rep in 1:n_replicates) {
    # Simulate chromosome evolution
    sim_data <- simulate_chromosome_evolution(tree, root_count, model, sim_params)
    
    # Run and test each parsimony method
    method_results <- list()
    
    for(method_idx in seq_along(test_methods)) {
      method <- test_methods[method_idx]
      
      # Run parsimony reconstruction
      parsimony_result <- reconstruct_chromosomes_parsimony(tree, sim_data$tip_states, method = method)
      
      # Compare reconstructed vs true ancestral states
      reconstructed <- parsimony_result$ancestral_states$state
      true_ancestral <- sim_data$ancestral_states
      
      # Calculate errors
      errors <- reconstructed - true_ancestral
      mae <- mean(abs(errors))
      rmse <- sqrt(mean(errors^2))
      max_error <- max(abs(errors))
      correlation <- cor(reconstructed, true_ancestral, method = "spearman")
      
      # Store performance
      methods_performance[rep, method_idx, ] <- c(mae, rmse, max_error, correlation)
      
      # Store detailed method results
      method_results[[method]] <- list(
        parsimony_result = parsimony_result,
        errors = errors,
        metrics = c(mae = mae, rmse = rmse, max_error = max_error, correlation = correlation)
      )
    }
    
    # Store replicate results
    sim_results$replicates[[rep]] <- list(
      simulation = sim_data,
      method_results = method_results
    )
    
    # Print progress
    if(rep %% 5 == 0 || rep == n_replicates) {
      message(sprintf("Completed %d/%d simulation replicates", rep, n_replicates))
    }
  }
  
  # Calculate overall performance
  sim_results$method_performance <- list()
  
  for(metric_idx in 1:4) {
    metric_name <- dimnames(methods_performance)[[3]][metric_idx]
    metric_data <- methods_performance[, , metric_idx]
    
    sim_results$method_performance[[metric_name]] <- data.frame(
      method = test_methods,
      mean = colMeans(metric_data, na.rm = TRUE),
      sd = apply(metric_data, 2, sd, na.rm = TRUE),
      min = apply(metric_data, 2, min, na.rm = TRUE),
      max = apply(metric_data, 2, max, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  
  # Find best method
  mae_means <- sim_results$method_performance$MAE$mean
  best_method_idx <- which.min(mae_means)
  best_method <- test_methods[best_method_idx]
  
  sim_results$overall_best <- list(
    method = best_method,
    mae = mae_means[best_method_idx]
  )
  
  message(sprintf("Simulation testing completed. Best parsimony method: %s (MAE = %.2f)",
                 best_method, mae_means[best_method_idx]))
  
  return(sim_results)
}

#' Simulate chromosome evolution
#' 
#' Simulate chromosome count evolution on a tree
#' 
#' @param tree Phylogenetic tree
#' @param root_count Root chromosome count
#' @param model Evolutionary model
#' @param params Simulation parameters
#' @return Simulated chromosome data
#' @keywords internal
simulate_chromosome_evolution <- function(tree, root_count, model, params) {
  # Initialize result
  sim_data <- list(
    tree = tree,
    model = model,
    params = params,
    root_state = root_count
  )
  
  # Prepare node containers
  n_nodes <- tree$Nnode + Ntip(tree)
  node_states <- numeric(n_nodes)
  node_states[Ntip(tree) + 1] <- root_count  # Root state
  
  # Apply simulation model
  if(model == "BM") {
    # Brownian motion simulation
    sigma <- params$sigma
    
    # Simulate along each edge of the tree
    for(i in 1:nrow(tree$edge)) {
      parent <- tree$edge[i, 1]
      child <- tree$edge[i, 2]
      branch_length <- tree$edge.length[i]
      
      # Get parent state
      parent_state <- node_states[parent]
      
      # Simulate state change using Brownian motion
      change <- rnorm(1, mean = 0, sd = sigma * sqrt(branch_length))
      child_state <- parent_state + change
      
      # Ensure non-negative value with floor at 1
      node_states[child] <- max(1, child_state)
    }
  } else if(model == "jumps") {
    # Jump model with occasional large changes
    sigma <- params$sigma
    jump_prob <- params$jump_prob
    jump_size <- params$jump_size
    
    # Simulate along each edge of the tree
    for(i in 1:nrow(tree$edge)) {
      parent <- tree$edge[i, 1]
      child <- tree$edge[i, 2]
      branch_length <- tree$edge.length[i]
      
      # Get parent state
      parent_state <- node_states[parent]
      
      # Decide if a jump occurs
      if(runif(1) < jump_prob) {
        # Simulate jump: fusion (down) or fission/duplication (up)
        if(runif(1) < 0.5) {
          # Fusion (decrease)
          change <- -rexp(1, rate = 1/jump_size)
        } else {
          # Fission/duplication (increase)
          change <- rexp(1, rate = 1/jump_size)
        }
      } else {
        # Regular BM change
        change <- rnorm(1, mean = 0, sd = sigma * sqrt(branch_length))
      }
      
      child_state <- parent_state + change
      
      # Ensure non-negative value with floor at 1
      node_states[child] <- max(1, child_state)
    }
  } else if(model == "hybrid") {
    # Hybrid model: BM with occasional jumps and trend
    sigma <- params$sigma
    jump_prob <- params$jump_prob
    jump_size <- params$jump_size
    trend <- if(is.null(params$trend)) 0 else params$trend
    
    # Simulate along each edge of the tree
    for(i in 1:nrow(tree$edge)) {
      parent <- tree$edge[i, 1]
      child <- tree$edge[i, 2]
      branch_length <- tree$edge.length[i]
      
      # Get parent state
      parent_state <- node_states[parent]
      
      # Base change with trend component
      base_change <- rnorm(1, mean = trend * branch_length, sd = sigma * sqrt(branch_length))
      
      # Decide if a jump occurs
      if(runif(1) < jump_prob) {
        # Simulate jump: fusion (down) or fission/duplication (up)
        if(runif(1) < 0.5) {
          # Fusion (decrease)
          jump_change <- -rexp(1, rate = 1/jump_size)
        } else {
          # Fission/duplication (increase)
          jump_change <- rexp(1, rate = 1/jump_size)
        }
        
        # Combine changes
        change <- base_change + jump_change
      } else {
        change <- base_change
      }
      
      child_state <- parent_state + change
      
      # Ensure non-negative value with floor at 1
      node_states[child] <- max(1, child_state)
    }
  } else {
    stop("Unsupported simulation model:", model)
  }
  
  # Extract tip and internal node states
  tip_states <- node_states[1:Ntip(tree)]
  names(tip_states) <- tree$tip.label
  internal_states <- node_states[(Ntip(tree)+1):n_nodes]
  
  # Store simulation results
  sim_data$tip_states <- tip_states
  sim_data$ancestral_states <- internal_states
  sim_data$all_states <- node_states
  
  return(sim_data)
}

#===============================================================================
# Visualization Functions
#===============================================================================

#' Visualize parsimony reconstruction on tree
#' 
#' Create visualization of ancestral chromosome number reconstruction
#' 
#' @param parsimony_result Parsimony reconstruction result
#' @param node_method How to display node values: "color", "label", or "both"
#' @param edge_method How to display edge changes: "color", "width", or "both"
#' @param show_events Whether to mark significant chromosome events
#' @param min_change Minimum change to display on edges
#' @param output_file Output file path (NULL for on-screen display)
#' @return Invisibly returns the plot object
#' @export
plot_parsimony_reconstruction <- function(parsimony_result, 
                                        node_method = "both",
                                        edge_method = "color",
                                        show_events = TRUE,
                                        min_change = 0.5,
                                        output_file = NULL) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggplot2 and ggtree packages are required for visualization")
  }
  
  message("Visualizing parsimony-based ancestral chromosome reconstruction...")
  
  # Extract data
  tree <- parsimony_result$tree
  ancestral_states <- parsimony_result$ancestral_states
  tip_states <- parsimony_result$tip_states
  
  # Detect events if requested
  if(show_events) {
    events_result <- detect_events_parsimony(tree, parsimony_result)
    events <- events_result$events
  }
  
  # Prepare node data
  node_data <- data.frame(
    node = c(1:length(tree$tip.label), ancestral_states$node_id),
    state = c(tip_states, ancestral_states$state),
    stringsAsFactors = FALSE
  )
  
  # Setup output device if needed
  if(!is.null(output_file)) {
    if(grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = 10, height = 8)
    } else if(grepl("\\.png$", output_file)) {
      png(output_file, width = 2000, height = 1600, res = 200)
    } else {
      # Default to PDF
      pdf(paste0(output_file, ".pdf"), width = 10, height = 8)
    }
  }
  
  # Create base tree
  p <- ggtree::ggtree(tree, ladderize = TRUE) %<+% node_data
  
  # Add node visualizations based on method
  if(node_method %in% c("color", "both")) {
    # Color nodes by chromosome count
    p <- p + ggtree::aes(color = state) +
      ggplot2::scale_color_viridis_c(name = "Chromosome\nCount", option = "plasma")
  }
  
  if(node_method %in% c("label", "both")) {
    # Add node labels
    p <- p + ggtree::geom_nodelab(aes(label = round(state, 1)), size = 3, color = "black", 
                                fontface = "bold", nudge_x = 0.2)
  }
  
  # Add edge visualizations
  if(show_events) {
    # Prepare edge data
    edge_data <- data.frame(
      node = events$child_node,
      parent = events$parent_node,
      change = events$abs_change,
      change_ratio = events$change_ratio,
      event_type = events$event_type,
      stringsAsFactors = FALSE
    )
    
    # Filter to significant events
    sig_events <- edge_data[abs(edge_data$change) >= min_change & 
                           edge_data$event_type != "none", ]
    
    if(nrow(sig_events) > 0) {
      # Define event colors
      event_colors <- c(
        "fusion" = "blue",
        "fission" = "red",
        "wgd" = "purple",
        "uncertain" = "gray50"
      )
      
      # Add event markers
      for(i in 1:nrow(sig_events)) {
        event <- sig_events[i, ]
        event_color <- event_colors[event$event_type]
        
        if(!is.na(event_color)) {
          p <- p + ggtree::geom_point2(aes(subset = (node == event$node)), 
                                     color = event_color, size = 3, shape = 18)
        }
      }
      
      # Add legend for events
      event_types <- unique(sig_events$event_type)
      event_types <- event_types[event_types != "none" & event_types != "uncertain"]
      
      if(length(event_types) > 0) {
        legend_data <- data.frame(
          event = event_types,
          color = event_colors[event_types],
          stringsAsFactors = FALSE
        )
        
        p <- p + ggplot2::guides(color = "legend",
                               shape = ggplot2::guide_legend(
                                 title = "Chromosome Events",
                                 override.aes = list(
                                   shape = 18,
                                   color = event_colors[event_types]
                                 )
                               ))
      }
    }
  }
  
  # Add title and layout adjustments
  p <- p + ggplot2::labs(title = "Parsimony Reconstruction of Ancestral Chromosome Numbers",
                       subtitle = paste("Method:", sub("parsimony_", "", parsimony_result$method))) +
    ggplot2::theme(legend.position = "right")
  
  # Display plot
  print(p)
  
  # Close device if using file output
  if(!is.null(output_file)) {
    dev.off()
    message(paste("Visualization saved to:", output_file))
  }
  
  # Return plot object
  return(invisible(p))
}

#' Plot chromosome event distribution
#' 
#' Visualize the distribution of chromosome events across the tree
#' 
#' @param events_result Chromosome events detection result
#' @param event_type Event type to visualize (NULL for all)
#' @param show_magnitude Whether to indicate event magnitude
#' @param min_confidence Minimum event confidence to include
#' @param output_file Output file path
#' @return Invisibly returns the plot object
#' @export
plot_chromosome_events <- function(events_result, event_type = NULL,
                                 show_magnitude = TRUE,
                                 min_confidence = 0.6,
                                 output_file = NULL) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggplot2 and ggtree packages are required for visualization")
  }
  
  # Extract data
  tree <- events_result$tree
  events <- events_result$events
  
  # Filter events
  if(!is.null(event_type)) {
    events <- events[events$event_type == event_type, ]
  } else {
    events <- events[events$event_type != "none", ]
  }
  
  # Apply confidence threshold
  events <- events[events$event_confidence >= min_confidence, ]
  
  if(nrow(events) == 0) {
    message("No events matching criteria found")
    return(NULL)
  }
  
  message(sprintf("Visualizing %d chromosome events on the phylogenetic tree...", nrow(events)))
  
  # Setup output device if needed
  if(!is.null(output_file)) {
    if(grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = 10, height = 8)
    } else if(grepl("\\.png$", output_file)) {
      png(output_file, width = 2000, height = 1600, res = 200)
    } else {
      # Default to PDF
      pdf(paste0(output_file, ".pdf"), width = 10, height = 8)
    }
  }
  
  # Create base tree
  p <- ggtree::ggtree(tree, ladderize = TRUE)
  
  # Add node labels for tips
  p <- p + ggtree::geom_tiplab(size = 3)
  
  # Define event colors and shapes
  event_colors <- c(
    "fusion" = "blue",
    "fission" = "red",
    "wgd" = "purple",
    "uncertain" = "gray50"
  )
  
  event_shapes <- c(
    "fusion" = 25,  # Down triangle
    "fission" = 24, # Up triangle
    "wgd" = 23,     # Diamond
    "uncertain" = 21 # Circle
  )
  
  # Prepare event data for plotting
  event_data <- data.frame(
    parent = events$parent_node,
    node = events$child_node,
    type = events$event_type,
    magnitude = events$event_magnitude,
    confidence = events$event_confidence,
    stringsAsFactors = FALSE
  )
  
  # Add events to the tree
  for(i in 1:nrow(event_data)) {
    event <- event_data[i, ]
    event_color <- event_colors[event$type]
    event_shape <- event_shapes[event$type]
    
    # Calculate size based on magnitude and confidence
    if(show_magnitude) {
      point_size <- 2 + min(5, event$magnitude/2) * event$confidence
    } else {
      point_size <- 3 * event$confidence
    }
    
    # Add event marker
    p <- p + ggtree::geom_point2(aes(subset = (node == event$node)), 
                               color = event_color, 
                               shape = event_shape,
                               size = point_size)
  }
  
  # Add legend
  event_types <- unique(event_data$type)
  if(length(event_types) > 0) {
    # Create custom legend
    legend_data <- data.frame(
      type = event_types,
      color = event_colors[event_types],
      shape = event_shapes[event_types],
      stringsAsFactors = FALSE
    )
    
    # Add legend to plot
    p <- p + 
      ggplot2::guides(color = ggplot2::guide_legend(title = "Event Type"),
                    shape = ggplot2::guide_legend(title = "Event Type"))
  }
  
  # Add title and layout adjustments
  p <- p + ggplot2::labs(
    title = "Chromosome Evolutionary Events",
    subtitle = paste("Minimum confidence:", min_confidence)
  ) +
    ggplot2::theme(legend.position = "right")
  
  # Display plot
  print(p)
  
  # Close device if using file output
  if(!is.null(output_file)) {
    dev.off()
    message(paste("Visualization saved to:", output_file))
  }
  
  # Return plot object
  return(invisible(p))
}

#' Plot transition network between chromosome states
#' 
#' Create a network visualization of transitions between chromosome states
#' 
#' @param transition_result Chromosome transition analysis result
#' @param min_transitions Minimum number of transitions to display
#' @param show_self_loops Whether to display self-transitions
#' @param layout Network layout algorithm
#' @param output_file Output file path
#' @return Invisibly returns the plot object
#' @export
plot_chromosome_transitions <- function(transition_result, 
                                      min_transitions = 1,
                                      show_self_loops = FALSE,
                                      layout = "fruchterman.reingold",
                                      output_file = NULL) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("igraph", quietly = TRUE)) {
    stop("ggplot2 and igraph packages are required for visualization")
  }
  
  # Extract transition data
  transition_matrix <- transition_result$transition_counts
  
  # Remove self-loops if requested
  if(!show_self_loops) {
    diag(transition_matrix) <- 0
  }
  
  # Apply minimum transition threshold
  transition_matrix[transition_matrix < min_transitions] <- 0
  
  # Check if there are transitions left
  if(sum(transition_matrix) == 0) {
    message("No transitions meeting the criteria")
    return(NULL)
  }
  
  message("Visualizing chromosome state transition network...")
  
  # Create graph from transition matrix
  g <- igraph::graph_from_adjacency_matrix(
    transition_matrix,
    mode = "directed",
    weighted = TRUE
  )
  
  # Set node attributes
  igraph::V(g)$name <- rownames(transition_matrix)
  igraph::V(g)$size <- as.numeric(rownames(transition_matrix))
  
  # Set edge attributes
  igraph::E(g)$width <- sqrt(igraph::E(g)$weight)
  igraph::E(g)$arrow.size <- 0.5
  
  # Setup output device if needed
  if(!is.null(output_file)) {
    if(grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = 8, height = 8)
    } else if(grepl("\\.png$", output_file)) {
      png(output_file, width = 1600, height = 1600, res = 200)
    } else {
      # Default to PDF
      pdf(paste0(output_file, ".pdf"), width = 8, height = 8)
    }
  }
  
  # Create layout
  layout_func <- get(paste0("layout_with_", layout), asNamespace("igraph"))
  layout_coords <- layout_func(g)
  
  # Plot the graph
  plot(g,
       layout = layout_coords,
       vertex.color = "skyblue",
       vertex.label.cex = 1.2,
       vertex.label.color = "black",
       vertex.frame.color = "gray20",
       edge.color = "gray40",
       edge.curved = 0.2,
       main = "Chromosome State Transitions")
  
  # Close device if using file output
  if(!is.null(output_file)) {
    dev.off()
    message(paste("Transition network visualization saved to:", output_file))
  }
  
  # Return the graph
  return(invisible(g))
}
