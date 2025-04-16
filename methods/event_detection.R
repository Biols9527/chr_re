#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Event Detection Module
# Author: Bioinformatics Team
# Date: 2025-04-25
# Description: Implements methods for detecting and analyzing chromosome 
#              number change events along a phylogeny
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggplot2)
  library(ggtree)
  library(parallel)
})

#===============================================================================
# Event Detection Core Functions
#===============================================================================

#' Detect chromosome number change events along a phylogeny
#' 
#' Identifies putative events such as whole genome duplications (WGD),
#' chromosome fusions, fissions, and other changes based on ancestral state
#' reconstruction results
#' 
#' @param tree Phylogenetic tree
#' @param reconstruction Ancestral state reconstruction results
#' @param event_types Types of events to detect: "wgd", "fusion", "fission", "gain", "loss", or "all"
#' @param threshold_ratio For WGD detection, minimum ratio of chromosome number increase
#' @param min_change Minimum absolute change to consider an event
#' @param method Detection method: "discrete" or "continuous"
#' @param genome_data Optional genome structure data for refined event detection
#' @param probability_threshold Minimum posterior probability to consider an event
#' @return List with detected events and statistics
#' @export
detect_chromosome_events <- function(tree, 
                                   reconstruction, 
                                   event_types = "all",
                                   threshold_ratio = 1.8, 
                                   min_change = 1,
                                   method = "discrete",
                                   genome_data = NULL,
                                   probability_threshold = 0.8) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Extract reconstructed states
  if(is.data.frame(reconstruction$ancestral_states)) {
    anc_states <- reconstruction$ancestral_states$state
    node_ids <- reconstruction$ancestral_states$node_id
  } else if(is.vector(reconstruction$ancestral_states)) {
    anc_states <- reconstruction$ancestral_states
    node_ids <- names(anc_states)
    if(is.null(node_ids)) {
      node_ids <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)
    }
  } else {
    stop("Unrecognized format for ancestral state reconstruction")
  }
  
  # Determine event types to detect
  all_event_types <- c("wgd", "fusion", "fission", "gain", "loss")
  if(length(event_types) == 1 && event_types == "all") {
    event_types <- all_event_types
  } else {
    event_types <- match.arg(event_types, all_event_types, several.ok = TRUE)
  }
  
  # Handle uncertainties if available
  has_uncertainty <- FALSE
  if(!is.null(reconstruction$state_probabilities) || 
     !is.null(reconstruction$confidence_intervals)) {
    has_uncertainty <- TRUE
    
    if(!is.null(reconstruction$state_probabilities)) {
      # Discrete state probabilities
      state_probs <- reconstruction$state_probabilities
    } else {
      # Continuous confidence intervals
      ci_lower <- reconstruction$confidence_intervals$lower
      ci_upper <- reconstruction$confidence_intervals$upper
    }
  }
  
  # Initialize event data frame
  events <- data.frame(
    Edge = integer(0),
    Parent_Node = integer(0),
    Child_Node = integer(0),
    Parent_State = numeric(0),
    Child_State = numeric(0),
    Change = numeric(0),
    Rel_Change = numeric(0),
    Event_Type = character(0),
    Probability = numeric(0),
    Branch_Length = numeric(0),
    Description = character(0),
    stringsAsFactors = FALSE
  )
  
  # For each edge in the tree
  n_edges <- nrow(tree$edge)
  for(i in 1:n_edges) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get reconstructed states
    if(child <= length(tree$tip.label)) {
      # Child is a tip
      child_state <- get_tip_state(tree, reconstruction, child)
      child_is_tip <- TRUE
    } else {
      # Child is an internal node
      child_node_idx <- which(node_ids == child)
      child_state <- anc_states[child_node_idx]
      child_is_tip <- FALSE
    }
    
    # Get parent state
    parent_node_idx <- which(node_ids == parent)
    parent_state <- anc_states[parent_node_idx]
    
    # Calculate change
    change <- child_state - parent_state
    rel_change <- child_state / parent_state
    
    # Determine event probability
    if(has_uncertainty) {
      if(!is.null(reconstruction$state_probabilities)) {
        # Use state probabilities
        if(child_is_tip) {
          # For tips, we're certain about the state
          prob <- 1.0
        } else {
          # For internal nodes, get probability of the reconstructed state
          prob <- state_probs[child_node_idx, as.character(round(child_state))]
          # If NA, use a default (might not have probabilities for all states)
          if(is.na(prob)) prob <- probability_threshold
        }
      } else {
        # Use confidence intervals - probability based on how far state is from CI bounds
        ci_width <- ci_upper[child_node_idx] - ci_lower[child_node_idx]
        if(child_state >= ci_lower[child_node_idx] && 
           child_state <= ci_upper[child_node_idx]) {
          # State is within CI
          prob <- 1.0
        } else {
          # State is outside CI - calculate based on distance from boundary
          dist_from_boundary <- min(abs(child_state - ci_lower[child_node_idx]),
                                   abs(child_state - ci_upper[child_node_idx]))
          prob <- 1.0 - min(1.0, dist_from_boundary / ci_width)
        }
      }
    } else {
      # No uncertainty information - use default high probability
      prob <- 1.0
    }
    
    # Skip changes below the threshold or probability
    if(abs(change) < min_change || prob < probability_threshold) {
      next
    }
    
    # Determine event type
    event_type <- NA
    description <- NA
    
    if(change > 0) {
      # Chromosome gain
      if(rel_change >= threshold_ratio && "wgd" %in% event_types) {
        event_type <- "wgd"
        wgd_level <- round(rel_change)
        description <- sprintf("Putative WGD (%.1f-fold increase, ~%dx)", 
                             rel_change, wgd_level)
      } else if(change >= min_change) {
        if("gain" %in% event_types) {
          event_type <- "gain"
          description <- sprintf("Chromosome gain (+%d)", change)
        }
      }
    } else if(change < 0) {
      abs_change <- abs(change)
      
      if(abs_change <= 3 && "fusion" %in% event_types) {
        event_type <- "fusion"
        description <- sprintf("Putative chromosome fusion(s) (-%d)", abs_change)
      } else if("loss" %in% event_types) {
        event_type <- "loss"
        description <- sprintf("Chromosome loss (-%d)", abs_change)
      }
    }
    
    # If we have genome data, refine event classification
    if(!is.na(event_type) && !is.null(genome_data)) {
      # TO BE IMPLEMENTED: Use genome structure data to refine event detection
      # This would use synteny, chromosome length, or other genomic features
    }
    
    # Add event to result if identified
    if(!is.na(event_type)) {
      events <- rbind(events, data.frame(
        Edge = i,
        Parent_Node = parent,
        Child_Node = child,
        Parent_State = parent_state,
        Child_State = child_state,
        Change = change,
        Rel_Change = rel_change,
        Event_Type = event_type,
        Probability = prob,
        Branch_Length = branch_length,
        Description = description,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Create result structure
  result <- list(
    tree = tree,
    reconstruction = reconstruction,
    events = events,
    parameters = list(
      event_types = event_types,
      threshold_ratio = threshold_ratio,
      min_change = min_change,
      method = method,
      probability_threshold = probability_threshold
    )
  )
  
  # Calculate event statistics
  result$statistics <- calculate_event_statistics(result)
  
  # Generate visualizations
  result$plots <- create_event_visualizations(result)
  
  return(result)
}

#' Get state for a tip from reconstruction result
#' 
#' @param tree Phylogenetic tree
#' @param reconstruction Ancestral state reconstruction
#' @param tip_idx Tip index
#' @return Tip state
#' @keywords internal
get_tip_state <- function(tree, reconstruction, tip_idx) {
  # Try different ways the tip state might be stored
  if(!is.null(reconstruction$tip_states)) {
    # Directly from stored tip states
    return(reconstruction$tip_states[tip_idx])
  } else if(!is.null(reconstruction$data)) {
    # From the input data
    tip_name <- tree$tip.label[tip_idx]
    return(reconstruction$data[tip_name])
  } else {
    # Check if there's a chr_counts attribute on the tree
    if(!is.null(tree$chr_counts)) {
      tip_name <- tree$tip.label[tip_idx]
      return(tree$chr_counts[tip_name])
    }
    
    # Last resort - return NA
    warning(sprintf("Could not find state for tip %d", tip_idx))
    return(NA)
  }
}

#' Calculate statistics about detected events
#' 
#' @param event_results Event detection results
#' @return List of event statistics
#' @keywords internal
calculate_event_statistics <- function(event_results) {
  events <- event_results$events
  
  # If no events detected, return empty statistics
  if(nrow(events) == 0) {
    return(list(
      event_count = 0,
      event_types = character(0),
      events_per_type = data.frame(
        Event_Type = character(0),
        Count = integer(0),
        Proportion = numeric(0),
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # Count events by type
  event_counts <- table(events$Event_Type)
  
  # Create event statistics
  stats <- list(
    event_count = nrow(events),
    event_types = unique(events$Event_Type),
    events_per_type = data.frame(
      Event_Type = names(event_counts),
      Count = as.integer(event_counts),
      Proportion = as.numeric(event_counts) / sum(event_counts),
      stringsAsFactors = FALSE
    )
  )
  
  # Calculate chromosome change statistics
  stats$mean_abs_change <- mean(abs(events$Change))
  stats$max_change <- max(events$Change)
  stats$min_change <- min(events$Change)
  
  # Calculate branch length statistics
  stats$mean_branch_length <- mean(events$Branch_Length)
  stats$events_per_branch_length <- nrow(events) / sum(event_results$tree$edge.length)
  
  return(stats)
}

#' Create visualizations of detected events
#' 
#' @param event_results Event detection results
#' @return List of plot objects
#' @keywords internal
create_event_visualizations <- function(event_results) {
  # Initialize plots list
  plots <- list()
  
  # If no events detected, return empty plots
  if(nrow(event_results$events) == 0) {
    return(plots)
  }
  
  # Create event count barplot
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    plots$event_counts <- ggplot2::ggplot(
      event_results$statistics$events_per_type,
      ggplot2::aes(x = Event_Type, y = Count, fill = Event_Type)
    ) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::labs(
        title = "Chromosome Evolution Events",
        x = "Event Type",
        y = "Number of Events"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::scale_fill_brewer(palette = "Set1")
  }
  
  # Create phylogenetic tree with events
  if(requireNamespace("ggtree", quietly = TRUE)) {
    # Get data
    tree <- event_results$tree
    events <- event_results$events
    
    # Create basic tree
    p <- ggtree::ggtree(tree, ladderize = TRUE) + 
      ggtree::geom_tiplab(size = 3)
    
    # Define event symbols and colors
    event_shapes <- c(
      "wgd" = 24,    # Up-pointing triangle
      "fusion" = 21, # Circle
      "fission" = 22, # Square
      "gain" = 25,   # Up-pointing triangle
      "loss" = 23    # Diamond
    )
    
    event_colors <- c(
      "wgd" = "red",
      "fusion" = "blue",
      "fission" = "green",
      "gain" = "purple",
      "loss" = "orange"
    )
    
    # Add symbols for events
    if(nrow(events) > 0) {
      # Create event data for plotting
      event_data <- data.frame(
        x = numeric(nrow(events)),
        y = numeric(nrow(events)),
        event_type = events$Event_Type,
        description = events$Description,
        edge = events$Edge,
        stringsAsFactors = FALSE
      )
      
      # Get edge data for placement
      edge_data <- ggtree::get.edge.data(p$data)
      
      # Populate coordinates
      for(i in 1:nrow(event_data)) {
        edge_idx <- event_data$edge[i]
        edge_row <- which(edge_data$.edge == edge_idx)
        
        if(length(edge_row) > 0) {
          # Place symbol in the middle of the edge
          event_data$x[i] <- (edge_data[edge_row, "x"] + edge_data[edge_row, "xend"]) / 2
          event_data$y[i] <- (edge_data[edge_row, "y"] + edge_data[edge_row, "yend"]) / 2
        }
      }
      
      # Add event symbols to tree
      for(event_type in unique(event_data$event_type)) {
        type_data <- event_data[event_data$event_type == event_type, ]
        
        shape <- event_shapes[event_type]
        color <- event_colors[event_type]
        
        p <- p + ggtree::geom_point(
          data = type_data,
          ggplot2::aes(x = x, y = y, shape = event_type, color = event_type),
          size = 3
        )
      }
      
      # Add shape and color scales
      p <- p + 
        ggplot2::scale_shape_manual(name = "Event Type", values = event_shapes) +
        ggplot2::scale_color_manual(name = "Event Type", values = event_colors) +
        ggplot2::labs(title = "Chromosome Evolution Events") +
        ggplot2::theme(legend.position = "right")
    }
    
    plots$event_tree <- p
  }
  
  return(plots)
}

#===============================================================================
# Event Analysis Functions
#===============================================================================

#' Analyze the distribution of chromosome evolution events
#' 
#' Examines how events are distributed across the phylogeny and tests for
#' significant patterns (clustering, associations with traits, etc.)
#' 
#' @param event_results Event detection results
#' @param trait_data Optional trait data for association testing
#' @param analysis_type Type of analysis: "clustering", "association", "temporal", or "all"
#' @param n_permutations Number of permutations for significance testing
#' @return List with analysis results
#' @export
analyze_event_distribution <- function(event_results,
                                     trait_data = NULL,
                                     analysis_type = "all",
                                     n_permutations = 1000) {
  # Initialize results
  results <- list(
    event_results = event_results,
    analyses = list(),
    statistics = list(),
    plots = list()
  )
  
  # Check if any events were detected
  if(nrow(event_results$events) == 0) {
    warning("No events detected to analyze")
    return(results)
  }
  
  # Determine which analyses to run
  valid_analyses <- c("clustering", "association", "temporal")
  if(analysis_type == "all") {
    analyses_to_run <- valid_analyses
  } else {
    analysis_type <- match.arg(analysis_type, valid_analyses, several.ok = TRUE)
    analyses_to_run <- analysis_type
  }
  
  # Run specified analyses
  if("clustering" %in% analyses_to_run) {
    # Test for phylogenetic clustering of events
    clustering_result <- test_event_clustering(event_results, n_permutations)
    results$analyses$clustering <- clustering_result
    results$statistics$clustering <- clustering_result$statistics
    results$plots$clustering <- clustering_result$plots
  }
  
  if("association" %in% analyses_to_run && !is.null(trait_data)) {
    # Test for association between events and traits
    association_result <- test_trait_association(event_results, trait_data, n_permutations)
    results$analyses$association <- association_result
    results$statistics$association <- association_result$statistics
    results$plots$association <- association_result$plots
  }
  
  if("temporal" %in% analyses_to_run) {
    # Test for temporal patterns in event occurrence
    temporal_result <- test_temporal_patterns(event_results, n_permutations)
    results$analyses$temporal <- temporal_result
    results$statistics$temporal <- temporal_result$statistics
    results$plots$temporal <- temporal_result$plots
  }
  
  # Create summary statistics across all analyses
  results$summary <- summarize_event_analyses(results)
  
  return(results)
}

#' Test for phylogenetic clustering of events
#' 
#' @param event_results Event detection results
#' @param n_permutations Number of permutations for significance testing
#' @return List with clustering analysis results
#' @keywords internal
test_event_clustering <- function(event_results, n_permutations = 1000) {
  # Get event data
  events <- event_results$events
  tree <- event_results$tree
  
  # Initialize result
  result <- list(
    statistics = list(),
    p_values = list(),
    permutations = list(),
    plots = list()
  )
  
  # Calculate mean pairwise distance between events
  event_distances <- calculate_event_distances(events, tree)
  mean_dist <- mean(event_distances)
  result$statistics$mean_pairwise_distance <- mean_dist
  
  # Calculate nearest neighbor distances for events
  nn_distances <- calculate_nearest_neighbor_distances(events, tree)
  mean_nn_dist <- mean(nn_distances)
  result$statistics$mean_nearest_neighbor_distance <- mean_nn_dist
  
  # Permutation test for clustering
  message("Running permutation test for event clustering...")
  
  # Initialize vectors to store permutation results
  perm_mean_dists <- numeric(n_permutations)
  perm_mean_nn_dists <- numeric(n_permutations)
  
  # Run permutations
  for(i in 1:n_permutations) {
    # Create random distribution of events
    random_events <- randomize_events(events, tree)
    
    # Calculate distances for randomized events
    random_dists <- calculate_event_distances(random_events, tree)
    random_nn_dists <- calculate_nearest_neighbor_distances(random_events, tree)
    
    # Store results
    perm_mean_dists[i] <- mean(random_dists)
    perm_mean_nn_dists[i] <- mean(random_nn_dists)
  }
  
  # Store permutation results
  result$permutations$mean_pairwise_distances <- perm_mean_dists
  result$permutations$mean_nearest_neighbor_distances <- perm_mean_nn_dists
  
  # Calculate p-values
  p_mean_dist <- sum(perm_mean_dists <= mean_dist) / n_permutations
  p_mean_nn_dist <- sum(perm_mean_nn_dists <= mean_nn_dist) / n_permutations
  
  result$p_values$mean_pairwise_distance <- p_mean_dist
  result$p_values$mean_nearest_neighbor_distance <- p_mean_nn_dist
  
  # Interpret results
  is_clustered <- (p_mean_dist < 0.05) && (p_mean_nn_dist < 0.05)
  result$interpretation <- ifelse(is_clustered,
                               "Events show significant phylogenetic clustering",
                               "No significant phylogenetic clustering of events detected")
  
  # Create visualization if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Create histogram of permutation results
    perm_df <- data.frame(
      Mean_Distance = perm_mean_dists,
      stringsAsFactors = FALSE
    )
    
    p <- ggplot2::ggplot(perm_df, ggplot2::aes(x = Mean_Distance)) +
      ggplot2::geom_histogram(binwidth = (max(perm_mean_dists) - min(perm_mean_dists))/30, 
                           fill = "skyblue", color = "black") +
      ggplot2::geom_vline(xintercept = mean_dist, color = "red", size = 1) +
      ggplot2::labs(
        title = "Permutation Test for Event Clustering",
        subtitle = paste("P-value =", round(p_mean_dist, 3)),
        x = "Mean Pairwise Distance",
        y = "Frequency",
        caption = "Red line = observed value"
      ) +
      ggplot2::theme_minimal()
    
    result$plots$permutation_histogram <- p
  }
  
  return(result)
}

#' Calculate phylogenetic distances between events
#' 
#' @param events Event data frame
#' @param tree Phylogenetic tree
#' @return Matrix of pairwise distances
#' @keywords internal
calculate_event_distances <- function(events, tree) {
  # Get number of events
  n_events <- nrow(events)
  
  # Initialize distance matrix
  distances <- matrix(NA, n_events, n_events)
  
  # For each pair of events
  for(i in 1:(n_events-1)) {
    for(j in (i+1):n_events) {
      # Get nodes for both events
      node_i <- events$Child_Node[i]
      node_j <- events$Child_Node[j]
      
      # Calculate phylogenetic distance
      path_dist <- ape::dist.nodes(tree)[node_i, node_j]
      
      # Store in matrix
      distances[i, j] <- distances[j, i] <- path_dist
    }
  }
  
  # Convert to vector of distances (excluding diagonal)
  dist_vector <- distances[upper.tri(distances)]
  
  return(dist_vector)
}

#' Calculate nearest neighbor distances for events
#' 
#' @param events Event data frame
#' @param tree Phylogenetic tree
#' @return Vector of nearest neighbor distances
#' @keywords internal
calculate_nearest_neighbor_distances <- function(events, tree) {
  # Get number of events
  n_events <- nrow(events)
  
  # Initialize vector for nearest neighbor distances
  nn_distances <- numeric(n_events)
  
  # Get full distance matrix between nodes
  node_distances <- ape::dist.nodes(tree)
  
  # For each event
  for(i in 1:n_events) {
    # Get node for this event
    node_i <- events$Child_Node[i]
    
    # Get nodes for all other events
    other_nodes <- events$Child_Node[-i]
    
    # Calculate distances to all other events
    distances_to_others <- node_distances[node_i, other_nodes]
    
    # Find minimum distance (nearest neighbor)
    if(length(distances_to_others) > 0) {
      nn_distances[i] <- min(distances_to_others)
    } else {
      nn_distances[i] <- NA
    }
  }
  
  return(nn_distances)
}

#' Randomize event positions for permutation test
#' 
#' @param events Event data frame
#' @param tree Phylogenetic tree
#' @return Randomized event data frame
#' @keywords internal
randomize_events <- function(events, tree) {
  # Get number of events
  n_events <- nrow(events)
  
  # Get all possible edges
  all_edges <- 1:nrow(tree$edge)
  
  # Randomly sample edges for events
  random_edges <- sample(all_edges, n_events)
  
  # Create randomized events
  random_events <- events
  random_events$Edge <- random_edges
  
  # Update node information
  for(i in 1:n_events) {
    edge_idx <- random_events$Edge[i]
    random_events$Parent_Node[i] <- tree$edge[edge_idx, 1]
    random_events$Child_Node[i] <- tree$edge[edge_idx, 2]
    random_events$Branch_Length[i] <- tree$edge.length[edge_idx]
  }
  
  return(random_events)
}

#' Test for association between events and traits
#' 
#' @param event_results Event detection results
#' @param trait_data Data frame of trait values
#' @param n_permutations Number of permutations for significance testing
#' @return List with association analysis results
#' @keywords internal
test_trait_association <- function(event_results, trait_data, n_permutations = 1000) {
  # Initialize result
  result <- list(
    statistics = list(),
    p_values = list(),
    plots = list()
  )
  
  # TO BE IMPLEMENTED: Association testing between events and traits
  
  return(result)
}

#' Test for temporal patterns in event occurrence
#' 
#' @param event_results Event detection results
#' @param n_permutations Number of permutations for significance testing
#' @return List with temporal analysis results
#' @keywords internal
test_temporal_patterns <- function(event_results, n_permutations = 1000) {
  # Initialize result
  result <- list(
    statistics = list(),
    p_values = list(),
    plots = list()
  )
  
  # TO BE IMPLEMENTED: Analysis of temporal patterns in event occurrence
  
  return(result)
}

#' Summarize event analysis results
#' 
#' @param analysis_results Event analysis results
#' @return Summary of analysis results
#' @keywords internal
summarize_event_analyses <- function(analysis_results) {
  # Initialize summary
  summary <- list(
    significant_findings = character(0)
  )
  
  # Check clustering results
  if(!is.null(analysis_results$analyses$clustering)) {
    clustering <- analysis_results$analyses$clustering
    
    if(clustering$p_values$mean_pairwise_distance < 0.05) {
      summary$significant_findings <- c(
        summary$significant_findings,
        sprintf("Events show significant phylogenetic clustering (p = %.3f)",
              clustering$p_values$mean_pairwise_distance)
      )
    }
  }
  
  # Check association results
  if(!is.null(analysis_results$analyses$association)) {
    # TO BE IMPLEMENTED: Extract significant associations
  }
  
  # Check temporal results
  if(!is.null(analysis_results$analyses$temporal)) {
    # TO BE IMPLEMENTED: Extract significant temporal patterns
  }
  
  # If no significant findings, note that
  if(length(summary$significant_findings) == 0) {
    summary$significant_findings <- "No significant patterns detected in analyses"
  }
  
  return(summary)
}

#===============================================================================
# Workflow Functions
#===============================================================================

#' Run complete event detection and analysis workflow
#' 
#' @param tree Phylogenetic tree
#' @param reconstruction Ancestral state reconstruction results
#' @param output_dir Directory to save results
#' @param event_types Types of events to detect
#' @param threshold_ratio Threshold ratio for WGD detection
#' @param trait_data Optional trait data for association testing
#' @param generate_report Whether to generate HTML report
#' @param run_analyses Whether to analyze event distributions
#' @return List with complete workflow results
#' @export
run_event_analysis_workflow <- function(tree,
                                      reconstruction,
                                      output_dir = "event_analysis",
                                      event_types = "all",
                                      threshold_ratio = 1.8,
                                      trait_data = NULL,
                                      generate_report = TRUE,
                                      run_analyses = TRUE) {
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 1: Detect events
  message("Step 1: Detecting chromosome number change events...")
  events <- detect_chromosome_events(
    tree = tree,
    reconstruction = reconstruction,
    event_types = event_types,
    threshold_ratio = threshold_ratio
  )
  
  # Save event data
  write.csv(events$events, file.path(output_dir, "detected_events.csv"), row.names = FALSE)
  
  # Save event visualizations
  if(length(events$plots) > 0) {
    plots_dir <- file.path(output_dir, "plots")
    if(!dir.exists(plots_dir)) {
      dir.create(plots_dir)
    }
    
    for(plot_name in names(events$plots)) {
      if(requireNamespace("ggplot2", quietly = TRUE)) {
        ggplot2::ggsave(
          filename = file.path(plots_dir, paste0("event_", plot_name, ".pdf")),
          plot = events$plots[[plot_name]],
          width = 10,
          height = 8
        )
      }
    }
  }
  
  # Initialize workflow results
  workflow_results <- list(
    event_detection = events
  )
  
  # Step 2: Run analyses if requested
  if(run_analyses) {
    message("Step 2: Analyzing event distributions...")
    analyses <- analyze_event_distribution(
      event_results = events,
      trait_data = trait_data
    )
    
    # Save analysis results
    saveRDS(analyses, file.path(output_dir, "event_analysis_results.rds"))
    
    # Save analysis plots
    if(length(analyses$plots) > 0) {
      plots_dir <- file.path(output_dir, "plots")
      if(!dir.exists(plots_dir)) {
        dir.create(plots_dir)
      }
      
      for(category in names(analyses$plots)) {
        category_plots <- analyses$plots[[category]]
        
        for(plot_name in names(category_plots)) {
          if(requireNamespace("ggplot2", quietly = TRUE)) {
            ggplot2::ggsave(
              filename = file.path(plots_dir, paste0(category, "_", plot_name, ".pdf")),
              plot = category_plots[[plot_name]],
              width = 10,
              height = 8
            )
          }
        }
      }
    }
    
    workflow_results$event_analysis <- analyses
  }
  
  # Step 3: Generate report if requested
  if(generate_report && requireNamespace("rmarkdown", quietly = TRUE)) {
    message("Step 3: Generating HTML report...")
    
    # Create report data
    report_data <- list(
      events = events,
      analyses = if(run_analyses) workflow_results$event_analysis else NULL
    )
    
    # Save report data
    saveRDS(report_data, file.path(output_dir, "report_data.rds"))
    
    # Generate report
    # TO BE IMPLEMENTED: Create Rmd template and render it
  }
  
  message("Event analysis workflow complete.")
  message(paste0("Results saved to: ", output_dir))
  
  return(workflow_results)
}
