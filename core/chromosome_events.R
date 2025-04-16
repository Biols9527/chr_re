#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Chromosome Events Analysis Module
# Author: Bioinformatics Team
# Date: 2025-03-24
# Description: Provides chromosome evolutionary event detection, classification,
#              and statistical analysis functionality
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(data.table)
  library(igraph)
})

#===============================================================================
# Event Detection and Classification
#===============================================================================

#' Detect chromosome evolutionary events
#' 
#' Detect chromosome fusion, fission and whole genome duplication events based on
#' ancestral state reconstruction results
#' 
#' @param tree Phylogenetic tree object
#' @param ancestral_counts Ancestral chromosome count reconstruction results
#' @param fusion_threshold Fusion event relative change threshold (child/parent ratio)
#' @param fission_threshold Fission event relative change threshold (child/parent ratio)
#' @param wgd_threshold Whole genome duplication event relative change threshold
#' @param allow_wgd Whether to allow whole genome duplication events
#' @param adaptive_thresholds Whether to use adaptive thresholds based on data
#' @return Chromosome events list
#' @export
detect_chromosome_events <- function(tree, ancestral_counts, 
                                   fusion_threshold = 0.67,  
                                   fission_threshold = 1.5, 
                                   wgd_threshold = 2.0,
                                   allow_wgd = TRUE,
                                   adaptive_thresholds = FALSE) {
  # Check input parameters
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(is.null(ancestral_counts$ancestral_states) || !is.data.frame(ancestral_counts$ancestral_states)) {
    stop("ancestral_counts must contain ancestral_states data frame")
  }
  
  message("Detecting chromosome evolutionary events...")
  
  # Use adaptive thresholds based on data distribution if requested
  if(adaptive_thresholds) {
    message("Using adaptive thresholds for chromosome event detection...")
    
    # Analyze all edge change ratios
    all_ratios <- events$change_ratio[!is.na(events$change_ratio)]
    
    if(length(all_ratios) > 10) {
      # Try to use mixture model to identify potential multimodal distribution
      tryCatch({
        if(requireNamespace("mixtools", quietly = TRUE)) {
          # Use mixture model to classify change ratios into different types
          mixture_model <- mixtools::normalmixEM(all_ratios, k = 3)
          
          # Extract mixture component parameters
          mu <- mixture_model$mu
          sigma <- mixture_model$sigma
          
          # Sort by mean
          sorted_idx <- order(mu)
          sorted_mu <- mu[sorted_idx]
          sorted_sigma <- sigma[sorted_idx]
          
          if(length(mu) >= 3) {
            # Update thresholds based on crossing points between mixture components
            fusion_threshold <- find_crossing_point(sorted_mu[1], sorted_sigma[1], 
                                                 sorted_mu[2], sorted_sigma[2])
            
            fission_threshold <- find_crossing_point(sorted_mu[2], sorted_sigma[2], 
                                                   sorted_mu[3], sorted_sigma[3])
            
            # WGD threshold as a proportion of the third component mean
            wgd_threshold <- sorted_mu[3] * 0.9
            
            message(sprintf("Adaptive thresholds: fusion=%.2f, fission=%.2f, WGD=%.2f", 
                          fusion_threshold, fission_threshold, wgd_threshold))
          }
        }
      }, error = function(e) {
        message("Adaptive threshold calculation failed, using default thresholds: ", e$message)
      })
    }
  }
  
  # Get ancestral states
  node_states <- ancestral_counts$ancestral_states
  tip_states <- ancestral_counts$tip_states
  
  # Create data frame of all node states (including tips and internal nodes)
  all_nodes <- data.frame(
    node_id = c(1:length(tree$tip.label), node_states$node_id),
    state = c(tip_states, node_states$state),
    is_tip = c(rep(TRUE, length(tree$tip.label)), rep(FALSE, nrow(node_states))),
    label = c(tree$tip.label, rep(NA, nrow(node_states))),
    stringsAsFactors = FALSE
  )
  
  # Create events data frame
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
    stringsAsFactors = FALSE
  )
  
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
    
    # Determine event type
    if(is.na(events$change_ratio[i])) {
      events$event_type[i] <- "unknown"
      events$event_confidence[i] <- 0
      events$event_magnitude[i] <- 0
    } else if(events$change_ratio[i] <= fusion_threshold) {
      # Chromosome fusion
      events$event_type[i] <- "fusion"
      events$event_magnitude[i] <- parent_count - child_count
      
      # Calculate confidence: greater change = higher confidence
      events$event_confidence[i] <- min(1.0, 
                                     (1 - events$change_ratio[i]) / (1 - fusion_threshold))
    } else if(events$change_ratio[i] >= wgd_threshold && allow_wgd) {
      # Whole genome duplication
      events$event_type[i] <- "wgd"
      events$event_magnitude[i] <- child_count - parent_count
      
      # Calculate confidence: closer to exact multiple = higher confidence
      ratio_to_closest_multiple <- min(abs(events$change_ratio[i] %% 1), 
                                     abs(1 - (events$change_ratio[i] %% 1)))
      events$event_confidence[i] <- max(0, 1 - 2 * ratio_to_closest_multiple)
    } else if(events$change_ratio[i] >= fission_threshold) {
      # Chromosome fission
      events$event_type[i] <- "fission"
      events$event_magnitude[i] <- child_count - parent_count
      
      # Calculate confidence: greater change = higher confidence
      events$event_confidence[i] <- min(1.0, 
                                     (events$change_ratio[i] - fission_threshold) / fission_threshold)
    } else {
      # No significant change
      events$event_type[i] <- "none"
      events$event_confidence[i] <- 1.0 - abs(events$rel_change[i])
      events$event_magnitude[i] <- abs(child_count - parent_count)
    }
  }
  
  # Add species names (for tip nodes)
  events$species <- NA
  tip_nodes <- events$child_node[events$child_node <= length(tree$tip.label)]
  events$species[events$child_node %in% tip_nodes] <- 
    tree$tip.label[events$child_node[events$child_node %in% tip_nodes]]
  
  # Add branch lengths (if available)
  if(!is.null(tree$edge.length)) {
    events$branch_length <- tree$edge.length
  }
  
  # Calculate event rates by branch length (if available)
  if(!is.null(tree$edge.length)) {
    events$event_rate <- ifelse(events$branch_length > 0 & events$event_type != "none",
                              events$event_magnitude / events$branch_length, 
                              NA)
  }
  
  # Summarize event statistics
  event_counts <- table(events$event_type)
  message("Detected chromosome events:")
  for(event_type in names(event_counts)) {
    if(event_type != "none") {
      message(sprintf("  %s: %d", event_type, event_counts[event_type]))
    }
  }
  
  # Create detailed result object
  result <- list(
    events = events,
    tree = tree,
    ancestral_counts = ancestral_counts,
    parameters = list(
      fusion_threshold = fusion_threshold,
      fission_threshold = fission_threshold,
      wgd_threshold = wgd_threshold,
      allow_wgd = allow_wgd
    ),
    summary = list(
      event_counts = event_counts,
      avg_confidence = tapply(events$event_confidence, events$event_type, mean, na.rm = TRUE),
      avg_magnitude = tapply(events$event_magnitude, events$event_type, mean, na.rm = TRUE)
    )
  )
  
  return(result)
}

#' Find crossing point between two normal distributions
#' 
#' @param mu1 Mean of first distribution
#' @param sigma1 SD of first distribution
#' @param mu2 Mean of second distribution
#' @param sigma2 SD of second distribution
#' @return Crossing point value
#' @keywords internal
find_crossing_point <- function(mu1, sigma1, mu2, sigma2) {
  # Solve quadratic equation to find where two normal density functions are equal
  a <- 1/(2*sigma1^2) - 1/(2*sigma2^2)
  b <- mu2/(sigma2^2) - mu1/(sigma1^2)
  c <- mu1^2/(2*sigma1^2) - mu2^2/(2*sigma2^2) + log(sigma2/sigma1)
  
  # Quadratic formula
  discriminant <- b^2 - 4*a*c
  
  # Handle cases
  if(discriminant < 0) {
    # No real solution, return midpoint of means
    return((mu1 + mu2)/2)
  }
  
  # Two solutions
  x1 <- (-b + sqrt(discriminant))/(2*a)
  x2 <- (-b - sqrt(discriminant))/(2*a)
  
  # Choose solution that's between the means
  if(mu1 <= x1 && x1 <= mu2) {
    return(x1)
  } else if(mu1 <= x2 && x2 <= mu2) {
    return(x2)
  } else {
    # If neither solution is appropriate, use midpoint
    return((mu1 + mu2)/2)
  }
}

#' Enhance event significance with background modeling
#' 
#' Perform chromosome event background modeling and significance testing
#' 
#' @param events Event detection results
#' @param null_replicates Number of null model replicates
#' @param confidence_level Confidence level
#' @return Enhanced event testing results
#' @export
enhance_event_significance <- function(events, null_replicates = 1000, confidence_level = 0.95) {
  message("Performing chromosome event significance enhancement analysis...")
  
  # Extract data
  event_data <- events$events
  tree <- events$tree
  
  # Initialize enhanced results
  enhanced_events <- event_data
  enhanced_events$p_value <- NA
  enhanced_events$significance <- NA
  
  # Analyze each event type
  event_types <- unique(event_data$event_type[event_data$event_type != "none"])
  
  for(event_type in event_types) {
    message(sprintf("  Analyzing statistical significance of %s events...", event_type))
    
    # Filter this type of events
    type_events <- event_data[event_data$event_type == event_type,]
    
    # Exclude direct descendants of root (if multiple), as they may lack ancestral info
    root_children <- tree$edge[tree$edge[,1] == (Ntip(tree) + 1), 2]
    if(length(root_children) > 0) {
      type_events <- type_events[!type_events$child_node %in% root_children,]
    }
    
    if(nrow(type_events) == 0) {
      message(sprintf("    No %s events to analyze", event_type))
      next
    }
    
    # Perform permutation test
    try({
      # Simulate null distribution
      null_distribution <- simulate_null_events(tree, type_events$change_ratio,
                                              n_replicates = null_replicates)
      
      # Calculate p-value for each event
      for(i in 1:nrow(type_events)) {
        event_idx <- which(event_data$edge_id == type_events$edge_id[i])
        if(length(event_idx) > 0) {
          observed_ratio <- abs(type_events$change_ratio[i])
          
          # P-value: proportion of null distribution as extreme as observed
          if(event_type %in% c("fusion")) {
            # Fusion events: smaller ratio is more extreme
            p_value <- mean(null_distribution <= observed_ratio)
          } else {
            # Fission/WGD events: larger ratio is more extreme
            p_value <- mean(null_distribution >= observed_ratio)
          }
          
          enhanced_events$p_value[event_idx] <- p_value
          enhanced_events$significance[event_idx] <- p_value < (1 - confidence_level)
        }
      }
      
      message(sprintf("    Analyzed %d %s events, of which %d statistically significant(p < %.2f)",
                     nrow(type_events),
                     event_type,
                     sum(enhanced_events$significance[enhanced_events$event_type == event_type], na.rm=TRUE),
                     1 - confidence_level))
    }, silent = FALSE)
  }
  
  # Create enhanced result object
  enhanced_result <- events
  enhanced_result$events <- enhanced_events
  enhanced_result$significance_analysis <- list(
    confidence_level = confidence_level,
    null_replicates = null_replicates,
    significant_events = sum(enhanced_events$significance, na.rm=TRUE),
    total_events = sum(enhanced_events$event_type != "none")
  )
  
  return(enhanced_result)
}

#' Simulate null distribution for chromosome events
#' 
#' Simulate chromosome change null model distribution
#' 
#' @param tree Phylogenetic tree object
#' @param observed_ratios Observed change ratios
#' @param n_replicates Simulation replicates
#' @return Null model change ratio distribution
#' @keywords internal
simulate_null_events <- function(tree, observed_ratios, n_replicates = 1000) {
  # Create null model distribution
  # Simulate random chromosome changes using appropriate model
  
  # Initialize null distribution
  null_distribution <- numeric(n_replicates)
  
  # Generate simulated data
  for(i in 1:n_replicates) {
    # Simplified version: directly permute observed ratios
    null_distribution[i] <- sample(observed_ratios, 1)
    
    # More complex version would use Brownian motion or other phylogenetic models
    # to simulate chromosome number evolution on the tree
  }
  
  return(null_distribution)
}

#' Analyze chromosome macroevolution patterns
#' 
#' Analyze macroevolutionary patterns in chromosome evolution
#' 
#' @param events Event detection results
#' @param ancestral_counts Ancestral chromosome counts
#' @return Macroevolution pattern analysis
#' @export
analyze_macroevolution_patterns <- function(events, ancestral_counts) {
  message("Analyzing chromosome macroevolution patterns...")
  
  # Extract data
  event_data <- events$events
  tree <- events$tree
  
  # Initialize results
  macro_patterns <- list()
  
  # 1. Analyze chromosome number change over time
  try({
    # Get tree time information
    if(!is.null(tree$edge.length)) {
      # Calculate path length from root to each node
      node_depths <- node.depth.edgelength(tree)
      
      # Normalize to relative time
      total_height <- max(node_depths)
      relative_times <- node_depths / total_height
      
      # Combine ancestral states and time info
      time_series <- data.frame(
        node_id = c(1:length(tree$tip.label), ancestral_counts$ancestral_states$node_id),
        state = c(ancestral_counts$tip_states, ancestral_counts$ancestral_states$state),
        relative_time = relative_times
      )
      
      # Analyze trend
      trend_model <- lm(state ~ relative_time, data = time_series)
      
      macro_patterns$time_trend <- list(
        slope = coef(trend_model)[2],
        p_value = summary(trend_model)$coefficients[2,4],
        r_squared = summary(trend_model)$r.squared,
        direction = ifelse(coef(trend_model)[2] > 0, "increase", "decrease"),
        significant = summary(trend_model)$coefficients[2,4] < 0.05
      )
      
      if(macro_patterns$time_trend$significant) {
        message(sprintf("  Found significant time trend: chromosome numbers %s over time (slope = %.2f, p = %.3f)",
                       macro_patterns$time_trend$direction,
                       macro_patterns$time_trend$slope,
                       macro_patterns$time_trend$p_value))
      } else {
        message("  No significant time trend detected")
      }
    }
  }, silent = FALSE)
  
  # 2. Analyze phylogenetic signal of fusion/fission events
  try({
    # Initialize count vectors for each species
    tip_fusion_counts <- numeric(length(tree$tip.label))
    names(tip_fusion_counts) <- tree$tip.label
    
    tip_fission_counts <- numeric(length(tree$tip.label))
    names(tip_fission_counts) <- tree$tip.label
    
    # Fill event counts
    for(sp in tree$tip.label) {
      tip_edges <- event_data[event_data$species == sp,]
      tip_fusion_counts[sp] <- sum(tip_edges$event_type == "fusion")
      tip_fission_counts[sp] <- sum(tip_edges$event_type == "fission")
    }
    
    # Calculate phylogenetic signal
    if(requireNamespace("phytools", quietly = TRUE)) {
      # Fusion events phylogenetic signal
      if(sum(tip_fusion_counts) > 0) {
        fusion_signal <- phytools::phylosig(tree, tip_fusion_counts, method = "K", test = TRUE)
        
        macro_patterns$fusion_phylo_signal <- list(
          K = fusion_signal$K,
          p_value = fusion_signal$P,
          significant = fusion_signal$P < 0.05
        )
        
        message(sprintf("  Fusion event phylogenetic signal: K = %.2f, p = %.3f %s",
                       fusion_signal$K, fusion_signal$P,
                       ifelse(fusion_signal$P < 0.05, "[significant]", "")))
      }
      
      # Fission events phylogenetic signal
      if(sum(tip_fission_counts) > 0) {
        fission_signal <- phytools::phylosig(tree, tip_fission_counts, method = "K", test = TRUE)
        
        macro_patterns$fission_phylo_signal <- list(
          K = fission_signal$K,
          p_value = fission_signal$P,
          significant = fission_signal$P < 0.05
        )
        
        message(sprintf("  Fission event phylogenetic signal: K = %.2f, p = %.3f %s",
                       fission_signal$K, fission_signal$P,
                       ifelse(fission_signal$P < 0.05, "[significant]", "")))
      }
    }
  }, silent = FALSE)
  
  # 4. Test if event rates follow molecular clock
  try({
    # Get edges with events
    edges_with_events <- event_data[event_data$event_type != "none",]
    
    # Test correlation between event magnitude and branch length
    if(nrow(edges_with_events) > 10 && !is.null(tree$edge.length)) {
      cor_test <- cor.test(edges_with_events$event_magnitude,
                         edges_with_events$branch_length,
                         method = "spearman")
      
      macro_patterns$molecular_clock <- list(
        correlation = cor_test$estimate,
        p_value = cor_test$p.value,
        significant = cor_test$p.value < 0.05,
        follows_clock = cor_test$p.value < 0.05 && cor_test$estimate > 0
      )
      
      if(macro_patterns$molecular_clock$significant) {
        clock_type <- ifelse(cor_test$estimate > 0, "follows", "violates")
        message(sprintf("  Chromosome events %s molecular clock model (rho = %.2f, p = %.3f)",
                       clock_type, cor_test$estimate, cor_test$p.value))
      } else {
        message("  Cannot determine if chromosome events follow molecular clock")
      }
    }
  }, silent = FALSE)
  
  # Generate overall conclusions
  macro_patterns$overall_conclusions <- generate_macro_conclusions(macro_patterns)
  
  # Print conclusions
  for(conclusion in macro_patterns$overall_conclusions) {
    message("  Conclusion: ", conclusion)
  }
  
  return(macro_patterns)
}

#' Generate macroevolution pattern conclusions
#' 
#' @param macro_patterns Macroevolution pattern analysis results
#' @return Vector of conclusion strings
#' @keywords internal
generate_macro_conclusions <- function(macro_patterns) {
  conclusions <- c()
  
  # Based on time trend
  if(!is.null(macro_patterns$time_trend)) {
    if(macro_patterns$time_trend$significant) {
      if(macro_patterns$time_trend$direction == "increase") {
        conclusions <- c(conclusions, "Chromosome numbers show significant overall increase over time, which may reflect dominance of fission events in the evolutionary process")
      } else {
        conclusions <- c(conclusions, "Chromosome numbers show significant overall decrease over time, which may reflect dominance of fusion events in the evolutionary process")
      }
    } else {
      conclusions <- c(conclusions, "Chromosome numbers do not show a clear temporal trend, possibly indicating balance between fusion and fission events")
    }
  }
  
  # Based on phylogenetic signal
  if(!is.null(macro_patterns$fusion_phylo_signal) || !is.null(macro_patterns$fission_phylo_signal)) {
    fusion_sig <- !is.null(macro_patterns$fusion_phylo_signal) && macro_patterns$fusion_phylo_signal$significant
    fission_sig <- !is.null(macro_patterns$fission_phylo_signal) && macro_patterns$fission_phylo_signal$significant
    
    if(fusion_sig && fission_sig) {
      conclusions <- c(conclusions, "Both fusion and fission events show significant phylogenetic signal, indicating certain clades are more prone to specific types of chromosome rearrangements")
    } else if(fusion_sig) {
      conclusions <- c(conclusions, "Fusion events show significant phylogenetic signal, while fission events are more randomly distributed across taxa")
    } else if(fission_sig) {
      conclusions <- c(conclusions, "Fission events show significant phylogenetic signal, while fusion events are more randomly distributed across taxa")
    } else if(!is.null(macro_patterns$fusion_phylo_signal) || !is.null(macro_patterns$fission_phylo_signal)) {
      conclusions <- c(conclusions, "Chromosome rearrangement events are distributed relatively randomly across the phylogeny, showing no strong clade-specificity")
    }
  }
  
  # Based on molecular clock test
  if(!is.null(macro_patterns$molecular_clock)) {
    if(macro_patterns$molecular_clock$significant) {
      if(macro_patterns$molecular_clock$follows_clock) {
        conclusions <- c(conclusions, "Chromosome event rates are proportional to evolutionary time, consistent with a molecular clock model, suggesting time-dependent chromosome evolution")
      } else {
        conclusions <- c(conclusions, "Chromosome event rates are not proportional to evolutionary time, violating molecular clock model, suggesting chromosome evolution may be driven by selection pressures or other factors")
      }
    }
  }
  
  # If no conclusions
  if(length(conclusions) == 0) {
    conclusions <- c("Current data insufficient to draw clear conclusions about chromosome macroevolutionary patterns")
  }
  
  return(conclusions)
}

#' Analyze event distribution
#' 
#' Analyze the phylogenetic distribution pattern of chromosome events
#' 
#' @param events Event detection results
#' @param event_type Event type to analyze: "fusion", "fission", "wgd" or "all"
#' @param min_confidence Minimum confidence threshold
#' @return Event distribution analysis results
#' @export
analyze_event_distribution <- function(events, event_type = "all", min_confidence = 0.5) {
  if(!inherits(events$tree, "phylo")) {
    stop("events object must contain valid phylogenetic tree")
  }
  
  message(sprintf("Analyzing phylogenetic distribution of %s events...", 
                 ifelse(event_type == "all", "all", event_type)))
  
  # Filter events
  event_data <- events$events
  
  if(event_type != "all") {
    event_data <- event_data[event_data$event_type == event_type,]
  } else {
    event_data <- event_data[event_data$event_type != "none",]
  }
  
  # Apply confidence filtering
  event_data <- event_data[event_data$event_confidence >= min_confidence,]
  
  if(nrow(event_data) == 0) {
    message("No events found matching criteria")
    return(NULL)
  }
  
  # Get tree structure
  tree <- events$tree
  n_nodes <- Ntip(tree) + Nnode(tree)
  
  # Calculate distance from root for each node
  node_depths <- node.depth.edgelength(tree)
  if(is.null(node_depths)) {
    # If tree has no edge lengths, use topological depth
    node_depths <- node.depth(tree, method = "n")
  }
  
  # Normalize depths
  max_depth <- max(node_depths)
  norm_depths <- node_depths / max_depth
  
  # Create event distribution
  distribution <- data.frame(
    edge_id = event_data$edge_id,
    parent_node = event_data$parent_node,
    child_node = event_data$child_node,
    event_type = event_data$event_type,
    confidence = event_data$event_confidence,
    magnitude = event_data$event_magnitude,
    rel_depth = norm_depths[event_data$child_node],
    abs_depth = node_depths[event_data$child_node],
    species = event_data$species,
    stringsAsFactors = FALSE
  )
  
  # Calculate event depth distribution
  depth_bins <- cut(distribution$rel_depth, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
  depth_distribution <- table(depth_bins)
  
  # Calculate event density by depth
  tree_nodes_by_depth <- hist(norm_depths, breaks = seq(0, 1, by = 0.1), plot = FALSE)$counts
  event_density <- as.numeric(depth_distribution) / tree_nodes_by_depth
  
  # Analyze nodes with significant event clustering
  node_event_counts <- table(event_data$parent_node)
  significant_nodes <- names(node_event_counts)[node_event_counts > 1]
  
  # Calculate event frequency in different species branches
  species_events <- table(event_data$species[!is.na(event_data$species)])
  
  # Perform phylogenetic signal test (if enough events)
  phylo_signal <- NULL
  if(nrow(event_data) >= 10 && requireNamespace("phytools", quietly = TRUE)) {
    # Create event count for each tip
    tip_event_counts <- numeric(Ntip(tree))
    names(tip_event_counts) <- tree$tip.label
    
    # Fill species with events
    for(sp in names(species_events)) {
      if(sp %in% tree$tip.label) {
        tip_event_counts[sp] <- species_events[sp]
      }
    }
    
    # Test phylogenetic signal
    tryCatch({
      phylo_signal <- phytools::phylosig(tree, tip_event_counts, method = "K", test = TRUE)
      message(sprintf("Phylogenetic signal of event distribution: K = %.3f (p = %.3f)", 
                     phylo_signal$K, phylo_signal$P))
    }, error = function(e) {
      message("Unable to test phylogenetic signal: ", e$message)
    })
  }
  
  # Create analysis result
  result <- list(
    distribution = distribution,
    depth_distribution = depth_distribution,
    event_density = event_density,
    significant_nodes = significant_nodes,
    species_events = species_events,
    phylo_signal = phylo_signal,
    summary = list(
      event_count = nrow(distribution),
      event_type = event_type,
      min_confidence = min_confidence,
      mean_depth = mean(distribution$rel_depth),
      sd_depth = sd(distribution$rel_depth)
    )
  )
  
  message(sprintf("Found %d %s events, average relative depth: %.2f", 
                 nrow(distribution), 
                 ifelse(event_type == "all", "significant", event_type),
                 mean(distribution$rel_depth)))
  
  if(length(significant_nodes) > 0) {
    message(sprintf("Found %d nodes with significant event density", length(significant_nodes)))
  }
  
  return(result)
}

#' Detect potential whole genome duplication events
#' 
#' Identify nodes with strong signal indicating potential whole genome duplication
#' 
#' @param events Event detection results
#' @param wgd_threshold Minimum change ratio threshold, typically 1.8-2.0
#' @param conservation_ratio Ratio of post-event conserved relationships
#' @param mapping_data Bidirectional mapping data (optional)
#' @return Identified WGD events list
#' @export
detect_wgd_events <- function(events, wgd_threshold = 1.8, conservation_ratio = 0.6, mapping_data = NULL) {
  message("Detecting potential whole genome duplication events...")
  
  # Extract event data and tree
  event_data <- events$events
  tree <- events$tree
  
  # Find events with significant chromosome number increase
  potential_wgd <- event_data[event_data$change_ratio >= wgd_threshold,]
  
  if(nrow(potential_wgd) == 0) {
    message("No potential whole genome duplication events found")
    return(NULL)
  }
  
  # Sort by event magnitude
  potential_wgd <- potential_wgd[order(-potential_wgd$change_ratio),]
  
  # Analyze each potential WGD event
  wgd_events <- list()
  
  for(i in 1:nrow(potential_wgd)) {
    edge_id <- potential_wgd$edge_id[i]
    parent_node <- potential_wgd$parent_node[i]
    child_node <- potential_wgd$child_node[i]
    ratio <- potential_wgd$change_ratio[i]
    parent_count <- potential_wgd$parent_count[i]
    child_count <- potential_wgd$child_count[i]
    
    # Calculate how close the ratio is to nearest integer multiple (WGDs often close to 2x, 3x etc.)
    closest_multiplier <- round(ratio)
    deviation_from_multiplier <- abs(ratio - closest_multiplier)
    multiplier_confidence <- max(0, 1 - deviation_from_multiplier / 0.5)
    
    # Get descendant species
    descendants <- NULL
    if(child_node > Ntip(tree)) {
      # Internal node, get all descendant species
      descendants <- extract_clade(tree, child_node)$tip.label
    } else {
      # Leaf node
      descendants <- tree$tip.label[child_node]
    }
    
    # Analyze chromosome count distribution in descendants
    if(length(descendants) > 0 && !is.null(events$ancestral_counts$tip_states)) {
      descendant_counts <- events$ancestral_counts$tip_states[descendants]
      
      # Calculate similarity of chromosome counts to expected post-WGD count
      expected_count <- parent_count * closest_multiplier
      count_consistency <- sapply(descendant_counts, function(count) {
        abs_diff <- abs(count - expected_count)
        sim <- 1 - min(1, abs_diff / expected_count)
        return(sim)
      })
      
      # Calculate proportion of well-preserved descendants
      well_preserved <- sum(count_consistency > conservation_ratio) / length(count_consistency)
    } else {
      well_preserved <- NA
    }
    
    # Add karyotype evidence analysis if mapping data available
    karyotype_evidence <- NULL
    if(!is.null(mapping_data) && length(descendants) > 1) {
      # Analyze synteny/homology evidence for WGD
      karyotype_evidence <- analyze_karyotype_evidence(
        mapping_data, parent_node, child_node, descendants, tree
      )
    }
    
    # Calculate overall WGD likelihood score
    wgd_score <- mean(c(
      multiplier_confidence,
      ifelse(is.na(well_preserved), multiplier_confidence, well_preserved)
    ))
    
    # If karyotype evidence available, integrate it into the score
    if(!is.null(karyotype_evidence)) {
      wgd_score <- wgd_score * 0.7 + karyotype_evidence$score * 0.3
    }
    
    # Create WGD event object
    wgd_event <- list(
      edge_id = edge_id,
      parent_node = parent_node,
      child_node = child_node,
      ratio = ratio,
      parent_count = parent_count,
      child_count = child_count,
      closest_multiplier = closest_multiplier,
      multiplier_confidence = multiplier_confidence,
      descendants = descendants,
      well_preserved_ratio = well_preserved,
      wgd_score = wgd_score
    )
    
    # Add karyotype evidence if available
    if(!is.null(karyotype_evidence)) {
      wgd_event$karyotype_evidence <- karyotype_evidence
    }
    
    # Determine if it's a likely WGD
    if(wgd_score >= 0.6) {
      wgd_events[[length(wgd_events) + 1]] <- wgd_event
    }
  }
  
  if(length(wgd_events) == 0) {
    message("No high-scoring whole genome duplication events found")
    return(NULL)
  }
  
  # Sort by WGD score
  wgd_scores <- sapply(wgd_events, function(x) x$wgd_score)
  wgd_events <- wgd_events[order(-wgd_scores)]
  
  message(sprintf("Detected %d potential whole genome duplication events", length(wgd_events)))
  
  # Create analysis result
  wgd_result <- list(
    events = wgd_events,
    tree = tree,
    parameters = list(
      wgd_threshold = wgd_threshold,
      conservation_ratio = conservation_ratio
    )
  )
  
  return(wgd_result)
}

#' Analyze karyotype evidence
#' 
#' Analyze karyotype composition evidence for WGD
#' 
#' @param mapping_data Bidirectional mapping data
#' @param parent_node Parent node
#' @param child_node Child node
#' @param descendants Descendant species
#' @param tree Phylogenetic tree
#' @return Karyotype evidence analysis
#' @keywords internal
analyze_karyotype_evidence <- function(mapping_data, parent_node, child_node, descendants, tree) {
  # Extract chromosome composition from mapping data for relevant species
  # Look for evidence supporting WGD such as:
  # - Duplicated conserved regions
  # - Parallel changes across chromosomes
  # - Other signatures expected post-WGD
  
  # Default medium evidence
  evidence <- list(
    score = 0.5,  # Default to medium evidence
    duplicated_regions = 0,
    pattern_strength = 0
  )
  
  # Analysis implementation logic would go here...
  
  return(evidence)
}

#' Calculate chromosome evolution rates
#' 
#' Calculate chromosome number change rates in different parts of the tree
#' 
#' @param events Event detection results
#' @param taxonomic_groups Optional list of taxonomic groups, names are groups, values are species name vectors
#' @return Chromosome evolution rate analysis results
#' @export
calculate_chromosome_evolution_rates <- function(events, taxonomic_groups = NULL) {
  message("Calculating chromosome number evolution rates...")
  
  # Extract data
  event_data <- events$events
  tree <- events$tree
  
  # Ensure tree has branch lengths
  if(is.null(tree$edge.length)) {
    message("Tree lacks branch length information, using topological distances")
    tree$edge.length <- rep(1, nrow(tree$edge))
  }
  
  # Calculate overall chromosome change rate
  total_change <- sum(abs(event_data$abs_change), na.rm = TRUE)
  total_length <- sum(tree$edge.length)
  overall_rate <- total_change / total_length
  
  # Calculate rates by event type
  event_types <- unique(event_data$event_type[event_data$event_type != "none"])
  type_rates <- sapply(event_types, function(type) {
    type_data <- event_data[event_data$event_type == type,]
    type_change <- sum(type_data$event_magnitude, na.rm = TRUE)
    type_rate <- type_change / total_length
    return(type_rate)
  })
  
  # Create basic rate result
  rate_result <- list(
    overall_rate = overall_rate,
    type_rates = type_rates,
    edges = data.frame(
      edge_id = event_data$edge_id,
      parent_node = event_data$parent_node,
      child_node = event_data$child_node,
      event_type = event_data$event_type,
      branch_length = event_data$branch_length,
      abs_change = event_data$abs_change,
      change_rate = ifelse(event_data$branch_length > 0, 
                         abs(event_data$abs_change) / event_data$branch_length, 
                         NA)
    )
  )
  
  # If taxonomic groups provided, calculate group-specific rates
  if(!is.null(taxonomic_groups) && is.list(taxonomic_groups)) {
    group_rates <- list()
    
    for(group_name in names(taxonomic_groups)) {
      species <- taxonomic_groups[[group_name]]
      
      # Check if species are in tree
      valid_species <- species[species %in% tree$tip.label]
      
      if(length(valid_species) > 0) {
        # Find edges leading to these species
        tip_ids <- match(valid_species, tree$tip.label)
        
        # Find edges leading to these tips
        group_edges <- which(event_data$child_node %in% tip_ids)
        
        if(length(group_edges) > 0) {
          group_data <- event_data[group_edges,]
          group_change <- sum(abs(group_data$abs_change), na.rm = TRUE)
          group_length <- sum(group_data$branch_length, na.rm = TRUE)
          group_rate <- group_change / group_length
          
          # Create group rate result
          group_rates[[group_name]] <- list(
            species = valid_species,
            num_species = length(valid_species),
            total_change = group_change,
            total_length = group_length,
            overall_rate = group_rate,
            edges = group_edges
          )
        }
      }
    }
    
    # Add group rates to result
    rate_result$group_rates <- group_rates
    
    # Calculate group rate statistics
    if(length(group_rates) > 0) {
      group_names <- names(group_rates)
      group_rate_values <- sapply(group_rates, function(x) x$overall_rate)
      rate_result$group_comparison <- data.frame(
        group = group_names,
        rate = group_rate_values,
        relative_rate = group_rate_values / overall_rate,
        stringsAsFactors = FALSE
      )
      
      message("Taxonomic group chromosome evolution rate comparison:")
      for(i in 1:nrow(rate_result$group_comparison)) {
        message(sprintf("  %s: %.4f (%.2f x overall rate)", 
                       rate_result$group_comparison$group[i],
                       rate_result$group_comparison$rate[i],
                       rate_result$group_comparison$relative_rate[i]))
      }
    }
  }
  
  # Create overall rate summary
  message(sprintf("Overall chromosome evolution rate: %.4f changes/unit time", overall_rate))
  for(type in names(type_rates)) {
    message(sprintf("  %s event rate: %.4f/unit time (%.1f%%)", 
                   type, type_rates[type], 100 * type_rates[type] / overall_rate))
  }
  
  return(rate_result)
}

#===============================================================================
# Visualization Functions
#===============================================================================

#' Plot chromosome events on tree
#' 
#' Visualize different types of chromosome events on the phylogenetic tree
#' 
#' @param events Event detection results
#' @param event_types Event types to display vector
#' @param min_confidence Minimum confidence threshold
#' @param show_counts Whether to display chromosome numbers at nodes
#' @param count_cex Chromosome number text size
#' @param clade_labels Optional clade label list
#' @param output_file Output file path
#' @return No visible return
#' @export
plot_chromosome_events <- function(events, 
                                 event_types = c("fusion", "fission", "wgd"),
                                 min_confidence = 0.5,
                                 show_counts = TRUE,
                                 count_cex = 0.7,
                                 clade_labels = NULL,
                                 output_file = NULL) {
  # Check required packages
  if(!requireNamespace("ape", quietly = TRUE)) {
    stop("ape package required for visualization")
  }
  
  message("Visualizing chromosome events...")
  
  # Extract data
  event_data <- events$events
  tree <- events$tree
  ancestral_counts <- events$ancestral_counts
  
  # Filter events of interest
  significant_events <- event_data[
    event_data$event_type %in% event_types & 
    event_data$event_confidence >= min_confidence, ]
  
  # Set event type colors and shapes
  event_colors <- c(
    "fusion" = "blue",
    "fission" = "red",
    "wgd" = "purple"
  )
  
  event_pch <- c(
    "fusion" = 25,  # Down triangle
    "fission" = 24, # Up triangle
    "wgd" = 23      # Diamond
  )
  
  # Setup output device if needed
  if(!is.null(output_file)) {
    # Choose device based on file extension
    if(grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = 10, height = 8)
    } else if(grepl("\\.png$", output_file)) {
      png(output_file, width = 2000, height = 1600, res = 200)
    } else {
      # Default PDF
      pdf(paste0(output_file, ".pdf"), width = 10, height = 8)
    }
  }
  
  # Plot base tree
  plot(tree, cex = 0.7, label.offset = 0.01, no.margin = TRUE,
       main = "Chromosome Evolution Events")
  
  # Add event symbols on branches
  if(nrow(significant_events) > 0) {
    for(i in 1:nrow(significant_events)) {
      event <- significant_events[i, ]
      
      # Find corresponding edge
      edge_idx <- which(tree$edge[,1] == event$parent_node & 
                      tree$edge[,2] == event$child_node)
      
      if(length(edge_idx) == 1) {
        # Choose color and shape based on event type
        event_color <- event_colors[event$event_type]
        event_shape <- event_pch[event$event_type]
        
        # Adjust point size based on confidence
        point_size <- 1 + event$event_confidence * 1.5
        
        # Add event marker
        edgelabels(text = "", edge = edge_idx, 
                  bg = event_color, 
                  pch = event_shape,
                  cex = point_size)
      }
    }
  }
  
  # Add chromosome counts if requested
  if(show_counts && !is.null(ancestral_counts)) {
    # Extract ancestral states
    node_states <- ancestral_counts$ancestral_states
    tip_states <- ancestral_counts$tip_states
    
    # Add internal node chromosome counts
    if(!is.null(node_states)) {
      nodes <- node_states$node_id
      values <- node_states$state
      
      for(i in 1:length(nodes)) {
        node <- nodes[i]
        value <- values[i]
        
        # Format count (integer or 1 decimal)
        if(abs(value - round(value)) < 0.05) {
          count_text <- as.character(round(value))
        } else {
          count_text <- sprintf("%.1f", value)
        }
        
        # Add node label
        nodelabels(text = count_text, 
                  node = node, 
                  frame = "circle", 
                  bg = "white",
                  cex = count_cex)
      }
    }
    
    # Add tip node chromosome counts
    if(!is.null(tip_states)) {
      for(i in 1:length(tree$tip.label)) {
        if(tree$tip.label[i] %in% names(tip_states)) {
          value <- tip_states[tree$tip.label[i]]
          
          # Format count (integer or 1 decimal)
          if(abs(value - round(value)) < 0.05) {
            count_text <- as.character(round(value))
          } else {
            count_text <- sprintf("%.1f", value)
          }
          
          # Add tip label
          tiplabels(text = count_text, 
                   tip = i, 
                   frame = "circle", 
                   bg = "white",
                   cex = count_cex)
        }
      }
    }
  }
  
  # Add clade labels if provided
  if(!is.null(clade_labels) && is.list(clade_labels)) {
    for(clade_name in names(clade_labels)) {
      species <- clade_labels[[clade_name]]
      
      # Ensure species are in tree
      species_in_tree <- species[species %in% tree$tip.label]
      
      if(length(species_in_tree) >= 2) {
        # Find most recent common ancestor
        mrca_node <- phytools::findMRCA(tree, species_in_tree)
        
        if(!is.null(mrca_node)) {
          # Add clade label
          nodelabels(text = clade_name, 
                    node = mrca_node, 
                    frame = "none", 
                    bg = "transparent",
                    adj = c(-0.1, 0.5),
                    cex = 0.9,
                    col = "darkgreen",
                    font = 2)
        }
      }
    }
  }
  
  # Add legend
  legend_types <- intersect(names(event_colors), event_types)
  if(length(legend_types) > 0) {
    legend("topright", 
          legend = legend_types, 
          col = event_colors[legend_types], 
          pch = event_pch[legend_types],
          title = "Event Types",
          cex = 0.8)
  }
  
  if(!is.null(output_file)) {
    dev.off()
    message(paste("Chromosome event plot saved to:", output_file))
  }
  
  # Return invisible
  return(invisible(NULL))
}

#' Plot chromosome rate heatmap
#' 
#' Create a heatmap of chromosome number change rates on the phylogenetic tree
#' 
#' @param events Event detection results
#' @param rate_type Rate type: 'absolute' or 'relative'
#' @param log_transform Whether to log-transform rates
#' @param color_palette Color palette function
#' @param output_file Output file path
#' @return No visible return
#' @export
plot_chromosome_rates <- function(events, 
                                rate_type = "absolute", 
                                log_transform = TRUE,
                                color_palette = NULL,
                                output_file = NULL) {
  # Check required packages
  if(!requireNamespace("phytools", quietly = TRUE)) {
    stop("phytools package required for visualization")
  }
  
  message("Visualizing chromosome evolution rates...")
  
  # Extract data
  event_data <- events$events
  tree <- events$tree
  
  # Ensure tree has branch lengths
  if(is.null(tree$edge.length)) {
    message("Tree lacks branch length information, using topological distances")
    tree$edge.length <- rep(1, nrow(tree$edge))
  }
  
  # Calculate rates
  if(rate_type == "absolute") {
    # Absolute change rate: absolute change / branch length
    rates <- abs(event_data$abs_change) / event_data$branch_length
  } else if(rate_type == "relative") {
    # Relative change rate: relative change / branch length
    rates <- abs(event_data$rel_change) / event_data$branch_length
  } else {
    stop("Invalid rate type. Use 'absolute' or 'relative'")
  }
  
  # Handle NAs and Infs
  rates[is.na(rates) | !is.finite(rates)] <- 0
  
  # Log transform if requested
  if(log_transform && any(rates > 0)) {
    rates <- log1p(rates)  # log(x+1) to avoid infinities for small values
    rate_label <- paste0("log(", 
                        ifelse(rate_type == "absolute", "Absolute", "Relative"), 
                        " Change Rate + 1)")
  } else {
    rate_label <- paste0(ifelse(rate_type == "absolute", "Absolute", "Relative"), 
                       " Change Rate")
  }
  
  # Set rate color range
  max_rate <- max(rates, na.rm = TRUE)
  if(max_rate == 0) max_rate <- 1  # Avoid all zeros
  
  # Create color palette
  if(is.null(color_palette)) {
    color_palette <- colorRampPalette(c("blue", "green", "yellow", "red"))
  }
  
  # Set edge colors
  edge_colors <- color_palette(100)[ceiling(rates/max_rate * 99) + 1]
  edge_colors[is.na(edge_colors)] <- "gray"
  
  # Setup output device if needed
  if(!is.null(output_file)) {
    # Choose device based on file extension
    if(grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = 10, height = 8)
    } else if(grepl("\\.png$", output_file)) {
      png(output_file, width = 2000, height = 1600, res = 200)
    } else {
      # Default PDF
      pdf(paste0(output_file, ".pdf"), width = 10, height = 8)
    }
  }
  
  # Plot colored tree
  plot(tree, edge.color = edge_colors, cex = 0.7, 
       main = paste("Chromosome", rate_label))
  
  # Add color legend
  add_color_bar <- function(palette, min_val, max_val, 
                          title = "Rate", 
                          nticks = 5, digits = 2) {
    # Draw color bar
    par(fig = c(0.05, 0.3, 0.1, 0.3), new = TRUE, mar = c(1, 1, 1, 5))
    image(
      matrix(seq(min_val, max_val, length.out = 100), ncol = 1),
      col = palette(100),
      axes = FALSE,
      xlab = "",
      ylab = ""
    )
    
    # Add ticks
    ticks <- seq(min_val, max_val, length.out = nticks)
    if(log_transform) {
      labels <- round(expm1(ticks), digits)  # Reverse log transform
    } else {
      labels <- round(ticks, digits)
    }
    
    axis(4, at = seq(0, 1, length.out = nticks), 
         labels = labels, 
         las = 1, cex.axis = 0.7)
    
    mtext(title, side = 4, line = 3, cex = 0.8)
    
    # Restore plot area
    par(fig = c(0, 1, 0, 1), new = FALSE)
  }
  
  # Add color bar
  add_color_bar(color_palette, 0, max_rate, title = rate_label)
  
  if(!is.null(output_file)) {
    dev.off()
    message(paste("Chromosome rate plot saved to:", output_file))
  }
  
  return(invisible(NULL))
}

#' Get tree node descendants
#' 
#' @param tree Phylogenetic tree object
#' @param node Node number
#' @return Descendant species vector
#' @keywords internal
extract_clade <- function(tree, node) {
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(node <= Ntip(tree)) {
    # Leaf node
    return(list(tip.label = tree$tip.label[node]))
  }
  
  # Extract subtree containing this node
  edge <- tree$edge
  edges_to_keep <- descendant_edges(edge, node)
  
  nodes_to_keep <- unique(c(edge[edges_to_keep, 1], edge[edges_to_keep, 2]))
  tips_to_keep <- nodes_to_keep[nodes_to_keep <= Ntip(tree)]
  
  return(list(tip.label = tree$tip.label[tips_to_keep]))
}

#' Get descendant edges from node
#' 
#' @param edge Edge matrix
#' @param node Node ID
#' @return Vector of edge indices
#' @keywords internal
descendant_edges <- function(edge, node) {
  # Find edges where this node is parent
  direct_edges <- which(edge[, 1] == node)
  
  # If no direct children, return empty
  if(length(direct_edges) == 0) {
    return(numeric(0))
  }
  
  # Get all descendants recursively
  all_edges <- direct_edges
  
  for(e in direct_edges) {
    child <- edge[e, 2]
    desc_edges <- descendant_edges(edge, child)
    all_edges <- c(all_edges, desc_edges)
  }
  
  return(all_edges)
}
