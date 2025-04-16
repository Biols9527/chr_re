#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Karyotype Evolution Module
# Author: Bioinformatics Team
# Date: 2023-06-10
# Description: Methods for analyzing karyotype evolution patterns, including
#              fusion/fission events detection, polyploidy inference, and
#              statistical testing of evolutionary hypotheses
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(geiger)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})

#===============================================================================
# Rearrangement Detection Functions
#===============================================================================

#' Detect chromosome rearrangement events along phylogeny branches
#' 
#' Analyzes changes in chromosome numbers between ancestral and descendant nodes
#' to infer fusion, fission, and polyploidy events
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts for extant species
#' @param ancestral_states Named vector with reconstructed ancestral states
#' @param polyploidy_threshold Numeric threshold for detecting polyploidy events (e.g., 1.8 = 80% increase)
#' @param consider_branch_lengths Whether to account for branch lengths in detection
#' @return List containing detected events for each branch
#' @export
detect_chromosome_events <- function(tree, 
                                   chr_counts, 
                                   ancestral_states,
                                   polyploidy_threshold = 1.8,
                                   consider_branch_lengths = TRUE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Combine tip and ancestral state data
  all_states <- c(chr_counts[tree$tip.label], ancestral_states)
  node_ids <- c(1:length(tree$tip.label), 
                (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
  names(all_states) <- node_ids
  
  # Initialize results
  results <- list(
    tree = tree,
    events = list(),
    summary = list(
      fusion_count = 0,
      fission_count = 0,
      polyploidy_count = 0,
      no_change_count = 0,
      total_branches = nrow(tree$edge)
    ),
    branch_data = data.frame(
      parent = tree$edge[, 1],
      child = tree$edge[, 2],
      parent_chr = NA,
      child_chr = NA,
      delta = NA,
      event_type = NA,
      event_count = NA,
      stringsAsFactors = FALSE
    )
  )
  
  # Process each branch in the tree
  for(i in 1:nrow(tree$edge)) {
    parent_node <- tree$edge[i, 1]
    child_node <- tree$edge[i, 2]
    
    # Get chromosome numbers
    parent_chr <- all_states[as.character(parent_node)]
    child_chr <- all_states[as.character(child_node)]
    
    # Update branch data
    results$branch_data$parent_chr[i] <- parent_chr
    results$branch_data$child_chr[i] <- child_chr
    
    # Skip branches with missing data
    if(is.na(parent_chr) || is.na(child_chr)) {
      results$branch_data$event_type[i] <- "unknown"
      next
    }
    
    # Calculate change in chromosome number
    delta <- child_chr - parent_chr
    results$branch_data$delta[i] <- delta
    
    # Determine event type
    if(delta == 0) {
      # No change
      event_type <- "no_change"
      event_count <- 0
      results$summary$no_change_count <- results$summary$no_change_count + 1
    } else if(child_chr / parent_chr >= polyploidy_threshold) {
      # Polyploidy event
      event_type <- "polyploidy"
      
      # Estimate level of polyploidy (diploid, triploid, etc.)
      ploidy_level <- round(child_chr / parent_chr)
      event_count <- 1  # Count as a single event regardless of level
      
      results$summary$polyploidy_count <- results$summary$polyploidy_count + 1
      results$events[[i]] <- list(
        branch = i,
        parent = parent_node,
        child = child_node,
        event_type = event_type,
        ploidy_level = ploidy_level
      )
    } else if(delta > 0) {
      # Chromosome fission (increase in number)
      event_type <- "fission"
      event_count <- delta  # One event per added chromosome
      
      results$summary$fission_count <- results$summary$fission_count + delta
      results$events[[i]] <- list(
        branch = i,
        parent = parent_node,
        child = child_node,
        event_type = event_type,
        event_count = event_count
      )
    } else {  # delta < 0
      # Chromosome fusion (decrease in number)
      event_type <- "fusion"
      event_count <- abs(delta)  # One event per lost chromosome
      
      results$summary$fusion_count <- results$summary$fusion_count + abs(delta)
      results$events[[i]] <- list(
        branch = i,
        parent = parent_node,
        child = child_node,
        event_type = event_type,
        event_count = event_count
      )
    }
    
    # Update branch data
    results$branch_data$event_type[i] <- event_type
    results$branch_data$event_count[i] <- event_count
  }
  
  # Add summary stats
  total_changes <- results$summary$fusion_count + 
                  results$summary$fission_count + 
                  results$summary$polyploidy_count
  
  results$summary$total_events <- total_changes
  results$summary$events_per_branch <- total_changes / nrow(tree$edge)
  
  # Create event proportions
  if(total_changes > 0) {
    results$summary$fusion_proportion <- results$summary$fusion_count / total_changes
    results$summary$fission_proportion <- results$summary$fission_count / total_changes
    results$summary$polyploidy_proportion <- results$summary$polyploidy_count / total_changes
  } else {
    results$summary$fusion_proportion <- 0
    results$summary$fission_proportion <- 0
    results$summary$polyploidy_proportion <- 0
  }
  
  # Add class for the results object
  class(results) <- c("chr_events", class(results))
  
  return(results)
}

#' Test directional bias in chromosome evolution
#' 
#' Tests whether there is a significant bias toward fusions or fissions
#' in chromosome number evolution
#' 
#' @param events Results from detect_chromosome_events function
#' @param method Testing method: "binomial", "permutation", "mcmc"
#' @param n_permutations Number of permutations for permutation test
#' @return List with test results and statistics
#' @export
test_fusion_fission_bias <- function(events,
                                   method = "binomial",
                                   n_permutations = 1000) {
  # Validate input
  if(!inherits(events, "chr_events")) {
    stop("events must be a chr_events object from detect_chromosome_events")
  }
  
  # Get counts
  fusion_count <- events$summary$fusion_count
  fission_count <- events$summary$fission_count
  total_count <- fusion_count + fission_count
  
  # Skip test if no events
  if(total_count == 0) {
    return(list(
      method = method,
      fusion_count = 0,
      fission_count = 0,
      p_value = NA,
      significant = FALSE,
      direction = "none"
    ))
  }
  
  # Calculate observed proportion
  observed_prop <- fusion_count / total_count
  
  # Initialize results
  results <- list(
    method = method,
    fusion_count = fusion_count,
    fission_count = fission_count,
    total_count = total_count,
    observed_proportion = observed_prop
  )
  
  # Conduct test based on method
  if(method == "binomial") {
    # Binomial test against null of equal probability
    binom_test <- binom.test(fusion_count, total_count, p = 0.5)
    
    results$p_value <- binom_test$p.value
    results$significant <- binom_test$p.value < 0.05
    results$test_result <- binom_test
    
    # Determine direction of bias
    if(results$significant) {
      results$direction <- ifelse(observed_prop > 0.5, "fusion", "fission")
    } else {
      results$direction <- "none"
    }
  } else if(method == "permutation") {
    # Permutation test
    # Get branch-specific events
    branch_data <- events$branch_data
    
    # Only include branches with fusion or fission events
    branch_data <- branch_data[branch_data$event_type %in% c("fusion", "fission"), ]
    
    if(nrow(branch_data) == 0) {
      return(list(
        method = method,
        fusion_count = 0,
        fission_count = 0,
        p_value = NA,
        significant = FALSE,
        direction = "none"
      ))
    }
    
    # Calculate observed test statistic (difference in counts)
    observed_diff <- fusion_count - fission_count
    
    # Run permutations
    perm_diffs <- numeric(n_permutations)
    
    for(i in 1:n_permutations) {
      # Randomly permute event types
      perm_types <- sample(
        branch_data$event_type, 
        size = nrow(branch_data),
        replace = FALSE
      )
      
      # Calculate counts
      perm_fusion <- sum(branch_data$event_count[perm_types == "fusion"])
      perm_fission <- sum(branch_data$event_count[perm_types == "fission"])
      
      # Calculate difference
      perm_diffs[i] <- perm_fusion - perm_fission
    }
    
    # Calculate p-value (two-sided test)
    if(observed_diff > 0) {
      p_value <- sum(perm_diffs >= observed_diff) / n_permutations
    } else {
      p_value <- sum(perm_diffs <= observed_diff) / n_permutations
    }
    
    results$p_value <- p_value
    results$significant <- p_value < 0.05
    results$observed_diff <- observed_diff
    results$perm_diffs <- perm_diffs
    
    # Determine direction of bias
    if(results$significant) {
      results$direction <- ifelse(observed_diff > 0, "fusion", "fission")
    } else {
      results$direction <- "none"
    }
  } else {
    stop(paste("Unsupported method:", method))
  }
  
  return(results)
}

#' Calculate rates of chromosome evolution
#' 
#' Calculates overall and branch-specific rates of chromosome evolution,
#' accounting for branch lengths
#' 
#' @param tree Phylogenetic tree
#' @param events Results from detect_chromosome_events function
#' @param rate_metric Which metric to use for rate calculation: "events_per_time", "events_per_branch"
#' @param log_transform Whether to log-transform branch lengths
#' @return List with calculated rates
#' @export
calculate_chromosome_rates <- function(tree, 
                                     events,
                                     rate_metric = "events_per_time",
                                     log_transform = FALSE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!inherits(events, "chr_events")) {
    stop("events must be a chr_events object from detect_chromosome_events")
  }
  
  # Check if the tree has branch lengths
  if(is.null(tree$edge.length) && rate_metric == "events_per_time") {
    warning("Tree has no branch lengths. Using 'events_per_branch' metric instead.")
    rate_metric <- "events_per_branch"
  }
  
  # Get branch data
  branch_data <- events$branch_data
  
  # Add branch lengths
  branch_data$branch_length <- tree$edge.length
  
  if(log_transform) {
    branch_data$branch_length <- log(branch_data$branch_length + 0.001)
  }
  
  # Calculate branch-specific rates
  if(rate_metric == "events_per_time") {
    # Handle zero-length branches to avoid division by zero
    branch_data$branch_length[branch_data$branch_length <= 0] <- min(branch_data$branch_length[branch_data$branch_length > 0]) / 10
    
    # Calculate rate as events per unit branch length
    branch_data$rate <- branch_data$event_count / branch_data$branch_length
  } else {
    # Simply use event count as the rate
    branch_data$rate <- branch_data$event_count
  }
  
  # Calculate overall rates by event type
  rate_summary <- list()
  
  for(event_type in c("fusion", "fission", "polyploidy")) {
    type_branches <- branch_data$event_type == event_type
    
    if(sum(type_branches) > 0) {
      if(rate_metric == "events_per_time") {
        total_events <- sum(branch_data$event_count[type_branches])
        total_time <- sum(branch_data$branch_length[type_branches])
        rate_summary[[event_type]] <- total_events / total_time
      } else {
        total_events <- sum(branch_data$event_count[type_branches])
        n_branches <- sum(type_branches)
        rate_summary[[event_type]] <- total_events / n_branches
      }
    } else {
      rate_summary[[event_type]] <- 0
    }
  }
  
  # Calculate overall rate
  if(rate_metric == "events_per_time") {
    total_events <- sum(branch_data$event_count, na.rm = TRUE)
    total_time <- sum(branch_data$branch_length)
    rate_summary$overall <- total_events / total_time
  } else {
    total_events <- sum(branch_data$event_count, na.rm = TRUE)
    n_branches <- nrow(branch_data)
    rate_summary$overall <- total_events / n_branches
  }
  
  # Prepare results
  results <- list(
    branch_rates = branch_data,
    rate_summary = rate_summary,
    rate_metric = rate_metric,
    log_transform = log_transform
  )
  
  class(results) <- c("chr_rates", class(results))
  
  return(results)
}

#' Test for significant differences in rates among clades
#' 
#' Tests whether chromosome evolution rates differ significantly among
#' specified clades in the phylogeny
#' 
#' @param tree Phylogenetic tree
#' @param rates Results from calculate_chromosome_rates function
#' @param clades List of clades, each specified as a vector of tip labels
#' @param method Testing method: "anova", "permutation"
#' @param n_permutations Number of permutations for permutation test
#' @return List with test results
#' @export
test_rate_differences <- function(tree,
                                rates,
                                clades,
                                method = "permutation",
                                n_permutations = 1000) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!inherits(rates, "chr_rates")) {
    stop("rates must be a chr_rates object from calculate_chromosome_rates")
  }
  
  # Check that clades list is properly formatted
  if(!is.list(clades) || length(clades) < 2) {
    stop("clades must be a list of at least two vectors of tip labels")
  }
  
  # Extract branch rate data
  branch_data <- rates$branch_rates
  
  # Assign branches to clades
  branch_data$clade <- NA
  
  for(i in 1:length(clades)) {
    clade_tips <- clades[[i]]
    
    # Validate that tips exist in tree
    if(!all(clade_tips %in% tree$tip.label)) {
      missing_tips <- clade_tips[!clade_tips %in% tree$tip.label]
      warning(paste("Some tips not found in tree:", paste(missing_tips, collapse = ", ")))
      clade_tips <- clade_tips[clade_tips %in% tree$tip.label]
    }
    
    # Skip empty clades
    if(length(clade_tips) == 0) next
    
    # Find MRCA of clade
    if(length(clade_tips) == 1) {
      mrca_node <- which(tree$tip.label == clade_tips)
    } else {
      mrca_node <- phytools::findMRCA(tree, clade_tips)
    }
    
    # Find all descendants of MRCA
    descendants <- phytools::getDescendants(tree, mrca_node)
    
    # Assign all branches descending from MRCA to this clade
    branch_idx <- which(branch_data$child %in% descendants)
    branch_data$clade[branch_idx] <- i
  }
  
  # Remove branches not assigned to any clade
  branch_data <- branch_data[!is.na(branch_data$clade), ]
  
  # Check if we have enough data
  if(nrow(branch_data) < 3) {
    warning("Not enough branches assigned to clades for meaningful testing")
    
    return(list(
      method = method,
      p_value = NA,
      significant = FALSE
    ))
  }
  
  # Calculate mean rate per clade
  clade_rates <- tapply(branch_data$rate, branch_data$clade, mean, na.rm = TRUE)
  
  # Test based on method
  if(method == "anova") {
    # ANOVA test
    anova_result <- aov(rate ~ as.factor(clade), data = branch_data)
    anova_summary <- summary(anova_result)
    
    p_value <- anova_summary[[1]]["as.factor(clade)", "Pr(>F)"]
    
    results <- list(
      method = "anova",
      p_value = p_value,
      significant = p_value < 0.05,
      clade_rates = clade_rates,
      anova_result = anova_summary
    )
  } else if(method == "permutation") {
    # Calculate observed F statistic
    observed_f <- summary(aov(rate ~ as.factor(clade), data = branch_data))[[1]]["as.factor(clade)", "F value"]
    
    # Permutation test
    perm_f <- numeric(n_permutations)
    
    for(i in 1:n_permutations) {
      # Permute clade assignments
      perm_data <- branch_data
      perm_data$clade <- sample(branch_data$clade)
      
      # Calculate F statistic
      perm_f[i] <- summary(aov(rate ~ as.factor(clade), data = perm_data))[[1]]["as.factor(clade)", "F value"]
    }
    
    # Calculate p-value
    p_value <- sum(perm_f >= observed_f) / n_permutations
    
    results <- list(
      method = "permutation",
      p_value = p_value,
      significant = p_value < 0.05,
      clade_rates = clade_rates,
      observed_f = observed_f,
      perm_f = perm_f
    )
  } else {
    stop(paste("Unsupported method:", method))
  }
  
  # Add ranks of clades by rate
  results$clade_rank <- rank(-clade_rates)
  names(results$clade_rank) <- names(clade_rates)
  
  # Determine which clade has highest and lowest rates
  results$highest_rate_clade <- as.numeric(names(which.max(clade_rates)))
  results$lowest_rate_clade <- as.numeric(names(which.min(clade_rates)))
  
  # Add ratio of highest to lowest rate
  results$max_min_ratio <- max(clade_rates) / min(clade_rates)
  
  return(results)
}

#===============================================================================
# Visualization Functions
#===============================================================================

#' Plot chromosome events mapped onto a phylogenetic tree
#' 
#' Visualizes fusion, fission, and polyploidy events detected along branches
#' 
#' @param tree Phylogenetic tree
#' @param events Results from detect_chromosome_events function
#' @param show_counts Whether to show chromosome counts at nodes
#' @param show_tips Whether to display tip labels
#' @param node_size Size of nodes relative to number of events
#' @param event_colors Custom color scheme for events
#' @param layout Tree layout: "rectangular", "circular", "fan"
#' @return ggplot object with chromosome event visualization
#' @export
plot_chromosome_events <- function(tree, 
                                 events, 
                                 show_counts = TRUE, 
                                 show_tips = TRUE,
                                 node_size = 3,
                                 event_colors = NULL,
                                 layout = "rectangular") {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!inherits(events, "chr_events")) {
    stop("events must be a chr_events object from detect_chromosome_events")
  }
  
  if(!requireNamespace("ggplot2", quietly = TRUE) || 
     !requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggplot2 and ggtree packages are required for visualization")
  }
  
  # Extract branch data
  branch_data <- events$branch_data
  
  # Create data frame for edge mapping
  edge_data <- data.frame(
    parent = branch_data$parent,
    node = branch_data$child,
    event_type = branch_data$event_type,
    event_count = branch_data$event_count,
    parent_chr = branch_data$parent_chr,
    child_chr = branch_data$child_chr,
    delta = branch_data$delta,
    stringsAsFactors = FALSE
  )
  
  # Replace NA event types with "none" for clean plotting
  edge_data$event_type[is.na(edge_data$event_type)] <- "none"
  
  # Extract chromosome counts for all nodes
  node_data <- data.frame(
    node = c(1:length(tree$tip.label), 
             (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)),
    is_tip = c(rep(TRUE, length(tree$tip.label)), rep(FALSE, tree$Nnode)),
    stringsAsFactors = FALSE
  )
  
  # Add chromosome counts to node data
  node_data$chr_count <- NA
  
  # First add tip values
  for(i in 1:length(tree$tip.label)) {
    node_data$chr_count[i] <- branch_data$child_chr[branch_data$child == i]
  }
  
  # Then add internal node values
  for(i in (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)) {
    node_data$chr_count[node_data$node == i] <- branch_data$parent_chr[branch_data$parent == i][1]
  }
  
  # Default color scheme if not provided
  if(is.null(event_colors)) {
    event_colors <- c(
      "fusion" = "blue",
      "fission" = "red",
      "polyploidy" = "purple",
      "no_change" = "gray80",
      "unknown" = "gray50",
      "none" = "gray80"
    )
  }
  
  # Create tree base with appropriate layout
  if(layout == "circular") {
    p <- ggtree::ggtree(tree, layout = "circular")
  } else if(layout == "fan") {
    p <- ggtree::ggtree(tree, layout = "fan")
  } else {
    p <- ggtree::ggtree(tree)
  }
  
  # Add edge coloring by event type
  p <- p %<+% edge_data + 
    aes(color = event_type) +
    scale_color_manual(values = event_colors, name = "Event Type")
  
  # Adjust edge width based on event count (optional)
  if(any(!is.na(edge_data$event_count) & edge_data$event_count > 0)) {
    p <- p + ggtree::geom_tree(aes(size = ifelse(is.na(event_count) | event_count == 0, 1, 1 + log(event_count)))) +
      scale_size_continuous(name = "Event Count", range = c(0.5, 3), guide = "none")
  }
  
  # Add node points and chromosome counts
  p <- p %<+% node_data
  
  # Add node points with size based on importance
  p <- p + geom_nodepoint(aes(subset = !is_tip), size = node_size, color = "black", alpha = 0.7)
  
  # Show chromosome counts at nodes if requested
  if(show_counts) {
    p <- p + geom_text(aes(label = chr_count, subset = !is_tip), 
                      hjust = -0.3, vjust = 0.5, size = 3)
  }
  
  # Show tip labels if requested
  if(show_tips) {
    p <- p + geom_tiplab(size = 3, hjust = -0.05)
  }
  
  # Add title and theme
  p <- p + ggtree::labs(title = "Chromosome Evolution Events") +
    theme(legend.position = "right")
  
  return(p)
}

#' Create heatmap visualization of chromosome number changes
#' 
#' Generates a heatmap showing chromosome number changes across the phylogeny
#' 
#' @param tree Phylogenetic tree
#' @param events Results from detect_chromosome_events function
#' @param ancestral_states Ancestral chromosome counts
#' @param tip_values Tip chromosome counts (if different from events data)
#' @param color_scheme Color scheme for heatmap
#' @param show_labels Whether to show node labels
#' @return ggplot object with chromosome number heatmap
#' @export
plot_chromosome_heatmap <- function(tree,
                                  events,
                                  ancestral_states = NULL,
                                  tip_values = NULL,
                                  color_scheme = "viridis",
                                  show_labels = TRUE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!inherits(events, "chr_events")) {
    stop("events must be a chr_events object from detect_chromosome_events")
  }
  
  if(!requireNamespace("ggplot2", quietly = TRUE) || 
     !requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggplot2 and ggtree packages are required for visualization")
  }
  
  # Extract branch data
  branch_data <- events$branch_data
  
  # Initialize chromosome data
  if(is.null(ancestral_states)) {
    # Use event data to fill chromosome counts
    all_nodes <- c(1:length(tree$tip.label), 
                  (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
    chr_counts <- numeric(length(all_nodes))
    names(chr_counts) <- all_nodes
    
    # Fill tip values
    if(is.null(tip_values)) {
      for(i in 1:length(tree$tip.label)) {
        chr_counts[i] <- branch_data$child_chr[branch_data$child == i]
      }
    } else {
      chr_counts[1:length(tree$tip.label)] <- tip_values[tree$tip.label]
    }
    
    # Fill internal nodes
    for(i in (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)) {
      matches <- branch_data$parent == i
      if(any(matches)) {
        chr_counts[as.character(i)] <- branch_data$parent_chr[matches][1]
      }
    }
  } else {
    # Use provided ancestral states
    chr_counts <- ancestral_states
  }
  
  # Create data frame for mapping
  node_data <- data.frame(
    node = as.numeric(names(chr_counts)),
    chr_num = as.numeric(chr_counts),
    stringsAsFactors = FALSE
  )
  
  # Set up color palette
  if(color_scheme == "viridis") {
    color_pal <- viridis::viridis(100)
  } else if(color_scheme == "magma") {
    color_pal <- viridis::magma(100)
  } else if(color_scheme == "inferno") {
    color_pal <- viridis::inferno(100)
  } else if(color_scheme == "plasma") {
    color_pal <- viridis::plasma(100)
  } else {
    color_pal <- viridis::viridis(100)  # Default
  }
  
  # Create tree with mapped chromosome numbers
  p <- ggtree::ggtree(tree, aes(color = chr_num)) %<+% node_data +
    scale_color_gradientn(colors = color_pal, name = "Chromosome\nNumber") +
    ggtree::geom_nodepoint(aes(color = chr_num), size = 3) +
    labs(title = "Chromosome Number Evolution") +
    theme_tree2()
  
  # Add tip labels if requested
  if(show_labels) {
    p <- p + geom_tiplab(offset = 0.1, hjust = 0, size = 3)
  }
  
  # Find min and max values for scaling
  min_chr <- min(node_data$chr_num, na.rm = TRUE)
  max_chr <- max(node_data$chr_num, na.rm = TRUE)
  
  # Create heatmap data
  heat_data <- data.frame(
    node = node_data$node,
    value = node_data$chr_num
  )
  
  return(p)
}

#' Plot fusion-fission ratio across a phylogeny
#' 
#' Creates a visualization of the ratio between fusion and fission events
#' across different parts of the phylogeny
#' 
#' @param tree Phylogenetic tree
#' @param events Results from detect_chromosome_events function
#' @param window_size Number of branches to include in sliding window
#' @param color_ratio Whether to color branches by fusion/fission ratio
#' @param highlight_bias Whether to highlight branches with significant bias
#' @return ggplot object with fusion-fission ratio visualization
#' @export
plot_fusion_fission_ratio <- function(tree,
                                    events,
                                    window_size = 5,
                                    color_ratio = TRUE,
                                    highlight_bias = TRUE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!inherits(events, "chr_events")) {
    stop("events must be a chr_events object from detect_chromosome_events")
  }
  
  # Extract branch data
  branch_data <- events$branch_data
  
  # Calculate fusion-fission ratio for each branch
  branch_data$fusion_count <- ifelse(branch_data$event_type == "fusion", branch_data$event_count, 0)
  branch_data$fission_count <- ifelse(branch_data$event_type == "fission", branch_data$event_count, 0)
  
  # Calculate ratio (handle division by zero)
  branch_data$ff_ratio <- ifelse(branch_data$fission_count > 0,
                              branch_data$fusion_count / branch_data$fission_count,
                              ifelse(branch_data$fusion_count > 0, Inf, NA))
  
  # Convert infinite values to large number for plotting
  branch_data$ff_ratio[is.infinite(branch_data$ff_ratio)] <- 10
  
  # Create log ratio for better visualization
  branch_data$log_ratio <- log(branch_data$ff_ratio + 0.1)
  branch_data$log_ratio[is.na(branch_data$log_ratio)] <- 0
  
  # Set up plot
  p <- ggtree::ggtree(tree)
  
  # Color branches by ratio if requested
  if(color_ratio) {
    p <- p %<+% branch_data + aes(color = log_ratio) +
      scale_color_gradient2(name = "Fusion/Fission\nRatio (log)",
                          low = "blue", mid = "gray", high = "red", midpoint = 0,
                          na.value = "gray")
  }
  
  # Add dots to highlight biased branches
  if(highlight_bias) {
    # Define biased branches (fusion-biased or fission-biased)
    branch_data$bias <- "none"
    
    # At least 3 events and 2:1 ratio for clear bias
    fusion_biased <- branch_data$fusion_count >= 3 & 
                    branch_data$fusion_count >= 2 * branch_data$fission_count
    fission_biased <- branch_data$fission_count >= 3 & 
                     branch_data$fission_count >= 2 * branch_data$fusion_count
    
    branch_data$bias[fusion_biased] <- "fusion"
    branch_data$bias[fission_biased] <- "fission"
    
    p <- p + geom_point2(data = branch_data,
                       aes(subset = (bias != "none"), color = bias),
                       size = 3, alpha = 0.7) +
      scale_color_manual(values = c("fusion" = "blue", 
                                   "fission" = "red", 
                                   "none" = "gray"))
  }
  
  # Add title and theme
  p <- p + labs(title = "Fusion-Fission Ratio Across Phylogeny") +
    theme_tree2() + theme(legend.position = "right")
  
  return(p)
}

#===============================================================================
# Advanced Analysis Functions
#===============================================================================

#' Analyze polyploidy patterns across a phylogeny
#' 
#' Identifies putative whole genome duplication (WGD) events and analyzes
#' their phylogenetic distribution and impact on diversification
#' 
#' @param tree Phylogenetic tree
#' @param events Results from detect_chromosome_events function
#' @param min_fold_increase Minimum fold increase to consider as polyploidy
#' @param min_branches Minimum number of descendant branches to analyze impact
#' @param test_diversification Whether to test impact on diversification
#' @return List with polyploidy analysis results
#' @export
analyze_polyploidy_patterns <- function(tree,
                                      events,
                                      min_fold_increase = 1.8,
                                      min_branches = 3,
                                      test_diversification = TRUE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!inherits(events, "chr_events")) {
    stop("events must be a chr_events object from detect_chromosome_events")
  }
  
  # Extract branch data
  branch_data <- events$branch_data
  
  # Find putative WGD events
  putative_wgd <- which(branch_data$event_type == "polyploidy")
  
  # Initialize results
  results <- list(
    n_polyploidy_events = length(putative_wgd),
    polyploidy_nodes = branch_data$child[putative_wgd],
    summary = list(),
    events = list()
  )
  
  if(length(putative_wgd) == 0) {
    results$summary$text <- "No polyploidy events detected"
    return(results)
  }
  
  # Process each polyploidy event
  event_details <- list()
  
  for(i in 1:length(putative_wgd)) {
    event_idx <- putative_wgd[i]
    child_node <- branch_data$child[event_idx]
    parent_node <- branch_data$parent[event_idx]
    
    # Get chromosome numbers
    parent_chr <- branch_data$parent_chr[event_idx]
    child_chr <- branch_data$child_chr[event_idx]
    
    # Calculate fold increase
    fold_increase <- child_chr / parent_chr
    
    # Determine ploidy level
    ploidy_level <- round(fold_increase)
    ploidy_type <- case_when(
      ploidy_level == 2 ~ "diploidization",
      ploidy_level == 3 ~ "triploidization",
      ploidy_level == 4 ~ "tetraploidization",
      TRUE ~ paste0(ploidy_level, "-ploidization")
    )
    
    # Get all descendants of the WGD node
    descendants <- NULL
    if(requireNamespace("phytools", quietly = TRUE)) {
      descendants <- phytools::getDescendants(tree, child_node)
    } else {
      # Fallback method if phytools not available
      descendants <- get_descendants_from_tree(tree, child_node)
    }
    
    # Get descendant tips
    n_tips <- length(tree$tip.label)
    descendant_tips <- descendants[descendants <= n_tips]
    
    # Create event details
    event_details[[i]] <- list(
      event_id = i,
      parent_node = parent_node,
      node = child_node,
      parent_chr = parent_chr,
      child_chr = child_chr,
      fold_increase = fold_increase,
      ploidy_level = ploidy_level,
      ploidy_type = ploidy_type,
      n_descendants = length(descendants),
      n_descendant_tips = length(descendant_tips),
      descendant_tips = descendant_tips
    )
  }
  
  results$events <- event_details
  
  # Test impact on diversification if requested
  if(test_diversification && length(putative_wgd) > 0 && 
     requireNamespace("phytools", quietly = TRUE)) {
    
    div_results <- list()
    
    for(i in 1:length(event_details)) {
      event <- event_details[[i]]
      
      # Only test if enough descendants
      if(event$n_descendant_tips >= min_branches) {
        # Create a binary trait: descendants of WGD vs non-descendants
        wgd_state <- rep(0, length(tree$tip.label))
        wgd_state[event$descendant_tips] <- 1
        names(wgd_state) <- tree$tip.label
        
        # Test for diversification rate differences
        div_test <- tryCatch({
          # Attempt to test for state-dependent diversification
          if(requireNamespace("diversitree", quietly = TRUE)) {
            # This would require a more complex implementation with diversitree
            # For now, we'll use a simple comparison of species richness
            list(
              n_with_wgd = sum(wgd_state),
              n_without_wgd = length(wgd_state) - sum(wgd_state),
              wgd_age = get_node_age(tree, event$node),
              wgd_rate = sum(wgd_state) / get_node_age(tree, event$node),
              test_result = "Species richness comparison only"
            )
          } else {
            list(
              n_with_wgd = sum(wgd_state),
              n_without_wgd = length(wgd_state) - sum(wgd_state),
              test_result = "diversitree package required for full testing"
            )
          }
        }, error = function(e) {
          list(
            error = e$message,
            test_result = "Error in diversification test"
          )
        })
        
        div_results[[i]] <- div_test
      } else {
        div_results[[i]] <- list(
          n_descendant_tips = event$n_descendant_tips,
          test_result = "Too few descendants for meaningful testing"
        )
      }
    }
    
    results$diversification <- div_results
  }
  
  # Generate summary
  results$summary$text <- paste0("Detected ", length(putative_wgd), " putative polyploidy events")
  
  # Calculate age distribution of WGD events if tree is dated
  if(!is.null(tree$edge.length)) {
    event_ages <- sapply(results$events, function(e) {
      get_node_age(tree, e$node)
    })
    
    results$summary$age_range <- range(event_ages)
    results$summary$mean_age <- mean(event_ages)
    results$summary$median_age <- median(event_ages)
  }
  
  return(results)
}

#' Helper function to get all descendants of a node
#' 
#' @param tree Phylogenetic tree
#' @param node Node ID
#' @return Vector of descendant node IDs
#' @keywords internal
get_descendants_from_tree <- function(tree, node) {
  # Initialize with direct descendants
  edge_idx <- which(tree$edge[, 1] == node)
  desc <- tree$edge[edge_idx, 2]
  
  # Recursively get descendants of each child
  for(child in desc) {
    child_desc <- get_descendants_from_tree(tree, child)
    desc <- c(desc, child_desc)
  }
  
  return(unique(desc))
}

#' Get age of a node in a phylogenetic tree
#' 
#' @param tree Phylogenetic tree
#' @param node Node ID
#' @return Age of the node
#' @keywords internal
get_node_age <- function(tree, node) {
  if(is.null(tree$edge.length)) {
    return(NA)
  }
  
  # Get all tip descendants
  n_tips <- length(tree$tip.label)
  if(node <= n_tips) {
    # It's a tip, so age is 0
    return(0)
  }
  
  # For internal nodes, calculate mean path length to tips
  descendants <- get_descendants_from_tree(tree, node)
  tip_descendants <- descendants[descendants <= n_tips]
  
  # If no tip descendants found, return NA
  if(length(tip_descendants) == 0) {
    return(NA)
  }
  
  # Calculate paths from node to each tip
  paths <- lapply(tip_descendants, function(tip) {
    get_path_lengths(tree, node, tip)
  })
  
  # Return mean path length
  return(mean(unlist(paths)))
}

#' Get path length between two nodes
#' 
#' @param tree Phylogenetic tree
#' @param start_node Starting node ID
#' @param end_node Ending node ID
#' @return Path length between nodes
#' @keywords internal
get_path_lengths <- function(tree, start_node, end_node) {
  # Find path between nodes
  current <- end_node
  path_length <- 0
  
  while(current != start_node) {
    # Find parent of current node
    edge_idx <- which(tree$edge[, 2] == current)
    
    # If no parent found, path doesn't exist
    if(length(edge_idx) == 0) {
      return(NA)
    }
    
    # Add branch length to path
    path_length <- path_length + tree$edge.length[edge_idx]
    
    # Move to parent
    current <- tree$edge[edge_idx, 1]
  }
  
  return(path_length)
}

#' Test for chromosome number evolutionary models
#' 
#' Tests different models of chromosome number evolution
#' (constant rates, varying rates, directional trends)
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param models Character vector of models to test
#' @param method Model fitting method: "ML", "REML", "Bayesian"
#' @param criterion Model selection criterion: "AIC", "BIC", "likelihood"
#' @return List with model fitting and comparison results
#' @export
test_chromosome_models <- function(tree,
                                 chr_counts,
                                 models = c("BM", "OU", "EB", "trend"),
                                 method = "ML",
                                 criterion = "AIC") {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Get matching taxa
  shared_taxa <- intersect(tree$tip.label, names(chr_counts))
  
  if(length(shared_taxa) < 4) {
    stop("At least 4 shared taxa required for model testing")
  }
  
  # Prune tree and data to shared taxa
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, shared_taxa))
  pruned_counts <- chr_counts[shared_taxa]
  
  # Initialize results
  results <- list(
    model_fits = list(),
    model_selection = NULL,
    best_model = NULL,
    summary = list(
      n_taxa = length(shared_taxa),
      models_tested = models,
      criterion = criterion
    )
  )
  
  # Test each model
  if(requireNamespace("geiger", quietly = TRUE)) {
    # Fit models with fitContinuous
    model_fits <- list()
    
    for(model in models) {
      if(model %in% c("BM", "OU", "EB")) {
        # These models are directly supported by geiger
        fit <- tryCatch({
          geiger::fitContinuous(pruned_tree, pruned_counts, model = model, 
                              method = method, control = list(maxit = 1000))
        }, error = function(e) {
          warning(paste("Error fitting", model, "model:", e$message))
          return(NULL)
        })
        
        model_fits[[model]] <- fit
      } else if(model == "trend") {
        # Trend model (BM with trend parameter)
        fit <- tryCatch({
          # This model may not be directly available in geiger
          # Could implement custom model or use other packages
          warning("Trend model not currently implemented")
          return(NULL)
        }, error = function(e) {
          warning(paste("Error fitting trend model:", e$message))
          return(NULL)
        })
        
        model_fits[[model]] <- fit
      } else {
        warning(paste("Model", model, "not supported"))
      }
    }
    
    # Filter out failed fits
    model_fits <- model_fits[!sapply(model_fits, is.null)]
    
    # Extract model selection criteria
    if(length(model_fits) > 0) {
      model_selection <- data.frame(
        Model = names(model_fits),
        logLik = sapply(model_fits, function(x) x$opt$lnL),
        AIC = sapply(model_fits, function(x) x$opt$aic),
        Parameters = sapply(model_fits, function(x) length(x$opt$npars))
      )
      
      # Calculate BIC
      model_selection$BIC <- model_selection$AIC + 
                           model_selection$Parameters * log(length(shared_taxa)) - 
                           model_selection$Parameters * 2
      
      # Sort by selected criterion
      if(criterion == "AIC") {
        model_selection <- model_selection[order(model_selection$AIC), ]
        best_model <- model_selection$Model[1]
      } else if(criterion == "BIC") {
        model_selection <- model_selection[order(model_selection$BIC), ]
        best_model <- model_selection$Model[1]
      } else {
        model_selection <- model_selection[order(-model_selection$logLik), ]
        best_model <- model_selection$Model[1]
      }
      
      # Calculate model weights
      if(criterion == "AIC") {
        delta_aic <- model_selection$AIC - min(model_selection$AIC)
        weight <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
        model_selection$AIC_weight <- weight
      } else if(criterion == "BIC") {
        delta_bic <- model_selection$BIC - min(model_selection$BIC)
        weight <- exp(-0.5 * delta_bic) / sum(exp(-0.5 * delta_bic))
        model_selection$BIC_weight <- weight
      }
      
      results$model_fits <- model_fits
      results$model_selection <- model_selection
      results$best_model <- best_model
    }
  } else {
    warning("geiger package is required for model testing")
  }
  
  # Create visualization if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE) && 
     !is.null(results$model_selection)) {
    
    # Create AIC/BIC plot
    criterion_col <- paste0(criterion, "_weight")
    if(criterion_col %in% colnames(results$model_selection)) {
      weight_plot <- ggplot2::ggplot(results$model_selection, 
                                  ggplot2::aes(x = reorder(Model, .data[[criterion_col]]), 
                                            y = .data[[criterion_col]])) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
        ggplot2::labs(title = paste(criterion, "Weights for Chromosome Evolution Models"),
                    x = "Model",
                    y = paste(criterion, "Weight")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      results$plots$weight_plot <- weight_plot
    }
  }
  
  # Interpret results
  if(!is.null(results$best_model)) {
    best_fit <- model_fits[[results$best_model]]
    
    interpretation <- switch(results$best_model,
      "BM" = "Chromosome evolution follows a Brownian Motion model, suggesting random drift without directional trends.",
      "OU" = paste0("Chromosome evolution follows an Ornstein-Uhlenbeck model with alpha = ", 
                  round(best_fit$opt$alpha, 4), 
                  ", suggesting selection toward an optimal chromosome number."),
      "EB" = paste0("Chromosome evolution follows an Early Burst model with parameter = ", 
                  round(best_fit$opt$a, 4), 
                  ", suggesting rapid early diversification followed by slowdown."),
      "trend" = "Chromosome evolution shows a directional trend.",
      "No clear interpretation available"
    )
    
    results$summary$interpretation <- interpretation
  }
  
  return(results)
}

#' Analyze chromosome evolution by taxonomic group
#' 
#' Compares chromosome evolution patterns across different taxonomic groups
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param group_data Named vector or data frame with group assignments
#' @param test_method Method for testing group differences
#' @param rates_method Method for calculating rates
#' @return List with group comparison results
#' @export
compare_taxonomic_groups <- function(tree,
                                   chr_counts,
                                   group_data,
                                   test_method = "anova",
                                   rates_method = "contrasts") {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Process group data
  if(is.vector(group_data)) {
    # It's a vector, check names
    if(is.null(names(group_data))) {
      stop("group_data must have names corresponding to species")
    }
    
    groups <- group_data
  } else if(is.data.frame(group_data)) {
    # It's a data frame, extract groups
    if(is.null(rownames(group_data))) {
      stop("group_data must have row names corresponding to species")
    }
    
    if("group" %in% colnames(group_data)) {
      groups <- group_data$group
      names(groups) <- rownames(group_data)
    } else {
      stop("group_data must have a 'group' column")
    }
  } else {
    stop("group_data must be a named vector or data frame")
  }
  
  # Get common taxa
  shared_taxa <- intersect(tree$tip.label, names(chr_counts))
  shared_taxa <- intersect(shared_taxa, names(groups))
  
  if(length(shared_taxa) < 4) {
    stop("At least 4 shared taxa required for group comparison")
  }
  
  # Prune tree and data to shared taxa
  pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, shared_taxa))
  pruned_counts <- chr_counts[shared_taxa]
  pruned_groups <- groups[shared_taxa]
  
  # Ensure groups is a factor
  pruned_groups <- as.factor(pruned_groups)
  
  # Count species per group
  group_counts <- table(pruned_groups)
  
  # Filter out groups with too few species
  valid_groups <- names(group_counts)[group_counts >= 3]
  
  if(length(valid_groups) < 2) {
    stop("At least 2 groups with 3+ species required for comparison")
  }
  
  # Filter data to valid groups
  in_valid_group <- pruned_groups %in% valid_groups
  filtered_taxa <- shared_taxa[in_valid_group]
  filtered_tree <- ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, filtered_taxa))
  filtered_counts <- pruned_counts[filtered_taxa]
  filtered_groups <- pruned_groups[filtered_taxa]
  
  # Calculate basic statistics by group
  group_stats <- data.frame(
    Group = valid_groups,
    N = as.numeric(group_counts[valid_groups]),
    Mean = tapply(filtered_counts, filtered_groups, mean),
    Median = tapply(filtered_counts, filtered_groups, median),
    Min = tapply(filtered_counts, filtered_groups, min),
    Max = tapply(filtered_counts, filtered_groups, max),
    Range = tapply(filtered_counts, filtered_groups, function(x) max(x) - min(x)),
    SD = tapply(filtered_counts, filtered_groups, sd)
  )
  
  # Initialize results
  results <- list(
    group_stats = group_stats,
    group_test = NULL,
    rate_test = NULL,
    summary = list(
      n_taxa = length(filtered_taxa),
      n_groups = length(valid_groups),
      groups = valid_groups
    )
  )
  
  # Test for differences in chromosome numbers by group
  if(test_method == "anova") {
    # Standard ANOVA
    anova_result <- aov(filtered_counts ~ filtered_groups)
    anova_summary <- summary(anova_result)
    
    # Post-hoc test
    if(requireNamespace("multcomp", quietly = TRUE) && length(valid_groups) > 2) {
      posthoc <- multcomp::glht(anova_result, linfct = multcomp::mcp(filtered_groups = "Tukey"))
      posthoc_summary <- summary(posthoc)
      
      results$group_test <- list(
        method = "ANOVA with Tukey post-hoc",
        anova = anova_summary,
        p_value = anova_summary[[1]]["filtered_groups", "Pr(>F)"],
        significant = anova_summary[[1]]["filtered_groups", "Pr(>F)"] < 0.05,
        posthoc = posthoc_summary
      )
    } else {
      results$group_test <- list(
        method = "ANOVA",
        anova = anova_summary,
        p_value = anova_summary[[1]]["filtered_groups", "Pr(>F)"],
        significant = anova_summary[[1]]["filtered_groups", "Pr(>F)"] < 0.05
      )
    }
  } else if(test_method == "phylo_anova") {
    # Phylogenetic ANOVA
    if(requireNamespace("phytools", quietly = TRUE)) {
      # Convert to named vector for phytools
      group_vec <- filtered_groups
      names(group_vec) <- filtered_taxa
      
      phyanova_result <- tryCatch({
        phytools::phylANOVA(filtered_tree, group_vec, filtered_counts, posthoc = TRUE)
      }, error = function(e) {
        warning(paste("phylANOVA error:", e$message))
        return(NULL)
      })
      
      if(!is.null(phyanova_result)) {
        results$group_test <- list(
          method = "Phylogenetic ANOVA",
          phyanova = phyanova_result,
          p_value = phyanova_result$Pf,
          significant = phyanova_result$Pf < 0.05
        )
      }
    } else {
      warning("phytools package required for phylogenetic ANOVA")
    }
  }
  
  # Test for differences in rates by group
  if(rates_method == "contrasts") {
    # Calculate PICs for each group
    if(requireNamespace("ape", quietly = TRUE)) {
      group_contrasts <- list()
      
      for(group in valid_groups) {
        # Get taxa in this group
        group_taxa <- filtered_taxa[filtered_groups == group]
        
        # Skip if too few species
        if(length(group_taxa) < 4) next
        
        # Prune tree and get counts
        group_tree <- ape::drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, group_taxa))
        group_counts <- filtered_counts[group_taxa]
        
        # Calculate contrasts
        group_pics <- ape::pic(group_counts, group_tree)
        
        # Store absolute contrasts (measure of rate)
        group_contrasts[[group]] <- abs(group_pics)
      }
      
      # Test for differences in contrast magnitudes
      if(length(group_contrasts) >= 2) {
        # Convert to data frame for testing
        contrast_data <- stack(lapply(group_contrasts, function(x) {
          as.numeric(x)
        }))
        colnames(contrast_data) <- c("contrast", "group")
        
        # Test with Kruskal-Wallis (non-parametric) due to potential non-normality
        kw_test <- kruskal.test(contrast ~ group, data = contrast_data)
        
        # Pairwise comparisons if >2 groups
        pairwise_tests <- NULL
        if(length(group_contrasts) > 2) {
          pairwise_tests <- pairwise.wilcox.test(contrast_data$contrast, 
                                              contrast_data$group, 
                                              p.adjust.method = "fdr")
        }
        
        # Calculate mean absolute contrasts per group
        mean_contrasts <- sapply(group_contrasts, mean)
        
        results$rate_test <- list(
          method = "PIC magnitude comparison",
          test = kw_test,
          p_value = kw_test$p.value,
          significant = kw_test$p.value < 0.05,
          pairwise = pairwise_tests,
          mean_contrasts = mean_contrasts,
          group_with_highest_rate = names(which.max(mean_contrasts))
        )
      }
    } else {
      warning("ape package required for PIC calculation")
    }
  }
  
  # Create visualization if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Create boxplot of chromosome numbers by group
    p1 <- ggplot2::ggplot(data.frame(
      Chromosome_Number = filtered_counts,
      Group = filtered_groups
    ), ggplot2::aes(x = Group, y = Chromosome_Number, fill = Group)) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(title = "Chromosome Numbers by Group",
                  x = "Taxonomic Group",
                  y = "Chromosome Number") +
      ggplot2::theme_minimal()
    
    results$plots$boxplot <- p1
    
    # Create rates comparison plot if available
    if(!is.null(results$rate_test) && !is.null(results$rate_test$mean_contrasts)) {
      p2 <- ggplot2::ggplot(data.frame(
        Group = names(results$rate_test$mean_contrasts),
        Rate = results$rate_test$mean_contrasts
      ), ggplot2::aes(x = reorder(Group, Rate), y = Rate, fill = Group)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = "Chromosome Evolution Rate by Group",
                    subtitle = "Based on mean absolute phylogenetic independent contrasts",
                    x = "Taxonomic Group",
                    y = "Evolution Rate") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      results$plots$rate_plot <- p2
    }
  }
  
  # Create summary interpretation
  if(!is.null(results$group_test) && results$group_test$significant) {
    results$summary$number_interpretation <- paste(
      "Significant differences in chromosome numbers were detected among taxonomic groups",
      "(p =", format(results$group_test$p_value, digits = 3), ")"
    )
  } else if(!is.null(results$group_test)) {
    results$summary$number_interpretation <- paste(
      "No significant differences in chromosome numbers were detected among taxonomic groups",
      "(p =", format(results$group_test$p_value, digits = 3), ")"
    )
  }
  
  if(!is.null(results$rate_test) && results$rate_test$significant) {
    results$summary$rate_interpretation <- paste0(
      "Significant differences in chromosome evolution rates were detected among taxonomic groups ",
      "(p = ", format(results$rate_test$p_value, digits = 3), "). ",
      "The highest rate was observed in the ", results$rate_test$group_with_highest_rate, " group."
    )
  } else if(!is.null(results$rate_test)) {
    results$summary$rate_interpretation <- paste(
      "No significant differences in chromosome evolution rates were detected among taxonomic groups",
      "(p =", format(results$rate_test$p_value, digits = 3), ")"
    )
  }
  
  return(results)
}

#===============================================================================
# Reporting Functions
#===============================================================================

#' Generate a report summarizing karyotype evolution analysis
#' 
#' Creates a comprehensive report of chromosome evolution patterns and events
#' 
#' @param tree Phylogenetic tree
#' @param events Chromosome events detected using detect_chromosome_events
#' @param chr_counts Named vector of chromosome counts
#' @param ancestral_states Ancestral state reconstructions
#' @param rates Optional results from calculate_chromosome_rates
#' @param polyploidy Optional results from analyze_polyploidy_patterns
#' @param model_fits Optional results from test_chromosome_models
#' @param output_file Path to save the report (NULL for no saving)
#' @param format Output format: "html", "pdf", "docx"
#' @param include_plots Whether to include visualization plots in report
#' @return HTML report content as a character vector
#' @export
generate_karyotype_report <- function(tree,
                                    events,
                                    chr_counts,
                                    ancestral_states,
                                    rates = NULL,
                                    polyploidy = NULL,
                                    model_fits = NULL,
                                    output_file = NULL,
                                    format = "html",
                                    include_plots = TRUE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!inherits(events, "chr_events")) {
    stop("events must be a chr_events object from detect_chromosome_events")
  }
  
  # Check if rmarkdown is available for report generation
  if(!is.null(output_file) && !requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("rmarkdown package is required to generate reports")
  }
  
  # Extract basic statistics
  n_tips <- length(tree$tip.label)
  chr_range <- range(chr_counts, na.rm = TRUE)
  chr_mean <- mean(chr_counts, na.rm = TRUE)
  chr_median <- median(chr_counts, na.rm = TRUE)
  
  # Event statistics
  n_fusion <- events$summary$fusion_count
  n_fission <- events$summary$fission_count
  n_polyploidy <- events$summary$polyploidy_count
  total_events <- events$summary$total_events
  
  # Create report content
  report <- c(
    "---",
    "title: \"Karyotype Evolution Analysis Report\"",
    paste0("date: \"", format(Sys.Date(), "%B %d, %Y"), "\""),
    "output: html_document",
    "---",
    "",
    "## Overview",
    "",
    paste("This report summarizes the analysis of chromosome evolution across",
          n_tips, "species."),
    "",
    "### Summary Statistics",
    "",
    "| Metric | Value |",
    "| --- | --- |",
    paste("| Number of species | ", n_tips, " |"),
    paste("| Chromosome number range | ", chr_range[1], " - ", chr_range[2], " |"),
    paste("| Mean chromosome number | ", round(chr_mean, 2), " |"),
    paste("| Median chromosome number | ", chr_median, " |"),
    paste("| Total evolutionary events | ", total_events, " |"),
    paste("| Fusion events | ", n_fusion, " |"),
    paste("| Fission events | ", n_fission, " |"),
    paste("| Polyploidy events | ", n_polyploidy, " |"),
    "",
    "## Chromosome Evolution Events",
    "",
    paste("A total of", total_events, "chromosome rearrangement events were",
          "detected across the phylogeny. This includes", n_fusion, "fusion events,",
          n_fission, "fission events, and", n_polyploidy, "polyploidy events."),
    ""
  )
  
  # Add fusion-fission ratio information
  fusion_fission_ratio <- ifelse(n_fission > 0, n_fusion / n_fission, Inf)
  if(is.finite(fusion_fission_ratio)) {
    report <- c(report,
              paste("The fusion to fission ratio is", round(fusion_fission_ratio, 2),
                    "indicating a", ifelse(fusion_fission_ratio > 1, 
                                        "bias toward chromosome number reduction.", 
                                        "bias toward chromosome number increase.")),
              "")
  }
  
  # Add model testing results if available
  if(!is.null(model_fits) && !is.null(model_fits$best_model)) {
    report <- c(report,
              "## Evolutionary Model Testing",
              "",
              paste("The best-fitting model for chromosome number evolution was",
                    model_fits$best_model, "based on", model_fits$summary$criterion, "."),
              "",
              "### Model Selection Table",
              "",
              "```{r echo=FALSE}",
              "knitr::kable(model_fits$model_selection)",
              "```",
              "",
              "### Interpretation",
              "",
              model_fits$summary$interpretation,
              "")
  }
  
  # Add polyploidy analysis if available
  if(!is.null(polyploidy) && polyploidy$n_polyploidy_events > 0) {
    report <- c(report,
              "## Polyploidy Analysis",
              "",
              paste("Detected", polyploidy$n_polyploidy_events, 
                    "putative polyploidy events across the phylogeny."),
              "",
              "### Polyploidy Events",
              "",
              "```{r echo=FALSE}",
              "event_df <- data.frame(",
              "  Node = sapply(polyploidy$events, function(e) e$node),",
              "  Parent_Chr = sapply(polyploidy$events, function(e) e$parent_chr),",
              "  Child_Chr = sapply(polyploidy$events, function(e) e$child_chr),",
              "  Fold_Increase = sapply(polyploidy$events, function(e) round(e$fold_increase, 2)),",
              "  Ploidy_Type = sapply(polyploidy$events, function(e) e$ploidy_type),",
              "  Descendant_Tips = sapply(polyploidy$events, function(e) e$n_descendant_tips)",
              ")",
              "knitr::kable(event_df)",
              "```",
              "")
    
    # Add diversification impact if available
    if(!is.null(polyploidy$diversification)) {
      report <- c(report,
                "### Impact on Diversification",
                "",
                "The following table shows the potential impact of polyploidy events on species diversification:",
                "",
                "```{r echo=FALSE}",
                "div_data <- data.frame(",
                "  Event = 1:length(polyploidy$diversification),",
                "  N_With_WGD = sapply(polyploidy$diversification, function(d) ifelse(is.null(d$n_with_wgd), NA, d$n_with_wgd)),",
                "  N_Without_WGD = sapply(polyploidy$diversification, function(d) ifelse(is.null(d$n_without_wgd), NA, d$n_without_wgd)),",
                "  Result = sapply(polyploidy$diversification, function(d) ifelse(is.null(d$test_result), NA, d$test_result))",
                ")",
                "knitr::kable(div_data)",
                "```",
                "")
    }
  }
  
  # Add rate analysis if available
  if(!is.null(rates)) {
    report <- c(report,
              "## Chromosome Evolution Rates",
              "",
              "The analysis of rates revealed the following patterns in chromosome number evolution:",
              "",
              "### Rate Summary",
              "",
              "| Event Type | Rate |",
              "| --- | --- |",
              paste("| Overall | ", round(rates$rate_summary$overall, 4), " |"),
              paste("| Fusion | ", round(rates$rate_summary$fusion, 4), " |"),
              paste("| Fission | ", round(rates$rate_summary$fission, 4), " |"),
              paste("| Polyploidy | ", round(rates$rate_summary$polyploidy, 4), " |"),
              "",
              paste("Note: Rates are measured as", 
                    ifelse(rates$rate_metric == "events_per_time", 
                          "events per unit branch length.", 
                          "events per branch.")),
              "")
  }
  
  # Add visualization section if plots included
  if(include_plots) {
    report <- c(report,
              "## Visualizations",
              "",
              "### Chromosome Events Mapped on Phylogeny",
              "",
              "```{r fig.width=10, fig.height=8, echo=FALSE}",
              "plot_chromosome_events(tree, events, show_counts = TRUE)",
              "```",
              "",
              "### Fusion-Fission Ratio",
              "",
              "```{r fig.width=10, fig.height=8, echo=FALSE}",
              "plot_fusion_fission_ratio(tree, events)",
              "```",
              "")
    
    # Add chromosome number heatmap
    if(!is.null(ancestral_states)) {
      report <- c(report,
                "### Chromosome Number Evolution",
                "",
                "```{r fig.width=10, fig.height=8, echo=FALSE}",
                "plot_chromosome_heatmap(tree, events, ancestral_states)",
                "```",
                "")
    }
    
    # Add model plots if available
    if(!is.null(model_fits) && !is.null(model_fits$plots$weight_plot)) {
      report <- c(report,
                "### Model Selection",
                "",
                "```{r fig.width=8, fig.height=6, echo=FALSE}",
                "print(model_fits$plots$weight_plot)",
                "```",
                "")
    }
  }
  
  # Add conclusions
  report <- c(report,
            "## Conclusions",
            "",
            "Based on the analysis of chromosome evolution patterns:")
  
  # Add fusion-fission bias conclusion
  if(n_fusion > 0 || n_fission > 0) {
    if(fusion_fission_ratio > 1.5) {
      report <- c(report,
                "- There is a strong bias toward chromosome fusions, suggesting evolutionary pressure for chromosome number reduction.")
    } else if(fusion_fission_ratio < 0.67) {
      report <- c(report,
                "- There is a bias toward chromosome fissions, suggesting evolutionary pressure for chromosome number increase.")
    } else {
      report <- c(report,
                "- There is a roughly equal balance between fusion and fission events, suggesting no strong directional bias in chromosome number evolution.")
    }
  }
  
  # Add model conclusion
  if(!is.null(model_fits) && !is.null(model_fits$best_model)) {
    if(model_fits$best_model == "BM") {
      report <- c(report,
                "- Chromosome numbers evolve in a manner consistent with random drift (Brownian motion model).")
    } else if(model_fits$best_model == "OU") {
      report <- c(report,
                "- Chromosome numbers evolve under stabilizing selection toward an optimal value (Ornstein-Uhlenbeck model).")
    } else if(model_fits$best_model == "EB") {
      report <- c(report,
                "- Chromosome evolution shows a pattern of early rapid change followed by slowdown (Early Burst model).")
    }
  }
  
  # Add polyploidy conclusion
  if(!is.null(polyploidy) && polyploidy$n_polyploidy_events > 0) {
    report <- c(report,
              paste0("- ", polyploidy$n_polyploidy_events, " polyploidy events were detected, ",
                   "suggesting whole genome duplication has played a role in karyotype evolution of this group."))
  } else {
    report <- c(report,
              "- No clear evidence of polyploidy was detected in the evolutionary history of this group.")
  }
  
  # Generate the report
  if(!is.null(output_file)) {
    # Create a temporary Rmd file
    temp_rmd <- tempfile(fileext = ".Rmd")
    writeLines(report, temp_rmd)
    
    # Render the report
    rmarkdown::render(temp_rmd, 
                    output_format = paste0(format, "_document"), 
                    output_file = basename(output_file), 
                    output_dir = dirname(output_file))
    
    # Clean up temp file
    unlink(temp_rmd)
  }
  
  # Return report text
  return(report)
}
