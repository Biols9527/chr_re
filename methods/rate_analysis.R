#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Rate Analysis Module
# Author: Bioinformatics Team
# Date: 2025-05-20
# Description: Specialized methods for analyzing rates of chromosome number
#              evolution across phylogenies, including rate variation detection,
#              shift point analysis, and heterogeneity testing
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(geiger)
  library(l1ou)  # For OU process shift detection
  library(RRphylo) # For relative rates
  library(bayou)  # For Bayesian shift detection
  library(parallel)
  library(ggplot2)
  library(viridis)
})

#===============================================================================
# Rate Variation Detection Functions
#===============================================================================

#' Analyze rate heterogeneity in chromosome number evolution
#' 
#' Tests for variation in rates of chromosome number evolution across a phylogeny
#' using multiple methods and identifies branches with exceptional rates
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param methods Vector of methods to use: "contrasts", "rphylo", "shifts", "bayou", "all"
#' @param model Base evolutionary model: "BM", "OU", "EB" (only for some methods)
#' @param n_simulations Number of simulations for null distribution (for permutation tests)
#' @param alpha Significance level for hypothesis tests
#' @param min_clade_size Minimum size for a clade to be considered for shift analysis
#' @param n_cores Number of CPU cores for parallel processing (NULL = auto-detect)
#' @return List with rate analysis results
#' @export
analyze_rate_heterogeneity <- function(tree, 
                                     chr_counts,
                                     methods = "all",
                                     model = "BM",
                                     n_simulations = 1000,
                                     alpha = 0.05,
                                     min_clade_size = 4,
                                     n_cores = NULL) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Match data to tree
  chr_data <- match_data_to_tree(tree, chr_counts)
  
  # Drop tips with missing data
  has_data <- !is.na(chr_data)
  if(sum(has_data) < 4) {
    stop("At least 4 tips must have chromosome count data")
  }
  
  pruned_tree <- ape::drop.tip(tree, tree$tip.label[!has_data])
  pruned_data <- chr_data[has_data]
  
  # Determine which methods to use
  all_methods <- c("contrasts", "rphylo", "shifts", "bayou")
  if(length(methods) == 1 && methods == "all") {
    methods <- all_methods
  } else {
    methods <- match.arg(methods, all_methods, several.ok = TRUE)
  }
  
  # Set up parallel processing
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Initialize results
  results <- list(
    tree = pruned_tree,
    chr_counts = pruned_data,
    methods = methods,
    model = model,
    summary = list(),
    rate_estimates = list(),
    shift_points = list(),
    significance_tests = list(),
    plots = list()
  )
  
  # Apply each method
  if("contrasts" %in% methods) {
    contrast_results <- analyze_rate_contrasts(pruned_tree, pruned_data, 
                                            n_simulations, alpha)
    results$rate_estimates$contrasts <- contrast_results$rates
    results$significance_tests$contrasts <- contrast_results$tests
    results$summary$contrasts <- contrast_results$summary
    results$plots$contrasts <- contrast_results$plots
  }
  
  if("rphylo" %in% methods) {
    # Only run if RRphylo is available
    if(requireNamespace("RRphylo", quietly = TRUE)) {
      rphylo_results <- analyze_rate_rphylo(pruned_tree, pruned_data, 
                                         n_cores, alpha)
      results$rate_estimates$rphylo <- rphylo_results$rates
      results$shift_points$rphylo <- rphylo_results$shifts
      results$summary$rphylo <- rphylo_results$summary
      results$plots$rphylo <- rphylo_results$plots
    } else {
      warning("RRphylo package not available. Skipping RRphylo analysis.")
    }
  }
  
  if("shifts" %in% methods) {
    # Only run if l1ou is available
    if(requireNamespace("l1ou", quietly = TRUE)) {
      l1ou_results <- analyze_rate_shifts(pruned_tree, pruned_data, 
                                       min_clade_size, alpha)
      results$shift_points$l1ou <- l1ou_results$shifts
      results$summary$l1ou <- l1ou_results$summary
      results$plots$l1ou <- l1ou_results$plots
    } else {
      warning("l1ou package not available. Skipping L1OU shift analysis.")
    }
  }
  
  if("bayou" %in% methods) {
    # Only run if bayou is available
    if(requireNamespace("bayou", quietly = TRUE)) {
      bayou_results <- analyze_rate_bayou(pruned_tree, pruned_data, 
                                       model, n_simulations)
      results$shift_points$bayou <- bayou_results$shifts
      results$summary$bayou <- bayou_results$summary
      results$plots$bayou <- bayou_results$plots
    } else {
      warning("bayou package not available. Skipping Bayesian shift analysis.")
    }
  }
  
  # Create overall summary
  results$summary$overall <- generate_rate_summary(results)
  
  # Generate consensus results by combining all methods
  if(length(methods) > 1) {
    consensus_results <- generate_rate_consensus(results)
    results$consensus <- consensus_results
  }
  
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

#===============================================================================
# Contrast-Based Rate Analysis
#===============================================================================

#' Analyze rate variation using phylogenetic independent contrasts
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param n_simulations Number of simulations for null distribution
#' @param alpha Significance level
#' @return List with contrast-based rate analysis results
#' @keywords internal
analyze_rate_contrasts <- function(tree, 
                                 chr_counts,
                                 n_simulations = 1000,
                                 alpha = 0.05) {
  # Calculate contrasts
  pic_result <- ape::pic(chr_counts, tree)
  
  # Square contrasts (proportional to evolutionary rate)
  squared_contrasts <- pic_result^2
  
  # Get nodal values for plotting
  node_indices <- attr(pic_result, "indices")
  n_tips <- length(tree$tip.label)
  
  # Create data frame of contrast rates
  contrast_data <- data.frame(
    node = node_indices,
    contrast = pic_result,
    squared_contrast = squared_contrasts
  )
  
  # Calculate empirical p-values through simulation
  p_values <- calculate_contrast_pvalues(tree, chr_counts, squared_contrasts, 
                                       n_simulations)
  
  # Identify significant rates
  significant <- p_values < alpha
  
  # Add p-values to results
  contrast_data$p_value <- p_values
  contrast_data$significant <- significant
  
  # Create summary of results
  n_significant <- sum(significant)
  summary_text <- paste0(n_significant, " out of ", length(p_values), 
                      " nodes (", round(100 * n_significant / length(p_values), 1), 
                      "%) show significant rate heterogeneity")
  
  # Create visualization
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Edge-based rate plot
    edge_data <- data.frame(
      edge = 1:nrow(tree$edge),
      parent = tree$edge[, 1],
      child = tree$edge[, 2],
      stringsAsFactors = FALSE
    )
    
    # Match contrast rates to edges
    edge_data$rate <- NA
    edge_data$significant <- FALSE
    
    for(i in 1:nrow(contrast_data)) {
      node_idx <- contrast_data$node[i]
      edge_idx <- which(edge_data$child == node_idx)
      
      if(length(edge_idx) == 1) {
        edge_data$rate[edge_idx] <- contrast_data$squared_contrast[i]
        edge_data$significant[edge_idx] <- contrast_data$significant[i]
      }
    }
    
    # Normalize rates for better visualization
    if(!all(is.na(edge_data$rate))) {
      max_rate <- max(edge_data$rate, na.rm = TRUE)
      if(max_rate > 0) {
        edge_data$scaled_rate <- edge_data$rate / max_rate
      } else {
        edge_data$scaled_rate <- edge_data$rate
      }
    } else {
      edge_data$scaled_rate <- edge_data$rate
    }
    
    # Create plot using ggtree if available, otherwise use base plotting
    if(requireNamespace("ggtree", quietly = TRUE)) {
      tree_plot <- ggtree::ggtree(tree) %<+% edge_data + 
        ggtree::aes(color = scaled_rate) +
        viridis::scale_color_viridis(name = "Relative Rate", na.value = "gray80") +
        ggtree::geom_tiplab() +
        ggplot2::labs(title = "Chromosome Evolution Rate Variation", 
                   subtitle = summary_text) +
        ggplot2::theme_minimal()
      
      # Highlight significant edges
      if(sum(edge_data$significant, na.rm = TRUE) > 0) {
        tree_plot <- tree_plot + 
          ggtree::geom_hilight(data = subset(edge_data, significant), 
                            aes(node = child), alpha = 0.3, extend = 0.5)
      }
    } else {
      # Base plot placeholder
      tree_plot <- NULL
    }
    
    plots <- list(tree_plot = tree_plot)
  } else {
    plots <- list()
  }
  
  # Return results
  return(list(
    rates = contrast_data,
    tests = list(p_values = p_values, alpha = alpha),
    summary = list(text = summary_text, n_significant = n_significant),
    plots = plots
  ))
}

#' Calculate p-values for contrasts through simulation
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param observed_contrasts Observed squared contrasts
#' @param n_simulations Number of simulations
#' @return Vector of p-values
#' @keywords internal
calculate_contrast_pvalues <- function(tree, 
                                     chr_counts,
                                     observed_contrasts,
                                     n_simulations) {
  # Fit BM model to data
  if(requireNamespace("geiger", quietly = TRUE)) {
    fit <- geiger::fitContinuous(tree, chr_counts, model = "BM")
    sigma2 <- fit$opt$sigsq
  } else {
    # Estimate sigma2 from contrasts if geiger not available
    sigma2 <- mean(observed_contrasts)
  }
  
  # Simulate null distribution
  null_contrasts <- matrix(0, nrow = n_simulations, ncol = length(observed_contrasts))
  
  for(i in 1:n_simulations) {
    # Simulate trait under BM
    sim_trait <- phytools::fastBM(tree, sig2 = sigma2)
    
    # Calculate contrasts
    sim_pic <- ape::pic(sim_trait, tree)
    
    # Store squared contrasts
    null_contrasts[i, ] <- sim_pic^2
  }
  
  # Calculate p-values
  p_values <- numeric(length(observed_contrasts))
  
  for(j in 1:length(observed_contrasts)) {
    # Two-tailed test
    p_values[j] <- sum(null_contrasts[, j] >= observed_contrasts[j]) / n_simulations
  }
  
  return(p_values)
}

#===============================================================================
# RRphylo-Based Rate Analysis
#===============================================================================

#' Analyze rate variation using RRphylo method
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param n_cores Number of CPU cores for parallel processing
#' @param alpha Significance level
#' @return List with RRphylo-based rate analysis results
#' @keywords internal
analyze_rate_rphylo <- function(tree, 
                              chr_counts,
                              n_cores = 1,
                              alpha = 0.05) {
  # Check if RRphylo is available
  if(!requireNamespace("RRphylo", quietly = TRUE)) {
    stop("RRphylo package is required for this analysis")
  }
  
  # Run RRphylo analysis
  rr <- RRphylo::RRphylo(tree = tree, y = as.matrix(chr_counts), 
                       cov = NULL, clus = n_cores)
  
  # Search for rate shifts
  search <- RRphylo::search.shift(RR = rr, status.type = "clade", node = NULL, 
                                calibration = NULL, nrep = 100)
  
  # Extract shift points
  if(!is.null(search$single.clades)) {
    shift_points <- search$single.clades$node
    shift_probs <- search$single.clades$probability
    significant_shifts <- shift_probs > (1 - alpha)
  } else {
    shift_points <- integer(0)
    shift_probs <- numeric(0)
    significant_shifts <- logical(0)
  }
  
  # Extract branch-specific rates
  branch_rates <- data.frame(
    node = c(1:length(tree$tip.label), 
             (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)),
    rate = c(rr$rates[, 1], rep(NA, tree$Nnode)),
    stringsAsFactors = FALSE
  )
  
  # Add shift information
  branch_rates$is_shift <- branch_rates$node %in% shift_points
  branch_rates$shift_prob <- NA
  branch_rates$shift_prob[branch_rates$is_shift] <- shift_probs
  
  # Create summary
  n_shifts <- sum(significant_shifts)
  summary_text <- paste0(n_shifts, " significant rate shift(s) detected")
  
  # Create visualization if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggtree", quietly = TRUE)) {
    # Rate plot
    rate_plot <- ggtree::ggtree(tree) %<+% branch_rates +
      ggtree::aes(color = rate) +
      viridis::scale_color_viridis(name = "Evolutionary Rate", na.value = "gray80") +
      ggplot2::labs(title = "Chromosome Evolution Rates (RRphylo)",
                 subtitle = summary_text) +
      ggtree::geom_tiplab() +
      ggplot2::theme_minimal()
    
    # Plot shift points
    if(n_shifts > 0) {
      significant_nodes <- shift_points[significant_shifts]
      rate_plot <- rate_plot +
        ggtree::geom_hilight(node = significant_nodes, alpha = 0.3, extend = 0.5)
    }
    
    plots <- list(rate_plot = rate_plot)
  } else {
    plots <- list()
  }
  
  # Return results
  return(list(
    rates = branch_rates,
    shifts = data.frame(
      node = shift_points,
      probability = shift_probs,
      significant = significant_shifts,
      stringsAsFactors = FALSE
    ),
    summary = list(
      text = summary_text,
      n_shifts = n_shifts,
      shift_nodes = shift_points[significant_shifts]
    ),
    plots = plots,
    rr_results = rr,
    shift_results = search
  ))
}

#===============================================================================
# L1OU Shift Analysis
#===============================================================================

#' Analyze rate variation using L1OU method for OU process
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param min_clade_size Minimum size for a clade to be considered
#' @param alpha Significance level
#' @return List with L1OU-based shift analysis results
#' @keywords internal
analyze_rate_shifts <- function(tree, 
                              chr_counts,
                              min_clade_size = 4,
                              alpha = 0.05) {
  # Check if l1ou is available
  if(!requireNamespace("l1ou", quietly = TRUE)) {
    stop("l1ou package is required for this analysis")
  }
  
  # Convert named vector to matrix with column name "trait"
  trait_mat <- matrix(chr_counts, ncol = 1)
  rownames(trait_mat) <- names(chr_counts)
  colnames(trait_mat) <- "trait"
  
  # Run L1OU analysis
  l1ou_results <- l1ou::estimate_shift_configuration(tree, trait_mat,
                                                  min_subclade_size = min_clade_size)
  
  # Use BIC to determine optimal number of shifts
  bic_results <- l1ou::estimate_convergent_regimes(l1ou_results,
                                                maxCorePerNode = 1,
                                                criterion = "bic")
  
  # Extract shift points
  shift_edges <- l1ou_results$shift.configuration[bic_results$khat]$edges
  
  # Convert shift edges to nodes (shift occurs on edge leading to node)
  shift_nodes <- tree$edge[shift_edges, 2]
  
  # Extract regimes
  regimes <- bic_results$convergent.regimes[[bic_results$khat]]
  
  # Calculate regime-specific parameters
  regime_params <- l1ou::get_regime_optima(l1ou_results,
                                        regimes,
                                        k = bic_results$khat)
  regime_optima <- regime_params$convergent.parameters[, 1]
  
  # Create summary
  n_shifts <- length(shift_nodes)
  n_regimes <- bic_results$khat
  summary_text <- paste0(n_shifts, " shift point(s) detected, with ", 
                      n_regimes, " distinct selective regimes")
  
  # Create visualization if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggtree", quietly = TRUE)) {
    # Create regime mapping
    node_regimes <- rep(1, length(tree$tip.label) + tree$Nnode)
    for(regime in 1:n_regimes) {
      regime_nodes <- which(regimes == regime)
      node_regimes[regime_nodes] <- regime
    }
    
    # Create data frame for visualization
    node_data <- data.frame(
      node = 1:(length(tree$tip.label) + tree$Nnode),
      regime = factor(node_regimes),
      optimum = regime_optima[node_regimes],
      is_shift = 1:(length(tree$tip.label) + tree$Nnode) %in% shift_nodes,
      stringsAsFactors = FALSE
    )
    
    # Create plot
    shift_plot <- ggtree::ggtree(tree) %<+% node_data +
      ggtree::aes(color = regime) +
      ggplot2::labs(title = "Chromosome Evolution Regime Shifts (L1OU)",
                 subtitle = summary_text) +
      ggtree::geom_tiplab() +
      ggplot2::theme_minimal()
    
    # Mark shift points
    if(n_shifts > 0) {
      shift_plot <- shift_plot +
        ggtree::geom_point(data = subset(node_data, is_shift),
                         aes(x = x, y = y), shape = 18, size = 3, color = "red")
    }
    
    plots <- list(shift_plot = shift_plot)
  } else {
    plots <- list()
  }
  
  # Return results
  return(list(
    shifts = data.frame(
      node = shift_nodes,
      edge = shift_edges,
      stringsAsFactors = FALSE
    ),
    regimes = data.frame(
      regime = 1:n_regimes,
      optimum = regime_optima,
      stringsAsFactors = FALSE
    ),
    summary = list(
      text = summary_text,
      n_shifts = n_shifts,
      n_regimes = n_regimes,
      shift_nodes = shift_nodes
    ),
    plots = plots,
    l1ou_results = l1ou_results,
    bic_results = bic_results
  ))
}

#===============================================================================
# Bayesian Shift Analysis
#===============================================================================

#' Analyze rate variation using Bayesian shift detection
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param model Evolutionary model: "BM" or "OU"
#' @param n_simulations Number of MCMC iterations
#' @return List with Bayesian shift analysis results
#' @keywords internal
analyze_rate_bayou <- function(tree, 
                             chr_counts,
                             model = "OU",
                             n_simulations = 5000) {
  # Check if bayou is available
  if(!requireNamespace("bayou", quietly = TRUE)) {
    stop("bayou package is required for this analysis")
  }
  
  # Prepare data
  chr_data <- as.matrix(chr_counts)
  
  # Set up priors based on model
  if(model == "OU") {
    prior <- list(
      alpha = list(
        type = "halfcauchy",
        scale = 0.1
      ),
      sig2 = list(
        type = "halfcauchy",
        scale = 0.1
      ),
      k = list(
        type = "fixed",
        value = 2
      ),
      theta = list(
        type = "normal",
        mean = mean(chr_counts, na.rm = TRUE),
        sd = 2 * sd(chr_counts, na.rm = TRUE)
      ),
      slide = list(
        type = "fixed",
        value = 0
      )
    )
    model_spec <- "OU"
  } else {
    # Default to BM model
    prior <- list(
      sig2 = list(
        type = "halfcauchy",
        scale = 0.1
      ),
      k = list(
        type = "fixed",
        value = 2
      ),
      theta = list(
        type = "normal",
        mean = mean(chr_counts, na.rm = TRUE),
        sd = 2 * sd(chr_counts, na.rm = TRUE)
      )
    )
    model_spec <- "BM"
  }
  
  # Set up MCMC parameters
  mcmc_params <- list(
    ngen = n_simulations,
    samp = 10,
    printfreq = n_simulations/10,
    plot = FALSE
  )
  
  # Run Bayesian MCMC
  bayou_run <- bayou::bayou.mcmc(tree, chr_data, prior = prior,
                               model = model_spec, mcmc = mcmc_params,
                               start = NULL)
  
  # Process results
  shift_results <- bayou::summary.bayouMCMC(bayou_run)
  
  # Extract shift points with highest posterior probabilities
  shift_pp_threshold <- 0.2  # Only consider shifts with PP > 0.2
  if(length(shift_results$branch.posteriors) > 0) {
    shift_branches <- names(shift_results$branch.posteriors[shift_results$branch.posteriors > shift_pp_threshold])
    shift_probs <- shift_results$branch.posteriors[shift_results$branch.posteriors > shift_pp_threshold]
    shift_branches <- as.integer(shift_branches)
  } else {
    shift_branches <- integer(0)
    shift_probs <- numeric(0)
  }
  
  # Convert branch indices to nodes
  if(length(shift_branches) > 0) {
    shift_nodes <- tree$edge[shift_branches, 2]
  } else {
    shift_nodes <- integer(0)
  }
  
  # Create summary
  n_shifts <- length(shift_nodes)
  summary_text <- paste0(n_shifts, " shift(s) detected with posterior probability > ", 
                      shift_pp_threshold)
  
  # Create visualization if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggtree", quietly = TRUE)) {
    # Create data frame for visualization
    node_data <- data.frame(
      node = 1:(length(tree$tip.label) + tree$Nnode),
      is_shift = 1:(length(tree$tip.label) + tree$Nnode) %in% shift_nodes,
      stringsAsFactors = FALSE
    )
    
    # Add shift probabilities
    node_data$pp <- NA
    if(n_shifts > 0) {
      for(i in 1:length(shift_nodes)) {
        node_idx <- which(node_data$node == shift_nodes[i])
        if(length(node_idx) == 1) {
          node_data$pp[node_idx] <- shift_probs[i]
        }
      }
    }
    
    # Create plot
    shift_plot <- ggtree::ggtree(tree) %<+% node_data +
      ggplot2::labs(title = "Bayesian Shift Detection",
                 subtitle = summary_text) +
      ggtree::geom_tiplab() +
      ggplot2::theme_minimal()
    
    # Mark shift points
    if(n_shifts > 0) {
      shift_plot <- shift_plot +
        ggtree::geom_point(data = subset(node_data, is_shift),
                         aes(x = x, y = y, size = pp), 
                         shape = 18, color = "red") +
        ggplot2::scale_size_continuous(name = "Posterior Probability",
                                    range = c(3, 6))
    }
    
    plots <- list(shift_plot = shift_plot)
  } else {
    plots <- list()
  }
  
  # Return results
  return(list(
    shifts = data.frame(
      node = shift_nodes,
      branch = shift_branches,
      probability = shift_probs,
      stringsAsFactors = FALSE
    ),
    summary = list(
      text = summary_text,
      n_shifts = n_shifts,
      shift_nodes = shift_nodes
    ),
    plots = plots,
    bayou_results = bayou_run,
    bayou_summary = shift_results
  ))
}

#===============================================================================
# Consensus and Summary Functions
#===============================================================================

#' Generate consensus of rate variation results across methods
#' 
#' @param rate_results Results from analyze_rate_heterogeneity
#' @return List with consensus results
#' @keywords internal
generate_rate_consensus <- function(rate_results) {
  # Extract shift points from all methods
  all_shifts <- list()
  
  if(!is.null(rate_results$shift_points$rphylo)) {
    rp_shifts <- rate_results$shift_points$rphylo
    if(nrow(rp_shifts) > 0) {
      all_shifts$rphylo <- rp_shifts$node[rp_shifts$significant]
    }
  }
  
  if(!is.null(rate_results$shift_points$l1ou)) {
    l1_shifts <- rate_results$shift_points$l1ou
    if(nrow(l1_shifts) > 0) {
      all_shifts$l1ou <- l1_shifts$node
    }
  }
  
  if(!is.null(rate_results$shift_points$bayou)) {
    by_shifts <- rate_results$shift_points$bayou
    if(nrow(by_shifts) > 0) {
      all_shifts$bayou <- by_shifts$node
    }
  }
  
  # Get all unique shift nodes
  if(length(all_shifts) > 0) {
    all_nodes <- unique(unlist(all_shifts))
    
    # Count how many methods detect each node as a shift
    node_counts <- sapply(all_nodes, function(node) {
      sum(sapply(all_shifts, function(shifts) node %in% shifts))
    })
    names(node_counts) <- all_nodes
    
    # Create consensus data frame
    consensus_shifts <- data.frame(
      node = as.integer(names(node_counts)),
      n_methods = node_counts,
      stringsAsFactors = FALSE
    )
    
    # Add information about which methods detected each shift
    for(method in names(all_shifts)) {
      consensus_shifts[[paste0(method, "_detected")]] <- 
        consensus_shifts$node %in% all_shifts[[method]]
    }
    
    # Order by number of methods
    consensus_shifts <- consensus_shifts[order(-consensus_shifts$n_methods), ]
    
    # Create summary text
    n_consensus <- sum(consensus_shifts$n_methods >= length(all_shifts)/2)
    consensus_text <- paste0(n_consensus, " consensus shift points detected by at least half of the methods")
    
    # Create visualization if ggplot2 is available
    if(requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggtree", quietly = TRUE)) {
      # Create data frame for visualization
      tree <- rate_results$tree
      node_data <- data.frame(
        node = 1:(length(tree$tip.label) + tree$Nnode),
        is_shift = 1:(length(tree$tip.label) + tree$Nnode) %in% consensus_shifts$node,
        n_methods = 0,
        stringsAsFactors = FALSE
      )
      
      # Add method counts
      for(i in 1:nrow(consensus_shifts)) {
        node_idx <- which(node_data$node == consensus_shifts$node[i])
        if(length(node_idx) == 1) {
          node_data$n_methods[node_idx] <- consensus_shifts$n_methods[i]
        }
      }
      
      # Create plot
      consensus_plot <- ggtree::ggtree(tree) %<+% node_data +
        ggplot2::labs(title = "Consensus Shift Detection",
                   subtitle = consensus_text) +
        ggtree::geom_tiplab() +
        ggplot2::theme_minimal()
      
      # Mark shift points
      if(nrow(consensus_shifts) > 0) {
        consensus_plot <- consensus_plot +
          ggtree::geom_point(data = subset(node_data, is_shift),
                           aes(x = x, y = y, size = n_methods, color = factor(n_methods)), 
                           shape = 18) +
          ggplot2::scale_size_continuous(name = "Number of Methods",
                                      range = c(3, 6)) +
          ggplot2::scale_color_viridis_d(name = "Number of Methods")
      }
      
      plots <- list(consensus_plot = consensus_plot)
    } else {
      plots <- list()
    }
    
    # Return consensus results
    return(list(
      shifts = consensus_shifts,
      summary = list(
        text = consensus_text,
        n_consensus = n_consensus,
        consensus_nodes = consensus_shifts$node[consensus_shifts$n_methods >= length(all_shifts)/2]
      ),
      plots = plots
    ))
  } else {
    # No shifts detected by any method
    return(list(
      shifts = data.frame(node = integer(0), n_methods = integer(0)),
      summary = list(
        text = "No shift points detected by any method",
        n_consensus = 0,
        consensus_nodes = integer(0)
      ),
      plots = list()
    ))
  }
}

#' Generate summary of rate analysis results
#' 
#' @param rate_results Results from analyze_rate_heterogeneity
#' @return List with summary information
#' @keywords internal
generate_rate_summary <- function(rate_results) {
  # Compile summary information from all methods
  method_summaries <- list()
  
  if(!is.null(rate_results$summary$contrasts)) {
    method_summaries$contrasts <- rate_results$summary$contrasts$text
  }
  
  if(!is.null(rate_results$summary$rphylo)) {
    method_summaries$rphylo <- rate_results$summary$rphylo$text
  }
  
  if(!is.null(rate_results$summary$l1ou)) {
    method_summaries$l1ou <- rate_results$summary$l1ou$text
  }
  
  if(!is.null(rate_results$summary$bayou)) {
    method_summaries$bayou <- rate_results$summary$bayou$text
  }
  
  # Create consensus summary if available
  if(!is.null(rate_results$consensus)) {
    method_summaries$consensus <- rate_results$consensus$summary$text
  }
  
  # Compile overall summary text
  summary_text <- c(
    "CHROMOSOME EVOLUTION RATE ANALYSIS SUMMARY",
    "==========================================",
    ""
  )
  
  for(method in names(method_summaries)) {
    summary_text <- c(summary_text,
                    paste0(toupper(method), " METHOD:"),
                    method_summaries[[method]],
                    "")
  }
  
  # Create overall conclusion
  has_rate_variation <- FALSE
  
  if(!is.null(rate_results$summary$contrasts) && 
     rate_results$summary$contrasts$n_significant > 0) {
    has_rate_variation <- TRUE
  }
  
  if(!is.null(rate_results$summary$rphylo) && 
     rate_results$summary$rphylo$n_shifts > 0) {
    has_rate_variation <- TRUE
  }
  
  if(!is.null(rate_results$summary$l1ou) && 
     rate_results$summary$l1ou$n_shifts > 0) {
    has_rate_variation <- TRUE
  }
  
  if(!is.null(rate_results$summary$bayou) && 
     rate_results$summary$bayou$n_shifts > 0) {
    has_rate_variation <- TRUE
  }
  
  conclusion <- if(has_rate_variation) {
    "CONCLUSION: Rate heterogeneity detected in chromosome number evolution."
  } else {
    "CONCLUSION: No significant rate heterogeneity detected in chromosome number evolution."
  }
  
  summary_text <- c(summary_text, conclusion)
  
  # Return summary
  return(list(
    text = summary_text,
    method_summaries = method_summaries,
    has_rate_variation = has_rate_variation
  ))
}
