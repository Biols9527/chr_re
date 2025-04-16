#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Simulation & Bootstrap Module
# Author: Bioinformatics Team
# Date: 2023-07-15
# Description: Methods for simulating chromosome evolution, performing bootstrap
#              analyses, and evaluating confidence in ancestral reconstructions
#              through resampling and simulation approaches
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(parallel)
  library(ggplot2)
  library(viridis)
  library(dplyr)
})

#===============================================================================
# Chromosome Evolution Simulation Functions
#===============================================================================

#' Simulate chromosome number evolution along a phylogeny
#' 
#' Simulates chromosome evolution using various models including
#' Brownian motion, discrete rates, and custom event-based models
#' 
#' @param tree Phylogenetic tree
#' @param model Simulation model: "BM", "discrete", "fusion_fission", "custom"
#' @param params List of model parameters
#' @param root_value Chromosome number at the root of the tree
#' @param constrain_values Whether to constrain values to realistic ranges
#' @param min_chr Minimum chromosome number allowed (if constrained)
#' @param max_chr Maximum chromosome number allowed (if constrained)
#' @param discretize Whether to round simulated values to integers
#' @param seed Random seed for reproducibility
#' @return Named vector of simulated chromosome numbers for all nodes
#' @export
simulate_chromosome_evolution <- function(tree, 
                                        model = "BM", 
                                        params = list(sig2 = 0.1), 
                                        root_value = 10,
                                        constrain_values = TRUE,
                                        min_chr = 1,
                                        max_chr = 100,
                                        discretize = TRUE,
                                        seed = NULL) {
  # Set random seed if provided
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validate tree
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Validate model
  supported_models <- c("BM", "discrete", "fusion_fission", "custom")
  if(!model %in% supported_models) {
    stop(paste("Unsupported model. Please use one of:", paste(supported_models, collapse = ", ")))
  }
  
  # Simulate based on model
  sim_values <- NULL
  
  if(model == "BM") {
    # Check required parameters
    if(!("sig2" %in% names(params))) {
      params$sig2 <- 0.1  # Default variance rate
    }
    
    # Simulate using Brownian Motion
    sim_values <- phytools::fastBM(tree, 
                                 sig2 = params$sig2, 
                                 a = root_value,
                                 internal = TRUE)
  }
  else if(model == "discrete") {
    # Discrete transition model
    # Extract or set default parameters
    fusion_rate <- if("fusion_rate" %in% names(params)) params$fusion_rate else 0.05
    fission_rate <- if("fission_rate" %in% names(params)) params$fission_rate else 0.05
    
    # Simulate using discrete model
    sim_values <- simulate_discrete_changes(tree, 
                                          root_value = root_value,
                                          fusion_rate = fusion_rate,
                                          fission_rate = fission_rate)
  }
  else if(model == "fusion_fission") {
    # Custom fusion-fission model with possible WGD events
    # Extract or set default parameters
    fusion_rate <- if("fusion_rate" %in% names(params)) params$fusion_rate else 0.05
    fission_rate <- if("fission_rate" %in% names(params)) params$fission_rate else 0.05
    wgd_rate <- if("wgd_rate" %in% names(params)) params$wgd_rate else 0.01
    
    # Simulate using fusion-fission model
    sim_values <- simulate_fusion_fission_events(tree, 
                                               root_value = root_value,
                                               fusion_rate = fusion_rate,
                                               fission_rate = fission_rate,
                                               wgd_rate = wgd_rate)
  }
  else if(model == "custom") {
    # Custom user-provided simulation function
    if(!("sim_func" %in% names(params))) {
      stop("When using 'custom' model, params must include 'sim_func' function")
    }
    
    custom_func <- params$sim_func
    sim_values <- custom_func(tree, root_value, params)
  }
  
  # Round to integers if requested
  if(discretize) {
    sim_values <- round(sim_values)
  }
  
  # Constrain values if requested
  if(constrain_values) {
    sim_values[sim_values < min_chr] <- min_chr
    sim_values[sim_values > max_chr] <- max_chr
  }
  
  return(sim_values)
}

#' Simulate chromosome evolution using discrete changes (fusion, fission)
#' 
#' @param tree Phylogenetic tree
#' @param root_value Starting chromosome number at root
#' @param fusion_rate Rate at which fusion events occur per unit branch length
#' @param fission_rate Rate at which fission events occur per unit branch length
#' @return Named vector of chromosome numbers for all nodes
#' @keywords internal
simulate_discrete_changes <- function(tree, root_value, fusion_rate, fission_rate) {
  # Check that tree has branch lengths
  if(is.null(tree$edge.length)) {
    stop("Tree must have branch lengths for discrete simulation")
  }
  
  # Get number of nodes
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  
  # Initialize chromosome counts
  chr_counts <- numeric(n_nodes)
  names(chr_counts) <- 1:n_nodes
  
  # Set root value
  root_node <- n_tips + 1  # Assuming root is first internal node
  chr_counts[root_node] <- root_value
  
  # Create post-order traversal sequence
  postorder <- rev(postorder.edges(tree))
  
  # Process each edge in post-order
  for(i in postorder) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent chromosome count
    parent_chr <- chr_counts[parent]
    
    # Calculate number of events on this branch
    fusion_events <- rpois(1, fusion_rate * branch_length)
    fission_events <- rpois(1, fission_rate * branch_length)
    
    # Calculate child chromosome count
    child_chr <- parent_chr - fusion_events + fission_events
    
    # Ensure count is at least 1
    child_chr <- max(child_chr, 1)
    
    # Store child chromosome count
    chr_counts[child] <- child_chr
  }
  
  return(chr_counts)
}

#' Simulate chromosome evolution with fusion, fission, and WGD events
#' 
#' @param tree Phylogenetic tree
#' @param root_value Starting chromosome number at root
#' @param fusion_rate Rate at which fusion events occur per unit branch length
#' @param fission_rate Rate at which fission events occur per unit branch length
#' @param wgd_rate Rate at which whole genome duplication events occur
#' @return Named vector of chromosome numbers for all nodes
#' @keywords internal
simulate_fusion_fission_events <- function(tree, 
                                         root_value, 
                                         fusion_rate, 
                                         fission_rate,
                                         wgd_rate) {
  # Check that tree has branch lengths
  if(is.null(tree$edge.length)) {
    stop("Tree must have branch lengths for discrete simulation")
  }
  
  # Get number of nodes
  n_tips <- length(tree$tip.label)
  n_nodes <- n_tips + tree$Nnode
  
  # Initialize chromosome counts
  chr_counts <- numeric(n_nodes)
  names(chr_counts) <- 1:n_nodes
  
  # Set root value
  root_node <- n_tips + 1  # Assuming root is first internal node
  chr_counts[root_node] <- root_value
  
  # Create post-order traversal sequence
  postorder <- rev(postorder.edges(tree))
  
  # Initialize event tracking
  events <- list()
  
  # Process each edge in post-order
  for(i in postorder) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent chromosome count
    parent_chr <- chr_counts[parent]
    
    # Initialize child chromosome count
    child_chr <- parent_chr
    
    # Track events for this branch
    branch_events <- list(
      fusion = 0,
      fission = 0,
      wgd = 0
    )
    
    # Simulate events along this branch
    remaining_length <- branch_length
    current_chr <- parent_chr
    
    while(remaining_length > 0) {
      # Calculate event rates
      fusion_prob <- fusion_rate * current_chr  # Higher for more chromosomes
      fission_prob <- fission_rate * current_chr  # Higher for more chromosomes
      wgd_prob <- wgd_rate
      
      # Total event rate
      total_rate <- fusion_prob + fission_prob + wgd_prob
      
      # Generate time until next event
      time_to_event <- rexp(1, rate = total_rate)
      
      # If no more events along this branch
      if(time_to_event > remaining_length) {
        break
      }
      
      # Update remaining branch length
      remaining_length <- remaining_length - time_to_event
      
      # Determine which event occurred
      event_probs <- c(fusion_prob, fission_prob, wgd_prob) / total_rate
      event_type <- sample(c("fusion", "fission", "wgd"), 1, prob = event_probs)
      
      # Apply the event
      if(event_type == "fusion" && current_chr > 1) {
        current_chr <- current_chr - 1
        branch_events$fusion <- branch_events$fusion + 1
      } else if(event_type == "fission") {
        current_chr <- current_chr + 1
        branch_events$fission <- branch_events$fission + 1
      } else if(event_type == "wgd") {
        current_chr <- current_chr * 2
        branch_events$wgd <- branch_events$wgd + 1
      }
    }
    
    # Update child chromosome count
    chr_counts[child] <- current_chr
    
    # Store events
    events[[i]] <- branch_events
  }
  
  # Add events as an attribute
  attr(chr_counts, "events") <- events
  
  return(chr_counts)
}

#' Generate a parametric bootstrap sample for chromosome count data
#' 
#' Creates a bootstrap sample by simulating chromosome evolution under
#' a fitted model, preserving the original tree topology
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Original chromosome counts (for tips)
#' @param model Model for simulation: "BM", "discrete", "fusion_fission"
#' @param fitted_params Parameters from fitted model (if NULL, estimated from data)
#' @param root_value Chromosome number at the root (if NULL, estimated from data)
#' @param constrain_values Whether to constrain values to realistic ranges
#' @return Named vector of simulated chromosome counts for tree tips
#' @export
generate_bootstrap_sample <- function(tree, 
                                    chr_counts,
                                    model = "BM",
                                    fitted_params = NULL,
                                    root_value = NULL,
                                    constrain_values = TRUE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Match data to tree
  matched_counts <- match_data_to_tree(tree, chr_counts)
  
  # Estimate parameters if not provided
  if(is.null(fitted_params)) {
    fitted_params <- estimate_model_parameters(tree, matched_counts, model)
  }
  
  # Estimate root value if not provided
  if(is.null(root_value)) {
    if(model == "BM") {
      # For BM, use ACE to estimate ancestral state
      ace_result <- ace(matched_counts, tree, type = "continuous")
      root_node <- length(matched_counts) + 1
      root_value <- ace_result$ace[1]  # First node is root
    } else {
      # For discrete models, use median
      root_value <- median(matched_counts, na.rm = TRUE)
    }
  }
  
  # Simulate new data
  simulated_values <- simulate_chromosome_evolution(
    tree = tree,
    model = model,
    params = fitted_params,
    root_value = root_value,
    constrain_values = constrain_values
  )
  
  # Extract only tip values
  tip_indices <- 1:length(tree$tip.label)
  bootstrap_sample <- simulated_values[tip_indices]
  names(bootstrap_sample) <- tree$tip.label
  
  return(bootstrap_sample)
}

#' Estimate parameters for chromosome evolution models
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param model Model to fit: "BM", "discrete", "fusion_fission"
#' @return List of estimated parameters for the model
#' @keywords internal
estimate_model_parameters <- function(tree, chr_counts, model) {
  params <- list()
  
  if(model == "BM") {
    # Estimate Brownian motion rate parameter
    # Use maximum likelihood estimation from geiger package
    if(requireNamespace("geiger", quietly = TRUE)) {
      bm_fit <- geiger::fitContinuous(tree, chr_counts, model = "BM")
      params$sig2 <- bm_fit$opt$sigsq
    } else {
      # Simple variance-based estimate if geiger not available
      pic_values <- pic(chr_counts, tree)
      params$sig2 <- var(pic_values) / 2
    }
  }
  else if(model == "discrete" || model == "fusion_fission") {
    # Estimate discrete model parameters
    # Use ancestral state reconstruction to infer events
    ancestral_states <- estimate_ancestral_states(tree, chr_counts)
    
    # Calculate rates by counting events
    n_fusion <- 0
    n_fission <- 0
    n_wgd <- 0
    total_branch_length <- sum(tree$edge.length)
    
    for(i in 1:nrow(tree$edge)) {
      parent <- tree$edge[i, 1]
      child <- tree$edge[i, 2]
      
      parent_chr <- ancestral_states[as.character(parent)]
      child_chr <- ancestral_states[as.character(child)]
      
      if(is.na(parent_chr) || is.na(child_chr)) next
      
      # Count events
      if(child_chr < parent_chr) {
        # Fusion events
        n_fusion <- n_fusion + (parent_chr - child_chr)
      } else if(child_chr > parent_chr) {
        if(child_chr >= parent_chr * 1.8) {
          # Potential WGD
          n_wgd <- n_wgd + 1
        } else {
          # Fission events
          n_fission <- n_fission + (child_chr - parent_chr)
        }
      }
    }
    
    # Calculate rates
    params$fusion_rate <- n_fusion / total_branch_length
    params$fission_rate <- n_fission / total_branch_length
    
    if(model == "fusion_fission") {
      params$wgd_rate <- n_wgd / total_branch_length
    }
  }
  
  return(params)
}

#' Estimate ancestral chromosome numbers for parameter estimation
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @return Named vector of ancestral state estimates
#' @keywords internal
estimate_ancestral_states <- function(tree, chr_counts) {
  # Use ML estimation for this purpose
  ace_result <- ace(chr_counts, tree, type = "continuous")
  
  # Extract ancestral state estimates
  n_tips <- length(tree$tip.label)
  internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  
  # Combine tip and node states
  all_states <- c(chr_counts, ace_result$ace)
  names(all_states) <- c(1:n_tips, internal_nodes)
  
  return(all_states)
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
# Bootstrap Analysis Functions
#===============================================================================

#' Perform parametric bootstrap analysis of ancestral chromosome reconstructions
#' 
#' Uses simulation-based bootstrapping to estimate confidence in
#' reconstructed ancestral chromosome numbers
#'
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts for extant species
#' @param reconstruction_method Method for ancestral reconstruction: "parsimony", "ML", "Bayesian"
#' @param simulation_model Model for bootstrap simulation: "BM", "discrete", "fusion_fission"
#' @param n_bootstraps Number of bootstrap replicates
#' @param fitted_params Optional model parameters (estimated from data if NULL)
#' @param nodes Vector of nodes to calculate confidence for (NULL for all internal nodes)
#' @param ci_level Confidence interval level (e.g., 0.95 for 95% CI)
#' @param n_cores Number of cores for parallel processing (NULL = auto-detect)
#' @return List with bootstrap results and confidence intervals
#' @export
parametric_bootstrap_analysis <- function(tree,
                                        chr_counts,
                                        reconstruction_method = "ML",
                                        simulation_model = "BM",
                                        n_bootstraps = 100,
                                        fitted_params = NULL,
                                        nodes = NULL,
                                        ci_level = 0.95,
                                        n_cores = NULL) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Match data to tree
  matched_counts <- match_data_to_tree(tree, chr_counts)
  
  # Set up nodes to analyze
  n_tips <- length(tree$tip.label)
  if(is.null(nodes)) {
    # Default to all internal nodes
    nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  }
  
  # Perform original reconstruction
  original_recon <- reconstruct_ancestral_states(tree, matched_counts, reconstruction_method)
  
  # Estimate model parameters if not provided
  if(is.null(fitted_params)) {
    fitted_params <- estimate_model_parameters(tree, matched_counts, simulation_model)
  }
  
  # Estimate root value
  root_node <- n_tips + 1
  root_value <- original_recon[as.character(root_node)]
  if(is.na(root_value)) {
    root_value <- median(matched_counts, na.rm = TRUE)
  }
  
  # Set up parallel processing if requested
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Function to run a single bootstrap replicate
  run_bootstrap_replicate <- function(i) {
    # Generate bootstrap sample
    bootstrap_sample <- generate_bootstrap_sample(
      tree = tree,
      chr_counts = matched_counts,
      model = simulation_model,
      fitted_params = fitted_params,
      root_value = root_value
    )
    
    # Reconstruct ancestral states for the bootstrap sample
    bootstrap_recon <- reconstruct_ancestral_states(
      tree = tree,
      chr_counts = bootstrap_sample,
      method = reconstruction_method
    )
    
    # Extract states for nodes of interest
    node_states <- bootstrap_recon[as.character(nodes)]
    
    return(node_states)
  }
  
  # Run bootstrap replicates (in parallel if requested)
  if(n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    bootstrap_results <- parallel::parLapply(cl, 1:n_bootstraps, function(i) {
      run_bootstrap_replicate(i)
    })
  } else {
    bootstrap_results <- lapply(1:n_bootstraps, run_bootstrap_replicate)
  }
  
  # Organize results by node
  node_bootstraps <- list()
  
  for(node in nodes) {
    node_str <- as.character(node)
    
    # Extract bootstrap values for this node
    node_values <- sapply(bootstrap_results, function(x) x[node_str])
    
    # Calculate confidence interval
    alpha <- 1 - ci_level
    ci_lower <- quantile(node_values, alpha/2, na.rm = TRUE)
    ci_upper <- quantile(node_values, 1 - alpha/2, na.rm = TRUE)
    
    # Store node results
    node_bootstraps[[node_str]] <- list(
      node_id = node,
      original_value = original_recon[node_str],
      bootstrap_values = node_values,
      mean = mean(node_values, na.rm = TRUE),
      median = median(node_values, na.rm = TRUE),
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_level = ci_level
    )
  }
  
  # Prepare final results
  results <- list(
    tree = tree,
    original_reconstruction = original_recon,
    reconstruction_method = reconstruction_method,
    simulation_model = simulation_model,
    n_bootstraps = n_bootstraps,
    fitted_params = fitted_params,
    node_results = node_bootstraps,
    summary = list(
      n_nodes = length(nodes),
      ci_level = ci_level,
      mean_ci_width = mean(sapply(node_bootstraps, function(x) x$ci_upper - x$ci_lower)),
      reconstruction_method = reconstruction_method,
      simulation_model = simulation_model
    )
  )
  
  class(results) <- c("chr_bootstrap", class(results))
  
  return(results)
}

#' Reconstruct ancestral states using specified method
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param method Reconstruction method: "parsimony", "ML", "Bayesian"
#' @return Named vector of ancestral state estimates
#' @keywords internal
reconstruct_ancestral_states <- function(tree, chr_counts, method) {
  # Different reconstruction methods
  if(method == "parsimony") {
    # Parsimony reconstruction
    if(requireNamespace("phangorn", quietly = TRUE)) {
      # Create phyDat object for parsimony
      chr_matrix <- matrix(as.character(chr_counts[tree$tip.label]), 
                         nrow = length(tree$tip.label), ncol = 1)
      rownames(chr_matrix) <- tree$tip.label
      
      # Get range of chromosome counts
      all_counts <- unique(chr_counts[!is.na(chr_counts)])
      min_chr <- min(all_counts, na.rm = TRUE)
      max_chr <- max(all_counts, na.rm = TRUE)
      
      # Create phyDat object
      phydat <- phangorn::phyDat(chr_matrix, type = "USER", 
                               levels = as.character(min_chr:max_chr))
      
      # Perform ancestral reconstruction
      anc <- phangorn::ancestral.pars(tree, phydat)
      
      # Extract most likely state for each node
      n_tips <- length(tree$tip.label)
      n_nodes <- n_tips + tree$Nnode
      
      # Get node states
      node_states <- apply(anc, 1, function(x) {
        as.numeric(names(which.max(x)))
      })
      
      # Combine tip and node states
      all_states <- c(chr_counts[tree$tip.label], node_states)
      names(all_states) <- c(1:n_tips, (n_tips + 1):n_nodes)
      
      return(all_states)
    } else {
      stop("phangorn package required for parsimony reconstruction")
    }
  } else if(method == "ML") {
    # ML reconstruction
    ace_result <- ace(chr_counts, tree, type = "continuous")
    
    # Extract ancestral state estimates
    n_tips <- length(tree$tip.label)
    internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
    
    # Combine tip and node states
    all_states <- c(chr_counts, ace_result$ace)
    names(all_states) <- c(1:n_tips, internal_nodes)
    
    return(all_states)
  } else if(method == "Bayesian") {
    # Bayesian reconstruction
    if(requireNamespace("phytools", quietly = TRUE)) {
      mcmc_samples <- phytools::anc.Bayes(tree, chr_counts, ngen = 10000, 
                                        control = list(sample.freq = 10))
      
      # Calculate mean ancestral states
      n_tips <- length(tree$tip.label)
      internal_nodes <- (n_tips + 1):(n_tips + tree$Nnode)
      
      # Extract posterior means
      post_means <- colMeans(mcmc_samples$mcmc[, -c(1:2)])
      
      # Combine tip and node states
      all_states <- c(chr_counts, post_means)
      names(all_states) <- c(1:n_tips, internal_nodes)
      
      return(all_states)
    } else {
      stop("phytools package required for Bayesian reconstruction")
    }
  } else {
    stop(paste("Unsupported reconstruction method:", method))
  }
}

#===============================================================================
# Visualization Functions for Bootstrap Results
#===============================================================================

#' Plot bootstrap confidence intervals for ancestral chromosome numbers
#' 
#' Creates a plot showing the original reconstructed values with bootstrap
#' confidence intervals for selected nodes
#' 
#' @param bootstrap_results Results from parametric_bootstrap_analysis
#' @param nodes Vector of nodes to plot (NULL for all nodes in results)
#' @param sort_by How to sort nodes: "node_id", "value", "interval_width"
#' @param decreasing Whether to sort in decreasing order
#' @param color_scheme Color scheme for plot elements
#' @param show_node_labels Whether to show node IDs on plot
#' @param add_tree Whether to add a small tree visualization
#' @return ggplot object with confidence interval visualization
#' @export
plot_bootstrap_intervals <- function(bootstrap_results,
                                   nodes = NULL,
                                   sort_by = "node_id",
                                   decreasing = FALSE,
                                   color_scheme = "viridis",
                                   show_node_labels = TRUE,
                                   add_tree = FALSE) {
  # Validate input
  if(!inherits(bootstrap_results, "chr_bootstrap")) {
    stop("bootstrap_results must be from parametric_bootstrap_analysis function")
  }
  
  # Extract node results
  node_results <- bootstrap_results$node_results
  
  # Filter to specified nodes if provided
  if(!is.null(nodes)) {
    node_results <- node_results[as.character(nodes)]
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    node_id = sapply(node_results, function(x) x$node_id),
    original = sapply(node_results, function(x) x$original_value),
    mean = sapply(node_results, function(x) x$mean),
    median = sapply(node_results, function(x) x$median),
    lower = sapply(node_results, function(x) x$ci_lower),
    upper = sapply(node_results, function(x) x$ci_upper),
    width = sapply(node_results, function(x) x$ci_upper - x$ci_lower),
    stringsAsFactors = FALSE
  )
  
  # Sort data
  if(sort_by == "node_id") {
    plot_data <- plot_data[order(plot_data$node_id, decreasing = decreasing), ]
  } else if(sort_by == "value") {
    plot_data <- plot_data[order(plot_data$original, decreasing = decreasing), ]
  } else if(sort_by == "interval_width") {
    plot_data <- plot_data[order(plot_data$width, decreasing = decreasing), ]
  }
  
  # Add node factor for plotting order
  plot_data$node_factor <- factor(plot_data$node_id, 
                                levels = plot_data$node_id[order(nrow(plot_data):1)])
  
  # Set up colors based on color scheme
  if(color_scheme == "viridis") {
    point_color <- viridis(1, begin = 0.2)[1]
    interval_color <- viridis(1, begin = 0.7)[1]
    line_color <- viridis(1, begin = 0.5)[1]
  } else {
    # Default colors
    point_color <- "steelblue"
    interval_color <- "skyblue"
    line_color <- "darkblue"
  }
  
  # Create plot
  p <- ggplot(plot_data, aes(x = node_factor, y = original)) +
    # Add confidence intervals
    geom_errorbar(aes(ymin = lower, ymax = upper), 
                width = 0.3, color = interval_color, size = 1.2) +
    # Add original values
    geom_point(size = 3, color = point_color) +
    # Add bootstrap mean/median as line
    geom_point(aes(y = median), shape = 3, size = 2, color = line_color) +
    # Customize appearance
    coord_flip() +
    labs(
      title = "Bootstrap Confidence Intervals for Ancestral Chromosome Numbers",
      subtitle = paste0(bootstrap_results$ci_level * 100, 
                      "% confidence intervals based on ", 
                      bootstrap_results$n_bootstraps, " bootstrap replicates"),
      x = "Node ID",
      y = "Chromosome Number"
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "gray80"),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  
  # Add node labels if requested
  if(show_node_labels) {
    p <- p + geom_text(aes(label = node_id), hjust = -0.5, size = 3)
  }
  
  # Add tree if requested
  if(add_tree && requireNamespace("ggtree", quietly = TRUE)) {
    tree <- bootstrap_results$tree
    
    # Create tree plot
    tree_plot <- ggtree::ggtree(tree) + 
      ggtree::geom_nodelab(aes(label = node)) + 
      ggtree::theme_tree2()
    
    # Combine plots using patchwork if available
    if(requireNamespace("patchwork", quietly = TRUE)) {
      p <- patchwork::wrap_plots(
        tree_plot, p, 
        widths = c(1, 2)
      ) + 
        patchwork::plot_annotation(
          title = "Ancestral Chromosome Number Bootstrap Analysis",
          theme = theme(plot.title = element_text(hjust = 0.5))
        )
    }
  }
  
  return(p)
}

#' Create density plots of bootstrap distributions for ancestral nodes
#' 
#' Visualizes the distributions of bootstrap estimates for selected nodes
#' 
#' @param bootstrap_results Results from parametric_bootstrap_analysis
#' @param nodes Vector of node IDs to visualize (NULL for all)
#' @param arrange_in_grid Whether to arrange plots in a grid
#' @param ncol Number of columns for grid arrangement
#' @param color_by_uncertainty Whether to color distributions by uncertainty width
#' @param add_original_values Whether to show original reconstruction values
#' @param add_ci Whether to show confidence interval
#' @return List of ggplot objects or a combined grid plot
#' @export
plot_bootstrap_distributions <- function(bootstrap_results,
                                       nodes = NULL,
                                       arrange_in_grid = TRUE,
                                       ncol = 3,
                                       color_by_uncertainty = TRUE,
                                       add_original_values = TRUE,
                                       add_ci = TRUE) {
  # Validate input
  if(!inherits(bootstrap_results, "chr_bootstrap")) {
    stop("bootstrap_results must be from parametric_bootstrap_analysis function")
  }
  
  # Extract node results
  node_results <- bootstrap_results$node_results
  
  # Filter to specified nodes if provided
  if(!is.null(nodes)) {
    node_results <- node_results[as.character(nodes)]
  }
  
  # Create a list to store individual plots
  plots <- list()
  
  # Create density plots for each node
  for(node_id in names(node_results)) {
    node_data <- node_results[[node_id]]
    
    # Create data frame for plotting
    plot_data <- data.frame(
      chromosome_number = node_data$bootstrap_values,
      stringsAsFactors = FALSE
    )
    
    # Calculate uncertainty (width of confidence interval)
    uncertainty <- node_data$ci_upper - node_data$ci_lower
    
    # Create base plot
    p <- ggplot(plot_data, aes(x = chromosome_number)) +
      # Add density
      geom_density(fill = ifelse(color_by_uncertainty, 
                               viridis(1, alpha = 0.7, begin = max(0, 1 - uncertainty/5))[1], 
                               "steelblue"), 
                 alpha = 0.7) +
      # Customize appearance
      labs(
        title = paste("Node", node_id),
        subtitle = paste("Mean:", round(node_data$mean, 2),
                       "CI:", round(node_data$ci_lower, 2), "to", 
                       round(node_data$ci_upper, 2)),
        x = "Chromosome Number",
        y = "Density"
      ) +
      theme_minimal() +
      theme(
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "gray80")
      )
    
    # Add original value if requested
    if(add_original_values) {
      p <- p + geom_vline(xintercept = node_data$original_value, 
                        color = "red", linetype = "dashed", size = 1) +
        annotate("text", x = node_data$original_value, y = 0, 
               label = "Original", color = "red", angle = 90, 
               vjust = -0.5, hjust = -0.2, size = 3)
    }
    
    # Add confidence interval if requested
    if(add_ci) {
      p <- p + 
        geom_vline(xintercept = node_data$ci_lower, 
                 color = "darkgreen", linetype = "dotted") +
        geom_vline(xintercept = node_data$ci_upper, 
                 color = "darkgreen", linetype = "dotted") +
        annotate("rect", xmin = node_data$ci_lower, xmax = node_data$ci_upper, 
               ymin = 0, ymax = Inf, fill = "darkgreen", alpha = 0.1)
    }
    
    plots[[node_id]] <- p
  }
  
  # Arrange in grid if requested
  if(arrange_in_grid && length(plots) > 1) {
    if(requireNamespace("patchwork", quietly = TRUE)) {
      # Use patchwork to combine plots
      combined_plot <- patchwork::wrap_plots(plots, ncol = ncol)
      return(combined_plot)
    } else if(requireNamespace("gridExtra", quietly = TRUE)) {
      # Use gridExtra as fallback
      combined_plot <- gridExtra::grid.arrange(grobs = plots, ncol = ncol)
      return(combined_plot)
    } else {
      warning("patchwork or gridExtra package required for grid arrangement. Returning list of plots.")
    }
  }
  
  # Return list of individual plots
  return(plots)
}

#' Plot chromosome evolution simulations on a phylogenetic tree
#' 
#' Visualizes simulated chromosome evolution for verification and testing
#' 
#' @param tree Phylogenetic tree
#' @param sim_values Simulated chromosome numbers from simulate_chromosome_evolution
#' @param layout Tree layout: "rectangular", "circular", "fan"
#' @param show_labels Whether to show node labels with chromosome numbers
#' @param edge_width Edge width for tree branches
#' @param node_size Size of nodes proportional to chromosome number
#' @param color_scheme Color scheme for chromosome number mapping
#' @param highlight_events Whether to highlight events (requires events attribute)
#' @return ggplot object with tree visualization
#' @export
plot_simulated_evolution <- function(tree,
                                   sim_values,
                                   layout = "rectangular",
                                   show_labels = TRUE,
                                   edge_width = 1,
                                   node_size = 3,
                                   color_scheme = "viridis",
                                   highlight_events = TRUE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggtree package required for tree visualization")
  }
  
  # Create data frame for plotting
  node_data <- data.frame(
    node = as.numeric(names(sim_values)),
    chr_num = as.numeric(sim_values),
    stringsAsFactors = FALSE
  )
  
  # Create tree plot with mapped chromosome numbers
  if(layout == "circular") {
    p <- ggtree::ggtree(tree, layout = "circular", size = edge_width)
  } else if(layout == "fan") {
    p <- ggtree::ggtree(tree, layout = "fan", size = edge_width)
  } else {
    p <- ggtree::ggtree(tree, size = edge_width)
  }
  
  # Add chromosome number coloring
  p <- p %<+% node_data + aes(color = chr_num)
  
  # Apply color scheme
  if(color_scheme == "viridis") {
    p <- p + scale_color_viridis_c(name = "Chromosome\nNumber")
  } else if(color_scheme == "magma") {
    p <- p + scale_color_viridis_c(option = "magma", name = "Chromosome\nNumber")
  } else {
    p <- p + scale_color_continuous(name = "Chromosome\nNumber")
  }
  
  # Add node points
  p <- p + ggtree::geom_nodepoint(aes(size = chr_num), alpha = 0.8)
  
  # Scale point sizes
  p <- p + scale_size_continuous(range = c(1, node_size * 2), guide = "none")
  
  # Add tip labels
  if(show_labels) {
    p <- p + ggtree::geom_tiplab(size = 3)
  }
  
  # Add node labels with chromosome numbers
  if(show_labels) {
    p <- p + ggtree::geom_text(aes(label = chr_num), 
                            color = "black", hjust = -0.3, size = 3)
  }
  
  # Highlight events if requested and available
  if(highlight_events && !is.null(attr(sim_values, "events"))) {
    events <- attr(sim_values, "events")
    
    # Create data for highlighting branches with events
    event_data <- data.frame(
      edge = as.numeric(names(events)),
      fusion = sapply(events, function(e) e$fusion > 0),
      fission = sapply(events, function(e) e$fission > 0),
      wgd = sapply(events, function(e) e$wgd > 0),
      stringsAsFactors = FALSE
    )
    
    # Extract edge information
    edge_data <- as.data.frame(tree$edge)
    names(edge_data) <- c("parent", "node")
    edge_data$edge <- 1:nrow(edge_data)
    
    # Add event information to edge data
    edge_data <- merge(edge_data, event_data, by = "edge", all.x = TRUE)
    
    # Extract node positions to add event markers
    node_positions <- ggplot2::ggplot_build(p)$data[[1]]
    
    # Add event markers to the plot (midpoint of branches)
    for(event_type in c("fusion", "fission", "wgd")) {
      if(any(edge_data[[event_type]], na.rm = TRUE)) {
        event_edges <- edge_data[edge_data[[event_type]], ]
        
        for(i in 1:nrow(event_edges)) {
          edge_idx <- event_edges$edge[i]
          parent <- event_edges$parent[i]
          node <- event_edges$node[i]
          
          # Get parent and node positions
          parent_pos <- node_positions[node_positions$node == parent, c("x", "y")]
          node_pos <- node_positions[node_positions$node == node, c("x", "y")]
          
          # Calculate midpoint
          mid_x <- (parent_pos$x + node_pos$x) / 2
          mid_y <- (parent_pos$y + node_pos$y) / 2
          
          # Get event count
          event_count <- events[[as.character(edge_idx)]][[event_type]]
          
          if(event_count > 0) {
            if(event_type == "fusion") {
              p <- p + geom_point(aes(x = mid_x, y = mid_y), 
                               color = "blue", shape = 25, size = 3, data = data.frame(mid_x, mid_y))
            } else if(event_type == "fission") {
              p <- p + geom_point(aes(x = mid_x, y = mid_y), 
                               color = "red", shape = 24, size = 3, data = data.frame(mid_x, mid_y))
            } else if(event_type == "wgd") {
              p <- p + geom_point(aes(x = mid_x, y = mid_y), 
                               color = "purple", shape = 23, size = 3, data = data.frame(mid_x, mid_y))
            }
          }
        }
      }
    }
  }
  
  # Add title and theme
  p <- p + labs(title = "Simulated Chromosome Evolution") +
    theme(legend.position = "right")
  
  return(p)
}

#===============================================================================
# Non-parametric Bootstrap Functions
#===============================================================================

#' Perform non-parametric bootstrap analysis for ancestral reconstruction
#' 
#' Uses empirical bootstrapping by sampling species with replacement to
#' assess confidence in ancestral state reconstructions
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param reconstruction_method Method for ancestral reconstruction: "parsimony", "ML", "Bayesian"
#' @param n_bootstraps Number of bootstrap replicates
#' @param nodes Vector of nodes to calculate confidence for (NULL for all)
#' @param ci_level Confidence interval level (e.g., 0.95 for 95% CI)
#' @param n_cores Number of cores for parallel processing (NULL = auto-detect)
#' @return List with bootstrap results and confidence intervals
#' @export
nonparametric_bootstrap_analysis <- function(tree,
                                           chr_counts,
                                           reconstruction_method = "ML",
                                           n_bootstraps = 100,
                                           nodes = NULL,
                                           ci_level = 0.95,
                                           n_cores = NULL) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Match data to tree
  matched_counts <- match_data_to_tree(tree, chr_counts)
  
  # Set up nodes to analyze
  n_tips <- length(tree$tip.label)
  if(is.null(nodes)) {
    # Default to all internal nodes
    nodes <- (n_tips + 1):(n_tips + tree$Nnode)
  }
  
  # Perform original reconstruction
  original_recon <- reconstruct_ancestral_states(tree, matched_counts, reconstruction_method)
  
  # Set up parallel processing if requested
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Function to run a single bootstrap replicate
  run_bootstrap_replicate <- function(i) {
    # Create a bootstrap sample by resampling tips with replacement
    bootstrap_tips <- sample(tree$tip.label, size = length(tree$tip.label), replace = TRUE)
    bootstrap_counts <- matched_counts[bootstrap_tips]
    
    # Since we're sampling with replacement, we need to handle duplicated tips
    # We'll create a modified dataset with unique tip names
    unique_tips <- make.unique(bootstrap_tips)
    names(bootstrap_counts) <- unique_tips
    
    # Create a modified tree with the bootstrapped tips
    bootstrap_tree <- tree
    bootstrap_tree$tip.label <- unique_tips
    
    # Reconstruct ancestral states for the bootstrap sample
    bootstrap_recon <- tryCatch({
      reconstruct_ancestral_states(bootstrap_tree, bootstrap_counts, reconstruction_method)
    }, error = function(e) {
      # If reconstruction fails, return NULL
      warning(paste("Bootstrap replicate", i, "failed:", e$message))
      return(NULL)
    })
    
    # If reconstruction succeeded, extract states for nodes of interest
    if(!is.null(bootstrap_recon)) {
      node_states <- bootstrap_recon[as.character(nodes)]
      return(node_states)
    } else {
      return(NULL)
    }
  }
  
  # Run bootstrap replicates (in parallel if requested)
  if(n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    bootstrap_results <- parallel::parLapply(cl, 1:n_bootstraps, function(i) {
      run_bootstrap_replicate(i)
    })
  } else {
    bootstrap_results <- lapply(1:n_bootstraps, run_bootstrap_replicate)
  }
  
  # Remove failed replicates
  bootstrap_results <- bootstrap_results[!sapply(bootstrap_results, is.null)]
  
  # If no successful replicates, return error
  if(length(bootstrap_results) == 0) {
    stop("All bootstrap replicates failed. Try a different reconstruction method or check your data.")
  }
  
  # Organize results by node
  node_bootstraps <- list()
  
  for(node in nodes) {
    node_str <- as.character(node)
    
    # Extract bootstrap values for this node
    node_values <- sapply(bootstrap_results, function(x) x[node_str])
    
    # Calculate confidence interval
    alpha <- 1 - ci_level
    ci_lower <- quantile(node_values, alpha/2, na.rm = TRUE)
    ci_upper <- quantile(node_values, 1 - alpha/2, na.rm = TRUE)
    
    # Store node results
    node_bootstraps[[node_str]] <- list(
      node_id = node,
      original_value = original_recon[node_str],
      bootstrap_values = node_values,
      mean = mean(node_values, na.rm = TRUE),
      median = median(node_values, na.rm = TRUE),
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_level = ci_level
    )
  }
  
  # Prepare final results
  results <- list(
    tree = tree,
    original_reconstruction = original_recon,
    reconstruction_method = reconstruction_method,
    n_bootstraps = n_bootstraps,
    successful_bootstraps = length(bootstrap_results),
    node_results = node_bootstraps,
    summary = list(
      n_nodes = length(nodes),
      ci_level = ci_level,
      mean_ci_width = mean(sapply(node_bootstraps, function(x) x$ci_upper - x$ci_lower)),
      reconstruction_method = reconstruction_method,
      bootstrap_type = "non-parametric"
    )
  )
  
  class(results) <- c("chr_bootstrap", class(results))
  
  return(results)
}

#' Compare bootstrap results from different reconstruction methods
#' 
#' Analyzes differences in confidence intervals and ancestral reconstructions
#' between different reconstruction methods
#' 
#' @param bootstrap_list List of bootstrap results from different methods
#' @param method_names Optional vector of method names to label results
#' @param nodes Vector of nodes to compare (NULL for all common nodes)
#' @param output_dir Optional directory to save comparison plots
#' @return List with method comparison results
#' @export
compare_bootstrap_methods <- function(bootstrap_list,
                                    method_names = NULL,
                                    nodes = NULL,
                                    output_dir = NULL) {
  # Validate inputs
  if(!is.list(bootstrap_list) || length(bootstrap_list) < 2) {
    stop("bootstrap_list must be a list of at least two bootstrap results")
  }
  
  # Check that all elements are bootstrap results
  if(!all(sapply(bootstrap_list, function(x) inherits(x, "chr_bootstrap")))) {
    stop("All elements in bootstrap_list must be results from bootstrap analysis")
  }
  
  # Set method names if not provided
  if(is.null(method_names)) {
    method_names <- sapply(bootstrap_list, function(x) x$reconstruction_method)
    # If names are still not unique, add numbers
    if(length(unique(method_names)) < length(method_names)) {
      method_names <- paste0(method_names, "_", 1:length(method_names))
    }
  } else {
    if(length(method_names) != length(bootstrap_list)) {
      stop("method_names must have the same length as bootstrap_list")
    }
  }
  
  # Find nodes common to all bootstrap results
  common_nodes <- Reduce(intersect, lapply(bootstrap_list, function(x) names(x$node_results)))
  
  # Filter to specified nodes if provided
  if(!is.null(nodes)) {
    nodes <- as.character(nodes)
    common_nodes <- intersect(common_nodes, nodes)
  }
  
  if(length(common_nodes) == 0) {
    stop("No common nodes found across all bootstrap results")
  }
  
  # Initialize results
  results <- list(
    method_names = method_names,
    n_methods = length(bootstrap_list),
    common_nodes = common_nodes,
    node_comparisons = list(),
    summary = list()
  )
  
  # Compare methods for each node
  for(node in common_nodes) {
    # Extract results for this node across all methods
    node_results <- lapply(bootstrap_list, function(x) x$node_results[[node]])
    
    # Create comparison data frame
    comparison <- data.frame(
      Method = method_names,
      Original_Value = sapply(node_results, function(x) x$original_value),
      Mean = sapply(node_results, function(x) x$mean),
      Median = sapply(node_results, function(x) x$median),
      CI_Lower = sapply(node_results, function(x) x$ci_lower),
      CI_Upper = sapply(node_results, function(x) x$ci_upper),
      CI_Width = sapply(node_results, function(x) x$ci_upper - x$ci_lower),
      stringsAsFactors = FALSE
    )
    
    # Calculate agreement metrics
    # 1. Calculate overlap between CIs
    ci_overlaps <- matrix(NA, nrow = length(method_names), ncol = length(method_names))
    rownames(ci_overlaps) <- colnames(ci_overlaps) <- method_names
    
    for(i in 1:(length(method_names) - 1)) {
      for(j in (i + 1):length(method_names)) {
        # Calculate overlap between CIs
        ci1 <- c(comparison$CI_Lower[i], comparison$CI_Upper[i])
        ci2 <- c(comparison$CI_Lower[j], comparison$CI_Upper[j])
        
        # Check if CIs overlap
        if(ci1[2] < ci2[1] || ci2[2] < ci1[1]) {
          # No overlap
          ci_overlaps[i, j] <- ci_overlaps[j, i] <- 0
        } else {
          # Calculate overlap percentage
          overlap <- min(ci1[2], ci2[2]) - max(ci1[1], ci2[1])
          total_range <- max(ci1[2], ci2[2]) - min(ci1[1], ci2[1])
          ci_overlaps[i, j] <- ci_overlaps[j, i] <- overlap / total_range
        }
      }
    }
    
    # 2. Calculate absolute differences in point estimates
    point_diffs <- matrix(NA, nrow = length(method_names), ncol = length(method_names))
    rownames(point_diffs) <- colnames(point_diffs) <- method_names
    
    for(i in 1:(length(method_names) - 1)) {
      for(j in (i + 1):length(method_names)) {
        # Calculate absolute difference
        point_diffs[i, j] <- point_diffs[j, i] <- abs(comparison$Original_Value[i] - comparison$Original_Value[j])
      }
    }
    
    # Store node comparison
    results$node_comparisons[[node]] <- list(
      comparison = comparison,
      ci_overlaps = ci_overlaps,
      point_diffs = point_diffs,
      agreement_summary = list(
        mean_ci_overlap = mean(ci_overlaps[upper.tri(ci_overlaps)]),
        mean_point_diff = mean(point_diffs[upper.tri(point_diffs)]),
        all_cis_overlap = all(ci_overlaps[upper.tri(ci_overlaps)] > 0)
      )
    )
    
    # Create comparison plot
    if(requireNamespace("ggplot2", quietly = TRUE)) {
      # Create plot data for CI visualization
      plot_data <- data.frame(
        Method = comparison$Method,
        Original = comparison$Original_Value,
        Lower = comparison$CI_Lower,
        Upper = comparison$CI_Upper,
        stringsAsFactors = FALSE
      )
      
      # Order methods by their point estimates
      plot_data$Method <- factor(plot_data$Method, 
                               levels = plot_data$Method[order(plot_data$Original)])
      
      # Create CI comparison plot
      p <- ggplot(plot_data, aes(x = Method, y = Original, color = Method)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, size = 1) +
        labs(
          title = paste("Node", node, "- Method Comparison"),
          subtitle = "Point estimates with confidence intervals",
          y = "Chromosome Number",
          x = "Reconstruction Method"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
      
      results$node_comparisons[[node]]$plot <- p
      
      # Save plot if output directory provided
      if(!is.null(output_dir)) {
        if(!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        
        # Save plot
        filename <- file.path(output_dir, paste0("node_", node, "_comparison.pdf"))
        ggsave(filename, p, width = 8, height = 6)
      }
    }
  }
  
  # Create overall summary
  # Calculate overall metrics
  ci_overlap_summary <- sapply(results$node_comparisons, function(x) x$agreement_summary$mean_ci_overlap)
  point_diff_summary <- sapply(results$node_comparisons, function(x) x$agreement_summary$mean_point_diff)
  ci_width_summary <- lapply(results$node_comparisons, function(x) {
    x$comparison$CI_Width
  })
  
  # Convert to data frame for easier analysis
  ci_width_df <- do.call(rbind, lapply(1:length(ci_width_summary), function(i) {
    data.frame(
      Node = common_nodes[i],
      Method = method_names,
      CI_Width = ci_width_summary[[i]],
      stringsAsFactors = FALSE
    )
  }))
  
  # Calculate method-specific metrics
  method_summary <- aggregate(CI_Width ~ Method, data = ci_width_df, FUN = mean)
  method_summary$CI_Width_Rank <- rank(method_summary$CI_Width)
  
  results$summary <- list(
    mean_ci_overlap = mean(ci_overlap_summary),
    mean_point_diff = mean(point_diff_summary),
    all_methods_agree = sum(sapply(results$node_comparisons, function(x) x$agreement_summary$all_cis_overlap)) / length(common_nodes),
    method_summary = method_summary,
    best_confidence = method_summary$Method[which.min(method_summary$CI_Width)],
    node_level_summary = data.frame(
      Node = common_nodes,
      Mean_CI_Overlap = ci_overlap_summary,
      Mean_Point_Diff = point_diff_summary,
      stringsAsFactors = FALSE
    )
  )
  
  # Create overall comparison plot
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    overall_plot <- ggplot(ci_width_df, aes(x = Method, y = CI_Width, fill = Method)) +
      geom_boxplot() +
      labs(
        title = "Comparison of Confidence Interval Widths Across Methods",
        x = "Reconstruction Method",
        y = "Confidence Interval Width"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
    
    results$summary$overall_plot <- overall_plot
    
    # Save overall plot if output directory provided
    if(!is.null(output_dir)) {
      filename <- file.path(output_dir, "method_comparison_overview.pdf")
      ggsave(filename, overall_plot, width = 10, height = 8)
    }
  }
  
  # Set class for results object
  class(results) <- c("chr_method_comparison", class(results))
  
  return(results)
}

#===============================================================================
# Simulation Accuracy Evaluation Functions
#===============================================================================

#' Evaluate accuracy of ancestral reconstruction methods using simulations
#' 
#' Simulates chromosome evolution with known ancestral states, then tests
#' reconstruction methods against the true values
#' 
#' @param tree Phylogenetic tree
#' @param simulation_model Model for simulation: "BM", "discrete", "fusion_fission"
#' @param simulation_params Parameters for simulation
#' @param root_value Root chromosome number for simulation
#' @param reconstruction_methods Vector of methods to test
#' @param n_simulations Number of simulation replicates
#' @param discretize Whether to round simulated values to integers
#' @param n_cores Number of cores for parallel processing
#' @param seed Random seed for reproducibility
#' @return List with method performance metrics
#' @export
evaluate_reconstruction_accuracy <- function(tree,
                                          simulation_model = "discrete",
                                          simulation_params = list(fusion_rate = 0.05, fission_rate = 0.05),
                                          root_value = 10,
                                          reconstruction_methods = c("parsimony", "ML", "Bayesian"),
                                          n_simulations = 50,
                                          discretize = TRUE,
                                          n_cores = NULL,
                                          seed = NULL) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Set random seed if provided
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # Set up parallel processing if requested
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Set up storage for results
  n_methods <- length(reconstruction_methods)
  n_internal_nodes <- tree$Nnode
  n_tips <- length(tree$tip.label)
  
  # Function to run a single simulation and test reconstruction methods
  run_simulation <- function(sim_idx) {
    # Simulate chromosome evolution
    sim_values <- simulate_chromosome_evolution(
      tree = tree,
      model = simulation_model,
      params = simulation_params,
      root_value = root_value,
      discretize = discretize,
      seed = if(!is.null(seed)) seed + sim_idx else NULL
    )
    
    # Extract true ancestral states
    true_states <- sim_values[(n_tips + 1):(n_tips + n_internal_nodes)]
    
    # Extract tip values to use for reconstruction
    tip_values <- sim_values[1:n_tips]
    names(tip_values) <- tree$tip.label
    
    # Initialize storage for method results
    method_results <- list()
    
    # Test each reconstruction method
    for(method in reconstruction_methods) {
      # Reconstruct ancestral states
      recon_states <- try({
        # Use custom function to get only internal node states
        reconstructed <- reconstruct_ancestral_states(tree, tip_values, method)
        reconstructed[(n_tips + 1):(n_tips + n_internal_nodes)]
      }, silent = TRUE)
      
      # If reconstruction failed, store NAs
      if(inherits(recon_states, "try-error")) {
        recon_states <- rep(NA, n_internal_nodes)
        method_results[[method]] <- list(
          success = FALSE,
          error = attr(recon_states, "condition")$message
        )
        next
      }
      
      # Calculate accuracy metrics
      errors <- true_states - recon_states
      abs_errors <- abs(errors)
      
      # Compile results
      method_results[[method]] <- list(
        success = TRUE,
        reconstructed = recon_states,
        true_values = true_states,
        errors = errors,
        abs_errors = abs_errors,
        mae = mean(abs_errors, na.rm = TRUE),  # Mean Absolute Error
        rmse = sqrt(mean(errors^2, na.rm = TRUE)),  # Root Mean Squared Error
        correlation = cor(true_states, recon_states, use = "pairwise.complete.obs")
      )
      
      # Calculate additional accuracy metrics
      # Proportion of correctly reconstructed values (exact match)
      method_results[[method]]$exact_match <- sum(abs_errors == 0, na.rm = TRUE) / sum(!is.na(abs_errors))
      
      # Proportion within 1 chromosome
      method_results[[method]]$within_1 <- sum(abs_errors <= 1, na.rm = TRUE) / sum(!is.na(abs_errors))
      
      # Proportion within 2 chromosomes
      method_results[[method]]$within_2 <- sum(abs_errors <= 2, na.rm = TRUE) / sum(!is.na(abs_errors))
      
      # Direction accuracy (increase vs. decrease from parent)
      # (This is more complex and would require analyzing pairs of nodes)
    }
    
    # Return simulation results
    return(list(
      sim_idx = sim_idx,
      true_values = sim_values,
      method_results = method_results
    ))
  }
  
  # Run all simulations (in parallel if requested)
  if(n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    # Export necessary functions and objects to cluster
    parallel::clusterExport(cl, c("simulate_chromosome_evolution", 
                               "reconstruct_ancestral_states",
                               "simulation_model", 
                               "simulation_params", 
                               "root_value",
                               "tree",
                               "n_tips",
                               "n_internal_nodes",
                               "discretize"))
    
    # Run simulations in parallel
    sim_results <- parallel::parLapply(cl, 1:n_simulations, run_simulation)
  } else {
    # Run simulations sequentially
    sim_results <- lapply(1:n_simulations, run_simulation)
  }
  
  # Compile results across simulations
  method_summary <- list()
  
  for(method in reconstruction_methods) {
    # Extract results for this method from all simulations
    method_data <- lapply(sim_results, function(x) x$method_results[[method]])
    
    # Count successful reconstructions
    successful <- sum(sapply(method_data, function(x) !is.null(x) && x$success))
    
    # Skip if no successful reconstructions
    if(successful == 0) {
      method_summary[[method]] <- list(
        success_rate = 0,
        error = "All reconstructions failed"
      )
      next
    }
    
    # Extract metrics from successful reconstructions
    successful_data <- method_data[sapply(method_data, function(x) !is.null(x) && x$success)]
    
    # Calculate mean metrics across simulations
    mae_values <- sapply(successful_data, function(x) x$mae)
    rmse_values <- sapply(successful_data, function(x) x$rmse)
    cor_values <- sapply(successful_data, function(x) x$correlation)
    exact_match_values <- sapply(successful_data, function(x) x$exact_match)
    within_1_values <- sapply(successful_data, function(x) x$within_1)
    within_2_values <- sapply(successful_data, function(x) x$within_2)
    
    # Compile method summary
    method_summary[[method]] <- list(
      success_rate = successful / n_simulations,
      mean_mae = mean(mae_values, na.rm = TRUE),
      se_mae = sd(mae_values, na.rm = TRUE) / sqrt(length(mae_values)),
      mean_rmse = mean(rmse_values, na.rm = TRUE),
      se_rmse = sd(rmse_values, na.rm = TRUE) / sqrt(length(rmse_values)),
      mean_correlation = mean(cor_values, na.rm = TRUE),
      se_correlation = sd(cor_values, na.rm = TRUE) / sqrt(length(cor_values)),
      mean_exact_match = mean(exact_match_values, na.rm = TRUE),
      se_exact_match = sd(exact_match_values, na.rm = TRUE) / sqrt(length(exact_match_values)),
      mean_within_1 = mean(within_1_values, na.rm = TRUE),
      se_within_1 = sd(within_1_values, na.rm = TRUE) / sqrt(length(within_1_values)),
      mean_within_2 = mean(within_2_values, na.rm = TRUE),
      se_within_2 = sd(within_2_values, na.rm = TRUE) / sqrt(length(within_2_values)),
      all_mae = mae_values,
      all_rmse = rmse_values,
      all_correlation = cor_values
    )
  }
  
  # Create method comparison table
  if(length(method_summary) > 0) {
    comparison_table <- data.frame(
      Method = names(method_summary),
      Success_Rate = sapply(method_summary, function(x) x$success_rate),
      MAE = sapply(method_summary, function(x) x$mean_mae),
      RMSE = sapply(method_summary, function(x) x$mean_rmse),
      Correlation = sapply(method_summary, function(x) x$mean_correlation),
      Exact_Match = sapply(method_summary, function(x) x$mean_exact_match),
      Within_1 = sapply(method_summary, function(x) x$mean_within_1),
      Within_2 = sapply(method_summary, function(x) x$mean_within_2),
      stringsAsFactors = FALSE
    )
    
    # Rank methods by different metrics
    comparison_table$MAE_Rank <- rank(comparison_table$MAE)
    comparison_table$RMSE_Rank <- rank(comparison_table$RMSE)
    comparison_table$Correlation_Rank <- rank(-comparison_table$Correlation)  # Higher correlation is better
    comparison_table$Exact_Match_Rank <- rank(-comparison_table$Exact_Match)  # Higher proportion is better
    comparison_table$Overall_Rank <- rowMeans(comparison_table[, c("MAE_Rank", "RMSE_Rank", "Correlation_Rank", "Exact_Match_Rank")])
  } else {
    comparison_table <- data.frame()
  }
  
  # Prepare overall results
  results <- list(
    n_simulations = n_simulations,
    n_successful = sapply(method_summary, function(x) x$success_rate * n_simulations),
    simulation_model = simulation_model,
    simulation_params = simulation_params,
    reconstruction_methods = reconstruction_methods,
    method_summary = method_summary,
    comparison_table = comparison_table,
    best_method = if(nrow(comparison_table) > 0) comparison_table$Method[which.min(comparison_table$Overall_Rank)] else NA,
    sim_results = sim_results
  )
  
  # Create visualization if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE) && length(method_summary) > 0) {
    # Create data for error boxplots
    error_data <- do.call(rbind, lapply(names(method_summary), function(method) {
      data.frame(
        Method = method,
        MAE = method_summary[[method]]$all_mae,
        RMSE = method_summary[[method]]$all_rmse,
        Correlation = method_summary[[method]]$all_correlation,
        stringsAsFactors = FALSE
      )
    }))
    
    # Create boxplot of MAE by method
    mae_plot <- ggplot(error_data, aes(x = Method, y = MAE, fill = Method)) +
      geom_boxplot() +
      labs(
        title = "Mean Absolute Error by Reconstruction Method",
        x = "Method",
        y = "MAE"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Create boxplot of correlation by method
    cor_plot <- ggplot(error_data, aes(x = Method, y = Correlation, fill = Method)) +
      geom_boxplot() +
      labs(
        title = "True vs. Reconstructed Correlation by Method",
        x = "Method",
        y = "Correlation"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Store plots in results
    results$plots <- list(
      mae_plot = mae_plot,
      cor_plot = cor_plot
    )
    
    # Combine plots if patchwork is available
    if(requireNamespace("patchwork", quietly = TRUE)) {
      results$plots$combined <- mae_plot + cor_plot +
        patchwork::plot_annotation(
          title = "Reconstruction Method Performance Comparison",
          subtitle = paste("Based on", n_simulations, "simulations using", simulation_model, "model")
        )
    }
  }
  
  # Set class for results object
  class(results) <- c("chr_accuracy_evaluation", class(results))
  
  return(results)
}
