#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Simulation Module
# Author: Bioinformatics Team
# Date: 2025-04-01
# Description: Provides functionality for simulating chromosome evolution under
#              different models for testing and benchmarking reconstruction methods
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(MASS)
  library(parallel)
})

#===============================================================================
# Simulation Core Functions
#===============================================================================

#' Simulate chromosome evolution on a phylogenetic tree
#' 
#' Simulates chromosome number evolution under various models
#' 
#' @param tree Phylogenetic tree
#' @param model Evolution model: "BM" (Brownian Motion), "OU" (Ornstein-Uhlenbeck),
#'              "ACDC" (Accelerating/Decelerating), "jumps" (punctuated evolution),
#'              "hybrid" (mixed continuous/punctuated evolution), or "custom"
#' @param params List of model parameters
#' @param root_value Starting chromosome number at root (ancestral state)
#' @param constraints List of evolutionary constraints (min_value, max_value, etc.)
#' @param discrete Whether to discretize (round) chromosome numbers
#' @return Simulated chromosome evolution data
#' @export
simulate_chromosome_evolution <- function(tree, model = "BM", 
                                        params = NULL, 
                                        root_value = 12,
                                        constraints = list(min_value = 1, 
                                                           max_value = 100),
                                        discrete = TRUE) {
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Check if tree has branch lengths
  if(is.null(tree$edge.length)) {
    stop("Tree must have branch lengths for simulation")
  }
  
  # Verify model
  supported_models <- c("BM", "OU", "ACDC", "jumps", "hybrid", "custom")
  if(!model %in% supported_models) {
    stop(paste("Unsupported model. Choose from:", paste(supported_models, collapse = ", ")))
  }
  
  # Set default parameters if not provided
  if(is.null(params)) {
    params <- get_default_parameters(model)
  } else {
    # Merge with defaults for any missing parameters
    default_params <- get_default_parameters(model)
    missing_params <- setdiff(names(default_params), names(params))
    for(param in missing_params) {
      params[[param]] <- default_params[[param]]
    }
  }
  
  message(sprintf("Simulating chromosome evolution under %s model...", model))
  
  # Apply simulation model
  if(model == "BM") {
    sim_result <- simulate_bm(tree, root_value, params, constraints)
  } else if(model == "OU") {
    sim_result <- simulate_ou(tree, root_value, params, constraints)
  } else if(model == "ACDC") {
    sim_result <- simulate_acdc(tree, root_value, params, constraints)
  } else if(model == "jumps") {
    sim_result <- simulate_jumps(tree, root_value, params, constraints)
  } else if(model == "hybrid") {
    sim_result <- simulate_hybrid(tree, root_value, params, constraints)
  } else if(model == "custom") {
    if(is.null(params$sim_func)) {
      stop("Custom model requires a simulation function in params$sim_func")
    }
    sim_result <- params$sim_func(tree, root_value, params, constraints)
  }
  
  # Discretize if requested
  if(discrete) {
    sim_result$tip_states <- round(sim_result$tip_states)
    sim_result$node_states <- round(sim_result$node_states)
    sim_result$discrete <- TRUE
  }
  
  # Return simulation result with metadata
  result <- list(
    tree = tree,
    tip_states = sim_result$tip_states,
    node_states = sim_result$node_states,
    model = model,
    params = params,
    root_value = root_value,
    constraints = constraints,
    events = sim_result$events,
    discrete = discrete
  )
  
  # Add some statistics
  result$statistics <- calculate_simulation_stats(result)
  
  return(result)
}

#' Get default parameters for simulation models
#' 
#' @param model Evolution model
#' @return Default parameter list
#' @keywords internal
get_default_parameters <- function(model) {
  if(model == "BM") {
    # Brownian Motion
    return(list(
      sigma = 1.0,   # Rate parameter (standard deviation of steps)
      trend = 0.0    # Optional directional trend
    ))
  } else if(model == "OU") {
    # Ornstein-Uhlenbeck
    return(list(
      sigma = 1.0,   # Rate parameter
      alpha = 0.1,   # Strength of selection
      theta = 12     # Optimum value
    ))
  } else if(model == "ACDC") {
    # Accelerating/Decelerating
    return(list(
      sigma = 1.0,     # Base rate parameter
      beta = 0.1,      # Rate change parameter (positive = accelerating, negative = decelerating)
      min_rate = 0.1   # Minimum rate to prevent numerical issues
    ))
  } else if(model == "jumps") {
    # Punctuated evolution
    return(list(
      background_rate = 0.1,  # Background evolution rate
      jump_rate = 0.1,        # Rate of jump events per unit branch length
      jump_size_mean = 0,     # Mean jump size (0 = symmetric)
      jump_size_sd = 2,       # Standard deviation of jump sizes
      jump_distribution = "normal"  # Distribution of jumps: "normal", "exponential", or "custom"
    ))
  } else if(model == "hybrid") {
    # Hybrid of continuous and punctuated evolution
    return(list(
      sigma = 0.5,            # Rate for continuous BM component
      jump_rate = 0.05,       # Rate of jump events per unit branch length
      jump_size_mean = 0,     # Mean jump size
      jump_size_sd = 3,       # Standard deviation of jump sizes
      fusion_prob = 0.5,      # Probability of fusion vs. fission when jumping
      wgd_prob = 0.1,         # Probability of whole genome duplication jumps
      trend = 0.0             # Directional trend in continuous component
    ))
  } else if(model == "custom") {
    # Custom model requires user to provide simulation function
    return(list(
      sim_func = NULL  # User must provide this
    ))
  } else {
    stop(paste("Unsupported model:", model))
  }
}

#' Simulate under Brownian Motion model
#' 
#' @param tree Phylogenetic tree
#' @param root_value Root state
#' @param params Model parameters
#' @param constraints Evolutionary constraints
#' @return Simulation results
#' @keywords internal
simulate_bm <- function(tree, root_value, params, constraints) {
  # Extract parameters
  sigma <- params$sigma
  trend <- params$trend
  
  # Setup simulation containers
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  
  # Initialize all states
  all_states <- numeric(n_total)
  all_states[n_tips + 1] <- root_value  # Set root state
  
  # Initialize events tracking
  events <- data.frame(
    edge = integer(0),
    parent = integer(0),
    child = integer(0),
    time = numeric(0),
    type = character(0),
    magnitude = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Simulate along each edge
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent state
    parent_state <- all_states[parent]
    
    # Calculate expected change due to trend
    expected_change <- trend * branch_length
    
    # Simulate state change using Brownian motion
    random_change <- rnorm(1, mean = 0, sd = sigma * sqrt(branch_length))
    total_change <- expected_change + random_change
    
    # Calculate child state
    child_state <- parent_state + total_change
    
    # Apply constraints
    if(!is.null(constraints$min_value)) {
      child_state <- max(constraints$min_value, child_state)
    }
    if(!is.null(constraints$max_value)) {
      child_state <- min(constraints$max_value, child_state)
    }
    
    # Store child state
    all_states[child] <- child_state
  }
  
  # Extract tip and node states
  tip_states <- all_states[1:n_tips]
  names(tip_states) <- tree$tip.label
  node_states <- all_states[(n_tips+1):n_total]
  
  return(list(
    tip_states = tip_states,
    node_states = node_states,
    events = events
  ))
}

#' Simulate under Ornstein-Uhlenbeck model
#' 
#' @param tree Phylogenetic tree
#' @param root_value Root state
#' @param params Model parameters
#' @param constraints Evolutionary constraints
#' @return Simulation results
#' @keywords internal
simulate_ou <- function(tree, root_value, params, constraints) {
  # Extract parameters
  sigma <- params$sigma
  alpha <- params$alpha
  theta <- params$theta
  
  # Setup simulation containers
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  
  # Initialize all states
  all_states <- numeric(n_total)
  all_states[n_tips + 1] <- root_value  # Set root state
  
  # Initialize events tracking
  events <- data.frame(
    edge = integer(0),
    parent = integer(0),
    child = integer(0),
    time = numeric(0),
    type = character(0),
    magnitude = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Simulate along each edge
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent state
    parent_state <- all_states[parent]
    
    # Calculate expected value under OU model
    expected_value <- theta + (parent_state - theta) * exp(-alpha * branch_length)
    
    # Calculate variance under OU model
    variance <- (sigma^2 / (2 * alpha)) * (1 - exp(-2 * alpha * branch_length))
    
    # Simulate state
    child_state <- rnorm(1, mean = expected_value, sd = sqrt(variance))
    
    # Apply constraints
    if(!is.null(constraints$min_value)) {
      child_state <- max(constraints$min_value, child_state)
    }
    if(!is.null(constraints$max_value)) {
      child_state <- min(constraints$max_value, child_state)
    }
    
    # Store child state
    all_states[child] <- child_state
  }
  
  # Extract tip and node states
  tip_states <- all_states[1:n_tips]
  names(tip_states) <- tree$tip.label
  node_states <- all_states[(n_tips+1):n_total]
  
  return(list(
    tip_states = tip_states,
    node_states = node_states,
    events = events
  ))
}

#' Simulate under Accelerating/Decelerating model
#' 
#' @param tree Phylogenetic tree
#' @param root_value Root state
#' @param params Model parameters
#' @param constraints Evolutionary constraints
#' @return Simulation results
#' @keywords internal
simulate_acdc <- function(tree, root_value, params, constraints) {
  # Extract parameters
  sigma <- params$sigma
  beta <- params$beta  # Rate change parameter
  min_rate <- params$min_rate
  
  # Setup simulation containers
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  
  # Initialize all states
  all_states <- numeric(n_total)
  all_states[n_tips + 1] <- root_value  # Set root state
  
  # Calculate node ages
  node_ages <- calculate_node_ages(tree)
  
  # Initialize events tracking
  events <- data.frame(
    edge = integer(0),
    parent = integer(0),
    child = integer(0),
    time = numeric(0),
    type = character(0),
    magnitude = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Simulate along each edge
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent state
    parent_state <- all_states[parent]
    
    # Calculate age at start of edge
    start_age <- node_ages[parent]
    
    # Determine rate modifiers based on age
    rate_modifier <- exp(beta * start_age)
    
    # Ensure minimum rate
    rate_modifier <- max(rate_modifier, min_rate / sigma)
    
    # Apply rate to sigma
    effective_sigma <- sigma * rate_modifier
    
    # Simulate state change using modified Brownian motion
    random_change <- rnorm(1, mean = 0, sd = effective_sigma * sqrt(branch_length))
    
    # Calculate child state
    child_state <- parent_state + random_change
    
    # Apply constraints
    if(!is.null(constraints$min_value)) {
      child_state <- max(constraints$min_value, child_state)
    }
    if(!is.null(constraints$max_value)) {
      child_state <- min(constraints$max_value, child_state)
    }
    
    # Store child state
    all_states[child] <- child_state
  }
  
  # Extract tip and node states
  tip_states <- all_states[1:n_tips]
  names(tip_states) <- tree$tip.label
  node_states <- all_states[(n_tips+1):n_total]
  
  return(list(
    tip_states = tip_states,
    node_states = node_states,
    events = events
  ))
}

#' Simulate under punctuated evolution model with jumps
#' 
#' @param tree Phylogenetic tree
#' @param root_value Root state
#' @param params Model parameters
#' @param constraints Evolutionary constraints
#' @return Simulation results
#' @keywords internal
simulate_jumps <- function(tree, root_value, params, constraints) {
  # Extract parameters
  background_rate <- params$background_rate
  jump_rate <- params$jump_rate
  jump_size_mean <- params$jump_size_mean
  jump_size_sd <- params$jump_size_sd
  jump_distribution <- params$jump_distribution
  
  # Setup simulation containers
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  
  # Initialize all states
  all_states <- numeric(n_total)
  all_states[n_tips + 1] <- root_value  # Set root state
  
  # Initialize events tracking
  events <- data.frame(
    edge = integer(0),
    parent = integer(0),
    child = integer(0),
    time = numeric(0),
    type = character(0),
    magnitude = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Simulate along each edge
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent state
    parent_state <- all_states[parent]
    child_state <- parent_state
    
    # First, add background evolution (small continuous changes)
    if(background_rate > 0) {
      bg_change <- rnorm(1, mean = 0, sd = background_rate * sqrt(branch_length))
      child_state <- child_state + bg_change
    }
    
    # Determine number of jumps on this branch
    n_jumps <- rpois(1, lambda = jump_rate * branch_length)
    
    # Simulate each jump
    if(n_jumps > 0) {
      for(j in 1:n_jumps) {
        # Generate jump time (uniformly distributed along branch)
        jump_time <- runif(1, 0, branch_length)
        
        # Generate jump size based on specified distribution
        if(jump_distribution == "normal") {
          jump_size <- rnorm(1, mean = jump_size_mean, sd = jump_size_sd)
        } else if(jump_distribution == "exponential") {
          # For exponential, we'll use the direction from jump_size_mean
          direction <- sign(jump_size_mean)
          if(direction == 0) direction <- sample(c(-1, 1), 1)
          jump_size <- direction * rexp(1, rate = 1/abs(jump_size_sd))
        } else {
          # Default to normal
          jump_size <- rnorm(1, mean = jump_size_mean, sd = jump_size_sd)
        }
        
        # Apply jump
        child_state <- child_state + jump_size
        
        # Determine event type
        event_type <- if(jump_size > 0) "fission" else "fusion"
        
        # Record event
        events <- rbind(events, data.frame(
          edge = i,
          parent = parent,
          child = child,
          time = jump_time,
          type = event_type,
          magnitude = abs(jump_size),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Apply constraints
    if(!is.null(constraints$min_value)) {
      child_state <- max(constraints$min_value, child_state)
    }
    if(!is.null(constraints$max_value)) {
      child_state <- min(constraints$max_value, child_state)
    }
    
    # Store child state
    all_states[child] <- child_state
  }
  
  # Extract tip and node states
  tip_states <- all_states[1:n_tips]
  names(tip_states) <- tree$tip.label
  node_states <- all_states[(n_tips+1):n_total]
  
  return(list(
    tip_states = tip_states,
    node_states = node_states,
    events = events
  ))
}

#' Simulate under hybrid model (continuous with discrete jumps)
#' 
#' @param tree Phylogenetic tree
#' @param root_value Root state
#' @param params Model parameters
#' @param constraints Evolutionary constraints
#' @return Simulation results
#' @keywords internal
simulate_hybrid <- function(tree, root_value, params, constraints) {
  # Extract parameters
  sigma <- params$sigma
  jump_rate <- params$jump_rate
  jump_size_mean <- params$jump_size_mean
  jump_size_sd <- params$jump_size_sd
  fusion_prob <- params$fusion_prob
  wgd_prob <- params$wgd_prob
  trend <- params$trend
  
  # Setup simulation containers
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  
  # Initialize all states
  all_states <- numeric(n_total)
  all_states[n_tips + 1] <- root_value  # Set root state
  
  # Initialize events tracking
  events <- data.frame(
    edge = integer(0),
    parent = integer(0),
    child = integer(0),
    time = numeric(0),
    type = character(0),
    magnitude = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Simulate along each edge
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent state
    parent_state <- all_states[parent]
    child_state <- parent_state
    
    # Add continuous background evolution with trend
    bg_change <- rnorm(1, mean = trend * branch_length, sd = sigma * sqrt(branch_length))
    child_state <- child_state + bg_change
    
    # Determine number of jumps on this branch
    n_jumps <- rpois(1, lambda = jump_rate * branch_length)
    
    # Simulate each jump
    if(n_jumps > 0) {
      for(j in 1:n_jumps) {
        # Generate jump time (uniformly distributed along branch)
        jump_time <- runif(1, 0, branch_length)
        
        # Determine jump type (WGD or fusion/fission)
        if(runif(1) < wgd_prob) {
          # Whole genome duplication - approximately doubles chromosome number
          jump_size <- child_state * (runif(1, 0.9, 1.1))  # Add some noise around exact doubling
          child_state <- child_state + jump_size
          event_type <- "wgd"
        } else {
          # Fusion or fission event
          if(runif(1) < fusion_prob) {
            # Fusion - chromosome number decreases
            jump_size <- -abs(rnorm(1, mean = jump_size_mean, sd = jump_size_sd))
            event_type <- "fusion"
          } else {
            # Fission - chromosome number increases
            jump_size <- abs(rnorm(1, mean = jump_size_mean, sd = jump_size_sd))
            event_type <- "fission"
          }
          child_state <- child_state + jump_size
        }
        
        # Record event
        events <- rbind(events, data.frame(
          edge = i,
          parent = parent,
          child = child,
          time = jump_time,
          type = event_type,
          magnitude = abs(jump_size),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Apply constraints
    if(!is.null(constraints$min_value)) {
      child_state <- max(constraints$min_value, child_state)
    }
    if(!is.null(constraints$max_value)) {
      child_state <- min(constraints$max_value, child_state)
    }
    
    # Store child state
    all_states[child] <- child_state
  }
  
  # Extract tip and node states
  tip_states <- all_states[1:n_tips]
  names(tip_states) <- tree$tip.label
  node_states <- all_states[(n_tips+1):n_total]
  
  return(list(
    tip_states = tip_states,
    node_states = node_states,
    events = events
  ))
}

#' Calculate node ages from a phylogenetic tree
#' 
#' @param tree Phylogenetic tree
#' @return Vector of node ages
#' @keywords internal
calculate_node_ages <- function(tree) {
  # Get max root-to-tip distance (tree height)
  heights <- node.depth.edgelength(tree)
  tree_height <- max(heights)
  
  # Convert from node depths to ages (time from present)
  ages <- tree_height - heights
  
  return(ages)
}

#' Calculate simulation statistics
#' 
#' @param sim_result Simulation result
#' @return List of summary statistics
#' @keywords internal
calculate_simulation_stats <- function(sim_result) {
  # Basic statistics
  tip_states <- sim_result$tip_states
  node_states <- sim_result$node_states
  
  stats <- list(
    n_tips = length(tip_states),
    n_nodes = length(node_states),
    tip_mean = mean(tip_states),
    tip_median = median(tip_states),
    tip_sd = sd(tip_states),
    tip_range = range(tip_states),
    node_mean = mean(node_states),
    node_median = median(node_states),
    node_sd = sd(node_states),
    node_range = range(node_states)
  )
  
  # Event statistics if events are present
  if(!is.null(sim_result$events) && nrow(sim_result$events) > 0) {
    events <- sim_result$events
    
    # Count by event type
    event_counts <- table(events$type)
    
    # Calculate average magnitude by type
    event_magnitude <- tapply(events$magnitude, events$type, mean)
    
    stats$events <- list(
      total_events = nrow(events),
      event_counts = event_counts,
      event_magnitude = event_magnitude
    )
  }
  
  return(stats)
}

#===============================================================================
# Batch Simulation Functions
#===============================================================================

#' Run multiple simulations under different conditions
#' 
#' @param tree Phylogenetic tree
#' @param models Character vector of models to test
#' @param n_replicates Number of replicates per model
#' @param root_values Vector of root values to test (NULL for default)
#' @param param_variations List of parameter variations to test (NULL for default)
#' @param constraints Evolutionary constraints
#' @param discrete Whether to discretize chromosome numbers
#' @param use_parallel Whether to use parallel processing
#' @return List of simulation results
#' @export
simulate_chromosome_batch <- function(tree, 
                                    models = c("BM", "OU", "jumps", "hybrid"),
                                    n_replicates = 10, 
                                    root_values = NULL,
                                    param_variations = NULL,
                                    constraints = list(min_value = 1),
                                    discrete = TRUE,
                                    use_parallel = TRUE) {
  # Default root values if not specified
  if(is.null(root_values)) {
    root_values <- 12
  }
  
  # Create parameter combinations
  if(is.null(param_variations)) {
    param_list <- lapply(models, function(model) {
      get_default_parameters(model)
    })
    names(param_list) <- models
  } else {
    # Use custom parameter variations
    param_list <- param_variations
  }
  
  # Build simulation configurations
  sim_configs <- list()
  config_id <- 1
  
  for(model in models) {
    for(root in root_values) {
      for(params in param_list[[model]]) {
        # Create config
        config <- list(
          id = config_id,
          model = model,
          root_value = root,
          params = params,
          constraints = constraints,
          discrete = discrete
        )
        
        sim_configs[[config_id]] <- config
        config_id <- config_id + 1
      }
    }
  }
  
  # Create replicates of each configuration
  all_configs <- list()
  config_id <- 1
  
  for(base_config in sim_configs) {
    for(rep in 1:n_replicates) {
      config <- base_config
      config$replicate <- rep
      all_configs[[config_id]] <- config
      config_id <- config_id + 1
    }
  }
  
  message(sprintf("Running %d total simulations across %d configurations...", 
                 length(all_configs), length(sim_configs)))
  
  # Run simulations
  if(use_parallel && requireNamespace("parallel", quietly = TRUE)) {
    # Parallel processing
    n_cores <- min(parallel::detectCores() - 1, length(all_configs))
    if(n_cores < 1) n_cores <- 1
    
    message(sprintf("Using %d cores for parallel simulation", n_cores))
    
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, c("simulate_chromosome_evolution", "get_default_parameters",
                               "simulate_bm", "simulate_ou", "simulate_acdc", 
                               "simulate_jumps", "simulate_hybrid",
                               "calculate_node_ages", "calculate_simulation_stats"), 
                          envir = environment())
    
    # Load required packages on each node
    parallel::clusterEvalQ(cl, {
      library(ape)
      library(phytools)
    })
    
    # Run simulations in parallel
    sim_results <- parallel::parLapply(cl, all_configs, function(config) {
      # Run simulation
      sim <- simulate_chromosome_evolution(
        tree = tree,
        model = config$model,
        params = config$params,
        root_value = config$root_value,
        constraints = config$constraints,
        discrete = config$discrete
      )
      
      # Add config info
      sim$config_id <- config$id
      sim$replicate <- config$replicate
      
      return(sim)
    })
    
    # Close cluster
    parallel::stopCluster(cl)
    
  } else {
    # Sequential processing
    sim_results <- lapply(all_configs, function(config) {
      # Run simulation
      sim <- simulate_chromosome_evolution(
        tree = tree,
        model = config$model,
        params = config$params,
        root_value = config$root_value,
        constraints = config$constraints,
        discrete = config$discrete
      )
      
      # Add config info
      sim$config_id <- config$id
      sim$replicate <- config$replicate
      
      return(sim)
    })
  }
  
  # Organize results
  batch_results <- list(
    simulations = sim_results,
    configs = sim_configs,
    models = models,
    n_replicates = n_replicates,
    tree = tree
  )
  
  # Add summary statistics
  batch_results$summary <- summarize_batch_results(batch_results)
  
  message(sprintf("Completed %d simulations", length(sim_results)))
  
  return(batch_results)
}

#' Summarize batch simulation results
#' 
#' @param batch_results Batch simulation results
#' @return Summary statistics
#' @keywords internal
summarize_batch_results <- function(batch_results) {
  # Extract data
  sims <- batch_results$simulations
  configs <- batch_results$configs
  
  # Group simulations by configuration
  grouped_sims <- list()
  
  for(sim in sims) {
    config_id <- sim$config_id
    if(is.null(grouped_sims[[as.character(config_id)]])) {
      grouped_sims[[as.character(config_id)]] <- list(sim)
    } else {
      grouped_sims[[as.character(config_id)]] <- c(grouped_sims[[as.character(config_id)]], list(sim))
    }
  }
  
  # Calculate summary statistics for each configuration
  config_summaries <- list()
  
  for(config_id in names(grouped_sims)) {
    # Get simulations for this config
    config_sims <- grouped_sims[[config_id]]
    
    # Extract tip states from all replicates
    tip_means <- sapply(config_sims, function(sim) mean(sim$tip_states))
    tip_sds <- sapply(config_sims, function(sim) sd(sim$tip_states))
    tip_ranges <- sapply(config_sims, function(sim) diff(range(sim$tip_states)))
    
    # Calculate summary statistics
    config_summary <- list(
      config_id = as.numeric(config_id),
      model = config_sims[[1]]$model,
      n_replicates = length(config_sims),
      tip_state_mean = mean(tip_means),
      tip_state_sd = mean(tip_sds),
      tip_state_range = mean(tip_ranges),
      variability_between_replicates = sd(tip_means)
    )
    
    config_summaries[[config_id]] <- config_summary
  }
  
  # Overall summary by model
  model_summaries <- list()
  
  for(model in batch_results$models) {
    # Find simulations for this model
    model_sims <- Filter(function(sim) sim$model == model, sims)
    
    if(length(model_sims) > 0) {
      # Extract tip states from all replicates
      tip_means <- sapply(model_sims, function(sim) mean(sim$tip_states))
      tip_sds <- sapply(model_sims, function(sim) sd(sim$tip_states))
      
      # Calculate summary statistics
      model_summary <- list(
        model = model,
        n_simulations = length(model_sims),
        tip_state_mean = mean(tip_means),
        tip_state_sd = mean(tip_sds),
        variability_between_simulations = sd(tip_means)
      )
      
      model_summaries[[model]] <- model_summary
    }
  }
  
  # Return all summaries
  return(list(
    by_config = config_summaries,
    by_model = model_summaries
  ))
}

// ...existing code...

#===============================================================================
# Simulation Visualization Functions
#===============================================================================

#' Visualize simulated chromosome evolution on a tree
#' 
#' Creates a visualization of simulated chromosome changes on a phylogenetic tree
#' 
#' @param sim_result Simulation result object
#' @param show_events Whether to display discrete events
#' @param event_types Which event types to show (NULL for all)
#' @param color_palette Color palette for chromosome counts
#' @param node_labels Whether to show node labels
#' @param output_file Output file path (NULL for on-screen display)
#' @return Invisibly returns the plot object
#' @export
plot_simulated_evolution <- function(sim_result, 
                                    show_events = TRUE,
                                    event_types = c("wgd", "fusion", "fission"),
                                    color_palette = "viridis",
                                    node_labels = TRUE,
                                    output_file = NULL) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggplot2 and ggtree packages are required for visualization")
  }
  
  message("Visualizing simulated chromosome evolution...")
  
  # Extract data
  tree <- sim_result$tree
  tip_states <- sim_result$tip_states
  node_states <- sim_result$node_states
  model <- sim_result$model
  events <- sim_result$events
  
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
  
  # Create node data
  n_tips <- length(tip_states)
  node_ids <- (n_tips + 1):(n_tips + length(node_states))
  
  node_data <- data.frame(
    node = c(1:n_tips, node_ids),
    state = c(tip_states, node_states),
    stringsAsFactors = FALSE
  )
  
  # Create base tree
  p <- ggtree::ggtree(tree, ladderize = TRUE) %<+% node_data
  
  # Add node colors by chromosome count
  if(color_palette == "viridis") {
    p <- p + ggtree::aes(color = state) +
      ggplot2::scale_color_viridis_c(name = "Chromosome\nCount", option = "plasma")
  } else if(color_palette == "rainbow") {
    p <- p + ggtree::aes(color = state) +
      ggplot2::scale_color_gradientn(name = "Chromosome\nCount", 
                                  colors = rainbow(10))
  } else {
    p <- p + ggtree::aes(color = state) +
      ggplot2::scale_color_gradient(name = "Chromosome\nCount", 
                                 low = "blue", high = "red")
  }
  
  # Add node labels if requested
  if(node_labels) {
    p <- p + ggtree::geom_nodelab(aes(label = round(state, 1)), size = 3, color = "black", 
                               fontface = "bold", nudge_x = 0.2)
  }
  
  # Add tip labels
  p <- p + ggtree::geom_tiplab(size = 3)
  
  # Add events if requested
  if(show_events && !is.null(events) && nrow(events) > 0) {
    # Filter by event types if specified
    if(!is.null(event_types)) {
      events <- events[events$type %in% event_types, ]
    }
    
    if(nrow(events) > 0) {
      # Define event colors and shapes
      event_colors <- c(
        "fusion" = "blue",
        "fission" = "red",
        "wgd" = "purple"
      )
      
      event_shapes <- c(
        "fusion" = 25,  # Down triangle
        "fission" = 24, # Up triangle
        "wgd" = 23      # Diamond
      )
      
      # Add event markers
      for(i in 1:nrow(events)) {
        event <- events[i, ]
        
        # Add event marker on the edge
        # This is a simplified approach - ideally we would place it at the exact time
        child_node <- event$child
        event_type <- event$type
        
        if(event_type %in% names(event_colors)) {
          # Add marker at the node
          p <- p + ggtree::geom_point2(aes(subset = (node == child_node)), 
                                     color = event_colors[event_type], 
                                     shape = event_shapes[event_type],
                                     size = 3)
        }
      }
      
      # Add legend for events
      used_event_types <- unique(events$type)
      used_event_types <- used_event_types[used_event_types %in% names(event_colors)]
      
      if(length(used_event_types) > 0) {
        # Create custom legend entries
        for(type in used_event_types) {
          p <- p + ggplot2::annotate(
            "point", 
            x = -Inf, y = -Inf, 
            color = event_colors[type], 
            shape = event_shapes[type],
            size = 3
          )
        }
        
        # Add legend
        p <- p + ggplot2::guides(
          color = ggplot2::guide_colorbar(title = "Chromosome Count"),
          shape = ggplot2::guide_legend(title = "Events", 
                                    override.aes = list(
                                      shape = unname(event_shapes[used_event_types]),
                                      color = unname(event_colors[used_event_types])
                                    ))
        )
      }
    }
  }
  
  # Add simulation info and title
  p <- p + ggplot2::labs(
    title = paste("Simulated Chromosome Evolution -", model, "Model"),
    subtitle = paste("Root value:", sim_result$root_value, 
                    "| Parameters:", 
                    paste(names(sim_result$params), sim_result$params, 
                          sep = "=", collapse = ", "))
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

#' Plot chromosome count distribution from simulations
#' 
#' @param sim_result Simulation result or batch simulation results
#' @param compare_to Optional real data to compare with
#' @param by_model Whether to group by model (for batch results)
#' @param output_file Output file path
#' @return Invisibly returns the plot object
#' @export
plot_chromosome_distribution <- function(sim_result, 
                                       compare_to = NULL,
                                       by_model = TRUE,
                                       output_file = NULL) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }
  
  # Setup output device if needed
  if(!is.null(output_file)) {
    if(grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = 8, height = 6)
    } else if(grepl("\\.png$", output_file)) {
      png(output_file, width = 1600, height = 1200, res = 200)
    } else {
      # Default to PDF
      pdf(paste0(output_file, ".pdf"), width = 8, height = 6)
    }
  }
  
  # Determine if this is a single simulation or batch
  is_batch <- FALSE
  if(!is.null(sim_result$simulations)) {
    is_batch <- TRUE
    message("Plotting chromosome distribution from batch simulations...")
  } else {
    message("Plotting chromosome distribution from simulation...")
  }
  
  if(is_batch) {
    # Prepare data from batch simulations
    if(by_model) {
      # Group by model
      models <- unique(sapply(sim_result$simulations, function(s) s$model))
      plot_data <- data.frame()
      
      for(model in models) {
        model_sims <- Filter(function(s) s$model == model, sim_result$simulations)
        
        # Extract chromosome counts from all simulations of this model
        for(sim in model_sims) {
          sim_data <- data.frame(
            model = model,
            counts = sim$tip_states,
            stringsAsFactors = FALSE
          )
          plot_data <- rbind(plot_data, sim_data)
        }
      }
    } else {
      # Combine all simulations
      plot_data <- data.frame()
      
      for(sim in sim_result$simulations) {
        sim_data <- data.frame(
          model = sim$model,
          counts = sim$tip_states,
          stringsAsFactors = FALSE
        )
        plot_data <- rbind(plot_data, sim_data)
      }
    }
    
    # Create density plot
    if(by_model) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = counts, fill = model)) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::labs(
          title = "Chromosome Count Distributions from Simulations",
          x = "Chromosome Count",
          y = "Density"
        ) +
        ggplot2::theme_minimal()
    } else {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = counts)) +
        ggplot2::geom_density(fill = "steelblue", alpha = 0.7) +
        ggplot2::labs(
          title = "Chromosome Count Distribution from All Simulations",
          x = "Chromosome Count",
          y = "Density"
        ) +
        ggplot2::theme_minimal()
    }
    
  } else {
    # Single simulation
    plot_data <- data.frame(
      counts = sim_result$tip_states,
      stringsAsFactors = FALSE
    )
    
    # Create histogram
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = counts)) +
      ggplot2::geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
      ggplot2::labs(
        title = paste("Chromosome Count Distribution -", sim_result$model, "Model"),
        x = "Chromosome Count",
        y = "Frequency"
      ) +
      ggplot2::theme_minimal()
  }
  
  # Add real data comparison if provided
  if(!is.null(compare_to)) {
    real_data <- data.frame(
      counts = compare_to,
      type = "observed",
      stringsAsFactors = FALSE
    )
    
    # Add density of real data
    p <- p + ggplot2::geom_density(data = real_data, ggplot2::aes(x = counts, fill = "Observed"), 
                                alpha = 0.5) +
      ggplot2::scale_fill_manual(values = c("Observed" = "red"), 
                              name = "Data Source")
  }
  
  # Display plot
  print(p)
  
  # Close device if using file output
  if(!is.null(output_file)) {
    dev.off()
    message(paste("Plot saved to:", output_file))
  }
  
  # Return plot object
  return(invisible(p))
}

#' Create animation of chromosome evolution through time
#' 
#' @param sim_result Simulation result object
#' @param n_frames Number of frames in animation
#' @param output_file Output GIF file path
#' @return Invisibly returns animation object
#' @export
animate_chromosome_evolution <- function(sim_result, n_frames = 30, output_file = "evolution.gif") {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE) || 
     !requireNamespace("ggtree", quietly = TRUE) || 
     !requireNamespace("gganimate", quietly = TRUE)) {
    stop("ggplot2, ggtree, and gganimate packages are required for animation")
  }
  
  message("Creating chromosome evolution animation...")
  
  # Extract data
  tree <- sim_result$tree
  tip_states <- sim_result$tip_states
  node_states <- sim_result$node_states
  events <- sim_result$events
  
  # Calculate tree height for time scaling
  tree_height <- max(phytools::nodeHeights(tree))
  
  # Create time slices
  time_points <- seq(0, tree_height, length.out = n_frames)
  
  # Create a list to store each frame
  frames <- list()
  
  # For each time point, create a frame showing the state of evolution
  for(t in 1:length(time_points)) {
    current_time <- time_points[t]
    
    # Calculate states at current time
    # This requires using edge data to interpolate states along branches
    current_states <- vector("numeric", length = Ntip(tree) + Nnode(tree))
    
    # Initialize with root state
    current_states[Ntip(tree) + 1] <- node_states[1]
    
    # Process each edge to update current states
    for(i in 1:nrow(tree$edge)) {
      parent <- tree$edge[i, 1]
      child <- tree$edge[i, 2]
      edge_length <- tree$edge.length[i]
      
      # Get parent and child states
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
      
      # Get edge height data
      heights <- phytools::nodeHeights(tree)
      edge_start <- heights[i, 1]
      edge_end <- heights[i, 2]
      
      # Check if this edge is active at current time
      if(current_time >= edge_start && current_time <= edge_end) {
        # Calculate position along edge
        pos <- (current_time - edge_start) / (edge_end - edge_start)
        
        # Interpolate state
        current_state <- parent_state + (child_state - parent_state) * pos
        
        # Update current state for child node
        current_states[child] <- current_state
      } else if(current_time > edge_end) {
        # Edge has been traversed, use child state
        current_states[child] <- child_state
      }
    }
    
    # Create data frame for this time point
    frame_data <- data.frame(
      node = 1:(Ntip(tree) + Nnode(tree)),
      state = current_states,
      time = current_time,
      stringsAsFactors = FALSE
    )
    
    frames[[t]] <- frame_data
  }
  
  # Combine all frames
  animation_data <- do.call(rbind, frames)
  
  # Create initial tree plot
  p <- ggtree::ggtree(tree, ladderize = TRUE)
  
  # Add time transition
  p <- p + ggtree::aes(color = state) + 
    gganimate::transition_time(time) +
    ggplot2::scale_color_viridis_c(name = "Chromosome\nCount", option = "plasma") +
    ggplot2::labs(title = "Chromosome Evolution: Time {frame_time}") +
    ggplot2::theme(legend.position = "right")
  
  # Create animation
  anim <- gganimate::animate(
    p + ggtree::geom_tiplab(), 
    nframes = n_frames, 
    fps = 10, 
    width = 800, 
    height = 600,
    renderer = gganimate::gifski_renderer(output_file)
  )
  
  message(paste("Animation saved to:", output_file))
  
  return(invisible(anim))
}

#===============================================================================
# Simulation Validation and Testing Functions
#===============================================================================

#' Compare simulated data to empirical data
#' 
#' Tests how well simulated distributions match empirical chromosome counts
#' 
#' @param sim_result Simulation result or batch simulation results
#' @param empirical_counts Empirical chromosome counts for comparison
#' @param test_method Statistical test method: "ks" (Kolmogorov-Smirnov), 
#'                    "chisq" (Chi-squared), or "mean_diff" (Mean difference)
#' @return Comparison results
#' @export
compare_to_empirical <- function(sim_result, 
                               empirical_counts, 
                               test_method = "ks") {
  message("Comparing simulated data to empirical chromosome counts...")
  
  # Check if this is a batch simulation
  is_batch <- FALSE
  if(!is.null(sim_result$simulations)) {
    is_batch <- TRUE
  }
  
  # Prepare output structure
  comparison <- list(
    method = test_method,
    empirical_stats = list(
      n = length(empirical_counts),
      mean = mean(empirical_counts),
      median = median(empirical_counts),
      sd = sd(empirical_counts),
      range = range(empirical_counts)
    )
  )
  
  if(is_batch) {
    # Process each simulation in batch
    models <- unique(sapply(sim_result$simulations, function(s) s$model))
    
    model_results <- list()
    for(model in models) {
      # Combine tip states from all simulations of this model
      model_sims <- Filter(function(s) s$model == model, sim_result$simulations)
      
      all_tip_states <- numeric(0)
      for(sim in model_sims) {
        all_tip_states <- c(all_tip_states, sim$tip_states)
      }
      
      # Run statistical test
      test_result <- compare_distributions(all_tip_states, empirical_counts, test_method)
      
      # Store results for this model
      model_results[[model]] <- list(
        stats = list(
          n = length(all_tip_states),
          mean = mean(all_tip_states),
          median = median(all_tip_states),
          sd = sd(all_tip_states),
          range = range(all_tip_states)
        ),
        test = test_result
      )
    }
    
    comparison$model_results <- model_results
    
    # Find best fitting model
    if(test_method == "ks") {
      # Lower D-statistic is better
      test_stats <- sapply(model_results, function(x) x$test$statistic)
      best_model <- names(which.min(test_stats))
    } else if(test_method == "chisq") {
      # Lower chi-squared is better
      test_stats <- sapply(model_results, function(x) x$test$statistic)
      best_model <- names(which.min(test_stats))
    } else if(test_method == "mean_diff") {
      # Lower mean difference is better
      test_stats <- sapply(model_results, function(x) x$test$mean_diff)
      best_model <- names(which.min(test_stats))
    } else {
      # Default to p-value
      p_values <- sapply(model_results, function(x) x$test$p.value)
      best_model <- names(which.max(p_values))
    }
    
    comparison$best_model <- best_model
    
  } else {
    # Single simulation
    # Run statistical test
    test_result <- compare_distributions(sim_result$tip_states, empirical_counts, test_method)
    
    # Store results
    comparison$simulation_stats <- list(
      n = length(sim_result$tip_states),
      mean = mean(sim_result$tip_states),
      median = median(sim_result$tip_states),
      sd = sd(sim_result$tip_states),
      range = range(sim_result$tip_states)
    )
    comparison$test <- test_result
  }
  
  return(comparison)
}

#' Compare two distributions using statistical tests
#' 
#' @param simulated Simulated values
#' @param observed Observed values
#' @param test_method Statistical test method
#' @return Test results
#' @keywords internal
compare_distributions <- function(simulated, observed, test_method = "ks") {
  if(test_method == "ks") {
    # Kolmogorov-Smirnov test
    test_result <- stats::ks.test(simulated, observed)
    
    return(list(
      method = "Kolmogorov-Smirnov",
      statistic = test_result$statistic,
      p.value = test_result$p.value
    ))
    
  } else if(test_method == "chisq") {
    # Chi-squared test on binned data
    # Create bins that cover both distributions
    all_data <- c(simulated, observed)
    bins <- seq(floor(min(all_data)), ceiling(max(all_data)) + 1)
    
    # Count observations in each bin
    simulated_counts <- hist(simulated, breaks = bins, plot = FALSE)$counts
    observed_counts <- hist(observed, breaks = bins, plot = FALSE)$counts
    
    # Add small constant to avoid zeros
    simulated_counts <- simulated_counts + 0.1
    observed_counts <- observed_counts + 0.1
    
    # Scale to match total counts
    scale_factor <- sum(observed_counts) / sum(simulated_counts)
    simulated_counts <- simulated_counts * scale_factor
    
    # Run chi-squared test
    test_result <- chisq.test(rbind(simulated_counts, observed_counts))
    
    return(list(
      method = "Chi-squared",
      statistic = test_result$statistic,
      p.value = test_result$p.value
    ))
    
  } else if(test_method == "mean_diff") {
    # Simple difference in means
    mean_diff <- abs(mean(simulated) - mean(observed))
    sd_diff <- abs(sd(simulated) - sd(observed))
    
    # Calculate p-value based on bootstrapping
    n_bootstrap <- 1000
    combined_data <- c(simulated, observed)
    
    bootstrap_diffs <- numeric(n_bootstrap)
    
    for(i in 1:n_bootstrap) {
      # Shuffle data and split into two groups
      shuffled <- sample(combined_data)
      group1 <- shuffled[1:length(simulated)]
      group2 <- shuffled[(length(simulated)+1):length(shuffled)]
      
      # Calculate mean difference for this bootstrap sample
      bootstrap_diffs[i] <- abs(mean(group1) - mean(group2))
    }
    
    # Calculate p-value as proportion of bootstrap differences >= observed difference
    p_value <- mean(bootstrap_diffs >= mean_diff)
    
    return(list(
      method = "Mean difference",
      mean_diff = mean_diff,
      sd_diff = sd_diff,
      p.value = p_value
    ))
  }
}

#' Generate test data for method validation
#' 
#' Creates test datasets with known properties for evaluating reconstruction methods
#' 
#' @param tree Phylogenetic tree
#' @param scenario Test scenario: "gradual", "jumps", "clade_specific", or "rate_shift"
#' @param params Scenario-specific parameters
#' @param root_value Starting chromosome number
#' @return Simulated test data
#' @export
generate_test_data <- function(tree, 
                             scenario = "gradual", 
                             params = NULL,
                             root_value = 12) {
  message(sprintf("Generating test data: %s scenario", scenario))
  
  # Set default parameters if not provided
  if(is.null(params)) {
    if(scenario == "gradual") {
      params <- list(sigma = 0.5, trend = 0.1)
    } else if(scenario == "jumps") {
      params <- list(jump_rate = 0.2, jump_size = 3, background_rate = 0.1)
    } else if(scenario == "clade_specific") {
      params <- list(clades = NULL, rates = NULL)
    } else if(scenario == "rate_shift") {
      params <- list(shift_node = NULL, rate_before = 0.5, rate_after = 2)
    }
  }
  
  # Apply appropriate simulation model
  if(scenario == "gradual") {
    # Gradual evolution - use BM with trend
    sim_result <- simulate_bm(tree, root_value, params, list(min_value = 1))
    sim_model <- "BM"
  
  } else if(scenario == "jumps") {
    # Punctuated evolution - use jumps model
    jump_params <- list(
      background_rate = params$background_rate,
      jump_rate = params$jump_rate,
      jump_size_mean = 0,
      jump_size_sd = params$jump_size,
      jump_distribution = "normal"
    )
    sim_result <- simulate_jumps(tree, root_value, jump_params, list(min_value = 1))
    sim_model <- "jumps"
    
  } else if(scenario == "clade_specific") {
    # Clade-specific rates - implement custom simulation
    sim_result <- simulate_clade_specific(tree, root_value, params)
    sim_model <- "clade_specific"
    
  } else if(scenario == "rate_shift") {
    # Rate shift at a specific node - implement custom simulation
    sim_result <- simulate_rate_shift(tree, root_value, params)
    sim_model <- "rate_shift"
    
  } else {
    stop("Unsupported test scenario:", scenario)
  }
  
  # Discretize chromosome counts
  sim_result$tip_states <- round(sim_result$tip_states)
  sim_result$node_states <- round(sim_result$node_states)
  
  # Compile result
  result <- list(
    tree = tree,
    tip_states = sim_result$tip_states,
    node_states = sim_result$node_states,
    events = sim_result$events,
    model = sim_model,
    scenario = scenario,
    params = params,
    root_value = root_value
  )
  
  # Add some statistics
  result$statistics <- calculate_simulation_stats(result)
  
  return(result)
}

#' Simulate chromosome evolution with clade-specific rates
#' 
#' @param tree Phylogenetic tree
#' @param root_value Root state
#' @param params Parameters including clades and rates
#' @return Simulation results
#' @keywords internal
simulate_clade_specific <- function(tree, root_value, params) {
  # Extract parameters
  clades <- params$clades
  rates <- params$rates
  
  # If clades not specified, create random clade assignments
  if(is.null(clades)) {
    # Randomly select internal nodes to define clades
    n_clades <- max(2, ceiling(tree$Nnode / 5))  # Ensure at least 2 clades
    clade_nodes <- sample((Ntip(tree) + 1):(Ntip(tree) + tree$Nnode), n_clades)
    
    # For each clade node, identify all descendant tips
    clades <- list()
    for(i in 1:length(clade_nodes)) {
      descendants <- phytools::getDescendants(tree, clade_nodes[i])
      # Filter to keep only tips
      tips <- descendants[descendants <= Ntip(tree)]
      
      if(length(tips) > 0) {
        clades[[i]] <- tree$tip.label[tips]
      }
    }
    
    # Remove empty clades
    clades <- clades[sapply(clades, length) > 0]
  }
  
  # If rates not specified, assign random rates to clades
  if(is.null(rates)) {
    rates <- runif(length(clades), 0.2, 2.0)
  }
  
  # Setup simulation containers
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  
  # Initialize all states
  all_states <- numeric(n_total)
  all_states[n_tips + 1] <- root_value  # Set root state
  
  # Initialize events tracking
  events <- data.frame(
    edge = integer(0),
    parent = integer(0),
    child = integer(0),
    time = numeric(0),
    type = character(0),
    magnitude = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Create mapping of tip indices to clade rates
  tip_rates <- rep(1.0, n_tips)  # Default rate
  
  for(i in 1:length(clades)) {
    clade_tips <- clades[[i]]
    clade_rate <- rates[i]
    
    # Find tip indices for this clade
    tip_indices <- match(clade_tips, tree$tip.label)
    tip_indices <- tip_indices[!is.na(tip_indices)]
    
    if(length(tip_indices) > 0) {
      tip_rates[tip_indices] <- clade_rate
    }
  }
  
  # Simulate along each edge
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent state
    parent_state <- all_states[parent]
    
    # Determine rate for this edge
    if(child <= n_tips) {
      # Edge leading to a tip - use tip rate
      edge_rate <- tip_rates[child]
    } else {
      # Internal edge - calculate based on descendant tips
      descendants <- phytools::getDescendants(tree, child)
      tip_descendants <- descendants[descendants <= Ntip(tree)]
      
      if(length(tip_descendants) > 0) {
        edge_rate <- mean(tip_rates[tip_descendants])
      } else {
        edge_rate <- 1.0  # Default rate
      }
    }
    
    # Simulate state change using clade-specific rate
    random_change <- rnorm(1, mean = 0, sd = edge_rate * sqrt(branch_length))
    
    # Calculate child state
    child_state <- parent_state + random_change
    
    # Ensure non-negative value with floor at 1
    child_state <- max(1, child_state)
    
    # Store child state
    all_states[child] <- child_state
  }
  
  # Extract tip and node states
  tip_states <- all_states[1:n_tips]
  names(tip_states) <- tree$tip.label
  node_states <- all_states[(n_tips+1):n_total]
  
  return(list(
    tip_states = tip_states,
    node_states = node_states,
    events = events
  ))
}

#' Simulate chromosome evolution with a rate shift at a specific node
#' 
#' @param tree Phylogenetic tree
#' @param root_value Root state
#' @param params Parameters including shift node and rates
#' @return Simulation results
#' @keywords internal
simulate_rate_shift <- function(tree, root_value, params) {
  # Extract parameters
  shift_node <- params$shift_node
  rate_before <- params$rate_before
  rate_after <- params$rate_after
  
  # If shift node not specified, select a random internal node
  if(is.null(shift_node)) {
    internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)
    
    # Prefer a node that divides the tree somewhat evenly
    # Get node depths
    node_depths <- node.depth(tree)
    
    # Filter to nodes in the middle third of the tree
    tree_height <- max(node_depths)
    mid_height <- tree_height / 2
    mid_range <- tree_height / 3
    
    candidate_nodes <- internal_nodes[
      node_depths[internal_nodes] > (mid_height - mid_range) & 
      node_depths[internal_nodes] < (mid_height + mid_range)
    ]
    
    if(length(candidate_nodes) > 0) {
      shift_node <- sample(candidate_nodes, 1)
    } else {
      # Fallback to any internal node except root
      shift_node <- sample(internal_nodes[-1], 1)
    }
  }
  
  # Setup simulation containers
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_total <- n_tips + n_nodes
  
  # Initialize all states
  all_states <- numeric(n_total)
  all_states[n_tips + 1] <- root_value  # Set root state
  
  # Initialize events tracking
  events <- data.frame(
    edge = integer(0),
    parent = integer(0),
    child = integer(0),
    time = numeric(0),
    type = character(0),
    magnitude = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Get all descendants of the shift node
  shift_descendants <- c(shift_node, phytools::getDescendants(tree, shift_node))
  
  # Simulate along each edge
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    branch_length <- tree$edge.length[i]
    
    # Get parent state
    parent_state <- all_states[parent]
    
    # Determine rate for this edge
    if(child %in% shift_descendants) {
      # After shift - use higher/lower rate
      edge_rate <- rate_after
    } else {
      # Before shift - use baseline rate
      edge_rate <- rate_before
    }
    
    # Simulate state change using appropriate rate
    random_change <- rnorm(1, mean = 0, sd = edge_rate * sqrt(branch_length))
    
    # Calculate child state
    child_state <- parent_state + random_change
    
    # Ensure non-negative value with floor at 1
    child_state <- max(1, child_state)
    
    # Store child state
    all_states[child] <- child_state
    
    # Add a rate shift "event" at the shift node
    if(child == shift_node) {
      events <- rbind(events, data.frame(
        edge = i,
        parent = parent,
        child = child,
        time = 0,  # Placeholder
        type = "rate_shift",
        magnitude = rate_after / rate_before,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Extract tip and node states
  tip_states <- all_states[1:n_tips]
  names(tip_states) <- tree$tip.label
  node_states <- all_states[(n_tips+1):n_total]
  
  return(list(
    tip_states = tip_states,
    node_states = node_states,
    events = events
  ))
}

#' Test reconstruction method against known simulated data
#' 
#' Evaluates reconstruction method performance on datasets with known properties
#' 
#' @param reconstruct_method Function that takes tree and tip_states and returns reconstructed states
#' @param method_name Name of method being tested
#' @param tree Phylogenetic tree (NULL to generate random tree)
#' @param scenarios Vector of test scenarios to run
#' @param n_replicates Number of replicates per scenario
#' @return Test results
#' @export
test_reconstruction_method <- function(reconstruct_method, 
                                    method_name = "custom",
                                    tree = NULL, 
                                    scenarios = c("gradual", "jumps", "clade_specific", "rate_shift"),
                                    n_replicates = 5) {
  message(sprintf("Testing %s reconstruction method against simulated data...", method_name))
  
  # Generate tree if not provided
  if(is.null(tree)) {
    tree <- ape::rtree(30)
  }
  
  # Initialize results container
  test_results <- list(
    method = method_name,
    scenarios = scenarios,
    n_replicates = n_replicates,
    replicates = list()
  )
  
  # Run tests for each scenario
  for(scenario in scenarios) {
    message(sprintf("  Testing %s scenario...", scenario))
    
    scenario_results <- list()
    
    for(rep in 1:n_replicates) {
      # Generate test data
      test_data <- generate_test_data(tree, scenario)
      
      # Apply reconstruction method
      start_time <- Sys.time()
      reconstruction <- tryCatch({
        reconstruct_method(tree, test_data$tip_states)
      }, error = function(e) {
        message(sprintf("    Error in reconstruction: %s", e$message))
        return(NULL)
      })
      end_time <- Sys.time()
      runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      # Skip if reconstruction failed
      if(is.null(reconstruction)) {
        next
      }
      
      # Extract true and reconstructed states for comparison
      true_states <- test_data$node_states
      
      # Extract reconstructed states - handle different result formats
      reconstructed_states <- NULL
      
      if(is.list(reconstruction) && !is.null(reconstruction$ancestral_states)) {
        # Standard format with ancestral_states data frame
        anc_states <- reconstruction$ancestral_states
        
        # Check if node_id column exists
        if("node_id" %in% colnames(anc_states)) {
          # Ensure ordered correctly
          node_ids <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
          reconstructed_states <- numeric(length(node_ids))
          
          for(i in 1:length(node_ids)) {
            row_idx <- which(anc_states$node_id == node_ids[i])
            if(length(row_idx) == 1) {
              reconstructed_states[i] <- anc_states$state[row_idx]
            }
          }
        } else {
          # Assume states are in the correct order
          reconstructed_states <- anc_states$state
        }
      } else if(is.numeric(reconstruction)) {
        # Direct vector of states
        reconstructed_states <- reconstruction
      } else {
        message("    Unsupported reconstruction result format")
        next
      }
      
      # Calculate error metrics
      errors <- true_states - reconstructed_states
      mae <- mean(abs(errors))
      rmse <- sqrt(mean(errors^2))
      correlation <- cor(true_states, reconstructed_states, method = "spearman")
      
      # Store results for this replicate
      scenario_results[[rep]] <- list(
        test_data = test_data,
        reconstruction = reconstruction,
        errors = errors,
        metrics = list(
          mae = mae,
          rmse = rmse,
          correlation = correlation,
          runtime = runtime
        )
      )
      
      message(sprintf("    Replicate %d: MAE = %.3f, RMSE = %.3f, Corr = %.3f", 
                     rep, mae, rmse, correlation))
    }
    
    # Calculate average metrics for this scenario
    if(length(scenario_results) > 0) {
      mae_values <- sapply(scenario_results, function(x) x$metrics$mae)
      rmse_values <- sapply(scenario_results, function(x) x$metrics$rmse)
      corr_values <- sapply(scenario_results, function(x) x$metrics$correlation)
      runtime_values <- sapply(scenario_results, function(x) x$metrics$runtime)
      
      scenario_avg <- list(
        mae = list(mean = mean(mae_values), sd = sd(mae_values)),
        rmse = list(mean = mean(rmse_values), sd = sd(rmse_values)),
        correlation = list(mean = mean(corr_values), sd = sd(corr_values)),
        runtime = list(mean = mean(runtime_values), sd = sd(runtime_values))
      )
      
      message(sprintf("  %s scenario results: Mean MAE = %.3f, Mean RMSE = %.3f, Mean Corr = %.3f",
                     scenario, scenario_avg$mae$mean, scenario_avg$rmse$mean, 
                     scenario_avg$correlation$mean))
    } else {
      scenario_avg <- NULL
      message(sprintf("  %s scenario: No successful reconstructions", scenario))
    }
    
    # Store results for this scenario
    test_results$replicates[[scenario]] <- list(
      results = scenario_results,
      summary = scenario_avg
    )
  }
  
  # Create overall summary
  all_mae <- numeric(0)
  all_rmse <- numeric(0)
  all_corr <- numeric(0)
  all_runtime <- numeric(0)
  
  for(scenario in scenarios) {
    if(!is.null(test_results$replicates[[scenario]]$summary)) {
      all_mae <- c(all_mae, sapply(test_results$replicates[[scenario]]$results, 
                                  function(x) x$metrics$mae))
      all_rmse <- c(all_rmse, sapply(test_results$replicates[[scenario]]$results, 
                                    function(x) x$metrics$rmse))
      all_corr <- c(all_corr, sapply(test_results$replicates[[scenario]]$results, 
                                    function(x) x$metrics$correlation))
      all_runtime <- c(all_runtime, sapply(test_results$replicates[[scenario]]$results, 
                                         function(x) x$metrics$runtime))
    }
  }
  
  test_results$overall <- list(
    mae = list(mean = mean(all_mae), sd = sd(all_mae)),
    rmse = list(mean = mean(all_rmse), sd = sd(all_rmse)),
    correlation = list(mean = mean(all_corr), sd = sd(all_corr)),
    runtime = list(mean = mean(all_runtime), sd = sd(all_runtime))
  )
  
  message(sprintf("Overall method performance: Mean MAE = %.3f, Mean RMSE = %.3f, Mean Corr = %.3f",
                 test_results$overall$mae$mean, test_results$overall$rmse$mean,
                 test_results$overall$correlation$mean))
  
  return(test_results)
}
