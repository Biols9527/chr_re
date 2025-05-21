#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Bayesian Reconstruction Module
# Author: Bioinformatics Team
# Date: 2025-03-22
# Description: Implements Bayesian MCMC methods for ancestral chromosome
#              number reconstruction with parameter uncertainty quantification
#===============================================================================

suppressPackageStartupMessages({
  library(ape)        # For phylogenetic operations
  library(phytools)   # For phylogenetic utilities  
  library(coda)       # For MCMC diagnostics
  library(parallel)   # For parallel computation
  library(rjags)      # For JAGS interface
  library(MCMCpack)   # For additional MCMC tools
})

#===============================================================================
# Bayesian Reconstruction Core Functions
#===============================================================================

#' Initialize Bayesian Run (Internal Helper)
#'
#' Performs initial setup, input validation, default settings, and tree pruning.
#'
#' @param tree Phylogenetic tree object
#' @param chr_counts Chromosome counts, named vector with species names
#' @param model Evolutionary model string
#' @param mcmc_settings List of MCMC settings or NULL
#' @param priors List of prior distributions or NULL
#' @param allow_jumps Logical, whether to allow jumps (for JumpBM model)
#' @return A list containing: `pruned_tree`, `ordered_counts`, 
#'         `effective_mcmc_settings`, `effective_priors`.
#' @keywords internal
initialize_bayesian_run_internal <- function(tree, chr_counts, model, 
                                           mcmc_settings, priors, allow_jumps) {
  # Check input
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  if (is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector with species names")
  }
  supported_models <- c("BM", "OU", "JumpBM", "RateShift")
  if (!model %in% supported_models) {
    stop(paste("Unsupported model. Choose from:", paste(supported_models, collapse = ", ")))
  }
  if (!requireNamespace("rjags", quietly = TRUE)) {
    stop("Package 'rjags' and JAGS software are required for Bayesian analysis")
  }

  message(sprintf("Initializing Bayesian %s model run...", model))

  # Set default MCMC settings if not provided
  effective_mcmc_settings <- mcmc_settings
  if (is.null(effective_mcmc_settings)) {
    effective_mcmc_settings <- list(
      iterations = 50000,    # Total MCMC iterations
      burnin = 10000,        # Burn-in iterations to discard
      thinning = 20,         # Thinning interval
      chains = 4             # Number of MCMC chains
    )
  }
  
  # Set default priors if not provided
  effective_priors <- priors
  if (is.null(effective_priors)) {
    effective_priors <- get_default_priors(model)
  }
  
  # Ensure tree and data contain same species
  common_species <- intersect(names(chr_counts), tree$tip.label)
  if (length(common_species) < 3) {
    stop("Fewer than 3 species in common, cannot perform Bayesian reconstruction")
  }
  
  # Prune tree to match data
  pruned_tree <- ape::keep.tip(tree, common_species)
  
  # Order data to match tree
  ordered_counts <- chr_counts[pruned_tree$tip.label]
  
  return(list(
    pruned_tree = pruned_tree,
    ordered_counts = ordered_counts,
    effective_mcmc_settings = effective_mcmc_settings,
    effective_priors = effective_priors
  ))
}

#' Prepare JAGS Input (Internal Helper)
#'
#' Prepares the JAGS model string, data list, and parameters to monitor.
#'
#' @param pruned_tree Pruned phylogenetic tree object
#' @param ordered_counts Ordered chromosome counts for tips
#' @param model Evolutionary model string
#' @param effective_priors List of effective prior distributions
#' @param allow_jumps Logical, whether to allow jumps
#' @return A list containing: `jags_model_text`, `jags_data_list`, `monitor_parameters`.
#' @keywords internal
prepare_jags_input_internal <- function(pruned_tree, ordered_counts, model, 
                                      effective_priors, allow_jumps) {
  
  # Create tree structure data for JAGS model
  tree_data <- prepare_tree_data(pruned_tree) # Existing helper
  
  # Prepare JAGS model string
  jags_model_text <- build_bayesian_model(model, tree_data, effective_priors, allow_jumps) # Existing helper
  
  # Prepare data for JAGS
  jags_data_list <- list(
    N_tips = length(ordered_counts),
    N_nodes = pruned_tree$Nnode,
    N_edges = nrow(pruned_tree$edge),
    y = ordered_counts,
    parent = tree_data$parent,
    child = tree_data$child,
    edge_length = tree_data$edge_length,
    is_tip = tree_data$is_tip 
  )
  
  # Add model-specific data to jags_data_list
  if (model == "OU") {
    jags_data_list$node_height <- tree_data$node_height
  } else if (model == "RateShift") {
    jags_data_list$edge_group <- rep(1, nrow(pruned_tree$edge)) # Default: all edges in one group
    jags_data_list$N_groups <- 1
  }
  
  # Parameters to monitor
  monitor_parameters <- c("states", "sigma2")
  if (model == "OU") {
    monitor_parameters <- c(monitor_parameters, "alpha", "theta")
  } else if (model == "JumpBM") {
    monitor_parameters <- c(monitor_parameters, "jump_prob", "jump_size")
  } else if (model == "RateShift") {
    monitor_parameters <- c(monitor_parameters, "rate_multiplier")
  }
  
  return(list(
    jags_model_text = jags_model_text,
    jags_data_list = jags_data_list,
    monitor_parameters = monitor_parameters
  ))
}

#' Run MCMC Sampling (Internal Helper)
#'
#' Executes the MCMC sampling using JAGS, supporting parallel or sequential chains.
#'
#' @param jags_model_text Text string of the JAGS model.
#' @param jags_data_list List of data for JAGS.
#' @param monitor_parameters Vector of parameters to monitor.
#' @param effective_mcmc_settings List of effective MCMC settings.
#' @param use_parallel Logical, whether to run chains in parallel.
#' @param model_name_for_messages Name of the model, for messages.
#' @return The `mcmc.list` object containing the samples.
#' @keywords internal
run_mcmc_sampling_internal <- function(jags_model_text, jags_data_list, monitor_parameters,
                                     effective_mcmc_settings, use_parallel, model_name_for_messages) {
  
  message("Starting MCMC sampling...")
  
  if (use_parallel && effective_mcmc_settings$chains > 1 && requireNamespace("parallel", quietly = TRUE)) {
    n_cores <- min(parallel::detectCores() - 1, effective_mcmc_settings$chains)
    if (n_cores < 1) n_cores <- 1
    
    message(sprintf("Running %d MCMC chains for %s model on %d cores...", 
                    effective_mcmc_settings$chains, model_name_for_messages, n_cores))
    
    cl <- parallel::makeCluster(n_cores)
    # Export all necessary variables and specifically use the environment of this function
    # for jags_model_text, jags_data_list etc. to ensure they are found.
    parallel::clusterExport(cl, c("jags_model_text", "jags_data_list", "monitor_parameters", "effective_mcmc_settings"),
                            envir = environment())
    parallel::clusterEvalQ(cl, {
      library(rjags)
    })
    
    chain_results <- parallel::parLapply(cl, 1:effective_mcmc_settings$chains, function(chain_num) {
      # Each chain needs its own model compilation and adaptation
      jags_instance <- rjags::jags.model(textConnection(jags_model_text), data = jags_data_list,
                                         n.chains = 1, n.adapt = 5000, quiet = TRUE)
      # Burn-in for this specific chain
      update(jags_instance, n.iter = effective_mcmc_settings$burnin, progress.bar = "none")
      # Sample from posterior for this specific chain
      samples <- rjags::coda.samples(jags_instance, variable.names = monitor_parameters,
                                     n.iter = effective_mcmc_settings$iterations,
                                     thin = effective_mcmc_settings$thinning,
                                     progress.bar = "none")
      return(samples[[1]]) # Return the mcmc object from the list
    })
    parallel::stopCluster(cl)
    mcmc_samples <- coda::as.mcmc.list(chain_results)
    
  } else {
    message(sprintf("Running %d MCMC chains for %s model sequentially...", 
                    effective_mcmc_settings$chains, model_name_for_messages))
    
    jags_instance <- rjags::jags.model(textConnection(jags_model_text), data = jags_data_list,
                                       n.chains = effective_mcmc_settings$chains, n.adapt = 5000, quiet = FALSE)
    update(jags_instance, n.iter = effective_mcmc_settings$burnin) # Progress bar will show by default
    mcmc_samples <- rjags::coda.samples(jags_instance, variable.names = monitor_parameters,
                                        n.iter = effective_mcmc_settings$iterations,
                                        thin = effective_mcmc_settings$thinning) # Progress bar
  }
  
  message("MCMC sampling completed.")
  return(mcmc_samples)
}

#' Compile Bayesian Results (Internal Helper)
#'
#' Processes MCMC samples and assembles the final results object.
#'
#' @param mcmc_samples The `mcmc.list` object from `run_mcmc_sampling_internal`.
#' @param pruned_tree Pruned phylogenetic tree object.
#' @param ordered_counts Ordered chromosome counts for tips.
#' @param model Evolutionary model string.
#' @param effective_mcmc_settings List of effective MCMC settings.
#' @param effective_priors List of effective prior distributions.
#' @param original_chr_counts Original chromosome counts (before pruning/ordering).
#' @return The final result list for `reconstruct_chromosomes_bayesian`.
#' @keywords internal
compile_bayesian_results_internal <- function(mcmc_samples, pruned_tree, ordered_counts, model,
                                            effective_mcmc_settings, effective_priors, original_chr_counts) {
  
  # Process MCMC samples using the existing helper
  processed_results <- process_mcmc_samples(mcmc_samples, pruned_tree, ordered_counts, model) # Existing helper
  
  # Assemble the final result list
  final_result <- list(
    ancestral_states = processed_results$ancestral_states,
    tree = pruned_tree,
    method = "Bayesian",
    model = model,
    mcmc_samples = mcmc_samples,
    model_parameters = processed_results$model_parameters,
    convergence = processed_results$convergence,
    tip_states = ordered_counts,
    mcmc_settings = effective_mcmc_settings,
    priors = effective_priors,
    original_data = list(
      chr_counts = original_chr_counts # Use the original, full chr_counts here
    )
  )
  
  # Calculate reconstruction quality metrics using the existing helper
  final_result$quality <- calculate_bayesian_quality(final_result) # Existing helper
  
  return(final_result)
}

#' Reconstruct ancestral chromosome numbers using Bayesian MCMC
#' 
#' Applies Bayesian inference via MCMC to reconstruct ancestral chromosome numbers
#' with full posterior distributions
#' 
#' @param tree Phylogenetic tree object
#' @param chr_counts Chromosome counts, named vector with species names
#' @param model Evolutionary model: "BM" (Brownian Motion), "OU" (Ornstein-Uhlenbeck),
#'              "JumpBM" (Brownian Motion with jumps), or "RateShift" (variable rates)
#' @param mcmc_settings List of MCMC settings (iterations, burnin, thinning, chains)
#' @param priors List of prior distributions for model parameters
#' @param allow_jumps Whether to allow discrete jumps (chromosome events)
#' @param use_parallel Whether to run chains in parallel
#' @return Ancestral chromosome reconstruction results with posterior distributions
#' @export
reconstruct_chromosomes_bayesian <- function(tree, chr_counts, 
                                           model = "BM", 
                                           mcmc_settings = NULL,
                                           priors = NULL,
                                           allow_jumps = TRUE,
                                           use_parallel = TRUE) {

  # Step 1: Initialization (validation, defaults, tree pruning)
  # The message "Using Bayesian %s model..." is part of initialize_bayesian_run_internal
  init_results <- initialize_bayesian_run_internal(
    tree = tree, 
    chr_counts = chr_counts, 
    model = model, 
    mcmc_settings = mcmc_settings, 
    priors = priors, 
    allow_jumps = allow_jumps
  )
  
  # Step 2: Prepare JAGS model string, data, and parameters to monitor
  jags_input <- prepare_jags_input_internal(
    pruned_tree = init_results$pruned_tree,
    ordered_counts = init_results$ordered_counts,
    model = model,
    effective_priors = init_results$effective_priors,
    allow_jumps = allow_jumps
  )
  
  # Step 3: Run MCMC Sampling
  # Messages related to MCMC start and chain execution are in run_mcmc_sampling_internal
  mcmc_samples <- run_mcmc_sampling_internal(
    jags_model_text = jags_input$jags_model_text,
    jags_data_list = jags_input$jags_data_list,
    monitor_parameters = jags_input$monitor_parameters,
    effective_mcmc_settings = init_results$effective_mcmc_settings,
    use_parallel = use_parallel,
    model_name_for_messages = model # Pass model name for more descriptive messages
  )
  # The message "MCMC sampling completed" is part of run_mcmc_sampling_internal
  
  # Step 4: Process MCMC samples and compile final results
  final_bayesian_result <- compile_bayesian_results_internal(
    mcmc_samples = mcmc_samples,
    pruned_tree = init_results$pruned_tree,
    ordered_counts = init_results$ordered_counts,
    model = model,
    effective_mcmc_settings = init_results$effective_mcmc_settings,
    effective_priors = init_results$effective_priors,
    original_chr_counts = chr_counts # Pass original chr_counts for the final result object
  )
  
  return(final_bayesian_result)
}

#' Prepare tree structure data for JAGS model
#' 
#' @param tree Phylogenetic tree
#' @return List of tree structure data
#' @keywords internal
prepare_tree_data <- function(tree) {
  n_tips <- length(tree$tip.label)
  n_nodes <- tree$Nnode
  n_edges <- nrow(tree$edge)
  
  # Create tip indicator
  is_tip <- rep(FALSE, n_tips + n_nodes)
  is_tip[1:n_tips] <- TRUE
  
  # Get node heights if tree has branch lengths
  node_height <- numeric(n_tips + n_nodes)
  if(!is.null(tree$edge.length)) {
    # Calculate node heights (distance from root)
    heights <- phytools::nodeHeights(tree)
    for(i in 1:n_edges) {
      node_height[tree$edge[i, 2]] <- heights[i, 2]
    }
    # Root height is 0
    node_height[n_tips + 1] <- 0
  }
  
  # Prepare parent-child relationships
  parent <- tree$edge[, 1]
  child <- tree$edge[, 2]
  
  # Use edge lengths if available, otherwise use unit lengths
  edge_length <- if(!is.null(tree$edge.length)) tree$edge.length else rep(1, n_edges)
  
  return(list(
    parent = parent,
    child = child,
    edge_length = edge_length,
    is_tip = is_tip,
    node_height = node_height
  ))
}

#' Build Bayesian model for JAGS
#' 
#' @param model Model type
#' @param tree_data Tree structure data
#' @param priors Prior distributions
#' @param allow_jumps Whether to allow jumps
#' @return JAGS model code as string
#' @keywords internal
build_bayesian_model <- function(model, tree_data, priors, allow_jumps) {
  if(model == "BM") {
    return(build_bm_model(priors))
  } else if(model == "OU") {
    return(build_ou_model(priors))
  } else if(model == "JumpBM") {
    return(build_jumpbm_model(priors, allow_jumps))
  } else if(model == "RateShift") {
    return(build_rateshift_model(priors))
  } else {
    stop(paste("Unsupported model for JAGS:", model))
  }
}

#' Build Brownian Motion model for JAGS
#' 
#' @param priors Prior distributions
#' @return JAGS model code as string
#' @keywords internal
build_bm_model <- function(priors) {
  model_code <- "
  model {
    # Priors for model parameters
    sigma2 ~ dgamma(PRIOR_SIGMA_SHAPE, PRIOR_SIGMA_RATE)  # Evolutionary rate
    root_state ~ dnorm(PRIOR_ROOT_MEAN, PRIOR_ROOT_PREC)  # Root state
    
    # Assign root state
    states[N_tips + 1] <- root_state
    
    # Loop through edges to define tree structure
    for(i in 1:N_edges) {
      # Define state changes via Brownian motion
      states[child[i]] ~ dnorm(states[parent[i]], 1/(sigma2 * edge_length[i]))
    }
    
    # Data likelihood for tips
    for(i in 1:N_tips) {
      y[i] ~ dnorm(states[i], 1/0.01)  # Small observational variance
    }
  }
  "
  
  # Replace prior placeholders with actual values
  model_code <- gsub("PRIOR_SIGMA_SHAPE", priors$sigma_shape, model_code)
  model_code <- gsub("PRIOR_SIGMA_RATE", priors$sigma_rate, model_code)
  model_code <- gsub("PRIOR_ROOT_MEAN", priors$root_mean, model_code)
  model_code <- gsub("PRIOR_ROOT_PREC", 1/priors$root_var, model_code)
  
  return(model_code)
}

#' Build Ornstein-Uhlenbeck model for JAGS
#' 
#' @param priors Prior distributions
#' @return JAGS model code as string
#' @keywords internal
build_ou_model <- function(priors) {
  model_code <- "
  model {
    # Priors for model parameters
    sigma2 ~ dgamma(PRIOR_SIGMA_SHAPE, PRIOR_SIGMA_RATE)  # Evolutionary rate
    alpha ~ dgamma(PRIOR_ALPHA_SHAPE, PRIOR_ALPHA_RATE)   # Strength of selection
    theta ~ dnorm(PRIOR_THETA_MEAN, PRIOR_THETA_PREC)     # Optimum value
    root_state ~ dnorm(PRIOR_ROOT_MEAN, PRIOR_ROOT_PREC)  # Root state
    
    # Assign root state
    states[N_tips + 1] <- root_state
    
    # Loop through edges to define tree structure
    for(i in 1:N_edges) {
      # OU parameters
      expected <- theta + (states[parent[i]] - theta) * exp(-alpha * edge_length[i])
      precision <- 1/(sigma2/(2*alpha) * (1 - exp(-2*alpha * edge_length[i])))
      
      # Define state changes via OU process
      states[child[i]] ~ dnorm(expected, precision)
    }
    
    # Data likelihood for tips
    for(i in 1:N_tips) {
      y[i] ~ dnorm(states[i], 1/0.01)  # Small observational variance
    }
  }
  "
  
  # Replace prior placeholders with actual values
  model_code <- gsub("PRIOR_SIGMA_SHAPE", priors$sigma_shape, model_code)
  model_code <- gsub("PRIOR_SIGMA_RATE", priors$sigma_rate, model_code)
  model_code <- gsub("PRIOR_ALPHA_SHAPE", priors$alpha_shape, model_code)
  model_code <- gsub("PRIOR_ALPHA_RATE", priors$alpha_rate, model_code)
  model_code <- gsub("PRIOR_THETA_MEAN", priors$theta_mean, model_code)
  model_code <- gsub("PRIOR_THETA_PREC", 1/priors$theta_var, model_code)
  model_code <- gsub("PRIOR_ROOT_MEAN", priors$root_mean, model_code)
  model_code <- gsub("PRIOR_ROOT_PREC", 1/priors$root_var, model_code)
  
  return(model_code)
}

#' Build Brownian Motion with jumps model for JAGS
#' 
#' @param priors Prior distributions
#' @param allow_jumps Whether to allow jumps
#' @return JAGS model code as string
#' @keywords internal
build_jumpbm_model <- function(priors, allow_jumps) {
  model_code <- "
  model {
    # Priors for model parameters
    sigma2 ~ dgamma(PRIOR_SIGMA_SHAPE, PRIOR_SIGMA_RATE)               # Background evolutionary rate
    jump_prob ~ dbeta(PRIOR_JUMP_PROB_ALPHA, PRIOR_JUMP_PROB_BETA)     # Probability of jump event
    jump_size_mean ~ dnorm(PRIOR_JUMP_MEAN, PRIOR_JUMP_PREC)           # Mean jump size
    jump_size_sd ~ dgamma(PRIOR_JUMP_SD_SHAPE, PRIOR_JUMP_SD_RATE)     # SD of jump size
    root_state ~ dnorm(PRIOR_ROOT_MEAN, PRIOR_ROOT_PREC)               # Root state
    
    # Assign root state
    states[N_tips + 1] <- root_state
    
    # Define jump parameters for reporting
    jump_size <- jump_size_mean
    
    # Loop through edges to define tree structure
    for(i in 1:N_edges) {
      # Determine if edge has a jump (chromosome rearrangement)
      jump[i] ~ dbern(jump_prob)
      
      # Generate jump contribution
      jump_contrib[i] ~ dnorm(jump_size_mean, 1/jump_size_sd^2)
      
      # Expected value combines BM and potential jump
      expected[i] <- states[parent[i]] + jump[i] * jump_contrib[i]
      
      # Define state changes via BM process with possible jumps
      states[child[i]] ~ dnorm(expected[i], 1/(sigma2 * edge_length[i]))
    }
    
    # Data likelihood for tips
    for(i in 1:N_tips) {
      y[i] ~ dnorm(states[i], 1/0.01)  # Small observational variance
    }
  }
  "
  
  # Replace prior placeholders with actual values
  model_code <- gsub("PRIOR_SIGMA_SHAPE", priors$sigma_shape, model_code)
  model_code <- gsub("PRIOR_SIGMA_RATE", priors$sigma_rate, model_code)
  model_code <- gsub("PRIOR_JUMP_PROB_ALPHA", priors$jump_prob_alpha, model_code)
  model_code <- gsub("PRIOR_JUMP_PROB_BETA", priors$jump_prob_beta, model_code)
  model_code <- gsub("PRIOR_JUMP_MEAN", priors$jump_mean, model_code)
  model_code <- gsub("PRIOR_JUMP_PREC", 1/priors$jump_var, model_code)
  model_code <- gsub("PRIOR_JUMP_SD_SHAPE", priors$jump_sd_shape, model_code)
  model_code <- gsub("PRIOR_JUMP_SD_RATE", priors$jump_sd_rate, model_code)
  model_code <- gsub("PRIOR_ROOT_MEAN", priors$root_mean, model_code)
  model_code <- gsub("PRIOR_ROOT_PREC", 1/priors$root_var, model_code)
  
  return(model_code)
}

#' Build variable rate model for JAGS
#' 
#' @param priors Prior distributions
#' @return JAGS model code as string
#' @keywords internal
build_rateshift_model <- function(priors) {
  model_code <- "
  model {
    # Priors for model parameters
    sigma2 ~ dgamma(PRIOR_SIGMA_SHAPE, PRIOR_SIGMA_RATE)    # Base evolutionary rate
    root_state ~ dnorm(PRIOR_ROOT_MEAN, PRIOR_ROOT_PREC)    # Root state
    
    # Rate multipliers for different edge groups
    for(i in 1:N_groups) {
      rate_multiplier[i] ~ dgamma(PRIOR_RATE_MULT_SHAPE, PRIOR_RATE_MULT_RATE)
    }
    
    # Assign root state
    states[N_tips + 1] <- root_state
    
    # Loop through edges to define tree structure
    for(i in 1:N_edges) {
      # Calculate edge-specific rate
      edge_rate[i] <- sigma2 * rate_multiplier[edge_group[i]]
      
      # Define state changes via BM with variable rates
      states[child[i]] ~ dnorm(states[parent[i]], 1/(edge_rate[i] * edge_length[i]))
    }
    
    # Data likelihood for tips
    for(i in 1:N_tips) {
      y[i] ~ dnorm(states[i], 1/0.01)  # Small observational variance
    }
  }
  "
  
  # Replace prior placeholders with actual values
  model_code <- gsub("PRIOR_SIGMA_SHAPE", priors$sigma_shape, model_code)
  model_code <- gsub("PRIOR_SIGMA_RATE", priors$sigma_rate, model_code)
  model_code <- gsub("PRIOR_RATE_MULT_SHAPE", priors$rate_mult_shape, model_code)
  model_code <- gsub("PRIOR_RATE_MULT_RATE", priors$rate_mult_rate, model_code)
  model_code <- gsub("PRIOR_ROOT_MEAN", priors$root_mean, model_code)
  model_code <- gsub("PRIOR_ROOT_PREC", 1/priors$root_var, model_code)
  
  return(model_code)
}

#' Get default priors for Bayesian models
#' 
#' @param model Model type
#' @return List of prior values
#' @keywords internal
get_default_priors <- function(model) {
  # Common priors
  base_priors <- list(
    sigma_shape = 0.5,   # Shape for sigma2 gamma prior
    sigma_rate = 0.05,   # Rate for sigma2 gamma prior
    root_mean = 12,      # Mean for root state normal prior
    root_var = 100       # Variance for root state normal prior
  )
  
  # Model-specific priors
  if(model == "BM") {
    return(base_priors)
  } else if(model == "OU") {
    ou_priors <- c(base_priors, list(
      alpha_shape = 1,    # Shape for alpha gamma prior
      alpha_rate = 1,     # Rate for alpha gamma prior
      theta_mean = 12,    # Mean for theta normal prior
      theta_var = 100     # Variance for theta normal prior
    ))
    return(ou_priors)
  } else if(model == "JumpBM") {
    jump_priors <- c(base_priors, list(
      jump_prob_alpha = 1,   # Alpha for jump probability beta prior
      jump_prob_beta = 9,    # Beta for jump probability beta prior (favors rare jumps)
      jump_mean = 0,         # Mean for jump size normal prior
      jump_var = 25,         # Variance for jump size normal prior
      jump_sd_shape = 1,     # Shape for jump size SD gamma prior
      jump_sd_rate = 0.2     # Rate for jump size SD gamma prior
    ))
    return(jump_priors)
  } else if(model == "RateShift") {
    rateshift_priors <- c(base_priors, list(
      rate_mult_shape = 1,   # Shape for rate multiplier gamma prior
      rate_mult_rate = 1     # Rate for rate multiplier gamma prior
    ))
    return(rateshift_priors)
  } else {
    stop(paste("Unsupported model for priors:", model))
  }
}

#' Process MCMC samples to extract ancestral state estimates
#' 
#' @param mcmc_samples MCMC sample object
#' @param tree Phylogenetic tree
#' @param tip_states Observed tip states
#' @param model Model type
#' @return Processed results including ancestral states and model parameters
#' @keywords internal
process_mcmc_samples <- function(mcmc_samples, tree, tip_states, model) {
  # Convert to coda mcmc object if not already
  if(!inherits(mcmc_samples, "mcmc") && !inherits(mcmc_samples, "mcmc.list")) {
    mcmc_samples <- coda::as.mcmc(mcmc_samples)
  }
  
  # Check convergence
  convergence <- list(
    gelman_diag = NA,
    geweke_diag = NA,
    effective_size = NA
  )
  
  # Get Gelman-Rubin diagnostic if multiple chains
  if(inherits(mcmc_samples, "mcmc.list") && length(mcmc_samples) > 1) {
    convergence$gelman_diag <- coda::gelman.diag(mcmc_samples)
  }
  
  # Geweke diagnostic (for first chain if multiple)
  if(inherits(mcmc_samples, "mcmc.list")) {
    chain1 <- mcmc_samples[[1]]
  } else {
    chain1 <- mcmc_samples
  }
  
  # Get sample columns for ancestral states only
  n_tips <- length(tip_states)
  n_nodes <- tree$Nnode
  node_ids <- (n_tips + 1):(n_tips + n_nodes)
  
  # Extract ancestral state samples
  anc_col_pattern <- "^states\\["
  anc_cols <- grep(anc_col_pattern, colnames(chain1))
  
  # Create empty matrix for ancestral states
  anc_states_matrix <- matrix(NA, nrow = nrow(chain1), ncol = length(node_ids))
  colnames(anc_states_matrix) <- paste0("node_", node_ids)
  
  # Extract node samples
  for(i in 1:length(node_ids)) {
    node_id <- node_ids[i]
    col_idx <- which(grepl(paste0("states\\[", node_id, "\\]"), colnames(chain1)))
    if(length(col_idx) > 0) {
      anc_states_matrix[, i] <- chain1[, col_idx]
    }
  }
  
  # Create ancestral states data frame
  ancestral_states <- data.frame(
    node_id = node_ids,
    state = colMeans(anc_states_matrix),
    ci_lower = apply(anc_states_matrix, 2, quantile, probs = 0.025),
    ci_upper = apply(anc_states_matrix, 2, quantile, probs = 0.975),
    stringsAsFactors = FALSE
  )
  
  # Extract model parameters
  model_params <- list()
  
  # Common parameters
  sigma2_cols <- grep("^sigma2", colnames(chain1))
  if(length(sigma2_cols) > 0) {
    model_params$sigma2 <- mean(chain1[, sigma2_cols])
    model_params$sigma2_ci <- quantile(chain1[, sigma2_cols], probs = c(0.025, 0.975))
  }
  
  # Model-specific parameters
  if(model == "OU") {
    alpha_cols <- grep("^alpha", colnames(chain1))
    if(length(alpha_cols) > 0) {
      model_params$alpha <- mean(chain1[, alpha_cols])
      model_params$alpha_ci <- quantile(chain1[, alpha_cols], probs = c(0.025, 0.975))
    }
    
    theta_cols <- grep("^theta", colnames(chain1))
    if(length(theta_cols) > 0) {
      model_params$theta <- mean(chain1[, theta_cols])
      model_params$theta_ci <- quantile(chain1[, theta_cols], probs = c(0.025, 0.975))
    }
  } else if(model == "JumpBM") {
    jump_prob_cols <- grep("^jump_prob", colnames(chain1))
    if(length(jump_prob_cols) > 0) {
      model_params$jump_prob <- mean(chain1[, jump_prob_cols])
      model_params$jump_prob_ci <- quantile(chain1[, jump_prob_cols], probs = c(0.025, 0.975))
    }
    
    jump_size_cols <- grep("^jump_size", colnames(chain1))
    if(length(jump_size_cols) > 0) {
      model_params$jump_size <- mean(chain1[, jump_size_cols])
      model_params$jump_size_ci <- quantile(chain1[, jump_size_cols], probs = c(0.025, 0.975))
    }
  } else if(model == "RateShift") {
    rate_mult_cols <- grep("^rate_multiplier", colnames(chain1))
    if(length(rate_mult_cols) > 0) {
      model_params$rate_multiplier <- colMeans(chain1[, rate_mult_cols])
      model_params$rate_multiplier_ci <- apply(chain1[, rate_mult_cols], 2, 
                                             quantile, probs = c(0.025, 0.975))
    }
  }
  
  # Calculate effective sample size
  convergence$effective_size <- coda::effectiveSize(chain1[, sigma2_cols])
  
  # Geweke diagnostic for sigma2
  try({
    convergence$geweke_diag <- coda::geweke.diag(chain1[, sigma2_cols])
  }, silent = TRUE)
  
  return(list(
    ancestral_states = ancestral_states,
    model_parameters = model_params,
    convergence = convergence
  ))
}

#' Calculate Bayesian reconstruction quality metrics
#' 
#' @param bayesian_result Bayesian reconstruction result
#' @return Quality assessment metrics
#' @keywords internal
calculate_bayesian_quality <- function(bayesian_result) {
  # Extract parameters
  tree <- bayesian_result$tree
  ancestors <- bayesian_result$ancestral_states
  mcmc_samples <- bayesian_result$mcmc_samples
  convergence <- bayesian_result$convergence
  
  # Calculate confidence interval widths
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
  
  # Calculate MCMC diagnostics if available
  gelman_stat <- NA
  if(!is.null(convergence$gelman_diag) && !is.null(convergence$gelman_diag$psrf)) {
    gelman_stat <- max(convergence$gelman_diag$psrf[, 1])
  }
  
  geweke_stat <- NA
  if(!is.null(convergence$geweke_diag) && !is.null(convergence$geweke_diag$z)) {
    geweke_stat <- convergence$geweke_diag$z
  }
  
  # Effective sample size
  ess <- NA
  if(!is.null(convergence$effective_size)) {
    ess <- convergence$effective_size
  }
  
  # Return quality metrics
  quality <- list(
    mean_ci_width = mean_ci_width,
    max_ci_width = max(ci_widths),
    min_ci_width = min(ci_widths),
    root_uncertainty = root_uncertainty,
    gelman_stat = gelman_stat,
    geweke_stat = geweke_stat,
    effective_sample_size = ess,
    mcmc_quality = ifelse(is.na(gelman_stat) || gelman_stat < 1.1, "good", "poor")
  )
  
  return(quality)
}
