#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Simulation Tools
# Author: Bioinformatics Team
# Date: 2025-03-24 
# Description: Provides unified functions for simulating chromosome evolution 
#              under various models.
#===============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(ape)
  library(geiger)
  library(phytools) # For discrete character evolution if needed later
})

#===============================================================================
# Internal Helper Functions
#===============================================================================

#' Transform phylogenetic tree according to model parameters
#' 
#' @param tree Phylogenetic tree
#' @param model Model type: "EB", "ACDC", "lambda"
#' @param parameter Model parameter value
#' @return Transformed tree
#' @keywords internal
transform_tree_by_model <- function(tree, model, parameter) {
  # Create a copy of the tree to modify
  transformed_tree <- tree
  
  # Apply appropriate transformation
  if(model == "EB") {
    # Early Burst: exponential rate change over time
    # Branch lengths are transformed as: b' = b * (exp(r*T) - exp(r*(T-b)))/r
    # where T is the tree height, b is the branch length, r is the rate parameter
    
    tree_height <- max(node.depth.edgelength(tree))
    node_depths <- node.depth.edgelength(tree)
    
    for(i in 1:nrow(tree$edge)) {
      parent <- tree$edge[i, 1]
      child <- tree$edge[i, 2]
      
      # Calculate depth at start and end of branch
      if(child <= Ntip(tree)) {
        end_depth <- node_depths[child]
      } else {
        end_depth <- node_depths[child - Ntip(tree) + Ntip(tree)]
      }
      
      if(parent <= Ntip(tree)) {
        start_depth <- node_depths[parent]
      } else {
        start_depth <- node_depths[parent - Ntip(tree) + Ntip(tree)]
      }
      
      branch_length <- tree$edge.length[i]
      
      # Apply EB transformation
      if(abs(parameter) < 1e-6) {
        # For very small r, use linear approximation
        transformed_length <- branch_length
      } else {
        # Full transformation
        transformed_length <- (exp(parameter * start_depth) - 
                             exp(parameter * end_depth)) / parameter
      }
      
      transformed_tree$edge.length[i] <- transformed_length
    }
    
  } else if(model == "ACDC") {
    # ACDC: accelerating/decelerating rates
    # Branch lengths are transformed as: b' = b * exp(beta * d)
    # where d is the node depth, beta is the parameter
    
    node_depths <- node.depth.edgelength(tree)
    
    for(i in 1:nrow(tree$edge)) {
      parent <- tree$edge[i, 1]
      child <- tree$edge[i, 2]
      
      # Calculate depth at start of branch
      if(parent <= Ntip(tree)) {
        depth <- node_depths[parent]
      } else {
        depth <- node_depths[parent - Ntip(tree) + Ntip(tree)]
      }
      
      # Apply ACDC transformation
      transformed_tree$edge.length[i] <- tree$edge.length[i] * exp(parameter * depth)
    }
    
  } else if(model == "lambda") {
    # Lambda: scales internal branches
    # For internal branches: b' = b * lambda
    # For terminal branches: unchanged
    
    for(i in 1:nrow(tree$edge)) {
      child <- tree$edge[i, 2]
      
      # Check if branch leads to a tip
      if(child > Ntip(tree)) {
        # Internal branch
        transformed_tree$edge.length[i] <- tree$edge.length[i] * parameter
      }
    }
  }
  
  return(transformed_tree)
}

#' Simulate chromosome evolution under a specified continuous model
#' (Internal helper function)
#' @param tree Phylogenetic tree
#' @param model Evolutionary model name (e.g., "BM", "OU", "EB", "ACDC", "lambda")
#' @param root_state Starting state at root
#' @param sigma_sq Rate parameter (variance)
#' @param extra_params A list of additional model-specific parameters 
#'                     (e.g., alpha for OU, r for EB)
#' @return Simulated chromosome counts at tips and nodes
#' @keywords internal
simulate_continuous_model_internal <- function(tree, model, root_state, sigma_sq, extra_params = list()) {
  # Default to Brownian Motion
  if(model == "BM") {
    # Use standard BM simulation
    sim_all <- geiger::sim.char(tree, par = sigma_sq, model = "BM", root = root_state, nsim = 1)
    sim_tips <- sim_all[,,1][tree$tip.label]
    sim_nodes <- sim_all[,,1][(Ntip(tree)+1):nrow(sim_all)] # Geiger sim.char output needs careful indexing for nodes
  } else if(model == "OU") {
    # Extract OU parameters
    alpha <- extra_params$alpha
    theta <- extra_params$theta # Note: geiger::sim.char uses 'alpha' for strength of selection, 'theta' (optional) not directly used by sim.char's OU.
                               # Assuming theta is the ancestral optimum if provided, else root_state is used as such.
    
    # OU simulation
    sim_all <- geiger::sim.char(tree, par = list(sigsq = sigma_sq, alpha = alpha, theta = ifelse(is.null(theta), root_state, theta)), model = "OU", root = root_state, nsim = 1)
    sim_tips <- sim_all[,,1][tree$tip.label]
    sim_nodes <- sim_all[,,1][(Ntip(tree)+1):nrow(sim_all)]
  } else if(model %in% c("EB", "ACDC", "lambda")) {
    # Extract model-specific parameter
    rate_param <- NULL
    if(model == "EB") {
      rate_param <- extra_params$r
    } else if(model == "ACDC") {
      rate_param <- extra_params$beta
    } else if(model == "lambda") {
      rate_param <- extra_params$lambda
    }
    
    if(is.null(rate_param) && model != "lambda"){ # lambda can be 1, meaning no transformation
        stop(paste("Parameter for model", model, "not provided in extra_params"))
    }
    if(model == "lambda" && is.null(rate_param)) rate_param <- 1 # Default lambda to 1 (BM)

    # Transform tree according to the model
    transformed_tree <- transform_tree_by_model(tree, model, rate_param)
    
    # Simulate BM on transformed tree
    sim_all <- geiger::sim.char(transformed_tree, par = sigma_sq, model = "BM", root = root_state, nsim = 1)
    sim_tips <- sim_all[,,1][tree$tip.label]
    # Node states from sim.char on a transformed tree are relative to the transformed tree.
    # For simplicity, returning node states based on the transformed tree.
    # A more rigorous approach might map these back or re-estimate on original tree.
    sim_nodes <- sim_all[,,1][(Ntip(transformed_tree)+1):nrow(sim_all)]
    names(sim_nodes) <- (Ntip(tree)+1):(Ntip(tree)+Nnode(tree)) # Ensure correct node naming

  } else {
    stop(paste("Unsupported continuous model for simulation:", model))
  }
  
  # Ensure non-negative counts (impose floor at 1 if integer-like data is expected)
  # For now, returning raw continuous values as this is a general continuous trait simulator
  # sim_tips <- pmax(1, round(sim_tips))
  # sim_nodes <- pmax(1, round(sim_nodes))

  return(list(tip_states = sim_tips, ancestral_states = sim_nodes))
}

#' Simulate chromosome evolution with discrete events (Internal helper function)
#'
#' Simulate chromosome count evolution on a tree using discrete event models
#' like "BM" (simple Brownian on discrete states), "jumps", or "hybrid".
#'
#' @param tree Phylogenetic tree object.
#' @param root_state Root chromosome count.
#' @param model Evolutionary model: "BM", "jumps", or "hybrid".
#' @param params Simulation parameters list, specific to the model:
#'   - For "BM": `sigma` (standard deviation of change per unit branch length).
#'   - For "jumps": `sigma` (background BM rate), `jump_prob` (probability of a jump per unit branch length), 
#'                  `jump_size` (mean size of a jump, exponential distribution).
#'   - For "hybrid": `sigma`, `jump_prob`, `jump_size`, `trend` (drift parameter).
#' @return A list containing `tip_states` (named vector of simulated counts for tips)
#'         and `ancestral_states` (vector of simulated counts for internal nodes).
#' @keywords internal
simulate_discrete_event_model_internal <- function(tree, root_state, model, params) {
  # Initialize result list (though not strictly used here, kept for consistency with prior structure)
  # sim_data <- list(
  #   tree = tree,
  #   model = model,
  #   params = params,
  #   root_state = root_state
  # )

  # Prepare node containers
  n_total_nodes <- Ntip(tree) + Nnode(tree) # Total number of nodes in the tree
  node_states <- numeric(n_total_nodes)
  
  # Assign root state. The root is usually Ntip(tree) + 1.
  # Need to map this to the correct index in node_states if tree$edge node numbers are not sequential from 1 to n_total_nodes.
  # Assuming ape tree structure where tips are 1:Ntip(tree) and nodes are (Ntip(tree)+1):(Ntip(tree)+Nnode(tree))
  root_node_index <- Ntip(tree) + 1 
  node_states[root_node_index] <- root_state

  # Ensure edges are processed in an order that parents are typically calculated before children
  # (e.g. postorder traversal for assigning states from root to tips).
  # ape's tree$edge is usually ordered correctly or can be reordered.
  # Here, we iterate through edges assuming parent state is already set when we reach an edge.
  
  # For simulation from root to tips, a pre-order (root-first) traversal logic is better.
  # However, the provided code iterates edges and assumes parent state is known.
  # This works if node_states is globally updated.
  
  # Get edges in a way that ensures parent is simulated before child
  # This is implicitly handled by iterating tree$edge if root is set first.
  
  # Apply simulation model
  if(model == "BM") {
    # Brownian motion simulation (discretized)
    sigma <- params$sigma
    
    for(i in 1:nrow(tree$edge)) {
      parent_node_id <- tree$edge[i, 1]
      child_node_id <- tree$edge[i, 2]
      branch_length <- tree$edge.length[i]
      
      parent_state <- node_states[parent_node_id]
      
      change <- rnorm(1, mean = 0, sd = sigma * sqrt(branch_length))
      child_state <- parent_state + change
      
      node_states[child_node_id] <- max(1, round(child_state)) # Round for discrete-like behavior
    }
  } else if(model == "jumps") {
    sigma <- params$sigma
    jump_prob_per_unit_branch <- params$jump_prob 
    jump_size_mean <- params$jump_size # This was rexp rate = 1/jump_size, so jump_size is mean
    
    for(i in 1:nrow(tree$edge)) {
      parent_node_id <- tree$edge[i, 1]
      child_node_id <- tree$edge[i, 2]
      branch_length <- tree$edge.length[i]
      parent_state <- node_states[parent_node_id]
      
      # Number of jumps on a branch could be Poisson distributed if jump_prob is a rate
      # Original code implies jump_prob is per unit time, so for a branch, prob is jump_prob * branch_length
      # For simplicity, let's assume at most one jump event per edge for now, or interpret jump_prob differently.
      # The original code had if(runif(1) < jump_prob) which is not scaled by branch length.
      # Let's use a rate: number of jumps ~ Poisson(jump_prob_per_unit_branch * branch_length)
      
      num_jumps <- rpois(1, jump_prob_per_unit_branch * branch_length)
      total_jump_change <- 0
      if (num_jumps > 0) {
        for (j in 1:num_jumps) {
          if(runif(1) < 0.5) { # Fusion
            total_jump_change <- total_jump_change - rexp(1, rate = 1/jump_size_mean)
          } else { # Fission/duplication
            total_jump_change <- total_jump_change + rexp(1, rate = 1/jump_size_mean)
          }
        }
      }
      
      # Background BM change
      bm_change <- rnorm(1, mean = 0, sd = sigma * sqrt(branch_length))
      child_state <- parent_state + bm_change + total_jump_change
      node_states[child_node_id] <- max(1, round(child_state))
    }
  } else if(model == "hybrid") {
    sigma <- params$sigma
    jump_prob_per_unit_branch <- params$jump_prob
    jump_size_mean <- params$jump_size
    trend <- if(is.null(params$trend)) 0 else params$trend
    
    for(i in 1:nrow(tree$edge)) {
      parent_node_id <- tree$edge[i, 1]
      child_node_id <- tree$edge[i, 2]
      branch_length <- tree$edge.length[i]
      parent_state <- node_states[parent_node_id]
      
      base_change <- rnorm(1, mean = trend * branch_length, sd = sigma * sqrt(branch_length))
      
      num_jumps <- rpois(1, jump_prob_per_unit_branch * branch_length)
      total_jump_change <- 0
      if (num_jumps > 0) {
         for (j in 1:num_jumps) {
            if(runif(1) < 0.5) { # Fusion
              total_jump_change <- total_jump_change - rexp(1, rate = 1/jump_size_mean)
            } else { # Fission/duplication
              total_jump_change <- total_jump_change + rexp(1, rate = 1/jump_size_mean)
            }
        }
      }
      child_state <- parent_state + base_change + total_jump_change
      node_states[child_node_id] <- max(1, round(child_state))
    }
  } else {
    stop("Unsupported discrete event simulation model: ", model)
  }
  
  tip_states <- node_states[1:Ntip(tree)]
  names(tip_states) <- tree$tip.label
  # Internal nodes are (Ntip(tree)+1) to (Ntip(tree)+Nnode(tree))
  ancestral_states <- node_states[(Ntip(tree)+1):n_total_nodes]
  # Ensure ancestral_states are named if downstream code expects it, e.g. by node ID
  # names(ancestral_states) <- (Ntip(tree)+1):n_total_nodes
  
  return(list(tip_states = tip_states, ancestral_states = ancestral_states))
}

#===============================================================================
# Main Exported Simulation Function
#===============================================================================

#' Unified Chromosome Evolution Simulation Wrapper
#'
#' Runs a single simulation of chromosome evolution based on the specified type and model.
#'
#' @param tree Phylogenetic tree object.
#' @param simulation_type Type of simulation: "continuous_fitted" for models like BM, OU, EB, etc.,
#'                        or "discrete_event" for models like parsimony's BM, jumps, hybrid.
#' @param model_name Name of the evolutionary model (e.g., "BM", "OU", "EB", "ACDC", "lambda", "jumps", "hybrid").
#' @param root_state The ancestral state at the root of the tree.
#' @param params A list of parameters required for the chosen model.
#'   - For "continuous_fitted" models (delegated to `simulate_continuous_model_internal`):
#'     - `sigma_sq`: Rate parameter (variance or sigsq).
#'     - `extra_params`: A list of additional model-specific parameters (e.g., `alpha` for OU, `theta` for OU, `r` for EB, `beta` for ACDC, `lambda` for Lambda).
#'   - For "discrete_event" models (delegated to `simulate_discrete_event_model_internal`):
#'     - These are passed directly as `params` to the internal function. Examples:
#'       - For "BM": `sigma`.
#'       - For "jumps": `sigma`, `jump_prob`, `jump_size`.
#'       - For "hybrid": `sigma`, `jump_prob`, `jump_size`, `trend`.
#' @return A list containing:
#'   - `tip_states`: Named vector of simulated chromosome counts for tip species.
#'   - `ancestral_states`: Vector of simulated chromosome counts for internal nodes.
#'         (Order corresponds to internal nodes of the tree, typically (Ntip+1) to (Ntip+Nnode)).
#' @export
#' @examples
#' \dontrun{
#'   # Example for continuous model
#'   # library(ape)
#'   # tree <- rphylo(5, birth=1, death=0)
#'   # params_cont <- list(sigma_sq = 0.1, extra_params = list(alpha=0.5, theta=10))
#'   # sim_data_cont <- run_chromosome_simulation(tree, "continuous_fitted", "OU", 
#'   #                                          root_state = 10, params = params_cont)
#'   # print(sim_data_cont$tip_states)
#'
#'   # Example for discrete event model
#'   # params_disc <- list(sigma = 0.5, jump_prob = 0.1, jump_size = 2)
#'   # sim_data_disc <- run_chromosome_simulation(tree, "discrete_event", "jumps",
#'   #                                           root_state = 10, params = params_disc)
#'   # print(sim_data_disc$tip_states)
#' }
run_chromosome_simulation <- function(tree, simulation_type, model_name, root_state, params) {
  
  if (!inherits(tree, "phylo")) {
    stop("tree must be a phylo object.")
  }
  if (!is.character(simulation_type) || length(simulation_type) != 1) {
    stop("simulation_type must be a single string.")
  }
  if (!is.character(model_name) || length(model_name) != 1) {
    stop("model_name must be a single string.")
  }
  if (!is.numeric(root_state) || length(root_state) != 1) {
    stop("root_state must be a single numeric value.")
  }
  if (!is.list(params)) {
    stop("params must be a list.")
  }

  message(sprintf("Running simulation: Type=%s, Model=%s, Root=%s", 
                  simulation_type, model_name, root_state))

  if (simulation_type == "continuous_fitted") {
    # Check for necessary params for continuous models
    if (is.null(params$sigma_sq)) {
      stop("For continuous_fitted simulation, params must contain 'sigma_sq'.")
    }
    # extra_params is optional, defaults to empty list if not provided by simulate_continuous_model_internal
    extra_p <- if (is.null(params$extra_params)) list() else params$extra_params
    
    sim_results <- simulate_continuous_model_internal(
      tree = tree,
      model = model_name,
      root_state = root_state,
      sigma_sq = params$sigma_sq,
      extra_params = extra_p
    )
  } else if (simulation_type == "discrete_event") {
    # For discrete_event, params are passed directly
    sim_results <- simulate_discrete_event_model_internal(
      tree = tree,
      root_state = root_state,
      model = model_name, # e.g., "BM", "jumps", "hybrid" from parsimony context
      params = params      # e.g., list(sigma=1.0, jump_prob=0.05, jump_size=3)
    )
  } else {
    stop(paste("Unsupported simulation_type:", simulation_type, 
               ". Choose 'continuous_fitted' or 'discrete_event'."))
  }
  
  # Ensure output structure is consistent: list with tip_states and ancestral_states
  if (!is.list(sim_results) || !all(c("tip_states", "ancestral_states") %in% names(sim_results))) {
    stop("Internal simulation function did not return the expected structure (list with tip_states and ancestral_states).")
  }
  
  return(sim_results)
}
