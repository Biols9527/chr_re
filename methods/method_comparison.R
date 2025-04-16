#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Method Comparison Module
# Author: Bioinformatics Team
# Date: 2025-04-15
# Description: Implements tools for comparing different ancestral chromosome
#              reconstruction methods and evaluating their performance
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggplot2)
  library(reshape2)
  library(parallel)
})

#===============================================================================
# Method Comparison Core Functions
#===============================================================================

#' Compare multiple ancestral reconstruction methods
#' 
#' Systematically compares different methods for ancestral chromosome reconstruction
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts, named vector with species names
#' @param methods Vector of method names to compare: "parsimony", "ML", "bayesian", "ensemble"
#' @param method_params List of parameters for each method
#' @param eval_approach Evaluation approach: "cv" (cross-validation), "simulation", or "both"
#' @param selection_criterion Criterion for ranking methods: "rmse", "mae", "correlation"
#' @param cv_folds Number of cross-validation folds
#' @param create_plots Whether to create comparison plots
#' @param n_cores Number of cores for parallel processing (NULL = auto-detect)
#' @return List with method comparison results
#' @export
compare_reconstruction_methods <- function(tree, chr_counts, 
                                         methods = c("parsimony", "ML", "bayesian"),
                                         method_params = NULL,
                                         eval_approach = "cv", 
                                         selection_criterion = "rmse",
                                         cv_folds = 5,
                                         create_plots = TRUE,
                                         n_cores = NULL) {
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Ensure chr_counts is a named vector
  if(is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector with species names")
  }
  
  # Check if methods are supported
  supported_methods <- c("parsimony", "ML", "bayesian", "ensemble")
  unsupported <- setdiff(methods, supported_methods)
  if(length(unsupported) > 0) {
    stop(paste("Unsupported methods:", paste(unsupported, collapse = ", ")))
  }
  
  # Validate selection criterion
  valid_criteria <- c("rmse", "mae", "correlation")
  if(!selection_criterion %in% valid_criteria) {
    stop(paste0("Invalid selection criterion. Choose one of: ", paste(valid_criteria, collapse = ", ")))
  }
  
  # Ensure tree and data contain same species
  common_species <- intersect(names(chr_counts), tree$tip.label)
  if(length(common_species) < 5) {
    stop("At least 5 species in common required for method comparison")
  }
  
  # Prune tree and data to match
  pruned_tree <- ape::keep.tip(tree, common_species)
  pruned_counts <- chr_counts[common_species]
  
  # Setup default method parameters if not provided
  if(is.null(method_params)) {
    method_params <- list(
      parsimony = list(method = "wagner"),
      ML = list(model = "ER"),
      bayesian = list(model = "BM", mcmc_settings = list(
        iterations = 20000, burnin = 5000, thinning = 10, chains = 2
      )),
      ensemble = list(weights = "equal")
    )
  } else {
    # Fill in missing parameters with defaults
    default_params <- list(
      parsimony = list(method = "wagner"),
      ML = list(model = "ER"),
      bayesian = list(model = "BM", mcmc_settings = list(
        iterations = 20000, burnin = 5000, thinning = 10, chains = 2
      )),
      ensemble = list(weights = "equal")
    )
    
    for(method in methods) {
      if(is.null(method_params[[method]])) {
        method_params[[method]] <- default_params[[method]]
      }
    }
  }
  
  # Initialize results structure
  results <- list(
    tree = pruned_tree,
    data = pruned_counts,
    methods = methods,
    method_params = method_params,
    eval_approach = eval_approach,
    selection_criterion = selection_criterion,
    cv_results = NULL,
    simulation_results = NULL,
    overall_ranking = NULL,
    best_method = NULL,
    plots = list()
  )
  
  # Set up parallel processing
  parallel_setup <- setup_parallel(n_cores)
  use_parallel <- parallel_setup$use_parallel
  cl <- parallel_setup$cl
  
  # Run evaluations based on selected approach
  if(eval_approach %in% c("cv", "both")) {
    message("Performing cross-validation evaluation of reconstruction methods...")
    cv_results <- evaluate_methods_cv(pruned_tree, pruned_counts, methods, method_params, 
                                    cv_folds, use_parallel, cl)
    results$cv_results <- cv_results
  }
  
  if(eval_approach %in% c("simulation", "both")) {
    message("Performing simulation-based evaluation of reconstruction methods...")
    sim_results <- evaluate_methods_simulation(pruned_tree, pruned_counts, methods, method_params, 
                                             use_parallel, cl)
    results$simulation_results <- sim_results
  }
  
  # Clean up parallel cluster if used
  if(use_parallel && !is.null(cl)) {
    parallel::stopCluster(cl)
  }
  
  # Determine overall ranking and best method
  results$overall_ranking <- rank_methods(results, selection_criterion)
  results$best_method <- results$overall_ranking$Method[1]
  
  # Create comparison plots if requested
  if(create_plots) {
    results$plots <- create_comparison_plots(results)
  }
  
  message(paste("Method comparison complete. Best method based on", 
                selection_criterion, "is:", results$best_method))
  
  return(results)
}

#' Set up parallel processing environment
#' 
#' @param n_cores Number of cores to use (NULL = auto-detect)
#' @return List with parallel settings
#' @keywords internal
setup_parallel <- function(n_cores) {
  use_parallel <- FALSE
  cl <- NULL
  
  # Check if parallel processing is possible
  if(requireNamespace("parallel", quietly = TRUE)) {
    # Determine number of cores to use
    if(is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
      if(n_cores < 1) n_cores <- 1
    }
    
    # Only use parallel if we have more than 1 core
    if(n_cores > 1) {
      use_parallel <- TRUE
      cl <- parallel::makeCluster(n_cores)
      
      # Load required packages on each node
      parallel::clusterEvalQ(cl, {
        library(ape)
        library(phytools)
      })
      
      message(sprintf("Using %d cores for parallel processing", n_cores))
    }
  }
  
  return(list(use_parallel = use_parallel, cl = cl, n_cores = n_cores))
}

#' Evaluate methods using cross-validation
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param methods Methods to evaluate
#' @param method_params Method parameters
#' @param cv_folds Number of cross-validation folds
#' @param use_parallel Whether to use parallel processing
#' @param cl Parallel cluster
#' @return Cross-validation results
#' @keywords internal
evaluate_methods_cv <- function(tree, chr_counts, methods, method_params, 
                              cv_folds = 5, use_parallel = FALSE, cl = NULL) {
  # Ensure methods and method_params match
  methods_to_evaluate <- intersect(methods, names(method_params))
  
  if(length(methods_to_evaluate) == 0) {
    stop("No valid methods to evaluate")
  }
  
  message(sprintf("Evaluating %d methods using %d-fold cross-validation...", 
                 length(methods_to_evaluate), cv_folds))
  
  # Setup result containers
  cv_results <- list(
    method_performance = data.frame(
      Method = methods_to_evaluate,
      MAE = NA_real_,
      RMSE = NA_real_,
      Correlation = NA_real_,
      Success_Rate = NA_real_,
      stringsAsFactors = FALSE
    ),
    detailed_results = list(),
    folds = list()
  )
  
  # Create CV folds
  species <- names(chr_counts)
  n_species <- length(species)
  fold_size <- ceiling(n_species / cv_folds)
  
  # Shuffle species
  set.seed(42)  # For reproducibility
  shuffled_species <- sample(species)
  
  # Assign to folds
  folds <- list()
  for(i in 1:cv_folds) {
    start_idx <- (i - 1) * fold_size + 1
    end_idx <- min(i * fold_size, n_species)
    
    if(start_idx <= n_species) {
      test_species <- shuffled_species[start_idx:end_idx]
      train_species <- setdiff(species, test_species)
      
      folds[[i]] <- list(
        train = train_species,
        test = test_species
      )
    }
  }
  
  cv_results$folds <- folds
  
  # Run CV for each method
  for(method in methods_to_evaluate) {
    message(sprintf("  Cross-validating method: %s", method))
    
    # Get method parameters
    params <- method_params[[method]]
    
    # Setup result container for this method
    method_results <- list(
      predictions = data.frame(
        Species = character(0),
        Actual = numeric(0),
        Predicted = numeric(0),
        Fold = integer(0),
        stringsAsFactors = FALSE
      ),
      errors = data.frame(
        Species = character(0),
        Error = numeric(0),
        AbsError = numeric(0),
        SquaredError = numeric(0),
        Fold = integer(0),
        stringsAsFactors = FALSE
      ),
      fold_metrics = data.frame(
        Fold = integer(0),
        MAE = numeric(0),
        RMSE = numeric(0),
        Correlation = numeric(0),
        Success = logical(0),
        stringsAsFactors = FALSE
      )
    )
    
    # Process each fold
    for(fold_idx in 1:length(folds)) {
      fold <- folds[[fold_idx]]
      
      # Create training tree and data
      train_tree <- ape::keep.tip(tree, fold$train)
      train_counts <- chr_counts[fold$train]
      
      # Run reconstruction on training data
      reconstruction_result <- try({
        run_reconstruction_method(train_tree, train_counts, method, params)
      }, silent = TRUE)
      
      # Check if reconstruction was successful
      fold_success <- !inherits(reconstruction_result, "try-error") && !is.null(reconstruction_result)
      
      if(fold_success) {
        # Make predictions for test species
        predictions <- predict_test_species(tree, reconstruction_result, fold$test)
        
        # Calculate errors
        errors <- predictions$actual - predictions$predicted
        abs_errors <- abs(errors)
        squared_errors <- errors^2
        
        # Store predictions
        pred_df <- data.frame(
          Species = names(predictions$actual),
          Actual = predictions$actual,
          Predicted = predictions$predicted,
          Fold = fold_idx,
          stringsAsFactors = FALSE
        )
        method_results$predictions <- rbind(method_results$predictions, pred_df)
        
        # Store errors
        error_df <- data.frame(
          Species = names(predictions$actual),
          Error = errors,
          AbsError = abs_errors,
          SquaredError = squared_errors,
          Fold = fold_idx,
          stringsAsFactors = FALSE
        )
        method_results$errors <- rbind(method_results$errors, error_df)
        
        # Calculate fold metrics
        fold_mae <- mean(abs_errors)
        fold_rmse <- sqrt(mean(squared_errors))
        
        # Calculate correlation if we have enough points
        if(length(errors) >= 3) {
          fold_corr <- cor(predictions$actual, predictions$predicted, method = "spearman")
        } else {
          fold_corr <- NA
        }
        
        # Store fold metrics
        fold_metric <- data.frame(
          Fold = fold_idx,
          MAE = fold_mae,
          RMSE = fold_rmse,
          Correlation = fold_corr,
          Success = TRUE,
          stringsAsFactors = FALSE
        )
        method_results$fold_metrics <- rbind(method_results$fold_metrics, fold_metric)
        
      } else {
        # Record failure
        fold_metric <- data.frame(
          Fold = fold_idx,
          MAE = NA,
          RMSE = NA,
          Correlation = NA,
          Success = FALSE,
          stringsAsFactors = FALSE
        )
        method_results$fold_metrics <- rbind(method_results$fold_metrics, fold_metric)
      }
    }
    
    # Calculate overall metrics
    success_rate <- sum(method_results$fold_metrics$Success) / length(folds)
    
    if(nrow(method_results$errors) > 0) {
      overall_mae <- mean(method_results$errors$AbsError)
      overall_rmse <- sqrt(mean(method_results$errors$SquaredError))
      
      # Calculate overall correlation if we have enough points
      if(nrow(method_results$predictions) >= 5) {
        overall_corr <- cor(method_results$predictions$Actual, 
                           method_results$predictions$Predicted, 
                           method = "spearman")
      } else {
        overall_corr <- NA
      }
    } else {
      overall_mae <- NA
      overall_rmse <- NA
      overall_corr <- NA
    }
    
    # Store results in summary table
    method_idx <- which(cv_results$method_performance$Method == method)
    cv_results$method_performance$MAE[method_idx] <- overall_mae
    cv_results$method_performance$RMSE[method_idx] <- overall_rmse
    cv_results$method_performance$Correlation[method_idx] <- overall_corr
    cv_results$method_performance$Success_Rate[method_idx] <- success_rate
    
    # Store detailed results
    cv_results$detailed_results[[method]] <- method_results
    
    message(sprintf("    CV results: MAE = %.2f, RMSE = %.2f, Success Rate = %.1f%%", 
                   overall_mae, overall_rmse, success_rate * 100))
  }
  
  return(cv_results)
}

#' Evaluate methods using simulation
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param methods Methods to evaluate
#' @param method_params Method parameters
#' @param use_parallel Whether to use parallel processing
#' @param cl Parallel cluster
#' @return Simulation-based evaluation results
#' @keywords internal
evaluate_methods_simulation <- function(tree, chr_counts, methods, method_params, 
                                     use_parallel = FALSE, cl = NULL) {
  # Ensure methods and method_params match
  methods_to_evaluate <- intersect(methods, names(method_params))
  
  if(length(methods_to_evaluate) == 0) {
    stop("No valid methods to evaluate")
  }
  
  # Setup result containers
  sim_results <- list(
    method_performance = data.frame(
      Method = methods_to_evaluate,
      MAE = NA_real_,
      RMSE = NA_real_,
      Correlation = NA_real_,
      Success_Rate = NA_real_,
      stringsAsFactors = FALSE
    ),
    detailed_results = list(),
    simulations = list()
  )
  
  # Generate simulation data
  message("Generating simulation data for method testing...")
  
  # Fit different evolutionary models
  if(requireNamespace("geiger", quietly = TRUE)) {
    # Try to fit BM and OU models
    bm_fit <- try(geiger::fitContinuous(tree, chr_counts, model = "BM"), silent = TRUE)
    ou_fit <- try(geiger::fitContinuous(tree, chr_counts, model = "OU"), silent = TRUE)
    
    # Choose best-fitting model
    if(!inherits(bm_fit, "try-error") && !inherits(ou_fit, "try-error")) {
      if(bm_fit$opt$aicc < ou_fit$opt$aicc) {
        sim_model <- "BM"
        sim_params <- list(sigma = sqrt(bm_fit$opt$sigsq), trend = 0)
        root_value <- bm_fit$opt$z0
      } else {
        sim_model <- "OU"
        sim_params <- list(sigma = sqrt(ou_fit$opt$sigsq), alpha = ou_fit$opt$alpha, theta = ou_fit$opt$z0)
        root_value <- ou_fit$opt$z0
      }
    } else if(!inherits(bm_fit, "try-error")) {
      sim_model <- "BM"
      sim_params <- list(sigma = sqrt(bm_fit$opt$sigsq), trend = 0)
      root_value <- bm_fit$opt$z0
    } else {
      # Default to simple BM
      sim_model <- "BM"
      sim_params <- list(sigma = 1.0, trend = 0)
      root_value <- mean(chr_counts)
    }
  } else {
    # Default to simple BM
    sim_model <- "BM"
    sim_params <- list(sigma = 1.0, trend = 0)
    root_value <- mean(chr_counts)
  }
  
  # Generate simulations
  n_sims <- 10 # Number of simulation replicates
  simulations <- list()
  
  # Check if simulation module is available
  if(exists("simulate_chromosome_evolution", mode = "function")) {
    for(i in 1:n_sims) {
      sim <- simulate_chromosome_evolution(
        tree = tree,
        model = sim_model,
        params = sim_params,
        root_value = root_value,
        constraints = list(min_value = 1),
        discrete = TRUE
      )
      simulations[[i]] <- sim
    }
  } else {
    # Fallback to basic simulation
    warning("Simulation module not loaded. Using basic simulation.")
    for(i in 1:n_sims) {
      sim <- list(
        tree = tree,
        tip_states = NULL,
        node_states = NULL
      )
      
      # Simple BM simulation
      sim$node_states <- numeric(tree$Nnode)
      sim$node_states[1] <- root_value  # Root value
      
      # Simulate down the tree
      sim$tip_states <- phytools::fastBM(tree, sig2 = sim_params$sigma^2, a = root_value)
      sim$tip_states <- pmax(round(sim$tip_states), 1)  # Ensure positive integers
      
      simulations[[i]] <- sim
    }
  }
  
  sim_results$simulations <- simulations
  
  # Test each method on simulated data
  message(sprintf("Testing %d methods on %d simulated datasets...", 
                 length(methods_to_evaluate), length(simulations)))
  
  for(method in methods_to_evaluate) {
    message(sprintf("  Testing method: %s", method))
    
    # Get method parameters
    params <- method_params[[method]]
    
    # Setup result container for this method
    method_results <- list(
      sim_errors = list(),
      sim_metrics = data.frame(
        Simulation = integer(0),
        MAE = numeric(0),
        RMSE = numeric(0),
        Correlation = numeric(0),
        Success = logical(0),
        stringsAsFactors = FALSE
      )
    )
    
    # Process each simulation
    for(sim_idx in 1:length(simulations)) {
      sim <- simulations[[sim_idx]]
      
      # Run reconstruction on simulated tip data
      reconstruction_result <- try({
        run_reconstruction_method(sim$tree, sim$tip_states, method, params)
      }, silent = TRUE)
      
      # Check if reconstruction was successful
      sim_success <- !inherits(reconstruction_result, "try-error") && !is.null(reconstruction_result)
      
      if(sim_success) {
        # Compare reconstructed states with true simulated states
        comparison <- compare_reconstructed_states(reconstruction_result, sim$node_states)
        
        # Store errors
        method_results$sim_errors[[sim_idx]] <- comparison$errors
        
        # Store metrics
        sim_metric <- data.frame(
          Simulation = sim_idx,
          MAE = comparison$mae,
          RMSE = comparison$rmse,
          Correlation = comparison$correlation,
          Success = TRUE,
          stringsAsFactors = FALSE
        )
        method_results$sim_metrics <- rbind(method_results$sim_metrics, sim_metric)
        
      } else {
        # Record failure
        sim_metric <- data.frame(
          Simulation = sim_idx,
          MAE = NA,
          RMSE = NA,
          Correlation = NA,
          Success = FALSE,
          stringsAsFactors = FALSE
        )
        method_results$sim_metrics <- rbind(method_results$sim_metrics, sim_metric)
      }
    }
    
    # Calculate overall metrics
    success_rate <- sum(method_results$sim_metrics$Success) / length(simulations)
    
    successful_metrics <- method_results$sim_metrics[method_results$sim_metrics$Success, ]
    if(nrow(successful_metrics) > 0) {
      overall_mae <- mean(successful_metrics$MAE)
      overall_rmse <- mean(successful_metrics$RMSE)
      overall_corr <- mean(successful_metrics$Correlation)
    } else {
      overall_mae <- NA
      overall_rmse <- NA
      overall_corr <- NA
    }
    
    # Store results in summary table
    method_idx <- which(sim_results$method_performance$Method == method)
    sim_results$method_performance$MAE[method_idx] <- overall_mae
    sim_results$method_performance$RMSE[method_idx] <- overall_rmse
    sim_results$method_performance$Correlation[method_idx] <- overall_corr
    sim_results$method_performance$Success_Rate[method_idx] <- success_rate
    
    # Store detailed results
    sim_results$detailed_results[[method]] <- method_results
    
    message(sprintf("    Simulation results: MAE = %.2f, RMSE = %.2f, Success Rate = %.1f%%", 
                   overall_mae, overall_rmse, success_rate * 100))
  }
  
  return(sim_results)
}

#' Run a specific reconstruction method
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param method Method name
#' @param params Method parameters
#' @return Reconstruction result
#' @keywords internal
run_reconstruction_method <- function(tree, chr_counts, method, params) {
  if(method == "parsimony") {
    if(exists("reconstruct_chromosomes_parsimony", mode = "function")) {
      # Use framework's parsimony function if available
      return(reconstruct_chromosomes_parsimony(
        tree = tree,
        chr_counts = chr_counts,
        method = params$method,
        cost_matrix = params$cost_matrix,
        weight_function = params$weight_function,
        discrete = TRUE
      ))
    } else {
      # Fallback to basic parsimony
      warning("Framework's parsimony function not available. Using basic implementation.")
      anc <- ape::ace(chr_counts, tree, type = "discrete", method = "ML", model = "ER")
      return(list(
        ancestral_states = data.frame(
          node_id = (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode),
          state = anc$ace
        ),
        tree = tree,
        method = "parsimony_basic"
      ))
    }
  } else if(method == "ML") {
    if(exists("reconstruct_chromosomes_ml", mode = "function")) {
      # Use framework's ML function if available
      return(reconstruct_chromosomes_ml(
        tree = tree,
        chr_counts = chr_counts,
        model = params$model,
        discrete = params$discrete
      ))
    } else {
      # Fallback to basic ML
      warning("Framework's ML function not available. Using basic implementation.")
      anc <- ape::ace(chr_counts, tree, type = "continuous", method = "ML", model = "BM")
      return(list(
        ancestral_states = data.frame(
          node_id = (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode),
          state = anc$ace
        ),
        tree = tree,
        method = "ML_basic"
      ))
    }
  } else if(method == "bayesian") {
    if(exists("reconstruct_chromosomes_bayesian", mode = "function")) {
      # Use framework's Bayesian function if available
      return(reconstruct_chromosomes_bayesian(
        tree = tree,
        chr_counts = chr_counts,
        model = params$model,
        mcmc_settings = params$mcmc_settings
      ))
    } else {
      # Fallback to basic Bayesian
      warning("Framework's Bayesian function not available. Using ML as fallback.")
      anc <- ape::ace(chr_counts, tree, type = "continuous", method = "ML", model = "BM")
      return(list(
        ancestral_states = data.frame(
          node_id = (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode),
          state = anc$ace
        ),
        tree = tree,
        method = "bayesian_fallback"
      ))
    }
  } else if(method == "ensemble") {
    if(exists("reconstruct_chromosomes_ensemble", mode = "function")) {
      # Use framework's ensemble function if available
      return(reconstruct_chromosomes_ensemble(
        tree = tree,
        chr_counts = chr_counts,
        weights = params$weights,
        methods = params$methods
      ))
    } else {
      # Fallback to basic ensemble
      warning("Framework's ensemble function not available. Using ML as fallback.")
      anc <- ape::ace(chr_counts, tree, type = "continuous", method = "ML", model = "BM")
      return(list(
        ancestral_states = data.frame(
          node_id = (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode),
          state = anc$ace
        ),
        tree = tree,
        method = "ensemble_fallback"
      ))
    }
  } else {
    stop(paste("Unsupported method:", method))
  }
}

#' Predict chromosome counts for test species
#' 
#' @param tree Full phylogenetic tree
#' @param reconstruction Reconstruction result for training set
#' @param test_species Vector of test species names
#' @return List with actual and predicted counts
#' @keywords internal
predict_test_species <- function(tree, reconstruction, test_species) {
  # Get ancestral states
  anc_states <- reconstruction$ancestral_states
  
  # Get actual values for test species
  test_idx <- match(test_species, tree$tip.label)
  
  # Initialize predictions
  predicted <- numeric(length(test_species))
  names(predicted) <- test_species
  
  # For each test species, find its closest reconstructed ancestor
  for(i in seq_along(test_species)) {
    species <- test_species[i]
    sp_idx <- which(tree$tip.label == species)
    
    if(length(sp_idx) == 1) {
      # Find closest reconstructed ancestor
      ancestors <- find_ancestors(tree, sp_idx)
      
      for(ancestor in ancestors) {
        # Check if this ancestor was reconstructed
        anc_row <- which(anc_states$node_id == ancestor)
        
        if(length(anc_row) == 1) {
          # Use this ancestor's state as prediction
          predicted[i] <- anc_states$state[anc_row]
          break
        }
      }
    }
  }
  
  # Get actual values
  actual <- tree$chr_counts[test_species]
  
  return(list(
    actual = actual,
    predicted = predicted
  ))
}

#' Find ancestors of a tip in a tree
#' 
#' @param tree Phylogenetic tree
#' @param tip_idx Tip index
#' @return Vector of ancestor node indices, starting from most recent
#' @keywords internal
find_ancestors <- function(tree, tip_idx) {
  # Find edges leading to this tip
  edges <- which(tree$edge[, 2] == tip_idx)
  
  if(length(edges) == 0) {
    return(numeric(0))
  }
  
  # Get parent of the tip
  parent <- tree$edge[edges, 1]
  
  # Find all ancestors recursively
  ancestors <- parent
  current <- parent
  
  while(length(current) > 0) {
    # Find edges leading to the current node
    edges <- which(tree$edge[, 2] == current)
    
    if(length(edges) == 0) {
      break
    }
    
    # Get parent of the current node
    current <- tree$edge[edges, 1]
    ancestors <- c(ancestors, current)
  }
  
  return(ancestors)
}

#' Compare reconstructed states with true simulated states
#' 
#' @param reconstruction Reconstruction result
#' @param true_states True simulated node states
#' @return Comparison metrics
#' @keywords internal
compare_reconstructed_states <- function(reconstruction, true_states) {
  # Extract reconstructed states
  anc_states <- reconstruction$ancestral_states
  
  # Ensure node IDs match
  n_nodes <- length(true_states)
  node_ids <- anc_states$node_id
  
  # Extract reconstructed values in same order as true states
  reconstructed <- numeric(n_nodes)
  
  for(i in 1:n_nodes) {
    node_id <- node_ids[i]
    node_idx <- i
    
    if(node_idx <= n_nodes) {
      reconstructed[node_idx] <- anc_states$state[i]
    }
  }
  
  # Calculate error metrics
  errors <- true_states - reconstructed
  mae <- mean(abs(errors))
  rmse <- sqrt(mean(errors^2))
  
  # Calculate correlation if possible
  if(length(true_states) > 2) {
    correlation <- cor(true_states, reconstructed, method = "spearman")
  } else {
    correlation <- NA
  }
  
  return(list(
    errors = errors,
    mae = mae,
    rmse = rmse,
    correlation = correlation
  ))
}

#' Rank methods based on evaluation results
#' 
#' @param results Method comparison results
#' @param criterion Selection criterion
#' @return Ranked data frame of methods
#' @keywords internal
rank_methods <- function(results, criterion) {
  # Create ranking table
  methods <- results$methods
  
  ranking <- data.frame(
    Method = methods,
    CV_Score = NA_real_,
    Sim_Score = NA_real_,
    Overall_Score = NA_real_,
    stringsAsFactors = FALSE
  )
  
  # Extract scores from CV results if available
  if(!is.null(results$cv_results)) {
    cv_perf <- results$cv_results$method_performance
    
    if(criterion == "rmse") {
      cv_scores <- cv_perf$RMSE
      # Lower is better, so invert for scoring
      cv_scores <- max(cv_scores, na.rm = TRUE) - cv_scores
    } else if(criterion == "mae") {
      cv_scores <- cv_perf$MAE
      # Lower is better, so invert for scoring
      cv_scores <- max(cv_scores, na.rm = TRUE) - cv_scores
    } else if(criterion == "correlation") {
      cv_scores <- cv_perf$Correlation
    }
    
    # Normalize to 0-1 range
    cv_range <- range(cv_scores, na.rm = TRUE)
    if(diff(cv_range) > 0) {
      cv_scores_norm <- (cv_scores - cv_range[1]) / diff(cv_range)
    } else {
      cv_scores_norm <- rep(0.5, length(cv_scores))
    }
    
    # Apply success rate penalty
    cv_scores_norm <- cv_scores_norm * cv_perf$Success_Rate
    
    # Store in ranking table
    for(i in 1:length(methods)) {
      method <- methods[i]
      method_idx <- which(cv_perf$Method == method)
      
      if(length(method_idx) == 1) {
        ranking$CV_Score[i] <- cv_scores_norm[method_idx]
      }
    }
  }
  
  # Extract scores from simulation results if available
  if(!is.null(results$simulation_results)) {
    sim_perf <- results$simulation_results$method_performance
    
    if(criterion == "rmse") {
      sim_scores <- sim_perf$RMSE
      # Lower is better, so invert for scoring
      sim_scores <- max(sim_scores, na.rm = TRUE) - sim_scores
    } else if(criterion == "mae") {
      sim_scores <- sim_perf$MAE
      # Lower is better, so invert for scoring
      sim_scores <- max(sim_scores, na.rm = TRUE) - sim_scores
    } else if(criterion == "correlation") {
      sim_scores <- sim_perf$Correlation
    }
    
    # Normalize to 0-1 range
    sim_range <- range(sim_scores, na.rm = TRUE)
    if(diff(sim_range) > 0) {
      sim_scores_norm <- (sim_scores - sim_range[1]) / diff(sim_range)
    } else {
      sim_scores_norm <- rep(0.5, length(sim_scores))
    }
    
    # Apply success rate penalty
    sim_scores_norm <- sim_scores_norm * sim_perf$Success_Rate
    
    # Store in ranking table
    for(i in 1:length(methods)) {
      method <- methods[i]
      method_idx <- which(sim_perf$Method == method)
      
      if(length(method_idx) == 1) {
        ranking$Sim_Score[i] <- sim_scores_norm[method_idx]
      }
    }
  }
  
  # Calculate overall score as average of CV and simulation scores
  for(i in 1:nrow(ranking)) {
    scores <- c(ranking$CV_Score[i], ranking$Sim_Score[i])
    scores <- scores[!is.na(scores)]
    
    if(length(scores) > 0) {
      ranking$Overall_Score[i] <- mean(scores)
    }
  }
  
  # Sort by overall score
  ranking <- ranking[order(ranking$Overall_Score, decreasing = TRUE), ]
  
  return(ranking)
}

#' Create comparison plots
#' 
#' @param results Method comparison results
#' @return List of plot objects
#' @keywords internal
create_comparison_plots <- function(results) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package is required for creating comparison plots")
    return(list())
  }
  
  plots <- list()
  
  # Create overall ranking plot
  if(!is.null(results$overall_ranking)) {
    ranking <- results$overall_ranking
    
    # Melt the data for plotting
    plot_data <- reshape2::melt(
      ranking[, c("Method", "CV_Score", "Sim_Score", "Overall_Score")],
      id.vars = "Method",
      variable.name = "Metric",
      value.name = "Score"
    )
    
    # Clean up metric names
    plot_data$Metric <- gsub("_Score", "", plot_data$Metric)
    
    # Create barplot
    plots$ranking <- ggplot2::ggplot(
      plot_data, 
      ggplot2::aes(x = reorder(Method, Score), y = Score, fill = Metric)
    ) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::labs(
        title = "Method Performance Comparison",
        subtitle = "Higher scores indicate better performance",
        x = "Method",
        y = "Performance Score"
      ) +
      ggplot2::scale_fill_brewer(palette = "Set1") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  # Create cross-validation error plot
  if(!is.null(results$cv_results)) {
    # Combine all predictions
    all_predictions <- data.frame()
    
    for(method in names(results$cv_results$detailed_results)) {
      method_preds <- results$cv_results$detailed_results[[method]]$predictions
      if(nrow(method_preds) > 0) {
        method_preds$Method <- method
        all_predictions <- rbind(all_predictions, method_preds)
      }
    }
    
    if(nrow(all_predictions) > 0) {
      # Create scatterplot of predicted vs actual
      plots$cv_predictions <- ggplot2::ggplot(
        all_predictions, 
        ggplot2::aes(x = Actual, y = Predicted, color = Method)
      ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dashed") +
        ggplot2::labs(
          title = "Cross-Validation: Predicted vs Actual",
          x = "Actual Chromosome Count",
          y = "Predicted Chromosome Count"
        ) +
        ggplot2::facet_wrap(~Method) +
        ggplot2::theme_minimal()
      
      # Create boxplot of errors
      all_errors <- data.frame()
      
      for(method in names(results$cv_results$detailed_results)) {
        method_errors <- results$cv_results$detailed_results[[method]]$errors
        if(nrow(method_errors) > 0) {
          method_errors$Method <- method
          all_errors <- rbind(all_errors, method_errors)
        }
      }
      
      if(nrow(all_errors) > 0) {
        plots$cv_errors <- ggplot2::ggplot(
          all_errors, 
          ggplot2::aes(x = Method, y = AbsError, fill = Method)
        ) +
          ggplot2::geom_boxplot() +
          ggplot2::labs(
            title = "Cross-Validation: Absolute Errors",
            x = "Method",
            y = "Absolute Error"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      }
    }
  }
  
  # Create simulation error plot
  if(!is.null(results$simulation_results)) {
    # Combine all simulation metrics
    all_metrics <- data.frame()
    
    for(method in names(results$simulation_results$detailed_results)) {
      method_metrics <- results$simulation_results$detailed_results[[method]]$sim_metrics
      if(nrow(method_metrics) > 0) {
        method_metrics$Method <- method
        all_metrics <- rbind(all_metrics, method_metrics)
      }
    }
    
    if(nrow(all_metrics) > 0) {
      # Create boxplot of simulation errors
      plots$sim_errors <- ggplot2::ggplot(
        all_metrics, 
        ggplot2::aes(x = Method, y = MAE, fill = Method)
      ) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(
          title = "Simulation Testing: Mean Absolute Errors",
          x = "Method",
          y = "Mean Absolute Error"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  }
  
  return(plots)
}

#===============================================================================
# High-Level Method Workflow Functions
#===============================================================================

#' Run full method comparison workflow
#' 
#' Comprehensive workflow to compare methods and select the best option
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param methods Methods to compare
#' @param output_dir Output directory for results
#' @param run_cv Whether to run cross-validation
#' @param run_simulation Whether to run simulation-based testing
#' @param create_plots Whether to create comparison plots
#' @param selection_criterion Criterion for selecting best method
#' @param export_results Whether to export results to files
#' @return Method comparison results
#' @export
run_method_comparison_workflow <- function(tree, chr_counts,
                                         methods = c("parsimony", "ML", "bayesian"),
                                         output_dir = "method_comparison",
                                         run_cv = TRUE,
                                         run_simulation = TRUE,
                                         create_plots = TRUE,
                                         selection_criterion = "rmse",
                                         export_results = TRUE) {
  # Create output directory if needed
  if(export_results && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set evaluation approach based on parameters
  if(run_cv && run_simulation) {
    eval_approach <- "both"
  } else if(run_cv) {
    eval_approach <- "cv"
  } else if(run_simulation) {
    eval_approach <- "simulation"
  } else {
    stop("At least one of run_cv or run_simulation must be TRUE")
  }
  
  # Run method comparison
  results <- compare_reconstruction_methods(
    tree = tree,
    chr_counts = chr_counts,
    methods = methods,
    eval_approach = eval_approach,
    selection_criterion = selection_criterion,
    create_plots = create_plots
  )
  
  # Export results if requested
  if(export_results) {
    # Save RDS of full results
    saveRDS(results, file.path(output_dir, "method_comparison_results.rds"))
    
    # Export tables
    if(!is.null(results$overall_ranking)) {
      write.csv(results$overall_ranking, file.path(output_dir, "method_ranking.csv"), 
                row.names = FALSE)
    }
    
    if(!is.null(results$cv_results)) {
      write.csv(results$cv_results$method_performance, 
                file.path(output_dir, "cv_performance.csv"), row.names = FALSE)
    }
    
    if(!is.null(results$simulation_results)) {
      write.csv(results$simulation_results$method_performance, 
                file.path(output_dir, "simulation_performance.csv"), row.names = FALSE)
    }
    
    # Export plots
    if(create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
      plots_dir <- file.path(output_dir, "plots")
      if(!dir.exists(plots_dir)) {
        dir.create(plots_dir)
      }
      
      for(plot_name in names(results$plots)) {
        ggplot2::ggsave(
          filename = file.path(plots_dir, paste0(plot_name, ".pdf")),
          plot = results$plots[[plot_name]],
          width = 8,
          height = 6
        )
      }
    }
    
    # Create a summary report
    report_file <- file.path(output_dir, "method_comparison_report.txt")
    cat("Method Comparison Report\n", file = report_file)
    cat("========================\n\n", file = report_file, append = TRUE)
    
    cat("Best Method: ", results$best_method, "\n\n", file = report_file, append = TRUE)
    
    if(!is.null(results$overall_ranking)) {
      cat("Method Ranking:\n", file = report_file, append = TRUE)
      capture.output(results$overall_ranking, file = report_file, append = TRUE)
      cat("\n\n", file = report_file, append = TRUE)
    }
    
    if(!is.null(results$cv_results)) {
      cat("Cross-Validation Performance:\n", file = report_file, append = TRUE)
      capture.output(results$cv_results$method_performance, file = report_file, append = TRUE)
      cat("\n\n", file = report_file, append = TRUE)
    }
    
    if(!is.null(results$simulation_results)) {
      cat("Simulation Performance:\n", file = report_file, append = TRUE)
      capture.output(results$simulation_results$method_performance, file = report_file, append = TRUE)
    }
  }
  
  # Return the results
  return(results)
}

#' Apply the best reconstruction method based on comparison
#' 
#' @param method_comparison Method comparison results
#' @param tree Phylogenetic tree (NULL to use tree from comparison)
#' @param chr_counts Chromosome counts (NULL to use data from comparison)
#' @param custom_params Custom parameters for selected method
#' @return Reconstruction results from best method
#' @export
apply_best_method <- function(method_comparison, tree = NULL, chr_counts = NULL, custom_params = NULL) {
  # Extract best method
  best_method <- method_comparison$best_method
  
  if(is.null(best_method) || is.na(best_method)) {
    stop("No best method identified in method comparison results")
  }
  
  # Use tree and data from comparison if not provided
  if(is.null(tree)) {
    tree <- method_comparison$tree
  }
  
  if(is.null(chr_counts)) {
    chr_counts <- method_comparison$data
  }
  
  # Get default parameters for the best method
  if(is.null(custom_params)) {
    params <- method_comparison$method_params[[best_method]]
  } else {
    params <- custom_params
  }
  
  message(sprintf("Applying best method (%s) to reconstruct ancestral chromosome numbers...", best_method))
  
  # Run the selected method
  reconstruction <- run_reconstruction_method(tree, chr_counts, best_method, params)
  
  # Add method selection information
  reconstruction$selection_info <- list(
    comparison_summary = method_comparison$overall_ranking,
    cv_results = if(!is.null(method_comparison$cv_results)) method_comparison$cv_results$method_performance else NULL,
    sim_results = if(!is.null(method_comparison$simulation_results)) method_comparison$simulation_results$method_performance else NULL
  )
  
  return(reconstruction)
}
