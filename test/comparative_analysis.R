#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Comparative Analysis Module
# Author: Bioinformatics Team
# Date: 2025-05-15
# Description: Implements methods for comparative analysis of chromosome 
#              number evolution across clades, testing evolutionary hypotheses,
#              and analyzing correlations with other traits
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(geiger)
  library(picante)
  library(caper)
  library(nlme)
  library(diversitree)
  library(parallel)
  library(ggplot2)
  library(RColorBrewer)
})

#===============================================================================
# Clade Comparison Functions
#===============================================================================

#' Compare chromosome number evolution between clades
#' 
#' Performs statistical comparisons of chromosome number patterns between
#' specified clades in a phylogeny
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param clades List of clades, each a character vector of taxon names
#' @param clade_colors Optional named vector of colors for each clade
#' @param metrics Metrics to compare: "mean", "variance", "range", "signal", "rate", "all"
#' @param rate_method Method for estimating evolutionary rate: "contrasts", "model" 
#' @param signal_method Method for estimating phylogenetic signal: "K", "lambda"
#' @param n_replicates Number of replicates for randomization tests
#' @param n_cores Number of CPU cores for parallel processing (NULL = auto-detect)
#' @return List with comparison results and statistics
#' @export
compare_chromosome_evolution <- function(tree, 
                                       chr_counts, 
                                       clades,
                                       clade_colors = NULL,
                                       metrics = "all",
                                       rate_method = "contrasts",
                                       signal_method = "K",
                                       n_replicates = 1000,
                                       n_cores = NULL) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Validate clades
  if(!is.list(clades)) {
    stop("clades must be a list of taxon name vectors")
  }
  
  if(is.null(names(clades))) {
    names(clades) <- paste0("Clade_", seq_along(clades))
  }
  
  # Check that all clades have at least some representatives in the tree
  for(clade_name in names(clades)) {
    clade_taxa <- clades[[clade_name]]
    tree_match <- intersect(clade_taxa, tree$tip.label)
    
    if(length(tree_match) == 0) {
      stop(paste("No taxa from clade", clade_name, "found in the phylogeny"))
    }
  }
  
  # Determine which metrics to compute
  all_metrics <- c("mean", "variance", "range", "signal", "rate")
  if(length(metrics) == 1 && metrics == "all") {
    metrics <- all_metrics
  } else {
    metrics <- match.arg(metrics, all_metrics, several.ok = TRUE)
  }
  
  # Set up parallel processing
  if(is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  # Initialize results
  results <- list(
    summary = list(
      clades = names(clades),
      n_clades = length(clades),
      metrics = metrics,
      has_significant_differences = FALSE
    ),
    clade_data = list(),
    statistics = list(),
    comparisons = list(),
    permutations = list(),
    plots = list(),
    tree = tree
  )
  
  # Define colors for each clade if not provided
  if(is.null(clade_colors)) {
    n_clades <- length(clades)
    if(n_clades <= 8) {
      clade_colors <- RColorBrewer::brewer.pal(max(3, n_clades), "Dark2")[1:n_clades]
    } else {
      clade_colors <- rainbow(n_clades)
    }
    names(clade_colors) <- names(clades)
  }
  results$clade_colors <- clade_colors
  
  # Extract data for each clade
  clade_trees <- list()
  clade_counts <- list()
  
  for(clade_name in names(clades)) {
    clade_taxa <- clades[[clade_name]]
    
    # Get intersection of clade taxa and tree tips
    shared_tips <- intersect(clade_taxa, tree$tip.label)
    
    # Extract subtree for this clade
    clade_tree <- ape::keep.tip(tree, shared_tips)
    clade_trees[[clade_name]] <- clade_tree
    
    # Extract chr counts for this clade
    shared_taxa_with_data <- intersect(shared_tips, names(chr_counts))
    clade_chr_counts <- chr_counts[shared_taxa_with_data]
    clade_counts[[clade_name]] <- clade_chr_counts
    
    # Store basic information
    results$clade_data[[clade_name]] <- list(
      n_taxa = length(shared_tips),
      n_with_data = length(shared_taxa_with_data),
      percent_complete = 100 * length(shared_taxa_with_data) / length(shared_tips),
      min_chr = min(clade_chr_counts, na.rm = TRUE),
      max_chr = max(clade_chr_counts, na.rm = TRUE),
      mean_chr = mean(clade_chr_counts, na.rm = TRUE),
      median_chr = median(clade_chr_counts, na.rm = TRUE),
      sd_chr = sd(clade_chr_counts, na.rm = TRUE),
      tree = clade_tree,
      chr_counts = clade_chr_counts
    )
  }
  
  # Calculate statistics for each metric
  metric_results <- list()
  
  # Calculate mean chromosome numbers
  if("mean" %in% metrics) {
    mean_by_clade <- sapply(clade_counts, function(x) mean(x, na.rm = TRUE))
    
    # ANOVA test for differences in means
    means_data <- data.frame(
      Chromosome_Count = unlist(clade_counts),
      Clade = rep(names(clade_counts), sapply(clade_counts, length)),
      stringsAsFactors = FALSE
    )
    
    means_anova <- anova(lm(Chromosome_Count ~ Clade, data = means_data))
    
    # Store results
    metric_results$mean <- list(
      values = mean_by_clade,
      anova = means_anova,
      p_value = means_anova$`Pr(>F)`[1],
      significant = means_anova$`Pr(>F)`[1] < 0.05
    )
    
    # Create means plot
    if(requireNamespace("ggplot2", quietly = TRUE)) {
      mean_plot <- ggplot2::ggplot(means_data, ggplot2::aes(x = Clade, y = Chromosome_Count, fill = Clade)) +
        ggplot2::geom_boxplot() +
        ggplot2::scale_fill_manual(values = clade_colors) +
        ggplot2::labs(
          title = "Chromosome Number Distribution by Clade",
          subtitle = paste("ANOVA p-value =", format(means_anova$`Pr(>F)`[1], digits = 3)),
          x = "Clade",
          y = "Chromosome Number"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")
      
      results$plots$mean <- mean_plot
    }
  }
  
  # Calculate variance in chromosome numbers
  if("variance" %in% metrics) {
    var_by_clade <- sapply(clade_counts, function(x) var(x, na.rm = TRUE))
    
    # Bartlett's test for homogeneity of variance
    var_test <- bartlett.test(Chromosome_Count ~ Clade, data = means_data)
    
    # Store results
    metric_results$variance <- list(
      values = var_by_clade,
      bartlett = var_test,
      p_value = var_test$p.value,
      significant = var_test$p.value < 0.05
    )
    
    # Create variance plot
    if(requireNamespace("ggplot2", quietly = TRUE)) {
      var_data <- data.frame(
        Clade = names(var_by_clade),
        Variance = var_by_clade,
        stringsAsFactors = FALSE
      )
      
      var_plot <- ggplot2::ggplot(var_data, ggplot2::aes(x = Clade, y = Variance, fill = Clade)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = clade_colors) +
        ggplot2::labs(
          title = "Variance in Chromosome Numbers by Clade",
          subtitle = paste("Bartlett's test p-value =", format(var_test$p.value, digits = 3)),
          x = "Clade",
          y = "Variance"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")
      
      results$plots$variance <- var_plot
    }
  }
  
  # Calculate range in chromosome numbers
  if("range" %in% metrics) {
    range_by_clade <- sapply(clade_counts, function(x) {
      r <- range(x, na.rm = TRUE)
      r[2] - r[1]
    })
    
    # Randomization test for differences in range
    range_pvalue <- perform_range_randomization_test(clade_counts, n_replicates, n_cores)
    
    # Store results
    metric_results$range <- list(
      values = range_by_clade,
      p_value = range_pvalue,
      significant = range_pvalue < 0.05
    )
    
    # Create range plot
    if(requireNamespace("ggplot2", quietly = TRUE)) {
      range_data <- data.frame(
        Clade = names(range_by_clade),
        Range = range_by_clade,
        stringsAsFactors = FALSE
      )
      
      range_plot <- ggplot2::ggplot(range_data, ggplot2::aes(x = Clade, y = Range, fill = Clade)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = clade_colors) +
        ggplot2::labs(
          title = "Range of Chromosome Numbers by Clade",
          subtitle = paste("Permutation test p-value =", format(range_pvalue, digits = 3)),
          x = "Clade",
          y = "Range"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                     legend.position = "none")
      
      results$plots$range <- range_plot
    }
  }
  
  # Store metric results
  results$statistics <- metric_results
  
  # Create summary of results
  results$summary$has_significant_differences <- any(sapply(metric_results, function(x) 
    !is.null(x$significant) && !is.na(x$significant) && x$significant))
  
  return(results)
}

#' Perform randomization test for difference in chromosome range between clades
#' 
#' @param clade_counts List of chromosome counts by clade
#' @param n_replicates Number of replicates for randomization test
#' @param n_cores Number of CPU cores for parallel processing
#' @return P-value from randomization test
#' @keywords internal
perform_range_randomization_test <- function(clade_counts, n_replicates, n_cores) {
  # Calculate observed range differences
  ranges <- sapply(clade_counts, function(x) {
    r <- range(x, na.rm = TRUE)
    r[2] - r[1]
  })
  
  # Calculate test statistic: variance in ranges between clades
  observed_stat <- var(ranges, na.rm = TRUE)
  
  # All data combined
  all_data <- unlist(clade_counts)
  clade_sizes <- sapply(clade_counts, length)
  
  # Parallel function for randomization
  random_replicate <- function(i) {
    # Randomly shuffle all data
    shuffled <- sample(all_data)
    
    # Split into original clade sizes
    start <- 1
    shuffled_clades <- list()
    for(j in seq_along(clade_sizes)) {
      end <- start + clade_sizes[j] - 1
      shuffled_clades[[j]] <- shuffled[start:end]
      start <- end + 1
    }
    
    # Calculate ranges of shuffled data
    shuffled_ranges <- sapply(shuffled_clades, function(x) {
      r <- range(x, na.rm = TRUE)
      r[2] - r[1]
    })
    
    # Return test statistic
    var(shuffled_ranges, na.rm = TRUE)
  }
  
  # Run randomization in parallel
  if(n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    random_stats <- parallel::parSapply(cl, 1:n_replicates, random_replicate)
  } else {
    random_stats <- sapply(1:n_replicates, random_replicate)
  }
  
  # Calculate p-value
  p_value <- sum(random_stats >= observed_stat) / n_replicates
  
  return(p_value)
}

#' Calculate mean squared contrasts as a measure of evolutionary rate
#' 
#' @param tree Phylogenetic tree
#' @param trait Trait values for tips
#' @return Mean squared contrasts value
#' @keywords internal
calculate_mean_squared_contrasts <- function(tree, trait) {
  # Calculate contrasts
  contrasts <- ape::pic(trait, tree)
  
  # Return mean squared contrasts
  mean(contrasts^2, na.rm = TRUE)
}

#===============================================================================
# Trait Correlation Functions
#===============================================================================

#' Analyze correlation between chromosome numbers and other traits
#' 
#' Tests for evolutionary correlation between chromosome numbers and 
#' other continuous or discrete traits, accounting for phylogeny
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param traits Data frame or matrix of other traits, row names should match tree tip labels
#' @param method Correlation method: "PGLS", "pic", "gls", or "phyloANOVA"
#' @param transform Transform chromosome counts: "none", "log", "sqrt"
#' @param discrete_traits Character vector of trait names that should be treated as discrete
#' @param control_variables Character vector of traits to use as control variables
#' @param alpha Significance level for p-value
#' @param adjust_method Method for adjusting p-values for multiple comparisons
#' @return List with correlation analysis results
#' @export
analyze_trait_correlations <- function(tree, 
                                     chr_counts, 
                                     traits,
                                     method = "PGLS",
                                     transform = "none",
                                     discrete_traits = NULL,
                                     control_variables = NULL,
                                     alpha = 0.05,
                                     adjust_method = "fdr") {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  if(!is.data.frame(traits) && !is.matrix(traits)) {
    stop("traits must be a data frame or matrix")
  }
  
  # Check for matching taxa in tree, chromosome data, and trait data
  tree_tips <- tree$tip.label
  chr_taxa <- names(chr_counts)
  trait_taxa <- rownames(traits)
  
  if(is.null(trait_taxa)) {
    stop("traits must have row names matching species names")
  }
  
  # Get common taxa across all datasets
  common_taxa <- Reduce(intersect, list(tree_tips, chr_taxa, trait_taxa))
  
  if(length(common_taxa) < 4) {
    stop("At least 4 taxa must have chromosome counts, trait data, and be present in the tree")
  }
  
  # Subset data to common taxa
  pruned_tree <- ape::keep.tip(tree, common_taxa)
  pruned_chr <- chr_counts[common_taxa]
  pruned_traits <- traits[common_taxa, , drop = FALSE]
  
  # Apply transformation to chromosome counts if requested
  if(transform == "log") {
    pruned_chr <- log(pruned_chr)
  } else if(transform == "sqrt") {
    pruned_chr <- sqrt(pruned_chr)
  }
  
  # Initialize results
  results <- list(
    summary = list(
      method = method,
      transform = transform,
      n_taxa = length(common_taxa),
      n_traits = ncol(pruned_traits),
      n_significant = 0,
      significant_traits = character(0),
      alpha = alpha
    ),
    correlations = list(),
    plots = list()
  )
  
  # Analyze correlation for each trait
  trait_names <- colnames(pruned_traits)
  p_values <- numeric(length(trait_names))
  estimates <- numeric(length(trait_names))
  
  for(i in seq_along(trait_names)) {
    trait_name <- trait_names[i]
    trait_values <- pruned_traits[, trait_name]
    
    # Check if trait is discrete
    is_discrete <- !is.null(discrete_traits) && trait_name %in% discrete_traits
    
    # Calculate correlation based on method and whether trait is discrete
    if(is_discrete) {
      # For discrete traits
      if(method %in% c("PGLS", "gls", "pic")) {
        # Convert to factor
        trait_factor <- as.factor(trait_values)
        
        if(method == "PGLS") {
          # Use PGLS with factor as predictor
          if(requireNamespace("caper", quietly = TRUE)) {
            comp_data <- caper::comparative.data(pruned_tree, 
                                               data.frame(species = common_taxa, 
                                                         chr = pruned_chr, 
                                                         trait = trait_factor, 
                                                         stringsAsFactors = FALSE), 
                                               "species")
            pgls_model <- tryCatch({
              caper::pgls(chr ~ trait, data = comp_data)
            }, error = function(e) {
              warning(paste("PGLS error for trait", trait_name, ":", e$message))
              return(NULL)
            })
            
            if(!is.null(pgls_model)) {
              p_values[i] <- anova(pgls_model)$"Pr(>F)"[2]
              estimates[i] <- NA  # No single estimate for factor
            }
          }
        } else if(method == "gls") {
          # Use GLS with phylogenetic correlation structure
          if(requireNamespace("nlme", quietly = TRUE) && requireNamespace("ape", quietly = TRUE)) {
            vcv_matrix <- ape::vcv(pruned_tree)
            gls_model <- tryCatch({
              nlme::gls(chr ~ trait_factor, 
                       correlation = nlme::corSymm(ape::vcv2cor(vcv_matrix)[lower.tri(vcv_matrix)]))
            }, error = function(e) {
              warning(paste("GLS error for trait", trait_name, ":", e$message))
              return(NULL)
            })
            
            if(!is.null(gls_model)) {
              p_values[i] <- anova(gls_model)$"p-value"[2]
              estimates[i] <- NA  # No single estimate for factor
            }
          }
        } else {  # method == "pic"
          warning("PIC method not suitable for discrete traits")
        }
      } else if(method == "phyloANOVA") {
        # Phylogenetic ANOVA for discrete traits
        if(requireNamespace("phytools", quietly = TRUE)) {
          phyanova_result <- tryCatch({
            phytools::phylANOVA(pruned_tree, trait_factor, pruned_chr)
          }, error = function(e) {
            warning(paste("phylANOVA error for trait", trait_name, ":", e$message))
            return(NULL)
          })
          
          if(!is.null(phyanova_result)) {
            p_values[i] <- phyanova_result$Pf
            estimates[i] <- NA  # No single estimate for ANOVA
          }
        }
      }
    } else {
      # For continuous traits
      if(method == "PGLS") {
        # PGLS for continuous traits
        if(requireNamespace("caper", quietly = TRUE)) {
          comp_data <- caper::comparative.data(pruned_tree, 
                                             data.frame(species = common_taxa, 
                                                       chr = pruned_chr, 
                                                       trait = trait_values, 
                                                       stringsAsFactors = FALSE), 
                                             "species")
          pgls_model <- tryCatch({
            caper::pgls(chr ~ trait, data = comp_data)
          }, error = function(e) {
            warning(paste("PGLS error for trait", trait_name, ":", e$message))
            return(NULL)
          })
          
          if(!is.null(pgls_model)) {
            p_values[i] <- summary(pgls_model)$coefficients["trait", "Pr(>|t|)"]
            estimates[i] <- summary(pgls_model)$coefficients["trait", "Estimate"]
          }
        }
      } else if(method == "pic") {
        # Phylogenetic independent contrasts
        if(requireNamespace("ape", quietly = TRUE)) {
          chr_pics <- ape::pic(pruned_chr, pruned_tree)
          trait_pics <- ape::pic(trait_values, pruned_tree)
          
          pic_cor <- tryCatch({
            cor.test(chr_pics, trait_pics)
          }, error = function(e) {
            warning(paste("PIC error for trait", trait_name, ":", e$message))
            return(NULL)
          })
          
          if(!is.null(pic_cor)) {
            p_values[i] <- pic_cor$p.value
            estimates[i] <- pic_cor$estimate
          }
        }
      } else if(method == "gls") {
        # GLS with phylogenetic correlation structure
        if(requireNamespace("nlme", quietly = TRUE) && requireNamespace("ape", quietly = TRUE)) {
          vcv_matrix <- ape::vcv(pruned_tree)
          gls_model <- tryCatch({
            nlme::gls(chr ~ trait_values, 
                     correlation = nlme::corSymm(ape::vcv2cor(vcv_matrix)[lower.tri(vcv_matrix)]))
          }, error = function(e) {
            warning(paste("GLS error for trait", trait_name, ":", e$message))
            return(NULL)
          })
          
          if(!is.null(gls_model)) {
            p_values[i] <- summary(gls_model)$tTable["trait_values", "p-value"]
            estimates[i] <- summary(gls_model)$tTable["trait_values", "Value"]
          }
        }
      }
    }
    
    # Store correlation results for this trait
    results$correlations[[trait_name]] <- list(
      trait_name = trait_name,
      is_discrete = is_discrete,
      p_value = p_values[i],
      estimate = estimates[i]
    )
  }
  
  # Adjust p-values for multiple comparisons
  adjusted_p <- p.adjust(p_values, method = adjust_method)
  
  # Update results with adjusted p-values
  for(i in seq_along(trait_names)) {
    trait_name <- trait_names[i]
    results$correlations[[trait_name]]$adjusted_p_value <- adjusted_p[i]
    results$correlations[[trait_name]]$significant <- adjusted_p[i] < alpha
  }
  
  # Update summary
  significant_traits <- trait_names[adjusted_p < alpha]
  results$summary$n_significant <- length(significant_traits)
  results$summary$significant_traits <- significant_traits
  results$summary$adjustment_method <- adjust_method
  
  return(results)
}

#===============================================================================
# Model-Based Comparative Analysis
#===============================================================================

#' Test and compare evolutionary models for chromosome number evolution
#' 
#' Fits and compares different evolutionary models to test hypotheses about
#' chromosome number evolution
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param models Vector of models to fit: "BM", "OU", "EB", "trend", "white", "lambda", "kappa", "delta"
#' @param discrete Whether to treat chromosome counts as discrete
#' @param constrain_params List of parameter constraints for models
#' @param criterion Model selection criterion: "AIC", "AICc", "BIC"
#' @param transform Transform chromosome counts: "none", "log", "sqrt"
#' @param clades Optional list of clades for clade-specific models
#' @return List with model fitting and comparison results
#' @export
test_chromosome_evolution_models <- function(tree, 
                                          chr_counts, 
                                          models = c("BM", "OU", "EB", "trend"),
                                          discrete = TRUE,
                                          constrain_params = NULL,
                                          criterion = "AICc",
                                          transform = "none",
                                          clades = NULL) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Prune tree to match available data
  common_taxa <- intersect(tree$tip.label, names(chr_counts))
  
  if(length(common_taxa) < 4) {
    stop("At least 4 taxa must have chromosome counts and be present in the tree")
  }
  
  pruned_tree <- ape::keep.tip(tree, common_taxa)
  pruned_chr <- chr_counts[common_taxa]
  
  # Apply transformation if requested
  if(transform == "log") {
    pruned_chr <- log(pruned_chr)
  } else if(transform == "sqrt") {
    pruned_chr <- sqrt(pruned_chr)
  }
  
  # Initialize results
  results <- list(
    model_fits = list(),
    model_comparison = NULL,
    best_model = NULL,
    discrete = discrete,
    transform = transform,
    criterion = criterion,
    summary = NULL,
    plots = list()
  )
  
  # Fit continuous trait models if not discrete
  if(!discrete || transform != "none") {
    for(model in models) {
      # Different models require different packages
      model_fit <- NULL
      
      if(model %in% c("BM", "OU", "EB", "trend")) {
        if(requireNamespace("geiger", quietly = TRUE)) {
          # Map model names to geiger model names
          geiger_model <- switch(model,
                               "BM" = "BM",
                               "OU" = "OU",
                               "EB" = "EB",
                               "trend" = "trend")
          
          # Fit model
          model_fit <- tryCatch({
            geiger::fitContinuous(pruned_tree, pruned_chr, model = geiger_model)
          }, error = function(e) {
            warning(paste("Error fitting", model, "model:", e$message))
            return(NULL)
          })
        }
      } else if(model %in% c("lambda", "kappa", "delta")) {
        if(requireNamespace("phytools", quietly = TRUE)) {
          # Fit phylogenetic signal models
          model_fit <- tryCatch({
            phytools::phylosig(pruned_tree, pruned_chr, method = model, test = TRUE)
          }, error = function(e) {
            warning(paste("Error fitting", model, "model:", e$message))
            return(NULL)
          })
        }
      } else if(model == "white") {
        if(requireNamespace("geiger", quietly = TRUE)) {
          # White noise (no phylogenetic structure)
          model_fit <- tryCatch({
            geiger::fitContinuous(pruned_tree, pruned_chr, model = "white")
          }, error = function(e) {
            warning(paste("Error fitting white noise model:", e$message))
            return(NULL)
          })
        }
      }
      
      results$model_fits[[model]] <- model_fit
    }
    
    # Compare models
    model_comparison <- compare_fitted_models(results$model_fits, criterion = criterion)
    results$model_comparison <- model_comparison
    
    # Identify best model
    if(!is.null(model_comparison) && nrow(model_comparison) > 0) {
      best_idx <- which.min(model_comparison[[criterion]])
      results$best_model <- model_comparison$Model[best_idx]
    }
  } else {
    # TO BE IMPLEMENTED: Discrete trait evolution models
    warning("Discrete trait models not yet implemented")
  }
  
  # Create summary of results
  results$summary <- summarize_model_results(results)
  
  return(results)
}

#' Compare fitted models based on AIC/BIC
#' 
#' @param model_fits List of fitted models
#' @param criterion Model selection criterion: "AIC", "AICc", "BIC"
#' @return Data frame with model comparison statistics
#' @keywords internal
compare_fitted_models <- function(model_fits, criterion = "AICc") {
  # Initialize result data frame
  model_comparison <- data.frame(
    Model = character(0),
    LogLik = numeric(0),
    Parameters = numeric(0),
    AIC = numeric(0),
    AICc = numeric(0),
    BIC = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Extract information from each fitted model
  for(model_name in names(model_fits)) {
    model_fit <- model_fits[[model_name]]
    
    if(is.null(model_fit)) {
      next
    }
    
    # Handle different model fit structures
    if(inherits(model_fit, "fitContinuous")) {
      # geiger model
      loglik <- model_fit$opt$lnL
      params <- model_fit$opt$k
      sample_size <- length(model_fit$res)
      
      aic <- -2 * loglik + 2 * params
      aicc <- aic + (2 * params * (params + 1)) / (sample_size - params - 1)
      bic <- -2 * loglik + params * log(sample_size)
      
      model_comparison <- rbind(model_comparison, data.frame(
        Model = model_name,
        LogLik = loglik,
        Parameters = params,
        AIC = aic,
        AICc = aicc,
        BIC = bic,
        stringsAsFactors = FALSE
      ))
    } else if(inherits(model_fit, "phylosig")) {
      # phytools phylogenetic signal model
      loglik <- model_fit$logL
      params <- 2  # Assuming 2 parameters (signal parameter + root state)
      sample_size <- length(model_fit$x)
      
      aic <- -2 * loglik + 2 * params
      aicc <- aic + (2 * params * (params + 1)) / (sample_size - params - 1)
      bic <- -2 * loglik + params * log(sample_size)
      
      model_comparison <- rbind(model_comparison, data.frame(
        Model = model_name,
        LogLik = loglik,
        Parameters = params,
        AIC = aic,
        AICc = aicc,
        BIC = bic,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Sort by criterion
  if(nrow(model_comparison) > 0) {
    model_comparison <- model_comparison[order(model_comparison[[criterion]]), ]
    
    # Calculate delta and weight
    model_comparison[[paste0("d", criterion)]] <- model_comparison[[criterion]] - min(model_comparison[[criterion]])
    model_comparison[[paste0(criterion, "_weight")]] <- exp(-0.5 * model_comparison[[paste0("d", criterion)]]) / 
      sum(exp(-0.5 * model_comparison[[paste0("d", criterion)]]))
  }
  
  return(model_comparison)
}

#' Summarize model comparison results
#' 
#' @param model_results Model fitting results
#' @return Summary of model comparison
#' @keywords internal
summarize_model_results <- function(model_results) {
  summary <- list()
  
  # Get best model
  best_model <- model_results$best_model
  
  if(!is.null(best_model)) {
    summary$best_model <- best_model
    
    # Get model weights
    if(!is.null(model_results$model_comparison)) {
      comp <- model_results$model_comparison
      criterion <- model_results$criterion
      weight_col <- paste0(criterion, "_weight")
      
      if(weight_col %in% colnames(comp)) {
        best_idx <- which(comp$Model == best_model)
        best_weight <- comp[[weight_col]][best_idx]
        summary$best_model_weight <- best_weight
        
        # Interpret support
        if(best_weight > 0.9) {
          summary$support_strength <- "Very strong support"
        } else if(best_weight > 0.8) {
          summary$support_strength <- "Strong support"
        } else if(best_weight > 0.6) {
          summary$support_strength <- "Moderate support"
        } else {
          summary$support_strength <- "Weak support"
        }
      }
    }
    
    # Interpret best model
    summary$model_interpretation <- interpret_model(best_model, model_results$model_fits[[best_model]])
  } else {
    summary$best_model <- NA
    summary$model_interpretation <- "No models successfully fit the data"
  }
  
  return(summary)
}

#' Interpret evolutionary model results
#' 
#' @param model_name Name of the model
#' @param model_fit Fitted model object
#' @return Interpretation of model results
#' @keywords internal
interpret_model <- function(model_name, model_fit) {
  if(is.null(model_fit)) {
    return("No model fit available")
  }
  
  interpretation <- ""
  
  if(model_name == "BM") {
    # Brownian motion - random evolution
    sigma2 <- extract_param(model_fit, "sigma2")
    interpretation <- paste0(
      "Chromosome numbers evolve following a random walk (Brownian motion) ",
      "with evolutionary rate sigma² = ", format(sigma2, digits = 3), ". ",
      "This suggests chromosome changes occur randomly with no particular directionality."
    )
  } else if(model_name == "OU") {
    # Ornstein-Uhlenbeck - stabilizing selection
    alpha <- extract_param(model_fit, "alpha")
    theta <- extract_param(model_fit, "theta")
    
    interpretation <- paste0(
      "Chromosome numbers evolve under stabilizing selection (Ornstein-Uhlenbeck) ",
      "with strength alpha = ", format(alpha, digits = 3), " ",
      "toward an optimum value of ", format(theta, digits = 3), ". ",
      "This suggests natural selection favors a particular chromosome number."
    )
  } else if(model_name == "EB") {
    # Early burst - evolution slowing over time
    r <- extract_param(model_fit, "r")
    
    if(r < 0) {
      interpretation <- paste0(
        "Chromosome numbers evolved rapidly early in the tree, with evolutionary rate ",
        "declining over time (r = ", format(r, digits = 3), "). ",
        "This pattern is consistent with adaptive radiation scenarios."
      )
    } else {
      interpretation <- paste0(
        "Chromosome numbers evolved with an accelerating rate over time (r = ", 
        format(r, digits = 3), "). ",
        "This unusual pattern may indicate recent diversification or changing selection pressures."
      )
    }
  } else if(model_name == "trend") {
    # Directional trend
    mu <- extract_param(model_fit, "mu")
    
    if(mu > 0) {
      interpretation <- paste0(
        "Chromosome numbers evolve with a positive directional trend (mu = ", 
        format(mu, digits = 3), "), ",
        "suggesting systematic increase over evolutionary time."
      )
    } else {
      interpretation <- paste0(
        "Chromosome numbers evolve with a negative directional trend (mu = ", 
        format(mu, digits = 3), "), ",
        "suggesting systematic decrease over evolutionary time."
      )
    }
  } else if(model_name == "white") {
    # White noise - no phylogenetic signal
    interpretation <- paste0(
      "Chromosome numbers show no phylogenetic signal (white noise), ",
      "suggesting changes are random with respect to phylogeny. ",
      "This might indicate frequent convergent evolution or rampant hybridization."
    )
  } else if(model_name == "lambda") {
    # Pagel's lambda - phylogenetic signal
    lambda <- extract_param(model_fit, "lambda")
    
    if(lambda > 0.8) {
      interpretation <- paste0(
        "Chromosome numbers show strong phylogenetic signal (lambda = ", 
        format(lambda, digits = 3), "), ",
        "consistent with Brownian motion evolution."
      )
    } else if(lambda > 0.4) {
      interpretation <- paste0(
        "Chromosome numbers show moderate phylogenetic signal (lambda = ", 
        format(lambda, digits = 3), "), ",
        "suggesting partial phylogenetic constraint on evolution."
      )
    } else {
      interpretation <- paste0(
        "Chromosome numbers show weak phylogenetic signal (lambda = ", 
        format(lambda, digits = 3), "), ",
        "suggesting evolution occurs with little regard to phylogeny."
      )
    }
  } else if(model_name == "delta") {
    # Pagel's delta - tempo of trait evolution
    delta <- extract_param(model_fit, "delta")
    
    if(delta > 1) {
      interpretation <- paste0(
        "Chromosome evolution accelerated over time (delta = ", 
        format(delta, digits = 3), "), ",
        "with recent changes being more important than early changes."
      )
    } else {
      interpretation <- paste0(
        "Chromosome evolution was concentrated early in the phylogeny (delta = ", 
        format(delta, digits = 3), "), ",
        "with less change in recent branches."
      )
    }
  } else if(model_name == "kappa") {
    # Pagel's kappa - punctuational vs. gradual evolution
    kappa <- extract_param(model_fit, "kappa")
    
    if(kappa < 0.5) {
      interpretation <- paste0(
        "Chromosome evolution shows punctuational pattern (kappa = ", 
        format(kappa, digits = 3), "), ",
        "with changes concentrated at speciation events rather than proportional to branch lengths."
      )
    } else {
      interpretation <- paste0(
        "Chromosome evolution follows a gradual pattern (kappa = ", 
        format(kappa, digits = 3), "), ",
        "with changes occurring proportionally to evolutionary time."
      )
    }
  }
  
  return(interpretation)
}

#' Extract parameter from model fit object
#' 
#' @param model_fit Fitted model object
#' @param param_name Name of parameter to extract
#' @return Parameter value or NA if not found
#' @keywords internal
extract_param <- function(model_fit, param_name) {
  if(is.null(model_fit)) {
    return(NA)
  }
  
  # Handle different model fit structures
  if(inherits(model_fit, "fitContinuous")) {
    # geiger model
    if(param_name %in% names(model_fit$opt)) {
      return(model_fit$opt[[param_name]])
    } else if(param_name %in% names(model_fit$opt$param)) {
      return(model_fit$opt$param[[param_name]])
    }
  } else if(inherits(model_fit, "phylosig")) {
    # phytools phylogenetic signal model
    if(param_name == "lambda" || param_name == "K" || param_name == "delta" || param_name == "kappa") {
      return(model_fit[[param_name]])
    }
  }
  
  return(NA)
}

#===============================================================================
# Disparity Analysis Functions
#===============================================================================

#' Analyze disparity in chromosome numbers across a phylogeny
#' 
#' Computes metrics of morphological disparity (diversity) in chromosome numbers
#' and tests for patterns in disparity across time and clades
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param method Disparity calculation method: "variance", "range", "quantile", "pairwise"
#' @param clades Optional list of clades for clade-specific disparity analysis
#' @param time_slices Number of time slices for disparity-through-time analysis
#' @param n_replicates Number of simulation replicates for significance testing
#' @param plot_results Whether to create plots of disparity patterns
#' @return List with disparity analysis results
#' @export
analyze_chromosome_disparity <- function(tree,
                                      chr_counts,
                                      method = "variance",
                                      clades = NULL,
                                      time_slices = 100,
                                      n_replicates = 1000,
                                      plot_results = TRUE) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Match chromosome counts to tree tips
  matched_counts <- match_taxon_data(tree, chr_counts)
  
  # Initialize results
  results <- list(
    overall_disparity = NA,
    method = method,
    summary = list(),
    plots = list()
  )
  
  # Calculate overall disparity based on method
  if(method == "variance") {
    results$overall_disparity <- var(matched_counts, na.rm = TRUE)
  } else if(method == "range") {
    range_val <- range(matched_counts, na.rm = TRUE)
    results$overall_disparity <- range_val[2] - range_val[1]
  } else if(method == "quantile") {
    q_range <- quantile(matched_counts, probs = c(0.25, 0.75), na.rm = TRUE)
    results$overall_disparity <- q_range[2] - q_range[1]  # Interquartile range
  } else if(method == "pairwise") {
    # Mean pairwise difference
    n <- length(matched_counts)
    if(n > 1) {
      all_pairs <- combn(n, 2)
      diffs <- abs(matched_counts[all_pairs[1,]] - matched_counts[all_pairs[2,]])
      results$overall_disparity <- mean(diffs, na.rm = TRUE)
    } else {
      results$overall_disparity <- 0
    }
  } else {
    stop(paste("Unsupported disparity method:", method))
  }
  
  # Disparity through time (DTT) analysis if geiger is available
  if(requireNamespace("geiger", quietly = TRUE)) {
    # Prepare data for geiger
    taxa_names <- tree$tip.label
    matched_matrix <- matrix(matched_counts, ncol = 1, 
                           dimnames = list(taxa_names, "chr_count"))
    
    # Run DTT analysis
    dtt_result <- tryCatch({
      geiger::dtt(tree, matched_matrix, nsim = n_replicates, 
                 calcMDIp = TRUE, plot = FALSE)
    }, error = function(e) {
      warning(paste("Error in DTT analysis:", e$message))
      return(NULL)
    })
    
    if(!is.null(dtt_result)) {
      results$dtt <- dtt_result
      
      # Add MDI statistic and p-value
      results$summary$mdi <- dtt_result$MDI
      results$summary$mdi_p <- dtt_result$MDIp
      
      # Interpret DTT results
      results$summary$dtt_pattern <- interpret_dtt_results(dtt_result)
      
      # Create DTT plot if requested
      if(plot_results) {
        results$plots$dtt <- create_dtt_plot(dtt_result)
      }
    }
  }
  
  # Analyze disparity within clades if provided
  if(!is.null(clades) && is.list(clades)) {
    clade_disparities <- list()
    
    for(clade_name in names(clades)) {
      clade_tips <- clades[[clade_name]]
      
      # Check if tips are in tree
      valid_tips <- intersect(clade_tips, tree$tip.label)
      
      if(length(valid_tips) >= 2) {
        # Extract chromosome counts for this clade
        clade_counts <- matched_counts[valid_tips]
        
        # Calculate disparity based on method
        if(method == "variance") {
          clade_disp <- var(clade_counts, na.rm = TRUE)
        } else if(method == "range") {
          clade_range <- range(clade_counts, na.rm = TRUE)
          clade_disp <- clade_range[2] - clade_range[1]
        } else if(method == "quantile") {
          clade_q_range <- quantile(clade_counts, probs = c(0.25, 0.75), na.rm = TRUE)
          clade_disp <- clade_q_range[2] - clade_q_range[1]
        } else if(method == "pairwise") {
          n_clade <- length(clade_counts)
          if(n_clade > 1) {
            clade_pairs <- combn(n_clade, 2)
            clade_diffs <- abs(clade_counts[clade_pairs[1,]] - clade_counts[clade_pairs[2,]])
            clade_disp <- mean(clade_diffs, na.rm = TRUE)
          } else {
            clade_disp <- 0
          }
        }
        
        # Store clade disparity
        clade_disparities[[clade_name]] <- list(
          disparity = clade_disp,
          variance = var(clade_counts, na.rm = TRUE),
          range = diff(range(clade_counts, na.rm = TRUE)),
          n_taxa = length(valid_tips)
        )
      }
    }
    
    results$clade_disparities <- clade_disparities
    
    # Test for significant differences between clades
    if(length(clade_disparities) >= 2) {
      # Extract disparities for each clade
      clade_disp_values <- sapply(clade_disparities, function(x) x$disparity)
      
      # Test for significance using randomization
      sig_results <- test_disparity_differences(tree, matched_counts, 
                                             clades, method, n_replicates)
      
      results$clade_comparisons <- sig_results$comparisons
      results$summary$n_significant_comparisons <- sig_results$n_significant
      results$summary$significant_comparisons <- sig_results$significant_pairs
      
      # Create clade comparison plot if requested
      if(plot_results) {
        results$plots$clade_comparison <- create_clade_disparity_plot(clade_disparities)
      }
    }
  }
  
  # Set up summary statistics
  results$summary$overall_disparity <- results$overall_disparity
  
  return(results)
}

#' Match taxon data to tree tips
#' 
#' @param tree Phylogenetic tree
#' @param data Named vector of data
#' @return Named vector of data matched to tree tips
#' @keywords internal
match_taxon_data <- function(tree, data) {
  # Match data to tree tips
  matched_data <- numeric(length(tree$tip.label))
  names(matched_data) <- tree$tip.label
  
  for(tip in tree$tip.label) {
    if(tip %in% names(data)) {
      matched_data[tip] <- data[tip]
    } else {
      matched_data[tip] <- NA
    }
  }
  
  return(matched_data)
}

#' Interpret disparity-through-time results
#' 
#' @param dtt_result Results from a DTT analysis
#' @return Textual interpretation of results
#' @keywords internal
interpret_dtt_results <- function(dtt_result) {
  # Extract MDI statistic and p-value
  mdi <- dtt_result$MDI
  p_value <- dtt_result$MDIp
  
  # Interpret MDI value
  if(is.na(mdi) || is.null(mdi)) {
    return("Could not calculate MDI statistic")
  }
  
  if(mdi > 0 && p_value <= 0.05) {
    return(paste0("Significant positive MDI (", round(mdi, 3), 
                ", p = ", round(p_value, 3), 
                "): Disparity is higher than expected within subclades, suggesting niche-filling or adaptive radiation"))
  } else if(mdi < 0 && p_value <= 0.05) {
    return(paste0("Significant negative MDI (", round(mdi, 3), 
                ", p = ", round(p_value, 3), 
                "): Disparity is lower than expected within subclades, suggesting niche conservatism"))
  } else {
    return(paste0("Non-significant MDI (", round(mdi, 3), 
                ", p = ", round(p_value, 3),
                "): Disparity through time follows expectations under Brownian motion evolution"))
  }
}

#' Test for significant differences in disparity between clades
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts
#' @param clades List of clades
#' @param method Disparity calculation method
#' @param n_replicates Number of randomization replicates
#' @return List with significance test results
#' @keywords internal
test_disparity_differences <- function(tree, chr_counts, clades, method, n_replicates) {
  clade_names <- names(clades)
  n_clades <- length(clade_names)
  
  # Calculate observed disparities
  observed_disp <- numeric(n_clades)
  names(observed_disp) <- clade_names
  
  for(i in 1:n_clades) {
    clade_tips <- intersect(clades[[i]], tree$tip.label)
    clade_counts <- chr_counts[clade_tips]
    
    if(length(clade_tips) >= 2) {
      # Calculate disparity based on method
      if(method == "variance") {
        observed_disp[i] <- var(clade_counts, na.rm = TRUE)
      } else if(method == "range") {
        observed_disp[i] <- diff(range(clade_counts, na.rm = TRUE))
      } else if(method == "quantile") {
        q_range <- quantile(clade_counts, probs = c(0.25, 0.75), na.rm = TRUE)
        observed_disp[i] <- q_range[2] - q_range[1]
      } else if(method == "pairwise") {
        n <- length(clade_counts)
        if(n > 1) {
          all_pairs <- combn(n, 2)
          diffs <- abs(clade_counts[all_pairs[1,]] - clade_counts[all_pairs[2,]])
          observed_disp[i] <- mean(diffs, na.rm = TRUE)
        } else {
          observed_disp[i] <- 0
        }
      }
    } else {
      observed_disp[i] <- NA
    }
  }
  
  # Create all pairwise comparisons
  pairs <- combn(clade_names, 2)
  n_pairs <- ncol(pairs)
  
  # Initialize results
  p_values <- numeric(n_pairs)
  
  # Perform randomization test for each pair
  for(j in 1:n_pairs) {
    clade1 <- pairs[1, j]
    clade2 <- pairs[2, j]
    
    # Skip if either clade has NA disparity
    if(is.na(observed_disp[clade1]) || is.na(observed_disp[clade2])) {
      p_values[j] <- NA
      next
    }
    
    # Get tips for both clades
    tips1 <- intersect(clades[[clade1]], tree$tip.label)
    tips2 <- intersect(clades[[clade2]], tree$tip.label)
    
    # Calculate observed disparity difference
    obs_diff <- abs(observed_disp[clade1] - observed_disp[clade2])
    
    # Combine tips
    all_tips <- c(tips1, tips2)
    all_counts <- chr_counts[all_tips]
    
    # Randomization test
    random_diffs <- numeric(n_replicates)
    
    for(rep in 1:n_replicates) {
      # Randomly shuffle tips between clades
      shuffled <- sample(all_tips)
      random_tips1 <- shuffled[1:length(tips1)]
      random_tips2 <- shuffled[(length(tips1)+1):length(all_tips)]
      
      random_counts1 <- chr_counts[random_tips1]
      random_counts2 <- chr_counts[random_tips2]
      
      # Calculate disparities based on method
      if(method == "variance") {
        disp1 <- var(random_counts1, na.rm = TRUE)
        disp2 <- var(random_counts2, na.rm = TRUE)
      } else if(method == "range") {
        disp1 <- diff(range(random_counts1, na.rm = TRUE))
        disp2 <- diff(range(random_counts2, na.rm = TRUE))
      } else if(method == "quantile") {
        q_range1 <- quantile(random_counts1, probs = c(0.25, 0.75), na.rm = TRUE)
        q_range2 <- quantile(random_counts2, probs = c(0.25, 0.75), na.rm = TRUE)
        disp1 <- q_range1[2] - q_range1[1]
        disp2 <- q_range2[2] - q_range2[1]
      } else if(method == "pairwise") {
        n1 <- length(random_counts1)
        n2 <- length(random_counts2)
        
        if(n1 > 1) {
          pairs1 <- combn(n1, 2)
          diffs1 <- abs(random_counts1[pairs1[1,]] - random_counts1[pairs1[2,]])
          disp1 <- mean(diffs1, na.rm = TRUE)
        } else {
          disp1 <- 0
        }
        
        if(n2 > 1) {
          pairs2 <- combn(n2, 2)
          diffs2 <- abs(random_counts2[pairs2[1,]] - random_counts2[pairs2[2,]])
          disp2 <- mean(diffs2, na.rm = TRUE)
        } else {
          disp2 <- 0
        }
      }
      
      # Store randomized difference
      random_diffs[rep] <- abs(disp1 - disp2)
    }
    
    # Calculate p-value
    p_values[j] <- sum(random_diffs >= obs_diff) / n_replicates
  }
  
  # Create comparison results
  comparison_list <- list()
  
  for(j in 1:n_pairs) {
    clade1 <- pairs[1, j]
    clade2 <- pairs[2, j]
    
    comparison_list[[paste(clade1, "vs", clade2)]] <- list(
      disp1 = observed_disp[clade1],
      disp2 = observed_disp[clade2],
      diff = abs(observed_disp[clade1] - observed_disp[clade2]),
      p_value = p_values[j],
      significant = !is.na(p_values[j]) && p_values[j] <= 0.05
    )
  }
  
  # Find significant comparisons
  significant_indices <- which(p_values <= 0.05 & !is.na(p_values))
  n_significant <- length(significant_indices)
  
  significant_pairs <- character(0)
  if(n_significant > 0) {
    for(idx in significant_indices) {
      significant_pairs <- c(significant_pairs,
                           paste(pairs[1, idx], "vs", pairs[2, idx]))
    }
  }
  
  return(list(
    comparisons = comparison_list,
    n_significant = n_significant,
    significant_pairs = significant_pairs,
    p_values = p_values
  ))
}

#' Create disparity-through-time plot
#' 
#' @param dtt_result Results from a DTT analysis
#' @return ggplot object with DTT plot
#' @keywords internal
create_dtt_plot <- function(dtt_result) {
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package required for plotting")
    return(NULL)
  }
  
  # Extract DTT data
  times <- dtt_result$times
  dtt <- dtt_result$dtt
  sim <- dtt_result$sim
  mdi <- dtt_result$MDI
  mdip <- dtt_result$MDIp
  
  # Create data frames for plotting
  dtt_data <- data.frame(
    Time = times,
    DTT = dtt,
    stringsAsFactors = FALSE
  )
  
  sim_data <- data.frame()
  
  if(!is.null(sim)) {
    sim_means <- apply(sim, 1, mean)
    sim_lower <- apply(sim, 1, quantile, probs = 0.025)
    sim_upper <- apply(sim, 1, quantile, probs = 0.975)
    
    sim_data <- data.frame(
      Time = times,
      Mean = sim_means,
      Lower = sim_lower,
      Upper = sim_upper,
      stringsAsFactors = FALSE
    )
  }
  
  # Create the plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = dtt_data, 
                      ggplot2::aes(x = Time, y = DTT), 
                      color = "blue", size = 1) +
    ggplot2::labs(
      title = "Disparity Through Time (DTT)",
      subtitle = paste0("MDI = ", round(mdi, 3), 
                      ifelse(!is.null(mdip), paste0(", p = ", round(mdip, 3)), "")),
      x = "Relative Time",
      y = "Relative Disparity"
    ) +
    ggplot2::theme_minimal()
  
  # Add simulation envelope if available
  if(nrow(sim_data) > 0) {
    p <- p +
      ggplot2::geom_ribbon(data = sim_data, 
                         ggplot2::aes(x = Time, ymin = Lower, ymax = Upper),
                         fill = "gray80", alpha = 0.5) +
      ggplot2::geom_line(data = sim_data,
                       ggplot2::aes(x = Time, y = Mean),
                       color = "black", linetype = "dashed")
  }
  
  return(p)
}

#' Create clade disparity comparison plot
#' 
#' @param clade_disparities List of disparity values for each clade
#' @return ggplot object with clade comparison plot
#' @keywords internal
create_clade_disparity_plot <- function(clade_disparities) {
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package required for plotting")
    return(NULL)
  }
  
  # Extract data for plotting
  clade_names <- names(clade_disparities)
  disp_values <- sapply(clade_disparities, function(x) x$disparity)
  
  # Create data frame
  plot_data <- data.frame(
    Clade = clade_names,
    Disparity = disp_values,
    stringsAsFactors = FALSE
  )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, 
                      ggplot2::aes(x = reorder(Clade, -Disparity), y = Disparity)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(
      title = "Chromosome Number Disparity by Clade",
      x = "Clade",
      y = "Disparity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}

#===============================================================================
# Ancestral State Hypothesis Testing Functions
#===============================================================================

#' Test hypotheses about ancestral chromosome numbers
#' 
#' Tests specific hypotheses about chromosome numbers at nodes in the phylogeny
#' 
#' @param tree Phylogenetic tree
#' @param reconstruction Ancestral state reconstruction results
#' @param hypotheses List of hypotheses to test, each a list with node_id and expected_state
#' @param test_type Type of test: "direct" or "simulation"
#' @param alpha Significance level
#' @param n_simulations Number of simulations for simulation-based tests
#' @param confidence_method Confidence interval method: "ml", "bayesian", or "bootstrap"
#' @param confidence_level Confidence level for intervals (e.g., 0.95)
#' @return List with hypothesis test results
#' @export
test_ancestral_chromosome_hypotheses <- function(tree,
                                               reconstruction,
                                               hypotheses,
                                               test_type = "direct",
                                               alpha = 0.05,
                                               n_simulations = 1000,
                                               confidence_method = "ml",
                                               confidence_level = 0.95) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(is.null(hypotheses) || !is.list(hypotheses)) {
    stop("hypotheses must be a list of hypothesis specifications")
  }
  
  # Extract ancestral states from reconstruction
  if(is.list(reconstruction)) {
    if(!is.null(reconstruction$ancestral_states)) {
      states <- reconstruction$ancestral_states
    } else if(!is.null(reconstruction$states)) {
      states <- reconstruction$states
    } else {
      stop("Could not find ancestral states in reconstruction object")
    }
  } else if(is.vector(reconstruction) && !is.null(names(reconstruction))) {
    states <- reconstruction
  } else {
    stop("Invalid reconstruction format")
  }
  
  # Convert to named vector if necessary
  if(is.data.frame(states)) {
    node_ids <- states$node_id
    state_values <- states$state
    states <- state_values
    names(states) <- node_ids
  }
  
  # Extract confidence intervals if available
  conf_intervals <- NULL
  if(is.list(reconstruction) && !is.null(reconstruction$confidence_intervals)) {
    conf_intervals <- reconstruction$confidence_intervals
  }
  
  # Check if state distributions are available (for Bayesian reconstructions)
  state_distributions <- NULL
  if(is.list(reconstruction) && !is.null(reconstruction$state_distributions)) {
    state_distributions <- reconstruction$state_distributions
  }
  
  # Initialize results
  results <- list(
    hypotheses = hypotheses,
    tests = list(),
    summary = list(
      n_tested = length(hypotheses),
      n_rejected = 0,
      n_accepted = 0,
      significance_level = alpha,
      test_type = test_type,
      confidence_method = confidence_method
    )
  )
  
  # Test each hypothesis
  for(i in seq_along(hypotheses)) {
    hypothesis <- hypotheses[[i]]
    
    # Check if hypothesis has required fields
    if(!is.list(hypothesis) || 
       !all(c("node_id", "expected_state") %in% names(hypothesis))) {
      warning(paste("Hypothesis", i, "is missing required fields. Skipping."))
      next
    }
    
    node_id <- as.character(hypothesis$node_id)
    expected_state <- hypothesis$expected_state
    
    # Check if node exists in the reconstruction
    if(!node_id %in% names(states)) {
      warning(paste("Node", node_id, "not found in reconstruction. Skipping hypothesis", i))
      next
    }
    
    # Get reconstructed state
    reconstructed_state <- states[node_id]
    
    # Initialize test result
    test_result <- list(
      node_id = node_id,
      expected_state = expected_state,
      reconstructed_state = reconstructed_state,
      difference = reconstructed_state - expected_state,
      rejected = FALSE,
      p_value = NA,
      description = hypothesis$description
    )
    
    # Perform test based on test_type
    if(test_type == "direct") {
      # Direct test using confidence intervals
      if(!is.null(conf_intervals)) {
        lower <- conf_intervals$lower[node_id]
        upper <- conf_intervals$upper[node_id]
        
        test_result$conf_lower <- lower
        test_result$conf_upper <- upper
        
        # Reject if expected value falls outside confidence interval
        test_result$rejected <- expected_state < lower || expected_state > upper
        
        # Calculate p-value based on confidence level
        # For ML estimates, assuming normal distribution of error
        # Compute Z-score based on estimated standard error
        if(!is.null(conf_intervals$se) && !is.na(conf_intervals$se[node_id])) {
          se <- conf_intervals$se[node_id]
          z_score <- abs(reconstructed_state - expected_state) / se
          test_result$p_value <- 2 * pnorm(-abs(z_score))  # Two-tailed test
        }
      } else if(!is.null(state_distributions)) {
        # For Bayesian reconstruction, use posterior distribution
        if(is.list(state_distributions) && !is.null(state_distributions[[node_id]])) {
          posterior <- state_distributions[[node_id]]
          
          # Calculate credible interval
          cred_interval <- quantile(posterior, probs = c((1-confidence_level)/2, 
                                                       1-(1-confidence_level)/2))
          test_result$conf_lower <- cred_interval[1]
          test_result$conf_upper <- cred_interval[2]
          
          # Reject if expected value falls outside credible interval
          test_result$rejected <- expected_state < cred_interval[1] || 
                                expected_state > cred_interval[2]
          
          # Calculate p-value as proportion of posterior outside expected value
          if(expected_state <= median(posterior)) {
            test_result$p_value <- 2 * mean(posterior <= expected_state)
          } else {
            test_result$p_value <- 2 * mean(posterior >= expected_state)
          }
        }
      } else {
        warning("No confidence intervals or state distributions available for direct test")
      }
    } else if(test_type == "simulation") {
      # Simulation-based test
      if(requireNamespace("phytools", quietly = TRUE)) {
        # Extract model from reconstruction if available
        model <- NULL
        if(is.list(reconstruction) && !is.null(reconstruction$model)) {
          model <- reconstruction$model
        } else {
          # Default to Brownian motion
          model <- "BM"
        }
        
        # Extract model parameters
        params <- list()
        if(is.list(reconstruction) && !is.null(reconstruction$model_parameters)) {
          params <- reconstruction$model_parameters
        } else {
          # Estimate sigma^2 from the data
          if(model == "BM") {
            chr_tips <- states[1:length(tree$tip.label)]
            chr_tips <- chr_tips[!is.na(chr_tips)]
            
            if(length(chr_tips) >= 4) {
              # Estimate rate parameter from tip data
              pic_values <- ape::pic(chr_tips, tree)
              params$sigma2 <- mean(pic_values^2)
            }
          }
        }
        
        # Simulate under the model
        root_state <- NA
        if(is.list(reconstruction) && !is.null(reconstruction$root_state)) {
          root_state <- reconstruction$root_state
        } else {
          root_node <- length(tree$tip.label) + 1  # Standard root node numbering
          if(as.character(root_node) %in% names(states)) {
            root_state <- states[as.character(root_node)]
          } else {
            # Use mean of tip states as fallback
            root_state <- mean(states[1:length(tree$tip.label)], na.rm = TRUE)
          }
        }
        
        # Perform simulations
        simulated_states <- matrix(NA, nrow = n_simulations, 
                                 ncol = tree$Nnode + length(tree$tip.label))
        
        for(sim in 1:n_simulations) {
          if(model == "BM") {
            # Simulate under Brownian motion
            sim_result <- phytools::fastBM(tree, 
                                         sig2 = params$sigma2, 
                                         a = root_state,
                                         internal = TRUE)
          } else if(model == "OU") {
            # Simulate under OU process
            sim_result <- phytools::fastBM(tree, 
                                         sig2 = params$sigma2, 
                                         alpha = params$alpha,
                                         a = root_state,
                                         theta = params$theta,
                                         internal = TRUE)
          } else {
            warning("Unsupported model for simulation")
            sim_result <- rep(NA, tree$Nnode + length(tree$tip.label))
          }
          
          simulated_states[sim, ] <- sim_result
        }
        
        # Extract simulated values for the node of interest
        node_index <- as.integer(node_id)
        if(node_index <= length(tree$tip.label)) {
          sim_node_states <- simulated_states[, node_index]
        } else {
          # Adjust index for internal nodes
          internal_idx <- node_index - length(tree$tip.label)
          sim_node_states <- simulated_states[, length(tree$tip.label) + internal_idx]
        }
        
        # Calculate p-value from simulations
        if(expected_state <= reconstructed_state) {
          test_result$p_value <- mean(sim_node_states <= expected_state)
        } else {
          test_result$p_value <- mean(sim_node_states >= expected_state)
        }
        
        # Calculate confidence interval from simulations
        test_result$conf_lower <- quantile(sim_node_states, probs = (1-confidence_level)/2)
        test_result$conf_upper <- quantile(sim_node_states, probs = 1-(1-confidence_level)/2)
        
        # Determine rejection
        test_result$rejected <- test_result$p_value <= alpha
      } else {
        warning("phytools package required for simulation-based tests")
      }
    } else {
      warning(paste("Unsupported test type:", test_type))
    }
    
    # Store test result
    results$tests[[i]] <- test_result
  }
  
  # Update summary statistics
  results$summary$n_rejected <- sum(sapply(results$tests, function(x) isTRUE(x$rejected)))
  results$summary$n_accepted <- length(results$tests) - results$summary$n_rejected
  
  return(results)
}

#===============================================================================
# Workflow Functions
#===============================================================================

#' Run a complete comparative analysis workflow for chromosome evolution
#' 
#' Performs a series of comparative analyses to test evolutionary patterns
#' and hypotheses about chromosome number evolution
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param traits Optional data frame of traits for correlation analysis
#' @param clades Optional list of clades for clade comparison
#' @param output_dir Directory to save results
#' @param analyses Which analyses to run: "models", "traits", "disparity", "clades", "all"
#' @param generate_report Whether to generate an HTML report of results
#' @param plot_results Whether to create plots of results
#' @return List with results from all analyses
#' @export
run_comparative_analysis <- function(tree,
                                   chr_counts,
                                   traits = NULL,
                                   clades = NULL,
                                   output_dir = "comparative_results",
                                   analyses = "all",
                                   generate_report = TRUE,
                                   plot_results = TRUE) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Determine which analyses to run
  available_analyses <- c("models", "traits", "disparity", "clades")
  if(length(analyses) == 1 && analyses == "all") {
    analyses_to_run <- available_analyses
  } else {
    analyses_to_run <- match.arg(analyses, available_analyses, several.ok = TRUE)
  }
  
  # Initialize results
  results <- list(
    tree = tree,
    chr_counts = chr_counts,
    analyses_run = analyses_to_run,
    timestamp = Sys.time()
  )
  
  # Run evolutionary models analysis if requested
  if("models" %in% analyses_to_run) {
    message("Analyzing evolutionary models...")
    model_results <- test_chromosome_evolution_models(
      tree = tree,
      chr_counts = chr_counts,
      models = c("BM", "OU", "EB", "trend", "lambda", "kappa", "delta", "white"),
      discrete = FALSE
    )
    
    results$model_analysis <- model_results
    
    # Save model results
    saveRDS(model_results, file.path(output_dir, "model_results.rds"))
    
    # Create model comparison plot if requested
    if(plot_results && !is.null(model_results$model_comparison) && 
       requireNamespace("ggplot2", quietly = TRUE)) {
      model_plot <- create_model_comparison_plot(model_results)
      results$plots$model_comparison <- model_plot
      
      if(!is.null(model_plot)) {
        ggplot2::ggsave(
          file.path(output_dir, "model_comparison.pdf"),
          plot = model_plot,
          width = 8,
          height = 6
        )
      }
    }
  }
  
  # Run trait correlation analysis if requested and traits provided
  if("traits" %in% analyses_to_run && !is.null(traits)) {
    message("Analyzing trait correlations...")
    trait_results <- analyze_trait_correlations(
      tree = tree,
      chr_counts = chr_counts,
      traits = traits,
      method = "PGLS"
    )
    
    results$trait_analysis <- trait_results
    
    # Save trait results
    saveRDS(trait_results, file.path(output_dir, "trait_results.rds"))
    
    # Create correlation plot if requested
    if(plot_results && trait_results$summary$n_significant > 0 && 
       requireNamespace("ggplot2", quietly = TRUE)) {
      corr_plot <- create_trait_correlation_plot(trait_results)
      results$plots$trait_correlation <- corr_plot
      
      if(!is.null(corr_plot)) {
        ggplot2::ggsave(
          file.path(output_dir, "trait_correlations.pdf"),
          plot = corr_plot,
          width = 8,
          height = 6
        )
      }
    }
  }
  
  # Run disparity analysis if requested
  if("disparity" %in% analyses_to_run) {
    message("Analyzing chromosome number disparity...")
    disparity_results <- analyze_chromosome_disparity(
      tree = tree,
      chr_counts = chr_counts,
      plot_results = FALSE  # We'll handle plotting separately
    )
    
    results$disparity_analysis <- disparity_results
    
    # Save disparity results
    saveRDS(disparity_results, file.path(output_dir, "disparity_results.rds"))
    
    # Create disparity plots if requested
    if(plot_results && requireNamespace("ggplot2", quietly = TRUE)) {
      # Save DTT plot if available
      if(!is.null(disparity_results$plots$dtt)) {
        ggplot2::ggsave(
          file.path(output_dir, "disparity_through_time.pdf"),
          plot = disparity_results$plots$dtt,
          width = 8,
          height = 6
        )
      }
      
      # Save clade comparison plot if available
      if(!is.null(disparity_results$plots$clade_comparison)) {
        ggplot2::ggsave(
          file.path(output_dir, "clade_disparity.pdf"),
          plot = disparity_results$plots$clade_comparison,
          width = 8,
          height = 6
        )
      }
    }
  }
  
  # Run clade comparison analysis if requested and clades provided
  if("clades" %in% analyses_to_run && !is.null(clades)) {
    message("Comparing chromosome evolution across clades...")
    clade_results <- compare_chromosome_evolution(
      tree = tree,
      chr_counts = chr_counts,
      clades = clades,
      metrics = c("mean", "variance", "range")
    )
    
    results$clade_analysis <- clade_results
    
    # Save clade results
    saveRDS(clade_results, file.path(output_dir, "clade_results.rds"))
    
    # Create clade comparison plots if requested
    if(plot_results && requireNamespace("ggplot2", quietly = TRUE)) {
      # Save mean comparison plot if available
      if(!is.null(clade_results$plots$mean)) {
        ggplot2::ggsave(
          file.path(output_dir, "clade_mean_comparison.pdf"),
          plot = clade_results$plots$mean,
          width = 8,
          height = 6
        )
      }
      
      # Save variance comparison plot if available
      if(!is.null(clade_results$plots$variance)) {
        ggplot2::ggsave(
          file.path(output_dir, "clade_variance_comparison.pdf"),
          plot = clade_results$plots$variance,
          width = 8,
          height = 6
        )
      }
      
      # Save range comparison plot if available
      if(!is.null(clade_results$plots$range)) {
        ggplot2::ggsave(
          file.path(output_dir, "clade_range_comparison.pdf"),
          plot = clade_results$plots$range,
          width = 8,
          height = 6
        )
      }
    }
  }
  
  # Generate summary of all analyses
  results$summary <- generate_comparative_summary(results)
  
  # Write summary to text file
  writeLines(
    results$summary$text,
    file.path(output_dir, "comparative_analysis_summary.txt")
  )
  
  # Generate HTML report if requested
  if(generate_report && requireNamespace("rmarkdown", quietly = TRUE)) {
    message("Generating HTML report...")
    report_template <- generate_report_template(results)
    
    # Write temporary Rmd file
    report_path <- file.path(output_dir, "temp_report.Rmd")
    writeLines(report_template, report_path)
    
    # Render report
    rmarkdown::render(
      report_path,
      output_file = "comparative_analysis_report.html",
      output_dir = output_dir,
      quiet = TRUE
    )
    
    # Remove temporary Rmd file
    file.remove(report_path)
  }
  
  message("Comparative analysis complete. Results saved to: ", output_dir)
  
  return(results)
}

#' Create model comparison plot
#' 
#' @param model_results Results from test_chromosome_evolution_models
#' @return ggplot object
#' @keywords internal
create_model_comparison_plot <- function(model_results) {
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package required for plotting")
    return(NULL)
  }
  
  # Check if model comparison is available
  if(is.null(model_results$model_comparison) || 
     nrow(model_results$model_comparison) == 0) {
    warning("No model comparison results available")
    return(NULL)
  }
  
  comp <- model_results$model_comparison
  
  # Get AIC weights
  criterion <- model_results$criterion
  weight_col <- paste0(criterion, "_weight")
  
  if(!weight_col %in% colnames(comp)) {
    warning("Model weights not found in comparison results")
    return(NULL)
  }
  
  # Create plot
  p <- ggplot2::ggplot(comp, ggplot2::aes(x = reorder(Model, -get(weight_col)), 
                                        y = .data[[weight_col]])) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(
      title = paste("Model Comparison (", criterion, ")", sep = ""),
      subtitle = paste("Best model:", model_results$best_model),
      x = "Model",
      y = paste(criterion, "weight")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Create trait correlation plot
#' 
#' @param trait_results Results from analyze_trait_correlations
#' @return ggplot object
#' @keywords internal
create_trait_correlation_plot <- function(trait_results) {
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package required for plotting")
    return(NULL)
  }
  
  # Check if there are significant correlations
  if(trait_results$summary$n_significant == 0) {
    warning("No significant trait correlations found")
    return(NULL)
  }
  
  # Create data frame of significant correlations
  sig_traits <- trait_results$summary$significant_traits
  
  corr_data <- data.frame(
    Trait = character(0),
    Estimate = numeric(0),
    P_Value = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for(trait in sig_traits) {
    corr_data <- rbind(corr_data, data.frame(
      Trait = trait,
      Estimate = trait_results$correlations[[trait]]$estimate,
      P_Value = trait_results$correlations[[trait]]$p_value,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create plot
  p <- ggplot2::ggplot(corr_data, 
                      ggplot2::aes(x = reorder(Trait, -abs(Estimate)), y = Estimate)) +
    ggplot2::geom_bar(stat = "identity", 
                    fill = ifelse(corr_data$Estimate > 0, "darkblue", "darkred")) +
    ggplot2::labs(
      title = "Significant Trait Correlations with Chromosome Number",
      subtitle = paste("Method:", trait_results$summary$method),
      x = "Trait",
      y = "Correlation Coefficient"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' Generate summary of comparative analysis results
#' 
#' @param results Results from run_comparative_analysis
#' @return List with summary information
#' @keywords internal
generate_comparative_summary <- function(results) {
  # Initialize text summary
  summary_text <- c(
    "COMPARATIVE ANALYSIS OF CHROMOSOME EVOLUTION",
    "===========================================",
    "",
    paste("Analysis date:", format(results$timestamp, "%Y-%m-%d %H:%M:%S")),
    paste("Analyses performed:", paste(results$analyses_run, collapse = ", ")),
    paste("Number of taxa in tree:", length(results$tree$tip.label)),
    paste("Number of taxa with chromosome data:", sum(!is.na(results$chr_counts))),
    "",
    "KEY FINDINGS:",
    "--------------",
    ""
  )
  
  # Add model analysis summary
  if(!is.null(results$model_analysis)) {
    if(!is.null(results$model_analysis$summary)) {
      summary_text <- c(summary_text,
                      "Evolutionary Models:",
                      paste("  Best model:", results$model_analysis$summary$best_model),
                      paste("  Support strength:", results$model_analysis$summary$support_strength),
                      paste("  Interpretation:", results$model_analysis$summary$model_interpretation),
                      "")
    }
  }
  
  # Add trait correlation summary
  if(!is.null(results$trait_analysis)) {
    if(!is.null(results$trait_analysis$summary)) {
      summary_text <- c(summary_text,
                      "Trait Correlations:",
                      paste("  Number of traits tested:", results$trait_analysis$summary$n_traits),
                      paste("  Number of significant correlations:", 
                          results$trait_analysis$summary$n_significant))
      
      if(results$trait_analysis$summary$n_significant > 0) {
        summary_text <- c(summary_text,
                        "  Significant traits:")
        
        for(trait in results$trait_analysis$summary$significant_traits) {
          corr <- results$trait_analysis$correlations[[trait]]
          summary_text <- c(summary_text,
                          paste("    -", trait, ":", 
                              "estimate =", round(corr$estimate, 3),
                              ", p-value =", round(corr$p_value, 3)))
        }
      }
      
      summary_text <- c(summary_text, "")
    }
  }
  
  # Add disparity analysis summary
  if(!is.null(results$disparity_analysis)) {
    if(!is.null(results$disparity_analysis$summary)) {
      summary_text <- c(summary_text,
                      "Disparity Analysis:",
                      paste("  Overall disparity:", 
                          round(results$disparity_analysis$overall_disparity, 3)))
      
      if(!is.null(results$disparity_analysis$summary$dtt_pattern)) {
        summary_text <- c(summary_text,
                        paste("  Disparity through time:", 
                            results$disparity_analysis$summary$dtt_pattern))
      }
      
      if(!is.null(results$disparity_analysis$summary$n_significant_comparisons) &&
         results$disparity_analysis$summary$n_significant_comparisons > 0) {
        summary_text <- c(summary_text,
                        paste("  Significant clade disparity differences:", 
                            results$disparity_analysis$summary$n_significant_comparisons))
      }
      
      summary_text <- c(summary_text, "")
    }
  }
  
  # Add clade comparison summary
  if(!is.null(results$clade_analysis)) {
    if(!is.null(results$clade_analysis$summary)) {
      summary_text <- c(summary_text,
                      "Clade Comparisons:",
                      paste("  Number of clades compared:", 
                          results$clade_analysis$summary$n_clades))
      
      if(results$clade_analysis$summary$has_significant_differences) {
        summary_text <- c(summary_text,
                        "  Significant differences detected between clades")
        
        # Add specific metric differences
        for(metric in names(results$clade_analysis$statistics)) {
          metric_result <- results$clade_analysis$statistics[[metric]]
          if(!is.null(metric_result$significant) && metric_result$significant) {
            summary_text <- c(summary_text,
                            paste("  -", metric, "differs significantly (p =", 
                                round(metric_result$p_value, 3), ")"))
          }
        }
      } else {
        summary_text <- c(summary_text,
                        "  No significant differences detected between clades")
      }
      
      summary_text <- c(summary_text, "")
    }
  }
  
  # Return summary as list
  return(list(
    text = summary_text,
    key_findings = summary_text
  ))
}

#' Generate R Markdown report template
#' 
#' @param results Results from run_comparative_analysis
#' @return Character vector with R Markdown content
#' @keywords internal
generate_report_template <- function(results) {
  # Create report template
  report <- c(
    "---",
    "title: \"Comparative Analysis of Chromosome Evolution\"",
    paste0("date: \"", format(Sys.time(), "%Y-%m-%d"), "\""),
    "output: html_document",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = TRUE)",
    "```",
    "",
    "# Summary",
    "",
    "```{r}",
    "cat(paste(results$summary$text, collapse = '\\n'))",
    "```",
    "",
    "# Evolutionary Models",
    "",
    "```{r}",
    "if(!is.null(results$model_analysis)) {",
    "  cat('## Best Model\\n')",
    "  cat(paste('Best model:', results$model_analysis$summary$best_model, '\\n'))",
    "  cat(paste('Support strength:', results$model_analysis$summary$support_strength, '\\n'))",
    "  cat(paste('Interpretation:', results$model_analysis$summary$model_interpretation, '\\n'))",
    "  if(!is.null(results$plots$model_comparison)) {",
    "    print(results$plots$model_comparison)",
    "  }",
    "}",
    "```",
    "",
    "# Trait Correlations",
    "",
    "```{r}",
    "if(!is.null(results$trait_analysis)) {",
    "  cat('## Significant Trait Correlations\\n')",
    "  cat(paste('Number of traits tested:', results$trait_analysis$summary$n_traits, '\\n'))",
    "  cat(paste('Number of significant correlations:', results$trait_analysis$summary$n_significant, '\\n'))",
    "  if(results$trait_analysis$summary$n_significant > 0) {",
    "    cat('### Significant Traits\\n')",
    "    for(trait in results$trait_analysis$summary$significant_traits) {",
    "      corr <- results$trait_analysis$correlations[[trait]]",
    "      cat(paste('- ', trait, ': estimate =', round(corr$estimate, 3), ', p-value =', round(corr$p_value, 3), '\\n'))",
    "    }",
    "    if(!is.null(results$plots$trait_correlation)) {",
    "      print(results$plots$trait_correlation)",
    "    }",
    "  }",
    "}",
    "```",
    "",
    "# Disparity Analysis",
    "",
    "```{r}",
    "if(!is.null(results$disparity_analysis)) {",
    "  cat('## Overall Disparity\\n')",
    "  cat(paste('Overall disparity:', round(results$disparity_analysis$overall_disparity, 3), '\\n'))",
    "  if(!is.null(results$disparity_analysis$summary$dtt_pattern)) {",
    "    cat(paste('Disparity through time:', results$disparity_analysis$summary$dtt_pattern, '\\n'))",
    "  }",
    "  if(!is.null(results$disparity_analysis$summary$n_significant_comparisons) && results$disparity_analysis$summary$n_significant_comparisons > 0) {",
    "    cat(paste('Significant clade disparity differences:', results$disparity_analysis$summary$n_significant_comparisons, '\\n'))",
    "  }",
    "  if(!is.null(results$plots$dtt)) {",
    "    print(results$plots$dtt)",
    "  }",
    "  if(!is.null(results$plots$clade_comparison)) {",
    "    print(results$plots$clade_comparison)",
    "  }",
    "}",
    "```",
    "",
    "# Clade Comparisons",
    "",
    "```{r}",
    "if(!is.null(results$clade_analysis)) {",
    "  cat('## Clade Comparisons\\n')",
    "  cat(paste('Number of clades compared:', results$clade_analysis$summary$n_clades, '\\n'))",
    "  if(results$clade_analysis$summary$has_significant_differences) {",
    "    cat('Significant differences detected between clades\\n')",
    "    for(metric in names(results$clade_analysis$statistics)) {",
    "      metric_result <- results$clade_analysis$statistics[[metric]]",
    "      if(!is.null(metric_result$significant) && metric_result$significant) {",
    "        cat(paste('- ', metric, 'differs significantly (p =', round(metric_result$p_value, 3), ')\\n'))",
    "      }",
    "    }",
    "  } else {",
    "    cat('No significant differences detected between clades\\n')",
    "  }",
    "  if(!is.null(results$plots$mean)) {",
    "    print(results$plots$mean)",
    "  }",
    "  if(!is.null(results$plots$variance)) {",
    "    print(results$plots$variance)",
    "  }",
    "  if(!is.null(results$plots$range)) {",
    "    print(results$plots$range)",
    "  }",
    "}",
    "```"
  )
  
  return(report)
}

