#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Genome Integration Module
# Author: Bioinformatics Team
# Date: 2025-05-30
# Description: Methods for integrating genomic data with chromosome evolution 
#              analyses, including correlations between chromosome changes and 
#              genomic features, synteny preservation, and sequence-based insights
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(geiger)
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})

#===============================================================================
# Data Integration Functions
#===============================================================================

#' Integrate chromosome counts with genomic feature data
#' 
#' Combines chromosome count data with genomic features like genome size,
#' gene counts, repeat content, etc. to enable integrated analyses
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param genome_features Data frame of genomic features with row names matching species names
#' @param feature_types Optional vector specifying feature data types (numeric, categorical)
#' @param check_tree_matching Whether to check and prune data to match tree tips
#' @param impute_missing Method for imputing missing data: "none", "mean", "phylo" (phylogenetic imputation)
#' @param scale_features Whether to scale numeric features (mean = 0, sd = 1)
#' @return Integrated data object with matched chromosome and genomic data
#' @export
integrate_genomic_data <- function(tree, 
                                 chr_counts, 
                                 genome_features, 
                                 feature_types = NULL,
                                 check_tree_matching = TRUE,
                                 impute_missing = "none",
                                 scale_features = FALSE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  if(!is.data.frame(genome_features)) {
    stop("genome_features must be a data frame")
  }
  
  # Get species information
  tree_taxa <- tree$tip.label
  chr_taxa <- names(chr_counts)
  
  # Check if genome_features has row names
  if(is.null(rownames(genome_features))) {
    stop("genome_features must have row names corresponding to species names")
  }
  
  feature_taxa <- rownames(genome_features)
  
  # Find taxa present in all datasets
  common_taxa <- Reduce(intersect, list(tree_taxa, chr_taxa, feature_taxa))
  
  if(length(common_taxa) < 3) {
    stop("At least 3 taxa must be present in the tree, chromosome data, and genome features")
  }
  
  # Match data to common taxa
  if(check_tree_matching) {
    # Prune tree to match common taxa
    pruned_tree <- ape::keep.tip(tree, common_taxa)
    
    # Subset chromosome counts
    pruned_chr <- chr_counts[common_taxa]
    
    # Subset genome features
    pruned_features <- genome_features[common_taxa, , drop = FALSE]
  } else {
    # Keep original data but warn about mismatches
    pruned_tree <- tree
    pruned_chr <- chr_counts
    pruned_features <- genome_features
    
    tree_only <- setdiff(tree_taxa, union(chr_taxa, feature_taxa))
    chr_only <- setdiff(chr_taxa, union(tree_taxa, feature_taxa))
    feature_only <- setdiff(feature_taxa, union(tree_taxa, chr_taxa))
    
    if(length(tree_only) > 0 || length(chr_only) > 0 || length(feature_only) > 0) {
      warning(paste("Data mismatch:", 
                    length(tree_only), "taxa in tree only,",
                    length(chr_only), "taxa in chromosome data only,",
                    length(feature_only), "taxa in genome features only"))
    }
  }
  
  # Determine feature types if not provided
  if(is.null(feature_types)) {
    feature_types <- sapply(pruned_features, function(x) {
      if(is.numeric(x)) {
        return("numeric")
      } else if(is.factor(x) || is.character(x)) {
        return("categorical")
      } else if(is.logical(x)) {
        return("binary")
      } else {
        return("unknown")
      }
    })
  }
  
  # Handle missing values through imputation if requested
  if(impute_missing != "none") {
    pruned_features <- impute_missing_data(pruned_tree, pruned_features, 
                                         feature_types, method = impute_missing)
  }
  
  # Scale numeric features if requested
  if(scale_features) {
    for(col in colnames(pruned_features)) {
      if(feature_types[col] == "numeric") {
        pruned_features[, col] <- scale(pruned_features[, col])
      }
    }
  }
  
  # Create integrated data object
  integrated_data <- list(
    tree = pruned_tree,
    chr_counts = pruned_chr,
    genome_features = pruned_features,
    feature_types = feature_types,
    taxa = common_taxa,
    n_taxa = length(common_taxa),
    summary = list(
      n_features = ncol(pruned_features),
      feature_names = colnames(pruned_features)
    )
  )
  
  class(integrated_data) <- c("acr_integrated", class(integrated_data))
  
  return(integrated_data)
}

#' Impute missing values in genomic feature data
#' 
#' @param tree Phylogenetic tree
#' @param features Data frame of features with missing values
#' @param feature_types Vector specifying feature data types
#' @param method Imputation method: "mean", "median", "phylo"
#' @return Data frame with imputed values
#' @keywords internal
impute_missing_data <- function(tree, features, feature_types, method = "mean") {
  # Create a copy of the data
  imputed <- features
  
  # Process each column
  for(col in colnames(features)) {
    # Skip columns with no missing values
    if(!any(is.na(features[, col]))) {
      next
    }
    
    if(feature_types[col] == "numeric") {
      if(method == "mean") {
        # Simple mean imputation
        col_mean <- mean(features[, col], na.rm = TRUE)
        imputed[is.na(features[, col]), col] <- col_mean
      } else if(method == "median") {
        # Median imputation
        col_median <- median(features[, col], na.rm = TRUE)
        imputed[is.na(features[, col]), col] <- col_median
      } else if(method == "phylo") {
        # Phylogenetic imputation (ancestral state reconstruction)
        if(requireNamespace("phytools", quietly = TRUE)) {
          # Extract non-missing data
          non_missing <- features[!is.na(features[, col]), col, drop = FALSE]
          species_with_data <- rownames(non_missing)
          trait_data <- non_missing[, 1]
          names(trait_data) <- species_with_data
          
          # Perform ancestral state reconstruction
          ace_result <- phytools::fastAnc(tree, trait_data)
          
          # For each missing value, find nearest relative with a value
          for(species in rownames(features)[is.na(features[, col])]) {
            if(!(species %in% tree$tip.label)) next
            
            # Find distance to all other species
            distances <- sapply(species_with_data, function(other_sp) {
              if(other_sp == species) return(Inf)
              mrca_node <- ape::getMRCA(tree, c(species, other_sp))
              dist_to_mrca <- sum(tree$edge.length[get_path_to_root(tree, species, mrca_node)])
              dist_from_mrca <- sum(tree$edge.length[get_path_to_root(tree, other_sp, mrca_node)])
              return(dist_to_mrca + dist_from_mrca)
            })
            
            # Use value from nearest relative
            nearest <- species_with_data[which.min(distances)]
            imputed[species, col] <- features[nearest, col]
          }
        } else {
          # Fallback to mean imputation if phytools not available
          warning("phytools package not available for phylogenetic imputation. Using mean instead.")
          col_mean <- mean(features[, col], na.rm = TRUE)
          imputed[is.na(features[, col]), col] <- col_mean
        }
      }
    } else if(feature_types[col] == "categorical" || feature_types[col] == "binary") {
      if(method == "mean" || method == "median") {
        # Mode imputation for categorical data
        tab <- table(features[, col])
        mode_val <- names(tab)[which.max(tab)]
        imputed[is.na(features[, col]), col] <- mode_val
      } else if(method == "phylo") {
        # For categorical data, use nearest neighbor
        for(species in rownames(features)[is.na(features[, col])]) {
          if(!(species %in% tree$tip.label)) next
          
          # Find non-missing species
          species_with_data <- rownames(features)[!is.na(features[, col])]
          
          # Find distance to all other species
          distances <- sapply(species_with_data, function(other_sp) {
            if(other_sp == species) return(Inf)
            mrca_node <- ape::getMRCA(tree, c(species, other_sp))
            dist_to_mrca <- sum(tree$edge.length[get_path_to_root(tree, species, mrca_node)])
            dist_from_mrca <- sum(tree$edge.length[get_path_to_root(tree, other_sp, mrca_node)])
            return(dist_to_mrca + dist_from_mrca)
          })
          
          # Use value from nearest relative
          nearest <- species_with_data[which.min(distances)]
          imputed[species, col] <- features[nearest, col]
        }
      }
    }
  }
  
  return(imputed)
}

#' Get path from tip to root or specified node
#' 
#' @param tree Phylogenetic tree
#' @param tip_name Name of the tip
#' @param stop_node Node to stop at (default is root)
#' @return Vector of edge indices forming the path
#' @keywords internal
get_path_to_root <- function(tree, tip_name, stop_node = NULL) {
  # Find tip index
  tip_idx <- which(tree$tip.label == tip_name)
  if(length(tip_idx) == 0) {
    stop(paste("Tip not found:", tip_name))
  }
  
  # Set default stop node to root if not specified
  if(is.null(stop_node)) {
    stop_node <- length(tree$tip.label) + 1  # Root node
  }
  
  # Initialize path
  path <- c()
  current <- tip_idx
  
  # Traverse from tip to root
  while(current != stop_node) {
    # Find the edge connecting current node to its parent
    edge_idx <- which(tree$edge[, 2] == current)
    if(length(edge_idx) == 0) {
      break  # Reached root or error
    }
    
    # Add edge to path
    path <- c(path, edge_idx)
    
    # Move to parent
    current <- tree$edge[edge_idx, 1]
    
    # Check if we've reached the stop node
    if(current == stop_node) {
      break
    }
  }
  
  return(path)
}

#===============================================================================
# Correlation Functions
#===============================================================================

#' Analyze correlations between chromosome numbers and genomic features
#' 
#' Tests for phylogenetically-corrected correlations between chromosome 
#' counts and various genomic features
#' 
#' @param integrated_data Integrated data from integrate_genomic_data
#' @param method Correlation method: "PGLS", "pic", "standard"
#' @param features Features to test (NULL for all features)
#' @param model Evolutionary model for PGLS: "BM", "OU", "lambda"
#' @param p_adjustment Method for p-value adjustment: "fdr", "bonferroni", etc.
#' @param ci_level Confidence level for correlation intervals
#' @return List with correlation analysis results
#' @export
analyze_genomic_correlations <- function(integrated_data, 
                                       method = "PGLS",
                                       features = NULL,
                                       model = "BM",
                                       p_adjustment = "fdr",
                                       ci_level = 0.95) {
  # Validate input
  if(!inherits(integrated_data, "acr_integrated")) {
    stop("integrated_data must be created by integrate_genomic_data function")
  }
  
  # Extract data
  tree <- integrated_data$tree
  chr_counts <- integrated_data$chr_counts
  genome_features <- integrated_data$genome_features
  feature_types <- integrated_data$feature_types
  
  # Determine which features to analyze
  if(is.null(features)) {
    features <- colnames(genome_features)
  } else {
    # Check if specified features exist
    missing_features <- setdiff(features, colnames(genome_features))
    if(length(missing_features) > 0) {
      warning(paste("Some features not found:", paste(missing_features, collapse = ", ")))
    }
    features <- intersect(features, colnames(genome_features))
  }
  
  # Filter to numeric features only for correlation
  numeric_features <- features[feature_types[features] == "numeric"]
  categorical_features <- features[feature_types[features] %in% c("categorical", "binary")]
  
  # Initialize results
  results <- list(
    correlations = list(),
    categorical_tests = list(),
    summary = list(
      n_numeric = length(numeric_features),
      n_categorical = length(categorical_features),
      n_significant = 0,
      method = method,
      model = model
    ),
    plots = list()
  )
  
  # Analyze correlations for numeric features
  if(length(numeric_features) > 0) {
    correlation_table <- data.frame(
      Feature = character(0),
      Coefficient = numeric(0),
      p_value = numeric(0),
      CI_lower = numeric(0),
      CI_upper = numeric(0),
      stringsAsFactors = FALSE
    )
    
    for(feature in numeric_features) {
      # Extract feature values
      feature_values <- genome_features[, feature]
      
      # Check for sufficient non-NA data
      valid <- !is.na(feature_values) & !is.na(chr_counts)
      if(sum(valid) < 4) {
        warning(paste("Not enough valid data points for feature:", feature))
        next
      }
      
      # Calculate correlation based on method
      result <- NULL
      
      if(method == "PGLS") {
        # Phylogenetic Generalized Least Squares
        if(requireNamespace("caper", quietly = TRUE)) {
          # Create comparative data
          species <- names(chr_counts)
          data <- data.frame(
            species = species,
            chr = chr_counts,
            feature = feature_values
          )
          
          comp_data <- caper::comparative.data(tree, data, "species")
          
          # Run PGLS
          pgls_model <- tryCatch({
            caper::pgls(chr ~ feature, data = comp_data, lambda = "ML")
          }, error = function(e) {
            warning(paste("PGLS error for feature", feature, ":", e$message))
            return(NULL)
          })
          
          if(!is.null(pgls_model)) {
            coef <- summary(pgls_model)$coefficients["feature", "Estimate"]
            p_value <- summary(pgls_model)$coefficients["feature", "Pr(>|t|)"]
            
            # Calculate confidence intervals
            ci <- confint(pgls_model)["feature", ]
            
            result <- list(
              feature = feature,
              coefficient = coef,
              p_value = p_value,
              ci_lower = ci[1],
              ci_upper = ci[2],
              model = pgls_model
            )
          }
        } else {
          warning("caper package required for PGLS analysis. Using standard correlation instead.")
          method <- "standard"
        }
      }
      
      if(method == "pic") {
        # Phylogenetic Independent Contrasts
        if(requireNamespace("ape", quietly = TRUE)) {
          # Calculate PICs for both variables
          chr_pic <- ape::pic(chr_counts, tree)
          feat_pic <- ape::pic(feature_values, tree)
          
          # Test correlation of contrasts
          cor_test <- cor.test(chr_pic, feat_pic, method = "pearson")
          
          result <- list(
            feature = feature,
            coefficient = cor_test$estimate,
            p_value = cor_test$p.value,
            ci_lower = cor_test$conf.int[1],
            ci_upper = cor_test$conf.int[2],
            test = cor_test
          )
        } else {
          warning("ape package required for PIC analysis. Using standard correlation instead.")
          method <- "standard"
        }
      }
      
      if(method == "standard") {
        # Standard Pearson correlation
        cor_test <- cor.test(chr_counts, feature_values, method = "pearson")
        
        result <- list(
          feature = feature,
          coefficient = cor_test$estimate,
          p_value = cor_test$p.value,
          ci_lower = cor_test$conf.int[1],
          ci_upper = cor_test$conf.int[2],
          test = cor_test
        )
      }
      
      # Store results
      if(!is.null(result)) {
        results$correlations[[feature]] <- result
        
        correlation_table <- rbind(correlation_table, data.frame(
          Feature = feature,
          Coefficient = result$coefficient,
          p_value = result$p_value,
          CI_lower = result$ci_lower,
          CI_upper = result$ci_upper,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Adjust p-values for multiple testing
    if(nrow(correlation_table) > 0) {
      correlation_table$adjusted_p <- p.adjust(correlation_table$p_value, method = p_adjustment)
      correlation_table$significant <- correlation_table$adjusted_p < 0.05
      
      # Sort by absolute correlation strength
      correlation_table <- correlation_table[order(-abs(correlation_table$Coefficient)), ]
      
      # Count significant correlations
      results$summary$n_significant_numeric <- sum(correlation_table$significant)
      
      # Store correlation table
      results$correlation_table <- correlation_table
    }
  }
  
  # Analyze relationships for categorical features
  if(length(categorical_features) > 0) {
    categorical_table <- data.frame(
      Feature = character(0),
      Test = character(0),
      Statistic = numeric(0),
      p_value = numeric(0),
      stringsAsFactors = FALSE
    )
    
    for(feature in categorical_features) {
      # Extract feature values
      feature_values <- genome_features[, feature]
      
      # Check for sufficient non-NA data
      valid <- !is.na(feature_values) & !is.na(chr_counts)
      if(sum(valid) < 4) {
        warning(paste("Not enough valid data points for feature:", feature))
        next
      }
      
      # Use appropriate test based on number of categories
      n_categories <- length(unique(feature_values))
      
      if(n_categories == 2) {
        # Binary feature - t-test or phylogenetic t-test
        if(method == "PGLS" || method == "pic") {
          # Phylogenetic t-test
          if(requireNamespace("phytools", quietly = TRUE)) {
            # Convert chr_counts to named vector
            chr_vector <- chr_counts
            names(chr_vector) <- names(chr_counts)
            
            # Convert feature to factor
            feature_factor <- as.factor(feature_values)
            names(feature_factor) <- names(feature_values)
            
            # Run phylogenetic ANOVA (effectively t-test for binary factor)
            phyanova_result <- tryCatch({
              phytools::phylANOVA(tree, feature_factor, chr_vector)
            }, error = function(e) {
              warning(paste("phylANOVA error for feature", feature, ":", e$message))
              return(NULL)
            })
            
            if(!is.null(phyanova_result)) {
              result <- list(
                feature = feature,
                test = "phyloANOVA",
                statistic = phyanova_result$F,
                p_value = phyanova_result$Pf,
                full_result = phyanova_result
              )
              
              results$categorical_tests[[feature]] <- result
              
              categorical_table <- rbind(categorical_table, data.frame(
                Feature = feature,
                Test = "phyloANOVA",
                Statistic = result$statistic,
                p_value = result$p_value,
                stringsAsFactors = FALSE
              ))
            }
          } else {
            warning("phytools required for phylogenetic t-test. Using standard t-test instead.")
            method = "standard"
          }
        }
        
        if(method == "standard") {
          # Regular t-test
          group1 <- chr_counts[feature_values == unique(feature_values)[1]]
          group2 <- chr_counts[feature_values == unique(feature_values)[2]]
          
          t_test <- t.test(group1, group2)
          
          result <- list(
            feature = feature,
            test = "t-test",
            statistic = t_test$statistic,
            p_value = t_test$p.value,
            full_result = t_test
          )
          
          results$categorical_tests[[feature]] <- result
          
          categorical_table <- rbind(categorical_table, data.frame(
            Feature = feature,
            Test = "t-test",
            Statistic = result$statistic,
            p_value = result$p_value,
            stringsAsFactors = FALSE
          ))
        }
      } else if(n_categories > 2) {
        # Multi-category feature - ANOVA or phylogenetic ANOVA
        if(method == "PGLS" || method == "pic") {
          # Phylogenetic ANOVA
          if(requireNamespace("phytools", quietly = TRUE)) {
            # Convert chr_counts to named vector
            chr_vector <- chr_counts
            names(chr_vector) <- names(chr_counts)
            
            # Convert feature to factor
            feature_factor <- as.factor(feature_values)
            names(feature_factor) <- names(feature_values)
            
            # Run phylogenetic ANOVA
            phyanova_result <- tryCatch({
              phytools::phylANOVA(tree, feature_factor, chr_vector)
            }, error = function(e) {
              warning(paste("phylANOVA error for feature", feature, ":", e$message))
              return(NULL)
            })
            
            if(!is.null(phyanova_result)) {
              result <- list(
                feature = feature,
                test = "phyloANOVA",
                statistic = phyanova_result$F,
                p_value = phyanova_result$Pf,
                full_result = phyanova_result
              )
              
              results$categorical_tests[[feature]] <- result
              
              categorical_table <- rbind(categorical_table, data.frame(
                Feature = feature,
                Test = "phyloANOVA",
                Statistic = result$statistic,
                p_value = result$p_value,
                stringsAsFactors = FALSE
              ))
            }
          } else {
            warning("phytools required for phylogenetic ANOVA. Using standard ANOVA instead.")
            method = "standard"
          }
        }
        
        if(method == "standard") {
          # Regular ANOVA
          anova_model <- aov(chr_counts ~ as.factor(feature_values))
          anova_summary <- summary(anova_model)
          
          result <- list(
            feature = feature,
            test = "ANOVA",
            statistic = anova_summary[[1]]$F[1],
            p_value = anova_summary[[1]]$Pr[1],
            full_result = anova_summary
          )
          
          results$categorical_tests[[feature]] <- result
          
          categorical_table <- rbind(categorical_table, data.frame(
            Feature = feature,
            Test = "ANOVA",
            Statistic = result$statistic,
            p_value = result$p_value,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Adjust p-values for multiple testing
    if(nrow(categorical_table) > 0) {
      categorical_table$adjusted_p <- p.adjust(categorical_table$p_value, method = p_adjustment)
      categorical_table$significant <- categorical_table$adjusted_p < 0.05
      
      # Sort by p-value
      categorical_table <- categorical_table[order(categorical_table$p_value), ]
      
      # Count significant tests
      results$summary$n_significant_categorical <- sum(categorical_table$significant)
      
      # Store categorical test table
      results$categorical_table <- categorical_table
    }
  }
  
  # Update overall count of significant features
  results$summary$n_significant <- 
    ifelse(is.null(results$summary$n_significant_numeric), 0, results$summary$n_significant_numeric) + 
    ifelse(is.null(results$summary$n_significant_categorical), 0, results$summary$n_significant_categorical)
  
  # Create visualization of significant correlations
  if(requireNamespace("ggplot2", quietly = TRUE) && 
     !is.null(results$correlation_table) && 
     sum(results$correlation_table$significant) > 0) {
    
    # Filter to significant correlations
    sig_correlations <- results$correlation_table[results$correlation_table$significant, ]
    
    # Create correlation plot
    corr_plot <- ggplot2::ggplot(sig_correlations, 
                               ggplot2::aes(y = reorder(Feature, Coefficient), 
                                          x = Coefficient, 
                                          xmin = CI_lower, 
                                          xmax = CI_upper)) +
      ggplot2::geom_point(size = 3, color = "steelblue") +
      ggplot2::geom_errorbarh(height = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = "Significant Correlations with Chromosome Number",
                  subtitle = paste("Method:", method),
                  x = "Correlation Coefficient",
                  y = "Genomic Feature") +
      ggplot2::theme_minimal()
    
    results$plots$correlation_plot <- corr_plot
    
    # Create scatterplots for top significant correlations
    top_n <- min(3, nrow(sig_correlations))
    top_features <- sig_correlations$Feature[1:top_n]
    
    scatter_plots <- list()
    
    for(feature in top_features) {
      feature_data <- genome_features[, feature]
      
      scatter_plot <- ggplot2::ggplot(data = data.frame(
        chr = chr_counts,
        feature = feature_data
      ), ggplot2::aes(x = feature, y = chr)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "lm", color = "blue") +
        ggplot2::labs(
          title = paste("Chromosome Number vs.", feature),
          x = feature,
          y = "Chromosome Number"
        ) +
        ggplot2::theme_minimal()
      
      scatter_plots[[feature]] <- scatter_plot
    }
    
    results$plots$scatter_plots <- scatter_plots
  }
  
  # Create visualization for categorical features
  if(requireNamespace("ggplot2", quietly = TRUE) && 
     !is.null(results$categorical_table) && 
     sum(results$categorical_table$significant) > 0) {
    
    # Filter to significant tests
    sig_categorical <- results$categorical_table[results$categorical_table$significant, ]
    
    # Create boxplots for significant categorical features
    boxplot_list <- list()
    
    for(i in 1:nrow(sig_categorical)) {
      feature <- sig_categorical$Feature[i]
      feature_values <- genome_features[, feature]
      
      box_data <- data.frame(
        chr = chr_counts,
        factor = as.factor(feature_values)
      )
      
      boxplot <- ggplot2::ggplot(box_data, ggplot2::aes(x = factor, y = chr, fill = factor)) +
        ggplot2::geom_boxplot() +
        ggplot2::labs(
          title = paste("Chromosome Number by", feature),
          subtitle = paste("p =", format(sig_categorical$adjusted_p[i], digits = 3)),
          x = feature,
          y = "Chromosome Number"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
      
      boxplot_list[[feature]] <- boxplot
    }
    
    results$plots$categorical_plots <- boxplot_list
  }
  
  return(results)
}

#===============================================================================
# Synteny Analysis Functions
#===============================================================================

#' Analyze synteny conservation and chromsome evolution
#' 
#' Correlates synteny conservation between species with chromosome 
#' rearrangements and number changes
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param synteny_data Data frame with pairwise synteny data
#' @param synteny_metric Column in synteny_data to use as conservation metric
#' @param metric_type Type of synteny metric: "similarity" or "distance"
#' @param ancestral_recon Optional ancestral reconstruction results
#' @param rate_results Optional rate analysis results
#' @param permutations Number of permutations for significance testing
#' @return List with synteny analysis results
#' @export
analyze_synteny_conservation <- function(tree,
                                      chr_counts,
                                      synteny_data,
                                      synteny_metric = "conserved_blocks",
                                      metric_type = "similarity",
                                      ancestral_recon = NULL,
                                      rate_results = NULL,
                                      permutations = 1000) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  if(!is.data.frame(synteny_data)) {
    stop("synteny_data must be a data frame")
  }
  
  # Check required columns in synteny data
  required_cols <- c("species1", "species2", synteny_metric)
  if(!all(required_cols %in% colnames(synteny_data))) {
    stop(paste("synteny_data must contain columns:", paste(required_cols, collapse = ", ")))
  }
  
  # Initialize results
  results <- list(
    synteny_vs_distance = list(),
    synteny_vs_changes = list(),
    node_results = list(),
    summary = list(
      n_pairs = nrow(synteny_data),
      synteny_metric = synteny_metric,
      metric_type = metric_type
    ),
    plots = list()
  )
  
  # Calculate phylogenetic distances between species pairs
  if(requireNamespace("ape", quietly = TRUE)) {
    # Get species in tree
    tree_species <- tree$tip.label
    
    # Filter synteny data to species in tree
    valid_pairs <- synteny_data$species1 %in% tree_species & synteny_data$species2 %in% tree_species
    if(sum(valid_pairs) == 0) {
      stop("No synteny pairs found with both species in the tree")
    }
    
    filtered_synteny <- synteny_data[valid_pairs, ]
    
    # Calculate phylogenetic distance for each pair
    distances <- numeric(nrow(filtered_synteny))
    
    for(i in 1:nrow(filtered_synteny)) {
      sp1 <- filtered_synteny$species1[i]
      sp2 <- filtered_synteny$species2[i]
      
      # Calculate patristic distance (sum of branch lengths)
      distances[i] <- ape::cophenetic.phylo(tree)[sp1, sp2]
    }
    
    # Add distances to data
    filtered_synteny$phylo_distance <- distances
    
    # Calculate absolute chromosome count differences
    chr_diffs <- numeric(nrow(filtered_synteny))
    
    for(i in 1:nrow(filtered_synteny)) {
      sp1 <- filtered_synteny$species1[i]
      sp2 <- filtered_synteny$species2[i]
      
      # Get chromosome counts (if available)
      chr1 <- ifelse(sp1 %in% names(chr_counts), chr_counts[sp1], NA)
      chr2 <- ifelse(sp2 %in% names(chr_counts), chr_counts[sp2], NA)
      
      # Calculate difference
      chr_diffs[i] <- ifelse(!is.na(chr1) && !is.na(chr2), abs(chr1 - chr2), NA)
    }
    
    # Add chromosome differences to data
    filtered_synteny$chr_diff <- chr_diffs
    
    # Analyze relationship between synteny and phylogenetic distance
    valid_dist <- !is.na(filtered_synteny$phylo_distance) & 
                !is.na(filtered_synteny[[synteny_metric]])
    
    if(sum(valid_dist) >= 4) {
      # For distance metrics, we expect positive correlation with phylo distance
      # For similarity metrics, we expect negative correlation with phylo distance
      cor_method <- "spearman"  # Non-parametric correlation
      
      dist_cor <- cor.test(
        filtered_synteny$phylo_distance[valid_dist],
        filtered_synteny[[synteny_metric]][valid_dist],
        method = cor_method
      )
      
      # Store results
      results$synteny_vs_distance <- list(
        correlation = dist_cor$estimate,
        p_value = dist_cor$p.value,
        expected_sign = ifelse(metric_type == "similarity", "negative", "positive"),
        observed_sign = ifelse(dist_cor$estimate < 0, "negative", "positive"),
        matches_expectation = (metric_type == "similarity" && dist_cor$estimate < 0) ||
                              (metric_type == "distance" && dist_cor$estimate > 0),
        test = dist_cor
      )
      
      # Create scatter plot
      if(requireNamespace("ggplot2", quietly = TRUE)) {
        dist_plot <- ggplot2::ggplot(filtered_synteny[valid_dist, ], 
                                  ggplot2::aes(x = phylo_distance, y = .data[[synteny_metric]])) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "lm", formula = 'y ~ x') +
          ggplot2::labs(
            title = "Synteny Conservation vs. Phylogenetic Distance",
            subtitle = paste("Correlation =", round(dist_cor$estimate, 3), 
                           "(p =", format(dist_cor$p.value, digits = 3), ")"),
            x = "Phylogenetic Distance",
            y = synteny_metric
          ) +
          ggplot2::theme_minimal()
        
        results$plots$distance_plot <- dist_plot
      }
    }
    
    # Analyze relationship between synteny and chromosome number differences
    valid_chr <- !is.na(filtered_synteny$chr_diff) & 
               !is.na(filtered_synteny[[synteny_metric]])
    
    if(sum(valid_chr) >= 4) {
      # For similarity metrics, we expect negative correlation with chr differences
      # For distance metrics, we expect positive correlation with chr differences
      chr_cor <- cor.test(
        filtered_synteny$chr_diff[valid_chr],
        filtered_synteny[[synteny_metric]][valid_chr],
        method = cor_method
      )
      
      # Run permutation test
      observed_cor <- chr_cor$estimate
      permutation_results <- numeric(permutations)
      
      for(i in 1:permutations) {
        # Shuffle chromosome differences
        shuffled_diffs <- sample(filtered_synteny$chr_diff[valid_chr])
        
        # Calculate correlation
        permutation_results[i] <- cor(
          shuffled_diffs,
          filtered_synteny[[synteny_metric]][valid_chr],
          method = cor_method
        )
      }
      
      # Calculate empirical p-value
      if(observed_cor > 0) {
        perm_p_value <- sum(permutation_results >= observed_cor) / permutations
      } else {
        perm_p_value <- sum(permutation_results <= observed_cor) / permutations
      }
      
      # Store results
      results$synteny_vs_changes <- list(
        correlation = chr_cor$estimate,
        p_value = chr_cor$p.value,
        perm_p_value = perm_p_value,
        expected_sign = ifelse(metric_type == "similarity", "negative", "positive"),
        observed_sign = ifelse(chr_cor$estimate < 0, "negative", "positive"),
        matches_expectation = (metric_type == "similarity" && chr_cor$estimate < 0) ||
                              (metric_type == "distance" && chr_cor$estimate > 0),
        test = chr_cor,
        permutation_dist = permutation_results
      )
      
      # Create scatter plot
      if(requireNamespace("ggplot2", quietly = TRUE)) {
        chr_plot <- ggplot2::ggplot(filtered_synteny[valid_chr, ], 
                                 ggplot2::aes(x = chr_diff, y = .data[[synteny_metric]])) +
          ggplot2::geom_point() +
          ggplot2::geom_smooth(method = "lm", formula = 'y ~ x') +
          ggplot2::labs(
            title = "Synteny Conservation vs. Chromosome Number Differences",
            subtitle = paste("Correlation =", round(chr_cor$estimate, 3), 
                           "(p =", format(perm_p_value, digits = 3), ", permutation)"),
            x = "Chromosome Number Difference",
            y = synteny_metric
          ) +
          ggplot2::theme_minimal()
        
        results$plots$chr_diff_plot <- chr_plot
      }
    }
    
    # Store processed synteny data
    results$synteny_data <- filtered_synteny
  } else {
    warning("ape package required for phylogenetic distance calculation")
  }
  
  # If ancestral reconstruction and rate results are provided, perform additional analyses
  if(!is.null(ancestral_recon) && !is.null(rate_results)) {
    # Add analyses of synteny conservation vs. evolutionary rates
    # This is a placeholder for a more complex analysis that could be implemented
    results$summary$additional_analyses <- "Additional analyses with ancestral reconstruction not yet implemented"
  }
  
  # Create summary
  results$summary$synteny_distance_correlation <- results$synteny_vs_distance$correlation
  results$summary$synteny_distance_p <- results$synteny_vs_distance$p_value
  results$summary$synteny_chr_correlation <- results$synteny_vs_changes$correlation
  results$summary$synteny_chr_p <- results$synteny_vs_changes$perm_p_value
  
  return(results)
}

#===============================================================================
# Genome Size Analysis Functions
#===============================================================================

#' Analyze relationship between genome size and chromosome evolution
#' 
#' Tests for correlations between genome size, chromosome numbers,
#' and chromosome evolution rates
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param genome_sizes Named vector of genome sizes (in Mb)
#' @param ancestral_recon Optional ancestral reconstruction results
#' @param rate_results Optional rate analysis results
#' @param method Analysis method: "PGLS", "pic", "standard"
#' @param transform Whether to log-transform genome sizes
#' @param normalize Whether to normalize data before analysis
#' @return List with genome size analysis results
#' @export
analyze_genome_size_relationship <- function(tree,
                                          chr_counts,
                                          genome_sizes,
                                          ancestral_recon = NULL,
                                          rate_results = NULL,
                                          method = "PGLS",
                                          transform = TRUE,
                                          normalize = FALSE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  if(!is.vector(genome_sizes) || is.null(names(genome_sizes))) {
    stop("genome_sizes must be a named vector")
  }
  
  # Find species with both chromosome counts and genome sizes
  common_species <- intersect(names(chr_counts), names(genome_sizes))
  common_species <- intersect(common_species, tree$tip.label)
  
  if(length(common_species) < 4) {
    stop("At least 4 species must have both chromosome counts and genome sizes")
  }
  
  # Extract matched data
  matched_chr <- chr_counts[common_species]
  matched_sizes <- genome_sizes[common_species]
  
  # Apply transformations if requested
  if(transform) {
    matched_sizes <- log10(matched_sizes)
  }
  
  if(normalize) {
    matched_chr <- scale(matched_chr)
    matched_sizes <- scale(matched_sizes)
  }
  
  # Prune tree to common species
  pruned_tree <- ape::keep.tip(tree, common_species)
  
  # Initialize results
  results <- list(
    direct_correlation = NULL,
    chromosome_evolution = NULL,
    summary = list(
      n_species = length(common_species),
      transform = transform,
      normalize = normalize,
      method = method
    ),
    plots = list()
  )
  
  # Analyze direct correlation between genome size and chromosome number
  correlation_result <- NULL
  
  if(method == "PGLS") {
    # Phylogenetic Generalized Least Squares
    if(requireNamespace("caper", quietly = TRUE)) {
      # Create comparative data
      data <- data.frame(
        species = common_species,
        chr = matched_chr,
        size = matched_sizes
      )
      
      comp_data <- caper::comparative.data(pruned_tree, data, "species")
      
      # Run PGLS
      pgls_model <- tryCatch({
        caper::pgls(chr ~ size, data = comp_data, lambda = "ML")
      }, error = function(e) {
        warning(paste("PGLS error:", e$message))
        return(NULL)
      })
      
      if(!is.null(pgls_model)) {
        coef <- summary(pgls_model)$coefficients["size", "Estimate"]
        p_value <- summary(pgls_model)$coefficients["size", "Pr(>|t|)"]
        
        # Calculate confidence intervals
        ci <- confint(pgls_model)["size", ]
        
        correlation_result <- list(
          coefficient = coef,
          p_value = p_value,
          ci_lower = ci[1],
          ci_upper = ci[2],
          lambda = pgls_model$param$lambda,
          model = pgls_model
        )
      }
    } else {
      warning("caper package required for PGLS analysis. Using standard correlation instead.")
      method <- "standard"
    }
  } else if(method == "pic") {
    # Phylogenetic Independent Contrasts
    if(requireNamespace("ape", quietly = TRUE)) {
      # Calculate PICs for both variables
      chr_pic <- ape::pic(matched_chr, pruned_tree)
      size_pic <- ape::pic(matched_sizes, pruned_tree)
      
      # Test correlation of contrasts
      cor_test <- cor.test(chr_pic, size_pic, method = "pearson")
      
      correlation_result <- list(
        coefficient = cor_test$estimate,
        p_value = cor_test$p.value,
        ci_lower = cor_test$conf.int[1],
        ci_upper = cor_test$conf.int[2],
        test = cor_test
      )
    } else {
      warning("ape package required for PIC analysis. Using standard correlation instead.")
      method <- "standard"
    }
  }
  
  if(method == "standard" || is.null(correlation_result)) {
    # Standard Pearson correlation
    cor_test <- cor.test(matched_chr, matched_sizes, method = "pearson")
    
    correlation_result <- list(
      coefficient = cor_test$estimate,
      p_value = cor_test$p.value,
      ci_lower = cor_test$conf.int[1],
      ci_upper = cor_test$conf.int[2],
      test = cor_test
    )
  }
  
  # Store correlation results
  results$direct_correlation <- correlation_result
  
  # Create scatter plot
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Basic scatter plot
    scatter_plot <- ggplot2::ggplot(data.frame(
      Chr = matched_chr,
      Size = matched_sizes,
      Species = common_species
    ), ggplot2::aes(x = Size, y = Chr)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", formula = 'y ~ x') +
      ggplot2::labs(
        title = "Chromosome Number vs. Genome Size",
        subtitle = paste(
          "Correlation =", round(correlation_result$coefficient, 3),
          "(p =", format(correlation_result$p_value, digits = 3), ")"
        ),
        x = ifelse(transform, "log10(Genome Size) (Mb)", "Genome Size (Mb)"),
        y = "Chromosome Number"
      ) +
      ggplot2::theme_minimal()
    
    # Add species labels if not too many
    if(length(common_species) <= 20) {
      scatter_plot <- scatter_plot +
        ggplot2::geom_text(ggplot2::aes(label = Species), 
                         nudge_x = 0.1, nudge_y = 0.1, size = 3)
    }
    
    results$plots$scatter_plot <- scatter_plot
  }
  
  # Additional analyses if rate results provided
  if(!is.null(rate_results) && !is.null(ancestral_recon)) {
    # Placeholder for additional analyses
    results$summary$additional_analyses <- "Additional analyses with rates not yet implemented"
  }
  
  # Create summary
  results$summary$correlation <- correlation_result$coefficient
  results$summary$p_value <- correlation_result$p_value
  results$summary$significant <- correlation_result$p_value < 0.05
  
  # Interpret results
  if(results$summary$significant) {
    if(correlation_result$coefficient > 0) {
      results$summary$interpretation <- "Significant positive correlation between genome size and chromosome number"
    } else {
      results$summary$interpretation <- "Significant negative correlation between genome size and chromosome number"
    }
  } else {
    results$summary$interpretation <- "No significant correlation between genome size and chromosome number"
  }
  
  return(results)
}

#===============================================================================
# Visualization Functions
#===============================================================================

#' Create integrated visualization of chromosome evolution and genomic features
#' 
#' Generates a multi-panel visualization showing chromosome number evolution
#' alongside changes in genomic features
#' 
#' @param integrated_data Integrated data from integrate_genomic_data
#' @param ancestral_states Ancestral chromosome reconstruction results
#' @param features Features to visualize (NULL for automatic selection of top correlated features)
#' @param highlight_nodes Optional vector of nodes to highlight
#' @param node_labels Whether to show node labels
#' @param color_palette Color palette for feature heatmap
#' @param use_patchwork Whether to use patchwork package for combining plots
#' @param output_file Optional path to save the visualization
#' @return ggplot or patchwork object with visualization
#' @export
visualize_integrated_evolution <- function(integrated_data,
                                        ancestral_states,
                                        features = NULL,
                                        highlight_nodes = NULL,
                                        node_labels = FALSE,
                                        color_palette = viridis::viridis(100),
                                        use_patchwork = TRUE,
                                        output_file = NULL) {
  # Validate input
  if(!inherits(integrated_data, "acr_integrated")) {
    stop("integrated_data must be created by integrate_genomic_data function")
  }
  
  # Extract data
  tree <- integrated_data$tree
  chr_counts <- integrated_data$chr_counts
  genome_features <- integrated_data$genome_features
  feature_types <- integrated_data$feature_types
  
  # Select features to visualize
  if(is.null(features)) {
    # Auto-select features most correlated with chromosome number
    feature_cors <- sapply(colnames(genome_features), function(feat) {
      if(feature_types[feat] == "numeric") {
        valid <- !is.na(genome_features[, feat]) & !is.na(chr_counts)
        if(sum(valid) >= 3) {
          return(abs(cor(chr_counts[valid], genome_features[valid, feat], 
                        method = "spearman", use = "pairwise.complete.obs")))
        }
      }
      return(0)
    })
    
    # Select top 3 numeric features
    features <- names(sort(feature_cors, decreasing = TRUE))[1:min(3, sum(feature_cors > 0))]
  }
  
  # Check if required packages are available
  if(!all(c("ggplot2", "ggtree") %in% installed.packages()[,"Package"])) {
    stop("ggplot2 and ggtree packages are required for visualization")
  }
  
  # Create chromosome mapping plot
  chr_mapping <- NULL
  if(!is.null(ancestral_states)) {
    # Combine tip and ancestral states
    all_states <- ancestral_states
    
    # Create data for tree nodes
    node_data <- data.frame(
      node = 1:(length(tree$tip.label) + tree$Nnode),
      chr_num = as.numeric(all_states[as.character(1:(length(tree$tip.label) + tree$Nnode))]),
      is_tip = 1:length(all_states) <= length(tree$tip.label),
      is_highlighted = 1:length(all_states) %in% highlight_nodes,
      stringsAsFactors = FALSE
    )
    
    # Create tree with chromosome mapping
    chr_mapping <- ggtree::ggtree(tree, aes(color = chr_num), size = 1.2) %<+% node_data +
      scale_color_viridis_c(name = "Chromosome\nNumber", option = "magma") +
      ggtree::geom_tiplab(aes(label = label), offset = 0.1, size = 3) +
      theme_tree2() +
      labs(title = "Chromosome Number Evolution")
    
    # Add node labels if requested
    if(node_labels) {
      chr_mapping <- chr_mapping + 
        ggtree::geom_nodelab(aes(label = round(chr_num, 1)), size = 2.5, hjust = -0.1)
    }
    
    # Highlight specified nodes if any
    if(!is.null(highlight_nodes) && any(node_data$is_highlighted)) {
      chr_mapping <- chr_mapping +
        ggtree::geom_point2(aes(subset = is_highlighted), size = 4, shape = 21, fill = "yellow", alpha = 0.8)
    }
  }
  
  # Create feature heatmap plots
  feature_plots <- list()
  for(feature in features) {
    if(feature_types[feature] == "numeric") {
      # Create data frame for gheatmap
      feature_df <- data.frame(feature_val = genome_features[, feature])
      rownames(feature_df) <- rownames(genome_features)
      
      # Create feature plot
      feat_plot <- ggtree::gheatmap(ggtree::ggtree(tree), 
                                  feature_df, 
                                  offset = 0.2, 
                                  width = 0.2, 
                                  colnames = FALSE) +
        scale_fill_viridis_c(name = feature, na.value = "grey80", option = "inferno") +
        labs(title = feature)
      
      feature_plots[[feature]] <- feat_plot
    } else {
      # Handle categorical features differently
      feature_df <- data.frame(feature_val = as.factor(genome_features[, feature]))
      rownames(feature_df) <- rownames(genome_features)
      
      # Create feature plot for categorical data
      n_levels <- length(unique(genome_features[, feature]))
      color_set <- viridis::viridis(max(3, n_levels))[1:n_levels]
      
      feat_plot <- ggtree::gheatmap(ggtree::ggtree(tree), 
                                  feature_df, 
                                  offset = 0.2, 
                                  width = 0.2, 
                                  colnames = FALSE) +
        scale_fill_manual(name = feature, values = color_set, na.value = "grey80") +
        labs(title = feature)
      
      feature_plots[[feature]] <- feat_plot
    }
  }
  
  # Combine plots
  if(use_patchwork && requireNamespace("patchwork", quietly = TRUE)) {
    if(!is.null(chr_mapping)) {
      # Combine chromosome mapping with feature plots
      combined_plot <- chr_mapping
      
      for(p in feature_plots) {
        combined_plot <- combined_plot + p
      }
      
      # Arrange in a grid with chromosome mapping larger
      final_viz <- combined_plot + 
        patchwork::plot_layout(ncol = 1, heights = c(3, rep(1, length(feature_plots))))
    } else {
      # Only feature plots
      combined_plot <- feature_plots[[1]]
      
      for(i in 2:length(feature_plots)) {
        combined_plot <- combined_plot + feature_plots[[i]]
      }
      
      final_viz <- combined_plot + 
        patchwork::plot_layout(ncol = 1)
    }
  } else {
    # Return list of plots if patchwork not available
    final_viz <- list(
      chr_mapping = chr_mapping,
      feature_plots = feature_plots
    )
  }
  
  # Save plot if requested
  if(!is.null(output_file) && requireNamespace("ggplot2", quietly = TRUE)) {
    ggsave(output_file, plot = final_viz, width = 12, height = 8 + 2*length(features), dpi = 300)
  }
  
  return(final_viz)
}

#' Create faceted scatterplots of chromosome counts vs genomic features
#' 
#' @param integrated_data Integrated data from integrate_genomic_data
#' @param features Features to plot (NULL for all numeric features)
#' @param color_by Optional factor to use for coloring points
#' @param add_trendlines Whether to add linear regression lines
#' @param facet_ncol Number of columns for the faceted plot
#' @param output_file Optional path to save the plot
#' @return ggplot object with faceted scatterplots
#' @export
plot_chr_feature_relationships <- function(integrated_data,
                                         features = NULL,
                                         color_by = NULL,
                                         add_trendlines = TRUE,
                                         facet_ncol = 2,
                                         output_file = NULL) {
  # Validate input
  if(!inherits(integrated_data, "acr_integrated")) {
    stop("integrated_data must be created by integrate_genomic_data function")
  }
  
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for this function")
  }
  
  # Extract data
  chr_counts <- integrated_data$chr_counts
  genome_features <- integrated_data$genome_features
  feature_types <- integrated_data$feature_types
  
  # Determine features to plot
  if(is.null(features)) {
    # Use all numeric features
    features <- names(feature_types)[feature_types == "numeric"]
    
    if(length(features) == 0) {
      stop("No numeric features found in the integrated data")
    }
  } else {
    # Check that specified features exist and are numeric
    non_existing <- setdiff(features, colnames(genome_features))
    if(length(non_existing) > 0) {
      warning(paste("Features not found:", paste(non_existing, collapse = ", ")))
    }
    
    non_numeric <- features[feature_types[features] != "numeric"]
    if(length(non_numeric) > 0) {
      warning(paste("Non-numeric features will be excluded:", paste(non_numeric, collapse = ", ")))
    }
    
    # Filter to existing numeric features
    features <- intersect(features, colnames(genome_features)[feature_types == "numeric"])
  }
  
  # Create long-format data frame for plotting
  plot_data <- data.frame(
    Species = rep(names(chr_counts), length(features)),
    Chromosome_Number = rep(chr_counts, length(features)),
    Feature = rep(features, each = length(chr_counts)),
    Value = unlist(lapply(features, function(f) genome_features[names(chr_counts), f])),
    stringsAsFactors = FALSE
  )
  
  # Add color variable if provided
  if(!is.null(color_by)) {
    if(color_by %in% colnames(genome_features)) {
      plot_data$Color <- rep(genome_features[names(chr_counts), color_by], length(features))
      color_aes <- ggplot2::aes(color = Color)
    } else {
      warning(paste("Color variable", color_by, "not found in genome_features"))
      color_aes <- NULL
    }
  } else {
    color_aes <- NULL
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Value, y = Chromosome_Number)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::facet_wrap(~ Feature, scales = "free_x", ncol = facet_ncol) +
    ggplot2::labs(
      title = "Chromosome Number vs. Genomic Features",
      x = "Feature Value",
      y = "Chromosome Number"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "lightblue", colour = "darkblue"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, color = "grey80")
    )
  
  # Add color if provided
  if(!is.null(color_aes)) {
    p <- p + color_aes + 
      ggplot2::guides(color = ggplot2::guide_legend(title = color_by))
  }
  
  # Add trend lines if requested
  if(add_trendlines) {
    p <- p + ggplot2::geom_smooth(method = "lm", formula = 'y ~ x', se = TRUE,
                               color = "blue", linetype = "dashed", size = 1)
  }
  
  # Add species labels if not too many
  if(length(chr_counts) <= 20) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = Species), 
                            hjust = -0.1, vjust = 0.5, size = 3)
  }
  
  # Save plot if requested
  if(!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, 
                  width = min(4 * facet_ncol, 12), 
                  height = 3 * ceiling(length(features) / facet_ncol),
                  dpi = 300)
  }
  
  return(p)
}

#===============================================================================
# Repeat Content Analysis Functions
#===============================================================================

#' Analyze relationship between repeat content and chromosome evolution
#' 
#' Tests for correlations between various classes of repetitive elements
#' and chromosome numbers, accounting for phylogenetic relationships
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param repeat_data Data frame with repeat content data (rows = species, columns = repeat classes)
#' @param genome_sizes Optional named vector of genome sizes for normalization
#' @param method Analysis method: "PGLS", "pic", "standard"
#' @param normalize Whether to normalize repeat content by genome size
#' @param percent Whether repeat data is already in percentages
#' @param p_adjustment Method for p-value adjustment: "fdr", "bonferroni", etc.
#' @param permutations Number of permutations for significance testing
#' @return List with repeat content analysis results
#' @export
analyze_repeat_content <- function(tree,
                                 chr_counts,
                                 repeat_data,
                                 genome_sizes = NULL,
                                 method = "PGLS",
                                 normalize = TRUE,
                                 percent = FALSE,
                                 p_adjustment = "fdr",
                                 permutations = 1000) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  if(!is.data.frame(repeat_data)) {
    stop("repeat_data must be a data frame")
  }
  
  # Check if repeat_data has row names
  if(is.null(rownames(repeat_data))) {
    stop("repeat_data must have row names corresponding to species names")
  }
  
  # Get common taxa across datasets
  tree_tips <- tree$tip.label
  chr_taxa <- names(chr_counts)
  repeat_taxa <- rownames(repeat_data)
  
  # Find species present in all datasets
  common_taxa <- Reduce(intersect, list(tree_tips, chr_taxa, repeat_taxa))
  
  if(length(common_taxa) < 4) {
    stop("At least 4 taxa must be present in the tree, chromosome data, and repeat data")
  }
  
  # Prune tree and data to common taxa
  pruned_tree <- ape::keep.tip(tree, common_taxa)
  pruned_chr <- chr_counts[common_taxa]
  pruned_repeat <- repeat_data[common_taxa, , drop = FALSE]
  
  # Normalize repeat data by genome size if requested
  if(normalize && !percent && !is.null(genome_sizes)) {
    # Check if genome sizes are available for all species
    common_with_size <- intersect(common_taxa, names(genome_sizes))
    
    if(length(common_with_size) < length(common_taxa)) {
      warning(paste(length(common_taxa) - length(common_with_size), 
                  "species lack genome size data; normalization may be incomplete"))
    }
    
    # Normalize by genome size (convert to percentage)
    for(i in 1:nrow(pruned_repeat)) {
      species <- rownames(pruned_repeat)[i]
      if(species %in% names(genome_sizes)) {
        pruned_repeat[i, ] <- 100 * pruned_repeat[i, ] / genome_sizes[species]
      }
    }
  }
  
  # Initialize results
  results <- list(
    correlations = list(),
    summary = list(
      n_taxa = length(common_taxa),
      n_repeat_classes = ncol(pruned_repeat),
      method = method,
      normalized = normalize && !percent && !is.null(genome_sizes),
      percent = percent || (normalize && !is.null(genome_sizes))
    ),
    plots = list()
  )
  
  # Test correlation for each repeat class
  correlation_table <- data.frame(
    Repeat_Class = character(0),
    Coefficient = numeric(0),
    p_value = numeric(0),
    CI_lower = numeric(0),
    CI_upper = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for(repeat_class in colnames(pruned_repeat)) {
    # Extract repeat data
    repeat_values <- pruned_repeat[, repeat_class]
    
    # Check for sufficient non-NA data
    valid <- !is.na(repeat_values) & !is.na(pruned_chr)
    if(sum(valid) < 4) {
      warning(paste("Not enough valid data points for repeat class:", repeat_class))
      next
    }
    
    # Calculate correlation based on method
    result <- NULL
    
    if(method == "PGLS") {
      # Phylogenetic Generalized Least Squares
      if(requireNamespace("caper", quietly = TRUE)) {
        # Create comparative data
        data <- data.frame(
          species = common_taxa,
          chr = pruned_chr,
          repeat_val = repeat_values
        )
        
        comp_data <- caper::comparative.data(pruned_tree, data, "species")
        
        # Run PGLS
        pgls_model <- tryCatch({
          caper::pgls(chr ~ repeat_val, data = comp_data, lambda = "ML")
        }, error = function(e) {
          warning(paste("PGLS error for repeat class", repeat_class, ":", e$message))
          return(NULL)
        })
        
        if(!is.null(pgls_model)) {
          coef <- summary(pgls_model)$coefficients["repeat_val", "Estimate"]
          p_value <- summary(pgls_model)$coefficients["repeat_val", "Pr(>|t|)"]
          
          # Calculate confidence intervals
          ci <- confint(pgls_model)["repeat_val", ]
          
          result <- list(
            repeat_class = repeat_class,
            coefficient = coef,
            p_value = p_value,
            ci_lower = ci[1],
            ci_upper = ci[2],
            lambda = pgls_model$param$lambda,
            model = pgls_model
          )
        }
      } else {
        warning("caper package required for PGLS analysis. Using standard correlation instead.")
        method <- "standard"
      }
    }
    
    if(method == "pic") {
      # Phylogenetic Independent Contrasts
      if(requireNamespace("ape", quietly = TRUE)) {
        # Calculate PICs for both variables
        chr_pic <- ape::pic(pruned_chr, pruned_tree)
        repeat_pic <- ape::pic(repeat_values, pruned_tree)
        
        # Test correlation of contrasts
        cor_test <- cor.test(chr_pic, repeat_pic, method = "pearson")
        
        result <- list(
          repeat_class = repeat_class,
          coefficient = cor_test$estimate,
          p_value = cor_test$p.value,
          ci_lower = cor_test$conf.int[1],
          ci_upper = cor_test$conf.int[2],
          test = cor_test
        )
      } else {
        warning("ape package required for PIC analysis. Using standard correlation instead.")
        method <- "standard"
      }
    }
    
    if(method == "standard" || is.null(result)) {
      # Standard Pearson correlation
      cor_test <- cor.test(pruned_chr, repeat_values, method = "pearson")
      
      result <- list(
        repeat_class = repeat_class,
        coefficient = cor_test$estimate,
        p_value = cor_test$p.value,
        ci_lower = cor_test$conf.int[1],
        ci_upper = cor_test$conf.int[2],
        test = cor_test
      )
    }
    
    # Store results
    results$correlations[[repeat_class]] <- result
    
    correlation_table <- rbind(correlation_table, data.frame(
      Repeat_Class = repeat_class,
      Coefficient = result$coefficient,
      p_value = result$p_value,
      CI_lower = result$ci_lower,
      CI_upper = result$ci_upper,
      stringsAsFactors = FALSE
    ))
  }
  
  # Adjust p-values for multiple testing
  if(nrow(correlation_table) > 0) {
    correlation_table$adjusted_p <- p.adjust(correlation_table$p_value, method = p_adjustment)
    correlation_table$significant <- correlation_table$adjusted_p < 0.05
    
    # Sort by absolute correlation strength
    correlation_table <- correlation_table[order(-abs(correlation_table$Coefficient)), ]
    
    # Count significant correlations
    results$summary$n_significant <- sum(correlation_table$significant)
    
    # Store correlation table
    results$correlation_table <- correlation_table
  }
  
  # Create visualization of significant correlations
  if(requireNamespace("ggplot2", quietly = TRUE) && 
     !is.null(results$correlation_table) && 
     sum(results$correlation_table$significant) > 0) {
    
    # Filter to significant correlations
    sig_correlations <- results$correlation_table[results$correlation_table$significant, ]
    
    # Create correlation plot
    corr_plot <- ggplot2::ggplot(sig_correlations, 
                               ggplot2::aes(y = reorder(Repeat_Class, Coefficient), 
                                          x = Coefficient, 
                                          xmin = CI_lower, 
                                          xmax = CI_upper)) +
      ggplot2::geom_point(size = 3, color = "steelblue") +
      ggplot2::geom_errorbarh(height = 0.2) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(title = "Significant Correlations Between Repeat Content and Chromosome Number",
                  subtitle = paste("Method:", method),
                  x = "Correlation Coefficient",
                  y = "Repeat Class") +
      ggplot2::theme_minimal()
    
    results$plots$correlation_plot <- corr_plot
    
    # Create scatter plots for top significant correlations
    top_n <- min(5, nrow(sig_correlations))
    top_repeats <- sig_correlations$Repeat_Class[1:top_n]
    
    scatter_data <- data.frame(
      Species = rep(common_taxa, top_n),
      Chromosome_Number = rep(pruned_chr, top_n),
      Repeat_Class = rep(top_repeats, each = length(common_taxa)),
      Repeat_Value = unlist(lapply(top_repeats, function(r) pruned_repeat[, r])),
      stringsAsFactors = FALSE
    )
    
    scatter_plot <- ggplot2::ggplot(scatter_data, 
                                  ggplot2::aes(x = Repeat_Value, y = Chromosome_Number)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::geom_smooth(method = "lm", formula = 'y ~ x', se = TRUE, color = "blue") +
      ggplot2::facet_wrap(~ Repeat_Class, scales = "free_x") +
      ggplot2::labs(
        title = "Chromosome Number vs. Top Correlated Repeat Classes",
        x = ifelse(results$summary$percent, "Repeat Content (%)", "Repeat Content"),
        y = "Chromosome Number"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "lightblue", colour = "darkblue"),
        strip.text = ggplot2::element_text(face = "bold")
      )
    
    results$plots$scatter_plot <- scatter_plot
  }
  
  return(results)
}

#===============================================================================
# Gene Density Analysis Functions
#===============================================================================

#' Analyze relationship between gene density and chromosome evolution
#' 
#' Tests how gene density correlates with chromosome numbers and evolution rates,
#' accounting for phylogenetic relationships
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param gene_counts Named vector of gene counts
#' @param genome_sizes Named vector of genome sizes (Mb)
#' @param ancestral_recon Optional ancestral reconstruction results
#' @param rate_results Optional rate analysis results
#' @param method Analysis method: "PGLS", "pic", "standard"
#' @return List with gene density analysis results
#' @export
analyze_gene_density <- function(tree,
                               chr_counts,
                               gene_counts,
                               genome_sizes,
                               ancestral_recon = NULL,
                               rate_results = NULL,
                               method = "PGLS") {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  if(!is.vector(gene_counts) || is.null(names(gene_counts))) {
    stop("gene_counts must be a named vector")
  }
  
  if(!is.vector(genome_sizes) || is.null(names(genome_sizes))) {
    stop("genome_sizes must be a named vector")
  }
  
  # Find species with all required data
  common_species <- Reduce(intersect, list(
    names(chr_counts),
    names(gene_counts),
    names(genome_sizes),
    tree$tip.label
  ))
  
  if(length(common_species) < 4) {
    stop("At least 4 species must have chromosome counts, gene counts, genome sizes, and be present in the tree")
  }
  
  # Extract matched data
  matched_chr <- chr_counts[common_species]
  matched_genes <- gene_counts[common_species]
  matched_sizes <- genome_sizes[common_species]
  
  # Calculate gene density (genes per Mb)
  gene_density <- matched_genes / matched_sizes
  names(gene_density) <- common_species
  
  # Calculate genes per chromosome
  genes_per_chr <- matched_genes / matched_chr
  names(genes_per_chr) <- common_species
  
  # Prune tree to common species
  pruned_tree <- ape::keep.tip(tree, common_species)
  
  # Initialize results
  results <- list(
    gene_density_correlation = NULL,
    genes_per_chr_correlation = NULL,
    summary = list(
      n_species = length(common_species),
      method = method
    ),
    plots = list()
  )
  
  # Analyze correlation between gene density and chromosome number
  density_result <- NULL
  
  if(method == "PGLS") {
    # Phylogenetic Generalized Least Squares
    if(requireNamespace("caper", quietly = TRUE)) {
      # Create comparative data
      data <- data.frame(
        species = common_species,
        chr = matched_chr,
        density = gene_density
      )
      
      comp_data <- caper::comparative.data(pruned_tree, data, "species")
      
      # Run PGLS
      pgls_model <- tryCatch({
        caper::pgls(chr ~ density, data = comp_data, lambda = "ML")
      }, error = function(e) {
        warning(paste("PGLS error for gene density:", e$message))
        return(NULL)
      })
      
      if(!is.null(pgls_model)) {
        coef <- summary(pgls_model)$coefficients["density", "Estimate"]
        p_value <- summary(pgls_model)$coefficients["density", "Pr(>|t|)"]
        
        # Calculate confidence intervals
        ci <- confint(pgls_model)["density", ]
        
        density_result <- list(
          coefficient = coef,
          p_value = p_value,
          ci_lower = ci[1],
          ci_upper = ci[2],
          lambda = pgls_model$param$lambda,
          model = pgls_model
        )
      }
    } else {
      warning("caper package required for PGLS analysis. Using standard correlation instead.")
      method <- "standard"
    }
  } else if(method == "pic") {
    # Phylogenetic Independent Contrasts
    if(requireNamespace("ape", quietly = TRUE)) {
      # Calculate PICs for both variables
      chr_pic <- ape::pic(matched_chr, pruned_tree)
      density_pic <- ape::pic(gene_density, pruned_tree)
      
      # Test correlation of contrasts
      cor_test <- cor.test(chr_pic, density_pic, method = "pearson")
      
      density_result <- list(
        coefficient = cor_test$estimate,
        p_value = cor_test$p.value,
        ci_lower = cor_test$conf.int[1],
        ci_upper = cor_test$conf.int[2],
        test = cor_test
      )
    } else {
      warning("ape package required for PIC analysis. Using standard correlation instead.")
      method <- "standard"
    }
  }
  
  if(method == "standard" || is.null(density_result)) {
    # Standard Pearson correlation
    cor_test <- cor.test(matched_chr, gene_density, method = "pearson")
    
    density_result <- list(
      coefficient = cor_test$estimate,
      p_value = cor_test$p.value,
      ci_lower = cor_test$conf.int[1],
      ci_upper = cor_test$conf.int[2],
      test = cor_test
    )
  }
  
  # Store gene density correlation results
  results$gene_density_correlation <- density_result
  
  # Analyze correlation between genes per chromosome and chromosome number
  per_chr_result <- NULL
  
  # Similar structure as above, but for genes per chromosome vs chromosome number
  if(method == "PGLS") {
    if(requireNamespace("caper", quietly = TRUE)) {
      data <- data.frame(
        species = common_species,
        chr = matched_chr,
        genes_per_chr = genes_per_chr
      )
      
      comp_data <- caper::comparative.data(pruned_tree, data, "species")
      
      pgls_model <- tryCatch({
        caper::pgls(chr ~ genes_per_chr, data = comp_data, lambda = "ML")
      }, error = function(e) {
        warning(paste("PGLS error for genes per chromosome:", e$message))
        return(NULL)
      })
      
      if(!is.null(pgls_model)) {
        coef <- summary(pgls_model)$coefficients["genes_per_chr", "Estimate"]
        p_value <- summary(pgls_model)$coefficients["genes_per_chr", "Pr(>|t|)"]
        
        ci <- confint(pgls_model)["genes_per_chr", ]
        
        per_chr_result <- list(
          coefficient = coef,
          p_value = p_value,
          ci_lower = ci[1],
          ci_upper = ci[2],
          lambda = pgls_model$param$lambda,
          model = pgls_model
        )
      }
    }
  } else if(method == "pic") {
    if(requireNamespace("ape", quietly = TRUE)) {
      chr_pic <- ape::pic(matched_chr, pruned_tree)
      per_chr_pic <- ape::pic(genes_per_chr, pruned_tree)
      
      cor_test <- cor.test(chr_pic, per_chr_pic, method = "pearson")
      
      per_chr_result <- list(
        coefficient = cor_test$estimate,
        p_value = cor_test$p.value,
        ci_lower = cor_test$conf.int[1],
        ci_upper = cor_test$conf.int[2],
        test = cor_test
      )
    }
  } else {
    cor_test <- cor.test(matched_chr, genes_per_chr, method = "pearson")
    
    per_chr_result <- list(
      coefficient = cor_test$estimate,
      p_value = cor_test$p.value,
      ci_lower = cor_test$conf.int[1],
      ci_upper = cor_test$conf.int[2],
      test = cor_test
    )
  }
  
  # Store genes per chromosome correlation results
  results$genes_per_chr_correlation <- per_chr_result
  
  # Create scatter plots if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Create data frame for plotting
    plot_data <- data.frame(
      Species = common_species,
      Chromosome_Number = matched_chr,
      Gene_Density = gene_density,
      Genes_Per_Chr = genes_per_chr,
      stringsAsFactors = FALSE
    )
    
    # Gene density plot
    density_plot <- ggplot2::ggplot(plot_data, 
                                  ggplot2::aes(x = Gene_Density, y = Chromosome_Number)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::geom_smooth(method = "lm", formula = 'y ~ x', se = TRUE, color = "blue") +
      ggplot2::labs(
        title = "Chromosome Number vs. Gene Density",
        subtitle = paste(
          "Correlation =", round(density_result$coefficient, 3),
          "(p =", format(density_result$p_value, digits = 3), ")"
        ),
        x = "Gene Density (genes/Mb)",
        y = "Chromosome Number"
      ) +
      ggplot2::theme_minimal()
    
    # Genes per chromosome plot
    per_chr_plot <- ggplot2::ggplot(plot_data, 
                                  ggplot2::aes(x = Genes_Per_Chr, y = Chromosome_Number)) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::geom_smooth(method = "lm", formula = 'y ~ x', se = TRUE, color = "blue") +
      ggplot2::labs(
        title = "Chromosome Number vs. Genes Per Chromosome",
        subtitle = paste(
          "Correlation =", round(per_chr_result$coefficient, 3),
          "(p =", format(per_chr_result$p_value, digits = 3), ")"
        ),
        x = "Genes Per Chromosome",
        y = "Chromosome Number"
      ) +
      ggplot2::theme_minimal()
    
    # Add species labels if not too many
    if(length(common_species) <= 20) {
      density_plot <- density_plot +
        ggplot2::geom_text(ggplot2::aes(label = Species), 
                         nudge_x = 0.05 * max(Gene_Density, na.rm = TRUE), 
                         nudge_y = 0.05 * max(Chromosome_Number, na.rm = TRUE), 
                         size = 3)
      
      per_chr_plot <- per_chr_plot +
        ggplot2::geom_text(ggplot2::aes(label = Species), 
                         nudge_x = 0.05 * max(Genes_Per_Chr, na.rm = TRUE), 
                         nudge_y = 0.05 * max(Chromosome_Number, na.rm = TRUE), 
                         size = 3)
    }
    
    results$plots$density_plot <- density_plot
    results$plots$per_chr_plot <- per_chr_plot
    
    # Combined plot if patchwork is available
    if(requireNamespace("patchwork", quietly = TRUE)) {
      results$plots$combined <- density_plot + per_chr_plot +
        patchwork::plot_layout(ncol = 2)
    }
  }
  
  # Create summary
  results$summary$gene_density_cor <- density_result$coefficient
  results$summary$gene_density_p <- density_result$p_value
  results$summary$gene_density_significant <- density_result$p_value < 0.05
  
  results$summary$genes_per_chr_cor <- per_chr_result$coefficient
  results$summary$genes_per_chr_p <- per_chr_result$p_value
  results$summary$genes_per_chr_significant <- per_chr_result$p_value < 0.05
  
  # Interpret results
  results$summary$gene_density_interpretation <- if(density_result$p_value < 0.05) {
    if(density_result$coefficient > 0) {
      "Significant positive correlation between gene density and chromosome number"
    } else {
      "Significant negative correlation between gene density and chromosome number"
    }
  } else {
    "No significant correlation between gene density and chromosome number"
  }
  
  results$summary$genes_per_chr_interpretation <- if(per_chr_result$p_value < 0.05) {
    if(per_chr_result$coefficient > 0) {
      "Significant positive correlation between genes per chromosome and chromosome number"
    } else {
      "Significant negative correlation between genes per chromosome and chromosome number"
    }
  } else {
    "No significant correlation between genes per chromosome and chromosome number"
  }
  
  return(results)
}

#===============================================================================
# Integrated Workflow Function
#===============================================================================

#' Run comprehensive genomic integration workflow
#' 
#' Performs a series of analyses exploring relationships between chromosome 
#' evolution and various genomic features
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param genome_features Data frame of genomic features
#' @param ancestral_recon Optional ancestral reconstruction results
#' @param repeat_data Optional data frame with repeat content data
#' @param synteny_data Optional data frame with synteny conservation data
#' @param output_dir Directory to save results and visualizations
#' @param analyses Vector of analyses to perform
#' @param method Correlation method: "PGLS", "pic", "standard"
#' @param create_plots Whether to generate and save plots
#' @param verbose Whether to print progress messages
#' @return List with results from all analyses
#' @export
run_genomic_integration_workflow <- function(tree,
                                          chr_counts,
                                          genome_features,
                                          ancestral_recon = NULL,
                                          repeat_data = NULL,
                                          synteny_data = NULL,
                                          output_dir = "genomic_integration_results",
                                          analyses = c("correlations", "synteny", "repeats", "gene_density"),
                                          method = "PGLS",
                                          create_plots = TRUE,
                                          verbose = TRUE) {
  # Validate inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  if(!is.data.frame(genome_features)) {
    stop("genome_features must be a data frame")
  }
  
  # Create output directory if it doesn't exist
  if(create_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    
    # Create subdirectory for plots
    plots_dir <- file.path(output_dir, "plots")
    if(!dir.exists(plots_dir)) {
      dir.create(plots_dir)
    }
  }
  
  # Initialize results
  results <- list(
    metadata = list(
      date = Sys.Date(),
      analyses_run = analyses,
      method = method
    ),
    data = list(
      tree = tree,
      chr_counts = chr_counts,
      genome_features = genome_features,
      ancestral_recon = ancestral_recon
    )
  )
  
  # Integrate data
  if(verbose) message("Integrating genomic data...")
  
  integrated_data <- integrate_genomic_data(
    tree = tree,
    chr_counts = chr_counts,
    genome_features = genome_features,
    check_tree_matching = TRUE,
    impute_missing = "phylo",
    scale_features = FALSE
  )
  
  results$integrated_data <- integrated_data
  
  # Run genomic feature correlations analysis
  if("correlations" %in% analyses) {
    if(verbose) message("Analyzing correlations with genomic features...")
    
    corr_results <- analyze_genomic_correlations(
      integrated_data = integrated_data,
      method = method,
      features = NULL,
      model = "BM",
      p_adjustment = "fdr"
    )
    
    results$feature_correlations <- corr_results
    
    # Save correlation results
    saveRDS(corr_results, file.path(output_dir, "feature_correlations.rds"))
    
    # Create and save plots
    if(create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
      if(!is.null(corr_results$plots$correlation_plot)) {
        ggplot2::ggsave(
          file.path(output_dir, "plots", "feature_correlations.pdf"),
          corr_results$plots$correlation_plot,
          width = 10,
          height = 6,
          dpi = 300
        )
      }
      
      # Save scatter plots for significant correlations
      if(!is.null(corr_results$plots$scatter_plots)) {
        for(feature in names(corr_results$plots$scatter_plots)) {
          ggplot2::ggsave(
            file.path(output_dir, "plots", paste0("scatter_", gsub("[^a-zA-Z0-9]", "_", feature), ".pdf")),
            corr_results$plots$scatter_plots[[feature]],
            width = 8,
            height = 6,
            dpi = 300
          )
        }
      }
    }
  }
  
  # Run synteny conservation analysis
  if("synteny" %in% analyses && !is.null(synteny_data)) {
    if(verbose) message("Analyzing synteny conservation...")
    
    synteny_results <- analyze_synteny_conservation(
      tree = tree,
      chr_counts = chr_counts,
      synteny_data = synteny_data,
      synteny_metric = colnames(synteny_data)[3], # Assumes the third column is the synteny metric
      metric_type = "similarity",
      ancestral_recon = ancestral_recon,
      permutations = 1000
    )
    
    results$synteny_analysis <- synteny_results
    
    # Save synteny results
    saveRDS(synteny_results, file.path(output_dir, "synteny_analysis.rds"))
    
    # Create and save plots
    if(create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
      if(!is.null(synteny_results$plots$distance_plot)) {
        ggplot2::ggsave(
          file.path(output_dir, "plots", "synteny_phylodist.pdf"),
          synteny_results$plots$distance_plot,
          width = 8,
          height = 6,
          dpi = 300
        )
      }
      
      if(!is.null(synteny_results$plots$chr_diff_plot)) {
        ggplot2::ggsave(
          file.path(output_dir, "plots", "synteny_chrdiff.pdf"),
          synteny_results$plots$chr_diff_plot,
          width = 8,
          height = 6,
          dpi = 300
        )
      }
    }
  }
  
  # Run repeat content analysis
  if("repeats" %in% analyses && !is.null(repeat_data)) {
    if(verbose) message("Analyzing repeat content relationships...")
    
    # Try to get genome sizes if available in genome_features
    genome_sizes <- NULL
    if("genome_size" %in% colnames(genome_features)) {
      genome_sizes <- genome_features[, "genome_size"]
      names(genome_sizes) <- rownames(genome_features)
    }
    
    repeat_results <- analyze_repeat_content(
      tree = tree,
      chr_counts = chr_counts,
      repeat_data = repeat_data,
      genome_sizes = genome_sizes,
      method = method,
      normalize = !is.null(genome_sizes),
      p_adjustment = "fdr"
    )
    
    results$repeat_analysis <- repeat_results
    
    # Save repeat analysis results
    saveRDS(repeat_results, file.path(output_dir, "repeat_analysis.rds"))
    
    # Create and save plots
    if(create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
      if(!is.null(repeat_results$plots$correlation_plot)) {
        ggplot2::ggsave(
          file.path(output_dir, "plots", "repeat_correlations.pdf"),
          repeat_results$plots$correlation_plot,
          width = 9,
          height = 6,
          dpi = 300
        )
      }
      
      if(!is.null(repeat_results$plots$scatter_plot)) {
        ggplot2::ggsave(
          file.path(output_dir, "plots", "repeat_scatterplots.pdf"),
          repeat_results$plots$scatter_plot,
          width = 10,
          height = 8,
          dpi = 300
        )
      }
    }
  }
  
  # Run gene density analysis
  if("gene_density" %in% analyses) {
    # Check if required data is available in genome_features
    has_gene_data <- all(c("gene_count", "genome_size") %in% colnames(genome_features))
    
    if(has_gene_data) {
      if(verbose) message("Analyzing gene density relationships...")
      
      gene_counts <- genome_features[, "gene_count"]
      names(gene_counts) <- rownames(genome_features)
      
      genome_sizes <- genome_features[, "genome_size"]
      names(genome_sizes) <- rownames(genome_features)
      
      gene_results <- analyze_gene_density(
        tree = tree,
        chr_counts = chr_counts,
        gene_counts = gene_counts,
        genome_sizes = genome_sizes,
        ancestral_recon = ancestral_recon,
        method = method
      )
      
      results$gene_density_analysis <- gene_results
      
      # Save gene density analysis results
      saveRDS(gene_results, file.path(output_dir, "gene_density_analysis.rds"))
      
      # Create and save plots
      if(create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
        if(!is.null(gene_results$plots$combined)) {
          ggplot2::ggsave(
            file.path(output_dir, "plots", "gene_density_analysis.pdf"),
            gene_results$plots$combined,
            width = 12,
            height = 6,
            dpi = 300
          )
        } else {
          if(!is.null(gene_results$plots$density_plot)) {
            ggplot2::ggsave(
              file.path(output_dir, "plots", "gene_density_plot.pdf"),
              gene_results$plots$density_plot,
              width = 8,
              height = 6,
              dpi = 300
            )
          }
          
          if(!is.null(gene_results$plots$per_chr_plot)) {
            ggplot2::ggsave(
              file.path(output_dir, "plots", "genes_per_chr_plot.pdf"),
              gene_results$plots$per_chr_plot,
              width = 8,
              height = 6,
              dpi = 300
            )
          }
        }
      }
    } else {
      warning("Gene density analysis requires 'gene_count' and 'genome_size' columns in genome_features")
    }
  }
  
  # Create integrated visualization
  if(create_plots && requireNamespace("ggplot2", quietly = TRUE) && 
     requireNamespace("ggtree", quietly = TRUE) && !is.null(ancestral_recon)) {
    
    if(verbose) message("Creating integrated visualization...")
    
    # Identify most correlated features
    feature_cors <- NULL
    if("correlations" %in% analyses && !is.null(results$feature_correlations$correlation_table)) {
      feature_cors <- results$feature_correlations$correlation_table
      selected_features <- feature_cors$Feature[order(-abs(feature_cors$Coefficient))][1:min(3, nrow(feature_cors))]
    } else {
      # Select first 3 numeric features
      numeric_features <- colnames(genome_features)[sapply(genome_features, is.numeric)]
      selected_features <- head(numeric_features, 3)
    }
    
    # Create integrated plot
    integrated_viz <- visualize_integrated_evolution(
      integrated_data = integrated_data,
      ancestral_states = ancestral_recon,
      features = selected_features,
      highlight_nodes = NULL,
      node_labels = TRUE,
      use_patchwork = TRUE,
      output_file = file.path(output_dir, "plots", "integrated_visualization.pdf")
    )
    
    results$integrated_visualization <- integrated_viz
  }
  
  # Generate faceted scatter plots of all genomic features
  if(create_plots && requireNamespace("ggplot2", quietly = TRUE)) {
    if(verbose) message("Creating feature relationship plots...")
    
    facet_plot <- plot_chr_feature_relationships(
      integrated_data = integrated_data,
      features = NULL,
      add_trendlines = TRUE,
      facet_ncol = 2,
      output_file = file.path(output_dir, "plots", "feature_relationships.pdf")
    )
    
    results$feature_relationship_plot <- facet_plot
  }
  
  # Create summary of all results
  results$summary <- create_genomic_integration_summary(results)
  
  # Save summary as a text file
  if(!is.null(results$summary) && !is.null(results$summary$text)) {
    writeLines(results$summary$text, file.path(output_dir, "genomic_integration_summary.txt"))
  }
  
  # Save complete results
  saveRDS(results, file.path(output_dir, "genomic_integration_workflow_results.rds"))
  
  return(results)
}

#' Create summary of genomic integration results
#' 
#' @param results Results from genomic integration analyses
#' @return List with summary information
#' @keywords internal
create_genomic_integration_summary <- function(results) {
  # Initialize text summary
  summary_text <- c(
    "GENOMIC INTEGRATION ANALYSIS SUMMARY",
    "=====================================",
    "",
    paste("Analysis date:", format(Sys.Date(), "%Y-%m-%d")),
    paste("Analyses performed:", paste(results$metadata$analyses_run, collapse = ", ")),
    paste("Correlation method used:", results$metadata$method),
    paste("Number of taxa analyzed:", results$integrated_data$n_taxa),
    paste("Number of genomic features:", ncol(results$integrated_data$genome_features)),
    "",
    "KEY FINDINGS:",
    "-------------",
    ""
  )
  
  # Add feature correlation summary
  if(!is.null(results$feature_correlations)) {
    if(!is.null(results$feature_correlations$summary)) {
      n_significant <- results$feature_correlations$summary$n_significant
      summary_text <- c(summary_text, 
                      paste("GENOMIC FEATURE CORRELATIONS:"),
                      paste("*", n_significant, "significant correlations found between chromosome numbers and genomic features"))
      
      # Add top correlations if available
      if(n_significant > 0 && !is.null(results$feature_correlations$correlation_table)) {
        top_corr <- head(results$feature_correlations$correlation_table[results$feature_correlations$correlation_table$significant, ], 3)
        for(i in 1:nrow(top_corr)) {
          dir <- ifelse(top_corr$Coefficient[i] > 0, "positive", "negative")
          summary_text <- c(summary_text, 
                          paste("  -", top_corr$Feature[i], "shows", dir, "correlation (r =", 
                              round(top_corr$Coefficient[i], 3), ", p =", 
                              format(top_corr$adjusted_p[i], digits = 3), ")"))
        }
      }
      summary_text <- c(summary_text, "")
    }
  }
  
  # Add synteny analysis summary
  if(!is.null(results$synteny_analysis)) {
    if(!is.null(results$synteny_analysis$summary)) {
      summary_text <- c(summary_text, 
                      paste("SYNTENY CONSERVATION:"),
                      paste("* Synteny metric:", results$synteny_analysis$summary$synteny_metric),
                      paste("* Synteny vs phylogenetic distance: correlation =", 
                          round(results$synteny_analysis$summary$synteny_distance_correlation, 3),
                          "(p =", format(results$synteny_analysis$summary$synteny_distance_p, digits = 3), ")"),
                      paste("* Synteny vs chromosome differences: correlation =", 
                          round(results$synteny_analysis$summary$synteny_chr_correlation, 3),
                          "(p =", format(results$synteny_analysis$summary$synteny_chr_p, digits = 3), ")"),
                      "")
    }
  }
  
  # Add repeat content analysis summary
  if(!is.null(results$repeat_analysis)) {
    if(!is.null(results$repeat_analysis$summary)) {
      n_significant <- results$repeat_analysis$summary$n_significant
      summary_text <- c(summary_text, 
                      paste("REPEAT CONTENT ANALYSIS:"),
                      paste("*", n_significant, "repeat classes significantly correlated with chromosome number"))
      
      # Add top correlations if available
      if(n_significant > 0 && !is.null(results$repeat_analysis$correlation_table)) {
        top_rep <- head(results$repeat_analysis$correlation_table[results$repeat_analysis$correlation_table$significant, ], 3)
        for(i in 1:nrow(top_rep)) {
          dir <- ifelse(top_rep$Coefficient[i] > 0, "positive", "negative")
          summary_text <- c(summary_text, 
                          paste("  -", top_rep$Repeat_Class[i], "shows", dir, "correlation (r =", 
                              round(top_rep$Coefficient[i], 3), ", p =", 
                              format(top_rep$adjusted_p[i], digits = 3), ")"))
        }
      }
      summary_text <- c(summary_text, "")
    }
  }
  
  # Add gene density analysis summary
  if(!is.null(results$gene_density_analysis)) {
    if(!is.null(results$gene_density_analysis$summary)) {
      summary_text <- c(summary_text, 
                      paste("GENE DENSITY ANALYSIS:"),
                      paste("* Gene density vs chromosome number: correlation =", 
                          round(results$gene_density_analysis$summary$gene_density_cor, 3),
                          "(p =", format(results$gene_density_analysis$summary$gene_density_p, digits = 3), ")"),
                      paste("* Genes per chromosome vs chromosome number: correlation =", 
                          round(results$gene_density_analysis$summary$genes_per_chr_cor, 3),
                          "(p =", format(results$gene_density_analysis$summary$genes_per_chr_p, digits = 3), ")"),
                      paste("* Interpretation:", results$gene_density_analysis$summary$gene_density_interpretation),
                      "")
    }
  }
  
  # Add overall conclusion
  summary_text <- c(summary_text,
                  "OVERALL CONCLUSIONS:",
                  "------------------")
  
  # Generate overall conclusion based on findings
  has_significant_findings <- FALSE
  
  if(!is.null(results$feature_correlations$summary$n_significant) && 
     results$feature_correlations$summary$n_significant > 0) {
    has_significant_findings <- TRUE
  }
  
  if(!is.null(results$synteny_analysis$summary$synteny_chr_p) && 
     results$synteny_analysis$summary$synteny_chr_p < 0.05) {
    has_significant_findings <- TRUE
  }
  
  if(!is.null(results$repeat_analysis$summary$n_significant) && 
     results$repeat_analysis$summary$n_significant > 0) {
    has_significant_findings <- TRUE
  }
  
  if(!is.null(results$gene_density_analysis$summary$gene_density_significant) && 
     results$gene_density_analysis$summary$gene_density_significant) {
    has_significant_findings <- TRUE
  }
  
  if(has_significant_findings) {
    summary_text <- c(summary_text,
                    "The analyses reveal significant relationships between chromosome numbers and various genomic features,",
                    "suggesting that chromosome evolution is linked to changes in genomic architecture in this group.",
                    "These findings support the hypothesis that chromosome restructuring is associated with broader",
                    "genomic changes rather than occurring in isolation.")
  } else {
    summary_text <- c(summary_text,
                    "The analyses did not find strong relationships between chromosome numbers and the examined genomic features.",
                    "This suggests that chromosome evolution in this group may be decoupled from changes in the examined",
                    "genomic features, possibly being driven by other factors not included in the current analysis.")
  }
  
  # Return summary
  return(list(
    text = summary_text,
    has_significant_findings = has_significant_findings
  ))
}
