#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Chromosome Data Handler Module
# Author: Bioinformatics Team
# Date: 2025-05-05
# Description: Functions for processing, validating, and preparing chromosome
#              count data for ancestral state reconstruction and analysis
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(parallel)
})

#===============================================================================
# Data Loading and Validation Functions
#===============================================================================

#' Load chromosome count data from different formats
#' 
#' Reads and parses chromosome count data from various file formats
#' 
#' @param file Path to the data file
#' @param format Format of the data file: "table", "nexus", "csv", or "auto" (detect)
#' @param species_col Column containing species names (for table format)
#' @param count_col Column containing chromosome counts (for table format)
#' @param sep Field separator for tabular data (default tab)
#' @param header Whether the file has a header row
#' @param clean_names Whether to clean species names (remove special characters)
#' @param verbose Whether to print messages during loading
#' @return Named vector of chromosome counts with species names
#' @export
load_chromosome_data <- function(file, 
                               format = "auto", 
                               species_col = 1, 
                               count_col = 2,
                               sep = "\t",
                               header = TRUE,
                               clean_names = TRUE,
                               verbose = TRUE) {
  # Check if file exists
  if(!file.exists(file)) {
    stop(paste("File not found:", file))
  }
  
  # Auto-detect format if needed
  if(format == "auto") {
    ext <- tolower(tools::file_ext(file))
    if(ext == "nex" || ext == "nexus") {
      format <- "nexus"
    } else if(ext == "csv") {
      format <- "csv"
      sep <- ","
    } else if(ext == "tsv" || ext == "txt") {
      format <- "table"
    } else {
      # Default to table format
      format <- "table"
    }
    if(verbose) message(paste("Auto-detected format:", format))
  }
  
  # Load data based on format
  chr_counts <- NULL
  
  if(format == "table" || format == "csv") {
    # Read tabular data
    data <- utils::read.table(file, 
                             sep = sep,
                             header = header,
                             stringsAsFactors = FALSE,
                             quote = "\"")
    
    # Extract species names and counts
    if(header) {
      if(is.character(species_col) && !species_col %in% colnames(data)) {
        stop(paste("Species column not found:", species_col))
      }
      if(is.character(count_col) && !count_col %in% colnames(data)) {
        stop(paste("Count column not found:", count_col))
      }
    } else {
      if(is.character(species_col)) {
        stop("Cannot use column names with header=FALSE")
      }
    }
    
    # Extract data
    species <- data[[species_col]]
    counts <- data[[count_col]]
    
    # Clean species names if requested
    if(clean_names) {
      species <- clean_species_names(species)
    }
    
    # Create named vector
    chr_counts <- counts
    names(chr_counts) <- species
    
  } else if(format == "nexus") {
    # Read Nexus format
    if(requireNamespace("ape", quietly = TRUE)) {
      # Try to read using APE
      nexus_data <- ape::read.nexus.data(file)
      
      # NEXUS typically stores discrete characters, extract first character set
      if(length(nexus_data) > 0) {
        # Convert character data to numeric
        chr_counts <- sapply(nexus_data, function(x) as.numeric(x[1]))
      } else {
        stop("No data found in NEXUS file")
      }
    } else {
      stop("ape package is required to read NEXUS files")
    }
  } else {
    stop(paste("Unsupported format:", format))
  }
  
  # Check for valid data
  if(is.null(chr_counts) || length(chr_counts) == 0) {
    stop("No chromosome count data loaded")
  }
  
  # Check for valid counts
  invalid_counts <- is.na(chr_counts) | !is.numeric(chr_counts) | chr_counts <= 0
  if(any(invalid_counts)) {
    if(verbose) {
      warning(paste("Found", sum(invalid_counts), "invalid chromosome counts. Converting to NA."))
    }
    chr_counts[invalid_counts] <- NA
  }
  
  # Print summary
  if(verbose) {
    message(sprintf("Loaded chromosome counts for %d species (%.1f%% complete)",
                  length(chr_counts),
                  100 * sum(!is.na(chr_counts)) / length(chr_counts)))
  }
  
  return(chr_counts)
}

#' Clean species names by standardizing format
#' 
#' @param species Vector of species names
#' @param to_title_case Convert to title case (first letter capitalized)
#' @param remove_special Remove special characters
#' @param hybrid_pattern Pattern to identify hybrid names
#' @return Vector of cleaned species names
#' @keywords internal
clean_species_names <- function(species, 
                               to_title_case = TRUE, 
                               remove_special = TRUE,
                               hybrid_pattern = " x ") {
  # Remove leading/trailing whitespace
  species <- trimws(species)
  
  # Convert to title case if requested
  if(to_title_case) {
    # Split genus and species
    parts <- strsplit(species, " ")
    
    # Process each name
    species <- sapply(parts, function(p) {
      # Capitalize genus
      if(length(p) >= 1) {
        p[1] <- paste0(toupper(substr(p[1], 1, 1)), substring(p[1], 2))
      }
      
      # Keep species epithet lowercase
      if(length(p) >= 2) {
        p[2] <- tolower(p[2])
      }
      
      # Join parts back together
      paste(p, collapse = " ")
    })
  }
  
  # Handle hybrid symbols
  if(!is.null(hybrid_pattern)) {
    # Replace various hybrid patterns with standard format
    species <- gsub("\\s*[×xX]\\s*", hybrid_pattern, species)
  }
  
  # Remove special characters if requested
  if(remove_special) {
    # Keep spaces, letters, numbers, periods, hyphens, and underscores
    # Replace other special characters with underscores
    species <- gsub("[^a-zA-Z0-9\\.\\-_ ]", "_", species)
  }
  
  return(species)
}

#' Validate chromosome count data against a tree
#' 
#' Checks compatibility between chromosome count data and a phylogenetic tree
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param strict Whether to error on missing data (TRUE) or just warn (FALSE)
#' @param match_tip_format Whether to try matching tip labels with species names
#' @param check_numeric Check if all counts are numeric
#' @param check_positive Check if all counts are positive
#' @param verbose Whether to print messages during validation
#' @return List with validation results and possibly corrected data
#' @export
validate_chromosome_data <- function(tree, 
                                   chr_counts, 
                                   strict = FALSE,
                                   match_tip_format = TRUE,
                                   check_numeric = TRUE,
                                   check_positive = TRUE,
                                   verbose = TRUE) {
  # Check input types
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts)) {
    stop("chr_counts must be a vector")
  }
  
  # Initialize results
  results <- list(
    valid = FALSE,
    species_in_tree = 0,
    species_in_data = length(chr_counts),
    species_match = 0,
    species_missing = 0,
    tips_missing_data = 0,
    data_not_in_tree = 0,
    corrected_data = NULL,
    tree = tree,
    issues = character(0)
  )
  
  # Get species lists
  tree_species <- tree$tip.label
  data_species <- names(chr_counts)
  
  # Count species in tree
  results$species_in_tree <- length(tree_species)
  
  # Check if data has names
  if(is.null(data_species)) {
    msg <- "Chromosome count data is not named with species labels"
    results$issues <- c(results$issues, msg)
    if(strict) {
      stop(msg)
    } else {
      warning(msg)
      return(results)
    }
  }
  
  # Find exact matches
  matches <- intersect(tree_species, data_species)
  results$species_match <- length(matches)
  
  # Calculate missing and extra species
  missing_tips <- setdiff(tree_species, data_species)
  results$species_missing <- length(missing_tips)
  
  extra_data <- setdiff(data_species, tree_species)
  results$data_not_in_tree <- length(extra_data)
  
  # Try to match species names if requested
  corrected_matches <- matches
  if(match_tip_format && results$species_missing > 0 && results$data_not_in_tree > 0) {
    # Clean both sets of names for matching
    clean_tree_species <- clean_species_names(tree_species)
    clean_data_species <- clean_species_names(data_species)
    
    # Create mapping from clean names to original names
    tree_name_map <- setNames(tree_species, clean_tree_species)
    data_name_map <- setNames(data_species, clean_data_species)
    
    # Find matches between clean names
    clean_missing_tips <- setdiff(clean_tree_species, clean_data_species)
    clean_extra_data <- setdiff(clean_data_species, clean_tree_species)
    
    # Update match count
    clean_matches <- intersect(clean_tree_species, clean_data_species)
    
    # Check if we found more matches
    if(length(clean_matches) > length(matches)) {
      # Create corrected data
      corrected_data <- chr_counts[data_name_map[clean_matches]]
      names(corrected_data) <- tree_name_map[clean_matches]
      
      # Add to results
      results$corrected_data <- corrected_data
      corrected_matches <- names(corrected_data)
      
      if(verbose) {
        message(sprintf("Found %d additional matches after cleaning species names",
                       length(clean_matches) - length(matches)))
      }
    }
  }
  
  # Count how many tips lack data
  missing_data_tips <- setdiff(tree_species, corrected_matches)
  results$tips_missing_data <- length(missing_data_tips)
  
  # Check data quality if requested
  data_subset <- chr_counts[corrected_matches]
  
  if(check_numeric) {
    non_numeric <- !is.numeric(data_subset)
    if(any(non_numeric)) {
      msg <- paste("Found", sum(non_numeric), "non-numeric chromosome counts")
      results$issues <- c(results$issues, msg)
      if(strict) {
        stop(msg)
      } else {
        warning(msg)
      }
    }
  }
  
  if(check_positive) {
    non_positive <- data_subset <= 0
    if(any(non_positive, na.rm = TRUE)) {
      msg <- paste("Found", sum(non_positive, na.rm = TRUE), "non-positive chromosome counts")
      results$issues <- c(results$issues, msg)
      if(strict) {
        stop(msg)
      } else {
        warning(msg)
      }
    }
  }
  
  # Check for NAs
  na_values <- is.na(data_subset)
  if(any(na_values)) {
    msg <- paste("Found", sum(na_values), "NA chromosome counts")
    results$issues <- c(results$issues, msg)
    if(verbose) {
      warning(msg)
    }
  }
  
  # Determine overall validity
  results$valid <- results$species_match > 0 && 
                  length(results$issues) == 0
  
  # Print summary if requested
  if(verbose) {
    message(sprintf("Validation of chromosome data against tree: %s",
                  ifelse(results$valid, "passed", "failed")))
    message(sprintf("  Species in tree: %d", results$species_in_tree))
    message(sprintf("  Species in data: %d", results$species_in_data))
    message(sprintf("  Species matched: %d", results$species_match))
    message(sprintf("  Tree tips missing data: %d", results$tips_missing_data))
    message(sprintf("  Data species not in tree: %d", results$data_not_in_tree))
    
    if(length(results$issues) > 0) {
      message("  Issues found:")
      for(issue in results$issues) {
        message(paste("   -", issue))
      }
    }
  }
  
  return(results)
}

#' Match chromosome count data to a tree
#' 
#' Creates a properly formatted named vector of chromosome counts
#' that matches the tip labels in a phylogenetic tree
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param impute Whether to impute missing values
#' @param impute_method Method for imputation: "mean", "median", "phylo", "kNN"
#' @param match_names Whether to attempt name matching between data and tree
#' @param verbose Whether to print messages during processing
#' @return Named vector of chromosome counts matching tree tip labels
#' @export
match_chromosome_data <- function(tree, 
                                chr_counts, 
                                impute = FALSE,
                                impute_method = "phylo",
                                match_names = TRUE,
                                verbose = TRUE) {
  # Validate data first
  validation <- validate_chromosome_data(tree, chr_counts, 
                                       strict = FALSE, 
                                       match_tip_format = match_names,
                                       verbose = verbose)
  
  # Initialize matched data
  matched_data <- rep(NA_real_, length(tree$tip.label))
  names(matched_data) <- tree$tip.label
  
  # If validation returned corrected data, use that
  if(!is.null(validation$corrected_data)) {
    # Fill in corrected matches
    matched_species <- names(validation$corrected_data)
    matched_data[matched_species] <- validation$corrected_data
  } else {
    # Use original matches
    matched_species <- intersect(names(chr_counts), tree$tip.label)
    matched_data[matched_species] <- chr_counts[matched_species]
  }
  
  # Check if we need to impute
  missing_count <- sum(is.na(matched_data))
  
  if(missing_count > 0) {
    if(verbose) {
      message(sprintf("Missing chromosome counts for %d out of %d species (%.1f%%)",
                    missing_count, length(matched_data),
                    100 * missing_count / length(matched_data)))
    }
    
    if(impute) {
      if(verbose) {
        message(sprintf("Imputing missing values using %s method...", impute_method))
      }
      
      # Impute missing values
      imputed_data <- impute_chromosome_counts(tree, matched_data, 
                                            method = impute_method, 
                                            verbose = verbose)
      
      # Update matched data with imputed values
      matched_data <- imputed_data
    } else if(verbose) {
      message("No imputation requested. Returning data with NA values.")
    }
  } else {
    if(verbose) {
      message("All tree tips have matching chromosome count data.")
    }
  }
  
  return(matched_data)
}

#===============================================================================
# Data Imputation Functions
#===============================================================================

#' Impute missing chromosome counts in a dataset
#' 
#' Fills in missing values using various methods
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts, with NAs for missing data
#' @param method Imputation method: "mean", "median", "mode", "phylo", "kNN"
#' @param k Number of neighbors for kNN imputation
#' @param weighted Whether to use distance-weighted imputation for kNN
#' @param round Whether to round imputed values to integers
#' @param verbose Whether to print messages during imputation
#' @return Named vector of chromosome counts with imputed values
#' @export
impute_chromosome_counts <- function(tree, 
                                   chr_counts,
                                   method = "phylo",
                                   k = 5,
                                   weighted = TRUE,
                                   round = TRUE,
                                   verbose = TRUE) {
  # Check input
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Make a copy of the data
  imputed_counts <- chr_counts
  
  # Get indices of missing values
  missing_idx <- which(is.na(imputed_counts))
  
  # If no missing values, return original data
  if(length(missing_idx) == 0) {
    if(verbose) {
      message("No missing values to impute")
    }
    return(imputed_counts)
  }
  
  # Get indices of non-missing values
  present_idx <- which(!is.na(imputed_counts))
  
  # Get species names
  all_species <- names(imputed_counts)
  missing_species <- all_species[missing_idx]
  present_species <- all_species[present_idx]
  
  # Perform imputation based on selected method
  if(method == "mean") {
    # Simple mean imputation
    mean_value <- mean(imputed_counts[present_idx], na.rm = TRUE)
    imputed_counts[missing_idx] <- mean_value
    
    if(verbose) {
      message(sprintf("Imputed %d missing values with mean: %.2f", 
                    length(missing_idx), mean_value))
    }
    
  } else if(method == "median") {
    # Median imputation
    median_value <- median(imputed_counts[present_idx], na.rm = TRUE)
    imputed_counts[missing_idx] <- median_value
    
    if(verbose) {
      message(sprintf("Imputed %d missing values with median: %.2f", 
                    length(missing_idx), median_value))
    }
    
  } else if(method == "mode") {
    # Mode imputation
    tab <- table(imputed_counts[present_idx])
    mode_value <- as.numeric(names(tab)[which.max(tab)])
    imputed_counts[missing_idx] <- mode_value
    
    if(verbose) {
      message(sprintf("Imputed %d missing values with mode: %.2f", 
                    length(missing_idx), mode_value))
    }
    
  } else if(method == "phylo") {
    # Phylogenetic imputation using ancestral state reconstruction
    if(requireNamespace("phytools", quietly = TRUE)) {
      # Extract subtree with non-missing data
      present_tree <- ape::keep.tip(tree, present_species)
      present_data <- imputed_counts[present_species]
      
      # Run ancestral state reconstruction
      ace_result <- ape::ace(present_data, present_tree, type = "continuous", method = "ML")
      
      # For each missing species, find its position relative to the reconstructed tree
      for(i in seq_along(missing_species)) {
        species <- missing_species[i]
        
        # Find closest tip or node
        closest <- find_closest_point(tree, species, present_tree, ace_result)
        
        # Use the reconstructed value
        imputed_counts[species] <- closest$state
        
        if(verbose && i %% 10 == 0) {
          message(sprintf("Imputed %d/%d missing values using phylogenetic method", 
                        i, length(missing_species)))
        }
      }
    } else {
      warning("phytools package required for phylogenetic imputation. Using median instead.")
      median_value <- median(imputed_counts[present_idx], na.rm = TRUE)
      imputed_counts[missing_idx] <- median_value
    }
    
  } else if(method == "kNN") {
    # k-nearest neighbors imputation
    if(requireNamespace("phytools", quietly = TRUE)) {
      # Calculate phylogenetic distance matrix
      dist_matrix <- ape::cophenetic.phylo(tree)
      
      # For each missing species, find k nearest neighbors with data
      for(i in seq_along(missing_species)) {
        species <- missing_species[i]
        
        # Get distances to all other species
        distances <- dist_matrix[species, present_species]
        
        # Find k nearest neighbors
        nearest_idx <- order(distances)[1:min(k, length(distances))]
        nearest_species <- present_species[nearest_idx]
        nearest_distances <- distances[nearest_idx]
        nearest_values <- imputed_counts[nearest_species]
        
        # Calculate imputed value
        if(weighted) {
          # Use inverse distance weighting
          weights <- 1 / (nearest_distances + 1e-10)  # Avoid division by zero
          weights <- weights / sum(weights)
          imputed_value <- sum(weights * nearest_values)
        } else {
          # Simple average
          imputed_value <- mean(nearest_values)
        }
        
        # Assign imputed value
        imputed_counts[species] <- imputed_value
        
        if(verbose && i %% 10 == 0) {
          message(sprintf("Imputed %d/%d missing values using kNN method", 
                        i, length(missing_species)))
        }
      }
    } else {
      warning("ape package required for kNN imputation. Using median instead.")
      median_value <- median(imputed_counts[present_idx], na.rm = TRUE)
      imputed_counts[missing_idx] <- median_value
    }
    
  } else {
    stop(paste("Unsupported imputation method:", method))
  }
  
  # Round imputed values if requested
  if(round) {
    imputed_counts[missing_idx] <- round(imputed_counts[missing_idx])
  }
  
  return(imputed_counts)
}

#' Find the closest point in a phylogeny with a known state
#' 
#' @param full_tree Complete phylogenetic tree
#' @param target_species Species for which to find closest point
#' @param subtree Tree containing only species with known states
#' @param ace_result Ancestral character estimation result for subtree
#' @return List with closest point information
#' @keywords internal
find_closest_point <- function(full_tree, target_species, subtree, ace_result) {
  # Get target tip index in full tree
  target_idx <- which(full_tree$tip.label == target_species)
  
  if(length(target_idx) == 0) {
    stop(paste("Species not found in tree:", target_species))
  }
  
  # Get all ancestors of target species
  ancestors <- get_all_ancestors(full_tree, target_idx)
  
  # Find which of these ancestors are in the subtree
  subtree_nodes <- c(subtree$tip.label, 
                     if(!is.null(subtree$node.label)) subtree$node.label else NULL)
  
  # Get node states from ACE result
  node_states <- ace_result$ace
  
  # Find the closest ancestor in the subtree
  min_dist <- Inf
  closest_node <- NULL
  closest_state <- NA
  
  for(anc in ancestors) {
    # Check if this ancestor is in the subtree
    if(anc %in% subtree_nodes) {
      # Get distance to this ancestor
      dist <- get_phylo_distance(full_tree, target_idx, which(subtree_nodes == anc))
      
      if(dist < min_dist) {
        min_dist <- dist
        closest_node <- anc
        
        # Get state for this node
        node_idx <- which(subtree_nodes == anc)
        if(node_idx <= length(subtree$tip.label)) {
          # It's a tip in the subtree
          closest_state <- ace_result$observed[node_idx]
        } else {
          # It's an internal node
          closest_state <- node_states[node_idx - length(subtree$tip.label)]
        }
      }
    }
  }
  
  # If no ancestor found in subtree, use closest tip
  if(is.null(closest_node)) {
    # Calculate distances to all tips in subtree
    distances <- numeric(length(subtree$tip.label))
    for(i in seq_along(subtree$tip.label)) {
      tip <- subtree$tip.label[i]
      tip_idx <- which(full_tree$tip.label == tip)
      distances[i] <- get_phylo_distance(full_tree, target_idx, tip_idx)
    }
    
    # Find closest tip
    closest_idx <- which.min(distances)
    closest_node <- subtree$tip.label[closest_idx]
    closest_state <- ace_result$observed[closest_idx]
    min_dist <- distances[closest_idx]
  }
  
  return(list(
    node = closest_node,
    state = closest_state,
    distance = min_dist
  ))
}

#' Get all ancestors of a node in a phylogenetic tree
#' 
#' @param tree Phylogenetic tree
#' @param node_idx Index of the node
#' @return Vector of ancestor indices
#' @keywords internal
get_all_ancestors <- function(tree, node_idx) {
  ancestors <- integer(0)
  
  # Find edges where this node is the child
  edge_idx <- which(tree$edge[, 2] == node_idx)
  
  if(length(edge_idx) > 0) {
    # Get parent of this node
    parent_idx <- tree$edge[edge_idx, 1]
    
    # Add parent to ancestors
    ancestors <- c(ancestors, parent_idx)
    
    # Recursively get ancestors of parent
    ancestors <- c(ancestors, get_all_ancestors(tree, parent_idx))
  }
  
  return(ancestors)
}

#' Calculate phylogenetic distance between two nodes
#' 
#' @param tree Phylogenetic tree
#' @param node1 Index of first node
#' @param node2 Index of second node
#' @return Phylogenetic distance
#' @keywords internal
get_phylo_distance <- function(tree, node1, node2) {
  # Use cophenetic distance matrix if available
  if(requireNamespace("ape", quietly = TRUE)) {
    # Calculate distance matrix
    dist_matrix <- ape::cophenetic.phylo(tree)
    
    # Get node labels
    node1_label <- if(node1 <= length(tree$tip.label)) tree$tip.label[node1] else node1
    node2_label <- if(node2 <= length(tree$tip.label)) tree$tip.label[node2] else node2
    
    # Return distance
    return(dist_matrix[node1_label, node2_label])
  } else {
    # Approximate with mean branch length
    return(mean(tree$edge.length) * 2)
  }
}

#===============================================================================
# Data Manipulation Functions
#===============================================================================

#' Standardize chromosome count data
#' 
#' Applies transformations and standardizations to chromosome count data
#' 
#' @param chr_counts Named vector of chromosome counts
#' @param transform Transformation to apply: "none", "log", "sqrt", "rank"
#' @param discrete Whether to treat values as discrete (integers)
#' @param standardize Whether to standardize values (mean=0, sd=1)
#' @param min_value Minimum valid chromosome count
#' @param verbose Whether to print messages during processing
#' @return Standardized chromosome count data
#' @export
standardize_chromosome_data <- function(chr_counts,
                                      transform = "none",
                                      discrete = TRUE,
                                      standardize = FALSE,
                                      min_value = 1,
                                      verbose = TRUE) {
  # Check input
  if(!is.vector(chr_counts)) {
    stop("chr_counts must be a vector")
  }
  
  # Make a copy of the data
  std_counts <- chr_counts
  
  # Apply minimum value constraint
  if(!is.null(min_value)) {
    below_min <- std_counts < min_value & !is.na(std_counts)
    if(any(below_min)) {
      if(verbose) {
        warning(sprintf("Found %d values below minimum threshold of %d. Setting to NA.",
                       sum(below_min), min_value))
      }
      std_counts[below_min] <- NA
    }
  }
  
  # Apply transformation
  if(transform != "none") {
    if(transform == "log") {
      std_counts <- log(std_counts)
      if(verbose) message("Applied log transformation")
    } else if(transform == "sqrt") {
      std_counts <- sqrt(std_counts)
      if(verbose) message("Applied square root transformation")
    } else if(transform == "rank") {
      std_counts <- rank(std_counts, na.last = "keep", ties.method = "average")
      if(verbose) message("Applied rank transformation")
    } else {
      stop(paste("Unsupported transformation:", transform))
    }
  }
  
  # Standardize if requested
  if(standardize) {
    mean_val <- mean(std_counts, na.rm = TRUE)
    sd_val <- sd(std_counts, na.rm = TRUE)
    
    std_counts <- (std_counts - mean_val) / sd_val
    
    if(verbose) {
      message(sprintf("Standardized data to mean=0, sd=1 (original: mean=%.2f, sd=%.2f)",
                    mean_val, sd_val))
    }
  }
  
  # Round if treating as discrete
  if(discrete && transform != "none") {
    std_counts <- round(std_counts)
    if(verbose) message("Rounded transformed values to integers")
  }
  
  return(std_counts)
}

#' Bin chromosome counts into categories
#' 
#' Creates categorical variables from continuous chromosome count data
#' 
#' @param chr_counts Named vector of chromosome counts
#' @param method Binning method: "equal_width", "equal_freq", "jenks", "custom"
#' @param n_bins Number of bins to create
#' @param custom_breaks Custom break points for binning
#' @param labels Optional labels for bins
#' @param include_counts Whether to include count ranges in labels
#' @param right Whether intervals are closed on the right
#' @param verbose Whether to print messages during processing
#' @return Factor vector of binned chromosome counts
#' @export
bin_chromosome_data <- function(chr_counts,
                              method = "equal_width",
                              n_bins = 5,
                              custom_breaks = NULL,
                              labels = NULL,
                              include_counts = TRUE,
                              right = FALSE,
                              verbose = TRUE) {
  # Check input
  if(!is.vector(chr_counts) || !is.numeric(chr_counts)) {
    stop("chr_counts must be a numeric vector")
  }
  
  # Remove NAs for break calculation
  values <- na.omit(chr_counts)
  
  # Calculate breaks based on method
  if(method == "custom") {
    if(is.null(custom_breaks)) {
      stop("custom_breaks must be provided when method='custom'")
    }
    breaks <- sort(unique(custom_breaks))
  } else if(method == "equal_width") {
    # Equal width bins
    min_val <- min(values)
    max_val <- max(values)
    breaks <- seq(min_val, max_val, length.out = n_bins + 1)
  } else if(method == "equal_freq") {
    # Equal frequency bins
    breaks <- stats::quantile(values, probs = seq(0, 1, length.out = n_bins + 1))
    breaks <- unique(breaks)  # Remove duplicates
  } else if(method == "jenks") {
    # Natural breaks (Jenks)
    if(requireNamespace("classInt", quietly = TRUE)) {
      cints <- classInt::classIntervals(values, n_bins, style = "jenks")
      breaks <- cints$brks
    } else {
      warning("classInt package not available. Using equal width bins instead.")
      min_val <- min(values)
      max_val <- max(values)
      breaks <- seq(min_val, max_val, length.out = n_bins + 1)
    }
  } else {
    stop(paste("Unsupported binning method:", method))
  }
  
  # Ensure breaks span the entire range
  if(min(breaks) > min(values, na.rm = TRUE)) {
    breaks <- c(min(values, na.rm = TRUE), breaks)
  }
  if(max(breaks) < max(values, na.rm = TRUE)) {
    breaks <- c(breaks, max(values, na.rm = TRUE))
  }
  
  # Create labels if not provided
  if(is.null(labels)) {
    if(include_counts) {
      # Create labels with count ranges
      labels <- character(length(breaks) - 1)
      for(i in 1:(length(breaks) - 1)) {
        if(right) {
          labels[i] <- sprintf("(%d,%d]", breaks[i], breaks[i + 1])
        } else {
          labels[i] <- sprintf("[%d,%d)", breaks[i], breaks[i + 1])
        }
      }
    } else {
      # Simple numeric labels
      labels <- as.character(1:(length(breaks) - 1))
    }
  }
  
  # Bin the data
  binned <- cut(chr_counts, breaks = breaks, labels = labels, 
               include.lowest = TRUE, right = right)
  
  # Print summary
  if(verbose) {
    # Count observations in each bin
    bin_counts <- table(binned)
    
    message(sprintf("Binned chromosome counts into %d categories using %s method:",
                  length(labels), method))
    for(i in 1:length(bin_counts)) {
      message(sprintf("  %s: %d observations (%.1f%%)",
                    names(bin_counts)[i], bin_counts[i],
                    100 * bin_counts[i] / sum(bin_counts)))
    }
  }
  
  return(binned)
}

#===============================================================================
# Visualization Functions
#===============================================================================

#' Create a basic histogram of chromosome count distribution
#' 
#' @param chr_counts Named vector of chromosome counts
#' @param binwidth Width of histogram bins
#' @param color Bar color
#' @param fill Bar fill color
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param add_mean Whether to add a vertical line for the mean
#' @param add_median Whether to add a vertical line for the median
#' @return ggplot object
#' @export
plot_chromosome_histogram <- function(chr_counts,
                                    binwidth = NULL,
                                    color = "black",
                                    fill = "steelblue",
                                    title = "Distribution of Chromosome Counts",
                                    xlab = "Chromosome Number",
                                    ylab = "Frequency",
                                    add_mean = TRUE,
                                    add_median = TRUE) {
  # Check for ggplot2
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for plotting")
  }
  
  # Remove NAs
  chr_counts <- chr_counts[!is.na(chr_counts)]
  
  # Create data frame
  plot_data <- data.frame(
    chr_num = chr_counts,
    stringsAsFactors = FALSE
  )
  
  # Calculate binwidth if not provided
  if(is.null(binwidth)) {
    range_val <- max(chr_counts) - min(chr_counts)
    binwidth <- max(1, round(range_val / 30))
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = chr_num)) +
    ggplot2::geom_histogram(binwidth = binwidth, color = color, fill = fill) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
    ggplot2::theme_minimal()
  
  # Add mean line if requested
  if(add_mean) {
    mean_val <- mean(chr_counts)
    p <- p + ggplot2::geom_vline(xintercept = mean_val, 
                               color = "red", linetype = "dashed", size = 1) +
      ggplot2::annotate("text", x = mean_val, y = 0, 
                      label = sprintf("Mean: %.1f", mean_val),
                      hjust = -0.1, vjust = -0.5, color = "red")
  }
  
  # Add median line if requested
  if(add_median) {
    median_val <- median(chr_counts)
    p <- p + ggplot2::geom_vline(xintercept = median_val, 
                               color = "blue", linetype = "dashed", size = 1) +
      ggplot2::annotate("text", x = median_val, y = 0, 
                      label = sprintf("Median: %d", median_val),
                      hjust = 1.1, vjust = -0.5, color = "blue")
  }
  
  return(p)
}

#' Map chromosome counts onto a phylogeny
#' 
#' Visualizes chromosome count data on a phylogenetic tree
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param method Visualization method: "color", "size", "labels", or "combined"
#' @param palette Color palette to use
#' @param reverse_palette Whether to reverse the color palette
#' @param show_legend Whether to show the legend
#' @param title Plot title
#' @param layout Tree layout: "rectangular", "circular", "fan", etc.
#' @param tip_labels Whether to show tip labels
#' @param label_offset Offset for tip labels
#' @return ggplot object
#' @export
plot_chromosome_phylogeny <- function(tree,
                                    chr_counts,
                                    method = "color",
                                    palette = "viridis",
                                    reverse_palette = FALSE,
                                    show_legend = TRUE,
                                    title = "Chromosome Numbers Across the Phylogeny",
                                    layout = "rectangular",
                                    tip_labels = TRUE,
                                    label_offset = 0.1) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggtree", quietly = TRUE)) {
    stop("ggplot2 and ggtree packages are required for phylogeny plotting")
  }
  
  # Ensure chr_counts match tree tips
  matched_counts <- match_chromosome_data(tree, chr_counts, impute = FALSE, verbose = FALSE)
  
  # Create base tree plot
  p <- ggtree::ggtree(tree, layout = layout)
  
  # Add chromosome data based on method
  if(method == "color" || method == "combined") {
    # Use color to represent chromosome counts
    p <- p %<+% data.frame(label = names(matched_counts), chr_num = matched_counts) +
      ggtree::aes(color = chr_num)
      
    # Add color scale
    if(palette == "viridis" && requireNamespace("viridis", quietly = TRUE)) {
      p <- p + viridis::scale_color_viridis(
        name = "Chromosome Number",
        na.value = "gray80",
        direction = ifelse(reverse_palette, -1, 1),
        option = "D"
      )
    } else {
      p <- p + ggplot2::scale_color_gradient(
        name = "Chromosome Number",
        low = ifelse(reverse_palette, "red", "blue"),
        high = ifelse(reverse_palette, "blue", "red"),
        na.value = "gray80"
      )
    }
  }
  
  if(method == "size" || method == "combined") {
    # Use point size to represent chromosome counts
    if(method == "size") {
      p <- p %<+% data.frame(label = names(matched_counts), chr_num = matched_counts)
    }
    
    # Add tip points
    p <- p + ggtree::geom_tippoint(aes(size = chr_num), alpha = 0.7)
    
    # Add size scale
    p <- p + ggplot2::scale_size_continuous(
      name = "Chromosome Number",
      range = c(1, 5)
    )
  }
  
  if(method == "labels") {
    # Add chromosome count labels directly
    p <- p + ggtree::geom_tiplab(aes(label = paste0(" ", tree$tip.label, " (", 
                                                 matched_counts, ")")),
                            align = TRUE, offset = label_offset)
  } else if(tip_labels) {
    # Add regular tip labels
    p <- p + ggtree::geom_tiplab(align = TRUE, offset = label_offset)
  }
  
  # Add title
  p <- p + ggplot2::labs(title = title)
  
  # Adjust legend
  if(!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  return(p)
}

#' Create a heatmap of chromosome number changes across the phylogeny
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts
#' @param ancestral_states Optional ancestral state reconstruction results
#' @param method Method for ancestral state reconstruction if not provided
#' @param palette Color palette to use
#' @param cluster_tips Whether to cluster species by chromosome number
#' @param show_counts Whether to show count values in cells
#' @param title Plot title
#' @return ggplot object
#' @export
plot_chromosome_heatmap <- function(tree,
                                  chr_counts,
                                  ancestral_states = NULL,
                                  method = "ML",
                                  palette = "viridis",
                                  cluster_tips = FALSE,
                                  show_counts = TRUE,
                                  title = "Chromosome Number Evolution") {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE) || 
     !requireNamespace("ggtree", quietly = TRUE) || 
     !requireNamespace("ape", quietly = TRUE)) {
    stop("ggplot2, ggtree, and ape packages are required for heatmap plotting")
  }
  
  # Ensure chr_counts match tree tips
  matched_counts <- match_chromosome_data(tree, chr_counts, impute = FALSE, verbose = FALSE)
  
  # Perform ancestral state reconstruction if not provided
  if(is.null(ancestral_states)) {
    if(method == "ML") {
      ace_result <- ape::ace(matched_counts, tree, type = "continuous", method = "ML")
      ancestral_states <- ace_result$ace
    } else if(method == "parsimony") {
      ace_result <- ape::ace(matched_counts, tree, type = "continuous", method = "REML")
      ancestral_states <- ace_result$ace
    } else {
      stop(paste("Unsupported ancestral state reconstruction method:", method))
    }
  }
  
  # Combine tip and internal node states
  all_states <- c(matched_counts, ancestral_states)
  
  # Create data frame for heatmap
  node_ids <- c(1:length(tree$tip.label), 
                (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode))
  
  heatmap_data <- data.frame(
    node = node_ids,
    chr_num = all_states,
    is_tip = c(rep(TRUE, length(tree$tip.label)), 
              rep(FALSE, tree$Nnode)),
    label = c(tree$tip.label, 
             (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode)),
    stringsAsFactors = FALSE
  )
  
  # Create tree with node positions
  tree_plot <- ggtree::ggtree(tree, layout = ifelse(cluster_tips, "circular", "rectangular"))
  
  # Extract node positions
  node_pos <- tree_plot$data[, c("node", "x", "y")]
  
  # Merge with heatmap data
  heatmap_data <- merge(heatmap_data, node_pos, by = "node")
  
  # Create heatmap plot
  p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = x, y = y, fill = chr_num)) +
    ggplot2::geom_tile(width = 0.8, height = 0.8) +
    ggplot2::labs(title = title, fill = "Chromosome Number") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  
  # Add color scale
  if(palette == "viridis" && requireNamespace("viridis", quietly = TRUE)) {
    p <- p + viridis::scale_fill_viridis(
      name = "Chromosome Number",
      na.value = "gray80",
      option = "D"
    )
  } else {
    p <- p + ggplot2::scale_fill_gradient(
      name = "Chromosome Number",
      low = "blue",
      high = "red",
      na.value = "gray80"
    )
  }
  
  # Add count values if requested
  if(show_counts) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = round(chr_num)),
      color = "white",
      size = 3
    )
  }
  
  # Add tree lines
  tree_data <- tree_plot$data
  edge_data <- tree_plot$data[tree_plot$data$.edge > 0, ]
  
  p <- p + ggplot2::geom_segment(
    data = edge_data,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
    color = "black",
    alpha = 0.5,
    size = 0.5
  )
  
  # Add tip labels
  p <- p + ggplot2::geom_text(
    data = heatmap_data[heatmap_data$is_tip, ],
    ggplot2::aes(label = label),
    hjust = -0.1,
    size = 3
  )
  
  return(p)
}

#===============================================================================
# Data Export Functions
#===============================================================================

#' Export processed chromosome data to various formats
#' 
#' @param chr_counts Named vector of chromosome counts
#' @param file Output file path
#' @param format Output format: "csv", "nexus", "phylip", or "json"
#' @param include_metadata Whether to include metadata in output
#' @param metadata Additional metadata to include
#' @param overwrite Whether to overwrite existing file
#' @param sep Field separator for tabular output
#' @return TRUE if export successful
#' @export
export_chromosome_data <- function(chr_counts,
                                 file,
                                 format = "csv",
                                 include_metadata = TRUE,
                                 metadata = NULL,
                                 overwrite = FALSE,
                                 sep = ",") {
  # Check if file exists
  if(file.exists(file) && !overwrite) {
    stop(paste("File already exists:", file, "(use overwrite=TRUE to replace)"))
  }
  
  # Check input
  if(!is.vector(chr_counts) || is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector")
  }
  
  # Create data frame
  df <- data.frame(
    Species = names(chr_counts),
    ChromosomeCount = chr_counts,
    stringsAsFactors = FALSE
  )
  
  # Add metadata if requested
  if(include_metadata) {
    # Add basic metadata
    df$TimeStamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    
    # Add custom metadata if provided
    if(!is.null(metadata) && is.list(metadata)) {
      for(field in names(metadata)) {
        df[[field]] <- metadata[[field]]
      }
    }
  }
  
  # Export based on format
  if(format == "csv") {
    utils::write.csv(df, file = file, row.names = FALSE)
  } else if(format == "tsv" || format == "txt") {
    utils::write.table(df, file = file, sep = sep, 
                      row.names = FALSE, quote = TRUE)
  } else if(format == "nexus") {
    if(requireNamespace("ape", quietly = TRUE)) {
      # Create NEXUS format
      nexus_data <- list()
      
      for(i in 1:nrow(df)) {
        species <- df$Species[i]
        count <- df$ChromosomeCount[i]
        
        # Handle missing values
        if(is.na(count)) {
          count_str <- "?"
        } else {
          count_str <- as.character(count)
        }
        
        nexus_data[[species]] <- count_str
      }
      
      # Write NEXUS file
      ape::write.nexus.data(nexus_data, file = file)
    } else {
      stop("ape package is required for NEXUS export")
    }
  } else if(format == "phylip") {
    # Create PHYLIP format
    phylip_lines <- c(
      paste(nrow(df), 1),  # Number of taxa and characters
      vapply(1:nrow(df), function(i) {
        # Format species name (padded to 10 characters)
        species <- df$Species[i]
        if(nchar(species) > 10) {
          species <- substr(species, 1, 10)
        } else {
          species <- sprintf("%-10s", species)
        }
        
        # Handle missing values
        count <- df$ChromosomeCount[i]
        if(is.na(count)) {
          count_str <- "?"
        } else {
          count_str <- as.character(count)
        }
        
        paste0(species, count_str)
      }, character(1))
    )
    
    # Write PHYLIP file
    writeLines(phylip_lines, file)
  } else if(format == "json") {
    if(requireNamespace("jsonlite", quietly = TRUE)) {
      # Create JSON format
      json_data <- df
      
      # Write JSON file
      jsonlite::write_json(json_data, file, pretty = TRUE)
    } else {
      stop("jsonlite package is required for JSON export")
    }
  } else {
    stop(paste("Unsupported export format:", format))
  }
  
  return(TRUE)
}

#===============================================================================
# Utility Functions
#===============================================================================

#' Summarize chromosome count data
#' 
#' Creates a comprehensive summary of chromosome count data
#' 
#' @param chr_counts Named vector of chromosome counts
#' @param tree Optional phylogenetic tree for taxonomic summary
#' @param by_clade Whether to summarize data by clade
#' @param quantiles Quantiles to include in summary
#' @param verbose Whether to print summary to console
#' @return List with summary statistics
#' @export
summarize_chromosome_data <- function(chr_counts,
                                    tree = NULL,
                                    by_clade = FALSE,
                                    quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),
                                    verbose = TRUE) {
  # Check input
  if(!is.vector(chr_counts)) {
    stop("chr_counts must be a vector")
  }
  
  # Initialize summary
  summary <- list(
    n_species = length(chr_counts),
    n_missing = sum(is.na(chr_counts)),
    percent_complete = 100 * (1 - sum(is.na(chr_counts)) / length(chr_counts)),
    min = min(chr_counts, na.rm = TRUE),
    max = max(chr_counts, na.rm = TRUE),
    range = diff(range(chr_counts, na.rm = TRUE)),
    mean = mean(chr_counts, na.rm = TRUE),
    median = median(chr_counts, na.rm = TRUE),
    mode = as.numeric(names(which.max(table(chr_counts)))),
    sd = sd(chr_counts, na.rm = TRUE),
    cv = 100 * sd(chr_counts, na.rm = TRUE) / mean(chr_counts, na.rm = TRUE),
    skewness = skewness(chr_counts),
    kurtosis = kurtosis(chr_counts)
  )
  
  # Add quantiles
  summary$quantiles <- stats::quantile(chr_counts, probs = quantiles, na.rm = TRUE)
  
  # Add frequency table
  freq_table <- table(chr_counts)
  summary$frequencies <- as.data.frame(freq_table)
  colnames(summary$frequencies) <- c("ChromosomeCount", "Frequency")
  summary$frequencies$Percentage <- 100 * summary$frequencies$Frequency / sum(summary$frequencies$Frequency)
  
  # Summarize by clade if requested and tree available
  if(by_clade && !is.null(tree)) {
    if(requireNamespace("ape", quietly = TRUE)) {
      # Identify major clades
      clades <- identify_major_clades(tree, min_size = 5)
      
      # Calculate summary for each clade
      clade_summaries <- list()
      
      for(clade_name in names(clades)) {
        # Extract species in this clade
        clade_species <- clades[[clade_name]]
        
        # Get chromosome counts for this clade
        clade_counts <- chr_counts[names(chr_counts) %in% clade_species]
        
        # Calculate summary statistics
        clade_summary <- list(
          n_species = length(clade_counts),
          n_missing = sum(is.na(clade_counts)),
          percent_complete = 100 * (1 - sum(is.na(clade_counts)) / length(clade_counts)),
          min = min(clade_counts, na.rm = TRUE),
          max = max(clade_counts, na.rm = TRUE),
          range = diff(range(clade_counts, na.rm = TRUE)),
          mean = mean(clade_counts, na.rm = TRUE),
          median = median(clade_counts, na.rm = TRUE),
          sd = sd(clade_counts, na.rm = TRUE),
          cv = 100 * sd(clade_counts, na.rm = TRUE) / mean(clade_counts, na.rm = TRUE)
        )
        
        clade_summaries[[clade_name]] <- clade_summary
      }
      
      summary$clade_summaries <- clade_summaries
    } else {
      warning("ape package is required for clade-based summary")
    }
  }
  
  # Print summary if requested
  if(verbose) {
    cat("Chromosome Count Data Summary\n")
    cat("============================\n\n")
    
    cat(sprintf("Number of species: %d (%d with data, %.1f%% complete)\n",
              summary$n_species, 
              summary$n_species - summary$n_missing,
              summary$percent_complete))
    
    cat(sprintf("\nRange: %d to %d (span: %d)\n",
              summary$min, summary$max, summary$range))
    
    cat(sprintf("Central tendency: mean = %.2f, median = %d, mode = %d\n",
              summary$mean, summary$median, summary$mode))
    
    cat(sprintf("Dispersion: SD = %.2f, CV = %.1f%%\n",
              summary$sd, summary$cv))
    
    cat(sprintf("Distribution shape: skewness = %.3f, kurtosis = %.3f\n",
              summary$skewness, summary$kurtosis))
    
    cat("\nQuantiles:\n")
    for(i in 1:length(summary$quantiles)) {
      cat(sprintf("  %s: %d\n", 
                names(summary$quantiles)[i], 
                summary$quantiles[i]))
    }
    
    cat("\nMost common chromosome numbers:\n")
    top_counts <- head(summary$frequencies[order(-summary$frequencies$Frequency), ], 5)
    for(i in 1:nrow(top_counts)) {
      cat(sprintf("  %d: %d species (%.1f%%)\n", 
                top_counts$ChromosomeCount[i],
                top_counts$Frequency[i],
                top_counts$Percentage[i]))
    }
    
    if(by_clade && !is.null(summary$clade_summaries)) {
      cat("\nSummary by clade:\n")
      for(clade_name in names(summary$clade_summaries)) {
        clade <- summary$clade_summaries[[clade_name]]
        cat(sprintf("\n  %s (%d species, %.1f%% complete):\n", 
                  clade_name, clade$n_species, clade$percent_complete))
        cat(sprintf("    Range: %d to %d (span: %d)\n",
                  clade$min, clade$max, clade$range))
        cat(sprintf("    Mean = %.2f, Median = %d, SD = %.2f, CV = %.1f%%\n",
                  clade$mean, clade$median, clade$sd, clade$cv))
      }
    }
  }
  
  return(summary)
}

#' Calculate skewness of a distribution
#' 
#' @param x Numeric vector
#' @return Skewness value
#' @keywords internal
skewness <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if(n < 3) return(NA)
  
  m3 <- sum((x - mean(x))^3) / n
  s3 <- (sum((x - mean(x))^2) / n)^(3/2)
  
  if(s3 == 0) return(0)
  return(m3 / s3)
}

#' Calculate kurtosis of a distribution
#' 
#' @param x Numeric vector
#' @return Kurtosis value
#' @keywords internal
kurtosis <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  if(n < 4) return(NA)
  
  m4 <- sum((x - mean(x))^4) / n
  s4 <- (sum((x - mean(x))^2) / n)^2
  
  if(s4 == 0) return(0)
  return(m4 / s4 - 3)  # Excess kurtosis (normal = 0)
}
