#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Data Integration Module
# Author: Bioinformatics Team
# Date: 2025-04-06
# Description: Integrates and processes phylogenetic trees, chromosome counts and
#              bidirectional mapping data, ensuring data consistency
#===============================================================================

#' Load and integrate multiple data sources
#' 
#' Integrate phylogenetic tree, chromosome counts and bidirectional mapping data,
#' ensuring consistent species are used in analysis
#' 
#' @param tree_file Phylogenetic tree file path (Newick format)
#' @param counts_file Chromosome counts file path (CSV format)
#' @param maps_dir Bidirectional mapping data directory (optional)
#' @param process_common_only Whether to only process species common to all datasets
#' @param min_common_species Minimum number of common species
#' @return List of integrated data
#' @export
integrate_data_sources <- function(tree_file, 
                                 counts_file, 
                                 maps_dir = NULL,
                                 process_common_only = TRUE,
                                 min_common_species = 3) {
  # Load phylogenetic tree
  tree <- load_phylogenetic_tree(tree_file)
  
  # Load chromosome counts
  chr_counts_data <- load_chromosome_counts(counts_file)
  chr_counts <- chr_counts_data$counts
  
  # Extract species names
  tree_species <- tree$tip.label
  count_species <- names(chr_counts)
  
  # Load mapping data if provided
  mapping_species <- character(0)
  mapping_data <- NULL
  if (!is.null(maps_dir) && dir.exists(maps_dir)) {
    map_data <- load_bidirectional_maps(maps_dir)
    mapping_species <- map_data$species
    mapping_data <- map_data$combined
  }
  
  # Find common species
  if (length(mapping_species) > 0) {
    common_species <- Reduce(intersect, list(tree_species, count_species, mapping_species))
  } else {
    common_species <- intersect(tree_species, count_species)
  }
  
  # Check if enough common species
  if (length(common_species) < min_common_species) {
    stop(paste0("Number of common species (", length(common_species), 
              ") less than minimum required (", min_common_species, ")"))
  }
  
  # If only processing common species, trim data
  if (process_common_only) {
    # Prune tree
    tree <- prune_phylogenetic_tree(tree, common_species)
    
    # Filter chromosome counts
    chr_counts <- filter_chromosome_counts(chr_counts, common_species)
    
    # Filter mapping data (if present)
    if (!is.null(mapping_data)) {
      mapping_data <- filter_bidirectional_maps(mapping_data, common_species)
    }
    
    message(paste("Restricted datasets to", length(common_species), "common species"))
  } else {
    message(paste("Found", length(common_species), "common species, but keeping all data"))
    # Logic for handling incomplete data could be added here
  }
  
  # Return integrated data
  return(list(
    tree = tree,
    chr_counts = chr_counts,
    mapping_data = mapping_data,
    common_species = common_species,
    all_tree_species = tree_species,
    all_count_species = count_species,
    all_mapping_species = mapping_species
  ))
}

#' Load phylogenetic tree
#' 
#' @param tree_file Tree file path
#' @return Phylogenetic tree object
#' @keywords internal
load_phylogenetic_tree <- function(tree_file) {
  message("Loading phylogenetic tree...")
  
  if (!file.exists(tree_file)) {
    stop("Phylogenetic tree file doesn't exist: ", tree_file)
  }
  
  tree <- try(ape::read.tree(tree_file), silent = TRUE)
  
  if (inherits(tree, "try-error")) {
    stop("Error reading phylogenetic tree file: ", tree_file)
  }
  
  # Basic validation
  if (is.null(tree) || length(tree$tip.label) < 3) {
    stop("Invalid tree file or insufficient species count (<3)")
  }
  
  message(paste("Loaded tree with", length(tree$tip.label), "species and", 
              tree$Nnode, "internal nodes"))
  
  return(tree)
}

#' Load chromosome count data
#' 
#' @param counts_file Chromosome counts file path
#' @return List of chromosome count data
#' @keywords internal
load_chromosome_counts <- function(counts_file) {
  message("Loading chromosome counts...")
  
  if (!file.exists(counts_file)) {
    stop("Chromosome counts file doesn't exist: ", counts_file)
  }
  
  counts_data <- try(read.csv(counts_file, stringsAsFactors = FALSE), silent = TRUE)
  
  if (inherits(counts_data, "try-error")) {
    stop("Error reading chromosome counts file")
  }
  
  # Check required columns
  col_species <- "species"
  col_count <- "count"
  
  # Try to auto-detect column names
  if (!col_species %in% colnames(counts_data)) {
    # Try to find possible species column
    possible_species_cols <- c("species", "Species", "taxon", "Taxon", "organism", "Organism")
    found_cols <- intersect(possible_species_cols, colnames(counts_data))
    
    if (length(found_cols) > 0) {
      col_species <- found_cols[1]
      message(paste("Using", col_species, "as species column"))
    } else {
      stop("Could not find species column. Ensure file contains a 'species' column")
    }
  }
  
  if (!col_count %in% colnames(counts_data)) {
    # Try to find possible count column
    possible_count_cols <- c("count", "Count", "chromosome_count", "ChromosomeCount",
                           "n", "2n", "haploid", "diploid", "chr_count")
    found_cols <- intersect(possible_count_cols, colnames(counts_data))
    
    if (length(found_cols) > 0) {
      col_count <- found_cols[1]
      message(paste("Using", col_count, "as chromosome count column"))
    } else {
      stop("Could not find chromosome count column. Ensure file contains a 'count' column")
    }
  }
  
  # Convert to named vector for ease of use
  chr_counts <- counts_data[[col_count]]
  names(chr_counts) <- counts_data[[col_species]]
  
  # Check and clean data
  chr_counts <- as.numeric(chr_counts)
  if (any(is.na(chr_counts))) {
    warning(paste("Found", sum(is.na(chr_counts)), "missing or non-numeric chromosome counts"))
  }
  
  message(paste("Loaded chromosome counts for", length(chr_counts), "species"))
  
  return(list(
    data_frame = counts_data,
    counts = chr_counts
  ))
}

#' Load bidirectional mapping data
#' 
#' @param maps_dir Mapping data directory
#' @return List of mapping data
#' @keywords internal
load_bidirectional_maps <- function(maps_dir) {
  message("Loading bidirectional mapping data...")
  
  if (!dir.exists(maps_dir)) {
    stop("Bidirectional mapping directory doesn't exist: ", maps_dir)
  }
  
  # Find all mapping files
  map_files <- list.files(maps_dir, pattern = "*_bidirectional\\.tsv$", full.names = TRUE)
  
  if (length(map_files) == 0) {
    # If no specific pattern found, try any TSV file
    map_files <- list.files(maps_dir, pattern = "\\.tsv$", full.names = TRUE)
    
    if (length(map_files) == 0) {
      stop("No mapping files found in directory: ", maps_dir)
    }
  }
  
  message(paste("Found", length(map_files), "bidirectional mapping files"))
  
  # Read all files
  has_data_table <- requireNamespace("data.table", quietly = TRUE)
  
  if (has_data_table) {
    # If data.table is installed, use more efficient fread
    all_maps <- lapply(map_files, function(file) {
      tryCatch({
        map_data <- data.table::fread(file)
        
        # Check required columns
        req_cols <- c("species_A", "chromosome_A", "species_B", "chromosome_B")
        if (!all(req_cols %in% names(map_data))) {
          warning("File missing required columns, skipping: ", file)
          return(NULL)
        }
        return(map_data)
      }, error = function(e) {
        warning("Error reading file: ", file, " - ", e$message)
        return(NULL)
      })
    })
    
    # Remove NULL items (failed files)
    all_maps <- all_maps[!sapply(all_maps, is.null)]
    
    if (length(all_maps) == 0) {
      stop("Could not load any valid mapping data")
    }
    
    # Combine all mapping data
    combined_maps <- data.table::rbindlist(all_maps, fill = TRUE)
  } else {
    # If no data.table, use base R functions
    all_maps <- lapply(map_files, function(file) {
      tryCatch({
        map_data <- read.delim(file, sep = "\t", stringsAsFactors = FALSE)
        
        # Check required columns
        req_cols <- c("species_A", "chromosome_A", "species_B", "chromosome_B")
        if (!all(req_cols %in% names(map_data))) {
          warning("File missing required columns, skipping: ", file)
          return(NULL)
        }
        return(map_data)
      }, error = function(e) {
        warning("Error reading file: ", file, " - ", e$message)
        return(NULL)
      })
    })
    
    # Remove NULL items (failed files)
    all_maps <- all_maps[!sapply(all_maps, is.null)]
    
    if (length(all_maps) == 0) {
      stop("Could not load any valid mapping data")
    }
    
    # Combine all mapping data
    combined_maps <- do.call(rbind, all_maps)
  }
  
  message(paste("Combined mapping data contains", nrow(combined_maps), "rows"))
  
  # Extract all species in mapping data
  mapping_species <- unique(c(combined_maps$species_A, combined_maps$species_B))
  message(paste("Mapping data contains", length(mapping_species), "species"))
  
  return(list(
    combined = combined_maps,
    species = mapping_species
  ))
}

#' Prune phylogenetic tree
#' 
#' @param tree Phylogenetic tree object
#' @param common_species Species to keep
#' @return Pruned tree
#' @keywords internal
prune_phylogenetic_tree <- function(tree, common_species) {
  message("Pruning phylogenetic tree to only include common species...")
  
  pruned_tree <- try(ape::keep.tip(tree, common_species), silent = TRUE)
  
  if (inherits(pruned_tree, "try-error")) {
    stop("Error pruning phylogenetic tree")
  }
  
  # Verify the pruning was correct
  if (length(pruned_tree$tip.label) != length(common_species)) {
    stop("Error: Pruned tree does not contain the exact number of common species")
  }
  
  # Check all leaves in pruned tree are in common_species
  if (!all(pruned_tree$tip.label %in% common_species)) {
    stop("Error: Pruned tree contains species not in common species list")
  }
  
  message(paste("Pruned tree contains", length(pruned_tree$tip.label), "species and", 
              pruned_tree$Nnode, "internal nodes"))
  
  return(pruned_tree)
}

#' Filter chromosome counts
#' 
#' @param chr_counts Chromosome counts named vector
#' @param common_species Species to keep
#' @return Filtered chromosome counts
#' @keywords internal
filter_chromosome_counts <- function(chr_counts, common_species) {
  message("Filtering chromosome counts to only common species...")
  
  filtered_counts <- chr_counts[common_species]
  
  # Verify filtering was correct
  if (length(filtered_counts) != length(common_species)) {
    stop("Error: Filtered counts does not contain the exact number of common species")
  }
  
  # Check all names in filtered counts are in common_species
  if (!all(names(filtered_counts) %in% common_species)) {
    stop("Error: Filtered counts contains species not in common species list")
  }
  
  message(paste("Filtered to", length(filtered_counts), "species' chromosome counts"))
  
  return(filtered_counts)
}

#' Filter bidirectional mapping data
#' 
#' @param mapping_data Bidirectional mapping data frame
#' @param common_species Species to keep
#' @return Filtered mapping data
#' @keywords internal
filter_bidirectional_maps <- function(mapping_data, common_species) {
  message("Filtering bidirectional mapping data to only common species...")
  
  if (requireNamespace("data.table", quietly = TRUE) && inherits(mapping_data, "data.table")) {
    filtered_maps <- mapping_data[species_A %in% common_species & 
                                species_B %in% common_species, ]
  } else {
    filtered_maps <- mapping_data[mapping_data$species_A %in% common_species & 
                                mapping_data$species_B %in% common_species, ]
  }
  
  if (nrow(filtered_maps) == 0) {
    warning("Filtered mapping data is empty, which may indicate no mapping relationships between common species")
  }
  
  message(paste("Filtered mapping data contains", nrow(filtered_maps), "rows"))
  
  return(filtered_maps)
}
