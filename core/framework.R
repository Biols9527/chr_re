#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Core Management Module
# Author: Bioinformatics Team
# Date: 2025-03-22
# Description: Provides modular architecture, unified data flow, config management
#              and coordinates between analysis modules
#===============================================================================

#' Load required packages
#' @importFrom magrittr %>%
#' @importFrom R6 R6Class
suppressPackageStartupMessages({
  library(R6)
  library(yaml)
  library(magrittr)
  library(ape)
  library(data.table)
  library(parallel)
})

#===============================================================================
# Configuration Management System
#===============================================================================

#' Load analysis configuration
#' 
#' Load configuration from YAML file, merging default config and user config
#' 
#' @param config_file Configuration file path, NULL to use default config
#' @param override Override configuration parameters list
#' @return Configuration list
#' @export
load_config <- function(config_file = NULL, override = NULL) {
  # Default configuration
  default_config <- list(
    # Phylogenetic parameters
    phylo = list(
      rooting_method = "midpoint",     # Tree rooting method: midpoint, outgroup
      branch_length_treatment = "keep", # Branch length treatment: keep, estimate, unit
      ultrametric = FALSE              # Whether to force ultrametric tree
    ),
    
    # Network analysis parameters
    network = list(
      weight_type = "combined",        # Edge weight type: combined, count, quality
      phylo_weighting = TRUE,          # Whether to use phylogenetic weighting
      phylo_weight_method = "distance", # Phylogenetic weighting method: distance, topology
      community_methods = c("louvain", "fast_greedy", "walktrap"), # Community detection algorithms
      use_ensemble = TRUE,             # Whether to use ensemble community detection
      resolution = 1.0                 # Community detection resolution
    ),
    
    # Event inference parameters
    events = list(
      fusion_threshold = 0.67,         # Chromosome fusion detection threshold
      fission_threshold = 1.5,         # Chromosome fission detection threshold
      wgd_threshold = 2.0,             # Whole genome duplication detection threshold
      allow_wgd = TRUE,                # Whether to allow whole genome duplication
      confidence_levels = list(        # Confidence thresholds
        low = 0.3, 
        medium = 0.6, 
        high = 0.8
      )
    ),
    
    # Ancestral reconstruction parameters
    ancestral = list(
      methods = c("parsimony", "bayesian", "ml"), # Methods to use
      compare_methods = TRUE,          # Whether to compare different methods
      generate_consensus = TRUE,       # Whether to generate consensus results
      bottom_up_weight = 0.6,          # Weight for bottom-up methods
      top_down_weight = 0.4,           # Weight for top-down methods
      max_chromosomes = 100            # Maximum allowed chromosome number
    ),
    
    # Computation parameters
    compute = list(
      use_parallel = TRUE,             # Whether to use parallel computation
      cores = 0,                       # Number of cores (0=auto)
      memory_limit = 8000,             # Memory limit (MB)
      temp_dir = tempdir()             # Temporary file directory
    ),
    
    # Output parameters
    output = list(
      base_dir = "results",            # Base results directory
      save_intermediates = TRUE,       # Whether to save intermediate results
      plots = TRUE,                    # Whether to generate plots
      plot_format = "pdf",             # Plot format: pdf, png
      report_type = "full",            # Report type: full, summary, minimal
      compress_results = FALSE         # Whether to compress results
    ),
    
    # Advanced parameters
    advanced = list(
      seed = 42,                       # Random seed
      debug = FALSE,                   # Debug mode
      verbose = TRUE                   # Verbose output
    )
  )
  
  # Load user configuration from YAML file
  user_config <- list()
  if(!is.null(config_file) && file.exists(config_file)) {
    tryCatch({
      user_config <- yaml::read_yaml(config_file)
      message("Configuration loaded from file: ", config_file)
    }, error = function(e) {
      warning("Could not load config file, using default config: ", e$message)
    })
  }
  
  # Merge configurations
  config <- merge_configs(default_config, user_config)
  
  # Apply override parameters
  if(!is.null(override) && is.list(override)) {
    config <- merge_configs(config, override)
  }
  
  # Validate configuration
  config <- validate_config(config)
  
  return(config)
}

#' Recursively merge configurations
#' 
#' @param default Default configuration list
#' @param user User configuration list
#' @return Merged configuration list
merge_configs <- function(default, user) {
  if(!is.list(user)) return(default)
  
  result <- default
  
  for(name in names(user)) {
    if(name %in% names(default) && is.list(default[[name]]) && is.list(user[[name]])) {
      result[[name]] <- merge_configs(default[[name]], user[[name]])
    } else {
      result[[name]] <- user[[name]]
    }
  }
  
  return(result)
}

#' Validate configuration validity
#' 
#' @param config Configuration list
#' @return Validated configuration list
validate_config <- function(config) {
  # Set random seed
  set.seed(config$advanced$seed)
  
  # If using parallel computation, set number of cores
  if(config$compute$use_parallel) {
    if(config$compute$cores <= 0) {
      config$compute$cores <- max(1, parallel::detectCores() - 1)
    }
    # Ensure reasonable core count
    config$compute$cores <- min(config$compute$cores, parallel::detectCores())
  }
  
  # Validate phylogenetic parameters
  if(!config$phylo$rooting_method %in% c("midpoint", "outgroup")) {
    warning("Invalid tree rooting method, using default 'midpoint'")
    config$phylo$rooting_method <- "midpoint"
  }
  
  # Validate event threshold parameters
  if(config$events$fusion_threshold <= 0 || config$events$fusion_threshold >= 1) {
    warning("Fusion threshold should be between 0-1, setting to default 0.67")
    config$events$fusion_threshold <- 0.67
  }
  
  if(config$events$fission_threshold <= 1) {
    warning("Fission threshold should be greater than 1, setting to default 1.5")
    config$events$fission_threshold <- 1.5
  }
  
  return(config)
}

#===============================================================================
# Data Flow Management Class
#===============================================================================

#' Data Flow Management Class
#' 
#' Provides state tracking and data exchange in the data processing flow
#' @export
DataFlow <- R6::R6Class(
  "DataFlow",
  public = list(
    #' @field data List for storing data
    data = list(),
    
    #' @field metadata Data metadata
    metadata = list(),
    
    #' @field history Operation history
    history = list(),
    
    #' @description Initialize data flow
    #' @return Data flow object
    initialize = function() {
      self$data <- list()
      self$metadata <- list(
        creation_time = Sys.time(),
        last_modified = Sys.time()
      )
      self$history <- list()
      self$log_event("initialize", "Created data flow object")
      return(invisible(self))
    },
    
    #' @description Add data
    #' @param name Data name
    #' @param value Data value
    #' @param stage Processing stage
    #' @param description Data description
    #' @return Data flow object
    add = function(name, value, stage = "unknown", description = "") {
      self$data[[name]] <- value
      
      # Record metadata
      self$metadata[[name]] <- list(
        added_time = Sys.time(),
        stage = stage,
        description = description,
        class = class(value)[1],
        size = object.size(value)
      )
      
      # Update last modified time
      self$metadata$last_modified <- Sys.time()
      
      # Record history
      self$log_event("add", sprintf("Added data '%s' (stage: %s)", name, stage))
      
      return(invisible(self))
    },
    
    #' @description Get data
    #' @param name Data name
    #' @return Data value
    get = function(name) {
      if(name %in% names(self$data)) {
        self$log_event("get", sprintf("Retrieved data '%s'", name))
        return(self$data[[name]])
      } else {
        warning(sprintf("Data item does not exist: '%s'", name))
        return(NULL)
      }
    },
    
    #' @description Update data
    #' @param name Data name
    #' @param value Data value
    #' @param stage Processing stage
    #' @param description Data description
    #' @return Data flow object
    update = function(name, value, stage = NULL, description = NULL) {
      if(name %in% names(self$data)) {
        old_value <- self$data[[name]]
        self$data[[name]] <- value
        
        # Update metadata
        md <- self$metadata[[name]]
        md$updated_time <- Sys.time()
        if(!is.null(stage)) md$stage <- stage
        if(!is.null(description)) md$description <- description
        md$class <- class(value)[1]
        md$size <- object.size(value)
        self$metadata[[name]] <- md
        
        # Update last modified time
        self$metadata$last_modified <- Sys.time()
        
        # Record history
        self$log_event("update", sprintf("Updated data '%s' (stage: %s)", 
                                        name, ifelse(is.null(stage), md$stage, stage)))
      } else {
        if(is.null(description)) description <- ""
        self$add(name, value, stage, description)
      }
      
      return(invisible(self))
    },
    
    #' @description Remove data
    #' @param name Data name
    #' @return Data flow object
    remove = function(name) {
      if(name %in% names(self$data)) {
        self$data[[name]] <- NULL
        self$metadata[[name]] <- NULL
        
        # Update last modified time
        self$metadata$last_modified <- Sys.time()
        
        # Record history
        self$log_event("remove", sprintf("Removed data '%s'", name))
      }
      
      return(invisible(self))
    },
    
    #' @description Save data flow to file
    #' @param file File path
    #' @return Data flow object
    save = function(file) {
      # Update metadata before saving
      self$metadata$saved_time <- Sys.time()
      self$metadata$saved_file <- file
      
      # Save data and metadata
      to_save <- list(
        data = self$data,
        metadata = self$metadata,
        history = self$history
      )
      
      saveRDS(to_save, file)
      self$log_event("save", sprintf("Saved data flow to file '%s'", file))
      
      return(invisible(self))
    },
    
    #' @description Load data flow from file
    #' @param file File path
    #' @return Data flow object
    load = function(file) {
      if(file.exists(file)) {
        loaded <- readRDS(file)
        
        # Load data and metadata
        self$data <- loaded$data
        self$metadata <- loaded$metadata
        old_history <- self$history  # Keep current history
        self$history <- loaded$history
        
        # Add load record to history
        load_entry <- list(
          time = Sys.time(),
          action = "load",
          description = sprintf("Loaded data flow from file '%s'", file)
        )
        self$history <- c(self$history, list(load_entry), old_history)
        
        # Update last modified time
        self$metadata$last_modified <- Sys.time()
        self$metadata$loaded_time <- Sys.time()
        self$metadata$loaded_file <- file
        
      } else {
        warning(sprintf("File does not exist: '%s'", file))
      }
      
      return(invisible(self))
    },
    
    #' @description Log event to history
    #' @param action Action
    #' @param description Description
    #' @return Data flow object
    log_event = function(action, description) {
      entry <- list(
        time = Sys.time(),
        action = action,
        description = description
      )
      
      self$history <- c(self$history, list(entry))
      
      return(invisible(self))
    },
    
    #' @description Print history
    #' @param n Number of recent events to display, NULL to display all
    #' @return Data flow object
    print_history = function(n = NULL) {
      history_length <- length(self$history)
      
      if(is.null(n)) {
        n <- history_length
      } else {
        n <- min(n, history_length)
      }
      
      if(n > 0) {
        start_idx <- history_length - n + 1
        for(i in start_idx:history_length) {
          entry <- self$history[[i]]
          timestamp <- format(entry$time, "%Y-%m-%d %H:%M:%S")
          cat(sprintf("[%s] %s - %s\n", timestamp, entry$action, entry$description))
        }
      }
      
      return(invisible(self))
    },
    
    #' @description Get data summary
    #' @return Data summary data frame
    summary = function() {
      if(length(self$data) == 0) {
        cat("Data flow is empty\n")
        return(invisible(data.frame()))
      }
      
      summary_df <- data.frame(
        Name = character(),
        Type = character(),
        Size = character(),
        Stage = character(),
        AddedTime = character(),
        stringsAsFactors = FALSE
      )
      
      for(name in names(self$data)) {
        md <- self$metadata[[name]]
        if(is.null(md)) next
        
        summary_df <- rbind(summary_df, data.frame(
          Name = name,
          Type = md$class,
          Size = format(md$size, units = "auto"),
          Stage = md$stage,
          AddedTime = format(md$added_time, "%Y-%m-%d %H:%M:%S"),
          stringsAsFactors = FALSE
        ))
      }
      
      return(summary_df)
    },
    
    #' @description Print object
    #' @return Data flow object
    print = function() {
      cat("Ancestral Chromosome Analysis Data Flow [", format(self$metadata$creation_time, "%Y-%m-%d %H:%M:%S"), "]\n")
      cat("Data items: ", length(self$data), "\n")
      cat("History records: ", length(self$history), "\n")
      
      if(length(self$data) > 0) {
        cat("\nData Summary:\n")
        print(self$summary())
      }
      
      return(invisible(self))
    }
  )
)

#===============================================================================
# Analysis Pipeline Management Class
#===============================================================================

#' Analysis Pipeline Management Class
#' 
#' Class for creating and executing analysis pipelines
#' @export
AnalysisPipeline <- R6::R6Class(
  "AnalysisPipeline",
  public = list(
    #' @field name Pipeline name
    name = NULL,
    
    #' @field config Configuration
    config = NULL,
    
    #' @field data_flow Data flow object
    data_flow = NULL,
    
    #' @field stages Analysis stages
    stages = list(),
    
    #' @field current_stage Current stage
    current_stage = NULL,
    
    #' @field status Pipeline status
    status = "initialized",
    
    #' @field start_time Start time
    start_time = NULL,
    
    #' @field end_time End time
    end_time = NULL,
    
    #' @description Initialize pipeline
    #' @param name Pipeline name
    #' @param config Configuration
    #' @param data_flow Data flow object
    #' @return Pipeline object
    initialize = function(name = "Ancestral Chromosome Analysis", config = NULL, data_flow = NULL) {
      self$name <- name
      
      # Initialize configuration
      self$config <- if(is.null(config)) load_config() else config
      
      # Initialize data flow
      self$data_flow <- if(is.null(data_flow)) DataFlow$new() else data_flow
      
      # Initialize stages list
      self$stages <- list()
      self$current_stage <- NULL
      self$status <- "initialized"
      
      message(sprintf("Analysis pipeline '%s' initialized", name))
      return(invisible(self))
    },
    
    #' @description Add analysis stage
    #' @param name Stage name
    #' @param func Stage function
    #' @param dependencies Dependencies stages
    #' @param description Stage description
    #' @return Pipeline object
    add_stage = function(name, func, dependencies = NULL, description = "") {
      if(name %in% names(self$stages)) {
        warning(sprintf("Stage '%s' already exists, will be overwritten", name))
      }
      
      self$stages[[name]] <- list(
        name = name,
        func = func,
        dependencies = dependencies,
        description = description,
        status = "pending",
        start_time = NULL,
        end_time = NULL,
        duration = NULL,
        error = NULL
      )
      
      message(sprintf("Added analysis stage: '%s'", name))
      if(!is.null(dependencies) && length(dependencies) > 0) {
        message(sprintf("  Dependencies: %s", paste(dependencies, collapse = ", ")))
      }
      
      return(invisible(self))
    },
    
    #' @description Run analysis stage
    #' @param name Stage name
    #' @param force Whether to force rerun
    #' @return Pipeline object
    run_stage = function(name, force = FALSE) {
      if(!name %in% names(self$stages)) {
        stop(sprintf("Stage does not exist: '%s'", name))
      }
      
      stage <- self$stages[[name]]
      
      # Skip if stage is already completed and not forcing rerun
      if(stage$status == "completed" && !force) {
        message(sprintf("Stage '%s' already completed, skipping (use force=TRUE to rerun)", name))
        return(invisible(self))
      }
      
      # Check dependencies
      if(!is.null(stage$dependencies)) {
        for(dep in stage$dependencies) {
          if(!dep %in% names(self$stages)) {
            stop(sprintf("Dependency stage does not exist: '%s'", dep))
          }
          
          dep_stage <- self$stages[[dep]]
          if(dep_stage$status != "completed") {
            message(sprintf("Running dependency stage: '%s'", dep))
            self$run_stage(dep, force = force)
          }
        }
      }
      
      # Run stage
      self$current_stage <- name
      message(sprintf("\n===== Running Stage: '%s' =====", name))
      if(stage$description != "") {
        message(sprintf("Description: %s", stage$description))
      }
      
      # Update status
      self$stages[[name]]$status <- "running"
      self$stages[[name]]$start_time <- Sys.time()
      self$status <- "running"
      
      # Execute stage function
      result <- tryCatch({
        # Call stage function, passing config and data flow
        stage$func(self$config, self$data_flow)
      }, error = function(e) {
        self$stages[[name]]$status <- "failed"
        self$stages[[name]]$error <- e$message
        self$stages[[name]]$end_time <- Sys.time()
        self$stages[[name]]$duration <- difftime(self$stages[[name]]$end_time, 
                                               self$stages[[name]]$start_time, 
                                               units = "secs")
        
        message(sprintf("Stage execution failed: '%s'", name))
        message(sprintf("Error message: %s", e$message))
        
        if(self$config$advanced$debug) {
          message("Stack trace:")
          print(sys.calls())
        }
        
        stop(sprintf("Stage '%s' execution failed: %s", name, e$message))
      })
      
      # Update stage status
      self$stages[[name]]$status <- "completed"
      self$stages[[name]]$end_time <- Sys.time()
      self$stages[[name]]$duration <- difftime(self$stages[[name]]$end_time, 
                                             self$stages[[name]]$start_time, 
                                             units = "secs")
      
      message(sprintf("Stage '%s' completed, duration: %s", name, 
                     format(self$stages[[name]]$duration, digits = 2)))
      
      return(invisible(self))
    },
    
    #' @description Run all analysis stages
    #' @param force Whether to force rerun
    #' @return Pipeline object
    run_all = function(force = FALSE) {
      message(sprintf("\n====== Starting Analysis Pipeline: '%s' ======\n", self$name))
      
      self$start_time <- Sys.time()
      self$status <- "running"
      
      # Run all stages in order
      stages_order <- self$get_execution_order()
      
      for(name in stages_order) {
        tryCatch({
          self$run_stage(name, force = force)
        }, error = function(e) {
          message(sprintf("Pipeline terminated at stage '%s': %s", name, e$message))
          self$status <- "failed"
        })
        
        # Stop execution if pipeline has failed
        if(self$status == "failed") {
          break
        }
      }
      
      # Update pipeline status
      self$end_time <- Sys.time()
      if(self$status == "running") {
        self$status <- "completed"
      }
      
      total_time <- difftime(self$end_time, self$start_time, units = "mins")
      
      message(sprintf("\n====== Analysis Pipeline %s, Total Duration: %.2f minutes ======\n", 
                     ifelse(self$status == "completed", "Completed", "Failed"), 
                     as.numeric(total_time)))
      
      return(invisible(self))
    },
    
    #' @description Get stage execution order
    #' @return Vector of stage names in execution order
    get_execution_order = function() {
      # Create directed graph to represent dependencies
      edges <- list()
      nodes <- names(self$stages)
      
      for(name in nodes) {
        deps <- self$stages[[name]]$dependencies
        if(!is.null(deps)) {
          for(dep in deps) {
            edges <- c(edges, list(c(dep, name)))
          }
        }
      }
      
      # Topological sort
      visited <- setNames(rep(FALSE, length(nodes)), nodes)
      temp_mark <- setNames(rep(FALSE, length(nodes)), nodes)
      ordered <- character(0)
      
      visit <- function(node) {
        if(temp_mark[node]) {
          stop("Circular dependency detected")
        }
        
        if(!visited[node]) {
          temp_mark[node] <- TRUE
          
          # Process dependencies
          deps <- self$stages[[node]]$dependencies
          if(!is.null(deps)) {
            for(dep in deps) {
              visit(dep)
            }
          }
          
          temp_mark[node] <- FALSE
          visited[node] <- TRUE
          ordered <- c(node, ordered)
        }
        
        return(ordered)
      }
      
      # Apply DFS to each unvisited node
      for(node in nodes) {
        if(!visited[node]) {
          ordered <- visit(node)
        }
      }
      
      return(rev(ordered))
    },
    
    #' @description Get pipeline status summary
    #' @return Status summary data frame
    get_status = function() {
      stages_df <- data.frame(
        Stage = character(),
        Status = character(),
        Dependencies = character(),
        StartTime = character(),
        EndTime = character(),
        Duration = character(),
        stringsAsFactors = FALSE
      )
      
      for(name in names(self$stages)) {
        stage <- self$stages[[name]]
        deps <- if(is.null(stage$dependencies)) "" else paste(stage$dependencies, collapse = ", ")
        
        start_time <- if(is.null(stage$start_time)) "" else format(stage$start_time, "%H:%M:%S")
        end_time <- if(is.null(stage$end_time)) "" else format(stage$end_time, "%H:%M:%S")
        duration <- if(is.null(stage$duration)) "" else sprintf("%.2f sec", as.numeric(stage$duration))
        
        stages_df <- rbind(stages_df, data.frame(
          Stage = name,
          Status = stage$status,
          Dependencies = deps,
          StartTime = start_time,
          EndTime = end_time,
          Duration = duration,
          stringsAsFactors = FALSE
        ))
      }
      
      return(stages_df)
    },
    
    #' @description Print pipeline status
    #' @return Pipeline object
    print_status = function() {
      cat(sprintf("Analysis Pipeline: %s\n", self$name))
      cat(sprintf("Status: %s\n", self$status))
      
      if(!is.null(self$start_time)) {
        cat(sprintf("Start Time: %s\n", format(self$start_time, "%Y-%m-%d %H:%M:%S")))
      }
      
      if(!is.null(self$end_time)) {
        cat(sprintf("End Time: %s\n", format(self$end_time, "%Y-%m-%d %H:%M:%S")))
        
        if(!is.null(self$start_time)) {
          total_time <- difftime(self$end_time, self$start_time, units = "mins")
          cat(sprintf("Total Duration: %.2f minutes\n", as.numeric(total_time)))
        }
      }
      
      cat("\nStage Status:\n")
      print(self$get_status())
      
      return(invisible(self))
    },
    
    #' @description Print object
    #' @return Pipeline object
    print = function() {
      cat(sprintf("Ancestral Chromosome Analysis Pipeline: %s\n", self$name))
      cat(sprintf("Status: %s\n", self$status))
      cat(sprintf("Number of stages: %d\n", length(self$stages)))
      
      if(self$status != "initialized") {
        if(!is.null(self$start_time)) {
          cat(sprintf("Start time: %s\n", format(self$start_time, "%Y-%m-%d %H:%M:%S")))
        }
        
        completed <- sum(sapply(self$stages, function(s) s$status == "completed"))
        cat(sprintf("Completed stages: %d/%d\n", completed, length(self$stages)))
      }
      
      return(invisible(self))
    }
  )
)

#===============================================================================
# Helper Functions
#===============================================================================

#' Check and install required packages
#'
#' @param packages Vector of packages to check
#' @param quiet Whether to install silently
#' @return Vector of successfully installed packages
#' @export
ensure_packages <- function(packages, quiet = FALSE) {
  # Check which packages need installation
  new_packages <- packages[!packages %in% installed.packages()[,"Package"]]
  
  # Install missing packages
  if(length(new_packages) > 0) {
    if(!quiet) {
      message("Installing missing R packages: ", paste(new_packages, collapse = ", "))
    }
    install.packages(new_packages, repos = "https://cloud.r-project.org")
  }
  
  # Check if installation was successful
  still_missing <- new_packages[!new_packages %in% installed.packages()[,"Package"]]
  if(length(still_missing) > 0) {
    warning("Failed to install the following packages: ", paste(still_missing, collapse = ", "))
  }
  
  # Load all packages
  for(pkg in packages) {
    if(pkg %in% installed.packages()[,"Package"]) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }
  }
  
  # Return successfully installed packages
  return(setdiff(packages, still_missing))
}

#' Create complete analysis output directory structure
#'
#' @param base_dir Base directory
#' @param create Whether to create directories
#' @return Directory structure list
#' @export
create_output_dirs <- function(base_dir, create = TRUE) {
  # Define directory structure
  dirs <- list(
    base = base_dir,
    data = file.path(base_dir, "data"),
    phases = list(
      phase1 = file.path(base_dir, "phase1_data_integration"),
      phase2 = file.path(base_dir, "phase2_conserved_groups"),
      phase3 = file.path(base_dir, "phase3_chromosome_events"),
      phase4 = file.path(base_dir, "phase4_ancestral_reconstruction")
    ),
    methods = list(
      parsimony = file.path(base_dir, "methods", "parsimony"),
      bayesian = file.path(base_dir, "methods", "bayesian"),
      ml = file.path(base_dir, "methods", "maximum_likelihood"),
      ensemble = file.path(base_dir, "methods", "ensemble")
    ),
    plots = file.path(base_dir, "plots"),
    reports = file.path(base_dir, "reports"),
    logs = file.path(base_dir, "logs")
  )
  
  # Create directories
  if(create) {
    # Create base directory
    if(!dir.exists(dirs$base)) {
      dir.create(dirs$base, recursive = TRUE)
    }
    
    # Create data directory
    if(!dir.exists(dirs$data)) {
      dir.create(dirs$data, recursive = TRUE)
    }
    
    # Create phase directories
    for(phase in names(dirs$phases)) {
      if(!dir.exists(dirs$phases[[phase]])) {
        dir.create(dirs$phases[[phase]], recursive = TRUE)
      }
    }
    
    # Create method directories
    for(method in names(dirs$methods)) {
      if(!dir.exists(dirs$methods[[method]])) {
        dir.create(dirs$methods[[method]], recursive = TRUE)
      }
    }
    
    # Create other directories
    for(dir_name in c("plots", "reports", "logs")) {
      if(!dir.exists(dirs[[dir_name]])) {
        dir.create(dirs[[dir_name]], recursive = TRUE)
      }
    }
    
    message("Created output directory structure: ", base_dir)
  }
  
  return(dirs)
}
