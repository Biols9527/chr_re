# ...existing code...

#' Summarize evolutionary model test results
#' 
#' @param model_results Model test analysis results
#' @return Summary of model results
#' @keywords internal
summarize_model_results <- function(model_results) {
  summary <- list()
  
  # Best model and support
  if(!is.null(model_results$model_comparison)) {
    comp <- model_results$model_comparison
    best_idx <- which.min(comp$AICc)
    
    summary$best_model <- comp$Model[best_idx]
    summary$best_model_weight <- comp$AICc_weight[best_idx]
    
    # Check if there's a clear winner
    if(comp$AICc_weight[best_idx] > 0.9) {
      summary$model_support <- "Strong support for best model"
    } else if(comp$AICc_weight[best_idx] > 0.7) {
      summary$model_support <- "Moderate support for best model"
    } else {
      summary$model_support <- "Weak support for best model"
    }
    
    # Model interpretation
    if(summary$best_model == "BM1") {
      summary$model_interpretation <- "Random-walk evolution (Brownian Motion)"
    } else if(summary$best_model == "OU1") {
      summary$model_interpretation <- "Stabilizing selection toward optimum (Ornstein-Uhlenbeck)"
      
      # Extract optimum value
      if(!is.null(model_results$model_fits$OU1$opt$z0)) {
        summary$ou_optimum <- model_results$model_fits$OU1$opt$z0
      }
    } else if(summary$best_model == "EB") {
      eb_rate <- model_results$model_fits$EB$opt$a
      if(eb_rate < 0) {
        summary$model_interpretation <- "Early burst of evolution followed by slowdown"
      } else {
        summary$model_interpretation <- "Accelerating evolution through time"
      }
    } else if(summary$best_model == "trend") {
      summary$model_interpretation <- "Directional trend in chromosome evolution"
      
      # Extract trend direction and magnitude
      if(!is.null(model_results$model_fits$trend$opt$drift)) {
        trend <- model_results$model_fits$trend$opt$drift
        summary$trend_coefficient <- trend
        summary$trend_direction <- ifelse(trend > 0, "increasing", "decreasing")
      }
    } else if(summary$best_model == "drift") {
      summary$model_interpretation <- "Directional selection with stabilizing component"
      
      # Extract drift direction
      if(!is.null(model_results$model_fits$drift$opt$drift)) {
        drift <- model_results$model_fits$drift$opt$drift
        summary$drift_coefficient <- drift
        summary$drift_direction <- ifelse(drift > 0, "increasing", "decreasing")
      }
    } else if(summary$best_model == "white") {
      summary$model_interpretation <- "No phylogenetic signal (white noise)"
    }
  }
  
  return(summary)
}

#' Summarize all test results
#' 
#' @param test_results List of test results from different analyses
#' @return Summary of all test results
#' @keywords internal
summarize_test_results <- function(test_results) {
  # Initialize summary
  summary <- list(
    key_findings = character(0)
  )
  
  # Add trait correlation findings
  if(!is.null(test_results$trait_correlation)) {
    if(!is.null(test_results$trait_correlation$summary)) {
      corr_summary <- test_results$trait_correlation$summary
      
      if(corr_summary$n_significant > 0) {
        top_traits <- corr_summary$significant_traits[1:min(3, length(corr_summary$significant_traits))]
        top_traits_str <- paste(top_traits, collapse = ", ")
        
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Found %d significant trait correlations (top traits: %s)",
                                       corr_summary$n_significant, top_traits_str))
      } else {
        summary$key_findings <- c(summary$key_findings,
                               "No significant trait correlations detected")
      }
    }
  }
  
  # Add disparity findings
  if(!is.null(test_results$disparity)) {
    if(!is.null(test_results$disparity$summary)) {
      disp_summary <- test_results$disparity$summary
      
      if(!is.null(disp_summary$dtt_pattern)) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Disparity through time: %s", disp_summary$dtt_pattern))
      }
      
      if(!is.null(disp_summary$n_significant_comparisons) && 
         disp_summary$n_significant_comparisons > 0) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Found %d significant clade disparity differences",
                                      disp_summary$n_significant_comparisons))
      }
    }
  }
  
  # Add tempo and mode findings
  if(!is.null(test_results$tempo_mode)) {
    if(!is.null(test_results$tempo_mode$summary)) {
      tm_summary <- test_results$tempo_mode$summary
      
      if(!is.null(tm_summary$tempo_pattern)) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Tempo: %s", tm_summary$tempo_pattern))
      }
      
      if(!is.null(tm_summary$mode_pattern)) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Mode: %s", tm_summary$mode_pattern))
      }
      
      if(!is.null(tm_summary$node_height_pattern)) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Node-height test: %s", tm_summary$node_height_pattern))
      }
    }
  }
  
  # Add phylogenetic signal findings
  if(!is.null(test_results$phylo_signal)) {
    if(!is.null(test_results$phylo_signal$summary)) {
      sig_summary <- test_results$phylo_signal$summary
      
      if(!is.null(sig_summary$overall_assessment)) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Phylogenetic signal: %s", sig_summary$overall_assessment))
      }
      
      # Add specific signal metrics if significant
      if(!is.null(sig_summary$lambda_significant) && sig_summary$lambda_significant) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Pagel's lambda = %.2f (p = %.3f)",
                                      sig_summary$lambda, sig_summary$lambda_p_value))
      }
      
      if(!is.null(sig_summary$k_significant) && sig_summary$k_significant) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Blomberg's K = %.2f (p = %.3f)",
                                      sig_summary$k, sig_summary$k_p_value))
      }
    }
  }
  
  # Add evolutionary model findings
  if(!is.null(test_results$models)) {
    if(!is.null(test_results$models$summary)) {
      model_summary <- test_results$models$summary
      
      if(!is.null(model_summary$best_model)) {
        summary$key_findings <- c(summary$key_findings,
                               sprintf("Best-fit model: %s (weight = %.2f)",
                                      model_summary$best_model, model_summary$best_model_weight))
        
        if(!is.null(model_summary$model_interpretation)) {
          summary$key_findings <- c(summary$key_findings,
                                 sprintf("Model interpretation: %s", model_summary$model_interpretation))
        }
      }
    }
  }
  
  return(summary)
}

#' Create visualizations of test results
#' 
#' @param results Results from test_chromosome_evolution
#' @return List of visualization plots
#' @keywords internal
create_test_visualizations <- function(results) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 package is required for visualizations")
    return(list())
  }
  
  plots <- list()
  
  # Create visualizations for trait correlations
  if(!is.null(results$test_results$trait_correlation)) {
    trait_results <- results$test_results$trait_correlation
    
    if(length(trait_results$correlations) > 0) {
      # Extract significant correlations
      sig_traits <- character(0)
      sig_estimates <- numeric(0)
      sig_pvals <- numeric(0)
      
      for(trait_name in names(trait_results$correlations)) {
        if(trait_results$correlations[[trait_name]]$p_value < 0.05) {
          sig_traits <- c(sig_traits, trait_name)
          sig_estimates <- c(sig_estimates, trait_results$correlations[[trait_name]]$estimate)
          sig_pvals <- c(sig_pvals, trait_results$correlations[[trait_name]]$p_value)
        }
      }
      
      if(length(sig_traits) > 0) {
        # Create correlation plot for significant traits
        corr_data <- data.frame(
          Trait = sig_traits,
          Estimate = sig_estimates,
          P_Value = sig_pvals,
          stringsAsFactors = FALSE
        )
        
        corr_data$Significance <- ifelse(corr_data$P_Value < 0.01, "**", 
                                       ifelse(corr_data$P_Value < 0.05, "*", ""))
        
        plots$trait_correlations <- ggplot2::ggplot(corr_data, 
                                                  ggplot2::aes(x = reorder(Trait, -abs(Estimate)), y = Estimate, fill = Estimate)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::geom_text(ggplot2::aes(label = Significance, y = Estimate + sign(Estimate) * 0.05), vjust = 0) +
          ggplot2::labs(title = "Significant Trait Correlations with Chromosome Number",
                     subtitle = "* p < 0.05, ** p < 0.01",
                     x = "Trait",
                     y = "Correlation Coefficient") +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      }
    }
  }
  
  # Create disparity visualization
  if(!is.null(results$test_results$disparity)) {
    disp_results <- results$test_results$disparity
    
    # Create disparity through time plot if available
    if(!is.null(disp_results$dtt)) {
      dtt <- disp_results$dtt
      
      # Convert to data frame for plotting
      dtt_data <- data.frame(
        Time = dtt$times,
        Disparity = dtt$dtt,
        stringsAsFactors = FALSE
      )
      
      # Create simulation data
      sim_data <- data.frame()
      
      if(!is.null(dtt$sim)) {
        for(i in 1:ncol(dtt$sim)) {
          sim_data <- rbind(sim_data, data.frame(
            Time = dtt$times,
            Disparity = dtt$sim[, i],
            Simulation = i,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Create plot
      if(nrow(sim_data) > 0) {
        # With simulation envelope
        plots$disparity_through_time <- ggplot2::ggplot() +
          ggplot2::geom_line(data = sim_data, ggplot2::aes(x = Time, y = Disparity, group = Simulation), 
                          color = "gray80", alpha = 0.3) +
          ggplot2::geom_line(data = dtt_data, ggplot2::aes(x = Time, y = Disparity), 
                          color = "blue", size = 1) +
          ggplot2::labs(title = "Disparity Through Time",
                     subtitle = sprintf("MDI = %.3f (p = %.3f)", dtt$MDI, dtt$MDIp),
                     x = "Relative Time",
                     y = "Relative Disparity") +
          ggplot2::theme_minimal()
      } else {
        # Without simulation envelope
        plots$disparity_through_time <- ggplot2::ggplot(dtt_data, ggplot2::aes(x = Time, y = Disparity)) +
          ggplot2::geom_line(color = "blue", size = 1) +
          ggplot2::labs(title = "Disparity Through Time",
                     subtitle = sprintf("MDI = %.3f", dtt$MDI),
                     x = "Relative Time",
                     y = "Relative Disparity") +
          ggplot2::theme_minimal()
      }
    }
    
    # Create clade disparity comparison plot if available
    if(!is.null(disp_results$clade_disparities) && length(disp_results$clade_disparities) > 1) {
      # Extract variance values from each clade
      clade_names <- names(disp_results$clade_disparities)
      clade_vars <- sapply(disp_results$clade_disparities, function(x) x$variance)
      
      # Create data frame
      clade_data <- data.frame(
        Clade = clade_names,
        Variance = clade_vars,
        stringsAsFactors = FALSE
      )
      
      # Create plot
      plots$clade_disparity <- ggplot2::ggplot(clade_data, 
                                            ggplot2::aes(x = reorder(Clade, -Variance), y = Variance, fill = Clade)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = "Chromosome Number Disparity by Clade",
                   x = "Clade",
                   y = "Variance") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  }
  
  # Create tempo and mode visualization
  if(!is.null(results$test_results$tempo_mode)) {
    tm_results <- results$test_results$tempo_mode
    
    # Create node height test plot if available
    if(!is.null(tm_results$node_height_test)) {
      nh_test <- tm_results$node_height_test
      
      # Create data frame for plotting
      nh_data <- nh_test$node_data
      
      # Create plot
      plots$node_height_test <- ggplot2::ggplot(nh_data, ggplot2::aes(x = height, y = contrast)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "blue") +
        ggplot2::labs(title = "Node-Height Test",
                   subtitle = sprintf("Slope = %.3f, p = %.3f", nh_test$slope, nh_test$p_value),
                   x = "Node Height",
                   y = "Absolute Contrast") +
        ggplot2::theme_minimal()
    }
  }
  
  # Create phylogenetic signal visualization
  if(!is.null(results$test_results$phylo_signal)) {
    sig_results <- results$test_results$phylo_signal
    
    # Create bar plot of signal metrics
    if(!is.null(sig_results$signal_metrics)) {
      metrics <- sig_results$signal_metrics
      
      # Extract relevant metrics
      metric_names <- c()
      metric_values <- c()
      
      if(!is.null(metrics$K)) {
        metric_names <- c(metric_names, "Blomberg's K")
        metric_values <- c(metric_values, metrics$K)
      }
      
      if(!is.null(metrics$lambda)) {
        metric_names <- c(metric_names, "Pagel's lambda")
        metric_values <- c(metric_values, metrics$lambda)
      }
      
      if(!is.null(metrics$morans_i)) {
        metric_names <- c(metric_names, "Moran's I")
        metric_values <- c(metric_values, metrics$morans_i)
      }
      
      if(length(metric_names) > 0) {
        # Create data frame
        signal_data <- data.frame(
          Metric = metric_names,
          Value = metric_values,
          stringsAsFactors = FALSE
        )
        
        # Create plot
        plots$signal_metrics <- ggplot2::ggplot(signal_data, 
                                             ggplot2::aes(x = Metric, y = Value, fill = Metric)) +
          ggplot2::geom_bar(stat = "identity") +
          ggplot2::labs(title = "Phylogenetic Signal Metrics",
                     subtitle = "Higher values indicate stronger signal",
                     x = "Metric",
                     y = "Value") +
          ggplot2::theme_minimal()
      }
    }
  }
  
  # Create model comparison visualization
  if(!is.null(results$test_results$models)) {
    model_results <- results$test_results$models
    
    # Create model weight plot if available
    if(!is.null(model_results$model_comparison)) {
      comp <- model_results$model_comparison
      
      # Create plot
      plots$model_comparison <- ggplot2::ggplot(comp, 
                                             ggplot2::aes(x = reorder(Model, -AICc_weight), y = AICc_weight, fill = Model)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::labs(title = "Model Comparison",
                   subtitle = "AICc Weights",
                   x = "Model",
                   y = "AICc Weight") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  }
  
  return(plots)
}

#===============================================================================
# Workflow Functions
#===============================================================================

#' Run complete statistical testing workflow for chromosome evolution
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Chromosome counts data
#' @param traits Optional trait data for correlation tests
#' @param output_dir Directory to save results
#' @param tests Tests to run: "trait_correlation", "disparity", "tempo_mode",
#'              "phylo_signal", "models", or "all"
#' @param generate_report Whether to generate HTML report
#' @param run_simulations Whether to run simulation tests for best model
#' @return Results of statistical testing workflow
#' @export
run_statistical_testing_workflow <- function(tree, 
                                          chr_counts, 
                                          traits = NULL,
                                          output_dir = "statistical_tests",
                                          tests = "all",
                                          generate_report = TRUE,
                                          run_simulations = TRUE) {
  # Create output directory if it doesn't exist
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Step 1: Run statistical tests
  message("Step 1: Running statistical tests for chromosome evolution...")
  test_results <- test_chromosome_evolution(
    tree = tree,
    chr_counts = chr_counts,
    traits = traits,
    test_type = tests
  )
  
  # Save test results
  saveRDS(test_results, file.path(output_dir, "chromosome_test_results.rds"))
  
  # Save data tables for key results
  save_data_tables(test_results, output_dir)
  
  # Step 2: Run simulations for the best model if requested
  if(run_simulations && !is.null(test_results$test_results$models)) {
    message("Step 2: Running simulations for the best model...")
    best_model <- test_results$test_results$models$best_model
    
    if(!is.na(best_model)) {
      model_sim_results <- NULL
      
      # Simulate under the best model if possible
      if(exists("simulate_chromosome_evolution", mode = "function")) {
        # Extract parameters from the best model
        model_params <- extract_model_parameters(test_results$test_results$models, best_model)
        
        # Run simulation
        model_sim_results <- simulate_chromosome_evolution(
          tree = tree,
          model = model_name_to_sim_model(best_model),
          params = model_params,
          root_value = mean(chr_counts),
          discrete = TRUE
        )
        
        # Save simulation results
        saveRDS(model_sim_results, file.path(output_dir, "best_model_simulation.rds"))
        
        # Create and save simulation visualizations
        if(requireNamespace("ggplot2", quietly = TRUE)) {
          sim_plots <- create_simulation_plots(model_sim_results, chr_counts)
          
          for(plot_name in names(sim_plots)) {
            ggplot2::ggsave(
              filename = file.path(output_dir, paste0("sim_", plot_name, ".pdf")),
              plot = sim_plots[[plot_name]],
              width = 8,
              height = 6
            )
          }
        }
      } else {
        message("simulate_chromosome_evolution function not available. Skipping simulation step.")
      }
      
      test_results$simulation_results <- model_sim_results
    } else {
      message("No best model identified. Skipping simulation step.")
    }
  }
  
  # Step 3: Create and save visualizations
  message("Step 3: Creating visualizations...")
  plots <- create_test_visualizations(test_results)
  
  # Save plots
  if(length(plots) > 0) {
    for(plot_name in names(plots)) {
      if(requireNamespace("ggplot2", quietly = TRUE)) {
        ggplot2::ggsave(
          filename = file.path(output_dir, paste0(plot_name, ".pdf")),
          plot = plots[[plot_name]],
          width = 8,
          height = 6
        )
      }
    }
  }
  
  # Step 4: Generate report if requested
  if(generate_report && requireNamespace("rmarkdown", quietly = TRUE)) {
    message("Step 4: Generating HTML report...")
    
    # Create temporary Rmd file
    rmd_template <- create_report_template(test_results)
    rmd_file <- file.path(output_dir, "report_template.Rmd")
    writeLines(rmd_template, rmd_file)
    
    # Render report
    rmarkdown::render(
      input = rmd_file,
      output_file = "chromosome_evolution_report.html",
      output_dir = output_dir,
      quiet = TRUE
    )
    
    # Clean up temporary Rmd file
    if(file.exists(rmd_file)) {
      file.remove(rmd_file)
    }
  }
  
  message("Statistical testing workflow complete.")
  message(paste0("Results saved to: ", output_dir))
  
  return(test_results)
}

#' Save data tables from test results
#' 
#' @param test_results Results from test_chromosome_evolution
#' @param output_dir Directory to save data tables
#' @keywords internal
save_data_tables <- function(test_results, output_dir) {
  # Create tables directory
  tables_dir <- file.path(output_dir, "tables")
  if(!dir.exists(tables_dir)) {
    dir.create(tables_dir)
  }
  
  # Save trait correlation tables
  if(!is.null(test_results$test_results$trait_correlation)) {
    if(!is.null(test_results$test_results$trait_correlation$summary)) {
      corr_summary <- test_results$test_results$trait_correlation$summary
      
      if(corr_summary$n_significant > 0) {
        # Create table of significant correlations
        sig_traits <- data.frame(
          Trait = character(0),
          Estimate = numeric(0),
          P_Value = numeric(0),
          stringsAsFactors = FALSE
        )
        
        for(trait in corr_summary$significant_traits) {
          corr <- test_results$test_results$trait_correlation$correlations[[trait]]
          sig_traits <- rbind(sig_traits, data.frame(
            Trait = trait,
            Estimate = corr$estimate,
            P_Value = corr$p_value,
            stringsAsFactors = FALSE
          ))
        }
        
        # Save table
        write.csv(sig_traits, file.path(tables_dir, "significant_trait_correlations.csv"), row.names = FALSE)
      }
    }
  }
  
  # Save model comparison table
  if(!is.null(test_results$test_results$models)) {
    if(!is.null(test_results$test_results$models$model_comparison)) {
      comp <- test_results$test_results$models$model_comparison
      write.csv(comp, file.path(tables_dir, "model_comparison.csv"), row.names = FALSE)
    }
  }
  
  # Save phylogenetic signal results
  if(!is.null(test_results$test_results$phylo_signal)) {
    if(!is.null(test_results$test_results$phylo_signal$signal_metrics)) {
      metrics <- test_results$test_results$phylo_signal$signal_metrics
      
      # Convert to data frame
      metric_names <- c()
      metric_values <- c()
      
      for(name in names(metrics)) {
        metric_names <- c(metric_names, name)
        metric_values <- c(metric_values, metrics[[name]])
      }
      
      signal_df <- data.frame(
        Metric = metric_names,
        Value = metric_values,
        stringsAsFactors = FALSE
      )
      
      write.csv(signal_df, file.path(tables_dir, "phylogenetic_signal_metrics.csv"), row.names = FALSE)
    }
  }
  
  # Create a summary text file
  summary_file <- file.path(output_dir, "statistical_test_summary.txt")
  
  cat("Statistical Test Summary\n", file = summary_file)
  cat("=======================\n\n", file = summary_file, append = TRUE)
  
  # Add key findings
  if(!is.null(test_results$summary$key_findings)) {
    cat("Key Findings:\n", file = summary_file, append = TRUE)
    for(finding in test_results$summary$key_findings) {
      cat(sprintf("- %s\n", finding), file = summary_file, append = TRUE)
    }
    cat("\n", file = summary_file, append = TRUE)
  }
  
  # Add details for each test type
  for(test_type in names(test_results$test_results)) {
    if(!is.null(test_results$test_results[[test_type]]$summary)) {
      cat(sprintf("%s Analysis:\n", gsub("_", " ", test_type, fixed = TRUE)), 
         file = summary_file, append = TRUE)
      
      # Extract summary text based on the test type
      if(test_type == "trait_correlation") {
        # Add trait correlation summary
        if(test_results$test_results[[test_type]]$summary$n_significant > 0) {
          cat(sprintf("  Found %d significant trait correlations\n", 
                    test_results$test_results[[test_type]]$summary$n_significant), 
             file = summary_file, append = TRUE)
        } else {
          cat("  No significant trait correlations detected\n", file = summary_file, append = TRUE)
        }
      } else if(test_type == "models") {
        # Add model selection summary
        cat(sprintf("  Best model: %s\n", test_results$test_results[[test_type]]$summary$best_model), 
           file = summary_file, append = TRUE)
        cat(sprintf("  Model support: %s\n", test_results$test_results[[test_type]]$summary$model_support), 
           file = summary_file, append = TRUE)
        cat(sprintf("  Interpretation: %s\n", test_results$test_results[[test_type]]$summary$model_interpretation), 
           file = summary_file, append = TRUE)
      } else if(test_type == "phylo_signal") {
        # Add phylogenetic signal summary
        cat(sprintf("  Overall assessment: %s\n", 
                  test_results$test_results[[test_type]]$summary$overall_assessment), 
           file = summary_file, append = TRUE)
      } else if(test_type == "tempo_mode") {
        # Add tempo and mode summary
        if(!is.null(test_results$test_results[[test_type]]$summary$tempo_pattern)) {
          cat(sprintf("  Tempo: %s\n", test_results$test_results[[test_type]]$summary$tempo_pattern), 
             file = summary_file, append = TRUE)
        }
        if(!is.null(test_results$test_results[[test_type]]$summary$mode_pattern)) {
          cat(sprintf("  Mode: %s\n", test_results$test_results[[test_type]]$summary$mode_pattern), 
             file = summary_file, append = TRUE)
        }
      } else if(test_type == "disparity") {
        # Add disparity summary
        if(!is.null(test_results$test_results[[test_type]]$summary$dtt_pattern)) {
          cat(sprintf("  Disparity through time: %s\n", 
                    test_results$test_results[[test_type]]$summary$dtt_pattern), 
             file = summary_file, append = TRUE)
        }
        if(!is.null(test_results$test_results[[test_type]]$summary$n_significant_comparisons)) {
          cat(sprintf("  Significant clade differences: %d\n", 
                    test_results$test_results[[test_type]]$summary$n_significant_comparisons), 
             file = summary_file, append = TRUE)
        }
      }
      
      cat("\n", file = summary_file, append = TRUE)
    }
  }
}

#' Extract parameters from fitted model for simulation
#' 
#' @param model_results Model test results
#' @param best_model Name of best model
#' @return List of parameters for simulation
#' @keywords internal
extract_model_parameters <- function(model_results, best_model) {
  params <- list()
  
  # Get the model fit object
  if(!is.null(model_results$model_fits[[best_model]])) {
    fit <- model_results$model_fits[[best_model]]
    
    if(best_model == "BM1") {
      params$sigma <- sqrt(fit$parameters$sigma2)
      params$drift <- 0  # No trend in regular BM
    } else if(best_model == "OU1") {
      params$sigma <- sqrt(fit$parameters$sigma2)
      params$alpha <- fit$parameters$alpha
      params$theta <- fit$parameters$theta
    } else if(best_model == "EB") {
      params$sigma <- sqrt(fit$parameters$sigma2)
      params$rate_change <- fit$parameters$rate  # Beta parameter for ACDC
    } else if(best_model == "trend") {
      params$sigma <- sqrt(fit$parameters$sigma2)
      params$drift <- fit$parameters$drift
    } else if(best_model == "JumpBM" || best_model == "JumpOU") {
      params$sigma <- sqrt(fit$parameters$sigma2)
      params$jump_rate <- fit$parameters$jump_prob
      params$jump_size_mean <- fit$parameters$jump_size
      params$jump_size_sd <- fit$parameters$jump_size * 0.3  # Approx. based on typical ratio
    }
  }
  
  return(params)
}

#' Convert model name to simulation model type
#' 
#' @param model_name Name of model from test results
#' @return Simulation model name
#' @keywords internal
model_name_to_sim_model <- function(model_name) {
  # Map model names from test results to simulation model types
  sim_model <- switch(model_name,
                     "BM1" = "BM",
                     "OU1" = "OU",
                     "EB" = "ACDC",
                     "trend" = "BM",  # Use BM with trend parameter
                     "JumpBM" = "jumps",
                     "JumpOU" = "hybrid",
                     "BM")  # Default to BM
  
  return(sim_model)
}

#' Create plots for simulation results
#' 
#' @param sim_results Simulation results
#' @param real_data Real chromosome count data
#' @return List of plot objects
#' @keywords internal
create_simulation_plots <- function(sim_results, real_data) {
  plots <- list()
  
  # Create histogram comparison of real vs. simulated data
  if(!is.null(sim_results$tip_states)) {
    # Create data frame for plotting
    real_df <- data.frame(
      Count = real_data,
      Type = "Real Data",
      stringsAsFactors = FALSE
    )
    
    sim_df <- data.frame(
      Count = sim_results$tip_states,
      Type = "Simulated",
      stringsAsFactors = FALSE
    )
    
    combined_df <- rbind(real_df, sim_df)
    
    # Create histogram
    plots$count_histogram <- ggplot2::ggplot(combined_df, ggplot2::aes(x = Count, fill = Type)) +
      ggplot2::geom_histogram(position = "dodge", bins = 20, alpha = 0.7) +
      ggplot2::labs(title = "Real vs. Simulated Chromosome Counts",
                 subtitle = paste("Simulation model:", sim_results$model),
                 x = "Chromosome Count",
                 y = "Frequency") +
      ggplot2::theme_minimal()
    
    # Create density plot
    plots$count_density <- ggplot2::ggplot(combined_df, ggplot2::aes(x = Count, fill = Type)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::labs(title = "Distribution of Chromosome Counts",
                 subtitle = paste("Simulation model:", sim_results$model),
                 x = "Chromosome Count",
                 y = "Density") +
      ggplot2::theme_minimal()
  }
  
  return(plots)
}

#' Create report template for chromosome evolution statistics
#' 
#' @param test_results Results from test_chromosome_evolution
#' @return Character string with Rmd template
#' @keywords internal
create_report_template <- function(test_results) {
  # Create template header
  template <- c(
    "---",
    "title: \"Chromosome Evolution Statistical Analysis Report\"",
    "date: \"`r Sys.Date()`\"",
    "output: html_document",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(ggplot2)",
    "library(knitr)",
    "library(kableExtra)",
    "# Load results",
    "results <- readRDS(\"chromosome_test_results.rds\")",
    "```",
    "",
    "## Summary of Findings",
    "",
    "```{r summary}",
    "if(!is.null(results$summary$key_findings)) {",
    "  cat(\"<ul>\")",
    "  for(finding in results$summary$key_findings) {",
    "    cat(sprintf(\"<li>%s</li>\", finding))",
    "  }",
    "  cat(\"</ul>\")",
    "}",
    "```",
    ""
  )
  
  # Add sections for each test type
  test_types <- names(test_results$test_results)
  
  for(test_type in test_types) {
    section_title <- switch(test_type,
                          "trait_correlation" = "Trait Correlation Analysis",
                          "disparity" = "Disparity Analysis",
                          "tempo_mode" = "Tempo and Mode of Evolution",
                          "phylo_signal" = "Phylogenetic Signal",
                          "models" = "Evolutionary Model Comparison",
                          gsub("_", " ", test_type, fixed = TRUE))
    
    # Add section header
    template <- c(template, 
                 paste0("## ", section_title),
                 "")
    
    # Add section content based on test type
    template <- c(template, create_report_section(test_type))
  }
  
  # Add references section
  template <- c(template,
              "## References",
              "",
              "1. Harmon, L. J., Weir, J. T., Brock, C. D., Glor, R. E., & Challenger, W. (2008). GEIGER: investigating evolutionary radiations. Bioinformatics, 24(1), 129-131.",
              "2. Paradis, E., Claude, J., & Strimmer, K. (2004). APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20(2), 289-290.",
              "3. Revell, L. J. (2012). phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3(2), 217-223."
              )
  
  return(paste(template, collapse = "\n"))
}

#' Create report section for specific test type
#' 
#' @param test_type Type of test
#' @return Character vector with section content
#' @keywords internal
create_report_section <- function(test_type) {
  if(test_type == "trait_correlation") {
    return(c(
      "```{r trait-correlation}",
      "# Check if results are available",
      "if(!is.null(results$test_results$trait_correlation) && ",
      "   !is.null(results$test_results$trait_correlation$summary)) {",
      "  summary <- results$test_results$trait_correlation$summary",
      "  ",
      "  # Add summary text",
      "  cat(sprintf(\"Analysis found %d significant correlations between traits and chromosome numbers (out of %d tested). ",
      "             The analysis was conducted using %s method with %s transformation.\", ",
      "             summary$n_significant, summary$n_traits_tested,",
      "             summary$method, summary$transformation))",
      "  ",
      "  # Create table of significant correlations if any",
      "  if(summary$n_significant > 0) {",
      "    sig_traits <- data.frame(",
      "      Trait = character(0),",
      "      Estimate = numeric(0),",
      "      P_Value = numeric(0),",
      "      stringsAsFactors = FALSE",
      "    )",
      "    ",
      "    for(trait in summary$significant_traits) {",
      "      corr <- results$test_results$trait_correlation$correlations[[trait]]",
      "      sig_traits <- rbind(sig_traits, data.frame(",
      "        Trait = trait,",
      "        Estimate = corr$estimate,",
      "        P_Value = corr$p_value,",
      "        stringsAsFactors = FALSE",
      "      ))",
      "    }",
      "    ",
      "    # Show table",
      "    kable(sig_traits, caption = \"Significant Trait Correlations\", digits = 3) %>%",
      "      kable_styling(bootstrap_options = c(\"striped\", \"hover\"))",
      "  }",
      "} else {",
      "  cat(\"No trait correlation analysis results available.\")",
      "}",
      "```",
      "",
      "### Visualization",
      "",
      "```{r trait-correlation-plot, fig.width=8, fig.height=6}",
      "# Show correlation plot if available",
      "if(file.exists(\"trait_correlations.pdf\")) {",
      "  knitr::include_graphics(\"trait_correlations.pdf\")",
      "} else {",
      "  cat(\"No correlation plot available.\")",
      "}",
      "```",
      ""
    ))
  } else if(test_type == "disparity") {
    return(c(
      "```{r disparity}",
      "# Check if results are available",
      "if(!is.null(results$test_results$disparity) && ",
      "   !is.null(results$test_results$disparity$summary)) {",
      "  summary <- results$test_results$disparity$summary",
      "  ",
      "  # Add summary text",
      "  cat(sprintf(\"Overall disparity in chromosome numbers: %.2f (variance)\\n\", ",
      "             summary$overall_disparity))",
      "  ",
      "  if(!is.null(summary$mdi)) {",
      "    cat(sprintf(\"\\nDisparity through time analysis: MDI = %.3f (p = %.3f)\\n\", ",
      "               summary$mdi, summary$mdi_p))",
      "    cat(sprintf(\"Pattern: %s\\n\", summary$dtt_pattern))",
      "  }",
      "  ",
      "  if(!is.null(summary$n_clades) && summary$n_clades > 0) {",
      "    cat(sprintf(\"\\nAnalyzed disparity across %d clades\\n\", summary$n_clades))",
      "    ",
      "    if(!is.null(summary$significant_comparisons) && length(summary$significant_comparisons) > 0) {",
      "      cat(\"\\nSignificant clade differences:\\n\")",
      "      for(comp in summary$significant_comparisons) {",
      "        cat(sprintf(\"- %s\\n\", comp))",
      "      }",
      "    }",
      "  }",
      "} else {",
      "  cat(\"No disparity analysis results available.\")",
      "}",
      "```",
      "",
      "### Visualization",
      "",
      "```{r disparity-plots, fig.width=8, fig.height=6}",
      "# Show disparity plots if available",
      "if(file.exists(\"disparity_through_time.pdf\")) {",
      "  knitr::include_graphics(\"disparity_through_time.pdf\")",
      "}",
      "",
      "if(file.exists(\"clade_disparity.pdf\")) {",
      "  knitr::include_graphics(\"clade_disparity.pdf\")",
      "}",
      "```",
      ""
    ))
  } else if(test_type == "tempo_mode") {
    return(c(
      "```{r tempo-mode}",
      "# Check if results are available",
      "if(!is.null(results$test_results$tempo_mode) && ",
      "   !is.null(results$test_results$tempo_mode$summary)) {",
      "  summary <- results$test_results$tempo_mode$summary",
      "  ",
      "  # Add summary text",
      "  if(!is.null(summary$tempo_pattern)) {",
      "    cat(sprintf(\"Tempo: %s\\n\", summary$tempo_pattern))",
      "    if(!is.null(summary$eb_vs_bm_delta_aicc)) {",
      "      cat(sprintf(\"Early Burst vs. Brownian Motion: Delta AICc = %.2f\\n\", ",
      "                 summary$eb_vs_bm_delta_aicc))",
      "    }",
      "  }",
      "  ",
      "  if(!is.null(summary$mode_pattern)) {",
      "    cat(sprintf(\"\\nMode: %s\\n\", summary$mode_pattern))",
      "    if(!is.null(summary$mode_p_value)) {",
      "      cat(sprintf(\"Directionality test p-value: %.3f\\n\", summary$mode_p_value))",
      "    }",
      "  }",
      "  ",
      "  if(!is.null(summary$node_height_pattern)) {",
      "    cat(sprintf(\"\\nNode-height test: %s\\n\", summary$node_height_pattern))",
      "    if(!is.null(summary$node_height_p_value)) {",
      "      cat(sprintf(\"Node-height test p-value: %.3f\\n\", summary$node_height_p_value))",
      "    }",
      "  }",
      "} else {",
      "  cat(\"No tempo and mode analysis results available.\")",
      "}",
      "```",
      "",
      "### Visualization",
      "",
      "```{r tempo-mode-plot, fig.width=8, fig.height=6}",
      "# Show node height test plot if available",
      "if(file.exists(\"node_height_test.pdf\")) {",
      "  knitr::include_graphics(\"node_height_test.pdf\")",
      "}",
      "```",
      ""
    ))
  } else if(test_type == "phylo_signal") {
    return(c(
      "```{r phylo-signal}",
      "# Check if results are available",
      "if(!is.null(results$test_results$phylo_signal) && ",
      "   !is.null(results$test_results$phylo_signal$summary)) {",
      "  summary <- results$test_results$phylo_signal$summary",
      "  ",
      "  # Add summary text",
      "  cat(sprintf(\"Overall assessment: %s\\n\\n\", summary$overall_assessment))",
      "  ",
      "  # Show detailed metrics if available",
      "  metrics_table <- data.frame(",
      "    Metric = character(0),",
      "    Value = numeric(0),",
      "    P_Value = numeric(0),",
      "    Interpretation = character(0),",
      "    stringsAsFactors = FALSE",
      "  )",
      "  ",
      "  if(!is.null(summary$k)) {",
      "    metrics_table <- rbind(metrics_table, data.frame(",
      "      Metric = \"Blomberg's K\",",
      "      Value = summary$k,",
      "      P_Value = summary$k_p_value,",
      "      Interpretation = summary$k_interpretation,",
      "      stringsAsFactors = FALSE",
      "    ))",
      "  }",
      "  ",
      "  if(!is.null(summary$lambda)) {",
      "    metrics_table <- rbind(metrics_table, data.frame(",
      "      Metric = \"Pagel's lambda\",",
      "      Value = summary$lambda,",
      "      P_Value = summary$lambda_p_value,",
      "      Interpretation = summary$lambda_interpretation,",
      "      stringsAsFactors = FALSE",
      "    ))",
      "  }",
      "  ",
      "  if(!is.null(summary$morans_i)) {",
      "    metrics_table <- rbind(metrics_table, data.frame(",
      "      Metric = \"Moran's I\",",
      "      Value = summary$morans_i,",
      "      P_Value = summary$morans_i_p_value,",
      "      Interpretation = summary$morans_i_interpretation,",
      "      stringsAsFactors = FALSE",
      "    ))",
      "  }",
      "  ",
      "  if(nrow(metrics_table) > 0) {",
      "    kable(metrics_table, caption = \"Phylogenetic Signal Metrics\", digits = 3) %>%",
      "      kable_styling(bootstrap_options = c(\"striped\", \"hover\"))",
      "  }",
      "} else {",
      "  cat(\"No phylogenetic signal analysis results available.\")",
      "}",
      "```",
      "",
      "### Visualization",
      "",
      "```{r phylo-signal-plot, fig.width=8, fig.height=6}",
      "# Show signal metrics plot if available",
      "if(file.exists(\"signal_metrics.pdf\")) {",
      "  knitr::include_graphics(\"signal_metrics.pdf\")",
      "}",
      "```",
      ""
    ))
  } else if(test_type == "models") {
    return(c(
      "```{r models}",
      "# Check if results are available",
      "if(!is.null(results$test_results$models) && ",
      "   !is.null(results$test_results$models$summary)) {",
      "  summary <- results$test_results$models$summary",
      "  ",
      "  # Add summary text",
      "  cat(sprintf(\"Best-fit model: %s\\n\", summary$best_model))",
      "  cat(sprintf(\"Model support: %s\\n\", summary$model_support))",
      "  cat(sprintf(\"Interpretation: %s\\n\\n\", summary$model_interpretation))",
      "  ",
      "  # Show model-specific details if available",
      "  if(summary$best_model == \"OU1\" && !is.null(summary$ou_optimum)) {",
      "    cat(sprintf(\"Estimated optimum (theta): %.2f\\n\", summary$ou_optimum))",
      "  } else if(summary$best_model == \"trend\" && !is.null(summary$trend_coefficient)) {",
      "    cat(sprintf(\"Trend coefficient: %.3f (%s)\\n\", ",
      "               summary$trend_coefficient, summary$trend_direction))",
      "  } else if(summary$best_model == \"drift\" && !is.null(summary$drift_coefficient)) {",
      "    cat(sprintf(\"Drift coefficient: %.3f (%s)\\n\", ",
      "               summary$drift_coefficient, summary$drift_direction))",
      "  }",
      "  ",
      "  # Show model comparison table if available",
      "  if(!is.null(results$test_results$models$model_comparison)) {",
      "    comp <- results$test_results$models$model_comparison",
      "    comp_table <- comp[, c(\"Model\", \"LogLik\", \"Parameters\", \"AICc\", \"dAICc\", \"AICc_weight\")]",
      "    kable(comp_table, caption = \"Model Comparison\", digits = 3) %>%",
      "      kable_styling(bootstrap_options = c(\"striped\", \"hover\"))",
      "  }",
      "} else {",
      "  cat(\"No model comparison results available.\")",
      "}",
      "```",
      "",
      "### Visualization",
      "",
      "```{r model-comparison-plot, fig.width=8, fig.height=6}",
      "# Show model comparison plot if available",
      "if(file.exists(\"model_comparison.pdf\")) {",
      "  knitr::include_graphics(\"model_comparison.pdf\")",
      "}",
      "```",
      "",
      "```{r model-simulation-plot, fig.width=8, fig.height=6}",
      "# Show simulation plots if available",
      "if(file.exists(\"sim_count_histogram.pdf\")) {",
      "  knitr::include_graphics(\"sim_count_histogram.pdf\")",
      "}",
      "```",
      ""
    ))
  } else {
    # Generic section for other test types
    return(c(
      sprintf("```{r %s}", gsub("_", "-", test_type, fixed = TRUE)),
      "# Generic display for this test type",
      sprintf("if(!is.null(results$test_results$%s)) {", test_type),
      sprintf("  cat(\"Results available for %s analysis. See full results for details.\")", 
            gsub("_", " ", test_type, fixed = TRUE)),
      "} else {",
      sprintf("  cat(\"No %s analysis results available.\")", 
            gsub("_", " ", test_type, fixed = TRUE)),
      "}",
      "```",
      ""
    ))
  }
}
