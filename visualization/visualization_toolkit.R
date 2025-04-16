#===============================================================================
# Ancestral Chromosome Reconstruction Framework: Advanced Visualization Toolkit
# Author: Bioinformatics Team
# Date: 2025-05-10
# Description: Specialized visualization functions for ancestral chromosome
#              reconstruction, including interactive plots, comparative views,
#              and publication-ready figures
#===============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(ggplot2)
  library(ggtree)
  library(RColorBrewer)
  library(viridis)
  library(plotly)
  library(gridExtra)
  library(patchwork)
})

#===============================================================================
# Core Visualization Functions
#===============================================================================

#' Create interactive ancestral chromosome visualization
#' 
#' Generates an interactive plot showing ancestral chromosome states on a phylogeny
#' with tooltips and interactive elements
#' 
#' @param tree Phylogenetic tree
#' @param ancestral_reconstruction Ancestral reconstruction results
#' @param display_options List of display options
#' @param highlight_nodes Vector of nodes to highlight
#' @param node_data Optional data frame with additional node information
#' @param color_scheme Color scheme for states: "viridis", "rainbow", "divergent", etc.
#' @param discrete_states Whether chromosome counts should be treated as discrete
#' @param width Width of the plot
#' @param height Height of the plot
#' @return Interactive plotly object
#' @export
create_interactive_chromosomes <- function(tree, 
                                         ancestral_reconstruction, 
                                         display_options = list(
                                           show_node_labels = FALSE,
                                           show_confidence = TRUE,
                                           show_events = TRUE
                                         ),
                                         highlight_nodes = NULL,
                                         node_data = NULL,
                                         color_scheme = "viridis",
                                         discrete_states = TRUE,
                                         width = 800,
                                         height = 600) {
  # Check for required packages
  if(!requireNamespace("plotly", quietly = TRUE)) {
    stop("The 'plotly' package is required for interactive visualizations")
  }
  
  # Extract reconstruction data
  if(is.list(ancestral_reconstruction) && 
     (!is.null(ancestral_reconstruction$ancestral_states) || 
      !is.null(ancestral_reconstruction$states))) {
    
    # Get states based on the structure of the reconstruction object
    if(!is.null(ancestral_reconstruction$ancestral_states)) {
      states <- ancestral_reconstruction$ancestral_states
      
      # Handle different structures of ancestral_states
      if(is.data.frame(states)) {
        node_ids <- states$node_id
        node_states <- states$state
      } else if(is.vector(states) && !is.null(names(states))) {
        node_ids <- as.integer(names(states))
        node_states <- states
      } else {
        stop("Unrecognized format for ancestral_states")
      }
      
    } else if(!is.null(ancestral_reconstruction$states)) {
      states <- ancestral_reconstruction$states
      
      # Handle different structures of states
      if(is.data.frame(states)) {
        node_ids <- states$node_id
        node_states <- states$state
      } else if(is.vector(states) && !is.null(names(states))) {
        node_ids <- as.integer(names(states))
        node_states <- states
      } else {
        stop("Unrecognized format for states")
      }
    } else {
      stop("Could not find state information in ancestral reconstruction")
    }
    
    # Handle confidence intervals or probability distributions if requested
    has_uncertainty <- FALSE
    if(display_options$show_confidence) {
      if(!is.null(ancestral_reconstruction$confidence_intervals)) {
        ci_lower <- ancestral_reconstruction$confidence_intervals$lower
        ci_upper <- ancestral_reconstruction$confidence_intervals$upper
        has_uncertainty <- TRUE
      } else if(!is.null(ancestral_reconstruction$state_probabilities)) {
        state_probs <- ancestral_reconstruction$state_probabilities
        has_uncertainty <- TRUE
      }
    }
    
    # Handle events if requested
    has_events <- FALSE
    if(display_options$show_events && !is.null(ancestral_reconstruction$events)) {
      events <- ancestral_reconstruction$events
      has_events <- TRUE
    }
    
  } else {
    stop("Invalid ancestral_reconstruction object")
  }
  
  # Prepare node data for visualization
  node_data_df <- data.frame(
    node_id = node_ids,
    state = node_states,
    stringsAsFactors = FALSE
  )
  
  # Add confidence intervals if available
  if(has_uncertainty) {
    if(exists("ci_lower") && exists("ci_upper")) {
      node_data_df$ci_lower <- ci_lower
      node_data_df$ci_upper <- ci_upper
      node_data_df$uncertainty_type <- "interval"
    } else if(exists("state_probs")) {
      # For discrete states, store the probability of the most likely state
      if(is.matrix(state_probs) || is.data.frame(state_probs)) {
        node_data_df$prob <- apply(state_probs, 1, max)
      } else {
        node_data_df$prob <- NA
      }
      node_data_df$uncertainty_type <- "probability"
    }
  }
  
  # Add custom node data if provided
  if(!is.null(node_data) && is.data.frame(node_data)) {
    if("node_id" %in% colnames(node_data)) {
      node_data_df <- merge(node_data_df, node_data, by = "node_id", all.x = TRUE)
    }
  }
  
  # Create base tree plot using ggtree
  base_tree <- ggtree::ggtree(tree)
  
  # Get node coordinates from ggtree
  tree_coords <- base_tree$data
  
  # Merge coordinates with node data
  node_data_df <- merge(node_data_df, tree_coords[, c("node", "x", "y", "label")],
                      by.x = "node_id", by.y = "node", all.x = TRUE)
  
  # Create node color mapping
  if(discrete_states && max(node_states, na.rm = TRUE) - min(node_states, na.rm = TRUE) < 30) {
    # For discrete states with limited range, use categorical colors
    unique_states <- sort(unique(node_states))
    n_states <- length(unique_states)
    
    if(color_scheme == "viridis") {
      node_colors <- viridis::viridis_pal()(n_states)
    } else if(color_scheme == "rainbow") {
      node_colors <- rainbow(n_states)
    } else if(color_scheme == "divergent") {
      node_colors <- RColorBrewer::brewer.pal(min(9, n_states), "RdYlBu")
      if(n_states > 9) {
        node_colors <- colorRampPalette(node_colors)(n_states)
      }
    } else {
      # Default color scheme
      node_colors <- colorRampPalette(c("blue", "green", "red"))(n_states)
    }
    
    # Map states to colors
    node_data_df$color <- node_colors[match(node_data_df$state, unique_states)]
    
  } else {
    # For continuous states, use gradient
    state_range <- range(node_states, na.rm = TRUE)
    
    if(color_scheme == "viridis") {
      color_func <- colorRamp(viridis::viridis(100))
    } else if(color_scheme == "rainbow") {
      color_func <- colorRamp(rainbow(100))
    } else if(color_scheme == "divergent") {
      color_func <- colorRamp(RColorBrewer::brewer.pal(11, "RdYlBu"))
    } else {
      # Default gradient
      color_func <- colorRamp(c("blue", "green", "red"))
    }
    
    # Normalize states to [0,1] and get RGB colors
    norm_states <- (node_data_df$state - state_range[1]) / diff(state_range)
    rgb_colors <- color_func(norm_states) / 255
    
    # Convert RGB to hex color codes
    node_data_df$color <- rgb(rgb_colors[,1], rgb_colors[,2], rgb_colors[,3])
  }
  
  # Highlight specific nodes if requested
  if(!is.null(highlight_nodes)) {
    node_data_df$highlight <- node_data_df$node_id %in% highlight_nodes
    node_data_df$size <- ifelse(node_data_df$highlight, 10, 5)
    node_data_df$opacity <- ifelse(node_data_df$highlight, 1, 0.7)
  } else {
    node_data_df$highlight <- FALSE
    node_data_df$size <- 5
    node_data_df$opacity <- 0.7
  }
  
  # Create tooltip text
  node_data_df$tooltip <- paste0(
    "Node: ", node_data_df$node_id, "<br>",
    "State: ", round(node_data_df$state, 2)
  )
  
  # Add confidence intervals to tooltip if available
  if(has_uncertainty) {
    if("ci_lower" %in% colnames(node_data_df) && "ci_upper" %in% colnames(node_data_df)) {
      node_data_df$tooltip <- paste0(
        node_data_df$tooltip, "<br>",
        "95% CI: [", round(node_data_df$ci_lower, 2), ", ", 
        round(node_data_df$ci_upper, 2), "]"
      )
    } else if("prob" %in% colnames(node_data_df)) {
      node_data_df$tooltip <- paste0(
        node_data_df$tooltip, "<br>",
        "Probability: ", round(node_data_df$prob, 3)
      )
    }
  }
  
  # Create interactive plot with plotly
  
  # First create segments for the tree edges
  edges <- data.frame()
  for(i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    
    parent_row <- tree_coords[tree_coords$node == parent, ]
    child_row <- tree_coords[tree_coords$node == child, ]
    
    # Create horizontal segment
    edges <- rbind(edges, data.frame(
      x = parent_row$x,
      y = parent_row$y,
      xend = child_row$x,
      yend = parent_row$y,
      edge_id = i,
      stringsAsFactors = FALSE
    ))
    
    # Create vertical segment
    edges <- rbind(edges, data.frame(
      x = child_row$x,
      y = parent_row$y,
      xend = child_row$x,
      yend = child_row$y,
      edge_id = i,
      stringsAsFactors = FALSE
    ))
  }
  
  # Create basic plotly object with tree edges
  p <- plotly::plot_ly(height = height, width = width) %>%
    plotly::add_segments(
      data = edges,
      x = ~x,
      y = ~y,
      xend = ~xend,
      yend = ~yend,
      color = I("black"),
      line = list(width = 1),
      showlegend = FALSE,
      hoverinfo = "none"
    )
  
  # Add nodes with states
  p <- p %>% 
    plotly::add_markers(
      data = node_data_df,
      x = ~x,
      y = ~y,
      color = I(node_data_df$color),
      size = ~size,
      opacity = ~opacity,
      text = ~tooltip,
      hoverinfo = "text",
      marker = list(
        line = list(
          width = ifelse(node_data_df$highlight, 2, 0.5),
          color = ifelse(node_data_df$highlight, "black", "white")
        )
      )
    )
  
  # Add tip labels
  if("label" %in% colnames(tree_coords)) {
    tip_labels <- tree_coords[!is.na(tree_coords$label), ]
    p <- p %>%
      plotly::add_text(
        data = tip_labels,
        x = ~x,
        y = ~y,
        text = ~label,
        textposition = "right",
        showlegend = FALSE
      )
  }
  
  # Add node labels if requested
  if(display_options$show_node_labels) {
    p <- p %>%
      plotly::add_text(
        data = node_data_df,
        x = ~x,
        y = ~y,
        text = ~node_id,
        textposition = "bottom center",
        showlegend = FALSE
      )
  }
  
  # Add events if requested and available
  if(display_options$show_events && has_events) {
    # TO BE IMPLEMENTED: Event visualization
  }
  
  # Add color scale
  if(discrete_states && max(node_states, na.rm = TRUE) - min(node_states, na.rm = TRUE) < 30) {
    # For discrete states, add a colorbar with tick values
    unique_states <- sort(unique(node_states))
    
    p <- p %>%
      plotly::colorbar(
        title = "Chromosome Number",
        tickvals = unique_states,
        ticktext = unique_states
      )
  } else {
    # For continuous states, add a continuous colorbar
    p <- p %>%
      plotly::colorbar(
        title = "Chromosome Number"
      )
  }
  
  # Configure layout
  p <- p %>%
    plotly::layout(
      title = "Interactive Ancestral Chromosome Visualization",
      xaxis = list(
        title = "",
        showticklabels = FALSE,
        zeroline = FALSE
      ),
      yaxis = list(
        title = "",
        showticklabels = FALSE,
        zeroline = FALSE
      ),
      showlegend = FALSE
    )
  
  return(p)
}

#' Create a comparative visualization of multiple reconstruction methods
#' 
#' Generates side-by-side visualizations of chromosome states reconstructed by
#' different methods for easy comparison
#' 
#' @param tree Phylogenetic tree
#' @param reconstructions List of reconstruction results from different methods
#' @param method_names Names of the methods (labels for each visualization)
#' @param layout Layout style: "side_by_side", "grid", or "overlay"
#' @param color_scheme Color scheme for states
#' @param highlight_diffs Whether to highlight differences between methods
#' @param diff_threshold Threshold for considering states as different
#' @param node_subset Optional subset of nodes to include (NULL for all)
#' @return ggplot or list of ggplot objects
#' @export
compare_reconstruction_methods_visually <- function(tree, 
                                                 reconstructions, 
                                                 method_names = NULL,
                                                 layout = "side_by_side",
                                                 color_scheme = "viridis",
                                                 highlight_diffs = TRUE,
                                                 diff_threshold = 1,
                                                 node_subset = NULL) {
  # Check for required packages
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }
  
  # Validate inputs
  if(!is.list(reconstructions)) {
    stop("reconstructions must be a list of reconstruction results")
  }
  
  # Assign method names if not provided
  if(is.null(method_names)) {
    method_names <- paste("Method", 1:length(reconstructions))
  }
  
  # Check that method_names match reconstructions length
  if(length(method_names) != length(reconstructions)) {
    stop("Number of method names must match number of reconstructions")
  }
  
  # Extract ancestral states from each reconstruction
  all_states <- list()
  for(i in seq_along(reconstructions)) {
    recon <- reconstructions[[i]]
    
    # Get the states based on structure
    if(!is.null(recon$ancestral_states)) {
      states <- recon$ancestral_states
    } else if(!is.null(recon$states)) {
      states <- recon$states
    } else {
      stop(paste("Could not find states in reconstruction", i))
    }
    
    # Convert to named vector if necessary
    if(is.data.frame(states)) {
      node_ids <- states$node_id
      state_values <- states$state
      states <- state_values
      names(states) <- node_ids
    }
    
    all_states[[i]] <- states
  }
  
  # Determine common nodes across all reconstructions
  common_nodes <- Reduce(intersect, lapply(all_states, names))
  
  # Filter to node_subset if provided
  if(!is.null(node_subset)) {
    common_nodes <- intersect(common_nodes, as.character(node_subset))
  }
  
  # Create data frame with all states
  state_df <- data.frame(node_id = common_nodes, stringsAsFactors = FALSE)
  
  for(i in seq_along(all_states)) {
    method <- method_names[i]
    state_df[[method]] <- as.numeric(all_states[[i]][common_nodes])
  }
  
  # Calculate differences if requested
  if(highlight_diffs && length(method_names) > 1) {
    for(i in 1:(length(method_names)-1)) {
      for(j in (i+1):length(method_names)) {
        diff_col <- paste0("diff_", method_names[i], "_", method_names[j])
        state_df[[diff_col]] <- abs(state_df[[method_names[i]]] - state_df[[method_names[j]]])
        state_df[[paste0(diff_col, "_significant")]] <- state_df[[diff_col]] > diff_threshold
      }
    }
  }
  
  # Create plots based on layout
  if(layout == "side_by_side") {
    plots <- list()
    
    for(i in seq_along(method_names)) {
      method <- method_names[i]
      
      # Create tree with states for this method
      states_for_plot <- all_states[[i]]
      
      # Create ggtree visualization
      p <- create_ggtree_with_states(tree, states_for_plot, 
                                    title = method,
                                    color_scheme = color_scheme)
      
      plots[[method]] <- p
    }
    
    # If we're using patchwork or gridExtra for arranging plots
    if(requireNamespace("patchwork", quietly = TRUE)) {
      # Convert list to patchwork object
      final_plot <- plots[[1]]
      for(i in 2:length(plots)) {
        final_plot <- final_plot + plots[[i]]
      }
      
      # Adjust layout if needed
      if(length(plots) > 2) {
        n_col <- min(3, length(plots))
        final_plot <- final_plot + patchwork::plot_layout(ncol = n_col)
      }
      
    } else if(requireNamespace("gridExtra", quietly = TRUE)) {
      # Arrange with gridExtra
      n_col <- min(3, length(plots))
      final_plot <- do.call(gridExtra::grid.arrange, 
                           c(plots, list(ncol = n_col)))
    } else {
      # Return list if neither package is available
      final_plot <- plots
    }
    
    return(final_plot)
    
  } else if(layout == "overlay") {
    # Create one tree with multiple layers for different methods
    
    # Create base tree
    p <- ggtree::ggtree(tree)
    
    # Set up color palettes for each method
    palettes <- list()
    
    for(i in seq_along(method_names)) {
      method <- method_names[i]
      states <- all_states[[i]]
      
      # Determine color range
      state_range <- range(as.numeric(states), na.rm = TRUE)
      state_values <- seq(state_range[1], state_range[2], length.out = 100)
      
      # Create color palette
      if(color_scheme == "viridis") {
        palette <- viridis::viridis(100)
      } else if(color_scheme == "rainbow") {
        palette <- rainbow(100)
      } else {
        palette <- colorRampPalette(c("blue", "green", "red"))(100)
      }
      
      palettes[[method]] <- palette
      
      # TO BE IMPLEMENTED: Add layers for each method
      # This would require a custom implementation beyond standard ggtree
    }
    
    # This layout approach needs custom implementation
    stop("Overlay layout not yet implemented")
    
  } else if(layout == "grid") {
    # Create a grid comparison with differences highlighted
    
    # TO BE IMPLEMENTED: Grid layout
    stop("Grid layout not yet implemented")
    
  } else {
    stop(paste("Unsupported layout:", layout))
  }
}

#' Create a ggtree visualization with ancestral states
#' 
#' Helper function to create consistent tree visualizations
#' 
#' @param tree Phylogenetic tree
#' @param states Named vector of states with node IDs as names
#' @param title Plot title
#' @param color_scheme Color scheme to use
#' @return ggplot object
#' @keywords internal
create_ggtree_with_states <- function(tree, states, title = "", color_scheme = "viridis") {
  # Create basic tree
  p <- ggtree::ggtree(tree)
  
  # Extract tree data
  tree_data <- p$data
  
  # Create data frame with node info
  node_data <- data.frame(
    node = as.integer(names(states)),
    state = as.numeric(states),
    stringsAsFactors = FALSE
  )
  
  # Merge with tree data
  tree_data <- merge(tree_data, node_data, by = "node", all.x = TRUE)
  
  # Create node plot
  if(color_scheme == "viridis") {
    p <- ggtree::ggtree(tree) +
      ggtree::geom_nodepoint(aes(color = states), size = 3) +
      viridis::scale_color_viridis(name = "Chromosome Number", option = "D") +
      ggtree::geom_tiplab() +
      ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else if(color_scheme == "rainbow") {
    p <- ggtree::ggtree(tree) +
      ggtree::geom_nodepoint(aes(color = states), size = 3) +
      ggplot2::scale_color_gradientn(
        name = "Chromosome Number",
        colors = rainbow(100)) +
      ggtree::geom_tiplab() +
      ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    # Default color scheme
    p <- ggtree::ggtree(tree) +
      ggtree::geom_nodepoint(aes(color = states), size = 3) +
      ggplot2::scale_color_gradientn(
        name = "Chromosome Number",
        colors = colorRampPalette(c("blue", "green", "red"))(100)) +
      ggtree::geom_tiplab() +
      ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  
  return(p)
}

#' Create a comprehensive chromosome evolution figure for publications
#' 
#' Generates publication-ready figures showing chromosome evolution and key findings
#' 
#' @param tree Phylogenetic tree
#' @param reconstruction Ancestral reconstruction results
#' @param events Optional chromosome evolution events
#' @param stats Optional statistical test results
#' @param clades Optional clade definitions for highlighting
#' @param figure_type Type of figure: "standard", "annotated", or "multi_panel"
#' @param color_scheme Color scheme for visualization
#' @param output_file Optional file path to save the figure
#' @param width Width of the figure in inches
#' @param height Height of the figure in inches
#' @param dpi DPI for saved figure
#' @param show_legend Whether to show legend
#' @return ggplot object with publication-ready figure
#' @export
create_publication_figure <- function(tree, 
                                    reconstruction, 
                                    events = NULL,
                                    stats = NULL,
                                    clades = NULL,
                                    figure_type = "standard",
                                    color_scheme = "viridis",
                                    output_file = NULL,
                                    width = 8,
                                    height = 10,
                                    dpi = 300,
                                    show_legend = TRUE) {
  # Check inputs
  if(!inherits(tree, "phylo")) {
    stop("tree must be a phylo object")
  }
  
  # Extract ancestral states
  if(is.list(reconstruction)) {
    if(!is.null(reconstruction$ancestral_states)) {
      states <- reconstruction$ancestral_states
    } else if(!is.null(reconstruction$states)) {
      states <- reconstruction$states
    } else {
      stop("Could not find states in reconstruction")
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
  
  # Create plot based on figure type
  if(figure_type == "standard") {
    # Create standard chromosome evolution visualization
    p <- create_standard_figure(tree, states, color_scheme, show_legend)
    
  } else if(figure_type == "annotated") {
    # Create annotated figure with events and information
    p <- create_annotated_figure(tree, states, events, color_scheme, show_legend)
    
  } else if(figure_type == "multi_panel") {
    # Create multi-panel figure with various information
    p <- create_multi_panel_figure(tree, states, events, stats, clades, color_scheme)
    
  } else {
    stop(paste("Unsupported figure type:", figure_type))
  }
  
  # Save figure if output_file provided
  if(!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
}

#' Create a standard chromosome evolution figure
#' 
#' @param tree Phylogenetic tree
#' @param states Ancestral states
#' @param color_scheme Color scheme
#' @param show_legend Whether to show legend
#' @return ggplot object
#' @keywords internal
create_standard_figure <- function(tree, states, color_scheme, show_legend) {
  # Calculate range of states for color scale
  state_range <- range(as.numeric(states), na.rm = TRUE)
  
  # Create basic tree layout
  p <- ggtree::ggtree(tree) +
    ggtree::geom_tiplab(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Chromosome Number Evolution",
      subtitle = paste("Range:", state_range[1], "to", state_range[2])
    )
  
  # Add tree scale if branch lengths are available
  if(!is.null(tree$edge.length)) {
    p <- p + ggtree::theme_tree2()
  }
  
  # Add node points with colors representing states
  if(color_scheme == "viridis" && requireNamespace("viridis", quietly = TRUE)) {
    p <- p + 
      ggtree::geom_nodepoint(aes(color = as.numeric(states[as.character(node)])), 
                           size = 3, na.rm = TRUE) +
      viridis::scale_color_viridis(name = "Chromosome Number", option = "D")
  } else {
    p <- p + 
      ggtree::geom_nodepoint(aes(color = as.numeric(states[as.character(node)])), 
                           size = 3, na.rm = TRUE) +
      ggplot2::scale_color_gradientn(
        name = "Chromosome Number",
        colors = colorRampPalette(c("blue", "green", "red"))(100)
      )
  }
  
  # Hide legend if requested
  if(!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  return(p)
}

#' Create an annotated chromosome evolution figure with events
#' 
#' @param tree Phylogenetic tree
#' @param states Ancestral states
#' @param events Chromosome evolution events
#' @param color_scheme Color scheme
#' @param show_legend Whether to show legend
#' @return ggplot object
#' @keywords internal
create_annotated_figure <- function(tree, states, events, color_scheme, show_legend) {
  # Start with standard figure
  p <- create_standard_figure(tree, states, color_scheme, show_legend)
  
  # Add event annotations if provided
  if(!is.null(events) && is.data.frame(events) && nrow(events) > 0) {
    # Check required columns
    required_cols <- c("Edge", "Event_Type", "Description")
    missing_cols <- setdiff(required_cols, colnames(events))
    
    if(length(missing_cols) == 0) {
      # Create event symbols and colors
      event_shapes <- c(
        "wgd" = 24,    # Up-pointing triangle
        "fusion" = 21, # Circle
        "fission" = 22, # Square
        "gain" = 25,   # Up-pointing triangle
        "loss" = 23    # Diamond
      )
      
      event_colors <- c(
        "wgd" = "red",
        "fusion" = "blue",
        "fission" = "green",
        "gain" = "purple",
        "loss" = "orange"
      )
      
      # Ensure Event_Type values are lowercase for matching
      events$Event_Type_lower <- tolower(events$Event_Type)
      
      # Get event coordinates
      tree_data <- p$data
      event_coords <- data.frame()
      
      for(i in 1:nrow(events)) {
        edge_idx <- events$Edge[i]
        if(edge_idx <= nrow(tree$edge)) {
          parent <- tree$edge[edge_idx, 1]
          child <- tree$edge[edge_idx, 2]
          
          parent_row <- tree_data[tree_data$node == parent, ]
          child_row <- tree_data[tree_data$node == child, ]
          
          # Place event in middle of edge
          x <- (parent_row$x + child_row$x) / 2
          y <- (parent_row$y + child_row$y) / 2
          
          event_coords <- rbind(event_coords, data.frame(
            x = x,
            y = y,
            edge = edge_idx,
            event_type = events$Event_Type[i],
            event_type_lower = events$Event_Type_lower[i],
            description = events$Description[i],
            stringsAsFactors = FALSE
          ))
        }
      }
      
      if(nrow(event_coords) > 0) {
        # Add events to plot
        for(type in unique(event_coords$event_type_lower)) {
          if(type %in% names(event_shapes)) {
            type_data <- event_coords[event_coords$event_type_lower == type, ]
            
            p <- p + 
              ggplot2::geom_point(
                data = type_data,
                ggplot2::aes(x = x, y = y),
                shape = event_shapes[type],
                color = "black",
                fill = event_colors[type],
                size = 4
              )
          }
        }
        
        # Add a legend for event types
        legend_data <- data.frame(
          event_type = names(event_shapes),
          stringsAsFactors = FALSE
        )
        
        # TO BE IMPLEMENTED: Add event type legend
      }
    } else {
      warning(paste("Events data frame missing required columns:", paste(missing_cols, collapse = ", ")))
    }
  }
  
  return(p)
}

#' Create a multi-panel chromosome evolution figure
#' 
#' @param tree Phylogenetic tree
#' @param states Ancestral states
#' @param events Chromosome evolution events
#' @param stats Statistical test results
#' @param clades Clade definitions
#' @param color_scheme Color scheme
#' @return ggplot object (multi-panel)
#' @keywords internal
create_multi_panel_figure <- function(tree, states, events, stats, clades, color_scheme) {
  # Check for required packages
  if(!requireNamespace("patchwork", quietly = TRUE) && !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Either 'patchwork' or 'gridExtra' package is required for multi-panel figures")
  }
  
  # Create main tree plot
  main_plot <- create_annotated_figure(tree, states, events, color_scheme, TRUE)
  
  # Initialize list of panels
  panels <- list(main_plot)
  
  # Add histogram of chromosome counts if we have tip states
  if(!is.null(tree$tip.label)) {
    tip_states <- extract_tip_states(tree, states)
    
    if(length(tip_states) > 0) {
      hist_plot <- ggplot2::ggplot(data.frame(count = tip_states), ggplot2::aes(x = count)) +
        ggplot2::geom_histogram(bins = 15, fill = "steelblue", color = "black") +
        ggplot2::labs(
          title = "Chromosome Number Distribution",
          x = "Chromosome Number",
          y = "Frequency"
        ) +
        ggplot2::theme_minimal()
      
      panels <- c(panels, list(hist_plot))
    }
  }
  
  # Add model comparison panel if stats provided
  if(!is.null(stats) && is.list(stats) && !is.null(stats$model_comparison)) {
    model_comp <- stats$model_comparison
    
    if(is.data.frame(model_comp) && nrow(model_comp) > 0 && 
       all(c("Model", "AICc_weight") %in% colnames(model_comp))) {
      
      model_plot <- ggplot2::ggplot(model_comp, 
                                    ggplot2::aes(x = reorder(Model, -AICc_weight), y = AICc_weight)) +
        ggplot2::geom_bar(stat = "identity", fill = "darkgreen", color = "black") +
        ggplot2::labs(
          title = "Model Comparison",
          x = "Model",
          y = "AICc Weight"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      panels <- c(panels, list(model_plot))
    }
  }
  
  # Add clade comparison panel if clades provided
  if(!is.null(clades) && is.list(clades) && length(clades) > 1) {
    # Calculate stats for each clade
    clade_stats <- list()
    
    for(clade_name in names(clades)) {
      clade_tips <- clades[[clade_name]]
      clade_states <- extract_tip_states(tree, states, clade_tips)
      
      if(length(clade_states) > 0) {
        clade_stats[[clade_name]] <- data.frame(
          clade = clade_name,
          mean = mean(clade_states, na.rm = TRUE),
          median = median(clade_states, na.rm = TRUE),
          min = min(clade_states, na.rm = TRUE),
          max = max(clade_states, na.rm = TRUE),
          n = length(clade_states),
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Create clade comparison plot
    if(length(clade_stats) > 1) {
      clade_data <- do.call(rbind, clade_stats)
      
      clade_plot <- ggplot2::ggplot(clade_data, 
                                    ggplot2::aes(x = reorder(clade, mean), y = mean)) +
        ggplot2::geom_bar(stat = "identity", fill = "orange", color = "black") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = min, ymax = max), width = 0.2) +
        ggplot2::labs(
          title = "Chromosome Numbers by Clade",
          x = "Clade",
          y = "Mean Chromosome Number"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      
      panels <- c(panels, list(clade_plot))
    }
  }
  
  # Arrange panels
  if(length(panels) == 1) {
    # Just return the main plot
    return(panels[[1]])
  } else {
    # Determine layout
    if(length(panels) == 2) {
      # 2 panels: 1 row, 2 columns
      n_rows <- 1
      n_cols <- 2
    } else if(length(panels) == 3 || length(panels) == 4) {
      # 3-4 panels: 2 rows, 2 columns
      n_rows <- 2
      n_cols <- 2
    } else {
      # 5+ panels: 3 rows, 2 columns
      n_rows <- 3
      n_cols <- 2
    }
    
    # Create multi-panel figure
    if(requireNamespace("patchwork", quietly = TRUE)) {
      # Use patchwork for layout
      final_plot <- panels[[1]]
      for(i in 2:length(panels)) {
        final_plot <- final_plot + panels[[i]]
      }
      
      final_plot <- final_plot + patchwork::plot_layout(ncol = n_cols)
      
    } else if(requireNamespace("gridExtra", quietly = TRUE)) {
      # Use gridExtra for layout
      final_plot <- do.call(gridExtra::grid.arrange, 
                           c(panels, list(ncol = n_cols)))
    }
    
    return(final_plot)
  }
}

#' Extract chromosome states for tree tips
#' 
#' @param tree Phylogenetic tree
#' @param states Named vector of states with node IDs as names
#' @param subset_tips Optional subset of tips to include
#' @return Named vector of tip states
#' @keywords internal
extract_tip_states <- function(tree, states, subset_tips = NULL) {
  # Get tip indices
  n_tips <- length(tree$tip.label)
  tip_ids <- 1:n_tips
  
  # Filter to subset if provided
  if(!is.null(subset_tips)) {
    subset_indices <- match(subset_tips, tree$tip.label)
    subset_indices <- subset_indices[!is.na(subset_indices)]
    
    if(length(subset_indices) == 0) {
      warning("No matching tips found in subset")
      return(numeric(0))
    }
    
    tip_ids <- subset_indices
  }
  
  # Extract states for tips
  tip_states <- numeric(length(tip_ids))
  names(tip_states) <- tree$tip.label[tip_ids]
  
  for(i in seq_along(tip_ids)) {
    tip_id <- tip_ids[i]
    
    # Try to get state from different possible structures
    if(as.character(tip_id) %in% names(states)) {
      tip_states[i] <- states[as.character(tip_id)]
    } else if(tip_id <= length(states) && !is.null(tree$tip.label)) {
      # If states is an unnamed vector matching tree tips
      tip_states[i] <- states[tip_id]
    } else {
      tip_states[i] <- NA
    }
  }
  
  return(tip_states)
}

#===============================================================================
# Specialized Visualization Functions
#===============================================================================

#' Create a chromosome number heatmap across phylogeny
#' 
#' Generates a heatmap visualization of chromosome numbers across a phylogeny,
#' allowing for visualization of patterns and trends
#' 
#' @param tree Phylogenetic tree
#' @param chr_counts Named vector of chromosome counts (tip states)
#' @param ancestral_states Optional ancestral state reconstruction results
#' @param cluster_method Clustering method: "none", "phylogenetic" or "hierarchical"
#' @param color_palette Color palette for the heatmap
#' @param include_dendrogram Whether to include dendrogram
#' @param row_labels Whether to show tip labels
#' @param node_labels Whether to show node labels
#' @param label_size Size of the labels
#' @return ggplot object with the heatmap
#' @export
create_chromosome_heatmap <- function(tree,
                                    chr_counts,
                                    ancestral_states = NULL,
                                    cluster_method = "phylogenetic",
                                    color_palette = "viridis",
                                    include_dendrogram = TRUE,
                                    row_labels = TRUE,
                                    node_labels = FALSE,
                                    label_size = 3) {
  # Check required packages
  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }
  
  # Ensure chr_counts is a named vector matching tree tips
  if(is.null(names(chr_counts))) {
    stop("chr_counts must be a named vector with species names")
  }
  
  # Match chromosome counts to tree tips
  tip_states <- numeric(length(tree$tip.label))
  names(tip_states) <- tree$tip.label
  
  for(tip in tree$tip.label) {
    if(tip %in% names(chr_counts)) {
      tip_states[tip] <- chr_counts[tip]
    } else {
      tip_states[tip] <- NA
    }
  }
  
  # Combine with ancestral states if provided
  if(!is.null(ancestral_states)) {
    if(is.list(ancestral_states) && !is.null(ancestral_states$ancestral_states)) {
      states <- ancestral_states$ancestral_states
    } else {
      states <- ancestral_states
    }
    
    # Convert to named vector if necessary
    if(is.data.frame(states)) {
      node_ids <- states$node_id
      state_values <- states$state
      states <- state_values
      names(states) <- node_ids
    }
    
    # Combine tip and ancestral states
    all_states <- c(tip_states, states)
  } else {
    all_states <- tip_states
  }
  
  # Order tips based on clustering method
  if(cluster_method == "none") {
    # Use default tree order
    tip_order <- tree$tip.label
  } else if(cluster_method == "phylogenetic") {
    # Order by tree structure
    if(requireNamespace("ape", quietly = TRUE)) {
      tree_ordered <- ape::ladderize(tree)
      tip_order <- tree_ordered$tip.label
    } else {
      warning("ape package required for phylogenetic clustering. Using default order.")
      tip_order <- tree$tip.label
    }
  } else if(cluster_method == "hierarchical") {
    # Hierarchical clustering by chromosome number
    if(!anyNA(tip_states)) {
      dist_mat <- dist(tip_states)
      hc <- hclust(dist_mat)
      tip_order <- tree$tip.label[hc$order]
    } else {
      warning("Missing values in chromosome counts. Using default order.")
      tip_order <- tree$tip.label
    }
  } else {
    stop(paste("Unsupported clustering method:", cluster_method))
  }
  
  # Prepare data for heatmap
  heatmap_data <- data.frame(
    Taxon = names(tip_states),
    Chromosome_Number = tip_states,
    stringsAsFactors = FALSE
  )
  
  # Order data by tip_order
  heatmap_data$Taxon <- factor(heatmap_data$Taxon, levels = tip_order)
  
  # Create the heatmap
  if(color_palette == "viridis" && requireNamespace("viridis", quietly = TRUE)) {
    p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = 1, y = Taxon, fill = Chromosome_Number)) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis(name = "Chromosome Number", na.value = "gray90") +
      ggplot2::labs(title = "Chromosome Number Heatmap", x = "", y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                  axis.text.x = ggplot2::element_blank())
  } else {
    p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = 1, y = Taxon, fill = Chromosome_Number)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradientn(
        name = "Chromosome Number",
        colors = colorRampPalette(c("blue", "green", "red"))(100),
        na.value = "gray90"
      ) +
      ggplot2::labs(title = "Chromosome Number Heatmap", x = "", y = "") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                  axis.text.x = ggplot2::element_blank())
  }
  
  # Hide row labels if requested
  if(!row_labels) {
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank())
  } else {
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_text(size = label_size))
  }
  
  # Add dendrogram if requested
  if(include_dendrogram && cluster_method != "none") {
    # TO BE IMPLEMENTED: Add dendrogram
    # This requires additional packages and custom handling
  }
  
  return(p)
}

#' Create animation of chromosome evolution
#' 
#' Generates an animation showing the progression of chromosome number changes
#' through time along the phylogeny
#' 
#' @param tree Phylogenetic tree
#' @param ancestral_states Ancestral state reconstruction results
#' @param n_frames Number of frames in animation
#' @param color_palette Color palette for states
#' @param animate_method Animation method: "sequential" or "all_at_once"
#' @param show_tips Whether to show tip labels
#' @param width Width of animation
#' @param height Height of animation
#' @param fps Frames per second
#' @param output_file Optional file path to save animation
#' @return gganim object
#' @export
animate_chromosome_evolution <- function(tree,
                                       ancestral_states,
                                       n_frames = 30,
                                       color_palette = "viridis",
                                       animate_method = "sequential",
                                       show_tips = TRUE,
                                       width = 800,
                                       height = 600,
                                       fps = 10,
                                       output_file = NULL) {
  # Check for gganimate package
  if(!requireNamespace("gganimate", quietly = TRUE)) {
    stop("gganimate package is required for animation")
  }
  
  # Extract ancestral states
  if(is.list(ancestral_states)) {
    if(!is.null(ancestral_states$ancestral_states)) {
      states <- ancestral_states$ancestral_states
    } else if(!is.null(ancestral_states$states)) {
      states <- ancestral_states$states
    } else {
      stop("Could not find states in ancestral_states")
    }
  } else {
    states <- ancestral_states
  }
  
  # Convert to named vector if necessary
  if(is.data.frame(states)) {
    node_ids <- states$node_id
    state_values <- states$state
    states <- state_values
    names(states) <- node_ids
  }
  
  # TO BE IMPLEMENTED: Create animation based on method
  if(animate_method == "sequential") {
    # Sequential animation showing states appearing through time
    stop("Sequential animation not yet implemented")
  } else if(animate_method == "all_at_once") {
    # All states shown at once, with transitions
    stop("All-at-once animation not yet implemented")
  } else {
    stop(paste("Unsupported animation method:", animate_method))
  }
  
  # This function needs to be implemented with gganimate
  stop("Animation functionality not yet implemented")
}

#===============================================================================
# Utility Functions
#===============================================================================

#' Create a color legend for chromosome states
#' 
#' @param states Vector of state values
#' @param color_scheme Color scheme
#' @param title Legend title
#' @param discrete Whether to treat states as discrete
#' @return ggplot object with legend only
#' @export
create_chromosome_color_legend <- function(states,
                                         color_scheme = "viridis",
                                         title = "Chromosome Number",
                                         discrete = FALSE) {
  # Create a dummy plot with the desired legend
  if(discrete) {
    unique_states <- sort(unique(states))
    
    if(color_scheme == "viridis" && requireNamespace("viridis", quietly = TRUE)) {
      colors <- viridis::viridis(length(unique_states))
    } else if(color_scheme == "rainbow") {
      colors <- rainbow(length(unique_states))
    } else {
      colors <- colorRampPalette(c("blue", "green", "red"))(length(unique_states))
    }
    
    p <- ggplot2::ggplot(data.frame(x = 1:length(unique_states), y = 1, 
                                   state = factor(unique_states))) +
      ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = state)) +
      ggplot2::scale_color_manual(name = title, values = colors) +
      ggplot2::theme_void()
  } else {
    state_range <- range(states, na.rm = TRUE)
    
    if(color_scheme == "viridis" && requireNamespace("viridis", quietly = TRUE)) {
      p <- ggplot2::ggplot(data.frame(x = 1, y = 1, state = median(states, na.rm = TRUE))) +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = state)) +
        viridis::scale_color_viridis(name = title, limits = state_range) +
        ggplot2::theme_void()
    } else if(color_scheme == "rainbow") {
      p <- ggplot2::ggplot(data.frame(x = 1, y = 1, state = median(states, na.rm = TRUE))) +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = state)) +
        ggplot2::scale_color_gradientn(name = title, limits = state_range, 
                                    colors = rainbow(100)) +
        ggplot2::theme_void()
    } else {
      p <- ggplot2::ggplot(data.frame(x = 1, y = 1, state = median(states, na.rm = TRUE))) +
        ggplot2::geom_point(ggplot2::aes(x = x, y = y, color = state)) +
        ggplot2::scale_color_gradientn(name = title, limits = state_range, 
                                    colors = colorRampPalette(c("blue", "green", "red"))(100)) +
        ggplot2::theme_void()
    }
  }
  
  # Extract just the legend
  legend <- ggplot2::cowplot::get_legend(p)
  
  # Create a plot with just the legend
  legend_plot <- ggplot2::ggplot() +
    ggplot2::cowplot::draw_grob(legend) +
    ggplot2::theme_void()
  
  return(legend_plot)
}
