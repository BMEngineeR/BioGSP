#' Simulate Multiple Center Patterns
#'
#' @description Generate spatial patterns with multiple circular centers at different scales.
#' Creates concentric circle patterns with inner circle A and outer ring B at various
#' radius combinations.
#'
#' @param grid_size Size of the spatial grid (default: 60)
#' @param n_centers Number of pattern centers to generate (default: 3)
#' @param Ra_seq Vector of inner circle radii (default: c(10, 5, 1))
#' @param Rb_seq Vector of outer ring radii (default: c(10, 5, 1))
#' @param seed Random seed for reproducible center placement (default: 123)
#' @param verbose Logical; if TRUE, show progress bar and messages (default: TRUE)
#'
#' @return List of data frames, each containing X, Y coordinates and circleA, circleB binary signals
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate multi-center patterns with default parameters
#' patterns <- simulate_multiscale()
#' 
#' # Custom parameters
#' Ra_seq <- seq(from = 10, to = 3, length.out = 6)
#' Rb_seq <- seq(from = 20, to = 3, length.out = 6)
#' patterns <- simulate_multiscale(Ra_seq = Ra_seq, Rb_seq = Rb_seq, n_centers = 3)
#' }
simulate_multiscale <- function(
    grid_size = 60,
    n_centers = 3,
    Ra_seq = c(10, 5, 1),   # Ra decreases left to right
    Rb_seq = c(10, 5, 1),   # Rb decreases top to bottom
    seed = 123,
    verbose = TRUE
) {
  grid <- expand.grid(X = 1:grid_size, Y = 1:grid_size)
  
  set.seed(seed)
  if (n_centers == 1) {
    # If only one center, place it at the center of the grid
    centers <- list(c(grid_size/2, grid_size/2))
  } else {
    # For multiple centers, place them randomly
    centers <- lapply(1:n_centers, function(i) {
      c(sample(10:(grid_size-10), 1), sample(10:(grid_size-10), 1))
    })
  }
  
  sim_list <- list()
  
  # Calculate total iterations for progress bar
  total_iterations <- length(Ra_seq) * length(Rb_seq)
  if (verbose) {
    cat("Generating", total_iterations, "multiscale patterns...\n")
  }
  
  # Initialize progress bar
  pb <- NULL
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)
  }
  iteration <- 0
  
  for (Ra in Ra_seq) {
    for (Rb in Rb_seq) {
      circleA <- rep(FALSE, nrow(grid))
      circleB <- rep(FALSE, nrow(grid))
      
      for (center in centers) {
        cx <- center[1]; cy <- center[2]
        dists <- sqrt((grid$X - cx)^2 + (grid$Y - cy)^2)
        
        circleA <- circleA | (dists <= Ra)
        circleB <- circleB | (dists > Ra & dists <= (Ra + Rb))
      }
      
      circleA <- ifelse(circleA, 1, 0)
      circleB <- ifelse(!circleA & circleB, 1, 0)
      
      df <- data.frame(
        X = grid$X,
        Y = grid$Y,
        signal_1 = circleA,
        signal_2 = circleB
      )
      name <- paste0("simulated_Ra_", Ra, "_Rb_", Rb)
      sim_list[[name]] <- df
      
      # Update progress bar
      iteration <- iteration + 1
      if (verbose && !is.null(pb)) {
        setTxtProgressBar(pb, iteration)
      }
    }
  }
  
  # Close progress bar
  if (verbose && !is.null(pb)) {
    close(pb)
    cat("\nMultiscale simulation completed!\n")
  }
  
  return(sim_list)
}

#' Simulate Stripe Patterns
#'
#' @description Generate stripe patterns with two parallel stripes separated by a gap.
#' Creates rotatable stripe patterns with configurable gap, width, and rotation angle.
#'
#' @param grid_size Size of the spatial grid (default: 100)
#' @param gap_seq Vector of gap distances between stripe centers (default: c(10))
#' @param width_seq Vector of stripe widths (default: c(5))
#' @param theta_seq Vector of rotation angles in degrees (default: c(0))
#' @param eps Small numeric value for open boundary conditions to avoid overlap at stripe edges (default: 1e-9)
#' @param verbose Logical; if TRUE, show progress messages (default: TRUE)
#'
#' @return List of data frames, each containing X, Y coordinates and signal_1, signal_2 binary signals
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate stripe patterns with default parameters
#' patterns <- simulate_stripe_patterns()
#' 
#' # Custom parameters
#' patterns <- simulate_stripe_patterns(
#'   grid_size = 80,
#'   gap_seq = c(10, 20),
#'   width_seq = c(5, 10, 20),
#'   theta_seq = c(0, 30, 60),
#'   eps = 1e-9,
#'   verbose = TRUE
#' )
#' }
simulate_stripe_patterns <- function(
    grid_size = 100,
    gap_seq   = c(10),
    width_seq = c(5),
    theta_seq = c(0),
    eps = 1e-9,
    verbose = TRUE
) {
  # Generate lattice grid
  grid <- expand.grid(X = 1:grid_size, Y = 1:grid_size)
  sim_list <- list()
  
  total_iterations <- length(gap_seq) * length(width_seq) * length(theta_seq)
  if (verbose) cat("Generating", total_iterations, "stripe patterns...\n")
  iteration <- 0
  
  for (gap in gap_seq) {
    for (width in width_seq) {
      for (theta in theta_seq) {
        
        # Baseline vertical lines at center ± gap/2
        x_left_line  <- grid_size/2 - gap/2
        x_right_line <- grid_size/2 + gap/2
        y_center     <- grid_size/2
        
        # Rotation function
        rotate <- function(x, y, theta, cx, cy) {
          th <- theta * pi/180
          x_shift <- as.numeric(x - cx)
          y_shift <- as.numeric(y - cy)
          x_rot <-  cos(th) * x_shift + sin(th) * y_shift
          y_rot <- -sin(th) * x_shift + cos(th) * y_shift
          data.frame(x = x_rot, y = y_rot)
        }
        
        if (gap == 0) {
          # Both stripes meet in the middle
          coords <- rotate(grid$X, grid$Y, theta, grid_size/2, y_center)
          signal_1 <- as.integer(coords$x >= -width & coords$x <= 0)
          signal_2 <- as.integer(coords$x > 0    & coords$x <= width)
        } else {
          # Normal case: two separated stripes
          coords_left  <- rotate(grid$X, grid$Y, theta, x_left_line, y_center)
          coords_right <- rotate(grid$X, grid$Y, theta, x_right_line, y_center)
          
          signal_1 <- as.integer(coords_left$x <= 0 & coords_left$x >= -width)
          signal_2 <- as.integer(coords_right$x >= 0 & coords_right$x <= width)
        }
        
        df <- data.frame(
          X = grid$X,
          Y = grid$Y,
          signal_1 = signal_1,
          signal_2 = signal_2,
          gap_param = gap,
          width_param = width,
          theta_param = theta
        )
        
        name <- paste0("gap_", gap, "_w_", width, "_th_", theta)
        sim_list[[name]] <- df
        
        iteration <- iteration + 1
        if (verbose && iteration %% max(1, floor(total_iterations/10)) == 0) {
          cat("Completed", iteration, "of", total_iterations, "patterns\n")
        }
      }
    }
  }
  
  if (verbose) cat("Stripe pattern simulation completed!\n")
  return(sim_list)
}

#' Visualize Multiple Center Simulation Results
#'
#' @description Create visualization plots for multiple center simulation patterns
#'
#' @param sim_data Output from simulate_multiscale function
#' @param Ra_seq Vector of Ra values used in simulation
#' @param Rb_seq Vector of Rb values used in simulation
#' @param bg_color Background color for plots (default: "grey")
#' @param signal1_color Color for signal 1 (default: "red")
#' @param signal2_color Color for signal 2 (default: "blue")
#' @param show_title Logical; if TRUE (default), add titles to plots with Ra and Rb values
#'
#' @return Combined ggplot object with all pattern visualizations
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate and visualize patterns
#' Ra_seq <- seq(from = 10, to = 3, length.out = 6)
#' Rb_seq <- seq(from = 20, to = 3, length.out = 6)
#' sim_data <- simulate_multiscale(Ra_seq = Ra_seq, Rb_seq = Rb_seq, n_centers = 3)
#' plot_grid <- visualize_multiscale(sim_data, Ra_seq, Rb_seq)
#' print(plot_grid)
#' }
visualize_multiscale <- function(sim_data, Ra_seq, Rb_seq, 
                                bg_color = "grey",      
                                signal1_color = "red",    
                                signal2_color = "blue",
                                show_title = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package is required for plot arrangement")
  }
  
  plot_list <- list()
  for (Rb in Rb_seq) {
    for (Ra in Ra_seq) {
      df <- sim_data[[paste0("simulated_Ra_", Ra, "_Rb_", Rb)]]
      df$label <- "Background"
      df$label[df$signal_1 == 1] <- "signal_1"
      df$label[df$signal_2 == 1] <- "signal_2"
      # print(table(df$label))
      p <- ggplot2::ggplot(df, ggplot2::aes_string(x = "X", y = "Y", fill = "label")) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_manual(values = c("Background" = bg_color,
                                        "signal_1" = signal1_color,
                                        "signal_2" = signal2_color)) +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "none")
      
      if (show_title) {
        p <- p + ggplot2::ggtitle(paste0("Ra=", round(Ra, 1), ", Rb=", round(Rb, 1)))
      }
      
      plot_list[[paste0("Ra_", Ra, "_Rb_", Rb)]] <- p
    }
  }
  
  # Arrange plots in grid (Rb = rows, Ra = cols)
  plot_grid <- patchwork::wrap_plots(plot_list, ncol = length(Ra_seq))
  return(plot_grid)
}

#' Simulate Multiple Center Patterns with Fixed Centers
#'
#' @description Generate spatial patterns with multiple circular centers at fixed positions.
#' Similar to simulate_multiscale but with centers placed at fixed locations for reproducible
#' pattern generation. Creates concentric circle patterns with inner circle A and outer ring B
#' at various radius combinations.
#'
#' @param grid_size Size of the spatial grid (default: 60)
#' @param n_centers Number of pattern centers to generate. If 1, center is placed at grid center.
#'   If > 1, centers are randomly placed but fixed by seed (default: 3)
#' @param Ra_seq Vector of inner circle radii (default: c(10, 5, 1))
#' @param Rb_seq Vector of outer ring radii (default: c(10, 5, 1))
#' @param seed Random seed for reproducible center placement (default: 123)
#' @param verbose Logical; if TRUE, show progress bar and messages (default: TRUE)
#'
#' @return List of data frames, each containing X, Y coordinates and signal_1, signal_2 binary signals
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate multi-center patterns with fixed centers
#' patterns <- simulate_multiscale_overlap()
#' 
#' # Single center at grid center
#' patterns_single <- simulate_multiscale_overlap(n_centers = 1)
#' 
#' # Custom parameters with multiple centers
#' Ra_seq <- seq(from = 10, to = 3, length.out = 4)
#' Rb_seq <- seq(from = 15, to = 2, length.out = 4)
#' patterns <- simulate_multiscale_overlap(
#'   Ra_seq = Ra_seq, 
#'   Rb_seq = Rb_seq, 
#'   n_centers = 2,
#'   seed = 456
#' )
#' }
simulate_multiscale_overlap <- function(
    grid_size = 60,
    n_centers = 3,
    Ra_seq = c(10, 5, 1),   # inner radius sequence
    Rb_seq = c(10, 5, 1),   # outer ring thickness sequence
    seed = 123,
    verbose = TRUE
) {
  # Create coordinate grid
  grid <- expand.grid(X = 1:grid_size, Y = 1:grid_size)
  
  set.seed(seed)
  
  # Define circle centers
  if (n_centers == 1) {
    # If only one center, fix at center of grid
    centers <- list(c(grid_size / 2, grid_size / 2))
  } else {
    # Otherwise randomly place multiple centers
    centers <- lapply(1:n_centers, function(i) {
      c(sample(10:(grid_size - 10), 1), sample(10:(grid_size - 10), 1))
    })
  }
  
  sim_list <- list()
  total_iterations <- length(Ra_seq) * length(Rb_seq)
  
  if (verbose) cat("Generating", total_iterations, "multiscale patterns...\n")
  
  pb <- if (verbose) txtProgressBar(min = 0, max = total_iterations, style = 3) else NULL
  iteration <- 0
  
  for (Ra in Ra_seq) {
    for (Rb in Rb_seq) {
      circleA <- rep(FALSE, nrow(grid))
      circleB <- rep(FALSE, nrow(grid))
      
      for (center in centers) {
        cx <- center[1]; cy <- center[2]
        dists <- sqrt((grid$X - cx)^2 + (grid$Y - cy)^2)
        
        # Inner circle (signal_1)
        circleA <- circleA | (dists <= Ra)
        # Outer ring (signal_2)
        circleB <- circleB | (dists > Ra & dists <= (Ra + Rb))
      }
      
      # Convert logicals to numeric signals
      signal_1 <- as.numeric(circleA)
      signal_2 <- as.numeric(circleB & !circleA)
      
      df <- data.frame(
        X = grid$X,
        Y = grid$Y,
        signal_1 = signal_1,
        signal_2 = signal_2
      )
      
      name <- paste0("simulated_Ra_", Ra, "_Rb_", Rb)
      sim_list[[name]] <- df
      
      iteration <- iteration + 1
      if (verbose && !is.null(pb)) setTxtProgressBar(pb, iteration)
    }
  }
  
  if (verbose && !is.null(pb)) {
    close(pb)
    cat("\nMultiscale (fixed-center) simulation completed!\n")
  }
  
  return(sim_list)
}


#' Visualize Stripe Pattern Simulation Results
#'
#' @description Create visualization plots for stripe pattern simulation results
#'
#' @param sim_data Output from simulate_stripe_patterns function
#' @param gap_seq Vector of gap values used in simulation
#' @param width_seq Vector of width values used in simulation
#' @param theta_seq Vector of theta (rotation angle) values used in simulation
#' @param bg_color Background color for plots (default: "grey")
#' @param signal1_color Color for signal 1 (default: "#1f6f8b")
#' @param signal2_color Color for signal 2 (default: "#e67e22")
#' @param overlap_color Color for overlapping regions (default: "#7a4dbf")
#' @param show_title Logical; if TRUE (default), add titles to plots with parameter values
#'
#' @return Combined ggplot object with all pattern visualizations
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate and visualize patterns
#' sim_data <- simulate_stripe_patterns(
#'   grid_size = 80,
#'   gap_seq = c(10, 20),
#'   width_seq = c(5, 10, 20),
#'   theta_seq = c(0, 30, 60)
#' )
#' plot_grid <- visualize_stripe_patterns(sim_data, 
#'                                        gap_seq = c(10, 20),
#'                                        width_seq = c(5, 10, 20),
#'                                        theta_seq = c(0, 30, 60))
#' print(plot_grid)
#' }
visualize_stripe_patterns <- function(
    sim_data, gap_seq, width_seq, theta_seq,
    bg_color = "grey",
    signal1_color = "#1f6f8b",
    signal2_color = "#e67e22",
    overlap_color = "#7a4dbf",
    show_title = TRUE
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package is required for plot arrangement")
  }
  
  plot_list <- list()
  
  for (gap in gap_seq) {
    for (width in width_seq) {
      for (theta in theta_seq) {
        df <- sim_data[[paste0("gap_", gap,
                               "_w_", width,
                               "_th_", theta)]]
        
        # Label background / stripes / overlap
        df$label <- "Background"
        df$label[df$signal_1 == 1 & df$signal_2 == 0] <- "signal_1"
        df$label[df$signal_2 == 1 & df$signal_1 == 0] <- "signal_2"
        df$label[df$signal_1 == 1 & df$signal_2 == 1] <- "overlap"
        
        # Raster-based plot (no seams)
        p <- ggplot2::ggplot(df, ggplot2::aes(x = X, y = Y, fill = label)) +
          ggplot2::geom_raster(interpolate = FALSE) +
          ggplot2::scale_fill_manual(values = c(
            "Background" = bg_color,
            "signal_1"   = signal1_color,
            "signal_2"   = signal2_color,
            "overlap"    = overlap_color
          )) +
          ggplot2::coord_fixed(expand = FALSE) +  # important for no gaps
          ggplot2::theme_void() +
          ggplot2::theme(legend.position = "none")
        
        if (show_title) {
          p <- p + ggplot2::ggtitle(paste0("g=", gap, ", w=", width, ", th=", theta))
        }
        
        plot_list[[paste0("gap_", gap, "_w_", width, "_th_", theta)]] <- p
      }
    }
  }
  
  patchwork::wrap_plots(plot_list, ncol = length(width_seq))
}
