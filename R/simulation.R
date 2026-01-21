#' Simulate Multi-center Multi-scale Concentric Ring Patterns
#'
#' @description Generate multi-center, multi-scale concentric ring simulation data.
#' Creates patterns with inner circles and outer rings where the outer radius shrinks
#' from a fixed starting point to a factor of the inner radius across multiple steps.
#'
#' @param grid_size Size of the spatial grid (default: 60)
#' @param Ra_seq Vector of inner circle radii (default: seq(2.5, 20, by = 2.5))
#' @param n_steps Number of outer radius shrinkage steps (default: 10)
#' @param n_centers Number of circle centers (default: 1)
#' @param outer_start Fixed starting outer radius (default: 40)
#' @param outer_end_factor Outer radius shrinks to this factor * Ra (default: 1.2)
#' @param seed Random seed for reproducible center placement (default: 123)
#' @param verbose Logical; if TRUE, show progress bar and messages (default: TRUE)
#'
#' @return List of data frames, each containing X, Y coordinates and signal_1, signal_2 binary signals
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate multi-center patterns with default parameters
#' patterns <- simulate_multiscale()
#' 
#' # Custom parameters
#' patterns <- simulate_multiscale(
#'   grid_size = 80,
#'   Ra_seq = seq(5, 25, by = 5),
#'   n_steps = 8,
#'   n_centers = 2,
#'   outer_start = 50
#' )
#' }
simulate_multiscale <- function(
    grid_size = 60,
    Ra_seq = seq(2.5, 20, by = 2.5), # inner radii sequence
    n_steps = 10,                    # number of outer radius shrinkage steps
    n_centers = 1,                   # number of circle centers
    outer_start = 40,                # fixed starting outer radius
    outer_end_factor = 1.2,          # outer radius shrinks to 1.2 * Ra
    seed = 123,
    verbose = TRUE
) {
  set.seed(seed)
  grid <- expand.grid(X = 1:grid_size, Y = 1:grid_size)
  
  # ---- Define circle centers ----
  if (n_centers == 1) {
    centers <- list(c(grid_size/2, grid_size/2))
  } else {
    centers <- lapply(1:n_centers, function(i) {
      c(sample(10:(grid_size - 10), 1), sample(10:(grid_size - 10), 1))
    })
  }
  
  sim_list <- list()
  total_iter <- length(Ra_seq) * n_steps
  if (verbose) cat("Generating", total_iter, "patterns...\n")
  pb <- if (verbose) txtProgressBar(min = 0, max = total_iter, style = 3) else NULL
  iter <- 0
  
  for (Ra in Ra_seq) {
    # Outer radius shrinks from fixed outer_start to Ra * outer_end_factor
    Rb_seq <- seq(from = outer_start,
                  to   = Ra * outer_end_factor,
                  length.out = n_steps)
    
    for (step_i in seq_along(Rb_seq)) {
      Rb <- Rb_seq[step_i]
      circleA <- rep(FALSE, nrow(grid))
      circleB <- rep(FALSE, nrow(grid))
      
      for (center in centers) {
        cx <- center[1]; cy <- center[2]
        dists <- sqrt((grid$X - cx)^2 + (grid$Y - cy)^2)
        circleA <- circleA | (dists <= Ra)
        circleB <- circleB | (dists > Ra & dists <= Rb)
      }
      
      df <- data.frame(
        X = grid$X,
        Y = grid$Y,
        signal_1 = as.numeric(circleA),
        signal_2 = as.numeric(circleB & !circleA)
      )
      
      sim_list[[paste0("Ra_", Ra, "_Step_", step_i)]] <- df
      
      iter <- iter + 1
      if (verbose && !is.null(pb)) setTxtProgressBar(pb, iter)
    }
  }
  
  if (verbose && !is.null(pb)) {
    close(pb)
    cat("\nMulticenter multiscale simulation completed!\n")
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

#' Visualize Multi-center Multi-scale Concentric Ring Patterns
#'
#' @description Visualize the simulated concentric ring patterns from simulate_multiscale
#'
#' @param sim_data Output from simulate_multiscale function
#' @param Ra_seq Vector of Ra values used in simulation
#' @param n_steps Number of steps used in simulation
#' @param bg_color Background color for plots (default: "grey90")
#' @param signal1_color Color for signal 1 (default: "#16964a")
#' @param signal2_color Color for signal 2 (default: "#2958a8")
#' @param show_subtitle Logical; if TRUE (default), show parameter values in facet labels
#' @param sort_order Order for sorting ("ascending" or "descending", default: "ascending")
#' @param panel_spacing Control spacing between panels in lines (default: 0.1)
#' @param title_size Size of title text (default: 12)
#'
#' @return ggplot object with faceted visualization
#' @importFrom dplyr bind_rows mutate case_when
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual coord_fixed facet_grid vars theme_void theme element_text element_blank unit
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate and visualize patterns
#' sim_data <- simulate_multiscale(
#'   Ra_seq = seq(2.5, 20, by = 2.5),
#'   n_steps = 10
#' )
#' plot_grid <- visualize_multiscale(sim_data, 
#'                                  Ra_seq = seq(2.5, 20, by = 2.5), 
#'                                  n_steps = 10)
#' print(plot_grid)
#' }
visualize_multiscale <- function(
    sim_data,
    Ra_seq,
    n_steps,
    bg_color = "grey90",
    signal1_color = "#16964a",
    signal2_color = "#2958a8",
    show_subtitle = TRUE,
    sort_order = c("ascending", "descending"),
    panel_spacing = 0.1,
    title_size = 12
) {
  sort_order <- match.arg(sort_order)
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for data manipulation")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }
  
  # ---- Merge all data and parse Ra / Step ----
  df_all <- dplyr::bind_rows(sim_data, .id = "name")
  df_all <- dplyr::mutate(df_all,
    Ra = as.numeric(sub(".*Ra_([0-9.]+)_Step_.*", "\\1", name)),
    Step = as.numeric(sub(".*_Step_([0-9]+)$", "\\1", name)),
    label = dplyr::case_when(
      signal_1 == 1 ~ "signal_1",
      signal_2 == 1 ~ "signal_2",
      TRUE ~ "Background"
    )
  )
  
  # ---- Sorting control ----
  Ra_vals <- sort(unique(Ra_seq), decreasing = (sort_order == "descending"))
  Step_vals <- sort(unique(df_all$Step), decreasing = (sort_order == "ascending"))
  
  df_all$Ra <- factor(df_all$Ra, levels = Ra_vals)
  df_all$Step <- factor(df_all$Step, levels = Step_vals)
  
  facet_labeller <- if (show_subtitle) {
    ggplot2::label_both
  } else {
    ggplot2::label_value   # still generate strips but without showing parameter values
  }
  
  # ---- Plot ----
  ggplot2::ggplot(df_all, ggplot2::aes(X, Y, fill = label)) +
    ggplot2::geom_tile(alpha = 0.95) +
    ggplot2::scale_fill_manual(values = c(
      "Background" = bg_color,
      "signal_1" = signal1_color,
      "signal_2" = signal2_color
    )) +
    ggplot2::coord_fixed() +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Ra),
      cols = ggplot2::vars(Step),
      labeller = facet_labeller,
      switch = "y"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      strip.text.x = if (show_subtitle) ggplot2::element_text(size = title_size - 2, face = "bold") else ggplot2::element_blank(),
      strip.text.y = if (show_subtitle) ggplot2::element_text(size = title_size - 2, face = "bold") else ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(panel_spacing, "lines")
    )
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

#' Simulate Moving Circles Pattern
#'
#' @description Generate patterns of two circles moving toward each other horizontally.
#' Creates mutually exclusive signals where overlapping pixels are assigned to signal_1 (circle 1).
#' The circles start at fixed horizontal distances from the midline and move toward the center.
#'
#' @param grid_size Size of the spatial grid (default: 60)
#' @param radius_seq Vector of radii for circle 1 (default: 6:14)
#' @param n_steps Number of movement steps (default: 10)
#' @param center_distance Initial horizontal distance from midline for both centers (default: 30)
#' @param radius2_factor Circle 2 radius = radius_seq * radius2_factor (default: 1.5)
#' @param seed Random seed for reproducibility (default: 123)
#' @param verbose Logical; if TRUE, show progress bar and messages (default: TRUE)
#'
#' @return List of data frames, each containing X, Y coordinates and signal_1, signal_2 binary signals
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate moving circles patterns with default parameters
#' patterns <- simulate_moving_circles()
#' 
#' # Custom parameters
#' patterns <- simulate_moving_circles(
#'   grid_size = 80,
#'   radius_seq = c(8, 12, 16),
#'   n_steps = 8,
#'   center_distance = 35,
#'   radius2_factor = 1.2
#' )
#' }
simulate_moving_circles <- function(
    grid_size       = 60,
    radius_seq      = 6:14,   # radii for circle 1
    n_steps         = 10,     # "Movement 1..n_steps"
    center_distance = 30,     # initial horizontal distance from midline for both centers
    radius2_factor  = 1.5,    # circle 2 radius = radius_seq * radius2_factor
    seed            = 123,
    verbose         = TRUE
) {
  stopifnot(length(radius_seq) >= 1, n_steps >= 1, grid_size >= 3)
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install 'dplyr'.")
  
  set.seed(seed)
  grid <- expand.grid(X = 1:grid_size, Y = 1:grid_size)
  
  # Deterministic movement: two circles move toward the midline horizontally
  movements <- data.frame(
    Step      = seq_len(n_steps),
    circle1_x = center_distance - seq(from = center_distance/2, to = 0, length.out = n_steps),
    circle2_x = center_distance + seq(from = center_distance/2, to = 0, length.out = n_steps),
    circle_y  = rep(grid_size / 2, n_steps)
  )
  
  sim_list <- list()
  total_iter <- length(radius_seq) * n_steps
  if (verbose) cat("Generating", total_iter, "patterns...\n")
  pb <- if (verbose) txtProgressBar(min = 0, max = total_iter, style = 3) else NULL
  iter <- 0
  
  for (r1 in radius_seq) {
    r2 <- r1 * radius2_factor
    for (k in seq_len(n_steps)) {
      cx1 <- movements$circle1_x[k]
      cx2 <- movements$circle2_x[k]
      cy  <- movements$circle_y[k]
      
      d1 <- sqrt((grid$X - cx1)^2 + (grid$Y - cy)^2)
      d2 <- sqrt((grid$X - cx2)^2 + (grid$Y - cy)^2)
      
      inside1 <- d1 <= r1
      inside2 <- d2 <= r2
      
      # ---- MUTUALLY EXCLUSIVE ASSIGNMENT ----
      # Priority rule: if a pixel falls in both, assign to circle 1 (signal_1)
      signal_1 <- as.integer(inside1)                   # 1 if in circle 1 (includes any overlap)
      signal_2 <- as.integer(inside2 & !inside1)        # 1 only if in circle 2 and NOT in circle 1
      
      df <- data.frame(
        X        = grid$X,
        Y        = grid$Y,
        signal_1 = signal_1,
        signal_2 = signal_2,
        Ra       = r1,
        Step     = k
      )
      
      sim_list[[paste0("Ra_", r1, "_Step_", k)]] <- df
      
      iter <- iter + 1
      if (!is.null(pb)) setTxtProgressBar(pb, iter)
    }
  }
  
  if (!is.null(pb)) close(pb)
  if (verbose) cat("\nMoving-circles simulation completed!\n")
  sim_list
}

#' Visualize Moving Circles Pattern
#'
#' @description Visualize the simulated moving circles patterns from simulate_moving_circles
#'
#' @param sim_data Output from simulate_moving_circles function
#' @param bg_color Background color for plots (default: "grey90")
#' @param signal1_color Color for signal 1 (default: "#16964a")
#' @param signal2_color Color for signal 2 (default: "#2958a8")
#' @param show_subtitle Logical; if TRUE (default), show parameter values in facet labels
#' @param sort_order Order for sorting ("ascending" or "descending", default: "ascending")
#' @param panel_spacing Control spacing between panels in lines (default: 0.1)
#' @param title_size Size of title text (default: 12)
#'
#' @return ggplot object with faceted visualization
#' @importFrom dplyr bind_rows mutate case_when
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_manual coord_fixed facet_grid vars theme_void theme element_text element_blank unit
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate and visualize patterns
#' sim_data <- simulate_moving_circles(
#'   radius_seq = 6:14,
#'   n_steps = 10
#' )
#' plot_grid <- visualize_moving_circles(sim_data)
#' print(plot_grid)
#' }
visualize_moving_circles <- function(
    sim_data,
    bg_color       = "grey90",
    signal1_color  = "#16964a",
    signal2_color  = "#2958a8",
    show_subtitle  = TRUE,
    sort_order     = c("ascending", "descending"),
    panel_spacing  = 0.1,
    title_size     = 12
) {
  sort_order <- match.arg(sort_order)
  if (!requireNamespace("dplyr", quietly = TRUE))   stop("Please install 'dplyr'.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  
  df_all <- dplyr::bind_rows(sim_data, .id = "name")
  
  # Parse Ra and Step from names
  df_all <- dplyr::mutate(df_all,
    Ra = as.numeric(gsub(".*Ra_([0-9.]+)_Step.*", "\\1", name)),
    Step = as.numeric(gsub(".*Step_([0-9]+)$", "\\1", name))
  )
  
  # Map to three mutually exclusive classes: signal_1, signal_2, Background
  df_all$label <- dplyr::case_when(
    df_all$signal_1 == 1 ~ "signal_1",
    df_all$signal_2 == 1 ~ "signal_2",
    TRUE                 ~ "Background"
  )
  
  fill_vals <- c(
    "Background" = bg_color,
    "signal_1"   = signal1_color,
    "signal_2"   = signal2_color
  )
  
  # ordering
  Ra_vals   <- sort(unique(df_all$Ra),  decreasing = (sort_order == "descending"))
  Step_vals <- sort(unique(df_all$Step), decreasing = FALSE)
  df_all$Ra   <- factor(df_all$Ra,   levels = Ra_vals)
  df_all$Step <- factor(df_all$Step, levels = Step_vals)
  
  facet_labeller <- if (show_subtitle) ggplot2::label_both else ggplot2::label_value
  
  ggplot2::ggplot(df_all, ggplot2::aes(X, Y, fill = label)) +
    ggplot2::geom_tile(alpha = 0.95) +
    ggplot2::scale_fill_manual(values = fill_vals) +
    ggplot2::coord_fixed() +
    ggplot2::facet_grid(
      rows     = ggplot2::vars(Ra),
      cols     = ggplot2::vars(Step),
      labeller = facet_labeller,
      switch   = "y"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      strip.text.x = if (show_subtitle) ggplot2::element_text(size = title_size - 2, face = "bold") else ggplot2::element_blank(),
      strip.text.y = if (show_subtitle) ggplot2::element_text(size = title_size - 2, face = "bold") else ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(panel_spacing, "lines")
    )
}
