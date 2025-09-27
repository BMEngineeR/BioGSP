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
#'
#' @return List of data frames, each containing X, Y coordinates and circleA, circleB binary signals
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
    seed = 123
) {
  grid <- expand.grid(X = 1:grid_size, Y = 1:grid_size)
  
  set.seed(seed)
  centers <- lapply(1:n_centers, function(i) {
    c(sample(10:(grid_size-10), 1), sample(10:(grid_size-10), 1))
  })
  
  sim_list <- list()
  
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
    }
  }
  
  return(sim_list)
}

#' Simulate Concentric Ring Patterns
#'
#' @description Generate concentric ring patterns with dynamic outer ring movement.
#' Creates a solid inner circle with a moving outer ring that closes in over time.
#'
#' @param grid_size Size of the spatial grid (default: 60)
#' @param radius_seq Vector of inner circle radii to simulate (default: seq(2.5, 20, by = 2.5))
#' @param n_movements Number of movement steps for the outer ring (default: 10)
#' @param center_x X coordinate of pattern center (default: grid_size/2)
#' @param center_y Y coordinate of pattern center (default: grid_size/2)
#'
#' @return List of data frames, each containing X, Y coordinates, movement indicators, 
#'         and signal_1 (solid circle), signal_2 (concentric ring) binary signals
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate concentric ring patterns with default parameters
#' patterns <- simulate_ringpattern()
#' 
#' # Custom parameters
#' radius_seq <- seq(2.5, 20, by = 2.5)
#' patterns <- simulate_ringpattern(radius_seq = radius_seq, n_movements = 5)
#' }
simulate_ringpattern <- function(
    grid_size = 60,
    radius_seq = seq(2.5, 20, by = 2.5),
    n_movements = 10,
    center_x = NULL,
    center_y = NULL
) {
  if (is.null(center_x)) center_x <- grid_size/2
  if (is.null(center_y)) center_y <- grid_size/2
  
  grid <- expand.grid(X = 1:grid_size, Y = 1:grid_size)
  sim_list <- list()
  
  for (r in radius_seq) {
    # movement trajectory (outer ring closing in)
    movements <- data.frame(
      r_outer = seq(from = 40, to = 1.2 * r, length.out = n_movements)
    )
    
    plot_data <- lapply(seq_len(n_movements), function(i) {
      r_out <- movements$r_outer[i]
      
      # Calculate distances from center
      distances <- sqrt((grid$X - center_x)^2 + (grid$Y - center_y)^2)
      
      # Define regions
      inside_circle1 <- distances <= r
      inside_circle2_outer <- distances <= r_out
      inside_ring <- inside_circle2_outer & !inside_circle1
      
      data.frame(
        X = grid$X,
        Y = grid$Y,
        signal_1 = as.numeric(inside_circle1),  # Solid Circle
        signal_2 = as.numeric(inside_ring),     # Concentric Ring
        movement = paste("Movement", i),
        movement_step = i,
        radius = r,
        outer_radius = r_out
      )
    })
    
    # Combine all movements for this radius
    combined_data <- do.call(rbind, plot_data)
    sim_list[[paste0("radius_", r)]] <- combined_data
  }
  
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
                                signal2_color = "blue") {
  
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
      print(table(df$label))
      p <- ggplot2::ggplot(df, ggplot2::aes_string(x = "X", y = "Y", fill = "label")) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_manual(values = c("Background" = bg_color,
                                        "signal_1" = signal1_color,
                                        "signal_2" = signal2_color)) +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ggtitle(paste0("Ra=", round(Ra, 1), ", Rb=", round(Rb, 1)))
      
      plot_list[[paste0("Ra_", Ra, "_Rb_", Rb)]] <- p
    }
  }
  
  # Arrange plots in grid (Rb = rows, Ra = cols)
  plot_grid <- patchwork::wrap_plots(plot_list, ncol = length(Ra_seq))
  return(plot_grid)
}

#' Visualize Concentric Ring Simulation Results
#'
#' @description Create visualization plots for concentric ring simulation patterns
#'
#' @param sim_data Output from simulate_ringpattern function
#' @param radius_seq Vector of radius values used in simulation
#' @param bg_color Background color for plots (default: "grey")
#' @param signal1_color Color for signal 1 (default: "#16964a")
#' @param signal2_color Color for signal 2 (default: "#2958a8")
#'
#' @return Combined ggplot object with all pattern visualizations
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate and visualize patterns
#' radius_seq <- seq(2.5, 20, by = 2.5)
#' sim_data <- simulate_ringpattern(radius_seq = radius_seq)
#' plot_grid <- visualize_ringpattern(sim_data, radius_seq)
#' print(plot_grid)
#' }
visualize_ringpattern <- function(sim_data, radius_seq, 
                                 bg_color = "grey",      
                                 signal1_color = "#16964a",    
                                 signal2_color = "#2958a8") {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package is required for plot arrangement")
  }
  
  plot_list <- list()
  for (r in radius_seq) {
    df <- sim_data[[paste0("radius_", r)]]
    
    # Create label column
    df$label <- "Outside"
    df$label[df$signal_1 == 1] <- "signal_1"
    df$label[df$signal_2 == 1] <- "signal_2"
    
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = "X", y = "Y", fill = "label")) +
      ggplot2::geom_tile(alpha = 0.95) +
      ggplot2::facet_wrap(~ movement, nrow = 1) +
      ggplot2::scale_fill_manual(values = c(
        "signal_1" =signal1_color,
        "signal_2" = signal2_color,
        "Outside" = bg_color
      )) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(paste("Radius =", r))
    
    plot_list[[paste0("radius_", r)]] <- p
  }
  
  # Arrange plots: each row = one radius
  plot_grid <- patchwork::wrap_plots(plot_list, ncol = 1)
  return(plot_grid)
}
