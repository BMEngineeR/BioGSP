#' Plot SGWT decomposition results
#'
#' @description Visualize SGWT decomposition components including original signal,
#' scaling function, wavelet coefficients, and reconstructed signal
#'
#' @param SG SGWT object with Forward and Inverse results computed
#' @param signal_name Name of signal to plot (default: first signal)
#' @param plot_scales Which wavelet scales to plot (default: first 4)
#' @param ncol Number of columns in the plot layout (default: 3)
#'
#' @return ggplot object with combined plots
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have SGWT object
#' plots <- plot_sgwt_decomposition(SG_object, signal_name = "signal1")
#' print(plots)
#' }
plot_sgwt_decomposition <- function(SG, signal_name = NULL, plot_scales = NULL, ncol = 3) {
  
  # Validate input
  if (!inherits(SG, "SGWT")) {
    stop("Input must be an SGWT object")
  }
  if (is.null(SG$Forward) || is.null(SG$Inverse)) {
    stop("SGWT object must have Forward and Inverse results computed")
  }
  
  # Default to first signal if not specified
  if (is.null(signal_name)) {
    signal_name <- names(SG$Forward)[1]
  }
  
  # Validate signal exists
  if (!signal_name %in% names(SG$Forward)) {
    stop(paste("Signal", signal_name, "not found in SGWT results"))
  }
  
  # Get decomposition and inverse results
  inverse_result <- SG$Inverse[[signal_name]]
  # Get coefficients from inverse results (vertex_approximations)
  coefficients <- inverse_result$vertex_approximations
  
  # Default scales to plot
  if (is.null(plot_scales)) {
    n_wavelets <- length(coefficients) - 1  # Exclude scaling
    plot_scales <- 1:min(4, n_wavelets)
  }
  
  # Prepare data for plotting
  data.in <- SG$Data$data
  x_col <- SG$Data$x_col
  y_col <- SG$Data$y_col
  
  # Create a helper function to create individual plots
  create_plot <- function(data, x_col, y_col, fill_var, title, subtitle = NULL) {
    # Use aes_string for compatibility and to avoid linting issues
    p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x_col, y = y_col, fill = fill_var)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(title = title, subtitle = subtitle) +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = 8, hjust = 0.5)
      )
    return(p)
  }
  
  plot_list <- list()
  
  # Plot original signal
  plot_data_orig <- data.in
  plot_data_orig$original <- as.numeric(data.in[[signal_name]])
  p_orig <- create_plot(plot_data_orig, x_col, y_col, "original", 
                       paste("Original Signal:", signal_name))
  plot_list[["original"]] <- p_orig
  
  # Plot scaling function coefficients (now called low_pass)
  plot_data_scaling <- data.in
  plot_data_scaling$scaling <- as.vector(Re(coefficients$low_pass))
  p_scaling <- create_plot(plot_data_scaling, x_col, y_col, "scaling", 
                          "Low-pass (Scaling)")
  plot_list[["scaling"]] <- p_scaling
  
  # Plot wavelet coefficients at selected scales (now called wavelet_1, wavelet_2, etc.)
  wavelet_names <- names(coefficients)[grep("^wavelet_", names(coefficients))]
  for (i in plot_scales) {
    wavelet_name <- paste0("wavelet_", i)
    if (wavelet_name %in% wavelet_names) {
      coeff_name <- paste0("wavelet_", i)
      plot_data_wavelet <- data.in
      plot_data_wavelet[[coeff_name]] <- as.vector(Re(coefficients[[wavelet_name]]))
      
      p_wavelet <- create_plot(plot_data_wavelet, x_col, y_col, coeff_name, 
                              paste("Band-pass Scale", i))
      plot_list[[coeff_name]] <- p_wavelet
    }
  }
  
  # Plot reconstructed signal
  plot_data_recon <- data.in
  plot_data_recon$reconstructed <- inverse_result$reconstructed_signal
  p_recon <- create_plot(plot_data_recon, x_col, y_col, "reconstructed", 
                        "Reconstructed", 
                        paste("RMSE:", round(inverse_result$reconstruction_error, 4)))
  plot_list[["reconstructed"]] <- p_recon
  
  # Validate that we have plots to combine
  if (length(plot_list) == 0) {
    stop("No plots were created. Check your SGWT object structure.")
  }
  
  # Combine plots using gridExtra (most reliable)
  n_plots <- length(plot_list)
  nrow <- ceiling(n_plots / ncol)
  
  # Use gridExtra::grid.arrange for reliable plot combination
  combined_plot <- gridExtra::grid.arrange(grobs = plot_list, ncol = ncol, nrow = nrow)
  
  # Return the combined plot
  return(combined_plot)
}

#' Analyze SGWT energy distribution across scales in Fourier domain
#'
#' @description Calculate and analyze energy distribution across different scales
#' using Fourier domain coefficients directly (consistent with Parseval's theorem).
#' Excludes DC component for more accurate energy analysis.
#'
#' @param SG SGWT object with Forward results computed
#' @param signal_name Name of signal to analyze (default: first signal)
#'
#' @return Data frame with energy analysis results computed in Fourier domain
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have SGWT object
#' energy_analysis <- sgwt_energy_analysis(SG_object, signal_name = "signal1")
#' print(energy_analysis)
#' }
sgwt_energy_analysis <- function(SG, signal_name = NULL) {
  
  # Validate input
  if (!inherits(SG, "SGWT")) {
    stop("Input must be an SGWT object")
  }
  if (is.null(SG$Forward)) {
    stop("SGWT object must have Forward results computed")
  }
  
  # Default to first signal if not specified
  if (is.null(signal_name)) {
    signal_name <- names(SG$Forward)[1]
  }
  
  # Validate signal exists
  if (!signal_name %in% names(SG$Forward)) {
    stop(paste("Signal", signal_name, "not found in SGWT Forward results"))
  }
  
  # Get Forward results and scales from Parameters
  forward_result <- SG$Forward[[signal_name]]
  fourier_coeffs <- forward_result$fourier_coefficients$filtered
  scales <- SG$Parameters$scales
  
  if (is.null(fourier_coeffs)) {
    stop("Fourier coefficients not found in Forward results")
  }
  
  # Calculate energies in Fourier domain (consistent with Parseval's theorem)
  energies <- numeric()
  scale_names <- character()
  scale_values <- numeric()
  
  # Scaling (low-pass) energy - exclude DC component
  if ("scaling" %in% names(fourier_coeffs)) {
    scaling_coeffs <- as.numeric(fourier_coeffs$scaling)
    # Exclude DC component (first coefficient)
    if (length(scaling_coeffs) > 1) {
      scaling_coeffs <- scaling_coeffs[-1]
    }
    scaling_energy <- sum(abs(scaling_coeffs)^2)
    
    energies <- c(energies, scaling_energy)
    scale_names <- c(scale_names, "low_pass")
    scale_values <- c(scale_values, scales[1])  # Use first scale for scaling function
  }
  
  # Wavelet energies - exclude DC components
  wavelet_names <- names(fourier_coeffs)[grep("^wavelet_scale_", names(fourier_coeffs))]
  if (length(wavelet_names) > 0) {
    # Order by scale index
    scale_indices <- as.integer(sub("^wavelet_scale_", "", wavelet_names))
    ord <- order(scale_indices)
    wavelet_names <- wavelet_names[ord]
    scale_indices <- scale_indices[ord]
    
    for (i in seq_along(wavelet_names)) {
      wavelet_coeffs <- as.numeric(fourier_coeffs[[wavelet_names[i]]])
      # Exclude DC component if present
      if (length(wavelet_coeffs) > 1) {
        wavelet_coeffs <- wavelet_coeffs[-1]
      }
      wavelet_energy <- sum(abs(wavelet_coeffs)^2)
      
      energies <- c(energies, wavelet_energy)
      scale_names <- c(scale_names, paste0("wavelet_", scale_indices[i]))
      scale_values <- c(scale_values, scales[scale_indices[i]])
    }
  }
  
  # Calculate energy ratios
  total_energy <- sum(energies)
  energy_ratios <- if (total_energy > 0) energies / total_energy else rep(0, length(energies))
  
  # Create results data frame
  energy_df <- data.frame(
    scale = scale_names,
    energy = energies,
    energy_ratio = energy_ratios,
    scale_value = scale_values,
    signal = signal_name,
    stringsAsFactors = FALSE
  )
  
  return(energy_df)
}

#' Plot Fourier modes (eigenvectors) from SGWT object
#'
#' @description Plot low-frequency and high-frequency Fourier modes (eigenvectors) 
#' from the graph Laplacian eigendecomposition in an SGWT object
#'
#' @param SG SGWT object with Graph slot computed (from runSpecGraph)
#' @param mode_type Type of modes to plot: "low", "high", or "both" (default: "both")
#' @param n_modes Number of modes to plot for each type (default: 6)
#' @param ncol Number of columns in plot layout (default: 3)
#' @param point_size Size of points in the plot (default: 1.5)
#'
#' @return Combined plot of Fourier modes
#' @export
#'
#' @examples
#' \dontrun{
#' # Plot both low and high frequency modes
#' SG <- initSGWT(data) %>% runSpecGraph()
#' plot_FM(SG, mode_type = "both", n_modes = 4)
#' 
#' # Plot only low frequency modes
#' plot_FM(SG, mode_type = "low", n_modes = 8)
#' }
plot_FM <- function(SG, mode_type = "both", n_modes = 6, ncol = 3, point_size = 1.5){
  
  # Validate input
  if (!inherits(SG, "SGWT")) {
    stop("Input must be an SGWT object")
  }
  if (is.null(SG$Graph)) {
    stop("SGWT object must have Graph slot computed. Run runSpecGraph() first.")
  }
  
  # Extract components
  data.in <- SG$Data$data
  x_col <- SG$Data$x_col
  y_col <- SG$Data$y_col
  eigenvalues <- SG$Graph$eigenvalues
  eigenvectors <- as.matrix(SG$Graph$eigenvectors)
  
  # Validate mode_type
  mode_type <- match.arg(mode_type, c("low", "high", "both"))
  
  # Determine which modes to plot based on eigenvalue spectrum
  n_total <- length(eigenvalues)
  n_modes <- min(n_modes, floor(n_total/2))  # Ensure we don't exceed available modes
  
  plot_list <- list()
  
  # Helper function to create individual Fourier mode plots
  create_mode_plot <- function(mode_data, mode_name, eigenval) {
    p <- ggplot2::ggplot(mode_data, ggplot2::aes_string(x = x_col, y = y_col, color = "mode_value")) +
      ggplot2::geom_point(size = point_size) +
      ggplot2::scale_color_viridis_c(option = "plasma") +
      ggplot2::labs(
        title = mode_name,
        subtitle = paste("\u03bb =", round(eigenval, 4))
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 8, hjust = 0.5)
      ) +
      ggplot2::coord_fixed()
    
    return(p)
  }
  
  # Plot low-frequency modes (smallest eigenvalues, skip DC component)
  if (mode_type %in% c("low", "both")) {
    low_indices <- 2:(n_modes + 1)  # Skip first mode (DC component), start from 2nd
    
    for (i in low_indices) {
      if (i <= n_total) {
        mode_data <- data.in
        mode_data$mode_value <- as.numeric(eigenvectors[, i])
        
        mode_name <- paste("Low Freq", i)
        eigenval <- eigenvalues[i]
        
        p <- create_mode_plot(mode_data, mode_name, eigenval)
        plot_list[[paste0("low_", i)]] <- p
      }
    }
  }
  
  # Plot high-frequency modes (largest eigenvalues)
  if (mode_type %in% c("high", "both")) {
    high_indices <- (n_total - n_modes + 1):n_total  # Last n_modes (highest frequencies)
    
    for (i in high_indices) {
      if (i >= 1) {
        mode_data <- data.in
        mode_data$mode_value <- as.numeric(eigenvectors[, i])
        
        mode_name <- paste("High Freq", i)
        eigenval <- eigenvalues[i]
        
        p <- create_mode_plot(mode_data, mode_name, eigenval)
        plot_list[[paste0("high_", i)]] <- p
      }
    }
  }
  
  # Validate that we have plots to combine
  if (length(plot_list) == 0) {
    stop("No plots were created. Check your SGWT object and parameters.")
  }
  
  # Create title based on mode_type
  main_title <- switch(mode_type,
                      "low" = paste("Low-Frequency Fourier Modes (n =", n_modes, ")"),
                      "high" = paste("High-Frequency Fourier Modes (n =", n_modes, ")"),
                      "both" = paste("Fourier Modes: Low &amp; High Frequency (n =", n_modes, "each)"))
  
  # Combine plots
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    # Calculate appropriate number of rows
    n_plots <- length(plot_list)
    nrow <- ceiling(n_plots / ncol)
    
    # Add main title
    title_grob <- grid::textGrob(main_title, 
                                gp = grid::gpar(fontsize = 14, fontface = "bold"))
    
    combined_plot <- gridExtra::grid.arrange(
      grobs = plot_list, 
      ncol = ncol, 
      nrow = nrow,
      top = title_grob
    )
  } else {
    stop("gridExtra package is required for plot combination. Please install it.")
  }
  
  return(combined_plot)
}

#' Visualize SGWT kernels and scaling functions
#'
#' @description Visualize the scaling function and wavelet kernels used in SGWT
#' based on the eigenvalue spectrum and selected parameters
#'
#' @param eigenvalues Vector of eigenvalues from graph Laplacian
#' @param scales Vector of scales for the wavelets (if NULL, auto-generated)
#' @param J Number of scales to generate if scales is NULL (default: 4)
#' @param scaling_factor Scaling factor between consecutive scales (default: 2)
#' @param kernel_type Type of wavelet kernel ("mexican_hat" or "meyer", default: "mexican_hat")
#' @param lmax Maximum eigenvalue (optional, computed if NULL)
#' @param eigenvalue_range Range of eigenvalues to plot (default: full range)
#' @param resolution Number of points for smooth curve plotting (default: 1000)
#'
#' @return List containing the filter visualization plot and filter values
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate some example eigenvalues
#' eigenvals <- seq(0, 2, length.out = 100)
#' 
#' # Visualize kernels with specific parameters
#' viz_result <- visualize_sgwt_kernels(
#'   eigenvalues = eigenvals,
#'   J = 4,
#'   scaling_factor = 2,
#'   kernel_type = "heat"
#' )
#' print(viz_result$plot)
#' }
visualize_sgwt_kernels <- function(eigenvalues, scales = NULL, J = 4, scaling_factor = 2,
                                   kernel_type = "heat", lmax = NULL,
                                   eigenvalue_range = NULL, resolution = 1000) {
  
  # Set lmax if not provided
  if (is.null(lmax)) {
    lmax <- max(eigenvalues) * 0.95
  }
  
  # Auto-generate scales if not provided
  if (is.null(scales)) {
    scales <- sgwt_auto_scales(lmax, J, scaling_factor)
  }
  
  # Set eigenvalue range for plotting
  if (is.null(eigenvalue_range)) {
    eigenvalue_range <- c(0, max(eigenvalues))
  }
  
  # Create smooth eigenvalue sequence for plotting
  lambda_smooth <- seq(eigenvalue_range[1], eigenvalue_range[2], length.out = resolution)
  
  # Compute filters for smooth sequence
  filters_smooth <- compute_sgwt_filters(lambda_smooth, scales, lmax)
  
  # Prepare data for plotting
  plot_data <- data.frame(
    eigenvalue = rep(lambda_smooth, length(filters_smooth)),
    filter_value = unlist(filters_smooth),
    filter_type = rep(c("Scaling Function", paste("Wavelet Scale", seq_along(scales))), 
                     each = length(lambda_smooth)),
    scale_param = rep(c(scales[1], scales), each = length(lambda_smooth))
  )
  
  # Create color palette
  n_filters <- length(filters_smooth)
  colors <- c("#E74C3C", viridis::viridis(n_filters - 1))
  
  # Create the plot
  p_kernels <- ggplot2::ggplot(plot_data, ggplot2::aes(x = eigenvalue, y = filter_value, 
                                                       color = filter_type)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = "SGWT Filter Bank: Scaling Function and Wavelet Kernels",
      subtitle = paste("Kernel Type:", kernel_type, "| J =", length(scales), 
                      "| Scaling Factor =", scaling_factor),
      x = "Eigenvalue (\u03bb)",
      y = "Filter Response",
      color = "Filter Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 12),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))
  
  # Add vertical lines for actual eigenvalues (sample)
  if (length(eigenvalues) <= 50) {
    eigenvalue_sample <- eigenvalues
  } else {
    eigenvalue_sample <- eigenvalues[seq(1, length(eigenvalues), length.out = 50)]
  }
  
  p_kernels <- p_kernels +
    ggplot2::geom_vline(xintercept = eigenvalue_sample, alpha = 0.3, color = "gray60", size = 0.3)
  
  # Create eigenvalue histogram subplot
  eigenval_data <- data.frame(eigenvalue = eigenvalues)
  p_eigenvals <- ggplot2::ggplot(eigenval_data, ggplot2::aes(x = eigenvalue)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
    ggplot2::labs(
      title = "Eigenvalue Distribution",
      x = "Eigenvalue (\u03bb)",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    )
  
  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined_plot <- p_kernels / p_eigenvals + patchwork::plot_layout(heights = c(3, 1))
  } else {
    combined_plot <- p_kernels
    cat("Note: Install 'patchwork' package to see eigenvalue distribution subplot\n")
  }
  
  # Compute filters for actual eigenvalues
  filters_actual <- compute_sgwt_filters(eigenvalues, scales, lmax)
  
  # Return results
  return(list(
    plot = combined_plot,
    filters_smooth = filters_smooth,
    filters_actual = filters_actual,
    lambda_smooth = lambda_smooth,
    scales = scales,
    parameters = list(
      J = length(scales),
      scaling_factor = scaling_factor,
      kernel_type = kernel_type,
      lmax = lmax
    )
  ))
}

#' Visualize similarity in low vs non-low frequency space
#'
#' @description Create a scatter plot with low-frequency similarity (c_low) on x-axis
#' and non-low-frequency similarity (c_nonlow) on y-axis from runSGCC results
#'
#' @param similarity_results List of similarity results from runSGCC function, or a single result
#' @param point_size Size of points in the plot (default: 2)
#' @param point_color Color of points (default: "steelblue")
#' @param add_diagonal Whether to add diagonal reference lines (default: TRUE)
#' @param add_axes_lines Whether to add x=0 and y=0 reference lines (default: TRUE)
#' @param title Plot title (default: "Low-frequency vs Non-low-frequency Similarity")
#' @param show_labels Whether to show point labels if names are available (default: FALSE)
#' @param show_names Whether to display data point names as text labels using ggrepel (default: FALSE).
#'   If more than 50 points, randomly samples 50 for labeling. Requires ggrepel package.
#'
#' @return ggplot object showing similarity space visualization
#' @export
#'
#' @examples
#' \dontrun{
#' # Single similarity result
#' sim_result <- runSGCC("signal1", "signal2", SG = SG_object)
#' plot <- visualize_similarity_xy(sim_result)
#' print(plot)
#' 
#' # Multiple similarity results
#' sim_results <- list(
#'   pair1 = runSGCC("signal1", "signal2", SG = SG_object1),
#'   pair2 = runSGCC("signal1", "signal2", SG = SG_object2)
#' )
#' plot <- visualize_similarity_xy(sim_results, show_names = TRUE)
#' print(plot)
#' 
#' # Show both labels and names (for comparison)
#' plot_both <- visualize_similarity_xy(sim_results, show_labels = TRUE, show_names = TRUE)
#' print(plot_both)
#' 
#' # With many data points (>50), names will be randomly sampled
#' # install.packages("ggrepel")  # Required for show_names = TRUE
#' plot_many <- visualize_similarity_xy(many_sim_results, show_names = TRUE)
#' print(plot_many)
#' }
visualize_similarity_xy <- function(similarity_results, 
                                   point_size = 2,
                                   point_color = "steelblue",
                                   add_diagonal = TRUE,
                                   add_axes_lines = TRUE,
                                   title = "Low-frequency vs Non-low-frequency Similarity",
                                   show_labels = FALSE,
                                   show_names = FALSE) {
  
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for visualization")
  }
  if (show_names && !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggrepel package is required when show_names = TRUE. Please install it with: install.packages('ggrepel')")
  }
  
  # Handle single result vs list of results
  if (is.list(similarity_results) && !is.null(similarity_results$c_low)) {
    # Single result - convert to list
    similarity_results <- list(result = similarity_results)
  }
  
  # Validate input structure
  if (!is.list(similarity_results)) {
    stop("similarity_results must be a list or a single runSGCC result")
  }
  
  # Extract c_low and c_nonlow values
  plot_data <- data.frame(
    c_low = numeric(0),
    c_nonlow = numeric(0),
    label = character(0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(similarity_results)) {
    result <- similarity_results[[i]]
    
    # Validate that result has required components
    if (is.null(result$c_low) || is.null(result$c_nonlow)) {
      warning(paste("Result", i, "missing c_low or c_nonlow components, skipping"))
      next
    }
    
    # Add to plot data
    plot_data <- rbind(plot_data, data.frame(
      c_low = result$c_low,
      c_nonlow = result$c_nonlow,
      label = if (is.null(names(similarity_results)[i])) paste("Point", i) else names(similarity_results)[i],
      stringsAsFactors = FALSE
    ))
  }
  
  # Check if we have data to plot
  if (nrow(plot_data) == 0) {
    stop("No valid similarity results found to plot")
  }
  
  # Create the base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "c_low", y = "c_nonlow")) +
    ggplot2::geom_point(size = point_size, color = point_color, alpha = 0.7) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(
      title = title,
      x = "Low-frequency Similarity (c_low)",
      y = "Non-low-frequency Similarity (c_nonlow)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Add reference lines if requested
  if (add_axes_lines) {
    p <- p + 
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray60", alpha = 0.7)
  }
  
  if (add_diagonal) {
    p <- p + 
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray40", alpha = 0.7) +
      ggplot2::geom_abline(slope = -1, intercept = 0, linetype = "dotted", color = "gray40", alpha = 0.7)
  }
  
  # Add labels if requested (backward compatibility)
  if (show_labels && nrow(plot_data) > 0) {
    p <- p + ggplot2::geom_text(ggplot2::aes_string(label = "label"), 
                               vjust = -0.5, hjust = 0.5, size = 3, color = "black")
  }
  
  # Add names if requested (new parameter with ggrepel)
  if (show_names && nrow(plot_data) > 0) {
    # Create a subset for labeling if there are too many points
    label_data <- plot_data
    n_points <- nrow(plot_data)
    
    if (n_points > 50) {
      # Random sample 50 points for labeling to avoid overcrowding
      set.seed(123)  # For reproducible sampling
      sample_indices <- sample(seq_len(n_points), size = 50, replace = FALSE)
      label_data <- plot_data[sample_indices, ]
      
      # Add a note about sampling
      subtitle_text <- paste("Showing", nrow(label_data), "of", n_points, "labels (random sample)")
      p <- p + ggplot2::labs(subtitle = subtitle_text)
    }
    
    # Use ggrepel for better text positioning
    p <- p + ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes_string(label = "label"),
      size = 2.5,
      color = "darkblue",
      fontface = "bold",
      box.padding = 0.35,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.2,
      max.overlaps = Inf,
      min.segment.length = 0.1
    )
  }
  
  return(p)
}

#' Demo function for SGWT
#'
#' @description Demonstration function showing basic SGWT usage with synthetic data
#' using the new workflow: initSGWT -> runSpecGraph -> runSGWT
#'
#' @return SGWT object with complete analysis
#' @export
#'
#' @examples
#' \dontrun{
#' SG <- demo_sgwt()
#' print(SG)
#' }
demo_sgwt <- function() {
  cat("=== SGWT Demo ===\n")
  
  # Generate synthetic spatial data
  set.seed(123)
  n_points <- 100
  
  # Create a simple 2D grid with some noise
  x_coords <- rep(1:10, each = 10) + stats::rnorm(n_points, 0, 0.1)
  y_coords <- rep(1:10, times = 10) + stats::rnorm(n_points, 0, 0.1)
  
  # Create synthetic signals
  signal1 <- sin(0.5 * x_coords) * cos(0.3 * y_coords) + stats::rnorm(n_points, 0, 0.1)
  signal2 <- 0.5 * sin(0.8 * x_coords * y_coords) + stats::rnorm(n_points, 0, 0.1)
  
  # Create data frame
  demo_data <- data.frame(
    x = x_coords,
    y = y_coords,
    signal1 = signal1,
    signal2 = signal2
  )
  
  cat("Generated synthetic data with", n_points, "points and", 2, "signals\n")
  
  # New SGWT workflow
  cat("Step 1: Initialize SGWT object\n")
  SG <- initSGWT(demo_data, signals = c("signal1", "signal2"), J = 4)
  
  cat("Step 2: Build spectral graph\n")
  SG <- runSpecGraph(SG, verbose = TRUE)
  
  cat("Step 3: Run SGWT analysis\n")
  SG <- runSGWT(SG, verbose = TRUE)
  
  cat("Step 4: Display results\n")
  print(SG)
  
  # Display energy analysis for first signal
  energy_analysis <- sgwt_energy_analysis(SG, "signal1")
  cat("\nEnergy analysis for signal1:\n")
  print(energy_analysis)
  
  cat("\n=== SGWT Demo Complete ===\n")
  
  return(SG)
} 