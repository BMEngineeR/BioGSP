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
  forward_result <- SG$Forward[[signal_name]]
  inverse_result <- SG$Inverse[[signal_name]]
  coefficients <- forward_result$coefficients
  
  # Default scales to plot
  if (is.null(plot_scales)) {
    n_wavelets <- length(coefficients) - 1  # Exclude scaling
    plot_scales <- 1:min(4, n_wavelets)
  }
  
  # Prepare data for plotting
  data.in <- SG$Data$data
  x_col <- SG$Data$x_col
  y_col <- SG$Data$y_col
  plot_data <- data.in
  
  plot_list <- list()
  
  # Plot original signal
  plot_data$original <- as.numeric(data.in[[signal_name]])
  p_orig <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = x_col, y = y_col, fill = "original")) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title = paste("Original Signal:", signal_name)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
  plot_list[["original"]] <- p_orig
  
  # Plot scaling function coefficients
  plot_data$scaling <- as.vector(Re(coefficients$scaling))
  p_scaling <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = x_col, y = y_col, fill = "scaling")) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title = "Low-pass (Scaling)") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
  plot_list[["scaling"]] <- p_scaling
  
  # Plot wavelet coefficients at selected scales
  wavelet_names <- names(coefficients)[grep("^wavelet_scale_", names(coefficients))]
  for (i in plot_scales) {
    wavelet_name <- paste0("wavelet_scale_", i)
    if (wavelet_name %in% wavelet_names) {
      coeff_name <- paste0("wavelet_", i)
      plot_data[[coeff_name]] <- as.vector(Re(coefficients[[wavelet_name]]))
      
      p_wavelet <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = x_col, y = y_col, fill = coeff_name)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::labs(title = paste("Band-pass Scale", i)) +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "none")
      
      plot_list[[coeff_name]] <- p_wavelet
    }
  }
  
  # Plot reconstructed signal
  plot_data$reconstructed <- inverse_result$reconstructed_signal
  p_recon <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = x_col, y = y_col, fill = "reconstructed")) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title = paste("Reconstructed (RMSE:", round(inverse_result$reconstruction_error, 4), ")")) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")
  plot_list[["reconstructed"]] <- p_recon
  
  # Combine plots
  combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol)
  return(combined_plot)
}

#' Analyze SGWT energy distribution across scales
#'
#' @description Calculate and analyze energy distribution across different scales
#' in the SGWT decomposition
#'
#' @param SG SGWT object with Forward results computed
#' @param signal_name Name of signal to analyze (default: first signal)
#'
#' @return Data frame with energy analysis results
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
  
  # Get forward result
  forward_result <- SG$Forward[[signal_name]]
  coefficients <- forward_result$coefficients
  scales <- forward_result$scales
  
  # Calculate energy at each scale
  energies <- sapply(coefficients, function(coeff) sum(abs(coeff)^2))
  total_energy <- sum(energies)
  energy_ratios <- energies / total_energy
  
  # Create results data frame
  scale_names <- names(coefficients)
  scale_values <- c(scales[1], scales)  # Use first scale for scaling function
  
  energy_df <- data.frame(
    scale = scale_names,
    energy = energies,
    energy_ratio = energy_ratios,
    scale_value = scale_values,
    signal = signal_name
  )
  
  return(energy_df)
}

#' Plot frequency modes
#'
#' @description Plot frequency modes from graph Fourier analysis
#'
#' @param input Input data (currently not used, for future compatibility)
#' @param FM_idx Indices of frequency modes to plot (default: 1:20)
#' @param ncol Number of columns in plot layout (default: 5)
#'
#' @return Combined plot of frequency modes
#' @export
#'
#' @examples
#' \dontrun{
#' # This function requires specific data structure (df_hex_combine)
#' # plot_FM(FM_idx = 1:10, ncol = 5)
#' }
plot_FM <- function(input = NULL, FM_idx = c(1:20), ncol = 5){
  
  # Create a vector of frequency column names
  freq_columns <- paste0("Freq", FM_idx)
  # Create a list to store the plots
  plot_list <- list()
  
  # Note: This function assumes existence of df_hex_combine
  # This is kept for compatibility with original code
  if (!exists("df_hex_combine")) {
    warning("df_hex_combine not found. This function requires specific data structure.")
    return(NULL)
  }
  
  # Loop through each frequency column and generate plots
  for (freq in freq_columns) {
    # Generate the scatter plot
    p <- ggplot2::ggplot(df_hex_combine, ggplot2::aes_string(x = "x", y = "y")) +
      ggplot2::geom_point(ggplot2::aes_string(color = freq), size = 1) +
      ggplot2::guides(color = "none") +
      ggplot2::theme_void() + 
      ggplot2::scale_color_viridis_c(option = "magma")
    
    # Add the plot to the list
    plot_list[[freq]] <- p
  }
  
  # Combine all plots into a single layout
  combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol)
  # Print the combined plot
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
  SG <- initSGWT(demo_data, signals = c("signal1", "signal2"), k = 8, J = 4)
  
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