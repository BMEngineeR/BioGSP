

#' Get a unified kernel family (low-pass and band-pass) by kernel_type
#'
#' @description Returns a pair of functions implementing the scaling (low-pass)
#' and wavelet (band-pass) kernels for a given kernel family. This enforces
#' consistency: a single kernel_type defines both filters.
#'
#' @param kernel_type Kernel family name ("mexican_hat", "meyer", or "heat")
#'
#' @return A list with two functions: list(scaling = function(x, scale_param), wavelet = function(x, scale_param))
#' @export
sgwt_get_kernels <- function(kernel_type = "mexican_hat") {
  if (kernel_type == "mexican_hat") {
    scaling_fun <- function(x, scale_param) {
      t <- x / scale_param
      exp(-0.5 * t^2)
    }
    wavelet_fun <- function(x, scale_param) {
      t <- x / scale_param
      (t^2) * exp(-0.5 * t^2)
    }
    return(list(scaling = scaling_fun, wavelet = wavelet_fun))
  } else if (kernel_type == "meyer") {
    scaling_fun <- function(x, scale_param) {
      a <- 0.5 * scale_param
      b <- 1.0 * scale_param
      if (x <= a) {
        1
      } else if (x >= b) {
        0
      } else {
        s <- (x - a) / (b - a)
        cos(pi/2 * s)
      }
    }
    wavelet_fun <- function(x, scale_param) {
      a <- 0.5 * scale_param
      b <- 1.0 * scale_param
      c <- 2.0 * scale_param
      if (x <= a || x >= c) {
        0
      } else if (x <= b) {
        s <- (x - a) / (b - a)
        sin(pi/2 * s)
      } else {
        s <- (x - b) / (c - b)
        cos(pi/2 * s)
      }
    }
    return(list(scaling = scaling_fun, wavelet = wavelet_fun))
  } else if (kernel_type == "heat") {
    scaling_fun <- function(x, scale_param) {
      # Heat kernel scaling: exponential decay h(t) = exp(-t)
      t <- x / scale_param
      exp(-t)
    }
    wavelet_fun <- function(x, scale_param) {
      # Heat kernel wavelet: derivative-like g(t) = t * exp(-t)
      t <- x / scale_param
      t * exp(-t)
    }
    return(list(scaling = scaling_fun, wavelet = wavelet_fun))
  } else {
    stop("Kernel type not supported. Use 'mexican_hat', 'meyer', or 'heat'.")
  }
}

#' Compute SGWT filters
#'
#' @description Compute wavelet and scaling function coefficients in the spectral domain
#'
#' @param eigenvalues Eigenvalues of the graph Laplacian
#' @param scales Vector of scales for the wavelets
#' @param lmax Maximum eigenvalue (optional)
#' @param kernel_type Kernel family that defines both scaling and wavelet filters (default: "mexican_hat", options: "mexican_hat", "meyer", "heat") 
#'
#' @return List of filters (scaling function + wavelets)
#' @export
#'
#' @examples
#' eigenvals <- c(0, 0.1, 0.5, 1.0, 1.5)
#' scales <- c(2, 1, 0.5)
#' filters <- compute_sgwt_filters(eigenvals, scales)
#' filters_meyer <- compute_sgwt_filters(eigenvals, scales, kernel_type = "meyer")
#' filters_heat <- compute_sgwt_filters(eigenvals, scales, kernel_type = "heat")
compute_sgwt_filters <- function(eigenvalues, scales, lmax = NULL, kernel_type = "mexican_hat") {
  if (is.null(lmax)) {
    lmax <- max(eigenvalues) * 0.95  # Avoid numerical issues at lambda_max
  }
  
  J <- length(scales)  # Number of scales
  
  # Initialize filter bank: J wavelets + 1 scaling function
  filters <- vector("list", J + 1)

  # Get unified kernel family
  kernels <- sgwt_get_kernels(kernel_type)
  
  # Scaling function (low-pass filter)
  filters[[1]] <- sapply(eigenvalues, function(lambda) {
    kernels$scaling(lambda, scale_param = scales[1])
  })
  
  # Wavelet functions at different scales
  for (j in seq_len(J)) {
    filters[[j + 1]] <- sapply(eigenvalues, function(lambda) {
      kernels$wavelet(lambda, scale_param = scales[j])
    })
  }
  
  return(filters)
}

#' Forward SGWT transform
#'
#' @description Decompose signal into wavelet coefficients using SGWT
#'
#' @param signal Input signal vector
#' @param eigenvectors Eigenvectors of the graph Laplacian
#' @param eigenvalues Eigenvalues of the graph Laplacian
#' @param scales Vector of scales for the wavelets
#' @param lmax Maximum eigenvalue (optional)
#' @param kernel_type Kernel family that defines both scaling and wavelet filters (default: "mexican_hat", options: "mexican_hat", "meyer", "heat")
#'
#' @return List containing coefficients, filters, scales, eigenvalues, and eigenvectors
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have eigenvalues, eigenvectors, and a signal
#' result <- sgwt_forward(signal, eigenvectors, eigenvalues, scales)
#' result_meyer <- sgwt_forward(signal, eigenvectors, eigenvalues, scales, kernel_type = "meyer")
#' result_heat <- sgwt_forward(signal, eigenvectors, eigenvalues, scales, kernel_type = "heat")
#' }
sgwt_forward <- function(signal, eigenvectors, eigenvalues, scales, lmax = NULL, kernel_type = "mexican_hat") {
  # Compute filters
  filters <- compute_sgwt_filters(eigenvalues, scales, lmax, kernel_type)
  
  # Transform signal to spectral domain
  signal_hat <- gft(signal, eigenvectors)
  
  # Apply filters and store coefficients
  coefficients <- vector("list", length(filters))
  
  for (i in seq_along(filters)) {
    # Multiply signal spectrum with filter
    filtered_spectrum <- signal_hat * filters[[i]]
    
    # Transform back to vertex domain (inverse GFT)
    coefficients[[i]] <- eigenvectors %*% filtered_spectrum
  }
  
  names(coefficients) <- c("scaling", paste0("wavelet_scale_", seq_along(scales)))
  
  return(list(
    coefficients = coefficients,
    filters = filters,
    scales = scales,
    eigenvalues = eigenvalues,
    eigenvectors = eigenvectors
  ))
}

#' Inverse SGWT transform
#'
#' @description Reconstruct signal from wavelet coefficients
#'
#' @param sgwt_decomp SGWT decomposition object from sgwt_forward
#'
#' @return Reconstructed signal vector
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have an SGWT decomposition
#' reconstructed <- sgwt_inverse(sgwt_decomp)
#' }
sgwt_inverse <- function(sgwt_decomp) {
  coefficients <- sgwt_decomp$coefficients
  filters <- sgwt_decomp$filters
  eigenvectors <- sgwt_decomp$eigenvectors
  
  # Initialize reconstructed signal
  reconstructed_signal <- rep(0, nrow(eigenvectors))
  
  # Reconstruct by summing all filtered components
  for (i in seq_along(coefficients)) {
    # Transform coefficient to spectral domain
    coeff_hat <- gft(coefficients[[i]], eigenvectors)
    
    # Apply filter
    filtered_coeff <- coeff_hat * filters[[i]]
    
    # Transform back and add to reconstruction
    reconstructed_signal <- reconstructed_signal + as.vector(eigenvectors %*% filtered_coeff)
  }
  
  return(reconstructed_signal)
}

#' Generate automatic scales for SGWT
#'
#' @description Generate logarithmically spaced scales for SGWT
#'
#' @param lmax Maximum eigenvalue
#' @param J Number of scales
#' @param scaling_factor Scaling factor between consecutive scales
#'
#' @return Vector of scales
#' @export
#'
#' @examples
#' scales <- sgwt_auto_scales(lmax = 2.0, J = 5, scaling_factor = 2)
sgwt_auto_scales <- function(lmax, J = 5, scaling_factor = 2) {
  # Generate logarithmically spaced scales
  scales <- lmax / (scaling_factor^(0:(J - 1)))
  return(scales)
}

#' Compare different kernel families
#'
#' @description Visualize and compare different kernel families (both scaling and wavelet filters)
#'
#' @param x_range Range of x values to evaluate (default: c(0, 3))
#' @param scale_param Scale parameter for all functions (default: 1)
#' @param plot_results Whether to plot the comparison (default: TRUE)
#'
#' @return Data frame with x values and kernel values for each family
#' @export
#'
#' @examples
#' comparison <- compare_kernel_families()
#' comparison <- compare_kernel_families(x_range = c(0, 5), scale_param = 1.5)
compare_kernel_families <- function(x_range = c(0, 3), scale_param = 1, plot_results = TRUE) {
  x_vals <- seq(x_range[1], x_range[2], length.out = 200)
  
  kernel_types <- c("mexican_hat", "meyer", "heat")
  
  results <- data.frame(x = x_vals)
  
  for (type in kernel_types) {
    kernels <- sgwt_get_kernels(type)
    results[[paste0(type, "_scaling")]] <- sapply(x_vals, function(x) {
      kernels$scaling(x, scale_param = scale_param)
    })
    results[[paste0(type, "_wavelet")]] <- sapply(x_vals, function(x) {
      kernels$wavelet(x, scale_param = scale_param)
    })
  }
  
  if (plot_results && requireNamespace("graphics", quietly = TRUE)) {
    graphics::par(mfrow = c(1, 2))
    
    # Plot scaling functions
    scaling_max <- max(results$mexican_hat_scaling, results$meyer_scaling, results$heat_scaling, na.rm = TRUE)
    graphics::plot(x_vals, results$mexican_hat_scaling, type = "l", col = 1, lwd = 2, 
                   xlab = "Eigenvalue", ylab = "Scaling Function Value",
                   main = "SGWT Scaling Functions",
                   ylim = c(0, scaling_max * 1.1))
    graphics::lines(x_vals, results$meyer_scaling, col = 2, lwd = 2)
    graphics::lines(x_vals, results$heat_scaling, col = 3, lwd = 2)
    graphics::legend("topright", legend = c("Mexican Hat", "Meyer", "Heat"), col = 1:3, lwd = 2)
    
    # Plot wavelet functions  
    wavelet_max <- max(results$mexican_hat_wavelet, results$meyer_wavelet, results$heat_wavelet, na.rm = TRUE)
    graphics::plot(x_vals, results$mexican_hat_wavelet, type = "l", col = 1, lwd = 2,
                   xlab = "Eigenvalue", ylab = "Wavelet Function Value", 
                   main = "SGWT Wavelet Functions",
                   ylim = c(0, wavelet_max * 1.1))
    graphics::lines(x_vals, results$meyer_wavelet, col = 2, lwd = 2)
    graphics::lines(x_vals, results$heat_wavelet, col = 3, lwd = 2)
    graphics::legend("topright", legend = c("Mexican Hat", "Meyer", "Heat"), col = 1:3, lwd = 2)
  }
  
  return(results)
} 