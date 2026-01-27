

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
sgwt_get_kernels <- function(kernel_type = "heat") {
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
compute_sgwt_filters <- function(eigenvalues, scales, lmax = NULL, kernel_type = "heat") {
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

#' Forward SGWT transform (single or batch)
#'
#' @description Transform signal(s) to spectral domain and apply SGWT filters.
#' Handles both single signals (vector) and multiple signals (matrix) efficiently.
#' Stores original and filtered Fourier coefficients for analysis.
#'
#' @param signal Input signal vector OR matrix where each column is a signal (n_vertices x n_signals)
#' @param eigenvectors Eigenvectors of the graph Laplacian
#' @param eigenvalues Eigenvalues of the graph Laplacian
#' @param scales Vector of scales for the wavelets
#' @param lmax Maximum eigenvalue (optional)
#' @param kernel_type Kernel family that defines both scaling and wavelet filters (default: "heat")
#'
#' @return List containing:
#'   \describe{
#'     \item{fourier_coefficients}{List with original and filtered Fourier coefficients}
#'     \item{filters}{Filter bank used}
#'   }
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data and compute graph
#' data <- data.frame(x = runif(50), y = runif(50), signal = rnorm(50))
#' SG <- initSGWT(data, signals = "signal", J = 3)
#' SG <- runSpecGraph(SG, k = 10)
#' eigenvectors <- SG$Graph$eigenvectors
#' eigenvalues <- SG$Graph$eigenvalues
#' scales <- SG$Parameters$scales
#' 
#' # Single signal
#' signal <- data$signal
#' result <- sgwt_forward(signal, eigenvectors, eigenvalues, scales)
#' 
#' # Multiple signals (batch processing)
#' signals_matrix <- cbind(data$signal, data$signal * 2, data$signal * 0.5)
#' result <- sgwt_forward(signals_matrix, eigenvectors, eigenvalues, scales)
#' }
sgwt_forward <- function(signal, eigenvectors, eigenvalues, scales, lmax = NULL, kernel_type = "heat") {
  
  # Validate dimensions
  if (is.matrix(signal)) {
    if (nrow(eigenvectors) != nrow(signal)) {
      stop("Number of vertices in signals matrix must match eigenvectors")
    }
  } else {
    if (nrow(eigenvectors) != length(signal)) {
      stop("Number of vertices in signal vector must match eigenvectors")
    }
  }
  
  # Compute filters
  filters <- compute_sgwt_filters(eigenvalues, scales, lmax, kernel_type)
  
  # Transform signal(s) to spectral domain using GFT
  # gft handles both vectors and matrices automatically
  signal_hat <- gft(signal, eigenvectors)
  
  # Store original and filtered Fourier coefficients
  fourier_coefficients <- list(
    original = signal_hat,
    filtered = vector("list", length(filters))
  )
  
  # Apply filters and store filtered Fourier coefficients
  for (i in seq_along(filters)) {
    # Apply filter - works for both single signals and matrices
    filtered_spectrum <- signal_hat * as.vector(filters[[i]])
    fourier_coefficients$filtered[[i]] <- filtered_spectrum
  }
  
  names(fourier_coefficients$filtered) <- c("scaling", paste0("wavelet_scale_", seq_along(scales)))
  
  return(list(
    fourier_coefficients = fourier_coefficients,
    filters = filters
  ))
}


#' Inverse SGWT transform (single or batch)
#'
#' @description Reconstruct signal(s) from filtered Fourier coefficients using inverse GFT.
#' Handles both single signals and multiple signals efficiently.
#' Returns detailed inverse transform results including low-pass, band-pass approximations,
#' reconstructed signal(s), and reconstruction error(s).
#'
#' @param sgwt_decomp SGWT decomposition object from sgwt_forward
#' @param eigenvectors Eigenvectors of the graph Laplacian (for inverse GFT)
#' @param original_signal Original signal vector OR matrix (n_vertices x n_signals) for error calculation (optional)
#'
#' @return List containing:
#'   \describe{
#'     \item{vertex_approximations}{Named list with inverse-transformed signals in vertex domain:}
#'       \itemize{
#'         \item{\code{low_pass}: Low-pass (scaling) approximation}
#'         \item{\code{wavelet_1}, \code{wavelet_2}, etc.: Band-pass (wavelet) approximations by scale}
#'       }
#'     \item{reconstructed_signal}{Full reconstructed signal (vector or matrix)}
#'     \item{reconstruction_error}{RMSE (scalar for single signal, vector for multiple signals)}
#'   }
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data and perform forward transform
#' data <- data.frame(x = runif(50), y = runif(50), signal = rnorm(50))
#' SG <- initSGWT(data, signals = "signal", J = 3)
#' SG <- runSpecGraph(SG, k = 10)
#' eigenvectors <- SG$Graph$eigenvectors
#' eigenvalues <- SG$Graph$eigenvalues
#' scales <- SG$Parameters$scales
#' 
#' # Single signal - forward transform first
#' original_signal <- data$signal
#' sgwt_decomp <- sgwt_forward(original_signal, eigenvectors, eigenvalues, scales)
#' inverse_result <- sgwt_inverse(sgwt_decomp, eigenvectors, original_signal)
#' 
#' # Multiple signals (batch processing)
#' original_signals_matrix <- cbind(data$signal, data$signal * 2)
#' sgwt_decomp <- sgwt_forward(original_signals_matrix, eigenvectors, eigenvalues, scales)
#' inverse_result <- sgwt_inverse(sgwt_decomp, eigenvectors, original_signals_matrix)
#' }
sgwt_inverse <- function(sgwt_decomp, eigenvectors, original_signal = NULL) {
  filtered_fourier <- sgwt_decomp$fourier_coefficients$filtered
  
  # Perform inverse GFT on all filtered coefficients using igft function
  vertex_approximations <- vector("list", length(filtered_fourier))
  names(vertex_approximations) <- names(filtered_fourier)
  
  for (i in seq_along(filtered_fourier)) {
    # Use igft function for inverse transform - handles both single signals and batches
    vertex_approximations[[i]] <- igft(filtered_fourier[[i]], eigenvectors)
  }
  
  # Rename components for clarity
  names(vertex_approximations)[names(vertex_approximations) == "scaling"] <- "low_pass"
  wavelet_names <- names(vertex_approximations)[grep("^wavelet_scale_", names(vertex_approximations))]
  for (name in wavelet_names) {
    scale_num <- sub("^wavelet_scale_", "", name)
    names(vertex_approximations)[names(vertex_approximations) == name] <- paste0("wavelet_", scale_num)
  }
  
  # Reconstruct signal(s) - sum of all components
  reconstructed <- Reduce("+", vertex_approximations)
  
  # Calculate reconstruction error(s) if original signal provided
  reconstruction_error <- NULL
  if (!is.null(original_signal)) {
    # Validate dimensions
    if (is.matrix(original_signal) && is.matrix(reconstructed)) {
      if (ncol(original_signal) != ncol(reconstructed) || nrow(original_signal) != nrow(reconstructed)) {
        stop("Dimensions of original_signal must match reconstructed signal")
      }
      # Calculate RMSE for each signal (column)
      reconstruction_error <- sqrt(colMeans((original_signal - reconstructed)^2))
    } else if (is.vector(original_signal) && (is.vector(reconstructed) || ncol(reconstructed) == 1)) {
      reconstructed_vec <- as.vector(reconstructed)
      if (length(original_signal) != length(reconstructed_vec)) {
        stop("Length of original_signal must match reconstructed signal")
      }
      # Calculate RMSE for single signal
      reconstruction_error <- sqrt(mean((original_signal - reconstructed_vec)^2))
    } else {
      stop("Type mismatch between original_signal and reconstructed signal")
    }
  }
  
  return(list(
    vertex_approximations = vertex_approximations,
    reconstructed_signal = reconstructed,
    reconstruction_error = reconstruction_error
  ))
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
#' @importFrom graphics par plot lines legend
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
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
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