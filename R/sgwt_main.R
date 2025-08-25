#' Spectral Graph Wavelet Transform (SGWT)
#'
#' @description Main function for performing Spectral Graph Wavelet Transform analysis
#' on spatial data. Provides a comprehensive interface for multi-scale analysis of
#' spatial signals using graph wavelets.
#'
#' @param data.in Data frame containing spatial coordinates and signal data.
#'   Must contain columns 'x' and 'y' for spatial coordinates.
#' @param signal Character string specifying the column name of the signal to analyze,
#'   or a numeric vector of signal values.
#' @param k Number of nearest neighbors for graph construction (default: 25)
#' @param scales Vector of scales for the wavelets. If NULL, scales are auto-generated.
#' @param J Number of scales to generate if scales is NULL (default: 5)
#' @param scaling_factor Scaling factor between consecutive scales (default: 2)
#' @param kernel_type Kernel family ("mexican_hat", "meyer", or "heat") that defines both scaling and wavelet filters
#' @param laplacian_type Type of graph Laplacian ("unnormalized", "normalized", or "randomwalk", default: "normalized")
#' @param k_fold Parameter for eigendecomposition (default: 15)
#' @param return_all Whether to return all analysis results (default: TRUE)
#'
#' @return If return_all = TRUE, returns a list containing:
#' \describe{
#'   \item{decomposition}{SGWT decomposition results}
#'   \item{reconstructed_signal}{Reconstructed signal for validation}
#'   \item{reconstruction_error}{RMSE between original and reconstructed signal}
#'   \item{original_signal}{Original input signal}
#'   \item{graph_info}{Graph construction information (adjacency matrix, Laplacian, eigenvalues, eigenvectors)}
#'   \item{data}{Original input data}
#'   \item{parameters}{Analysis parameters used}
#' }
#' If return_all = FALSE, returns only the SGWT decomposition.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate synthetic spatial data
#' set.seed(123)
#' n_points <- 100
#' x_coords <- rep(1:10, each = 10) + rnorm(n_points, 0, 0.1)
#' y_coords <- rep(1:10, times = 10) + rnorm(n_points, 0, 0.1)
#' signal_data <- sin(0.5 * x_coords) * cos(0.3 * y_coords) + rnorm(n_points, 0, 0.1)
#'
#' demo_data <- data.frame(x = x_coords, y = y_coords, signal = signal_data)
#'
#' # Apply SGWT
#' result <- SGWT(data.in = demo_data, signal = "signal", k = 8, J = 4)
#'
#' # View reconstruction error
#' print(result$reconstruction_error)
#' }
#'
#' @references
#' Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011).
#' Wavelets on graphs via spectral graph theory.
#' Applied and Computational Harmonic Analysis, 30(2), 129-150.
SGWT <- function(data.in = NULL, signal = NULL, k = 25,
                 scales = NULL, J = 5, scaling_factor = 2,
                 kernel_type = "mexican_hat",
                 laplacian_type = "normalized",
                 k_fold = 15, return_all = TRUE) {

  # Input validation
  if (is.null(data.in) || is.null(signal)) {
    stop("Both data.in and signal must be provided")
  }

  # Build graph from spatial coordinates
  if (all(c("x", "y") %in% colnames(data.in))) {
    cat("Building graph from spatial coordinates...\n")
    nn <- RANN::nn2(data.in[, c("x", "y")], k = k + 1)
    adj_list <- lapply(seq_len(nrow(data.in)), function(i) setdiff(nn$nn.idx[i, ], i))
    edges <- do.call(rbind, lapply(seq_along(adj_list), function(i) cbind(i, adj_list[[i]])))
    edges <- unique(t(apply(edges, 1, sort)))
    g <- igraph::graph_from_edgelist(edges, directed = FALSE)
    A <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  } else {
    stop("Data must contain 'x' and 'y' columns for spatial coordinates")
  }

  # Compute Laplacian and eigendecomposition
  cat("Computing Laplacian and eigendecomposition...\n")
  L <- cal_laplacian(A, laplacian_type)
  L_decompose_res <- FastDecompositionLap(L, k_fold = k_fold, which = "SM")
  eigenvalues <- L_decompose_res$evalues
  eigenvectors <- L_decompose_res$evectors

  # Auto-generate scales if not provided
  if (is.null(scales)) {
    lmax <- max(eigenvalues) * 0.95
    scales <- sgwt_auto_scales(lmax, J, scaling_factor)
    cat(paste("Auto-generated scales:", paste(round(scales, 4), collapse = ", "), "\n"))
  }

  # Extract signal
  if (is.character(signal)) {
    signal_data <- as.numeric(data.in[, signal])
  } else {
    signal_data <- as.numeric(signal)
  }

  # Perform SGWT decomposition
  eigenvectors <- as.matrix(eigenvectors)
  cat("Performing SGWT decomposition...\n")
  sgwt_result <- sgwt_forward(signal = signal_data,
                              eigenvectors = eigenvectors, eigenvalues = eigenvalues,
                              scales = scales,
                              kernel_type = kernel_type)

  # Reconstruction for validation
  reconstructed <- sgwt_inverse(sgwt_result)
  reconstruction_error <- sqrt(mean((signal_data - reconstructed)^2))

  cat(paste("Reconstruction RMSE:", round(reconstruction_error, 6), "\n"))

  if (return_all) {
    return(list(
      decomposition = sgwt_result,
      reconstructed_signal = reconstructed,
      reconstruction_error = reconstruction_error,
      original_signal = signal_data,
      graph_info = list(
        adjacency_matrix = A,
        laplacian_matrix = L,
        eigenvalues = eigenvalues,
        eigenvectors = eigenvectors
      ),
      data = data.in,
      parameters = list(
        k = k,
        scales = scales,
        J = J,
        kernel_type = kernel_type,
        laplacian_type = laplacian_type
      )
    ))
  } else {
    return(sgwt_result)
  }
}
