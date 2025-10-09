#' Calculate Graph Laplacian Matrix
#'
#' @description Compute unnormalized, normalized, or random-walk Laplacian from an adjacency matrix.
#'
#' @param W A square adjacency matrix (can be dense or sparse).
#' @param type Type of Laplacian to compute: "unnormalized", "normalized", or "randomwalk".
#'
#' @return Laplacian matrix of the same class as input.
#' @export
#'
#' @examples
#' \dontrun{
#' W <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
#' cal_laplacian(W, type = "normalized")
#' }
cal_laplacian <- function(W, type = c("unnormalized", "normalized", "randomwalk")) {
  type <- match.arg(type)

  is_sparse <- methods::is(W, "sparseMatrix")
  n <- nrow(W)

  # Compute row degrees
  deg <- if (is_sparse) Matrix::rowSums(W) else rowSums(W)
  deg[deg == 0] <- 1  # Avoid division by zero

  # Degree matrices
  D <- if (is_sparse) Matrix::Diagonal(n, deg) else diag(deg)
  D_half_inv <- if (is_sparse) Matrix::Diagonal(n, 1 / sqrt(deg)) else diag(1 / sqrt(deg))
  D_inv <- if (is_sparse) Matrix::Diagonal(n, 1 / deg) else diag(1 / deg)

  # Compute Laplacian
  L <- switch(type,
              unnormalized = D - W,
              normalized = diag(n) - D_half_inv %*% W %*% D_half_inv,
              randomwalk = diag(n) - D_inv %*% W)

  return(L)
}


#' Fast eigendecomposition of Laplacian matrix
#'
#' @description Perform fast eigendecomposition using RSpectra for large matrices
#'
#' @param laplacianMat Laplacian matrix
#' @param k_eigen Number of eigenvalues to compute (default: 25)
#' @param which Which eigenvalues to compute ("LM", "SM", etc.)
#' @param sigma Shift parameter for eigenvalue computation
#' @param opts Additional options for eigenvalue computation
#' @param lower Whether to compute from lower end of spectrum
#' @param ... Additional arguments
#'
#' @return List with eigenvalues (evalues) and eigenvectors (evectors)
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a Laplacian matrix and decompose
#' L <- matrix(c(2, -1, -1, -1, 2, -1, -1, -1, 2), nrow = 3)
#' decomp <- FastDecompositionLap(L, k_eigen = 25)
#' }
FastDecompositionLap <- function(laplacianMat = NULL, k_eigen = 25, which = "LM", sigma = NULL, opts = list(),
                                 lower = TRUE, ...) {
  res_decom <- RSpectra::eigs_sym(laplacianMat, k = k_eigen, which = which, sigma = sigma, opts = opts,
                        lower = lower)
  return(list(evalues = rev(res_decom$values),
              evectors = res_decom$vectors[, rev(seq_len(ncol(res_decom$vectors)))]))
}

#' Graph Fourier Transform
#'
#' @description Compute the Graph Fourier Transform (GFT) of a signal using Laplacian eigenvectors.
#'
#' @param signal Input signal (vector or matrix)
#' @param U Matrix of eigenvectors (dense matrix preferred)
#'
#' @return Transformed signal in the spectral domain (vector or matrix)
#' @export
gft <- function(signal, U) {
  # Convert signal to column matrix if it's a vector
  if (is.vector(signal)) {
    signal <- matrix(signal, ncol = 1)
  }

  # Ensure U is a matrix (convert if it's data.frame or other object)
  U <- as.matrix(U)

  # Compute GFT: transpose(U) %*% signal
  return(t(U) %*% signal)
}

#' Inverse Graph Fourier Transform
#'
#' @description Compute the Inverse Graph Fourier Transform (IGFT) of spectral coefficients using Laplacian eigenvectors.
#'
#' @param fourier_coeffs Input Fourier coefficients (vector or matrix)
#' @param U Matrix of eigenvectors (dense matrix preferred)
#'
#' @return Reconstructed signal in the vertex domain (vector or matrix)
#' @export
#'
#' @examples
#' \dontrun{
#' # Single signal
#' signal_reconstructed <- igft(fourier_coeffs, eigenvectors)
#' 
#' # Multiple signals (batch processing)
#' signals_reconstructed <- igft(fourier_coeffs_matrix, eigenvectors)
#' }
igft <- function(fourier_coeffs, U) {
  # Check if input was originally a vector
  was_vector <- is.vector(fourier_coeffs)
  
  # Convert fourier_coeffs to column matrix if it's a vector
  if (was_vector) {
    fourier_coeffs <- matrix(fourier_coeffs, ncol = 1)
  }

  # Ensure U is a matrix (convert if it's data.frame or other object)
  U <- as.matrix(U)

  # Compute IGFT: U %*% fourier_coeffs
  result <- U %*% fourier_coeffs
  
  # Return as vector if input was originally a vector (single signal)
  if (was_vector && ncol(result) == 1) {
    return(as.vector(result))
  }
  
  return(result)
}


#' Calculate cosine similarity between two vectors
#'
#' @description Calculate cosine similarity between two numeric vectors with numerical stability
#'
#' @param x First vector
#' @param y Second vector
#' @param eps Small numeric for numerical stability when norms are near zero (default 1e-12)
#'
#' @return Cosine similarity value (between -1 and 1)
#' @export
#'
#' @examples
#' x <- c(1, 2, 3)
#' y <- c(2, 3, 4)
#' similarity <- cosine_similarity(x, y)
#' # With custom eps for numerical stability
#' similarity2 <- cosine_similarity(x, y, eps = 1e-10)
cosine_similarity <- function(x, y, eps = 1e-12) {
  # Convert to numeric vectors
  x <- as.numeric(x)
  y <- as.numeric(y)

  # Check for same length
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }

  # Calculate magnitudes
  magnitude_x <- sqrt(sum(x^2))
  magnitude_y <- sqrt(sum(y^2))

  # Handle zero magnitude cases with numerical stability
  if (magnitude_x < eps && magnitude_y < eps) {
    return(0)
  }
  if (magnitude_x < eps || magnitude_y < eps) {
    return(0)
  }

  # Calculate dot product and cosine similarity
  dot_product <- sum(x * y)
  cosine_sim <- dot_product / max(magnitude_x * magnitude_y, eps)
  
  # Clamp to [-1, 1] range for numerical stability
  cosine_sim <- max(-1, min(1, cosine_sim))

  return(cosine_sim)
}


#' Find knee point in a curve
#'
#' @description Simple knee point detection using the maximum curvature method
#'
#' @param y Numeric vector of y values
#' @param sensitivity Sensitivity parameter (not used in this simple implementation)
#'
#' @return Index of the knee point
#' @export
#'
#' @examples
#' y <- c(1, 2, 3, 10, 11, 12)  # curve with a knee
#' knee_idx <- find_knee_point(y)
find_knee_point <- function(y, sensitivity = 1) {
  if (length(y) < 3) return(1)
  
  # Normalize the data to [0, 1]
  x <- seq_along(y)
  y_norm <- (y - min(y)) / (max(y) - min(y))
  x_norm <- (x - min(x)) / (max(x) - min(x))
  
  # Calculate distances from the line connecting first and last points
  distances <- rep(0, length(y))
  
  for (i in 2:(length(y) - 1)) {
    # Distance from point to line connecting first and last points
    x1 <- x_norm[1]
    y1 <- y_norm[1]
    x2 <- x_norm[length(x_norm)]
    y2 <- y_norm[length(y_norm)]
    xi <- x_norm[i]
    yi <- y_norm[i]
    
    # Point-to-line distance formula
    distances[i] <- abs((y2 - y1) * xi - (x2 - x1) * yi + x2 * y1 - y2 * x1) / 
                    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  }
  
  # Return the index with maximum distance (knee point)
  return(which.max(distances))
}



#' Install and load packages
#'
#' @description Utility function to install and load packages from CRAN or GitHub
#'
#' @param packages Named vector where names are package names and values are source URLs
#'
#' @return NULL (side effect: installs and loads packages)
#' @export
#' @importFrom utils install.packages
#'
#' @examples
#' \dontrun{
#' packages <- c("ggplot2" = "ggplot2", "devtools" = "r-lib/devtools")
#' install_and_load(packages)
#' }
install_and_load <- function(packages) {
  # Ensure that the 'remotes' package is available
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  # Loop through each package and install from the appropriate source
  for (pkg in names(packages)) {
    package_name <- pkg
    source_url <- packages[pkg]

    # Check if the package is already installed
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      # Install from GitHub if URL is provided
      if (grepl("github", source_url)) {
        cat(sprintf("Installing %s from GitHub repository %s\n", package_name, source_url))
        remotes::install_github(source_url)
      } else {  # Install from CRAN
        cat(sprintf("Installing package: %s from CRAN\n", package_name))
        install.packages(package_name)
      }

      # Load the package
      library(package_name, character.only = TRUE)
    } else {
      cat(sprintf("Package already installed: %s\n", package_name))
    }
  }
}

#' Hello function for SGWT package demonstration
#'
#' @description Simple hello function to demonstrate package loading
#'
#' @return Character string with greeting
#' @export
#'
#' @examples
#' hello_sgwt()
hello_sgwt <- function() {
  return("Hello from SGWT package! Ready for Spectral Graph Wavelet Transform analysis.")
}

# Private helper function: Extract decomposition from SGWT result
.get_decomp <- function(x) {
  if (!is.null(x$decomposition)) x$decomposition else x
}

#' Check K-band limited property of signals
#'
#' @description Analyze whether signals are k-band limited by comparing low-frequency 
#' and high-frequency Fourier coefficients using eigendecomposition and statistical testing.
#' Builds graph and computes Laplacian directly from SGWT data.
#'
#' @param SG SGWT object with Data slot (from initSGWT)
#' @param signals Character vector of signal names to analyze. If NULL, uses all signals from SG$Data$signals
#' @param alpha Significance level for Wilcoxon test (default: 0.05)
#' @param verbose Logical; if TRUE, print progress messages (default: TRUE)
#' @param k Number of nearest neighbors for graph construction (default: 25)
#' @param laplacian_type Type of Laplacian ("unnormalized", "normalized", or "randomwalk") (default: "normalized")
#' @return List containing:
#'   \describe{
#'     \item{is_kband_limited}{Logical; TRUE if all signals are k-band limited}
#'     \item{knee_point_low}{Integer; knee point index for low-frequency eigenvalues}
#'     \item{knee_point_high}{Integer; knee point index for high-frequency eigenvalues}
#'     \item{signal_results}{List with per-signal test results including p-values and Fourier coefficients}
#'   }
#' @export
#' @importFrom stats median wilcox.test
#' @importFrom igraph graph_from_edgelist as_adjacency_matrix
#'
#' @examples
#' \dontrun{
#' # Initialize SGWT object (no need to run runSpecGraph)
#' SG <- initSGWT(data, signals = c("signal1", "signal2"))
#' 
#' # Check k-band limited property
#' result <- checkKband(SG, signals = c("signal1", "signal2"), k = 30)
#' if (result$is_kband_limited) {
#'   cat("All signals are k-band limited")
#' }
#' }
checkKband <- function(SG, signals = NULL, alpha = 0.05, verbose = TRUE, k = 25, laplacian_type = "normalized") {
  
  # Input validation
  if (!inherits(SG, "SGWT")) {
    stop("SG must be an SGWT object")
  }
  
  if (is.null(SG$Data)) {
    stop("Data slot not found in SGWT object")
  }
  
  # Use all signals if not specified
  if (is.null(signals)) {
    signals <- SG$Data$signals
  }
  
  # Validate signal names
  missing_signals <- setdiff(signals, SG$Data$signals)
  if (length(missing_signals) > 0) {
    stop(paste("Signals not found in SGWT object:", paste(missing_signals, collapse = ", ")))
  }
  
  # Extract coordinates and build graph directly
  data.in <- SG$Data$data
  x_col <- SG$Data$x_col
  y_col <- SG$Data$y_col
  
  coords <- data.in[, c(x_col, y_col)]
  n_nodes <- nrow(coords)
  k_eigen <- floor(4*sqrt(n_nodes))
  
  if (verbose) {
    cat("Analyzing k-band limited property for", length(signals), "signals\n")
    cat("Number of nodes:", n_nodes, "\n")
    cat("Building graph with k =", k, "neighbors...\n")
  }
  
  # Build k-nearest neighbor graph (same as runSpecGraph)
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("RANN package is required for k-nearest neighbor graph construction")
  }
  
  # Build k-nearest neighbor graph (unweighted connectivity)
  nn <- RANN::nn2(coords, k = k + 1)
  adj_list <- lapply(seq_len(n_nodes), function(i) setdiff(nn$nn.idx[i, ], i))
  edges <- do.call(rbind, lapply(seq_along(adj_list), function(i) cbind(i, adj_list[[i]])))
  edges <- unique(t(apply(edges, 1, sort)))
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  adjacency_matrix <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  
  # Compute Laplacian matrix
  if (verbose) cat("Computing Laplacian matrix...\n")
  laplacian_matrix <- cal_laplacian(adjacency_matrix, type = laplacian_type)
  
  if (verbose) {
    cat("Computing", k_eigen, "low and high frequency eigenvalues/eigenvectors\n")
  }
  
  # Compute low-frequency eigendecomposition (smallest eigenvalues)
  if (verbose) cat("Computing low-frequency eigendecomposition...\n")
  low_decomp <- FastDecompositionLap(
    laplacianMat = laplacian_matrix,
    k_eigen = k_eigen,
    which = "SM",
    lower = TRUE
  )
  
  # Compute high-frequency eigendecomposition (largest eigenvalues)  
  if (verbose) cat("Computing high-frequency eigendecomposition...\n")
  high_decomp <- FastDecompositionLap(
    laplacianMat = laplacian_matrix,
    k_eigen = k_eigen,
    which = "LM",
    lower = FALSE
  )
  
  # Find knee points
  if (verbose) cat("Finding knee points in eigenvalue spectra...\n")
  knee_point_low <- find_knee_point(low_decomp$evalues)
  knee_point_high <- find_knee_point(high_decomp$evalues)
  
  if (verbose) {
    cat("Low-frequency knee point at index:", knee_point_low, "\n")
    cat("High-frequency knee point at index:", knee_point_high, "\n")
  }
  
  # Select eigenvectors up to knee points
  low_freq_eigenvecs <- low_decomp$evectors[, 1:knee_point_low, drop = FALSE]
  high_freq_eigenvecs <- high_decomp$evectors[, 1:knee_point_high, drop = FALSE]
  
  # Concatenate low and high frequency eigenvectors
  combined_eigenvecs <- cbind(low_freq_eigenvecs, high_freq_eigenvecs)
  
  if (verbose) {
    cat("Selected", knee_point_low, "low-frequency and", knee_point_high, "high-frequency components\n")
  }
  
  # Analyze each signal
  signal_results <- list()
  p_values <- numeric(length(signals))
  names(p_values) <- signals
  all_kband_limited <- TRUE
  
  for (sig in signals) {
    if (verbose) cat("Analyzing signal:", sig, "\n")
    
    # Extract signal vector
    signal_vec <- as.numeric(SG$Data$data[[sig]])
    
    # Compute Fourier coefficients using GFT
    fourier_coeffs <- gft(signal_vec, combined_eigenvecs)
    
    # Split into low and high frequency components
    low_freq_fc <- fourier_coeffs[1:knee_point_low]
    high_freq_fc <- fourier_coeffs[(knee_point_low + 1):(knee_point_low + knee_point_high)]
    
    # Remove DC component from low-frequency FC (first component)
    if (length(low_freq_fc) > 1) {
      low_freq_fc_no_dc <- low_freq_fc[-1]
    } else {
      warning(paste("Signal", sig, "has only DC component in low-frequency band"))
      low_freq_fc_no_dc <- numeric(0)
    }
    
    # Perform Wilcoxon test (one-sided: low > high)
    if (length(low_freq_fc_no_dc) > 0 && length(high_freq_fc) > 0) {
      # Test if low-frequency FC magnitudes are significantly greater than high-frequency FC magnitudes
      low_freq_fc_no_dc_filtered <- abs(low_freq_fc_no_dc)[abs(low_freq_fc_no_dc) > mean(abs(low_freq_fc_no_dc))]
      high_freq_fc_filtered <- abs(high_freq_fc)[abs(high_freq_fc) > mean(abs(high_freq_fc))]
      # wilcox test with Wilcoxon rank sum test with continuity correction
      wilcox_result <- wilcox.test(
        low_freq_fc_no_dc_filtered, 
        high_freq_fc_filtered, 
        alternative = "greater"
      )
      # adjust the p-value for multiple testing (Bonferroni correction)
      p_value_adjusted <- wilcox_result$p.value *length(signals)
      is_significant <- p_value_adjusted < alpha
    } else {
      p_value_adjusted <- NA
      is_significant <- FALSE
      warning(paste("Cannot perform Wilcoxon test for signal", sig, "due to insufficient data"))
    }
    
    p_values[sig] <- p_value_adjusted
    
    signal_results[[sig]] <- list(
      low_freq_fc = low_freq_fc,
      high_freq_fc = high_freq_fc,
      low_freq_fc_no_dc = low_freq_fc_no_dc,
      p_value_adjusted = p_value_adjusted,
      is_kband_limited = is_significant
    )
    
    if (!is_significant) {
      all_kband_limited <- FALSE
    }
    
    if (verbose) {
      cat("  Wilcoxon test p-value:", round(p_value_adjusted, 6), "\n")
      cat("  K-band limited:", is_significant, "\n")
    }
  }
  
  if (verbose) {
    cat("\n=== Summary ===\n")
    cat("All signals k-band limited:", all_kband_limited, "\n")
  }
  
  return(list(
    is_kband_limited = all_kband_limited,
    knee_point_low = knee_point_low,
    knee_point_high = knee_point_high,
    signal_results = signal_results
  ))
}



