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
#' @param k_fold Multiplier for number of eigenvalues to compute (default: 1.5)
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
#' decomp <- FastDecompositionLap(L, k_fold = 2)
#' }
FastDecompositionLap <- function(laplacianMat = NULL, k_fold = 1.5, which = "LM", sigma = NULL, opts = list(),
                                 lower = TRUE, ...) {
  res_decom <- RSpectra::eigs_sym(laplacianMat, k = k_fold * sqrt(ncol(laplacianMat)), which = which, sigma = sigma, opts = opts,
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

#' Calculate eigenvalues and eigenvectors with knee detection
#'
#' @description Calculate eigenvalues and eigenvectors of a spatial graph with automatic
#' detection of the low-frequency cutoff using knee detection
#'
#' @param data.in Data frame with spatial coordinates
#' @param x_col Character string specifying the column name for X coordinates (default: "x")
#' @param y_col Character string specifying the column name for Y coordinates (default: "y")
#' @param k Number of nearest neighbors (default: 25)
#' @param k_fold Eigendecomposition parameter (default: 15)
#' @param sensitivity Sensitivity parameter for knee detection (default: 2)
#'
#' @return List containing knee point, eigenvectors, and eigenvalues
#' @export
#' @importFrom graphics abline plot
#'
#' @examples
#' \dontrun{
#' # Create spatial data
#' data <- data.frame(x = runif(100), y = runif(100))
#' result <- Cal_Eigen(data, k = 10)
#' 
#' # With custom column names
#' data2 <- data.frame(X = runif(100), Y = runif(100))
#' result2 <- Cal_Eigen(data2, x_col = "X", y_col = "Y", k = 10)
#' }
Cal_Eigen <- function(data.in = NULL, x_col = "x", y_col = "y", k = 25, k_fold = 15, sensitivity = 2){
  data.in <- data.in
  k_fold <- k_fold # how many eigen values and vectors need from low frequency
  sensitivity <- sensitivity
  k <- k  # Number of nearest neighbors
  # Validate column names
  if (!all(c(x_col, y_col) %in% colnames(data.in))) {
    stop(paste("Data must contain columns:", x_col, "and", y_col, "for spatial coordinates"))
  }
  
  nn <- RANN::nn2(data.in[, c(x_col, y_col)], k = k + 1)  # Include the point itself in neighbors

  adj_list <- lapply(seq_len(nrow(data.in)), function(i) setdiff(nn$nn.idx[i, ], i))  # Remove self-loops
  # Convert adjacency list to edge list
  edges <- do.call(rbind, lapply(seq_along(adj_list), function(i) cbind(i, adj_list[[i]])))
  edges <- unique(t(apply(edges, 1, sort)))  # Remove duplicate edges
  # Create the graph
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  # Ensure vertex names are correctly set
  igraph::V(g)$name <- 1:igraph::vcount(g)

  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  #
  # Create the graph
  # Create a Graph object
  L <- cal_laplacian(A, "normalized")
  #
  L_decompose_res <- FastDecompositionLap(L, k_fold = k_fold, which = "SM")
  # calculate eigen value and vector
  eigenvalue <- L_decompose_res$evalues
  eigenvector <- L_decompose_res$evectors
  # cut low frequency using simple knee detection
  knee <- find_knee_point(eigenvalue, sensitivity = sensitivity)
  print(paste0("cutoff of low-frequency FM: ", knee))
  plot(eigenvalue, type = "o", col = "blue", xlab = "Index", ylab = "Value", main = "Line Plot of a Numeric Vector")
  abline(v = knee, col = "red", lwd = 2, lty = 2)
  return(list(knee, eigenvector, eigenvalue))
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

#' Calculate Graph Cross-Correlation (GCC) - DEPRECATED
#'
#' @description \strong{DEPRECATED:} This function is deprecated. Use \code{\link{sgwt_similarity}} instead 
#' for more comprehensive signal similarity analysis with energy normalization and advanced features.
#' 
#' Calculate Graph Cross-Correlation between two signals using Graph Fourier Transform.
#' This is a simplified approach that only considers low-frequency GFT components.
#'
#' @param data.in Data frame containing the signals
#' @param knee Knee point for frequency cutoff
#' @param signal1 Name of first signal column
#' @param signal2 Name of second signal column
#' @param eigenvector Matrix of eigenvectors
#'
#' @return Cosine similarity value
#' @export
#'
#' @examples
#' \dontrun{
#' # DEPRECATED - use sgwt_similarity instead
#' # gcc_value <- Cal_GCC(data, knee = 10, signal1 = "sig1", signal2 = "sig2", eigenvector = U)
#' 
#' # NEW RECOMMENDED APPROACH:
#' # sgwt1 <- SGWT(data, signal = "signal1", k = 25, J = 4)
#' # sgwt2 <- SGWT(data, signal = "signal2", k = 25, J = 4)  
#' # similarity <- sgwt_similarity(sgwt1, sgwt2)
#' }
Cal_GCC <- function(data.in = NULL, knee = NULL, signal1 = NULL, signal2 = NULL, eigenvector = NULL){
  .Deprecated("sgwt_similarity", package = "BioGSP", 
              msg = "Cal_GCC is deprecated. Use sgwt_similarity() for more comprehensive signal similarity analysis.")
  
  data.in <- as.data.frame(data.in)
  knee <- knee
  eigenvector <- eigenvector
  #
  # Apply GFT to both signals
  gft_signal1 <- gft(as.numeric(data.in[, signal1]), eigenvector)
  gft_signal2 <- gft(as.numeric(data.in[, signal2]), eigenvector)
  cosine_sim <- cosine_similarity(gft_signal1[2:knee[1]], gft_signal2[2:knee[1]])
  return(cosine_sim)
}

#' Comprehensive Signal Similarity Analysis
#'
#' @description Compute similarity between two signals using either raw signals (via SGWT) 
#' or pre-computed SGWT decompositions. This function provides a unified interface for 
#' signal similarity analysis with energy normalization and advanced features.
#'
#' @param signal1 Either a character string (column name in data.in), numeric vector, or SGWT result object
#' @param signal2 Either a character string (column name in data.in), numeric vector, or SGWT result object  
#' @param data.in Data frame containing signals (required if signal1/signal2 are column names)
#' @param x_col Column name for X coordinates (default: "x")
#' @param y_col Column name for Y coordinates (default: "y") 
#' @param k Number of nearest neighbors for graph construction (default: 25)
#' @param J Number of wavelet scales (default: 4)
#' @param kernel_type Wavelet kernel type (default: "mexican_hat")
#' @param eps Numerical stability parameter (default: 1e-12)
#' @param low_only If TRUE, use only low-frequency similarity (default: FALSE)
#' @param return_parts If TRUE, return detailed components; if FALSE, return scalar similarity (default: FALSE)
#' @param ... Additional arguments passed to SGWT()
#'
#' @return Similarity score (scalar) or detailed similarity analysis (list) depending on return_parts
#' @export
#'
#' @examples
#' \dontrun{
#' # Method 1: Direct signals in data frame
#' data <- data.frame(x = runif(100), y = runif(100), 
#'                   signal1 = rnorm(100), signal2 = rnorm(100))
#' sim1 <- sgwt_similarity("signal1", "signal2", data.in = data)
#' 
#' # Method 2: Pre-computed SGWT results
#' sgwt1 <- SGWT(data, signal = "signal1", k = 25, J = 4)
#' sgwt2 <- SGWT(data, signal = "signal2", k = 25, J = 4)
#' sim2 <- sgwt_similarity(sgwt1, sgwt2)
#' 
#' # Method 3: Mixed - one SGWT result, one raw signal
#' sim3 <- sgwt_similarity(sgwt1, "signal2", data.in = data)
#' }
sgwt_similarity <- function(signal1, signal2, data.in = NULL, 
                           x_col = "x", y_col = "y", k = 25, J = 4, 
                           kernel_type = "mexican_hat", eps = 1e-12, 
                           low_only = FALSE, return_parts = FALSE, ...) {
  
  # Helper function to process input signals
  .process_signal <- function(signal, data.in, x_col, y_col, k, J, kernel_type, ...) {
    if (is.character(signal) && length(signal) == 1) {
      # Signal is a column name - need to compute SGWT
      if (is.null(data.in)) {
        stop("data.in must be provided when signal is specified as column name")
      }
      return(SGWT(data.in = data.in, x_col = x_col, y_col = y_col, 
                  signal = signal, k = k, J = J, kernel_type = kernel_type, 
                  return_all = TRUE, ...))
    } else if (is.numeric(signal)) {
      # Signal is a numeric vector - need data.in for spatial coordinates
      if (is.null(data.in)) {
        stop("data.in must be provided when signal is a numeric vector")
      }
      temp_data <- data.in
      temp_data$temp_signal <- signal
      return(SGWT(data.in = temp_data, x_col = x_col, y_col = y_col,
                  signal = "temp_signal", k = k, J = J, kernel_type = kernel_type,
                  return_all = TRUE, ...))
    } else {
      # Assume it's already an SGWT result
      return(signal)
    }
  }
  
  # Process both signals
  sgwt1 <- .process_signal(signal1, data.in, x_col, y_col, k, J, kernel_type, ...)
  sgwt2 <- .process_signal(signal2, data.in, x_col, y_col, k, J, kernel_type, ...)
  
  # Use the comprehensive weighted similarity function
  return(sgwt_weighted_similarity(sgwt1, sgwt2, eps = eps, 
                                 low_only = low_only, return_parts = return_parts))
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

# Private helper function: Return matrix n×J of wavelet coefficients ordered by scale index
.wavelet_matrix <- function(decomp) {
  nm <- names(decomp$coefficients)
  idx <- grep("^wavelet_scale_", nm)
  if (length(idx) == 0L) stop("No wavelet coefficients found (expected names starting with 'wavelet_scale_').")
  # Order by numeric suffix
  ord <- order(as.integer(sub("^wavelet_scale_", "", nm[idx])))
  mats <- lapply(nm[idx][ord], function(k) as.numeric(decomp$coefficients[[k]]))
  W <- do.call(cbind, mats) # n × J (column = scale)
  # Replace non-finite with 0
  W[!is.finite(W)] <- 0
  return(W)
}


#' Energy-normalized weighted similarity between two SGWT results
#'
#' @description Compute low-frequency cosine similarity (scaling), non-low cosine similarity
#' (flattened wavelet coefficients), and an overall energy-normalized weighted score.
#' If `low_only = TRUE`, compute only the low-frequency cosine and set `S = c_low`.
#'
#' @param sgwt_a SGWT output for signal A. Either the full list returned by SGWT (with `$decomposition`)
#'   or a decomposition list as returned by `sgwt_forward()`.
#' @param sgwt_b SGWT output for signal B. Same format as `sgwt_a`.
#' @param eps Small numeric for numerical stability when norms are near zero (default 1e-12).
#' @param validate Logical; if TRUE, check consistency of dimensions, scale count/order, and kernel family.
#' @param return_parts Logical; if TRUE (default), return a list with components; if FALSE, return only the scalar S.
#' @param low_only Logical; if TRUE, compute **low-frequency-only** similarity (skip non-low and set `S = c_low`).
#'
#' @return If `return_parts=TRUE`, a list with:
#'   * `c_low`     — cosine on scaling coefficients
#'   * `c_nonlow`  — cosine on flattened wavelet coefficients (**NA if `low_only = TRUE`**)
#'   * `w_low`     — macro weight for the low-frequency part (in [0,1])
#'   * `w_NL`      — 1 - w_low (non-low weight)
#'   * `S`         — final weighted similarity in [-1,1]
#'   * `E_low_a`, `E_NL_a`, `E_low_b`, `E_NL_b` — energy diagnostics per signal (**E_NL_* = NA if `low_only`**)
#'   * `n`, `J`    — length of signal and number of wavelet scales (**J = NA if `low_only`**)
#' If `return_parts=FALSE`, returns the scalar `S`.
#' @export
#' @examples
#' \dontrun{
#' # Assume two SGWT results sgwt_res1 and sgwt_res2 from SGWT(..., return_all=TRUE)
#' sim <- sgwt_weighted_similarity(sgwt_res1, sgwt_res2)
#' sim_low <- sgwt_weighted_similarity(sgwt_res1, sgwt_res2, low_only = TRUE)
#' str(sim)
#' }
sgwt_weighted_similarity <- function(sgwt_a, sgwt_b, eps = 1e-12, validate = TRUE, return_parts = TRUE, low_only = FALSE) {
  
  # 1. Accept both wrapper or decomposition
  A <- .get_decomp(sgwt_a)
  B <- .get_decomp(sgwt_b)
  
  # 2. Extract scaling vectors
  y_low_a <- as.numeric(A$coefficients$scaling)
  y_low_b <- as.numeric(B$coefficients$scaling)
  
  # Replace non-finite entries with 0 and warn
  if (any(!is.finite(y_low_a))) {
    warning("Non-finite values found in scaling coefficients of sgwt_a, replacing with 0")
    y_low_a[!is.finite(y_low_a)] <- 0
  }
  if (any(!is.finite(y_low_b))) {
    warning("Non-finite values found in scaling coefficients of sgwt_b, replacing with 0")
    y_low_b[!is.finite(y_low_b)] <- 0
  }
  
  # Energies
  E_low_a <- sum(y_low_a^2)
  E_low_b <- sum(y_low_b^2)
  
  # 3. Compute low-frequency cosine similarity
  c_low <- cosine_similarity(y_low_a, y_low_b, eps)
  
  # Short-circuit for low-only
  if (isTRUE(low_only)) {
    return(if (isTRUE(return_parts)) list(
      c_low = c_low, c_nonlow = NA_real_, w_low = 1.0, w_NL = 0.0, S = c_low,
      E_low_a = E_low_a, E_NL_a = NA_real_, E_low_b = E_low_b, E_NL_b = NA_real_,
      n = length(y_low_a), J = NA_integer_
    ) else c_low)
  }
  
  # 4. Collect wavelet coefficients
  W_a <- .wavelet_matrix(A)
  W_b <- .wavelet_matrix(B)
  
  # Flatten column-major
  v_a <- as.vector(W_a)
  v_b <- as.vector(W_b)
  
  # Energies
  E_NL_a <- sum(v_a^2)
  E_NL_b <- sum(v_b^2)
  
  # Number of scales and signal length
  J <- ncol(W_a)
  n <- nrow(W_a)
  
  # 5. Compute non-low cosine similarity
  c_nonlow <- cosine_similarity(v_a, v_b, eps)
  
  # 6. Macro weights (energy normalization)
  w_low_a <- E_low_a / (E_low_a + E_NL_a + eps)
  w_low_b <- E_low_b / (E_low_b + E_NL_b + eps)
  w_low <- pmax(0, pmin(1, 0.5 * (w_low_a + w_low_b)))
  w_NL <- 1 - w_low
  S <- w_low * c_low + w_NL * c_nonlow
  
  # 7. Validation when validate=TRUE
  if (isTRUE(validate)) {
    # Check scaling coefficient lengths match
    if (length(y_low_a) != length(y_low_b)) {
      stop("Scaling coefficients must have the same length between sgwt_a and sgwt_b")
    }
    
    # Check wavelet matrix dimensions match
    if (nrow(W_a) != nrow(W_b)) {
      stop("Signal lengths must match between sgwt_a and sgwt_b")
    }
    if (ncol(W_a) != ncol(W_b)) {
      stop("Number of wavelet scales must match between sgwt_a and sgwt_b")
    }
    
    # Extract and compare scales if present
    scales_a <- A$scales
    scales_b <- B$scales
    if (!is.null(scales_a) && !is.null(scales_b)) {
      if (length(scales_a) != length(scales_b) || !all.equal(scales_a, scales_b)) {
        warning("Scales differ between sgwt_a and sgwt_b")
      }
    } else if (!is.null(scales_a) || !is.null(scales_b)) {
      warning("Scales present in only one of sgwt_a or sgwt_b")
    }
    
    # Check that each has ≥1 wavelet scale
    if (J < 1) {
      stop("Each SGWT decomposition must have at least 1 wavelet scale")
    }
  }
  
  # 8. Return structured list or scalar S per return_parts
  if (isTRUE(return_parts)) {
    return(list(
      c_low = c_low,
      c_nonlow = c_nonlow,
      w_low = w_low,
      w_NL = w_NL,
      S = S,
      E_low_a = E_low_a,
      E_NL_a = E_NL_a,
      E_low_b = E_low_b,
      E_NL_b = E_NL_b,
      n = n,
      J = J
    ))
  } else {
    return(S)
  }
}
