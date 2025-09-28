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
#' @param k_neighbor Number of neighbors for eigenvalue computation (default: 25)
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
#' decomp <- FastDecompositionLap(L, k_neighbor = 25)
#' }
FastDecompositionLap <- function(laplacianMat = NULL, k_neighbor = 25, which = "LM", sigma = NULL, opts = list(),
                                 lower = TRUE, ...) {
  res_decom <- RSpectra::eigs_sym(laplacianMat, k = k_neighbor, which = which, sigma = sigma, opts = opts,
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


