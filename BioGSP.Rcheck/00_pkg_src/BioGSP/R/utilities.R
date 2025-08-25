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
              evectors = res_decom$vectors[, ncol(res_decom$vectors):1]))
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
#' @description Calculate cosine similarity between two numeric vectors
#'
#' @param x First vector
#' @param y Second vector
#'
#' @return Cosine similarity value (between -1 and 1)
#' @export
#'
#' @examples
#' x <- c(1, 2, 3)
#' y <- c(2, 3, 4)
#' similarity <- cosine_similarity(x, y)
cosine_similarity <- function(x, y) {
  # Convert to numeric vectors
  x <- as.numeric(x)
  y <- as.numeric(y)

  # Check for same length
  if (length(x) != length(y)) {
    stop("Vectors must have the same length")
  }

  # Calculate dot product
  dot_product <- sum(x * y)

  # Calculate magnitudes
  magnitude_x <- sqrt(sum(x^2))
  magnitude_y <- sqrt(sum(y^2))

  # Handle zero magnitude case
  if (magnitude_x == 0 || magnitude_y == 0) {
    return(0)
  }

  # Calculate cosine similarity
  cosine_sim <- dot_product / (magnitude_x * magnitude_y)

  return(cosine_sim)
}

#' Calculate eigenvalues and eigenvectors with knee detection
#'
#' @description Calculate eigenvalues and eigenvectors of a spatial graph with automatic
#' detection of the low-frequency cutoff using knee detection
#'
#' @param data.in Data frame with x and y coordinates
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
#' }
Cal_Eigen <- function(data.in = NULL, k = 25, k_fold = 15, sensitivity = 2){
  data.in <- data.in
  k_fold <- k_fold # how many eigen values and vectors need from low frequency
  sensitivity <- sensitivity
  k <- k  # Number of nearest neighbors
  nn <- RANN::nn2(data.in[, c("x", "y")], k = k + 1)  # Include the point itself in neighbors

  adj_list <- lapply(1:nrow(data.in), function(i) setdiff(nn$nn.idx[i, ], i))  # Remove self-loops
  # Convert adjacency list to edge list
  edges <- do.call(rbind, lapply(1:length(adj_list), function(i) cbind(i, adj_list[[i]])))
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

#' Calculate Graph Cross-Correlation (GCC)
#'
#' @description Calculate Graph Cross-Correlation between two signals using Graph Fourier Transform
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
#' # Assuming you have data with two signals and eigenvectors
#' # gcc_value <- Cal_GCC(data, knee = 10, signal1 = "sig1", signal2 = "sig2", eigenvector = U)
#' }
Cal_GCC <- function(data.in = NULL, knee = NULL, signal1 = NULL, signal2 = NULL, eigenvector = NULL){
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
