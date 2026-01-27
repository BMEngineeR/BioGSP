#' Initialize SGWT object
#'
#' @description Build an SGWT object with Data and Parameters slots, validate inputs.
#' 
#' @param data.in Data frame containing spatial coordinates and signal data
#' @param x_col Character string specifying the column name for X coordinates (default: "x")
#' @param y_col Character string specifying the column name for Y coordinates (default: "y") 
#' @param signals Character vector of signal column names to analyze. If NULL, all non-coordinate columns are used.
#' @param scales Vector of scales for the wavelets. If NULL, scales are auto-generated.
#' @param J Number of scales to generate if scales is NULL (default: 5)
#' @param scaling_factor Scaling factor between consecutive scales (default: 2)
#' @param kernel_type Kernel family ("mexican_hat", "meyer", or "heat") (default: "heat")
#'
#' @return SGWT object with Data and Parameters slots initialized
#' @export
#'
#' @examples
#' \donttest{
#' # Initialize SGWT object
#' data <- data.frame(x = runif(100), y = runif(100), 
#'                   signal1 = rnorm(100), signal2 = rnorm(100))
#' SG <- initSGWT(data, signals = c("signal1", "signal2"))
#' }
initSGWT <- function(data.in, x_col = "x", y_col = "y", signals = NULL,
                     scales = NULL, J = 5, scaling_factor = 2,
                     kernel_type = "heat") {
  
  # Input validation
  if (is.null(data.in)) stop("data.in must be provided")
  if (!all(c(x_col, y_col) %in% colnames(data.in))) {
    stop(paste("Data must contain columns:", x_col, "and", y_col))
  }
  
  # Auto-detect signals if not provided
  if (is.null(signals)) {
    signals <- setdiff(colnames(data.in), c(x_col, y_col))
    if (length(signals) == 0) {
      stop("No signal columns found in data")
    }
  }
  
  # Validate signal columns exist
  missing_signals <- setdiff(signals, colnames(data.in))
  if (length(missing_signals) > 0) {
    stop(paste("Signal columns not found in data:", paste(missing_signals, collapse = ", ")))
  }
  
  # Create SGWT object
  SG <- list(
    Data = list(
      data = data.in,
      x_col = x_col,
      y_col = y_col,
      signals = signals
    ),
    Graph = NULL,
    Forward = NULL,
    Inverse = NULL,
    Parameters = list(
      scales = scales,
      J = J,
      scaling_factor = scaling_factor,
      kernel_type = kernel_type
    )
  )
  
  class(SG) <- "SGWT"
  return(SG)
}

#' Build spectral graph for SGWT object
#'
#' @description Generate Graph slot information including adjacency matrix, 
#' Laplacian matrix, eigenvalues, and eigenvectors.
#'
#' @param SG SGWT object from initSGWT()
#' @param k Number of nearest neighbors for graph construction (default: 25)
#' @param laplacian_type Type of graph Laplacian ("unnormalized", "normalized", or "randomwalk") (default: "normalized")
#' @param length_eigenvalue Number of eigenvalues/eigenvectors to compute (default: NULL, uses full length)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return Updated SGWT object with Graph slot populated
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data
#' data <- data.frame(x = runif(100), y = runif(100), signal = rnorm(100))
#' SG <- initSGWT(data, signals = "signal")
#' 
#' # Uses full length by default
#' SG <- runSpecGraph(SG, k = 30, laplacian_type = "normalized")  
#' 
#' # Or specify custom length
#' SG2 <- initSGWT(data, signals = "signal")
#' SG2 <- runSpecGraph(SG2, k = 30, laplacian_type = "normalized", 
#'                     length_eigenvalue = 30)  
#' }
runSpecGraph <- function(SG, k = 25, laplacian_type = "normalized", length_eigenvalue = NULL, verbose = TRUE) {
  
  # Validate input
  if (!inherits(SG, "SGWT")) {
    stop("Input must be an SGWT object from initSGWT()")
  }
  if (is.null(SG$Data)) {
    stop("SGWT object must have Data slot initialized")
  }
  
  # Extract data
  data.in <- SG$Data$data
  x_col <- SG$Data$x_col
  y_col <- SG$Data$y_col
  
  # Set default length_eigenvalue to full length if not specified
  if (is.null(length_eigenvalue)) {
    length_eigenvalue <- nrow(data.in)
  }
  
  if (verbose) cat("Building graph from spatial coordinates...\n")
  
  # Build k-nearest neighbor graph
  nn <- RANN::nn2(data.in[, c(x_col, y_col)], k = k + 1)
  adj_list <- lapply(seq_len(nrow(data.in)), function(i) setdiff(nn$nn.idx[i, ], i))
  edges <- do.call(rbind, lapply(seq_along(adj_list), function(i) cbind(i, adj_list[[i]])))
  edges <- unique(t(apply(edges, 1, sort)))
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  
  if (verbose) cat("Computing Laplacian and eigendecomposition...\n")
  
  # Compute Laplacian matrix
  L <- cal_laplacian(A, laplacian_type)
  
  # Eigendecomposition
  decomp <- FastDecompositionLap(L, k_eigen = length_eigenvalue, which = "SM")
  
  # Update SGWT object
  SG$Graph <- list(
    adjacency_matrix = A,
    laplacian_matrix = L,
    eigenvalues = decomp$evalues,
    eigenvectors = decomp$evectors
  )
  
  # Auto-generate scales if not provided (now that we have eigenvalues)
  if (is.null(SG$Parameters$scales)) {
    lmax <- max(decomp$evalues) * 0.95
    SG$Parameters$scales <- sgwt_auto_scales(lmax, SG$Parameters$J, SG$Parameters$scaling_factor)
    if (verbose) {
      cat(paste("Auto-generated scales:", paste(round(SG$Parameters$scales, 4), collapse = ", "), "\n"))
    }
  }
  
  if (verbose) cat("Graph construction completed.\n")
  
  return(SG)
}

#' Run SGWT forward and inverse transforms for all signals
#'
#' @description Perform SGWT analysis on all signals in the SGWT object.
#' Uses batch processing for multiple signals when possible for efficiency.
#' Assumes Graph slot is populated by runSpecGraph().
#'
#' @param SG SGWT object with Graph slot populated
#' @param use_batch Whether to use batch processing for multiple signals (default: TRUE)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return Updated SGWT object with Forward and Inverse slots populated
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data
#' data <- data.frame(x = runif(100), y = runif(100), signal = rnorm(100))
#' SG <- initSGWT(data, signals = "signal")
#' SG <- runSpecGraph(SG, k = 15)
#' 
#' # Uses batch processing by default
#' SG <- runSGWT(SG)
#' 
#' # Or force individual processing
#' SG2 <- initSGWT(data, signals = "signal")
#' SG2 <- runSpecGraph(SG2, k = 15)
#' SG2 <- runSGWT(SG2, use_batch = FALSE)
#' }
runSGWT <- function(SG, use_batch = TRUE, verbose = TRUE) {
  
  # Validate input
  if (!inherits(SG, "SGWT")) {
    stop("Input must be an SGWT object")
  }
  if (is.null(SG$Graph)) {
    stop("Graph slot is empty. Run runSpecGraph() first.")
  }
  
  # Extract components
  eigenvalues <- SG$Graph$eigenvalues
  eigenvectors <- as.matrix(SG$Graph$eigenvectors)
  params <- SG$Parameters
  signals <- SG$Data$signals
  data.in <- SG$Data$data
  
  # Scales should have been generated in runSpecGraph
  if (is.null(params$scales)) {
    stop("Scales not found. Make sure to run runSpecGraph() before runSGWT().")
  }
  
  n_signals <- length(signals)
  
  if (verbose) cat("Performing SGWT analysis for", n_signals, "signals...\n")
  
  if (use_batch && n_signals > 1) {
    # Batch processing for multiple signals
    if (verbose) cat("Using batch processing for efficiency...\n")
    
    # Create signal matrix (n_vertices x n_signals)
    signals_matrix <- as.matrix(data.in[, signals, drop = FALSE])
    
    # Batch forward transform
    batch_forward <- sgwt_forward(signals_matrix, eigenvectors, eigenvalues, params$scales, 
                                 kernel_type = params$kernel_type)
    
    # Batch inverse transform
    batch_inverse <- sgwt_inverse(batch_forward, eigenvectors, signals_matrix)
    
    # Split results back into individual signals
    forward_list <- list()
    inverse_list <- list()
    
    for (i in seq_along(signals)) {
      sig <- signals[i]
      
      # Extract individual signal results from batch
      forward_list[[sig]] <- list(
        fourier_coefficients = list(
          original = as.vector(batch_forward$fourier_coefficients$original[, i]),
          filtered = lapply(batch_forward$fourier_coefficients$filtered, function(x) as.vector(x[, i]))
        ),
        filters = batch_forward$filters
      )
      
      # Extract individual vertex approximations and create coefficients structure
      coefficients_individual <- list()
      for (comp_name in names(batch_inverse$vertex_approximations)) {
        if (is.matrix(batch_inverse$vertex_approximations[[comp_name]])) {
          coefficients_individual[[comp_name]] <- as.vector(batch_inverse$vertex_approximations[[comp_name]][, i])
        } else {
          coefficients_individual[[comp_name]] <- batch_inverse$vertex_approximations[[comp_name]]
        }
      }
      
      inverse_list[[sig]] <- list(
        vertex_approximations = coefficients_individual,
        reconstructed_signal = if (is.matrix(batch_inverse$reconstructed_signal)) {
          as.vector(batch_inverse$reconstructed_signal[, i])
        } else {
          batch_inverse$reconstructed_signal
        },
        reconstruction_error = if (is.vector(batch_inverse$reconstruction_error) && length(batch_inverse$reconstruction_error) > 1) {
          batch_inverse$reconstruction_error[i]
        } else {
          batch_inverse$reconstruction_error
        }
      )
    }
    
  } else {
    # Individual processing (original method)
    if (verbose && n_signals > 1) cat("Using individual processing...\n")
    
    forward_list <- list()
    inverse_list <- list()
    
    for (sig in signals) {
      if (verbose) cat("Processing signal:", sig, "\n")
      
      # Extract signal vector
      sig_vec <- as.numeric(data.in[[sig]])
      
      # Forward transform
      fwd <- sgwt_forward(sig_vec, eigenvectors, eigenvalues, params$scales, 
                         kernel_type = params$kernel_type)
      forward_list[[sig]] <- fwd
      
      # Inverse transform
      inv <- sgwt_inverse(fwd, eigenvectors, sig_vec)
      inverse_list[[sig]] <- inv
    }
  }
  
  # Update SGWT object
  SG$Forward <- forward_list
  SG$Inverse <- inverse_list
  
  if (verbose) cat("SGWT analysis completed.\n")
  
  return(SG)
}

#' Run SGCC weighted similarity analysis in Fourier domain
#'
#' @description Calculate energy-normalized weighted similarity between two signals
#' using Fourier domain coefficients directly (no vertex domain reconstruction).
#' Excludes DC component and uses energy-based weighting consistent with Parseval's theorem.
#'
#' @param signal1 Either a signal name (character) for SG object, or SGWT Forward result, or SGWT object
#' @param signal2 Either a signal name (character) for SG object, or SGWT Forward result, or SGWT object  
#' @param SG SGWT object (required if signal1/signal2 are signal names)
#' @param eps Small numeric for numerical stability (default: 1e-12)
#' @param validate Logical; if TRUE, check consistency (default: TRUE)
#' @param return_parts Logical; if TRUE, return detailed components (default: TRUE)
#' @param low_only Logical; if TRUE, compute only low-frequency similarity (default: FALSE)
#'
#' @return Similarity analysis results computed in Fourier domain
#' @export
#'
#' @examples
#' \donttest{
#' # Create example data and compute SGWT
#' data <- data.frame(x = runif(100), y = runif(100),
#'                   signal1 = rnorm(100), signal2 = rnorm(100))
#' SG <- initSGWT(data, signals = c("signal1", "signal2"))
#' SG <- runSpecGraph(SG, k = 15)
#' SG <- runSGWT(SG)
#' 
#' # Between two signals in same SGWT object
#' similarity <- runSGCC("signal1", "signal2", SG = SG)
#' print(similarity)
#' 
#' # Between two SGWT objects
#' data2 <- data.frame(x = runif(100), y = runif(100), signal = rnorm(100))
#' SG2 <- initSGWT(data2, signals = "signal")
#' SG2 <- runSpecGraph(SG2, k = 15)
#' SG2 <- runSGWT(SG2)
#' 
#' similarity2 <- runSGCC(SG, SG2)
#' print(similarity2)
#' }
runSGCC <- function(signal1, signal2, SG = NULL, eps = 1e-12, validate = TRUE, 
                    return_parts = TRUE, low_only = FALSE) {
  
  # Helper function to extract decomposition
  .get_decomp <- function(x, SG = NULL) {
    if (is.character(x) && length(x) == 1) {
      # x is a signal name
      if (is.null(SG) || is.null(SG$Forward)) {
        stop("SG object with Forward slot required when using signal names")
      }
      if (!x %in% names(SG$Forward)) {
        stop(paste("Signal", x, "not found in SGWT Forward results"))
      }
      return(SG$Forward[[x]])
    } else if (inherits(x, "SGWT")) {
      # x is an SGWT object - use first signal
      if (is.null(x$Forward) || length(x$Forward) == 0) {
        stop("SGWT object must have Forward results")
      }
      return(x$Forward[[1]])
    } else if (is.list(x) && !is.null(x$fourier_coefficients)) {
      # x is already a forward decomposition
      return(x)
    } else {
      stop("Invalid input type for signal")
    }
  }
  
  # Get decompositions
  A <- .get_decomp(signal1, SG)
  B <- .get_decomp(signal2, SG)
  
  # Extract filtered Fourier coefficients directly (no vertex domain reconstruction)
  fourier_a <- A$fourier_coefficients$filtered
  fourier_b <- B$fourier_coefficients$filtered
  
  if (is.null(fourier_a) || is.null(fourier_b)) {
    stop("Fourier coefficients not found in Forward results")
  }
  
  # Extract scaling (low-pass) Fourier coefficients, excluding DC component (first element)
  if (!"scaling" %in% names(fourier_a) || !"scaling" %in% names(fourier_b)) {
    stop("Scaling coefficients not found in filtered Fourier coefficients")
  }
  
  # Get scaling coefficients and exclude DC component (first element)
  f_low_a <- as.numeric(fourier_a$scaling)
  f_low_b <- as.numeric(fourier_b$scaling)
  
  # Exclude DC component (first coefficient, corresponding to λ = 0)
  if (length(f_low_a) > 1) f_low_a <- f_low_a[-1]
  if (length(f_low_b) > 1) f_low_b <- f_low_b[-1]
  
  # Handle non-finite values
  if (any(!is.finite(f_low_a))) {
    warning("Non-finite values found in scaling Fourier coefficients of signal1, replacing with 0")
    f_low_a[!is.finite(f_low_a)] <- 0
  }
  if (any(!is.finite(f_low_b))) {
    warning("Non-finite values found in scaling Fourier coefficients of signal2, replacing with 0")
    f_low_b[!is.finite(f_low_b)] <- 0
  }
  
  # Ensure both scaling coefficient vectors have the same length (for cross-object comparison)
  min_length <- min(length(f_low_a), length(f_low_b))
  if (length(f_low_a) != length(f_low_b)) {
    if (validate) {
      warning(paste("Scaling coefficients have different lengths:", length(f_low_a), "vs", length(f_low_b), 
                   ". Truncating to minimum length:", min_length))
    }
    f_low_a <- f_low_a[1:min_length]
    f_low_b <- f_low_b[1:min_length]
  }
  
  # Energies in Fourier domain (consistent with Parseval's theorem)
  E_low_a <- sum(f_low_a^2)
  E_low_b <- sum(f_low_b^2)
  
  # Low-frequency cosine similarity in Fourier domain
  c_low <- cosine_similarity(f_low_a, f_low_b, eps)
  
  # Short-circuit for low-only
  if (isTRUE(low_only)) {
    return(if (isTRUE(return_parts)) list(
      c_low = c_low, c_nonlow = NA_real_, w_low = 1.0, w_NL = 0.0, S = c_low,
      E_low_a = E_low_a, E_NL_a = NA_real_, E_low_b = E_low_b, E_NL_b = NA_real_,
      n = length(f_low_a), J = NA_integer_
    ) else c_low)
  }
  
  # Collect wavelet Fourier coefficients (non-low frequencies)
  wavelet_names_a <- names(fourier_a)[grep("^wavelet_scale_", names(fourier_a))]
  wavelet_names_b <- names(fourier_b)[grep("^wavelet_scale_", names(fourier_b))]
  
  if (length(wavelet_names_a) == 0 || length(wavelet_names_b) == 0) {
    stop("No wavelet Fourier coefficients found")
  }
  
  # Order wavelet coefficients by scale index
  ord_a <- order(as.integer(sub("^wavelet_scale_", "", wavelet_names_a)))
  ord_b <- order(as.integer(sub("^wavelet_scale_", "", wavelet_names_b)))
  wavelet_names_a <- wavelet_names_a[ord_a]
  wavelet_names_b <- wavelet_names_b[ord_b]
  
  # Handle different numbers of scales by using the minimum common scales
  min_scales <- min(length(wavelet_names_a), length(wavelet_names_b))
  if (length(wavelet_names_a) != length(wavelet_names_b)) {
    if (validate) {
      warning(paste("Different numbers of wavelet scales:", length(wavelet_names_a), "vs", length(wavelet_names_b), 
                   ". Using first", min_scales, "scales for comparison."))
    }
    wavelet_names_a <- wavelet_names_a[1:min_scales]
    wavelet_names_b <- wavelet_names_b[1:min_scales]
  }
  
  # Extract and flatten wavelet Fourier coefficients
  wavelet_coeffs_a <- lapply(wavelet_names_a, function(name) {
    as.numeric(fourier_a[[name]])
  })
  
  wavelet_coeffs_b <- lapply(wavelet_names_b, function(name) {
    as.numeric(fourier_b[[name]])
  })
  
  # Flatten all wavelet coefficients into single vectors
  f_wave_a <- unlist(wavelet_coeffs_a)
  f_wave_b <- unlist(wavelet_coeffs_b)
  
  # Handle non-finite values
  if (any(!is.finite(f_wave_a))) {
    warning("Non-finite values found in wavelet Fourier coefficients of signal1, replacing with 0")
    f_wave_a[!is.finite(f_wave_a)] <- 0
  }
  if (any(!is.finite(f_wave_b))) {
    warning("Non-finite values found in wavelet Fourier coefficients of signal2, replacing with 0")
    f_wave_b[!is.finite(f_wave_b)] <- 0
  }
  
  # Ensure both wavelet coefficient vectors have the same length (for cross-object comparison)
  min_wave_length <- min(length(f_wave_a), length(f_wave_b))
  if (length(f_wave_a) != length(f_wave_b)) {
    if (validate) {
      warning(paste("Wavelet coefficients have different lengths:", length(f_wave_a), "vs", length(f_wave_b), 
                   ". Truncating to minimum length:", min_wave_length))
    }
    f_wave_a <- f_wave_a[1:min_wave_length]
    f_wave_b <- f_wave_b[1:min_wave_length]
  }
  
  # Energies in Fourier domain for wavelet components
  E_NL_a <- sum(f_wave_a^2)
  E_NL_b <- sum(f_wave_b^2)
  
  # Number of scales and effective signal length (excluding DC)
  J <- length(wavelet_names_a)
  n <- length(f_low_a)  # Signal length after DC removal
  
  # Non-low cosine similarity in Fourier domain
  c_nonlow <- cosine_similarity(f_wave_a, f_wave_b, eps)
  
  # Energy-based macro weights (consistent with Parseval's theorem and Littlewood-Paley)
  w_low_a <- E_low_a / (E_low_a + E_NL_a + eps)
  w_low_b <- E_low_b / (E_low_b + E_NL_b + eps)
  w_low <- pmax(0, pmin(1, 0.5 * (w_low_a + w_low_b)))
  w_NL <- 1 - w_low
  S <- w_low * c_low + w_NL * c_nonlow
  
  # Validation
  if (isTRUE(validate)) {
    # Note: Length checks are now handled gracefully above with truncation
    # Only check for critical mismatches
    if (length(wavelet_names_a) != length(wavelet_names_b)) {
      stop("Number of wavelet scales must match")
    }
    
    # Check scales consistency if available
    if (!is.null(SG) && !is.null(SG$Parameters$scales)) {
      # Scales are consistent by design when using the same SGWT object
    }
    
    if (J < 1) {
      stop("Each SGWT decomposition must have at least 1 wavelet scale")
    }
  }
  
  # Return results
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


#' Print method for SGWT objects
#'
#' @param x SGWT object to print
#' @param ... Additional arguments passed to print methods
#' 
#' @return Invisibly returns the input SGWT object. Called for side effects (prints object summary to console).
#' @export
print.SGWT <- function(x, ...) {
  cat("SGWT Object\n")
  cat("===========\n")
  
  # Data information
  if (!is.null(x$Data)) {
    cat("Data:\n")
    cat("  Dimensions:", nrow(x$Data$data), "x", ncol(x$Data$data), "\n")
    cat("  Coordinates:", x$Data$x_col, ",", x$Data$y_col, "\n")
    cat("  Signals:", paste(x$Data$signals, collapse = ", "), "\n")
  }
  
  # Parameters
  if (!is.null(x$Parameters)) {
    cat("\nParameters:\n")
    cat("  k (neighbors):", x$Parameters$k, "\n")
    cat("  J (scales):", x$Parameters$J, "\n")
    cat("  Kernel type:", x$Parameters$kernel_type, "\n")
    cat("  Laplacian type:", x$Parameters$laplacian_type, "\n")
    if (!is.null(x$Parameters$scales)) {
      cat("  Scales:", paste(round(x$Parameters$scales, 4), collapse = ", "), "\n")
    }
  }
  
  # Status
  cat("\nStatus:\n")  
  cat("  Graph computed:", !is.null(x$Graph), "\n")
  cat("  Forward computed:", !is.null(x$Forward), "\n")
  cat("  Inverse computed:", !is.null(x$Inverse), "\n")
  
  if (!is.null(x$Inverse)) {
    cat("\nReconstruction Errors:\n")
    for (sig in names(x$Inverse)) {
      err <- x$Inverse[[sig]]$reconstruction_error
      if (!is.null(err)) {
        cat("  ", sig, ":", round(err, 6), "\n")
      }
    }
  }
  
  invisible(x)
}