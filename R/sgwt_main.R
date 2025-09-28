#' Initialize SGWT object
#'
#' @description Build an SGWT object with Data and Parameters slots, validate inputs.
#' 
#' @param data.in Data frame containing spatial coordinates and signal data
#' @param x_col Character string specifying the column name for X coordinates (default: "x")
#' @param y_col Character string specifying the column name for Y coordinates (default: "y") 
#' @param signals Character vector of signal column names to analyze. If NULL, all non-coordinate columns are used.
#' @param k Number of nearest neighbors for graph construction (default: 25)
#' @param scales Vector of scales for the wavelets. If NULL, scales are auto-generated.
#' @param J Number of scales to generate if scales is NULL (default: 5)
#' @param scaling_factor Scaling factor between consecutive scales (default: 2)
#' @param kernel_type Kernel family ("mexican_hat", "meyer", or "heat") (default: "heat")
#' @param laplacian_type Type of graph Laplacian ("unnormalized", "normalized", or "randomwalk") (default: "normalized")
#'
#' @return SGWT object with Data and Parameters slots initialized
#' @export
#'
#' @examples
#' \dontrun{
#' # Initialize SGWT object
#' data <- data.frame(x = runif(100), y = runif(100), 
#'                   signal1 = rnorm(100), signal2 = rnorm(100))
#' SG <- initSGWT(data, signals = c("signal1", "signal2"))
#' }
initSGWT <- function(data.in, x_col = "x", y_col = "y", signals = NULL,
                     k = 25, scales = NULL, J = 5, scaling_factor = 2,
                     kernel_type = "heat", laplacian_type = "normalized") {
  
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
      k = k,
      scales = scales,
      J = J,
      scaling_factor = scaling_factor,
      kernel_type = kernel_type,
      laplacian_type = laplacian_type
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
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return Updated SGWT object with Graph slot populated
#' @export
#'
#' @examples
#' \dontrun{
#' SG <- initSGWT(data)
#' SG <- runSpecGraph(SG)
#' }
runSpecGraph <- function(SG, verbose = TRUE) {
  
  # Validate input
  if (!inherits(SG, "SGWT")) {
    stop("Input must be an SGWT object from initSGWT()")
  }
  if (is.null(SG$Data)) {
    stop("SGWT object must have Data slot initialized")
  }
  
  # Extract data and parameters
  data.in <- SG$Data$data
  x_col <- SG$Data$x_col
  y_col <- SG$Data$y_col
  k <- SG$Parameters$k
  laplacian_type <- SG$Parameters$laplacian_type
  
  if (verbose) cat("Building graph from spatial coordinates...\\n")
  
  # Build k-nearest neighbor graph
  nn <- RANN::nn2(data.in[, c(x_col, y_col)], k = k + 1)
  adj_list <- lapply(seq_len(nrow(data.in)), function(i) setdiff(nn$nn.idx[i, ], i))
  edges <- do.call(rbind, lapply(seq_along(adj_list), function(i) cbind(i, adj_list[[i]])))
  edges <- unique(t(apply(edges, 1, sort)))
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  A <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  
  if (verbose) cat("Computing Laplacian and eigendecomposition...\\n")
  
  # Compute Laplacian matrix
  L <- cal_laplacian(A, laplacian_type)
  
  # Eigendecomposition
  decomp <- FastDecompositionLap(L, k_neighbor = 25, which = "SM")
  
  # Update SGWT object
  SG$Graph <- list(
    adjacency_matrix = A,
    laplacian_matrix = L,
    eigenvalues = decomp$evalues,
    eigenvectors = decomp$evectors
  )
  
  if (verbose) cat("Graph construction completed.\\n")
  
  return(SG)
}

#' Run SGWT forward and inverse transforms for all signals
#'
#' @description Perform SGWT analysis on all signals in the SGWT object.
#' Assumes Graph slot is populated by runSpecGraph().
#'
#' @param SG SGWT object with Graph slot populated
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return Updated SGWT object with Forward and Inverse slots populated
#' @export
#'
#' @examples
#' \dontrun{
#' SG <- initSGWT(data)
#' SG <- runSpecGraph(SG)
#' SG <- runSGWT(SG)
#' }
runSGWT <- function(SG, verbose = TRUE) {
  
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
  
  # Auto-generate scales if not provided
  if (is.null(params$scales)) {
    lmax <- max(eigenvalues) * 0.95
    params$scales <- sgwt_auto_scales(lmax, params$J, params$scaling_factor)
    if (verbose) {
      cat(paste("Auto-generated scales:", paste(round(params$scales, 4), collapse = ", "), "\\n"))
    }
  }
  
  if (verbose) cat("Performing SGWT analysis for", length(signals), "signals...\\n")
  
  # Process each signal
  forward_list <- list()
  inverse_list <- list()
  
  for (sig in signals) {
    if (verbose) cat("Processing signal:", sig, "\\n")
    
    # Extract signal vector
    sig_vec <- as.numeric(data.in[[sig]])
    
    # Forward transform
    fwd <- sgwt_forward(sig_vec, eigenvectors, eigenvalues, params$scales, 
                       kernel_type = params$kernel_type)
    forward_list[[sig]] <- fwd
    
    # Inverse transform
    inv <- sgwt_inverse(fwd, sig_vec)
    inverse_list[[sig]] <- inv
  }
  
  # Update SGWT object
  SG$Forward <- forward_list
  SG$Inverse <- inverse_list
  SG$Parameters <- params  # Update with auto-generated scales if applicable
  
  if (verbose) cat("SGWT analysis completed.\\n")
  
  return(SG)
}

#' Run SGCC weighted similarity analysis
#'
#' @description Calculate energy-normalized weighted similarity between two signals
#' from SGWT Forward results or between two SGWT objects.
#'
#' @param signal1 Either a signal name (character) for SG object, or SGWT Forward result, or SGWT object
#' @param signal2 Either a signal name (character) for SG object, or SGWT Forward result, or SGWT object  
#' @param SG SGWT object (required if signal1/signal2 are signal names)
#' @param eps Small numeric for numerical stability (default: 1e-12)
#' @param validate Logical; if TRUE, check consistency (default: TRUE)
#' @param return_parts Logical; if TRUE, return detailed components (default: TRUE)
#' @param low_only Logical; if TRUE, compute only low-frequency similarity (default: FALSE)
#'
#' @return Similarity analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # Between two signals in same SGWT object
#' similarity <- runSGCC("signal1", "signal2", SG = SG_object)
#' 
#' # Between two SGWT objects
#' similarity <- runSGCC(SG_object1, SG_object2)
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
    } else if (is.list(x) && !is.null(x$coefficients)) {
      # x is already a forward decomposition
      return(x)
    } else {
      stop("Invalid input type for signal")
    }
  }
  
  # Get decompositions
  A <- .get_decomp(signal1, SG)
  B <- .get_decomp(signal2, SG)
  
  # Extract scaling vectors
  y_low_a <- as.numeric(A$coefficients$scaling)
  y_low_b <- as.numeric(B$coefficients$scaling)
  
  # Handle non-finite values
  if (any(!is.finite(y_low_a))) {
    warning("Non-finite values found in scaling coefficients of signal1, replacing with 0")
    y_low_a[!is.finite(y_low_a)] <- 0
  }
  if (any(!is.finite(y_low_b))) {
    warning("Non-finite values found in scaling coefficients of signal2, replacing with 0")
    y_low_b[!is.finite(y_low_b)] <- 0
  }
  
  # Energies
  E_low_a <- sum(y_low_a^2)
  E_low_b <- sum(y_low_b^2)
  
  # Low-frequency cosine similarity
  c_low <- cosine_similarity(y_low_a, y_low_b, eps)
  
  # Short-circuit for low-only
  if (isTRUE(low_only)) {
    return(if (isTRUE(return_parts)) list(
      c_low = c_low, c_nonlow = NA_real_, w_low = 1.0, w_NL = 0.0, S = c_low,
      E_low_a = E_low_a, E_NL_a = NA_real_, E_low_b = E_low_b, E_NL_b = NA_real_,
      n = length(y_low_a), J = NA_integer_
    ) else c_low)
  }
  
  # Collect wavelet coefficients
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
  
  # Non-low cosine similarity
  c_nonlow <- cosine_similarity(v_a, v_b, eps)
  
  # Macro weights (energy normalization)
  w_low_a <- E_low_a / (E_low_a + E_NL_a + eps)
  w_low_b <- E_low_b / (E_low_b + E_NL_b + eps)
  w_low <- pmax(0, pmin(1, 0.5 * (w_low_a + w_low_b)))
  w_NL <- 1 - w_low
  S <- w_low * c_low + w_NL * c_nonlow
  
  # Validation
  if (isTRUE(validate)) {
    if (length(y_low_a) != length(y_low_b)) {
      stop("Scaling coefficients must have the same length")
    }
    if (nrow(W_a) != nrow(W_b)) {
      stop("Signal lengths must match")
    }
    if (ncol(W_a) != ncol(W_b)) {
      stop("Number of wavelet scales must match")
    }
    
    # Check scales consistency
    scales_a <- A$scales
    scales_b <- B$scales
    if (!is.null(scales_a) && !is.null(scales_b)) {
      if (length(scales_a) != length(scales_b) || !isTRUE(all.equal(scales_a, scales_b))) {
        warning("Scales differ between signals")
      }
    } else if (!is.null(scales_a) || !is.null(scales_b)) {
      warning("Scales present in only one signal")
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

#' Print method for SGWT objects
#'
#' @param x SGWT object to print
#' @param ... Additional arguments passed to print methods
#' @export
print.SGWT <- function(x, ...) {
  cat("SGWT Object\\n")
  cat("===========\\n")
  
  # Data information
  if (!is.null(x$Data)) {
    cat("Data:\\n")
    cat("  Dimensions:", nrow(x$Data$data), "x", ncol(x$Data$data), "\\n")
    cat("  Coordinates:", x$Data$x_col, ",", x$Data$y_col, "\\n")
    cat("  Signals:", paste(x$Data$signals, collapse = ", "), "\\n")
  }
  
  # Parameters
  if (!is.null(x$Parameters)) {
    cat("\\nParameters:\\n")
    cat("  k (neighbors):", x$Parameters$k, "\\n")
    cat("  J (scales):", x$Parameters$J, "\\n")
    cat("  Kernel type:", x$Parameters$kernel_type, "\\n")
    cat("  Laplacian type:", x$Parameters$laplacian_type, "\\n")
    if (!is.null(x$Parameters$scales)) {
      cat("  Scales:", paste(round(x$Parameters$scales, 4), collapse = ", "), "\\n")
    }
  }
  
  # Status
  cat("\\nStatus:\\n")
  cat("  Graph computed:", !is.null(x$Graph), "\\n")
  cat("  Forward computed:", !is.null(x$Forward), "\\n")
  cat("  Inverse computed:", !is.null(x$Inverse), "\\n")
  
  if (!is.null(x$Inverse)) {
    cat("\\nReconstruction Errors:\\n")
    for (sig in names(x$Inverse)) {
      err <- x$Inverse[[sig]]$reconstruction_error
      if (!is.null(err)) {
        cat("  ", sig, ":", round(err, 6), "\\n")
      }
    }
  }
}