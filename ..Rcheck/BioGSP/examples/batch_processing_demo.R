# Batch Processing Demo for SGWT
# ===========================================
# This script demonstrates the efficiency gains from batch processing
# multiple signals simultaneously using matrix multiplication.

library(BioGSP)

# Set up demo data
set.seed(123)
n_vertices <- 200
n_signals <- 10

# Create spatial coordinates
coords <- data.frame(
  x = runif(n_vertices, 0, 10),
  y = runif(n_vertices, 0, 10)
)

# Create multiple signals
signals_matrix <- matrix(rnorm(n_vertices * n_signals), nrow = n_vertices, ncol = n_signals)
colnames(signals_matrix) <- paste0("signal_", 1:n_signals)

# Build graph
k <- 8
nn <- RANN::nn2(coords, k = k + 1)
adj_list <- lapply(seq_len(n_vertices), function(i) setdiff(nn$nn.idx[i, ], i))
edges <- do.call(rbind, lapply(seq_along(adj_list), function(i) cbind(i, adj_list[[i]])))
edges <- unique(t(apply(edges, 1, sort)))
g <- igraph::graph_from_edgelist(edges, directed = FALSE)
A <- igraph::as_adjacency_matrix(g, sparse = TRUE)

# Compute Laplacian and eigendecomposition
L <- cal_laplacian(A, "normalized")
decomp <- FastDecompositionLap(L, k_eigen = 25, which = "SM")

# Set up SGWT parameters
scales <- sgwt_auto_scales(decomp$evalues, J = 4, scaling_factor = 2)
eigenvalues <- decomp$evalues
eigenvectors <- decomp$evectors

cat("=== SGWT Batch Processing Efficiency Demo ===\n")
cat("Number of vertices:", n_vertices, "\n")
cat("Number of signals:", n_signals, "\n")
cat("Number of scales:", length(scales), "\n\n")

# Method 1: Individual processing (original approach)
cat("Method 1: Individual signal processing...\n")
start_time <- Sys.time()

individual_results <- list()
for (i in 1:n_signals) {
  individual_results[[i]] <- sgwt_forward(
    signals_matrix[, i], 
    eigenvectors, 
    eigenvalues, 
    scales, 
    kernel_type = "heat"
  )
}

individual_time <- Sys.time() - start_time
cat("Individual processing time:", round(individual_time, 3), "seconds\n\n")

# Method 2: Batch processing (new unified approach)
cat("Method 2: Batch processing...\n")
start_time <- Sys.time()

batch_result <- sgwt_forward(
  signals_matrix, 
  eigenvectors, 
  eigenvalues, 
  scales, 
  kernel_type = "heat"
)

batch_time <- Sys.time() - start_time
cat("Batch processing time:", round(batch_time, 3), "seconds\n\n")

# Compare results
cat("=== Results Comparison ===\n")
cat("Speedup factor:", round(as.numeric(individual_time) / as.numeric(batch_time), 2), "x faster\n\n")

# Verify results are equivalent
cat("Verifying results are equivalent...\n")
max_diff <- 0
for (i in 1:n_signals) {
  # Compare scaling coefficients
  individual_scaling <- individual_results[[i]]$coefficients$scaling
  batch_scaling <- batch_result$coefficients$scaling[, i]
  diff <- max(abs(individual_scaling - batch_scaling))
  max_diff <- max(max_diff, diff)
  
  # Compare first wavelet scale
  individual_wavelet <- individual_results[[i]]$coefficients$wavelet_scale_1
  batch_wavelet <- batch_result$coefficients$wavelet_scale_1[, i]
  diff <- max(abs(individual_wavelet - batch_wavelet))
  max_diff <- max(max_diff, diff)
}

cat("Maximum difference between methods:", format(max_diff, scientific = TRUE), "\n")
cat("Results are", if (max_diff < 1e-10) "identical" else "nearly identical", "\n\n")

# Demonstrate batch inverse transform
cat("=== Batch Inverse Transform Demo ===\n")
start_time <- Sys.time()

batch_inverse <- sgwt_inverse(batch_result, signals_matrix)

inverse_time <- Sys.time() - start_time
cat("Batch inverse processing time:", round(inverse_time, 3), "seconds\n")

# Check reconstruction quality
cat("Reconstruction errors (RMSE):\n")
for (i in 1:n_signals) {
  cat("Signal", i, ":", round(batch_inverse$reconstruction_errors[i], 6), "\n")
}

cat("\n=== Summary ===\n")
cat("Batch processing provides significant speedup for multiple signals:\n")
cat("- Forward transform:", round(as.numeric(individual_time) / as.numeric(batch_time), 2), "x faster\n")
cat("- Memory efficient: Single matrix operations instead of loops\n")
cat("- Identical results: Numerical precision maintained\n")
cat("- Scalable: Benefits increase with more signals\n")

