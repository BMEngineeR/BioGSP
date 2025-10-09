# Test and demonstrate different SGWT scaling functions
# This script shows how to use the various scaling function options

# Load necessary libraries
source("../R/sgwt_core.R")
source("../R/utilities.R")

# Create sample eigenvalues for testing
eigenvals <- c(0, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
scales <- c(2, 1, 0.5)

# Test different scaling functions
cat("Testing different scaling function types:\n")

scaling_types <- c("gaussian", "exponential", "polynomial", "meyer", "spline", "hann", "tight_frame")

# Test each scaling function individually
for (scaling_type in scaling_types) {
  cat("\n--- Testing", scaling_type, "scaling function ---\n")
  
  # Test the scaling function directly
  test_values <- sapply(c(0, 0.5, 1.0, 1.5, 2.0), function(x) {
    sgwt_scaling_kernel(x, scale_param = 1, scaling_type = scaling_type)
  })
  
  cat("Values at x = [0, 0.5, 1.0, 1.5, 2.0]:", test_values, "\n")
  
  # Test filter computation
  tryCatch({
    filters <- compute_sgwt_filters(eigenvals, scales, scaling_type = scaling_type)
    cat("Successfully computed", length(filters), "filters\n")
    cat("Filter 1 (scaling) range: [", min(filters[[1]]), ",", max(filters[[1]]), "]\n")
  }, error = function(e) {
    cat("Error computing filters:", e$message, "\n")
  })
}

# Compare all scaling functions visually (if possible)
cat("\n--- Comparing all scaling functions ---\n")
tryCatch({
  comparison <- compare_scaling_functions(x_range = c(0, 3), scale_param = 1, plot_results = FALSE)
  cat("Comparison data generated with", nrow(comparison), "points\n")
  cat("Available scaling types:", names(comparison)[-1], "\n")
  
  # Print some sample values
  sample_indices <- c(1, 50, 100, 150, 200)
  cat("\nSample values at different x points:\n")
  print(comparison[sample_indices, ])
  
}, error = function(e) {
  cat("Error in comparison:", e$message, "\n")
})

# Test with different scale parameters
cat("\n--- Testing with different scale parameters ---\n")
scale_params <- c(0.5, 1.0, 2.0)
x_test <- 1.0

for (scale_param in scale_params) {
  cat("Scale parameter:", scale_param, "\n")
  values <- sapply(scaling_types, function(type) {
    sgwt_scaling_kernel(x_test, scale_param = scale_param, scaling_type = type)
  })
  names(values) <- scaling_types
  print(values)
  cat("\n")
}

cat("Testing completed successfully!\n")