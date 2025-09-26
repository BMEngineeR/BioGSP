#!/usr/bin/env Rscript

# BioGSP Simulation Pattern Demo
# This script demonstrates the new simulation functions and flexible column naming

# Load required libraries
library(BioGSP)
library(ggplot2)
library(patchwork)

# Set random seed for reproducibility
set.seed(123)

cat("=== BioGSP Simulation Pattern Demo ===\n\n")

# ============================================================================
# Part 1: Multiple Center Simulation
# ============================================================================

cat("1. Generating Multiple Center Patterns...\n")

# Define parameter sequences
Ra_seq <- c(8, 5, 3)   # Inner circle radii
Rb_seq <- c(12, 8, 4)  # Outer ring radii

# Generate simulation data
multicenter_data <- simulate_multiscale(
  grid_size = 40,
  n_centers = 2,
  Ra_seq = Ra_seq,
  Rb_seq = Rb_seq,
  seed = 123
)

cat("Generated", length(multicenter_data), "pattern combinations\n")

# Select one pattern for SGWT analysis
pattern_data <- multicenter_data[["simulated_Ra_8_Rb_12"]]
cat("Selected pattern dimensions:", nrow(pattern_data), "points\n")
cat("Column names:", paste(colnames(pattern_data), collapse = ", "), "\n\n")

# ============================================================================
# Part 2: SGWT Analysis with Flexible Column Naming
# ============================================================================

cat("2. Performing SGWT Analysis on signal_1 (Inner Circles)...\n")

# Perform SGWT analysis using custom column names (X, Y, signal_1)
sgwt_result <- SGWT(
  data.in = pattern_data,
  x_col = "X",           # Custom X coordinate column name
  y_col = "Y",           # Custom Y coordinate column name  
  signal = "signal_1",   # Analyze the inner circles
  k = 8,
  J = 3,
  scaling_factor = 2,
  kernel_type = "mexican_hat",
  k_fold = 6
)

cat("SGWT Analysis completed!\n")
cat("Reconstruction RMSE:", round(sgwt_result$reconstruction_error, 6), "\n")
cat("Number of scales:", length(sgwt_result$decomposition$scales), "\n")
cat("Scales:", paste(round(sgwt_result$decomposition$scales, 4), collapse = ", "), "\n\n")

# ============================================================================
# Part 3: Concentric Ring Simulation
# ============================================================================

cat("3. Generating Concentric Ring Patterns...\n")

# Generate ring patterns
ring_data <- simulate_ringpattern(
  grid_size = 40,
  radius_seq = c(5, 10),
  n_movements = 4
)

cat("Generated", length(ring_data), "ring patterns\n")

# Select ring data for analysis
selected_ring <- ring_data[["radius_10"]]
# Use movement step 2
movement_data <- selected_ring[selected_ring$movement_step == 2, ]

cat("Selected ring pattern dimensions:", nrow(movement_data), "points\n")
cat("Column names:", paste(colnames(movement_data), collapse = ", "), "\n\n")

# ============================================================================
# Part 4: SGWT Analysis on Ring Pattern
# ============================================================================

cat("4. Performing SGWT Analysis on signal_2 (Dynamic Ring)...\n")

# Analyze the dynamic ring using Meyer kernel
sgwt_ring_result <- SGWT(
  data.in = movement_data,
  x_col = "X",
  y_col = "Y",
  signal = "signal_2",   # Analyze the dynamic ring
  k = 8,
  J = 3,
  scaling_factor = 2,
  kernel_type = "meyer",
  k_fold = 6
)

cat("Ring SGWT Analysis completed!\n")
cat("Reconstruction RMSE:", round(sgwt_ring_result$reconstruction_error, 6), "\n\n")

# ============================================================================
# Part 5: Energy Analysis Comparison
# ============================================================================

cat("5. Comparing Energy Distributions...\n")

# Energy analysis for both patterns
energy_circles <- sgwt_energy_analysis(sgwt_result)
energy_rings <- sgwt_energy_analysis(sgwt_ring_result)

cat("Energy Distribution - Inner Circles:\n")
print(energy_circles)

cat("\nEnergy Distribution - Dynamic Rings:\n")
print(energy_rings)

# ============================================================================
# Part 6: Cross-Signal Analysis
# ============================================================================

cat("\n6. Performing Graph Cross-Correlation Analysis...\n")

# Calculate eigendecomposition for GCC
eigen_result <- Cal_Eigen(
  data.in = pattern_data,
  x_col = "X",
  y_col = "Y",
  k = 10,
  k_fold = 8,
  sensitivity = 2
)

knee_point <- eigen_result[[1]]
eigenvectors <- eigen_result[[2]]

# Calculate GCC between signal_1 and signal_2
gcc_value <- Cal_GCC(
  data.in = pattern_data,
  knee = knee_point,
  signal1 = "signal_1",  # Inner circles
  signal2 = "signal_2",  # Outer rings
  eigenvector = eigenvectors
)

cat("Graph Cross-Correlation between signals:", round(gcc_value, 4), "\n")
cat("Knee point (frequency cutoff):", knee_point, "\n\n")

# ============================================================================
# Part 7: Kernel Comparison
# ============================================================================

cat("7. Comparing Different Kernel Types...\n")

kernel_types <- c("mexican_hat", "meyer", "heat")
reconstruction_errors <- numeric(length(kernel_types))

for (i in seq_along(kernel_types)) {
  result <- SGWT(
    data.in = pattern_data,
    x_col = "X",
    y_col = "Y",
    signal = "signal_1",
    k = 8,
    J = 3,
    kernel_type = kernel_types[i],
    k_fold = 6
  )
  reconstruction_errors[i] <- result$reconstruction_error
}

comparison_df <- data.frame(
  Kernel = kernel_types,
  RMSE = reconstruction_errors
)

cat("Kernel Performance Comparison:\n")
print(comparison_df)

cat("\n=== Demo Complete ===\n")
cat("All functions successfully demonstrated flexible column naming support!\n")
cat("Key features tested:\n")
cat("- Multiple center pattern simulation\n")
cat("- Concentric ring pattern simulation\n")
cat("- SGWT analysis with custom column names (X, Y, signal_1, signal_2)\n")
cat("- Energy distribution analysis\n")
cat("- Graph Cross-Correlation analysis\n")
cat("- Multiple kernel type comparison\n")

