#!/usr/bin/env Rscript

# BioGSP Simulation Pattern Demo - Updated for New Workflow
# This script demonstrates the new SGWT workflow: initSGWT -> runSpecGraph -> runSGWT -> runSGCC

# Load required libraries
library(BioGSP)
library(ggplot2)
library(patchwork)

# Set random seed for reproducibility
set.seed(123)

cat("=== BioGSP Simulation Pattern Demo (New Workflow) ===\n\n")

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
# Part 2: NEW SGWT Workflow - Initialize and Build Graph
# ============================================================================

cat("2. NEW WORKFLOW - Step 1: Initialize SGWT Object...\n")

# Step 1: Initialize SGWT object with custom column names
SG <- initSGWT(
  data.in = pattern_data,
  x_col = "X",           # Custom X coordinate column name
  y_col = "Y",           # Custom Y coordinate column name  
  signals = c("signal_1", "signal_2"),  # Analyze both signals
  k = 8,
  J = 3,
  scaling_factor = 2,
  kernel_type = "heat"
)

cat("SGWT object initialized!\n")
print(SG)

cat("\nStep 2: Build Spectral Graph...\n")
# Step 2: Build spectral graph
SG <- runSpecGraph(SG, verbose = TRUE)

cat("\nStep 3: Run SGWT Analysis...\n")
# Step 3: Run SGWT forward and inverse transforms
SG <- runSGWT(SG, verbose = TRUE)

cat("SGWT Analysis completed!\n")
print(SG)

# ============================================================================
# Part 3: Concentric Ring Simulation
# ============================================================================

cat("\n3. Generating Concentric Ring Patterns...\n")

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
# Part 4: NEW WORKFLOW - Ring Pattern Analysis
# ============================================================================

cat("4. NEW WORKFLOW - Ring Pattern Analysis...\n")

# Initialize second SGWT object for ring pattern
SG_ring <- initSGWT(
  data.in = movement_data,
  x_col = "X",
  y_col = "Y",
  signals = c("signal_2"),   # Analyze the dynamic ring
  k = 8,
  J = 3,
  scaling_factor = 2,
  kernel_type = "meyer"
)

# Build graph and run analysis
SG_ring <- runSpecGraph(SG_ring, verbose = FALSE)
SG_ring <- runSGWT(SG_ring, verbose = FALSE)

cat("Ring SGWT Analysis completed!\n")

# ============================================================================
# Part 5: Energy Analysis Comparison
# ============================================================================

cat("\n5. Comparing Energy Distributions...\n")

# Energy analysis for both patterns
energy_circles <- sgwt_energy_analysis(SG, "signal_1")
energy_rings <- sgwt_energy_analysis(SG_ring, "signal_2")

cat("Energy Distribution - Inner Circles (signal_1):\n")
print(energy_circles)

cat("\nEnergy Distribution - Dynamic Rings (signal_2):\n")
print(energy_rings)

# ============================================================================
# Part 6: NEW WORKFLOW - Cross-Signal Similarity Analysis
# ============================================================================

cat("\n6. NEW WORKFLOW - Cross-Signal Similarity Analysis...\n")

# Calculate similarity between signal_1 and signal_2 in the same object
similarity_within <- runSGCC("signal_1", "signal_2", SG = SG, return_parts = TRUE)

cat("Similarity between signal_1 and signal_2 (within same graph):\n")
cat("Overall Similarity (S):", round(similarity_within$S, 4), "\n")
cat("Low-frequency similarity:", round(similarity_within$c_low, 4), "\n")
cat("Non-low-frequency similarity:", round(similarity_within$c_nonlow, 4), "\n")
cat("Energy weights - Low:", round(similarity_within$w_low, 3), ", Non-low:", round(similarity_within$w_NL, 3), "\n\n")

# Calculate similarity between different SGWT objects
similarity_cross <- runSGCC(SG, SG_ring, return_parts = FALSE)
cat("Cross-object similarity (circles vs rings):", round(similarity_cross, 4), "\n\n")

# ============================================================================
# Part 7: Kernel Comparison with New Workflow
# ============================================================================

cat("7. Comparing Different Kernel Types with New Workflow...\n")

kernel_types <- c("mexican_hat", "meyer", "heat")
reconstruction_errors <- numeric(length(kernel_types))

for (i in seq_along(kernel_types)) {
  # Create temporary SGWT object for each kernel
  SG_temp <- initSGWT(
    data.in = pattern_data,
    x_col = "X",
    y_col = "Y",
    signals = "signal_1",
    k = 8,
    J = 3,
    kernel_type = kernel_types[i]
  )
  
  SG_temp <- runSpecGraph(SG_temp, verbose = FALSE)
  SG_temp <- runSGWT(SG_temp, verbose = FALSE)
  
  reconstruction_errors[i] <- SG_temp$Inverse$signal_1$reconstruction_error
}

comparison_df <- data.frame(
  Kernel = kernel_types,
  RMSE = reconstruction_errors
)

cat("Kernel Performance Comparison:\n")
print(comparison_df)

# ============================================================================
# Part 8: Low-frequency Only Analysis
# ============================================================================

cat("\n8. Low-frequency Only Similarity Analysis...\n")

# Compare using only low-frequency components
similarity_low <- runSGCC("signal_1", "signal_2", SG = SG, low_only = TRUE, return_parts = TRUE)

cat("Low-frequency only similarity:\n")
cat("Similarity score:", round(similarity_low$S, 4), "\n")
cat("Note: c_nonlow is NA for low-only analysis\n\n")

# ============================================================================
# Part 9: Demo Visualization (if ggpubr is available)
# ============================================================================

cat("9. Demonstration Complete!\n")

# Try to create a simple visualization
if (requireNamespace("ggpubr", quietly = TRUE)) {
  cat("Creating visualization plots...\n")
  tryCatch({
    plots <- plot_sgwt_decomposition(SG, "signal_1")
    cat("Visualization plots created successfully!\n")
  }, error = function(e) {
    cat("Note: Visualization requires additional packages\n")
  })
} else {
  cat("Note: Install 'ggpubr' for visualization features\n")
}

cat("\n=== Demo Complete ===\n")
cat("NEW WORKFLOW successfully demonstrated!\n")
cat("Key features tested:\n")
cat("- initSGWT(): Initialize SGWT objects with flexible column naming\n")
cat("- runSpecGraph(): Build spectral graph structure\n")
cat("- runSGWT(): Perform forward and inverse SGWT transforms\n")
cat("- runSGCC(): Calculate energy-normalized weighted similarity\n")
cat("- sgwt_energy_analysis(): Analyze energy distribution across scales\n")
cat("- Multiple kernel type comparison\n")
cat("- Low-frequency only analysis\n")
cat("- Cross-object similarity comparison\n")
