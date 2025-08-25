## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load---------------------------------------------------------------------
library(BioGSP)

## ----synthetic_data-----------------------------------------------------------
# Generate synthetic spatial data
set.seed(123)
n_points <- 64
x_coords <- rep(1:8, each = 8) + rnorm(n_points, 0, 0.1)
y_coords <- rep(1:8, times = 8) + rnorm(n_points, 0, 0.1)

# Create a spatial pattern
signal_data <- sin(0.8 * x_coords) * cos(0.8 * y_coords) + 
               0.3 * sin(1.5 * x_coords * y_coords) + 
               rnorm(n_points, 0, 0.1)

demo_data <- data.frame(x = x_coords, y = y_coords, signal = signal_data)

# View the data structure
head(demo_data)

## ----sgwt_analysis------------------------------------------------------------
# Apply SGWT
result <- SGWT(data.in = demo_data, 
               signal = "signal",
               k = 4,               # number of nearest neighbors
               J = 3,               # number of scales
               scaling_factor = 2,
               k_fold = 5)          # k_fold * sqrt(64) = 5 * 8 = 40 < 64

# View reconstruction quality
cat("Reconstruction RMSE:", result$reconstruction_error, "\n")

## ----energy_analysis----------------------------------------------------------
# Analyze energy distribution
energy_df <- sgwt_energy_analysis(result)
print(energy_df)

# Plot energy distribution
library(ggplot2)
ggplot(energy_df, aes(x = scale, y = energy_ratio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Energy Distribution Across Scales",
       x = "Scale", y = "Energy Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----visualization, eval=FALSE------------------------------------------------
#  # Visualize SGWT decomposition (requires ggpubr)
#  if (requireNamespace("ggpubr", quietly = TRUE)) {
#    plots <- plot_sgwt_decomposition(result, demo_data, plot_scales = 1:3)
#    print(plots)
#  }

## ----custom_scales------------------------------------------------------------
# Define custom scales
custom_scales <- c(1.5, 0.8, 0.4, 0.2)

result_custom <- SGWT(data.in = demo_data,
                      signal = "signal",
                      scales = custom_scales,
                      k = 4,
                      k_fold = 5)

cat("Custom scales reconstruction RMSE:", result_custom$reconstruction_error, "\n")

## ----laplacian_types----------------------------------------------------------
# Normalized Laplacian (default)
result_norm <- SGWT(demo_data, "signal", k = 4, J = 3, laplacian_type = "normalized", k_fold = 5)

# Unnormalized Laplacian
result_unnorm <- SGWT(demo_data, "signal", k = 4, J = 3, laplacian_type = "unnormalized", k_fold = 5)

cat("Normalized Laplacian RMSE:", result_norm$reconstruction_error, "\n")
cat("Unnormalized Laplacian RMSE:", result_unnorm$reconstruction_error, "\n")

## ----gcc_example--------------------------------------------------------------
# Create a second signal
demo_data$signal2 <- cos(0.6 * x_coords) * sin(0.6 * y_coords) + rnorm(n_points, 0, 0.1)

# Calculate eigendecomposition
eigen_result <- Cal_Eigen(demo_data, k = 4, k_fold = 5, sensitivity = 2)

# Calculate cross-correlation
gcc_value <- Cal_GCC(data.in = demo_data,
                     knee = eigen_result[[1]],
                     signal1 = "signal",
                     signal2 = "signal2",
                     eigenvector = eigen_result[[2]])

cat("Graph Cross-Correlation:", gcc_value, "\n")

## ----demo---------------------------------------------------------------------
# Run the built-in demo
demo_result <- demo_sgwt()

# View demo results
print(demo_result$energy)

