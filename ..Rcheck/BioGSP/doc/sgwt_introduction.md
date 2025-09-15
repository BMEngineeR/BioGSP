---
title: "Introduction to SGWT: Spectral Graph Wavelet Transform"
author: "SGWT Development Team"
date: "2025-08-19"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to SGWT}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Introduction

The SGWT package provides tools for analyzing spatial patterns using Spectral Graph Wavelet Transform. This vignette introduces the main concepts and demonstrates basic usage.

## Loading the Package


```r
library(SGWT)
#> Registered S3 method overwritten by 'quantmod':
#>   method            from
#>   as.zoo.data.frame zoo
```

## Basic Concepts

The Spectral Graph Wavelet Transform extends classical wavelet analysis to irregular domains by:

1. **Graph Construction**: Building a graph from spatial coordinates
2. **Spectral Analysis**: Using graph Laplacian eigendecomposition
3. **Wavelet Design**: Creating wavelets in the spectral domain
4. **Multi-scale Decomposition**: Analyzing signals at multiple scales

## Simple Example

Let's start with a simple synthetic dataset:


```r
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
#>           x         y      signal
#> 1 0.9439524 0.8928209  0.70770135
#> 2 0.9769823 2.0303529  0.00483162
#> 3 1.1558708 3.0448210 -0.71636034
#> 4 1.0070508 4.0053004 -0.74394155
#> 5 1.0129288 5.0922267 -0.12979638
#> 6 1.1715065 6.2050085 -0.14034003
```

## SGWT Analysis

Now let's apply SGWT to analyze this spatial signal:


```r
# Apply SGWT
result <- SGWT(data.in = demo_data, 
               signal = "signal",
               k = 4,               # number of nearest neighbors
               J = 3,               # number of scales
               scaling_factor = 2,
               k_fold = 5)          # k_fold * sqrt(64) = 5 * 8 = 40 < 64
#> Building graph from spatial coordinates...
#> Computing Laplacian and eigendecomposition...
#> Auto-generated scales: 1.1682, 0.5841, 0.292 
#> Performing SGWT decomposition...
#> Reconstruction RMSE: 0.159434

# View reconstruction quality
cat("Reconstruction RMSE:", result$reconstruction_error, "\n")
#> Reconstruction RMSE: 0.1594344
```

## Energy Analysis

We can analyze how energy is distributed across different scales:


```r
# Analyze energy distribution
energy_df <- sgwt_energy_analysis(result)
print(energy_df)
#>                   scale    energy energy_ratio scale_value
#> scaling         scaling 12.470257  0.878613364   1.1681514
#> wavelet_scale_1 scale_1  1.146569  0.080783486   1.1681514
#> wavelet_scale_2 scale_2  0.452992  0.031916332   0.5840757
#> wavelet_scale_3 scale_3  0.123293  0.008686819   0.2920379

# Plot energy distribution
library(ggplot2)
ggplot(energy_df, aes(x = scale, y = energy_ratio)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Energy Distribution Across Scales",
       x = "Scale", y = "Energy Ratio") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

![plot of chunk energy_analysis](figure/energy_analysis-1.png)

## Visualization

The package provides visualization functions to explore the decomposition:


```r
# Visualize SGWT decomposition (requires ggpubr)
if (requireNamespace("ggpubr", quietly = TRUE)) {
  plots <- plot_sgwt_decomposition(result, demo_data, plot_scales = 1:3)
  print(plots)
}
```

## Advanced Usage

### Custom Scales

You can specify custom scales for analysis:


```r
# Define custom scales
custom_scales <- c(1.5, 0.8, 0.4, 0.2)

result_custom <- SGWT(data.in = demo_data,
                      signal = "signal",
                      scales = custom_scales,
                      k = 4,
                      k_fold = 5)
#> Building graph from spatial coordinates...
#> Computing Laplacian and eigendecomposition...
#> Performing SGWT decomposition...
#> Reconstruction RMSE: 0.151932

cat("Custom scales reconstruction RMSE:", result_custom$reconstruction_error, "\n")
#> Custom scales reconstruction RMSE: 0.151932
```

### Different Laplacian Types

You can use different types of graph Laplacians:


```r
# Normalized Laplacian (default)
result_norm <- SGWT(demo_data, "signal", k = 4, J = 3, laplacian_type = "normalized", k_fold = 5)
#> Building graph from spatial coordinates...
#> Computing Laplacian and eigendecomposition...
#> Auto-generated scales: 1.1682, 0.5841, 0.292 
#> Performing SGWT decomposition...
#> Reconstruction RMSE: 0.159434

# Unnormalized Laplacian
result_unnorm <- SGWT(demo_data, "signal", k = 4, J = 3, laplacian_type = "unnormalized", k_fold = 5)
#> Building graph from spatial coordinates...
#> Computing Laplacian and eigendecomposition...
#> Auto-generated scales: 5.1519, 2.5759, 1.288 
#> Performing SGWT decomposition...
#> Reconstruction RMSE: 1.12435

cat("Normalized Laplacian RMSE:", result_norm$reconstruction_error, "\n")
#> Normalized Laplacian RMSE: 0.1594344
cat("Unnormalized Laplacian RMSE:", result_unnorm$reconstruction_error, "\n")
#> Unnormalized Laplacian RMSE: 1.12435
```

## Graph Cross-Correlation

The package also provides tools for analyzing relationships between spatial signals:


```r
# Create a second signal
demo_data$signal2 <- cos(0.6 * x_coords) * sin(0.6 * y_coords) + rnorm(n_points, 0, 0.1)

# Calculate eigendecomposition
eigen_result <- Cal_Eigen(demo_data, k = 4, k_fold = 5, sensitivity = 2)
#> [1] "cutoff of low-frequency FM: 24"
```

![plot of chunk gcc_example](figure/gcc_example-1.png)

```r

# Calculate cross-correlation
gcc_value <- Cal_GCC(data.in = demo_data,
                     knee = eigen_result[[1]],
                     signal1 = "signal",
                     signal2 = "signal2",
                     eigenvector = eigen_result[[2]])

cat("Graph Cross-Correlation:", gcc_value, "\n")
#> Graph Cross-Correlation: -0.4267128
```

## Demo Function

The package includes a demo function for quick testing:


```r
# Run the built-in demo
demo_result <- demo_sgwt()
#> === SGWT Demo ===
#> Generated synthetic data with 100 points
#> Building graph from spatial coordinates...
#> Computing Laplacian and eigendecomposition...
#> Error in eigs_real_sym(A, nrow(A), k, which, sigma, opts, mattype = "sym_dgCMatrix", : 'k' must satisfy 0 < k < nrow(A)

# View demo results
print(demo_result$energy)
#> Error in eval(expr, envir, enclos): object 'demo_result' not found
```

## Conclusion

The SGWT package provides a comprehensive toolkit for multi-scale analysis of spatial data. Key features include:

- Flexible graph construction from spatial coordinates
- Multiple wavelet kernel options
- Energy analysis across scales
- Visualization tools
- Cross-correlation analysis

This makes it particularly useful for applications in spatial biology, image processing, and network analysis.

## References

Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011). Wavelets on graphs via spectral graph theory. Applied and Computational Harmonic Analysis, 30(2), 129-150. 
