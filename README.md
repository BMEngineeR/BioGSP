# SGWT: Spectral Graph Wavelet Transform for Spatial Data Analysis

## Overview

The SGWT package provides a comprehensive implementation of Spectral Graph Wavelet Transform (SGWT) for analyzing spatial patterns in biological and other spatial data. This implementation is based on the seminal work by Hammond, Vandergheynst, and Gribonval (2011) "Wavelets on Graphs via Spectral Graph Theory".

## Features

- **Multi-scale Analysis**: Decompose spatial signals into different frequency components
- **Graph Construction**: Automatic graph construction from spatial coordinates
- **Flexible Kernels**: Support for Mexican hat and Meyer wavelet kernels
- **Energy Analysis**: Analyze energy distribution across different scales
- **Visualization Tools**: Comprehensive plotting functions for results visualization
- **Reconstruction**: Perfect reconstruction capabilities with error analysis
- **Minimal Dependencies**: Custom implementations reduce external package requirements

## Installation

### From Source

```r
# Install required dependencies first
install.packages(c("Matrix", "igraph", "RANN", "RSpectra", "kneedle", 
                   "ggplot2", "patchwork", "ggpubr", "viridis"))

# Install the package from source
devtools::install_local("path/to/SGWT")
```

### Dependencies

The package requires the following R packages:

**Required:**
- Matrix
- igraph  
- RANN
- RSpectra
- kneedle
- ggplot2
- patchwork
- ggpubr
- viridis
- methods

**Suggested:**
- spdep
- spatialEco
- sp
- gasper
- egg
- knitr
- rmarkdown
- testthat

## Quick Start

```r
library(SGWT)

# Generate synthetic spatial data
set.seed(123)
n_points <- 100
x_coords <- rep(1:10, each = 10) + rnorm(n_points, 0, 0.1)
y_coords <- rep(1:10, times = 10) + rnorm(n_points, 0, 0.1)
signal_data <- sin(0.5 * x_coords) * cos(0.3 * y_coords) + rnorm(n_points, 0, 0.1)

demo_data <- data.frame(x = x_coords, y = y_coords, signal = signal_data)

# Apply SGWT
result <- SGWT(data.in = demo_data, signal = "signal", k = 8, J = 4)

# View reconstruction error
print(result$reconstruction_error)

# Analyze energy distribution
energy_analysis <- sgwt_energy_analysis(result)
print(energy_analysis)

# Visualize results
plots <- plot_sgwt_decomposition(result, demo_data)
print(plots)
```

## Main Functions

### Core Functions

- `SGWT()`: Main function for SGWT analysis
- `sgwt_forward()`: Forward SGWT transform
- `sgwt_inverse()`: Inverse SGWT transform
- `sgwt_energy_analysis()`: Energy distribution analysis

### Utility Functions

- `cal_laplacian()`: Calculate graph Laplacian matrices
- `FastDecompositionLap()`: Fast eigendecomposition
- `gft()`: Graph Fourier Transform
- `Cal_Eigen()`: Eigenvalue analysis with knee detection
- `Cal_GCC()`: Graph Cross-Correlation
- `cosine_similarity()`: Custom cosine similarity implementation

### Visualization Functions

- `plot_sgwt_decomposition()`: Visualize SGWT components
- `demo_sgwt()`: Demonstration with synthetic data

### Kernel Functions

- `sgwt_kernel()`: Wavelet kernel functions
- `sgwt_scaling_kernel()`: Scaling function kernel
- `compute_sgwt_filters()`: Compute filter bank
- `sgwt_auto_scales()`: Automatic scale generation

## Usage Examples

### Basic SGWT Analysis

```r
# Load your spatial data with x, y coordinates and a signal
# data should have columns: x, y, signal

result <- SGWT(data.in = your_data, 
               signal = "signal",
               k = 25,              # number of nearest neighbors
               J = 5,               # number of scales
               scaling_factor = 2,   # scale separation factor
               kernel_type = "mexican_hat",
               laplacian_type = "normalized")
```

### Custom Scale Analysis

```r
# Define custom scales
custom_scales <- c(2.0, 1.0, 0.5, 0.25, 0.125)

result <- SGWT(data.in = your_data,
               signal = "signal",
               scales = custom_scales,
               k = 20)
```

### Energy Analysis

```r
# Analyze energy distribution across scales
energy_df <- sgwt_energy_analysis(result)

# View energy ratios
print(energy_df)

# Plot energy distribution
library(ggplot2)
ggplot(energy_df, aes(x = scale, y = energy_ratio)) +
  geom_bar(stat = "identity") +
  theme_minimal()
```

### Cross-Correlation Analysis

```r
# Calculate eigendecomposition for your spatial data
eigen_result <- Cal_Eigen(your_data, k = 25, k_fold = 15)

# Calculate cross-correlation between two signals
gcc_value <- Cal_GCC(data.in = your_data,
                     knee = eigen_result[[1]],
                     signal1 = "signal1",
                     signal2 = "signal2", 
                     eigenvector = eigen_result[[2]])
```

### Custom Cosine Similarity

```r
# Use the built-in cosine similarity function
x <- c(1, 2, 3, 4)
y <- c(2, 3, 4, 5)
similarity <- cosine_similarity(x, y)
print(similarity)
```

## Applications

The SGWT package is particularly useful for:

- **Spatial Biology**: Analyzing cell distribution patterns in tissue imaging
- **Neuroscience**: Brain connectivity and signal analysis
- **Environmental Science**: Spatial pattern analysis in ecological data
- **Image Processing**: Multi-scale image analysis on irregular domains
- **Social Networks**: Analyzing signals on social network graphs

## Theory

The Spectral Graph Wavelet Transform extends classical wavelet analysis to graphs and irregular domains. Key concepts:

1. **Graph Construction**: Build a graph from spatial coordinates using k-nearest neighbors
2. **Spectral Domain**: Use graph Laplacian eigendecomposition to define frequency
3. **Wavelet Design**: Create wavelets in the spectral domain using kernel functions
4. **Multi-scale Analysis**: Analyze signals at different scales simultaneously

### Mathematical Foundation

For a signal f on a graph G with Laplacian L and eigenvectors {χₗ}, the SGWT coefficients are:

- Scaling coefficients: Wf(λ) = ⟨f, hχₗ⟩
- Wavelet coefficients: Wf(j,λ) = ⟨f, gⱼχₗ⟩

Where h is the scaling function and gⱼ are wavelets at different scales.

## Implementation Notes

### Reduced Dependencies
This package implements its own cosine similarity function to reduce dependencies. The custom implementation provides the same functionality as external libraries while minimizing package requirements.

### Performance Considerations
- Uses sparse matrix operations where possible
- Implements fast eigendecomposition for large matrices
- Optimized graph construction using k-nearest neighbors

## References

1. Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011). Wavelets on graphs via spectral graph theory. Applied and Computational Harmonic Analysis, 30(2), 129-150.

2. Shuman, D. I., Narang, S. K., Frossard, P., Ortega, A., & Vandergheynst, P. (2013). The emerging field of signal processing on graphs: Extending high-dimensional data analysis to networks and other irregular domains. IEEE signal processing magazine, 30(3), 83-98.

## License

GPL-3

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## Contact

For questions and support, please open an issue on the GitHub repository. 