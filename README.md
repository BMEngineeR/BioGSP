# BioGSP: Biological Graph Signal Processing for Spatial Data Analysis

## Overview

The BioGSP package provides a comprehensive implementation of Graph Signal Processing (GSP) methods including Spectral Graph Wavelet Transform (SGWT) for analyzing spatial patterns in biological data. This implementation is based on Hammond, Vandergheynst, and Gribonval (2011) "Wavelets on Graphs via Spectral Graph Theory" and provides tools for multi-scale analysis of spatial signals, including forward and inverse transforms, energy analysis, and visualization functions tailored for biological applications.

## Features

- **Multi-scale Analysis**: Decompose spatial signals into different frequency components
- **Graph Construction**: Automatic graph construction from spatial coordinates
- **Flexible Kernels**: Support for Mexican hat and Meyer wavelet kernels
- **Energy Analysis**: Analyze energy distribution across different scales
- **Visualization Tools**: Comprehensive plotting functions for results visualization
- **Reconstruction**: Perfect reconstruction capabilities with error analysis
- **SGCC Pipeline**: Spectral Graph Cross-Correlation for comparing graph signals
- **Minimal Dependencies**: Custom implementations reduce external package requirements

## SGCC Workflow Pipeline

The Spectral Graph Cross-Correlation (SGCC) pipeline enables comparison of two graph signals through a comprehensive multi-step process:

### Pipeline Overview

```mermaid
graph TD
    A[Input: Graph + Two Graph Signals] --> B[Cell Graph to Spot Graph Conversion]
    B --> C[Graph Fourier Transform - GFT]
    C --> D[Feature Extraction using Wavelet Kernel]
    D --> E{Branch Decision}
    E -->|k-bandlimited signals| F[Haar Mother Wavelet + Scaling Function]
    E -->|General signals| G[Other Kernels - Heat/Mexican Hat + Filter Bank]
    F --> H[Output: 1×N Vector]
    G --> I[Output: J+1×N Matrix]
    H --> J[Weighted Similarity Computation]
    I --> J
    J --> K[Energy Normalization]
    K --> L[SGCC Score]
```

### Detailed Workflow Steps

#### 1. **Input Graph and Two Graph Signals**
The pipeline starts with:
- A graph structure (e.g., cell graph, spot graph)
- Two associated graph signals to be compared
- Spatial coordinates defining the graph topology

#### 2. **Cell Graph to Spot Graph Conversion**
- Converts the graph representation to ensure both signals are in the same graph signal representation space
- Standardizes the graph structure for consistent analysis
- Maintains spatial relationships while enabling cross-modal comparison

#### 3. **Graph Fourier Transform (GFT)**
- **Input**: Node-domain features (signal values on graph nodes)
- **Output**: Spectral-domain features (frequency components of graph signals)
- Transforms spatial signals into the frequency domain using graph Laplacian eigendecomposition
- Enables frequency-based analysis of spatial patterns

#### 4. **Feature Extraction using Wavelet Kernel Function**
- Graph wavelet kernels are applied to extract multi-scale features from the spectral representation
- Captures both local and global signal characteristics
- Provides multi-resolution analysis of spatial patterns

#### 5. **Branch Decision**
The pipeline branches based on signal characteristics:

**Branch A: k-bandlimited Signals**
- **Condition**: If the graph signals are k-bandlimited
- **Method**: Haar mother wavelet with scaling function
- **Output**: 1×N vector (simple Haar case)
- **Use Case**: Signals with limited frequency content

**Branch B: General Signals**
- **Condition**: Non-bandlimited or complex signals
- **Method**: Other kernels (heat kernel or Mexican hat wavelet) with scaling + filter bank
- **Output**: (J+1)×N matrix (filter bank case)
- **Use Case**: Complex spatial patterns requiring multiple scales

#### 6. **Weighted Similarity Computation**
- Extracted multi-scale features are compared by computing similarity scores
- Similarity computed across graph nodes or spots
- Incorporates spatial weighting based on graph structure
- Accounts for both spectral and spatial relationships

#### 7. **Energy Normalization**
- Ensures similarity measures are normalized by the energy of the signals
- Provides fair comparison across different signal magnitudes
- Accounts for signal power variations
- Enables robust cross-correlation measurement

#### 8. **SGCC Score**
- **Final Output**: Spectral Graph Cross-Correlation score between the two signals
- Quantifies the similarity between graph signals in the spectral domain
- Ranges typically from -1 to 1 (depending on normalization)
- Higher scores indicate greater spectral similarity

## Installation

### From Source

```r
# Install required dependencies first
install.packages(c("Matrix", "igraph", "RANN", "RSpectra", "kneedle", 
                   "ggplot2", "patchwork", "ggpubr", "viridis"))

# Install the package from source
devtools::install_local("path/to/BioGSP")
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
- dplyr
- tidyr

## Quick Start

### Basic SGWT Analysis

```r
library(BioGSP)

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

### SGCC Pipeline Example

```r
library(BioGSP)

# Generate two related spatial signals for comparison
set.seed(123)
n_points <- 100
x_coords <- rep(1:10, each = 10) + rnorm(n_points, 0, 0.1)
y_coords <- rep(1:10, times = 10) + rnorm(n_points, 0, 0.1)

# Signal 1: Original pattern
signal1 <- sin(0.5 * x_coords) * cos(0.3 * y_coords) + rnorm(n_points, 0, 0.1)

# Signal 2: Related pattern with some variation
signal2 <- sin(0.4 * x_coords) * cos(0.35 * y_coords) + rnorm(n_points, 0, 0.15)

# Create data frame
comparison_data <- data.frame(
  x = x_coords, 
  y = y_coords, 
  signal1 = signal1,
  signal2 = signal2
)

# Step 1: Calculate eigendecomposition for the spatial data
eigen_result <- Cal_Eigen(comparison_data, k = 25, k_fold = 15)
knee_point <- eigen_result[[1]]
eigenvectors <- eigen_result[[2]]

# Step 2: Calculate SGCC score between the two signals
sgcc_score <- Cal_GCC(
  data.in = comparison_data,
  knee = knee_point,
  signal1 = "signal1",
  signal2 = "signal2", 
  eigenvector = eigenvectors
)

print(paste("SGCC Score:", round(sgcc_score, 4)))

# Step 3: Individual SGWT analysis for each signal
result1 <- SGWT(data.in = comparison_data, signal = "signal1", k = 25, J = 4)
result2 <- SGWT(data.in = comparison_data, signal = "signal2", k = 25, J = 4)

# Step 4: Energy analysis comparison
energy1 <- sgwt_energy_analysis(result1)
energy2 <- sgwt_energy_analysis(result2)

print("Energy Distribution Comparison:")
print(cbind(Signal1 = energy1$energy_ratio, Signal2 = energy2$energy_ratio))
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

The BioGSP package is particularly useful for:

### Biological Applications
- **Spatial Transcriptomics**: Comparing gene expression patterns across tissue regions
- **Single-cell Analysis**: Cross-correlation between cell populations and spatial locations
- **Tissue Architecture**: Analyzing spatial organization of cellular structures
- **Developmental Biology**: Tracking pattern formation and morphogenesis
- **Cancer Research**: Identifying spatial heterogeneity in tumor microenvironments

### General Applications
- **Neuroscience**: Brain connectivity and neural signal analysis
- **Environmental Science**: Spatial pattern analysis in ecological data
- **Image Processing**: Multi-scale image analysis on irregular domains
- **Social Networks**: Analyzing signals on social network graphs
- **Sensor Networks**: Spatial correlation analysis in distributed sensing

### SGCC-Specific Use Cases
- **Cross-modal Comparison**: Correlating different measurement modalities on the same spatial domain
- **Temporal Analysis**: Comparing spatial patterns across time points
- **Batch Effect Correction**: Identifying and correcting systematic differences between datasets
- **Quality Control**: Assessing spatial consistency in high-throughput experiments
- **Pattern Discovery**: Finding recurring spatial motifs across different samples

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

### SGCC Mathematical Framework

The Spectral Graph Cross-Correlation extends traditional cross-correlation to graph signals:

1. **Graph Fourier Transform**: For signals f₁, f₂ on graph G:
   - f̂₁(λᵢ) = ⟨f₁, χᵢ⟩ and f̂₂(λᵢ) = ⟨f₂, χᵢ⟩

2. **Spectral Cross-Correlation**: 
   - SGCC(f₁, f₂) = Σᵢ w(λᵢ) · f̂₁(λᵢ) · f̂₂*(λᵢ) / √(E₁ · E₂)
   
3. **Energy Normalization**:
   - E₁ = Σᵢ |f̂₁(λᵢ)|², E₂ = Σᵢ |f̂₂(λᵢ)|²
   
4. **Wavelet-based Features**:
   - Multi-scale decomposition: Wf(j,n) = Σᵢ gⱼ(λᵢ) f̂(λᵢ) χᵢ(n)

Where w(λᵢ) are frequency-dependent weights and gⱼ(λᵢ) are wavelet kernels at scale j.

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