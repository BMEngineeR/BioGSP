# BioGSP Developer Documentation

This document provides comprehensive technical documentation for the BioGSP package, including detailed function specifications, implementation details, and advanced usage examples.

## Table of Contents

- [Package Architecture](#package-architecture)
- [SGWT Object Structure](#sgwt-object-structure)  
- [Core Functions](#core-functions)
- [Utility Functions](#utility-functions)
- [Visualization Functions](#visualization-functions)
- [Kernel Functions](#kernel-functions)
- [Mathematical Framework](#mathematical-framework)
- [Implementation Details](#implementation-details)
- [Advanced Examples](#advanced-examples)
- [Performance Considerations](#performance-considerations)
- [Contributing](#contributing)

## Package Architecture

The BioGSP package is organized into several key components:

```
BioGSP/
├── R/
│   ├── sgwt_main.R          # Main SGWT functions
│   ├── sgwt_core.R          # Core transform functions  
│   ├── utilities.R          # Utility and similarity functions
│   ├── visualization.R      # Plotting functions
│   └── kernels.R           # Wavelet kernel implementations
├── data/                    # Example datasets
├── man/                     # Documentation files
├── vignettes/              # Package vignettes
└── inst/examples/          # Example scripts
```

## SGWT Object Structure

When you run `SGWT()`, it returns a comprehensive list containing all components of the spectral graph wavelet transform analysis.

### Main SGWT Result Object

```r
sgwt_result <- SGWT(data.in = your_data, signal = "signal", k = 25, J = 4)

# Structure (when return_all = TRUE, which is default):
sgwt_result$
├── decomposition        # SGWT decomposition results from sgwt_forward()
├── reconstructed_signal # Reconstructed signal from inverse transform
├── reconstruction_error # RMSE between original and reconstructed signal
├── original_signal     # Original input signal values
├── graph_info          # Graph construction information
├── data               # Original input data frame
└── parameters         # Analysis parameters used
```

### Detailed Component Descriptions

#### 1. **`decomposition`** (List)
- **Type**: List containing the core SGWT decomposition results
- **Content**: Output from `sgwt_forward()` function
- **Components**:
  ```r
  decomposition$
  ├── coefficients      # List of coefficient vectors for each scale
  │   ├── scaling       # Scaling function coefficients (low-pass/approximation)
  │   ├── wavelet_scale_1 # Wavelet coefficients at scale 1 (finest scale)
  │   ├── wavelet_scale_2 # Wavelet coefficients at scale 2
  │   └── ...           # Additional scales up to J
  ├── filters           # List of filter values in spectral domain
  ├── scales           # Vector of scale parameters used
  ├── eigenvalues      # Eigenvalues of the graph Laplacian
  └── eigenvectors     # Eigenvectors of the graph Laplacian
  ```
- **Meaning**: The actual multi-scale wavelet decomposition of your spatial signal
- **Interpretation**: 
  - `coefficients$scaling`: Smooth, low-frequency approximation of the signal
  - `coefficients$wavelet_scale_*`: Detail coefficients at different spatial scales
  - Higher scale numbers → coarser spatial details
  - Lower scale numbers → finer spatial details

#### 2. **`reconstructed_signal`** (Vector)
- **Type**: Numeric vector of length N (where N = number of spatial points)
- **Content**: Signal reconstructed from the SGWT decomposition using inverse transform
- **Meaning**: Result of `sgwt_inverse()` applied to the decomposition
- **Quality Check**: Should closely match `original_signal` (check `reconstruction_error`)
- **Interpretation**: Validates the perfect reconstruction property of the transform

#### 3. **`reconstruction_error`** (Scalar)
- **Type**: Single numeric value
- **Content**: Root Mean Square Error between original and reconstructed signal
- **Formula**: `sqrt(mean((original_signal - reconstructed_signal)^2))`
- **Interpretation**: 
  - Values close to 0 indicate perfect reconstruction
  - Higher values may indicate numerical issues or implementation problems
- **Usage**: Quality control metric for the transform

#### 4. **`original_signal`** (Vector)
- **Type**: Numeric vector of length N
- **Content**: The original signal values that were analyzed
- **Meaning**: Copy of the input signal for reference and validation
- **Interpretation**: Baseline for comparison with reconstructed signal

#### 5. **`graph_info`** (List)
- **Components**:
  ```r
  graph_info$
  ├── adjacency_matrix    # Sparse adjacency matrix of the spatial graph
  ├── laplacian_matrix    # Graph Laplacian matrix (normalized/unnormalized)
  ├── eigenvalues         # Eigenvalues of the Laplacian
  └── eigenvectors        # Eigenvectors of the Laplacian
  ```
- **Meaning**: Complete spectral information about the spatial graph structure
- **Interpretation**: Contains all graph-theoretic information needed for analysis and visualization

#### 6. **`data`** (Data Frame)
- **Type**: Data frame (copy of original input)
- **Content**: Original input data with spatial coordinates and signal
- **Meaning**: Preserves the original input for reference and further analysis
- **Usage**: Needed for visualization functions and further processing

#### 7. **`parameters`** (List)
- **Components**:
  ```r
  parameters$
  ├── k                   # Number of nearest neighbors used
  ├── scales             # Vector of scales used in the transform
  ├── J                  # Number of wavelet scales
  ├── kernel_type        # Wavelet kernel type used
  └── laplacian_type     # Type of Laplacian normalization
  ```
- **Meaning**: All parameters used in the analysis for reproducibility
- **Interpretation**: Enables exact reproduction of results and parameter tracking

### Usage Examples

```r
# Access the core decomposition
decomp <- sgwt_result$decomposition

# Access specific coefficient components
scaling_coeffs <- decomp$coefficients$scaling          # Low-frequency approximation
wavelet_1 <- decomp$coefficients$wavelet_scale_1       # Finest scale details
wavelet_2 <- decomp$coefficients$wavelet_scale_2       # Medium scale details
# ... additional scales as wavelet_scale_3, wavelet_scale_4, etc.

# Access filter information
filters <- decomp$filters                              # Spectral domain filters
scales_used <- decomp$scales                           # Scale parameters
eigenvals <- decomp$eigenvalues                        # Graph eigenvalues
eigenvecs <- decomp$eigenvectors                       # Graph eigenvectors

# Check reconstruction quality
error <- sgwt_result$reconstruction_error
print(paste("Reconstruction RMSE:", round(error, 10)))

# Get graph properties
eigenvalues <- sgwt_result$graph_info$eigenvalues
eigenvectors <- sgwt_result$graph_info$eigenvectors
adjacency <- sgwt_result$graph_info$adjacency_matrix

# Access parameters
scales_used <- sgwt_result$parameters$scales
k_neighbors <- sgwt_result$parameters$k
kernel_type <- sgwt_result$parameters$kernel_type

# Compare original vs reconstructed
original <- sgwt_result$original_signal
reconstructed <- sgwt_result$reconstructed_signal
plot(original, reconstructed)
abline(0, 1, col = "red")  # Perfect reconstruction line

# Analyze energy at each scale
scaling_energy <- sum(scaling_coeffs^2)
wavelet_energies <- sapply(seq_along(scales_used), function(i) {
  coeff_name <- paste0("wavelet_scale_", i)
  sum(decomp$coefficients[[coeff_name]]^2)
})
print("Energy distribution:")
print(c(scaling = scaling_energy, wavelet_energies))
```

## Core Functions

### `SGWT()`

Main function for Spectral Graph Wavelet Transform analysis.

```r
SGWT(data.in, signal, x_col = "x", y_col = "y", k = 25, 
     scales = NULL, J = 5, scaling_factor = 2, 
     kernel_type = "mexican_hat", laplacian_type = "normalized",
     return_all = TRUE, ...)
```

**Parameters:**
- `data.in`: Data frame with spatial coordinates and signal
- `signal`: Column name containing the signal to analyze
- `x_col`, `y_col`: Column names for spatial coordinates
- `k`: Number of nearest neighbors for graph construction
- `scales`: Custom scale vector (optional)
- `J`: Number of wavelet scales (if scales not provided)
- `scaling_factor`: Scale separation factor
- `kernel_type`: Wavelet kernel ("mexican_hat" or "meyer")
- `laplacian_type`: Laplacian normalization ("normalized" or "unnormalized")
- `return_all`: If TRUE, return complete analysis; if FALSE, return only decomposition

**Returns:**
- Complete SGWT analysis object (see [SGWT Object Structure](#sgwt-object-structure))

### `sgwt_forward()`

Forward SGWT transform - decomposes a graph signal into wavelet coefficients.

```r
sgwt_forward(signal, eigenvalues, eigenvectors, scales, 
             kernel_type = "mexican_hat")
```

**Parameters:**
- `signal`: Graph signal to decompose
- `eigenvalues`: Graph Laplacian eigenvalues
- `eigenvectors`: Graph Laplacian eigenvectors
- `scales`: Scale parameters for wavelet functions
- `kernel_type`: Type of wavelet kernel

**Returns:**
- List with coefficients, filters, scales, eigenvalues, and eigenvectors

### `sgwt_inverse()`

Inverse SGWT transform - reconstructs a graph signal from wavelet coefficients.

```r
sgwt_inverse(decomposition)
```

**Parameters:**
- `decomposition`: Output from `sgwt_forward()`

**Returns:**
- Reconstructed graph signal

### `sgwt_energy_analysis()`

Analyzes energy distribution across different scales in an SGWT decomposition.

```r
sgwt_energy_analysis(sgwt_result)
```

**Parameters:**
- `sgwt_result`: Output from `SGWT()` function

**Returns:**
- Data frame with scale-wise energy analysis

## Utility Functions

### `sgwt_similarity()`

Comprehensive function for computing similarity between two spatial signals.

```r
sgwt_similarity(signal1, signal2, data.in = NULL, 
                x_col = "x", y_col = "y", k = 25, J = 4, 
                kernel_type = "mexican_hat", eps = 1e-12, 
                low_only = FALSE, return_parts = FALSE, ...)
```

**Parameters:**
- `signal1`, `signal2`: Signal identifiers (column names) or pre-computed SGWT results
- `data.in`: Data frame containing the signals (if not pre-computed)
- `x_col`, `y_col`: Spatial coordinate column names
- `k`: Number of nearest neighbors for graph construction
- `J`: Number of wavelet scales
- `kernel_type`: Wavelet kernel type
- `eps`: Numerical stability parameter
- `low_only`: If TRUE, compute only low-frequency similarity
- `return_parts`: If TRUE, return detailed breakdown

**Returns:**
- Similarity score (scalar) or detailed analysis (list)

### `sgwt_weighted_similarity()`

Energy-normalized weighted similarity between two SGWT results.

```r
sgwt_weighted_similarity(sgwt_a, sgwt_b, eps = 1e-12, 
                        validate = TRUE, return_parts = TRUE, 
                        low_only = FALSE)
```

**Parameters:**
- `sgwt_a`, `sgwt_b`: SGWT results to compare
- `eps`: Numerical stability parameter
- `validate`: Whether to validate input consistency
- `return_parts`: Return detailed breakdown vs scalar result
- `low_only`: Compute only low-frequency similarity

**Returns:**
- If `return_parts=TRUE`: List with detailed similarity analysis
- If `return_parts=FALSE`: Scalar similarity score

### `cosine_similarity()`

Robust cosine similarity implementation with numerical stability.

```r
cosine_similarity(x, y, eps = 1e-12)
```

**Parameters:**
- `x`, `y`: Numeric vectors to compare
- `eps`: Numerical stability threshold

**Returns:**
- Cosine similarity value in [-1, 1]

### `Cal_Eigen()`

Eigenvalue analysis with automatic knee detection for determining the number of low-frequency components.

```r
Cal_Eigen(data.in, x_col = "x", y_col = "y", k = 25, k_fold = 15)
```

**Parameters:**
- `data.in`: Data frame with spatial coordinates
- `x_col`, `y_col`: Coordinate column names
- `k`: Number of nearest neighbors
- `k_fold`: Parameter for knee detection

**Returns:**
- List containing knee point and eigenvectors

### `cal_laplacian()`

Calculate graph Laplacian matrices from spatial data.

```r
cal_laplacian(data.in, x_col = "x", y_col = "y", k = 25, 
              laplacian_type = "normalized")
```

**Parameters:**
- `data.in`: Data frame with spatial coordinates
- `x_col`, `y_col`: Coordinate column names
- `k`: Number of nearest neighbors
- `laplacian_type`: "normalized" or "unnormalized"

**Returns:**
- List with adjacency matrix, Laplacian matrix, and spatial coordinates

### `gft()`

Graph Fourier Transform implementation.

```r
gft(signal, eigenvectors, inverse = FALSE)
```

**Parameters:**
- `signal`: Graph signal to transform
- `eigenvectors`: Graph Laplacian eigenvectors
- `inverse`: If TRUE, perform inverse GFT

**Returns:**
- Transformed signal

### `Cal_GCC()` (Deprecated)

**Note**: This function is deprecated. Use `sgwt_similarity()` instead.

```r
Cal_GCC(data.in = NULL, knee = NULL, signal1 = NULL, signal2 = NULL, 
        eigenvector = NULL)
```

## Visualization Functions

### `plot_sgwt_decomposition()`

Comprehensive visualization of SGWT decomposition results.

```r
plot_sgwt_decomposition(sgwt_result, data, point_size = 1.5, 
                       ncol = NULL, show_colorbar = TRUE)
```

**Parameters:**
- `sgwt_result`: Output from `SGWT()` function
- `data`: Original data frame with coordinates
- `point_size`: Size of points in plots
- `ncol`: Number of columns in plot layout
- `show_colorbar`: Whether to show color bars

**Returns:**
- Combined ggplot object with all decomposition components

### `demo_sgwt()`

Interactive demonstration function with synthetic data.

```r
demo_sgwt(n_points = 100, noise_level = 0.1, k = 15, J = 4)
```

**Parameters:**
- `n_points`: Number of spatial points to generate
- `noise_level`: Amount of noise to add to synthetic signal
- `k`: Number of nearest neighbors
- `J`: Number of wavelet scales

**Returns:**
- SGWT analysis results and generates visualization plots

## Kernel Functions

### `sgwt_kernel()`

Wavelet kernel functions for SGWT.

```r
sgwt_kernel(x, kernel_type = "mexican_hat")
```

**Parameters:**
- `x`: Input values (typically scaled eigenvalues)
- `kernel_type`: "mexican_hat" or "meyer"

**Returns:**
- Kernel values

### `sgwt_scaling_kernel()`

Scaling function kernel for low-pass filtering.

```r
sgwt_scaling_kernel(x, kernel_type = "mexican_hat")
```

**Parameters:**
- `x`: Input values
- `kernel_type`: Kernel type

**Returns:**
- Scaling kernel values

### `compute_sgwt_filters()`

Compute complete filter bank for SGWT analysis.

```r
compute_sgwt_filters(eigenvalues, scales, kernel_type = "mexican_hat")
```

**Parameters:**
- `eigenvalues`: Graph Laplacian eigenvalues
- `scales`: Scale parameters
- `kernel_type`: Wavelet kernel type

**Returns:**
- List of filter responses for each scale

### `sgwt_auto_scales()`

Automatic scale generation based on eigenvalue spectrum.

```r
sgwt_auto_scales(eigenvalues, J = 5, scaling_factor = 2)
```

**Parameters:**
- `eigenvalues`: Graph Laplacian eigenvalues
- `J`: Number of scales to generate
- `scaling_factor`: Geometric progression factor

**Returns:**
- Vector of automatically determined scales

## Mathematical Framework

### Graph Signal Processing Foundations

The mathematical foundation of BioGSP is based on extending classical signal processing to graph domains.

#### Graph Definition
A graph G = (V, E) consists of:
- **V**: Set of vertices (spatial locations)
- **E**: Set of edges (connections between nearby locations)
- **W**: Weighted adjacency matrix

#### Graph Laplacian
The normalized graph Laplacian is defined as:
```
L = I - D^(-1/2) W D^(-1/2)
```
where D is the degree matrix.

#### Graph Fourier Transform
For a signal f on graph G:
```
f̂(λᵢ) = ⟨f, χᵢ⟩ = Σⱼ f(j) χᵢ(j)
```
where χᵢ are the eigenvectors of L with eigenvalues λᵢ.

### Spectral Graph Wavelets

#### Wavelet Definition
Wavelets are defined in the spectral domain:
```
ψⱼ,ₙ = √(λₘₐₓ/s) g(λᵢ/s) χᵢ(n)
```
where g is the wavelet kernel and s is the scale parameter.

#### Multi-scale Decomposition
The SGWT decomposes a signal into:
```
f = Σᵢ ⟨f, φᵢ⟩ φᵢ + ΣⱼΣᵢ ⟨f, ψⱼ,ᵢ⟩ ψⱼ,ᵢ
```
where φ are scaling functions and ψ are wavelets.

### SGCC Mathematical Framework

#### Energy Normalization
The energy-normalized similarity between signals f₁ and f₂ is:
```
S = (w_low × c_low + w_NL × c_nonlow)
```
where:
- `w_low = 0.5 × (E_low_1/(E_low_1 + E_NL_1) + E_low_2/(E_low_2 + E_NL_2))`
- `c_low = cosine_similarity(scaling_coeffs_1, scaling_coeffs_2)`
- `c_nonlow = cosine_similarity(wavelet_coeffs_1, wavelet_coeffs_2)`

#### Robustness Features
- **Numerical Stability**: Epsilon-based thresholding prevents division by zero
- **Value Clamping**: Results bounded to [-1, 1] range
- **Energy Weighting**: Automatic balancing based on signal energy distribution

## Implementation Details

### Performance Optimizations

1. **Sparse Matrix Operations**: Uses `Matrix` package for efficient sparse computations
2. **Fast Eigendecomposition**: Leverages `RSpectra` for large-scale eigenvalue problems
3. **Vectorized Operations**: Minimizes loops in favor of vectorized R operations
4. **Memory Management**: Efficient handling of large spatial datasets

### Numerical Stability

1. **Epsilon Thresholding**: Prevents division by zero in similarity computations
2. **Condition Number Checking**: Validates matrix conditioning before eigendecomposition
3. **Finite Value Validation**: Replaces non-finite values with appropriate defaults
4. **Scale Normalization**: Ensures numerical stability across different scales

### Error Handling

1. **Input Validation**: Comprehensive checking of function parameters
2. **Dimension Consistency**: Validates matching dimensions across inputs
3. **Graceful Degradation**: Provides informative error messages
4. **Recovery Mechanisms**: Automatic parameter adjustment for edge cases

## Advanced Examples

### Custom Scale Analysis

```r
library(BioGSP)

# Load your spatial data
data <- read.csv("your_spatial_data.csv")

# Define custom scales based on your domain knowledge
# Smaller scales capture fine details, larger scales capture broad patterns
custom_scales <- c(0.1, 0.3, 0.7, 1.5, 3.0)

# Run SGWT with custom scales
result <- SGWT(data.in = data,
               signal = "expression_level",
               x_col = "x_coord",
               y_col = "y_coord",
               scales = custom_scales,
               k = 20,
               kernel_type = "mexican_hat")

# Analyze energy distribution
energy_df <- sgwt_energy_analysis(result)
print(energy_df)

# Visualize scale-specific patterns
plots <- plot_sgwt_decomposition(result, data)
print(plots)
```

### Batch Similarity Analysis

```r
# Compare multiple signals against a reference
reference_signal <- "reference_gene"
test_signals <- c("gene1", "gene2", "gene3", "gene4")

# Compute similarities
similarities <- sapply(test_signals, function(signal) {
  sgwt_similarity(reference_signal, signal, 
                 data.in = spatial_data,
                 k = 25, J = 4)
})

# Create results data frame
results_df <- data.frame(
  gene = test_signals,
  similarity = similarities,
  stringsAsFactors = FALSE
)

# Sort by similarity
results_df <- results_df[order(results_df$similarity, decreasing = TRUE), ]
print(results_df)
```

### Energy-based Feature Selection

```r
# Function to extract energy features from SGWT
extract_energy_features <- function(data, signal_cols, k = 25, J = 4) {
  features <- matrix(0, nrow = length(signal_cols), ncol = J + 1)
  colnames(features) <- c("scaling", paste0("wavelet_", 1:J))
  rownames(features) <- signal_cols
  
  for(i in seq_along(signal_cols)) {
    # Run SGWT
    result <- SGWT(data.in = data, signal = signal_cols[i], k = k, J = J)
    
    # Extract energy features
    energy_analysis <- sgwt_energy_analysis(result)
    features[i, ] <- energy_analysis$energy_ratio
  }
  
  return(features)
}

# Extract features for multiple genes
gene_cols <- c("gene1", "gene2", "gene3", "gene4", "gene5")
energy_features <- extract_energy_features(spatial_data, gene_cols)

# Perform clustering based on energy features
library(cluster)
clusters <- kmeans(energy_features, centers = 3)
print(clusters$cluster)
```

### Cross-Modal Comparison

```r
# Compare different measurement modalities on the same spatial domain
# Example: RNA-seq vs protein expression

# Prepare data with both modalities
combined_data <- merge(rna_data, protein_data, by = c("x", "y"))

# Compare corresponding genes/proteins
gene_protein_pairs <- list(
  c("GENE1_rna", "GENE1_protein"),
  c("GENE2_rna", "GENE2_protein"),
  c("GENE3_rna", "GENE3_protein")
)

cross_modal_similarities <- sapply(gene_protein_pairs, function(pair) {
  sgwt_similarity(pair[1], pair[2], 
                 data.in = combined_data,
                 k = 30, J = 5)
})

names(cross_modal_similarities) <- sapply(gene_protein_pairs, function(x) x[1])
print(cross_modal_similarities)
```

## Performance Considerations

### Memory Usage

- **Large Datasets**: For >10,000 points, consider reducing `k` or using sparse operations
- **Multiple Scales**: Memory usage increases linearly with number of scales `J`
- **Batch Processing**: Process multiple signals sequentially to manage memory

### Computational Complexity

- **Graph Construction**: O(N log N) for k-NN graph construction
- **Eigendecomposition**: O(N³) for dense matrices, O(N²) for sparse
- **SGWT Transform**: O(N × J) for forward/inverse transforms
- **Similarity Computation**: O(N × J) for coefficient comparison

### Optimization Tips

1. **Parameter Selection**:
   - Start with `k = 15-25` for most applications
   - Use `J = 3-5` scales for initial analysis
   - Increase parameters only if needed

2. **Data Preprocessing**:
   - Remove outliers that might affect graph construction
   - Normalize signals to similar scales
   - Consider spatial subsampling for very dense data

3. **Parallel Processing**:
   - Use `parallel` package for batch similarity computations
   - Consider GPU acceleration for large-scale problems

## Contributing

### Development Setup

```bash
# Clone the repository
git clone https://github.com/BMEngineeR/BioGSP.git
cd BioGSP

# Install development dependencies
R -e "install.packages(c('devtools', 'testthat', 'roxygen2'))"

# Load package for development
R -e "devtools::load_all()"
```

### Code Style

- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use roxygen2 for documentation
- Include unit tests for new functions
- Ensure all examples run without errors

### Testing

```r
# Run all tests
devtools::test()

# Check package
devtools::check()

# Build vignettes
devtools::build_vignettes()
```

### Submitting Changes

1. Fork the repository
2. Create a feature branch
3. Make your changes with tests
4. Update documentation
5. Submit a pull request

---

For basic usage information, see the main [README.md](README.md).

For questions and support, please open an issue on [GitHub](https://github.com/BMEngineeR/BioGSP/issues).
