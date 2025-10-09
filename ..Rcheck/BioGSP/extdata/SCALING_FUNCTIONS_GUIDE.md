# SGWT Scaling Functions Guide

This guide provides comprehensive information about the different scaling function options available in the Spectral Graph Wavelet Transform (SGWT) implementation.

## Overview

Scaling functions in SGWT serve as low-pass filters that capture the smooth, low-frequency components of signals on graphs. The choice of scaling function affects the localization properties, smoothness, and reconstruction quality of the wavelet transform.

## Available Scaling Functions

### 1. Gaussian (Default)
**Formula**: `h(λ) = exp(-λ²/(2σ²))`
- **Properties**: Smooth, infinite support, excellent localization in both vertex and spectral domains
- **Use case**: General-purpose applications, when smoothness is important
- **Characteristics**: Exponential decay, no compact support

```r
sgwt_scaling_kernel(x, scale_param = 1, scaling_type = "gaussian")
```

### 2. Exponential
**Formula**: `h(λ) = exp(-λ/σ)`
- **Properties**: Monotonic decay, infinite support, simple form
- **Use case**: When you want faster decay than Gaussian but still smooth
- **Characteristics**: Linear exponential decay, heavier tail than Gaussian

```r
sgwt_scaling_kernel(x, scale_param = 1, scaling_type = "exponential")
```

### 3. Polynomial
**Formula**: `h(λ) = (1 + (λ/σ)²)^(-p)` where p=2
- **Properties**: Algebraic decay, infinite support, heavy tails
- **Use case**: When you need polynomial decay behavior
- **Characteristics**: Slower decay than exponential functions

```r
sgwt_scaling_kernel(x, scale_param = 1, scaling_type = "polynomial")
```

### 4. Meyer
**Formula**: Piecewise function with smooth transitions
- **Properties**: Compact support, smooth transitions, bandlimited-like behavior
- **Use case**: When you need strict frequency band separation
- **Characteristics**: 
  - h(λ) = 1 for λ ≤ 0.5σ
  - h(λ) = 0 for λ ≥ σ
  - Smooth cosine transition in between

```r
sgwt_scaling_kernel(x, scale_param = 1, scaling_type = "meyer")
```

### 5. Spline
**Formula**: B-spline based function
- **Properties**: Compact support, piecewise polynomial, smooth
- **Use case**: When you need localized support with polynomial smoothness
- **Characteristics**: Cubic spline with support on [0, 2σ]

```r
sgwt_scaling_kernel(x, scale_param = 1, scaling_type = "spline")
```

### 6. Hann
**Formula**: `h(λ) = 0.5(1 + cos(πλ/σ))` for λ ≤ σ
- **Properties**: Compact support, cosine-based, smooth
- **Use case**: Signal processing applications, windowing
- **Characteristics**: Hann window shape, support on [0, σ]

```r
sgwt_scaling_kernel(x, scale_param = 1, scaling_type = "hann")
```

### 7. Tight Frame
**Formula**: Shifted Hann kernel designed for tight frame properties
- **Properties**: Ensures perfect reconstruction, tight frame property
- **Use case**: When perfect reconstruction is critical
- **Characteristics**: Designed to satisfy tight frame conditions

```r
sgwt_scaling_kernel(x, scale_param = 1, scaling_type = "tight_frame")
```

## Comparison of Properties

| Scaling Type | Support | Smoothness | Decay Rate | Tight Frame | Best For |
|--------------|---------|------------|------------|-------------|----------|
| Gaussian | Infinite | C∞ | Exponential | No | General use, smooth signals |
| Exponential | Infinite | C∞ | Exponential | No | Fast decay, simple form |
| Polynomial | Infinite | C∞ | Algebraic | No | Heavy-tailed signals |
| Meyer | Compact | C∞ | Sharp cutoff | No | Frequency separation |
| Spline | Compact | C² | Polynomial | No | Local support needed |
| Hann | Compact | C∞ | Cosine | No | Windowing applications |
| Tight Frame | Compact | C∞ | Cosine | Yes | Perfect reconstruction |

## Usage Examples

### Basic Usage
```r
# Load the SGWT functions
source("R/sgwt_core.R")

# Create sample data
eigenvals <- c(0, 0.1, 0.5, 1.0, 1.5, 2.0)
scales <- c(2, 1, 0.5)

# Use different scaling functions
filters_gaussian <- compute_sgwt_filters(eigenvals, scales, scaling_type = "gaussian")
filters_meyer <- compute_sgwt_filters(eigenvals, scales, scaling_type = "meyer")
filters_tight <- compute_sgwt_filters(eigenvals, scales, scaling_type = "tight_frame")
```

### Forward Transform with Different Scaling Functions
```r
# Assuming you have signal, eigenvectors, eigenvalues
result_gaussian <- sgwt_forward(signal, eigenvectors, eigenvalues, scales, 
                               scaling_type = "gaussian")
result_meyer <- sgwt_forward(signal, eigenvectors, eigenvalues, scales, 
                            scaling_type = "meyer")
```

### Comparing Scaling Functions
```r
# Visual comparison of all scaling functions
comparison <- compare_scaling_functions(x_range = c(0, 3), scale_param = 1)

# This will plot all scaling functions for comparison
# and return a data frame with the values
```

## Choosing the Right Scaling Function

### For General Applications
- **Gaussian**: Best all-around choice, smooth and well-localized
- **Exponential**: Good alternative when you want simpler form

### For Specific Requirements
- **Perfect Reconstruction**: Use `tight_frame`
- **Compact Support**: Use `meyer`, `spline`, or `hann`
- **Frequency Band Separation**: Use `meyer`
- **Heavy-tailed Signals**: Use `polynomial`

### For Signal Processing Applications
- **Windowing**: Use `hann`
- **Spectral Analysis**: Use `meyer` or `tight_frame`
- **Denoising**: Use `gaussian` or `exponential`

## Technical Notes

### Scale Parameter
The `scale_param` controls the width/spread of the scaling function:
- Larger values → wider support, more low-pass filtering
- Smaller values → narrower support, less smoothing

### Computational Considerations
- **Compact support functions** (`meyer`, `spline`, `hann`, `tight_frame`) are computationally more efficient
- **Infinite support functions** (`gaussian`, `exponential`, `polynomial`) provide better theoretical properties but may require truncation

### Reconstruction Quality
- **Tight frame** scaling functions guarantee perfect reconstruction
- Other functions may have small reconstruction errors depending on the graph and signal

## References

1. Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011). Wavelets on graphs via spectral graph theory. Applied and computational harmonic analysis, 30(2), 129-150.

2. Shuman, D. I., Wiesmeyr, C., Holighaus, N., & Vandergheynst, P. (2015). Spectrum-adapted tight graph wavelet and vertex-frequency frames. IEEE Transactions on Signal Processing, 63(16), 4223-4235.

3. Leonardi, N., & Van De Ville, D. (2013). Tight wavelet frames on multislice graphs. IEEE Transactions on Signal Processing, 61(13), 3357-3367.