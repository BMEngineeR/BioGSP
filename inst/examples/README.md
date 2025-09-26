# BioGSP Examples

This directory contains example scripts and demonstrations for the BioGSP package.

## Files

### simulation_demo.R
A comprehensive demonstration script showing:
- Multiple center pattern simulation
- Concentric ring pattern simulation  
- SGWT analysis with flexible column naming (X, Y, signal_1, signal_2)
- Energy distribution analysis
- Graph Cross-Correlation analysis
- Kernel type comparison

**Usage:**
```r
source("inst/examples/simulation_demo.R")
```

## Vignettes

For more detailed tutorials, see the vignettes:
- `vignettes/sgwt_simulation_demo.Rmd`: Complete RMarkdown tutorial with visualizations

## Key Features Demonstrated

1. **Flexible Column Naming**: All functions now support custom column names
   - X coordinates: any column name (e.g., "X", "x", "longitude")
   - Y coordinates: any column name (e.g., "Y", "y", "latitude") 
   - Signals: any column names (e.g., "signal_1", "signal_2", "gene_expression")

2. **Simulation Functions**:
   - `simulate_multiscale()`: Multiple center patterns with inner/outer regions
   - `simulate_ringpattern()`: Dynamic concentric ring patterns
   - `visualize_multiscale()`: Visualization for multiple center patterns
   - `visualize_ringpattern()`: Visualization for ring patterns

3. **Enhanced SGWT Functions**:
   - `SGWT()`: Main analysis function with x_col, y_col parameters
   - `Cal_Eigen()`: Eigendecomposition with flexible coordinate columns
   - `plot_sgwt_decomposition()`: Visualization with flexible column support

4. **Analysis Tools**:
   - Multiple kernel types: Mexican Hat, Meyer, Heat
   - Energy distribution analysis across scales
   - Graph Cross-Correlation between signals
   - Comprehensive visualization capabilities

