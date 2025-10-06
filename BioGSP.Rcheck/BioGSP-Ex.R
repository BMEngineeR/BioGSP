pkgname <- "BioGSP"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "BioGSP-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('BioGSP')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BioGSP-package")
### * BioGSP-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BioGSP-package
### Title: BioGSP: Biological Graph Signal Processing for Spatial Data
###   Analysis
### Aliases: BioGSP-package _PACKAGE
### Keywords: biological-data graph-theory package spatial-analysis
###   wavelets

### ** Examples

## Not run: 
##D # Load the package
##D library(BioGSP)
##D 
##D # Run a quick demo
##D demo_result <- demo_sgwt()
##D 
##D # Generate synthetic data
##D set.seed(123)
##D n <- 100
##D data <- data.frame(
##D   x = runif(n, 0, 10),
##D   y = runif(n, 0, 10),
##D   signal = sin(runif(n, 0, 2*pi))
##D )
##D 
##D # New workflow: Initialize -> Build Graph -> Run SGWT
##D SG <- initSGWT(data, signals = "signal", k = 8, J = 4, kernel_type = "heat")
##D SG <- runSpecGraph(SG)
##D SG <- runSGWT(SG)
##D 
##D # Analyze results
##D energy_analysis <- sgwt_energy_analysis(SG)
##D print(energy_analysis)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BioGSP-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("FastDecompositionLap")
### * FastDecompositionLap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FastDecompositionLap
### Title: Fast eigendecomposition of Laplacian matrix
### Aliases: FastDecompositionLap

### ** Examples

## Not run: 
##D # Create a Laplacian matrix and decompose
##D L <- matrix(c(2, -1, -1, -1, 2, -1, -1, -1, 2), nrow = 3)
##D decomp <- FastDecompositionLap(L, k_eigen = 25)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FastDecompositionLap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cal_laplacian")
### * cal_laplacian

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cal_laplacian
### Title: Calculate Graph Laplacian Matrix
### Aliases: cal_laplacian

### ** Examples

## Not run: 
##D W <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
##D cal_laplacian(W, type = "normalized")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cal_laplacian", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("codex_toy_data")
### * codex_toy_data

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: codex_toy_data
### Title: Toy CODEX Spatial Cell Type Data
### Aliases: codex_toy_data
### Keywords: CODEX SGWT datasets spatial

### ** Examples

# Load the toy dataset
data(codex_toy_data)

# Examine the structure
str(codex_toy_data)
head(codex_toy_data)

# Summary of cell types
table(codex_toy_data$Annotation5)

# Summary by ROI
table(codex_toy_data$ROI_num)
table(codex_toy_data$ROI_num, codex_toy_data$Annotation5)

# Quick visualization of spatial distribution
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  ggplot(codex_toy_data, aes(x = X_cent, y = Y_cent, color = Annotation5)) +
    geom_point(size = 0.8, alpha = 0.7) +
    facet_wrap(~ROI_num, scales = "free") +
    labs(title = "Toy CODEX Spatial Cell Distribution by ROI",
         x = "X Coordinate", y = "Y Coordinate") +
    theme_minimal() +
    scale_y_reverse()
}

# Basic SGWT analysis example
## Not run: 
##D # Focus on BCL6- B Cell cells in ROI_1 for SGWT analysis
##D bcl6nb_data <- codex_toy_data[codex_toy_data$Annotation5 == "BCL6- B Cell" & 
##D                               codex_toy_data$ROI_num == "ROI_1", ]
##D 
##D # Create binned representation
##D library(dplyr)
##D binned_data <- codex_toy_data %>%
##D   filter(Annotation5 == "BCL6- B Cell", ROI_num == "ROI_1") %>%
##D   mutate(
##D     x_bin = cut(X_cent, breaks = 20, labels = FALSE),
##D     y_bin = cut(Y_cent, breaks = 20, labels = FALSE)
##D   ) %>%
##D   group_by(x_bin, y_bin) %>%
##D   summarise(cell_count = n(), .groups = 'drop')
##D 
##D # Prepare for SGWT
##D complete_grid <- expand.grid(x_bin = 1:20, y_bin = 1:20)
##D sgwt_data <- complete_grid %>%
##D   left_join(binned_data, by = c("x_bin", "y_bin")) %>%
##D   mutate(
##D     cell_count = ifelse(is.na(cell_count), 0, cell_count),
##D     x = x_bin,
##D     y = y_bin,
##D     signal = cell_count / max(cell_count, na.rm = TRUE)
##D   ) %>%
##D   select(x, y, signal)
##D 
##D # Apply SGWT
##D sgwt_result <- SGWT(data.in = sgwt_data,
##D                     signal = "signal",
##D                     k = 8,
##D                     J = 3,
##D                     kernel_type = "heat")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("codex_toy_data", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compare_kernel_families")
### * compare_kernel_families

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compare_kernel_families
### Title: Compare different kernel families
### Aliases: compare_kernel_families

### ** Examples

comparison <- compare_kernel_families()
comparison <- compare_kernel_families(x_range = c(0, 5), scale_param = 1.5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compare_kernel_families", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_sgwt_filters")
### * compute_sgwt_filters

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_sgwt_filters
### Title: Compute SGWT filters
### Aliases: compute_sgwt_filters

### ** Examples

eigenvals <- c(0, 0.1, 0.5, 1.0, 1.5)
scales <- c(2, 1, 0.5)
filters <- compute_sgwt_filters(eigenvals, scales)
filters_meyer <- compute_sgwt_filters(eigenvals, scales, kernel_type = "meyer")
filters_heat <- compute_sgwt_filters(eigenvals, scales, kernel_type = "heat")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_sgwt_filters", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cosine_similarity")
### * cosine_similarity

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cosine_similarity
### Title: Calculate cosine similarity between two vectors
### Aliases: cosine_similarity

### ** Examples

x <- c(1, 2, 3)
y <- c(2, 3, 4)
similarity <- cosine_similarity(x, y)
# With custom eps for numerical stability
similarity2 <- cosine_similarity(x, y, eps = 1e-10)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cosine_similarity", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("demo_sgwt")
### * demo_sgwt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: demo_sgwt
### Title: Demo function for SGWT
### Aliases: demo_sgwt

### ** Examples

## Not run: 
##D SG <- demo_sgwt()
##D print(SG)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("demo_sgwt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("find_knee_point")
### * find_knee_point

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: find_knee_point
### Title: Find knee point in a curve
### Aliases: find_knee_point

### ** Examples

y <- c(1, 2, 3, 10, 11, 12)  # curve with a knee
knee_idx <- find_knee_point(y)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("find_knee_point", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("hello_sgwt")
### * hello_sgwt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: hello_sgwt
### Title: Hello function for SGWT package demonstration
### Aliases: hello_sgwt

### ** Examples

hello_sgwt()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("hello_sgwt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("igft")
### * igft

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: igft
### Title: Inverse Graph Fourier Transform
### Aliases: igft

### ** Examples

## Not run: 
##D # Single signal
##D signal_reconstructed <- igft(fourier_coeffs, eigenvectors)
##D 
##D # Multiple signals (batch processing)
##D signals_reconstructed <- igft(fourier_coeffs_matrix, eigenvectors)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("igft", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("initSGWT")
### * initSGWT

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: initSGWT
### Title: Initialize SGWT object
### Aliases: initSGWT

### ** Examples

## Not run: 
##D # Initialize SGWT object
##D data <- data.frame(x = runif(100), y = runif(100), 
##D                   signal1 = rnorm(100), signal2 = rnorm(100))
##D SG <- initSGWT(data, signals = c("signal1", "signal2"))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("initSGWT", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("install_and_load")
### * install_and_load

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: install_and_load
### Title: Install and load packages
### Aliases: install_and_load

### ** Examples

## Not run: 
##D packages <- c("ggplot2" = "ggplot2", "devtools" = "r-lib/devtools")
##D install_and_load(packages)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("install_and_load", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_FM")
### * plot_FM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_FM
### Title: Plot Fourier modes (eigenvectors) from SGWT object
### Aliases: plot_FM

### ** Examples

## Not run: 
##D # Plot both low and high frequency modes
##D SG <- initSGWT(data) %>% runSpecGraph()
##D plot_FM(SG, mode_type = "both", n_modes = 4)
##D 
##D # Plot only low frequency modes
##D plot_FM(SG, mode_type = "low", n_modes = 8)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_FM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_sgwt_decomposition")
### * plot_sgwt_decomposition

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_sgwt_decomposition
### Title: Plot SGWT decomposition results
### Aliases: plot_sgwt_decomposition

### ** Examples

## Not run: 
##D # Assuming you have SGWT object
##D plots <- plot_sgwt_decomposition(SG_object, signal_name = "signal1")
##D print(plots)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_sgwt_decomposition", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("runSGCC")
### * runSGCC

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: runSGCC
### Title: Run SGCC weighted similarity analysis in Fourier domain
### Aliases: runSGCC

### ** Examples

## Not run: 
##D # Between two signals in same SGWT object
##D similarity <- runSGCC("signal1", "signal2", SG = SG_object)
##D 
##D # Between two SGWT objects
##D similarity <- runSGCC(SG_object1, SG_object2)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("runSGCC", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("runSGWT")
### * runSGWT

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: runSGWT
### Title: Run SGWT forward and inverse transforms for all signals
### Aliases: runSGWT

### ** Examples

## Not run: 
##D SG <- initSGWT(data)
##D SG <- runSpecGraph(SG)
##D SG <- runSGWT(SG)  # Uses batch processing by default
##D SG <- runSGWT(SG, use_batch = FALSE)  # Force individual processing
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("runSGWT", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("runSpecGraph")
### * runSpecGraph

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: runSpecGraph
### Title: Build spectral graph for SGWT object
### Aliases: runSpecGraph

### ** Examples

## Not run: 
##D SG <- initSGWT(data)
##D # Uses full length by default
##D SG <- runSpecGraph(SG, k = 30, laplacian_type = "normalized")  
##D # Or specify custom length
##D SG <- runSpecGraph(SG, k = 30, laplacian_type = "normalized", 
##D                    length_eigenvalue = 30)  
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("runSpecGraph", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sgwt_auto_scales")
### * sgwt_auto_scales

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sgwt_auto_scales
### Title: Generate automatic scales for SGWT
### Aliases: sgwt_auto_scales

### ** Examples

scales <- sgwt_auto_scales(lmax = 2.0, J = 5, scaling_factor = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sgwt_auto_scales", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sgwt_energy_analysis")
### * sgwt_energy_analysis

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sgwt_energy_analysis
### Title: Analyze SGWT energy distribution across scales in Fourier domain
### Aliases: sgwt_energy_analysis

### ** Examples

## Not run: 
##D # Assuming you have SGWT object
##D energy_analysis <- sgwt_energy_analysis(SG_object, signal_name = "signal1")
##D print(energy_analysis)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sgwt_energy_analysis", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sgwt_forward")
### * sgwt_forward

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sgwt_forward
### Title: Forward SGWT transform (single or batch)
### Aliases: sgwt_forward

### ** Examples

## Not run: 
##D # Single signal
##D result <- sgwt_forward(signal, eigenvectors, eigenvalues, scales)
##D 
##D # Multiple signals (batch processing)
##D signals_matrix <- cbind(signal1, signal2, signal3)
##D result <- sgwt_forward(signals_matrix, eigenvectors, eigenvalues, scales)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sgwt_forward", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sgwt_inverse")
### * sgwt_inverse

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sgwt_inverse
### Title: Inverse SGWT transform (single or batch)
### Aliases: sgwt_inverse

### ** Examples

## Not run: 
##D # Single signal
##D inverse_result <- sgwt_inverse(sgwt_decomp, eigenvectors, original_signal)
##D 
##D # Multiple signals (batch processing)
##D inverse_result <- sgwt_inverse(sgwt_decomp, eigenvectors, original_signals_matrix)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sgwt_inverse", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_checkerboard")
### * simulate_checkerboard

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_checkerboard
### Title: Simulate checkerboard pattern
### Aliases: simulate_checkerboard

### ** Examples

## Not run: 
##D # Generate 8x8 checkerboard with 10x10 pixel tiles
##D df <- simulate_checkerboard(grid_size = 8, tile_size = 10)
##D p <- visualize_checkerboard(df)
##D print(p)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_checkerboard", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_multiscale")
### * simulate_multiscale

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_multiscale
### Title: Simulate Multiple Center Patterns
### Aliases: simulate_multiscale

### ** Examples

## Not run: 
##D # Generate multi-center patterns with default parameters
##D patterns <- simulate_multiscale()
##D 
##D # Custom parameters
##D Ra_seq <- seq(from = 10, to = 3, length.out = 6)
##D Rb_seq <- seq(from = 20, to = 3, length.out = 6)
##D patterns <- simulate_multiscale(Ra_seq = Ra_seq, Rb_seq = Rb_seq, n_centers = 3)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_multiscale", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_multiscale_overlap")
### * simulate_multiscale_overlap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_multiscale_overlap
### Title: Simulate Multiple Center Patterns with Fixed Centers
### Aliases: simulate_multiscale_overlap

### ** Examples

## Not run: 
##D # Generate multi-center patterns with fixed centers
##D patterns <- simulate_multiscale_overlap()
##D 
##D # Single center at grid center
##D patterns_single <- simulate_multiscale_overlap(n_centers = 1)
##D 
##D # Custom parameters with multiple centers
##D Ra_seq <- seq(from = 10, to = 3, length.out = 4)
##D Rb_seq <- seq(from = 15, to = 2, length.out = 4)
##D patterns <- simulate_multiscale_overlap(
##D   Ra_seq = Ra_seq, 
##D   Rb_seq = Rb_seq, 
##D   n_centers = 2,
##D   seed = 456
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_multiscale_overlap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("simulate_stripe_patterns")
### * simulate_stripe_patterns

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: simulate_stripe_patterns
### Title: Simulate Stripe Patterns
### Aliases: simulate_stripe_patterns

### ** Examples

## Not run: 
##D # Generate stripe patterns with default parameters
##D patterns <- simulate_stripe_patterns()
##D 
##D # Custom parameters
##D patterns <- simulate_stripe_patterns(
##D   grid_size = 80,
##D   gap_seq = c(10, 20),
##D   width_seq = c(5, 10, 20),
##D   theta_seq = c(0, 30, 60),
##D   eps = 1e-9,
##D   verbose = TRUE
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("simulate_stripe_patterns", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("visualize_checkerboard")
### * visualize_checkerboard

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: visualize_checkerboard
### Title: Visualize checkerboard pattern
### Aliases: visualize_checkerboard

### ** Examples

## Not run: 
##D df <- simulate_checkerboard(grid_size = 6, tile_size = 5)
##D p <- visualize_checkerboard(df, color1 = "darkblue", color2 = "lightgray")
##D print(p)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("visualize_checkerboard", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("visualize_multiscale")
### * visualize_multiscale

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: visualize_multiscale
### Title: Visualize Multiple Center Simulation Results
### Aliases: visualize_multiscale

### ** Examples

## Not run: 
##D # Generate and visualize patterns
##D Ra_seq <- seq(from = 10, to = 3, length.out = 6)
##D Rb_seq <- seq(from = 20, to = 3, length.out = 6)
##D sim_data <- simulate_multiscale(Ra_seq = Ra_seq, Rb_seq = Rb_seq, n_centers = 3)
##D plot_grid <- visualize_multiscale(sim_data, Ra_seq, Rb_seq)
##D print(plot_grid)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("visualize_multiscale", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("visualize_sgwt_kernels")
### * visualize_sgwt_kernels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: visualize_sgwt_kernels
### Title: Visualize SGWT kernels and scaling functions
### Aliases: visualize_sgwt_kernels

### ** Examples

## Not run: 
##D # Generate some example eigenvalues
##D eigenvals <- seq(0, 2, length.out = 100)
##D 
##D # Visualize kernels with specific parameters
##D viz_result <- visualize_sgwt_kernels(
##D   eigenvalues = eigenvals,
##D   J = 4,
##D   scaling_factor = 2,
##D   kernel_type = "heat"
##D )
##D print(viz_result$plot)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("visualize_sgwt_kernels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("visualize_similarity_xy")
### * visualize_similarity_xy

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: visualize_similarity_xy
### Title: Visualize similarity in low vs non-low frequency space
### Aliases: visualize_similarity_xy

### ** Examples

## Not run: 
##D # Single similarity result
##D sim_result <- runSGCC("signal1", "signal2", SG = SG_object)
##D plot <- visualize_similarity_xy(sim_result)
##D print(plot)
##D 
##D # Multiple similarity results
##D sim_results <- list(
##D   pair1 = runSGCC("signal1", "signal2", SG = SG_object1),
##D   pair2 = runSGCC("signal1", "signal2", SG = SG_object2)
##D )
##D plot <- visualize_similarity_xy(sim_results, show_names = TRUE)
##D print(plot)
##D 
##D # Show both labels and names (for comparison)
##D plot_both <- visualize_similarity_xy(sim_results, show_labels = TRUE, show_names = TRUE)
##D print(plot_both)
##D 
##D # With many data points (>50), names will be randomly sampled
##D # install.packages("ggrepel")  # Required for show_names = TRUE
##D plot_many <- visualize_similarity_xy(many_sim_results, show_names = TRUE)
##D print(plot_many)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("visualize_similarity_xy", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("visualize_stripe_patterns")
### * visualize_stripe_patterns

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: visualize_stripe_patterns
### Title: Visualize Stripe Pattern Simulation Results
### Aliases: visualize_stripe_patterns

### ** Examples

## Not run: 
##D # Generate and visualize patterns
##D sim_data <- simulate_stripe_patterns(
##D   grid_size = 80,
##D   gap_seq = c(10, 20),
##D   width_seq = c(5, 10, 20),
##D   theta_seq = c(0, 30, 60)
##D )
##D plot_grid <- visualize_stripe_patterns(sim_data, 
##D                                        gap_seq = c(10, 20),
##D                                        width_seq = c(5, 10, 20),
##D                                        theta_seq = c(0, 30, 60))
##D print(plot_grid)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("visualize_stripe_patterns", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
