#' BioGSP: Biological Graph Signal Processing for Spatial Data Analysis
#'
#' @description 
#' The BioGSP package provides a comprehensive implementation of Graph Signal Processing 
#' (GSP) methods including Spectral Graph Wavelet Transform (SGWT) for analyzing spatial 
#' patterns in biological data. This implementation is based on Hammond, Vandergheynst, 
#' and Gribonval (2011) "Wavelets on Graphs via Spectral Graph Theory".
#'
#' @details
#' The package enables multi-scale analysis of spatial signals by:
#' \itemize{
#'   \item Building graphs from spatial coordinates using k-nearest neighbors
#'   \item Computing graph Laplacian eigendecomposition for spectral analysis
#'   \item Designing wavelets in the spectral domain using various kernel functions
#'   \item Decomposing signals into scaling and wavelet components at multiple scales
#'   \item Providing reconstruction capabilities with error analysis
#'   \item Offering comprehensive visualization and analysis tools
#' }
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{initSGWT}}}{Initialize SGWT object with data and parameters}
#'   \item{\code{\link{runSpecGraph}}}{Build graph and compute eigendecomposition}
#'   \item{\code{\link{runSGWT}}}{Perform forward and inverse SGWT transforms}
#'   \item{\code{\link{runSGCC}}}{Calculate weighted similarity between signals}
#'   \item{\code{\link{sgwt_forward}}}{Forward SGWT transform}
#'   \item{\code{\link{sgwt_inverse}}}{Inverse SGWT transform}
#'   \item{\code{\link{sgwt_energy_analysis}}}{Energy distribution analysis}
#'   \item{\code{\link{plot_sgwt_decomposition}}}{Visualization of SGWT components}
#'   \item{\code{\link{demo_sgwt}}}{Demonstration with synthetic data}
#' }
#'
#' @section Applications:
#' The BioGSP package is particularly useful for:
#' \itemize{
#'   \item Spatial biology: Analyzing cell distribution patterns in tissue imaging (CODEX, Visium, etc.)
#'   \item Single-cell genomics: Spatial transcriptomics and proteomics analysis
#'   \item Neuroscience: Brain connectivity and signal analysis
#'   \item Pathology: Tumor microenvironment and tissue architecture analysis
#'   \item Developmental biology: Spatial pattern formation and cell fate mapping
#'   \item Immunology: Immune cell spatial organization and interactions
#' }
#'
#' @author BioGSP Development Team
#'
#' @references
#' Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011).
#' Wavelets on graphs via spectral graph theory.
#' Applied and Computational Harmonic Analysis, 30(2), 129-150.
#'
#' @keywords package spatial-analysis wavelets graph-theory biological-data
#' @aliases BioGSP-package
#'
#' @examples
#' \dontrun{
#' # Load the package
#' library(BioGSP)
#' 
#' # Run a quick demo
#' demo_result <- demo_sgwt()
#' 
#' # Generate synthetic data
#' set.seed(123)
#' n <- 100
#' data <- data.frame(
#'   x = runif(n, 0, 10),
#'   y = runif(n, 0, 10),
#'   signal = sin(runif(n, 0, 2*pi))
#' )
#' 
#' # New workflow: Initialize -> Build Graph -> Run SGWT
#' SG <- initSGWT(data, signals = "signal", k = 8, J = 4, kernel_type = "heat")
#' SG <- runSpecGraph(SG)
#' SG <- runSGWT(SG)
#' 
#' # Analyze results
#' energy_analysis <- sgwt_energy_analysis(SG)
#' print(energy_analysis)
#' }
#'
"_PACKAGE" 