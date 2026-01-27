#' Toy CODEX Spatial Cell Type Data
#'
#' @name codex_toy_data
#' @docType data
#'
#' @description A synthetic dataset mimicking CODEX multiplexed imaging data 
#' for demonstrating Spectral Graph Wavelet Transform (SGWT) analysis on 
#' spatial cell type distributions. The dataset contains spatial coordinates 
#' and cell type annotations for multiple immune cell populations arranged 
#' in realistic spatial clusters.
#'
#' @format A data frame with 18604 rows and 5 columns:
#' \describe{
#'   \item{cellLabel}{Character. Unique identifier for each cell}
#'   \item{Y_cent}{Numeric. Y coordinate of cell centroid (0-115 range)}
#'   \item{X_cent}{Numeric. X coordinate of cell centroid (0-116 range)}
#'   \item{Annotation5}{Character. Full descriptive cell type name}
#'   \item{ROI_num}{Character. Region of interest identifier ("ROI_0" through "ROI_15")}
#' }
#'
#' @details 
#' The dataset contains 16 regions of interest (ROI_0 through ROI_15) with different spatial patterns
#' and varying cell counts (945-1497 cells per ROI). Each ROI represents a distinct tissue region
#' with unique spatial arrangements of the same cell types.
#' 
#' ROI Distribution:
#' \itemize{
#'   \item \strong{ROI_0}: 952 cells
#'   \item \strong{ROI_1}: 945 cells  
#'   \item \strong{ROI_2}: 1155 cells
#'   \item \strong{ROI_3}: 1421 cells
#'   \item \strong{ROI_4}: 1096 cells
#'   \item \strong{ROI_5}: 1420 cells
#'   \item \strong{ROI_6-ROI_15}: 958-1497 cells each
#' }
#' 
#' Cell types across all ROIs include:
#' \itemize{
#'   \item \strong{BCL6- B Cell} (~3719 cells): Primary B cell population
#'   \item \strong{CD4 T} (~4092 cells): Helper T cells - largest population
#'   \item \strong{CD8 T} (~3346 cells): Cytotoxic T cells
#'   \item \strong{DC} (~2233 cells): Dendritic cells
#'   \item \strong{M1} (~1490 cells): M1 macrophages
#'   \item \strong{CD4 Treg} (~1490 cells): Regulatory T cells
#'   \item \strong{BCL6+ B Cell} (~931 cells): Activated B cells
#'   \item \strong{Endothelial} (~746 cells): Vascular cells
#'   \item \strong{M2} (~370 cells): M2 macrophages
#'   \item \strong{Myeloid} (~186 cells): Other myeloid cells
#'   \item \strong{Other} (~1 cells): Miscellaneous cell types
#' }
#'
#' This synthetic data is designed to demonstrate:
#' \itemize{
#'   \item Spatial clustering patterns of different cell types
#'   \item Multi-scale spatial analysis using SGWT
#'   \item Cross-cell type correlation analysis
#'   \item Graph construction and eigenvalue analysis
#'   \item Wavelet decomposition of spatial signals
#' }
#'
#' @source Generated synthetically using clustered normal distributions with 
#' realistic parameters based on real CODEX data characteristics.
#'
#' @usage data(codex_toy_data)
#'
#' @examples
#' # Load the toy dataset
#' data(codex_toy_data)
#' 
#' # Examine the structure
#' str(codex_toy_data)
#' head(codex_toy_data)
#' 
#' # Summary of cell types
#' table(codex_toy_data$Annotation5)
#' 
#' # Summary by ROI
#' table(codex_toy_data$ROI_num)
#' table(codex_toy_data$ROI_num, codex_toy_data$Annotation5)
#' 
#' # Quick visualization of spatial distribution
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   ggplot(codex_toy_data, aes(x = X_cent, y = Y_cent, color = Annotation5)) +
#'     geom_point(size = 0.8, alpha = 0.7) +
#'     facet_wrap(~ROI_num, scales = "free") +
#'     labs(title = "Toy CODEX Spatial Cell Distribution by ROI",
#'          x = "X Coordinate", y = "Y Coordinate") +
#'     theme_minimal() +
#'     scale_y_reverse()
#' }
#' 
#' # Basic SGWT analysis example
#' \donttest{
#' # Focus on BCL6- B Cell cells in ROI_1 for SGWT analysis
#' bcl6nb_data <- codex_toy_data[codex_toy_data$Annotation5 == "BCL6- B Cell" & 
#'                               codex_toy_data$ROI_num == "ROI_1", ]
#' 
#' # Create binned representation
#' library(dplyr)
#' binned_data <- codex_toy_data %>%
#'   filter(Annotation5 == "BCL6- B Cell", ROI_num == "ROI_1") %>%
#'   mutate(
#'     x_bin = cut(X_cent, breaks = 20, labels = FALSE),
#'     y_bin = cut(Y_cent, breaks = 20, labels = FALSE)
#'   ) %>%
#'   group_by(x_bin, y_bin) %>%
#'   summarise(cell_count = n(), .groups = 'drop')
#' 
#' # Prepare for SGWT
#' complete_grid <- expand.grid(x_bin = 1:20, y_bin = 1:20)
#' sgwt_data <- complete_grid %>%
#'   left_join(binned_data, by = c("x_bin", "y_bin")) %>%
#'   mutate(
#'     cell_count = ifelse(is.na(cell_count), 0, cell_count),
#'     x = x_bin,
#'     y = y_bin,
#'     signal = cell_count / max(cell_count, na.rm = TRUE)
#'   ) %>%
#'   select(x, y, signal)
#' 
#' # Apply SGWT using new workflow
#' SG <- initSGWT(sgwt_data, signals = "signal", J = 3, kernel_type = "heat")
#' SG <- runSpecGraph(SG, k = 8)
#' SG <- runSGWT(SG)
#' 
#' # View results
#' print(SG)
#' }
#'
#' @keywords datasets spatial CODEX SGWT
NULL