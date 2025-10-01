#' Global variables used in ggplot2 aesthetics
#' 
#' @name sgwt-globals
#' @description This file declares global variables used in ggplot2 aesthetics
#' to avoid R CMD check NOTEs about undefined global functions or variables.
#' 
#' @keywords internal
utils::globalVariables(c(
  "x", "y", "original", "scaling", "reconstructed", 
  "eigenvalue", "filter_value", "filter_type", "df_hex_combine",
  "X", "Y", "label"
))