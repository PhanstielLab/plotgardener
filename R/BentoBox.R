#' BentoBox: Coordinate-based Genomic Visualization Package for R
#'
#' BentoBox is a coordinate-based genomic visualization package for R. It grants users
#' the ability to programmatically produce complex, multi-paneled figures. Tailored for
#' genomics, BentoBox allows users to visualize large complex genomic datasets and
#' provides exquisite control over how plots are placed and arranged on a page.
#'
#'
#' @docType package
#' @name BentoBox
#' @useDynLib BentoBox
#' @importFrom grDevices colorRampPalette
#' @importFrom methods hasArg
#' @importFrom stats na.omit
#' @importFrom utils str
#' @importFrom Rcpp sourceCpp
#' @importFrom data.table fread
#' @importFrom data.table setDT
#' @importFrom Rmpfr mpfr
#' @importFrom ggplotify as.grob
#' @importFrom ggplotify base2grob
#' @importFrom RColorBrewer brewer.pal
#' @importFrom tools file_ext
#' @importFrom dplyr anti_join
#' @importFrom dplyr left_join
#' @importFrom dplyr inner_join
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom purrr map
#' @import grid
"_PACKAGE"
