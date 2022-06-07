#' plotgardener: Coordinate-based Genomic Visualization Package for R
#'
#' plotgardener is a coordinate-based genomic visualization package for R.
#' It grants users the ability to programmatically produce complex,
#' multi-paneled figures. Tailored for genomics, plotgardener allows users
#' to visualize large complex genomic datasets and
#' provides exquisite control over how plots are placed and arranged on
#' a page.
#' 
#' @author 
#' \strong{Maintainer}: Nicole Kramer \email{nekramer@live.unc.edu}
#' (\href{https://orcid.org/0000-0001-9617-9671}{ORCID})
#' 
#' Authors:
#' \itemize{
#'     \item Eric S. Davis \email{esdavis@live.unc.edu}
#'     (\href{https://orcid.org/0000-0003-4051-3217}{ORCID})
#'     \item Craig Wenger \email{craig.wenger@gmail.com}
#'     (\href{https://orcid.org/0000-0002-7361-8456}{ORCID})
#'     \item Douglas H. Phanstiel \email{douglas_phanstiel@med.unc.edu}
#'     [copyright holder]
#' }
#' 
#' Other contributors:
#' \itemize{
#'     \item Sarah Parker \email{sarmae@live.unc.edu} [contributor]
#'     \item Erika Deoudes \email{ed@erikadudes.com} [artist]
#'     \item Michael Love \email{milove@email.unc.edu} [contributor]
#' }
#' 
#' @seealso 
#' Useful links:
#' \itemize{
#'     \item \url{https://phanstiellab.github.io/plotgardener}
#'     \item \url{https://github.com/PhanstielLab/plotgardener}
#' }
#' @docType package
#' @name plotgardener
#' @useDynLib plotgardener
#' @importFrom grDevices colorRampPalette
#' @importFrom curl has_internet
#' @importFrom methods hasArg
#' @importFrom methods is
#' @importFrom stats na.omit
#' @importFrom utils str
#' @importFrom Rcpp sourceCpp
#' @importFrom data.table fread
#' @importFrom data.table setDT
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
#' @importFrom purrr pmap
#' @import grid
#' @importFrom strawr straw
#' @importFrom strawr readHicBpResolutions
#' @importFrom IRanges IRanges
#' @importFrom IRanges subsetByOverlaps
"_PACKAGE"
