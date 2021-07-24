#' Handle BentoBox color scaling parameters
#'
#' \code{colorby} should be used to create a set of parameters
#' that specify color scaling for the functions \code{bb_plotPairs},
#' \code{bb_plotPairsArches}, and \code{bb_plotRanges}.
#'
#' @param column String specifying name of data column to scale colors by.
#' @param palette (optional) A function describing the color palette to use for
#' color scaling.
#' @param range (optional) A numeric vector specifying the range of values to
#' apply a color scale to.
#' @param scalePerRegion A logical value indicating whether to adjust 
#' NULL range of numerical `colorby` values to subset of data in a plotted
#' genomic regoin. Default value is \code{scalePerRegion = FALSE}.
#'
#' @return Returns a "\code{bb_colorby}" object.
#'
#' @examples
#' ## Load paired ranges data in BEDPE format
#' library(BentoBoxData)
#' data("IMR90_DNAloops_pairs")
#'
#' ## Add a length column
#' IMR90_DNAloops_pairs$length <- 
#'         (IMR90_DNAloops_pairs$start2 - IMR90_DNAloops_pairs$start1) / 1000
#'
#' ## Plot pairs with colorby object set for `length` column
#' bedpePlot <- bb_plotPairs(
#'     data = IMR90_DNAloops_pairs,
#'     chrom = "chr21",
#'     chromstart = 27900000, chromend = 30700000,
#'     assembly = "hg19",
#'     fill = colorby("length", palette = 
#'                 colorRampPalette(c("dodgerblue2", "firebrick2"))),
#'     lwd = 2, spaceHeight = .7,
#' )
#' @export
colorby <- function(column, palette = NULL, range = NULL, 
                    scalePerRegion = FALSE) {
    colorbyObject <- structure(list(column = column, palette = palette, 
                                    range = range, 
                                    scalePerRegion = scalePerRegion),
        class = "bb_colorby"
    )
    return(colorbyObject)
}
