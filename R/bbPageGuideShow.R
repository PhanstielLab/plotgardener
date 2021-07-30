#' Reshow guides drawn with \code{bbPageCreate},
#' \code{bbPageGuideHorizontal}, and \code{bbPageGuideVertical}
#'
#' @seealso \link[BentoBox]{bbPageCreate},
#' \link[BentoBox]{bbPageGuideHorizontal},
#' \link[BentoBox]{bbPageGuideVertical}
#'
#' @return None.
#'
#' @examples
#' ## Load Hi-C data
#' library(BentoBoxData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create a page
#' bbPageCreate(width = 3, height = 3, default.units = "inches")
#'
#' ## Add a page guide
#' bbPageGuideHorizontal(y = 0.5, default.units = "inches")
#'
#' ## Plot and place Hi-C plot
#' hicPlot <- bbPlotHicSquare(
#'     data = IMR90_HiC_10kb, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     x = 0.5, y = 0.5, width = 2, height = 2,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Hide page guides
#' bbPageGuideHide()
#'
#' ## Re-show page guides
#' bbPageGuideShow()
#'
#' ## Annotate genome label
#' bbAnnoGenomeLabel(
#'     plot = hicPlot, scale = "Mb", axis = "x",
#'     x = 0.5, y = 2.53, just = c("left", "top")
#' )
#' @export
bbPageGuideShow <- function() {
    if (length(get("guide_grobs", envir = bbEnv)$children) == 0) {
        stop("No BentoBox page guides were previously drawn.")
    }

    add_guideGrobs <- function(child) {
        grid.add("guide_grobs", child = child, redraw = FALSE)
    }

    ## get the guide grobs
    grobList <- get("guide_grobs", envir = bbEnv)$children

    ## Get the last grob
    last_grob <- grobList[length(grobList)][[1]]

    ## Remove the last grob from the list of grobs
    grobList <- grobList[-length(grobList)]

    ## Re-add grobs to guide_grobs gTree without re-drawing
    invisible(lapply(grobList, add_guideGrobs))

    ## Re-add last grob, this time redrawing
    grid.add("guide_grobs", child = last_grob)
}
