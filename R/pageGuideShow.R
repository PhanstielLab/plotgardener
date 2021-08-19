#' Reshow guides drawn with \code{pageCreate},
#' \code{pageGuideHorizontal}, and \code{pageGuideVertical}
#'
#' @seealso \link[plotgardener]{pageCreate},
#' \link[plotgardener]{pageGuideHorizontal},
#' \link[plotgardener]{pageGuideVertical}
#'
#' @return None.
#'
#' @examples
#' ## Load Hi-C data
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create a page
#' pageCreate(width = 3, height = 3, default.units = "inches")
#'
#' ## Add a page guide
#' pageGuideHorizontal(y = 0.5, default.units = "inches")
#'
#' ## Plot and place Hi-C plot
#' hicPlot <- plotHicSquare(
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
#' pageGuideHide()
#'
#' ## Re-show page guides
#' pageGuideShow()
#'
#' ## Annotate genome label
#' annoGenomeLabel(
#'     plot = hicPlot, scale = "Mb", axis = "x",
#'     x = 0.5, y = 2.53, just = c("left", "top")
#' )
#' @export
pageGuideShow <- function() {
    if (length(get("guide_grobs", envir = pgEnv)$children) == 0) {
        stop("No `plotgardener` page guides were previously drawn.")
    }

    add_guideGrobs <- function(child) {
        grid.add("guide_grobs", child = child, redraw = FALSE)
    }

    ## get the guide grobs
    grobList <- get("guide_grobs", envir = pgEnv)$children

    ## Get the last grob
    last_grob <- grobList[length(grobList)][[1]]

    ## Remove the last grob from the list of grobs
    grobList <- grobList[-length(grobList)]

    ## Re-add grobs to guide_grobs gTree without re-drawing
    invisible(lapply(grobList, add_guideGrobs))

    ## Re-add last grob, this time redrawing
    grid.add("guide_grobs", child = last_grob)
}
