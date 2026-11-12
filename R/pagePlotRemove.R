#' Remove plotgardener plots and annotations
#'
#' Records a plot for deferred removal instead of calling `grid.remove()` on
#' the active device. The removal is applied by \code{pageExport()} or
#' \code{pageRedraw()}, avoiding extra pages in PDF output.
#'
#' @param plot Plot object to be removed, defined by the output
#' of a plotgardener plotting function.
#' @param redraw Logical value indicating whether to redraw immediately on a
#' screen device. File devices defer the change to \code{pageExport()}.
#' Default value is \code{redraw = TRUE}.
#'
#' @examples
#' ## Load Hi-C data
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create page
#' pageCreate(width = 5.5, height = 4, default.units = "inches")
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- plotHicSquare(
#'     data = IMR90_HiC_10kb, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     x = 0.5, y = 0.5, width = 2.5, height = 2.5,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Remove square Hi-C plot from page
#' pagePlotRemove(plot = hicPlot)
#' @return None.
#'
#' @export
pagePlotRemove <- function(plot, redraw = TRUE) {

    removed <- get("removed_plots", envir = pgEnv)

    if (is(plot, "yaxis") | is(plot, "xaxis")) {
        removed <- c(removed, plot$grobs$childrenOrder)
    } else {
        removed <- c(removed, plot$grobs$name)

        pg_vpTree <- get("pg_vpTree", envir = pgEnv)
        pg_vpTree[grep(plot$grobs$name, pg_vpTree)] <- NULL
        assign("pg_vpTree", pg_vpTree, envir = pgEnv)

        ## Keep viewport bookkeeping compatible with the previous behavior,
        ## but do not remove the plotted grob from the device.
        if (is(plot, "hicTriangle") |
            is(plot, "hicRectangle")) {
            vp_name <- plot$grobs$vp$name
            vp_name <- gsub("inside", "outside", vp_name)
            suppressMessages(seekViewport(vp_name))
            suppressMessages(popViewport())
        } else if (is(plot, "domain")) {
            if (plot$hicClass != "hicSquare") {
                vp_name <- plot$grobs$vp$name
                vp_name <- gsub("inside", "outside", vp_name)
                suppressMessages(seekViewport(vp_name))
                suppressMessages(popViewport())
            }
        } else if (is(plot, "signal")) {
            if (grepl("_v", plot$grobs$vp$name)) {
                vp_name <- plot$grobs$vp$name
                vp_name <- gsub("_v", "_vClip", vp_name)
                suppressMessages(seekViewport(vp_name))
                suppressMessages(popViewport())
            }
        }
    }

    assign("removed_plots", unique(removed), envir = pgEnv)
    if (isTRUE(redraw) && !.isFileDevice()) {
        pageRedraw()
    }
    invisible(NULL)
}
