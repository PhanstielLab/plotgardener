#' Remove BentoBox plots and annotations
#'
#' @param plot BentoBox plot object to be removed, defined by the output
#' of a BentoBox plotting function.
#'
#' @examples
#' ## Load Hi-C data
#' library(BentoBoxData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 5.5, height = 4, default.units = "inches")
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- bb_plotHicSquare(
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
#' bb_pagePlotRemove(plot = hicPlot)
#' @return None.
#'
#' @export
bb_pagePlotRemove <- function(plot) {
    grid.remove(gPath(plot$grobs$name))

    bb_vpTree <- get("bb_vpTree", envir = bbEnv)
    bb_vpTree[grep(plot$grobs$name, bb_vpTree)] <- NULL
    assign("bb_vpTree", bb_vpTree, envir = bbEnv)

    ## Need to remove outer viewport of bb_hicTriangle/bb_hicRectangle plot
    ## and domain Clips
    if (is(plot, "bb_hicTriangle") |
        is(plot, "bb_hicRectangle")) {
        vp_name <- plot$grobs$vp$name
        vp_name <- gsub("inside", "outside", vp_name)
        suppressMessages(seekViewport(vp_name))
        suppressMessages(popViewport())
    } else if (is(plot, "bb_domain")) {
        if (plot$hicClass != "bb_hicSquare") {
            vp_name <- plot$grobs$vp$name
            vp_name <- gsub("inside", "outside", vp_name)
            suppressMessages(seekViewport(vp_name))
            suppressMessages(popViewport())
        }
    } else if (is(plot, "bb_signal")){
        if (grepl("_v", plot$grobs$vp$name)){
            vp_name <- plot$grobs$vp$name
            vp_name <- gsub("_v", "_vClip", vp_name)
            suppressMessages(seekViewport(vp_name))
            suppressMessages(popViewport())
        }
    }
}
