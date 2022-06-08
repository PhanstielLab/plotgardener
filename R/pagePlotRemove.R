#' Remove plotgardener plots and annotations
#'
#' @param plot Plot object to be removed, defined by the output
#' of a plotgardener plotting function.
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
pagePlotRemove <- function(plot) {
    
    if (is(plot, "yaxis") | is(plot, "xaxis")){
        grid.remove(gPath(plot$grobs$childrenOrder))
    } else {
        grid.remove(gPath(plot$grobs$name))
        
        pg_vpTree <- get("pg_vpTree", envir = pgEnv)
        pg_vpTree[grep(plot$grobs$name, pg_vpTree)] <- NULL
        assign("pg_vpTree", pg_vpTree, envir = pgEnv)
        
        ## Need to remove outer viewport of bb_hicTriangle/bb_hicRectangle plot
        ## and domain Clips
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
        } else if (is(plot, "signal")){
            if (grepl("_v", plot$grobs$vp$name)){
                vp_name <- plot$grobs$vp$name
                vp_name <- gsub("_v", "_vClip", vp_name)
                suppressMessages(seekViewport(vp_name))
                suppressMessages(popViewport())
            }
        } 
    }

}
