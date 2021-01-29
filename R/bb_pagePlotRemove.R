#' Remove BentoBox plots and annotations
#'
#' @param plot BentoBox plot object to be removed, defined by the output of a BentoBox plotting function.
#'
#' @examples
#' ## Load Hi-C data
#' data("bb_hicData")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 4, height = 3.5, default.units = "inches", xgrid = 0, ygrid = 0)
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- bb_plotHicSquare(hicData = bb_hicData, resolution = 10000, zrange = c(0, 70),
#'                             chrom = "chr21", chromstart = 28000000, chromend = 30300000,
#'                             x = 0.5, y = 0.5, width = 2.5, height = 2.5,
#'                             just = c("left", "top"), default.units = "inches")
#'
#' ## Remove square Hi-C plot from page
#' bb_pagePlotRemove(plot = hicPlot)
#'
#' @export
bb_pagePlotRemove <- function(plot){

  ## error function: need a plot that's a valid BentoBox plot, make sure the plot is plotted
  grid.remove(gPath(plot$grobs$name))


  ## Need to remove outer viewport of bb_hicTriangle plot
  if (class(plot) == "bb_hicTriangle"){
    vp_name <- plot$grobs$vp$name
    vp_name <- gsub("inside", "outside", vp_name)
    seekViewport(vp_name)
    popViewport()

  }

}

