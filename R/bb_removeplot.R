#' Removes BentoBox plots and annotations
#'
#' @param plot BentoBox plot to remove, defined by the output of a BentoBox plotting function.
#'
#' @export
bb_removePlot <- function(plot){

  ## error function: need a plot that's a valid BentoBox plot, make sure the plot is plotted
  grid.remove(gPath(plot$grobs$name))


  ## Need to remove outer viewport of bb_trianglehic plot
  if (class(plot) == "bb_trianglehic"){
    vp_name <- plot$grobs$vp$name
    vp_name <- gsub("inside", "outside", vp_name)
    seekViewport(vp_name)
    popViewport()

  }

}

