#' Removes BentoBox plots and annotations
#'
#' @param plot BentoBox plot to remove, defined by the output of a BentoBox plotting function.
#'
#' @export
bb_removePlot <- function(plot){

  ## error function: need a plot that's a valid BentoBox plot, make sure the plot is plotted
  grid.remove(gPath(plot$grobs$name))

}

