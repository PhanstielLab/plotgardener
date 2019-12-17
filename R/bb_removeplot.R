#' removes BentoBox plots
#'
#' @param plot BentoBox plot to remove; defined by the output of a BentoBox plotting function
#'
#' @export
#'
bb_removePlot <- function(plot){

  ## error function: need a plot that's a valid BentoBox plot, make sure the plot is plotted

  ## get the grobs of the plot and convert to gpath
  grob_list <- lapply(plot$grobs$children, convert_gpath)

  ## get the last grob
  last_grob <- grob_list[length(grob_list)]

  ## remove the last grob from the list of grobs
  grob_list <- grob_list[- length(grob_list)]

  ## remove all grobs except last one, not redrawing each time (to save time)
  invisible(lapply(grob_list, grid.remove, redraw = FALSE))

  ## remove last grob with redrawing, now removing all the grobs
  invisible(lapply(last_grob, grid.remove))

  ## remove the viewport
  seekViewport(name = plot$grobs$vp$name)
  popViewport()


}
