
#' @export
#'
bb_removeplot <- function(plot){

  ## get the grobs of the plot
  grob_list <- plot$grobs

  ## get the last grob
  last_grob <- grob_list[length(grob_list)]

  ## remove the last grob from the list of grobs
  grob_list <- grob_list[- length(grob_list)]

  ## remove all grobs except last one, not redrawing each time (to save time)
  invisible(lapply(grob_list, grid.remove, redraw = FALSE))

  ## remove last grob with redrawing, now removing all the grobs
  invisible(lapply(last_grob, grid.remove))

}
