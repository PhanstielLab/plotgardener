#' Removes guides from a BentoBox page
#' @export
#'
bb_pageHideGuides <- function(){


  if (length(get("guide_grobs", envir = bbEnv)$children) == 0) {

    stop("No BentoBox page guides exist.")

  }


  ## Get the list of grobs from the guide grobs gTree and convert to gPaths
  grob_list <- lapply(get("guide_grobs", envir = bbEnv)$children, convert_gpath)

  ## Get the last grob
  last_grob <- grob_list[length(grob_list)]

  ## Remove the last grob from the list of grobs
  grob_list <- grob_list[- length(grob_list)]

  ## Remove all grobs except last one, not redrawing each time (to save time)
  invisible(lapply(grob_list, grid.remove, redraw = FALSE))

  ## Remove last grob with redrawing, now removing all the grobs
  invisible(lapply(last_grob, grid.remove))

}
