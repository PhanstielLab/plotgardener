#' draws a horizontal guideline at a specified y-coordinate
#'
#' @param y y-coordinate to draw guide
#' @param units units of y-coordinate
#' @param col color of guideline

#' @export
bb_hGuide <- function(y, units = "inches", col = "grey55"){

  ## Get the names of the current viewports
  current_viewports <- lapply(current.vpTree()$children, viewport_name)

  if (!"bb_page" %in% current_viewports){

    stop("Must make a BentoBox page with bb_makePage() before adding additional guidelines.")

  }

  y <- convertY(y, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)

  guide <- grid.segments(x0 = unit(0, units = "npc"), x1 = unit(1, units = "npc"), y0 = get("page_height", envir = bbEnv) - y, y1 = get("page_height", envir = bbEnv) - y,
                           default.units = get("page_units", envir = bbEnv), gp = gpar(col = col))
  assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = guide), envir = bbEnv)

}
