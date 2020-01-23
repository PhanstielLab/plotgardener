#' draws a vertical guideline at a specified x-coordinate
#'
#' @param x unit object specifying x-coordinate of guide
#' @param col color of guideline

#' @export
bb_vGuide <- function(x, col = "grey55"){

  ## Get the names of the current viewports
  # current_viewports <- lapply(current.vpTree()$children, viewport_name)
  #
  # if (!"bb_page" %in% current_viewports){
  #
  #   stop("Must make a BentoBox page with bb_makePage() before adding additional guidelines.")
  #
  # }

  guide <- grid.segments(x0 = x, x1 = x, y0 = unit(0, units = "npc"), y1 = unit(1, units = "npc"),
                           gp = gpar(col = col))
  assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = guide), envir = bbEnv)

}
