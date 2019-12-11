#' wrapper function for plotting an y axis in its own viewport
#' @param plot plot to add axis to
#' @param at a numeric vector of x-value locations for the tick marks
#' @param label a logical value indicating whether to draw the labels on the tick marks, or an expression or character vector which specify the labels to use. If not logical, must be the same length as the at argument
#' @param main A logical value indicating whether to draw the axis at the left (TRUE) or at the right (FALSE) of the plot
#' @export
bb_yAxis <- function(plot, at = NULL, label = TRUE, main = TRUE, gp = gpar()){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  yAxis <- structure(list(grobs = NULL, viewport = NULL), class = "bb_yaxis")

  # ======================================================================================================================================================================================
  # CREATE GROB WITHOUT DRAWING
  # ======================================================================================================================================================================================

  yGrob <- yaxisGrob(at = at, label = label, main = main, gp = gp, vp = plot$viewport)

  # ======================================================================================================================================================================================
  # GET CENTER OF INPUT PLOT VIEWPORT BASED ON JUSTIFICATION
  # ======================================================================================================================================================================================

  adjusted_vp <- adjust_vpCoords(viewport = plot$viewport)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Make viewport name
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_yaxis", length(grep(pattern = "bb_yaxis", x = current_viewports)) + 1)

  ## Define viewport
  if (main == TRUE){

    vp <- viewport(width = widthDetails(yGrob), height = plot$viewport$height,
                   x = adjusted_vp[[1]] - 0.5 * (plot$viewport$width), y = adjusted_vp[[2]] - 0.5 * (plot$viewport$height),
                   just = c("right", "bottom"), yscale = plot$viewport$yscale, name = vp_name)


  } else if (main == FALSE){

    vp <- viewport(width = widthDetails(yGrob), height = plot$viewport$height,
                   x = adjusted_vp[[1]] + 0.5 * (plot$viewport$width), y = adjusted_vp[[2]] - 0.5 * (plot$viewport$height),
                   just = c("left", "bottom"), yscale = plot$viewport$yscale, name = vp_name)

  }

  yAxis$viewport <- vp
  pushViewport(vp)
  upViewport()

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  yGrob <- grid.yaxis(at = at, label = label, main = main, gp = gp, vp = vp_name)
  yaxis_grobs <- gTree(name = "yaxis_grobs", children = gList(yGrob))
  yAxis$grobs <- yaxis_grobs$children

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(yAxis)

}
