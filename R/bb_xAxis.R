#' wrapper function for plotting an x axis in its own viewport
#' @param plot plot to add axis to
#' @param at a numeric vector of x-value locations for the tick marks
#' @param label a logical value indicating whether to draw the labels on the tick marks, or an expression or character vector which specify the labels to use. If not logical, must be the same length as the at argument
#' @param main A logical value indicating whether to draw the axis at the bottom (TRUE) or at the top (FALSE) of the plot
#' @export
bb_xAxis <- function(plot, at = NULL, label = TRUE, main = TRUE, gp = gpar()){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  xAxis <- structure(list(grobs = NULL, viewport = NULL), class = "bb_xaxis")

  # ======================================================================================================================================================================================
  # CREATE GROB WITHOUT DRAWING
  # ======================================================================================================================================================================================

  xGrob <- xaxisGrob(at = at, label = label, main = main, gp = gp, vp = plot$viewport)

  # ======================================================================================================================================================================================
  # GET CENTER OF INPUT PLOT VIEWPORT BASED ON JUSTIFICATION
  # ======================================================================================================================================================================================

  adjusted_vp <- adjust_vpCoords(viewport = plot$viewport)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================
  ## Make viewport name
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_xaxis", length(grep(pattern = "bb_xaxis", x = current_viewports)) + 1)

  ## Define viewport
  if (main == TRUE){

    vp <- viewport(width = plot$viewport$width, height = heightDetails(xGrob),
                   x = adjusted_vp[[1]] - 0.5 * (plot$viewport$width), y = adjusted_vp[[2]] - 0.5 * (plot$viewport$height),
                   just = c("left", "top"), xscale = plot$viewport$xscale, name = vp_name)


  } else if (main == FALSE){

    vp <- viewport(width = plot$viewport$width, height = heightDetails(xGrob),
                   x = adjusted_vp[[1]] - 0.5 * (plot$viewport$width), y = adjusted_vp[[2]] + 0.5 * (plot$viewport$height),
                   just = c("left", "bottom"), xscale = plot$viewport$xscale, name = vp_name)

  }

  xAxis$viewport <- vp
  pushViewport(vp)
  upViewport()

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  xGrob <- grid.xaxis(at = at, label = label, main = main, gp = gp, vp = vp_name)
  xaxis_grobs <- gTree(name = "xaxis_grobs", children = gList(xGrob))
  xAxis$grobs <- xaxis_grobs$children

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(xAxis)

}
