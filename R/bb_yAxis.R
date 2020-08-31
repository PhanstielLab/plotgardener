#' wrapper function for plotting an y axis in its own viewport
#'
#' @param plot plot to add axis to
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param at a numeric vector of x-value locations for the tick marks
#' @param label a logical value indicating whether to draw the labels on the tick marks, or an expression or character vector which specify the labels to use. If not logical, must be the same length as the at argument
#' @param main A logical value indicating whether to draw the axis at the left (TRUE) or at the right (FALSE) of the plot
#' @export
bb_yAxis <- function(plot, params = NULL, at = NULL, label = TRUE, main = TRUE, gp = gpar()){

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(label)) label <- NULL
  if(missing(main)) main <- NULL

  ## Check if plot argument is missing (could be in object)
  if(!hasArg(plot)) plot <- NULL

  ## Compile all parameters into an internal object
  bb_yInternal <- structure(list(plot = plot, at = at, label = label, main = main, gp = gp), class = "bb_yInternal")

  bb_yInternal <- parseParams(bb_params = params, object_params = bb_yInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_yInternal$label)) bb_yInternal$label <- TRUE
  if(is.null(bb_yInternal$main)) bb_yInternal$main <- TRUE

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_yInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  check_bbpage(error = "Cannot add an y-axis without a BentoBox page.")

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  yAxis <- structure(list(grobs = NULL), class = "bb_yaxis")

  # ======================================================================================================================================================================================
  # CREATE GROB WITHOUT DRAWING
  # ======================================================================================================================================================================================

  yGrob <- yaxisGrob(at = bb_yInternal$at, label = bb_yInternal$label, main = bb_yInternal$main, gp = bb_yInternal$gp, vp = bb_yInternal$plot$grobs$vp)

  # ======================================================================================================================================================================================
  # GET CENTER OF INPUT PLOT VIEWPORT BASED ON JUSTIFICATION
  # ======================================================================================================================================================================================

  adjusted_vp <- adjust_vpCoords(viewport = bb_yInternal$plot$grobs$vp)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Make viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_yaxis", length(grep(pattern = "bb_yaxis", x = currentViewports)) + 1)

  ## Define viewport
  if (bb_yInternal$main == TRUE){

    vp <- viewport(width = widthDetails(yGrob), height = bb_yInternal$plot$grobs$vp$height,
                   x = adjusted_vp[[1]] - 0.5 * (bb_yInternal$plot$grobs$vp$width), y = adjusted_vp[[2]] - 0.5 * (bb_yInternal$plot$grobs$vp$height),
                   just = c("right", "bottom"), yscale = bb_yInternal$plot$grobs$vp$yscale, name = vp_name)


  } else if (bb_yInternal$main == FALSE){

    vp <- viewport(width = widthDetails(yGrob), height = bb_yInternal$plot$grobs$vp$height,
                   x = adjusted_vp[[1]] + 0.5 * (bb_yInternal$plot$grobs$vp$width), y = adjusted_vp[[2]] - 0.5 * (bb_yInternal$plot$grobs$vp$height),
                   just = c("left", "bottom"), yscale = bb_yInternal$plot$grobs$vp$yscale, name = vp_name)

  }

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  yGrob <- grid.yaxis(at = bb_yInternal$at, label = bb_yInternal$label, main = bb_yInternal$main, gp = bb_yInternal$gp, vp = vp)
  yaxis_grobs <- gTree(vp = vp, children = gList(yGrob))
  yAxis$grobs <- yaxis_grobs

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(yAxis)

}
