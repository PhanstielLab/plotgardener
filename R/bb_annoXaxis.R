#' wrapper function for plotting an x axis in its own viewport
#'
#' @param plot plot to add axis to
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param at a numeric vector of x-value locations for the tick marks
#' @param label a logical value indicating whether to draw the labels on the tick marks, or an expression or character vector which specify the labels to use. If not logical, must be the same length as the at argument
#' @param main A logical value indicating whether to draw the axis at the bottom (TRUE) or at the top (FALSE) of the plot
#' @export
bb_annoXaxis <- function(plot, params = NULL, at = NULL, label = TRUE, main = TRUE, gp = gpar()){

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(label)) label <- NULL
  if(missing(main)) main <- NULL

  ## Check if plot argument is missing (could be in object)
  if(!hasArg(plot)) plot <- NULL

  ## Compile all parameters into an internal object
  bb_xInternal <- structure(list(plot = plot, at = at, label = label, main = main, gp = gp), class = "bb_xInternal")

  bb_xInternal <- parseParams(bb_params = params, object_params = bb_xInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_xInternal$label)) bb_xInternal$label <- TRUE
  if(is.null(bb_xInternal$main)) bb_xInternal$main <- TRUE

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_xInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  check_bbpage(error = "Cannot add an x-axis without a BentoBox page.")

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  xAxis <- structure(list(grobs = NULL), class = "bb_xaxis")

  # ======================================================================================================================================================================================
  # CREATE GROB WITHOUT DRAWING
  # ======================================================================================================================================================================================

  xGrob <- xaxisGrob(at = bb_xInternal$at, label = bb_xInternal$label, main = bb_xInternal$main, gp = bb_xInternal$gp, vp = bb_xInternal$plot$grobs$vp)

  # ======================================================================================================================================================================================
  # GET CENTER OF INPUT PLOT VIEWPORT BASED ON JUSTIFICATION
  # ======================================================================================================================================================================================

  adjusted_vp <- adjust_vpCoords(viewport = bb_xInternal$plot$grobs$vp)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Make viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_xaxis", length(grep(pattern = "bb_xaxis", x = currentViewports)) + 1)

  ## Define viewport
  if (bb_xInternal$main == TRUE){

    vp <- viewport(width = bb_xInternal$plot$grobs$vp$width, height = heightDetails(xGrob),
                   x = adjusted_vp[[1]] - 0.5 * (bb_xInternal$plot$grobs$vp$width), y = adjusted_vp[[2]] - 0.5 * (bb_xInternal$plot$grobs$vp$height),
                   just = c("left", "top"), xscale = bb_xInternal$plot$grobs$vp$xscale, name = vp_name)


  } else if (bb_xInternal$main == FALSE){

    vp <- viewport(width = bb_xInternal$plot$grobs$vp$width, height = heightDetails(xGrob),
                   x = adjusted_vp[[1]] - 0.5 * (bb_xInternal$plot$grobs$vp$width), y = adjusted_vp[[2]] + 0.5 * (bb_xInternal$plot$grobs$vp$height),
                   just = c("left", "bottom"), xscale = bb_xInternal$plot$grobs$vp$xscale, name = vp_name)

  }

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  xGrob <- grid.xaxis(at = bb_xInternal$at, label = bb_xInternal$label, main = bb_xInternal$main, gp = bb_xInternal$gp, vp = vp)
  xaxis_grobs <- gTree(vp = vp, children = gList(xGrob))
  xAxis$grobs <- xaxis_grobs

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(xAxis)

}
