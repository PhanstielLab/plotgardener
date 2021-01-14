#' adds a plot created with ggplot2 to a BentoBox layout
#'
#' @param plot ggplot object'
#' @param x A numeric or unit object specifying x-location.
#' @param y A numeric or unit object specifying y-location.
#' @param width A numeric or unit object specifying width.
#' @param height A numeric or unit object specifying height.
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param just A string or numeric vector specifying the justification of the plot relative to its (x, y) location
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#'
#' @export
bb_plotGG <- function(plot, x, y, width, height, params = NULL, just = c("left", "top"), default.units = "inches"){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if plot/x/y/width/height arguments are missing (could be in object)
  if(!hasArg(plot)) plot <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(width)) width <- NULL
  if(!hasArg(height)) height <- NULL

  ## Compile all parameters into an internal object
  bb_ggInternal <- structure(list(plot = plot, x = x, y = y, width = width, height = height,
                                     just = just, default.units = default.units), class = "bb_ggInternal")

  bb_ggInternal <- parseParams(bb_params = params, object_params = bb_ggInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_ggInternal$just)) bb_ggInternal$just <- c("left", "top")
  if(is.null(bb_ggInternal$default.units)) bb_ggInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # ERRORS
  # ======================================================================================================================================================================================

  if (is.null(bb_ggInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  if (is.null(bb_ggInternal$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if (is.null(bb_ggInternal$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if (is.null(bb_ggInternal$width)) stop("argument \"width\" is missing, with no default.", call. = FALSE)
  if (is.null(bb_ggInternal$height)) stop("argument \"height\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # INITIALIZE PLOT OBJECT
  # ======================================================================================================================================================================================

  gg_plot <- structure(list(width = bb_ggInternal$width, height = bb_ggInternal$height, x = bb_ggInternal$x, y = bb_ggInternal$y,
                            justification = bb_ggInternal$just), class = "bb_gg")

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  gg_plot <- defaultUnits(object = gg_plot, default.units = bb_ggInternal$default.units)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = gg_plot)

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_gg", length(grep(pattern = "bb_gg", x = currentViewports)) + 1)

  ## Make viewport for gene track
  vp <- viewport(height = page_coords$height, width = page_coords$width,
                 x = page_coords$x, y = page_coords$y,
                 just = bb_ggInternal$just, name = vp_name)

  # ======================================================================================================================================================================================
  # PRINT GGPLOT
  # ======================================================================================================================================================================================

  print(bb_ggInternal$plot, vp = vp)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(gg_plot)

}
