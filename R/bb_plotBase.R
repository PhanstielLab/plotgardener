#' adds a plot created with baseR to a BentoBox layout
#'
#' @param plot plot formula of base R plotting functions; ex: p <- ~plot(1:10) + abline(v = 2)
#' @param x A numeric or unit object specifying x-location.
#' @param y A numeric or unit object specifying y-location.
#' @param width A numeric or unit object specifying width.
#' @param height A numeric or unit object specifying height.
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param bg background color
#' @param just A string or numeric vector specifying the justification of the plot relative to its (x, y) location
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#'
#' @export
bb_plotBase <- function(plot, x, y, width, height, params = NULL, bg = NA, just = c("left", "top"), default.units = "inches"){

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(bg)) bg <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if plot/x/y/width/height arguments are missing (could be in object)
  if(!hasArg(plot)) plot <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(width)) width <- NULL
  if(!hasArg(height)) height <- NULL

  ## Compile all parameters into an internal object
  bb_baseInternal <- structure(list(plot = plot, x = x, y = y, width = width, height = height, bg = bg,
                                     just = just, default.units = default.units), class = "bb_baseInternal")

  bb_baseInternal <- parseParams(bb_params = params, object_params = bb_baseInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_baseInternal$bg)) bb_baseInternal$bg <- NA
  if(is.null(bb_baseInternal$just)) bb_baseInternal$just <- c("left", "top")
  if(is.null(bb_baseInternal$default.units)) bb_baseInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # INITIALIZE PLOT OBJECT
  # ======================================================================================================================================================================================

  bb_base <- structure(list(width = bb_baseInternal$width, height = bb_baseInternal$height, x = bb_baseInternal$x, y = bb_baseInternal$y,
                              justification = bb_baseInternal$just, grobs = NULL), class = "bb_base")

  # ======================================================================================================================================================================================
  # CALL ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_baseInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_baseInternal$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_baseInternal$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_baseInternal$width)) stop("argument \"width\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_baseInternal$height)) stop("argument \"height\" is missing, with no default.", call. = FALSE)

  check_bbpage(error = "Must have a BentoBox page before adding a base R plot.")

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  bb_base <- defaultUnits(object = bb_base, default.units = bb_baseInternal$default.units)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = bb_base)

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_base", length(grep(pattern = "bb_base", x = currentViewports)) + 1)

  ## Make viewport
  vp <- viewport(height = page_coords$height, width = page_coords$width,
                      x = page_coords$x, y = page_coords$y,
                      just = bb_baseInternal$just, name = vp_name)

  # ======================================================================================================================================================================================
  # CONVERT PLOT TO A GROB
  # ======================================================================================================================================================================================

  gtree <- ggplotify::base2grob(bb_baseInternal$plot)

  # ======================================================================================================================================================================================
  # ASSIGN VIEWPORT TO GTREE
  # ======================================================================================================================================================================================

  gtree$vp <- vp

  # ======================================================================================================================================================================================
  # BACKGROUND AND BORDER
  # ======================================================================================================================================================================================

  plotChildren <- gtree$childrenOrder
  backgroundGrob <- rectGrob(gp = gpar(fill = bb_baseInternal$bg, col = NA), name = "background")
  gtree <- addGrob(gTree = gtree, child = backgroundGrob)
  gtree$childrenOrder <- c("background", plotChildren)

  # ======================================================================================================================================================================================
  # DRAW GTREE
  # ======================================================================================================================================================================================

  grid.draw(gtree)

  # ======================================================================================================================================================================================
  # ADD GTREE TO OBJECT
  # ======================================================================================================================================================================================

  bb_base$grobs <- gtree

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_base)

}
