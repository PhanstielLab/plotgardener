#' wrapper to draw a textGrob based on bb_makePage coordinates and units
#'
#' @param label character or expression
#' @param x A numeric vector or unit object specifying x-location
#' @param y A numeric vector or unit object specifying y-location
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param just justification of text relative to its (x, y) location
#' @param fontcolor fontcolor
#' @param fontsize the size of text (in points)
#' @param rot the angle to rotate the song
#' @param check.overlap A logical value to indicate whether to check for and omit overlapping text.
#' @param default.units A string indicating the default units to use if x or y are only given as numeric vectors
#'
#' @export
bb_addText <- function(label, x, y, params = NULL, just = "center", fontcolor = "black", fontsize = 12, rot = 0, check.overlap = FALSE, default.units = "inches", ...){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(fontcolor)) fontcolor <- NULL
  if(missing(fontsize)) fontsize <- NULL
  if(missing(rot)) rot <- NULL
  if(missing(check.overlap)) check.overlap <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if label/x/y arguments are missing (could be in object)
  if(!hasArg(label)) label <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL

  ## Compile all parameters into an internal object
  bb_textInternal <- structure(list(label = label, x = x, y = y, just = just, fontcolor = fontcolor,
                                    fontsize = fontsize, rot = rot, check.overlap = check.overlap, default.units = default.units), class = "bb_textInternal")

  bb_textInternal <- parseParams(bb_params = params, object_params = bb_textInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_textInternal$just)) bb_textInternal$just <- "center"
  if(is.null(bb_textInternal$fontcolor)) bb_textInternal$fontcolor <- "black"
  if(is.null(bb_textInternal$fontsize)) bb_textInternal$fontsize <- 12
  if(is.null(bb_textInternal$rot)) bb_textInternal$rot <- 0
  if(is.null(bb_textInternal$check.overlap)) bb_textInternal$check.overlap <- FALSE
  if(is.null(bb_textInternal$default.units)) bb_textInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_text <- structure(list(label = bb_textInternal$label, x = bb_textInternal$x, y = bb_textInternal$y, just = bb_textInternal$just, rot = bb_textInternal$rot, grobs = NULL,
                            gp = gpar(col = bb_textInternal$fontcolor, fontsize = bb_textInternal$fontsize, ...)), class = "bb_text")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot add text without a BentoBox page.")
  if(is.null(bb_text$label)) stop("argument \"label\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_text$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_text$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  if (class(bb_text$x) != "unit"){

    if (!is.numeric(bb_text$x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_textInternal$default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_text$x <- unit(bb_text$x, bb_textInternal$default.units)

  }

  if (class(bb_text$y) != "unit"){

    if (!is.numeric(bb_text$y)){

      stop("y-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_textInternal$default.units)){

      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_text$y <- unit(bb_text$y, bb_textInternal$default.units)

  }

  ## Convert coordinates to page_units
  new_x <- convertX(bb_text$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(bb_text$y, unitTo = page_units, valueOnly = TRUE)

  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================

  text <- grid.text(label = bb_text$label, x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), just = bb_text$just,
            gp = bb_text$gp, rot = bb_text$rot, check.overlap = bb_textInternal$check.overlap)

  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================

  bb_text$grobs <- text

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_text)
}
