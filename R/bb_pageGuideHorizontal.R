#' Draw a horizontal guideline at a specified y-coordinate on a BentoBox page
#'
#' @usage bb_pageGuideHorizontal(y, default.units = "inches")
#'
#' @param y A numeric or unit object specifying y-coordinate of guide.
#' @param default.units A string indicating the default units to use if \code{y} is only given as a numeric. Default value is \code{default.units = "inches"}.
#' @param linecolor Character value indicating color of guideline. Default value is \code{linecolor = "grey55"}.
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @examples
#' ## Create a BentoBox page
#' bb_pageCreate(width = 4, height = 4, default.units = "inches")
#'
#' ## Add red horizontal guideline at y = 2.5 inches
#' bb_pageGuideHorizontal(y = 2.5, linecolor = "red")
#'
#' @export
bb_pageGuideHorizontal <- function(y, default.units = "inches", linecolor = "grey55", params = NULL, ...){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if y argument is missing (could be in object)
  if(!hasArg(y)) y <- NULL

  ## Compile all parameters into an internal object
  bb_hguide <- structure(list(y = y, linecolor = linecolor, default.units = default.units), class = "bb_hguide")

  bb_hguide <- parseParams(bb_params = params, object_params = bb_hguide)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_hguide$linecolor)) bb_hguide$linecolor <- "grey55"
  if(is.null(bb_hguide$default.units)) bb_hguide$default.units <- "inches"

  # ======================================================================================================================================================================================
  # ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_hguide$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # DEFAULT UNITS
  # ======================================================================================================================================================================================

  if (class(bb_hguide$y) != "unit"){

    if (!is.numeric(bb_hguide$y)){

      stop("y-coordinate is neither a unit object or a numeric value. Cannot place Hguide.", call. = FALSE)

    }

    if (is.null(bb_hguide$default.units)){

      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    y <- unit(bb_hguide$y, bb_hguide$default.units)
  }

  # ======================================================================================================================================================================================
  # MAKE GROB AND ASSIGN TO GTREE
  # ======================================================================================================================================================================================

  y <- convertY(y, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)

  guide <- grid.segments(x0 = unit(0, units = "npc"), x1 = unit(1, units = "npc"), y0 = get("page_height", envir = bbEnv) - y, y1 = get("page_height", envir = bbEnv) - y,
                           default.units = get("page_units", envir = bbEnv), gp = gpar(col = bb_hguide$linecolor, ...))
  assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = guide), envir = bbEnv)

}
