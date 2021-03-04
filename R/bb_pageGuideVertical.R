#' Draw a vertical guideline at a specified x-coordinate on a BentoBox page
#'
#' @param x A numeric or unit object specifying x-coordinate of guide.
#' @param default.units A string indicating the default units to use if \code{x} is only given as a numeric. Default value is \code{default.units = "inches"}.
#' @param linecolor Character value indicating color of guideline. Default value is \code{linecolor = "grey55"}.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @examples
#' ## Create a BentoBox page
#' bb_pageCreate(width = 6, height = 5, default.units = "inches")
#'
#' ## Add blue vertical guideline at x = 1.7 inches
#' bb_pageGuideVertical(x = 1.7, linecolor = "blue")
#'
#' @export
bb_pageGuideVertical <- function(x, default.units = "inches", linecolor = "grey55", params = NULL, ...){

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if x argument is missing (could be in object)
  if(!hasArg(x)) x <- NULL

  ## Compile all parameters into an internal object
  bb_vguide <- structure(list(x = x, linecolor = linecolor, default.units = default.units), class = "bb_vguide")

  bb_vguide <- parseParams(bb_params = params, object_params = bb_vguide)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_vguide$linecolor)) bb_vguide$linecolor <- "grey55"
  if(is.null(bb_vguide$default.units)) bb_vguide$default.units <- "inches"

  ## Set gp
  bb_vguide$gp <- gpar(col = bb_vguide$linecolor)
  bb_vguide$gp <- setGP(gpList = bb_vguide$gp, params = bb_vguide, ...)

  # ======================================================================================================================================================================================
  # ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_vguide$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)


  if (class(bb_vguide$x) != "unit"){

    if (!is.numeric(bb_vguide$x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_vguide$default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    x <- unit(bb_vguide$x, bb_vguide$default.units)
  }

  guide <- grid.segments(x0 = x, x1 = x, y0 = unit(0, units = "npc"), y1 = unit(1, units = "npc"),
                           gp = bb_vguide$gp)
  assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = guide), envir = bbEnv)

}
