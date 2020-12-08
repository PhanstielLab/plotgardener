#'draws a vertical guideline at a specified y-coordinate on a BentoBox page
#' @param x A numeric or unit object specifying x-coordinate of guide
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param col color of guideline
#' @param default.units A string indicating the default units to use if x is only given as numeric vectors

#' @export
bb_pageVerticalGuide <- function(x, params = NULL, col = "grey55", default.units = "inches", ...){

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(col)) col <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if x argument is missing (could be in object)
  if(!hasArg(x)) x <- NULL

  ## Compile all parameters into an internal object
  bb_vguide <- structure(list(x = x, col = col, default.units = default.units), class = "bb_vguide")

  bb_vguide <- parseParams(bb_params = params, object_params = bb_vguide)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_vguide$col)) bb_vguide$col <- "grey55"
  if(is.null(bb_vguide$default.units)) bb_vguide$default.units <- "inches"

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
                           gp = gpar(col = bb_vguide$col, ...))
  assign("guide_grobs", addGrob(gTree = get("guide_grobs", envir = bbEnv), child = guide), envir = bbEnv)

}
