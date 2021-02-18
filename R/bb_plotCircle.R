#' Plot a circle within a BentoBox layout
#'
#' @param x A numeric vector or unit object specifying circle x-locations relative to center.
#' @param y A numeric vector or unit object specifying circle y-locations relative to center.
#' @param r A numeric vector or unit object specifying radii.
#' @param default.units A string indicating the default units to use if \code{r}, \code{x}, or \code{y} are only given as numerics or numeric vectors. Default value is \code{default.units = "inches"}.
#' @param linecolor A character value specifying circle line color. Default value is \code{linecolor = "black"}.
#' @param lwd A numeric specifying circle line width. Default value is \code{lwd = 1}.
#' @param lty A numeric specifying circle line type. Default value is \code{lty = 1}.
#' @param fill A character value specifying circle fill color. Default value is \code{fill = NA}.
#' @param alpha Numeric value specifying color transparency. Default value is \code{alpha = 1}.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_circle} object containing relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a BentoBox page
#' bb_pageCreate(width = 2, height = 2, default.units = "inches")
#'
#' ## Plot two circles, one at a time
#' bb_plotCircle(x = 0.6, y = 0.5, r = 0.1, fill = "black", default.units = "inches")
#' bb_plotCircle(x = 1.4, y = 0.5, r = 0.1, fill = "black", default.units = "inches")
#'
#' ## Plot a vector of circles
#' xVals <- 1 + (0.5*cos(seq(0, pi, pi/8)))
#' yVals <- 1 + (0.5*sin(seq(0, pi, pi/8)))
#' bb_plotCircle(x = xVals, y = yVals, r = 0.05, default.units = "inches")
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @seealso \link[grid]{grid.circle}
#'
#' @export
bb_plotCircle <- function(x, y, r, default.units = "inches", linecolor = "black", lwd = 1, lty = 1, fill = NA, alpha = 1, params = NULL, ...){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(fill)) fill <- NULL
  if(missing(lwd)) lwd <- NULL
  if(missing(lty)) lty <- NULL
  if(missing(alpha)) alpha <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if label/x/y arguments are missing (could be in object)
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(r)) r <- NULL

  ## Compile all parameters into an internal object
  bb_circleInternal <- structure(list(x = x, y = y, r = r, linecolor = linecolor, fill = fill,
                                      lwd = lwd, lty = lty, alpha = alpha, default.units = default.units), class = "bb_circleInternal")

  bb_circleInternal <- parseParams(bb_params = params, object_params = bb_circleInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_circleInternal$linecolor)) bb_circleInternal$linecolor <- "black"
  if(is.null(bb_circleInternal$fill)) bb_circleInternal$fill <- NA
  if(is.null(bb_circleInternal$lwd)) bb_circleInternal$lwd <- 1
  if(is.null(bb_circleInternal$lty)) bb_circleInternal$lty <- 1
  if(is.null(bb_circleInternal$alpha)) bb_circleInternal$alpha <- 1
  if(is.null(bb_circleInternal$default.units)) bb_circleInternal$default.units <- "inches"

  ## Set gp
  bb_circleInternal$gp <- gpar(col = bb_circleInternal$linecolor, fill = bb_circleInternal$fill, lwd = bb_circleInternal$lwd,
                               lty = bb_circleInternal$lty, alpha = bb_circleInternal$alpha)
  bb_circleInternal$gp <- setGP(gpList = bb_circleInternal$gp, params = bb_circleInternal, ...)


  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_circle <- structure(list(x = bb_circleInternal$x, y = bb_circleInternal$y, r = bb_circleInternal$r, grobs = NULL,
                              gp = bb_circleInternal$gp), class = "bb_circle")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot plot circle without a BentoBox page.")
  if(is.null(bb_circle$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_circle$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_circle$r)) stop("argument \"r\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  if (!"unit" %in% class(bb_circle$x)){

    if (!is.numeric(bb_circle$x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot plot circle.", call. = FALSE)

    }

    if (is.null(bb_circleInternal$default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_circle$x <- unit(bb_circle$x, bb_circleInternal$default.units)

  }

  if (!"unit" %in% class(bb_circle$y)){

    if (!is.numeric(bb_circle$y)){

      stop("y-coordinate is neither a unit object or a numeric value. Cannot plot circle.", call. = FALSE)

    }

    if (is.null(bb_circleInternal$default.units)){

      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_circle$y <- unit(bb_circle$y, bb_circleInternal$default.units)

  }

  if (!"unit" %in% class(bb_circle$r)){

    if (!is.numeric(bb_circle$r)){

      stop("Radius is neither a unit object or a numeric value. Cannot plot circle.", call. = FALSE)

    }

    if (is.null(bb_circleInternal$default.units)){

      stop("Radius detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_circle$r <- unit(bb_circle$r, bb_circleInternal$default.units)

  }

  ## Convert coordinates to page_units
  new_x <- convertX(bb_circle$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(bb_circle$y, unitTo = page_units, valueOnly = TRUE)
  new_r <- convertUnit(bb_circle$r, unitTo = page_units, valueOnly = TRUE)

  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================

  circle <- grid.circle(x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), r = unit(new_r, page_units),
                        gp = bb_circle$gp)

  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================

  bb_circle$grobs <- circle

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_circle[", circle$name, "]"))
  invisible(bb_circle)
}
