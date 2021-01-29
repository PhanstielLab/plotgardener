#' Draw a line segment within a BentoBox layout
#'
#' @param x0 A numeric vector or unit object indicating the starting x-values of the line segments.
#' @param y0 A numeric vector or unit object indicating the starting y-values of the line segments.
#' @param x1 A numeric vector or unit object indicating the stopping x-values of the line segments.
#' @param y1 A numeric vector or unit object indicating the stopping y-values of the line segments.
#' @param default.units A string indicating the default units to use if \code{x0}, \code{y0}, \code{x1}, or \code{y1} are only given as numeric vectors. Default value is \code{default.units = "inches"}.
#' @param linecolor A character value specifying segment line color. Default value is \code{linecolor = "black"}.
#' @param lwd A numeric specifying segment line width. Default value is \code{lwd = 1}.
#' @param lty A numeric specifying segment line type. Default value is \code{lty = 1}.
#' @param lineend A character value specifying line end style. Default value is \code{lineend = "butt"}. Options are:
#' \itemize{
#' \item{\code{"round"}: Segment ends are rounded.}
#' \item{\code{"butt"}: Segment ends end exactly where ended.}
#' \item{\code{"square"}: Segment ends are squared.}
#' }
#' @param linejoin A character value specifying line join style. Default value is \code{linejoin = "mitre"}. Options are:
#' \itemize{
#' \item{\code{"round"}: }{Line joins are rounded.}
#' \item{\code{"mitre"}: }{Line joins are sharp corners.}
#' \item{\code{"bevel"}: }{Line joins are flattened corners.}
#' }
#' @param arrow A list describing arrow heads to place at either end of the line segments, as produced by the \link[grid]{arrow} function.
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_segments} object containing relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a BentoBox page
#' bb_pageCreate(width = 2, height = 2, default.units = "inches", xgrid = 0, ygrid = 0)
#'
#' ## Plot a line segment
#' bb_plotSegments(x0 = 0.5, y0 = 0.25, x1 = 1.5, y1 = 1.75, default.units = "inches")
#'
#' @seealso \link[grid]{grid.segments}
#'
#' @export
bb_plotSegments <- function(x0, y0, x1, y1, default.units = "inches", linecolor = "black", lwd = 1, lty = 1, lineend = "butt", linejoin = "mitre", arrow = NULL, params = NULL, ...){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(arrow)) arrow <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(lwd)) lwd <- NULL
  if(missing(lty)) lty <- NULL
  if(missing(lineend)) lineend <- NULL
  if(missing(linejoin)) linejoin <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if x0/y0/x1/y1 arguments are missing (could be in object)
  if(!hasArg(x0)) x0 <- NULL
  if(!hasArg(y0)) y0 <- NULL
  if(!hasArg(x1)) x1 <- NULL
  if(!hasArg(y1)) y1 <- NULL

  ## Compile all parameters into an internal object
  bb_segmentsInternal <- structure(list(x0 = x0, y0 = y0, x1 = x1, y1 = y1, arrow = arrow, linecolor = linecolor,
                                      lwd = lwd, lty = lty, lineend = lineend, linejoin = linejoin, default.units = default.units), class = "bb_segmentsInternal")

  bb_segmentsInternal <- parseParams(bb_params = params, object_params = bb_segmentsInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_segmentsInternal$arrow)) bb_segmentsInternal$arrow <- NULL
  if(is.null(bb_segmentsInternal$linecolor)) bb_segmentsInternal$linecolor <- "black"
  if(is.null(bb_segmentsInternal$lwd)) bb_segmentsInternal$lwd <- 1
  if(is.null(bb_segmentsInternal$lty)) bb_segmentsInternal$lty <- 1
  if(is.null(bb_segmentsInternal$lineend)) bb_segmentsInternal$lineend <- "butt"
  if(is.null(bb_segmentsInternal$linejoin)) bb_segmentsInternal$linejoin <- "mitre"
  if(is.null(bb_segmentsInternal$default.units)) bb_segmentsInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_segments <- structure(list(x0 = bb_segmentsInternal$x0, y0 = bb_segmentsInternal$y0, x1 = bb_segmentsInternal$x1, y1 = bb_segmentsInternal$y1,
                                arrow = bb_segmentsInternal$arrow, grobs = NULL,
                                gp = gpar(col = bb_segmentsInternal$linecolor, lwd = bb_segmentsInternal$lwd,lty = bb_segmentsInternal$lty, lineend = bb_segmentsInternal$lineend,
                                          linejoin = bb_segmentsInternal$linejoin, ...)), class = "bb_segmentsInternal")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot plot segment without a BentoBox page.")
  if(is.null(bb_segments$x0)) stop("argument \"x0\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_segments$y0)) stop("argument \"y0\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_segments$x1)) stop("argument \"x1\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_segments$y1)) stop("argument \"y1\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  if (!"unit" %in% class(bb_segments$x0)){

    if (!is.numeric(bb_segments$x0)){

      stop("Starting x-coordinate is neither a unit object or a numeric value. Cannot plot segment.", call. = FALSE)

    }

    if (is.null(bb_segmentsInternal$default.units)){

      stop("Starting x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_segments$x0 <- unit(bb_segments$x0, bb_segmentsInternal$default.units)

  }

  if (!"unit" %in% class(bb_segments$y0)){

    if (!is.numeric(bb_segments$y0)){

      stop("Starting y-coordinate is neither a unit object or a numeric value. Cannot plot segment.", call. = FALSE)

    }

    if (is.null(bb_segmentsInternal$default.units)){

      stop("Starting y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_segments$y0 <- unit(bb_segments$y0, bb_segmentsInternal$default.units)

  }

  if (!"unit" %in% class(bb_segments$x1)){

    if (!is.numeric(bb_segments$x1)){

      stop("Stopping x-coordinate is neither a unit object or a numeric value. Cannot plot segment.", call. = FALSE)

    }

    if (is.null(bb_segmentsInternal$default.units)){

      stop("Stopping x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_segments$x1 <- unit(bb_segments$x1, bb_segmentsInternal$default.units)

  }

  if (!"unit" %in% class(bb_segments$y1)){

    if (!is.numeric(bb_segments$y1)){

      stop("Stopping y-coordinate is neither a unit object or a numeric value. Cannot plot segment.", call. = FALSE)

    }

    if (is.null(bb_segmentsInternal$default.units)){

      stop("Stopping y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_segments$y1 <- unit(bb_segments$y1, bb_segmentsInternal$default.units)

  }
  ## Convert coordinates to page_units
  new_x0 <- convertX(bb_segments$x0, unitTo = page_units, valueOnly = TRUE)
  new_y0 <- convertY(bb_segments$y0, unitTo = page_units, valueOnly = TRUE)
  new_x1 <- convertX(bb_segments$x1, unitTo = page_units, valueOnly = TRUE)
  new_y1 <- convertY(bb_segments$y1, unitTo = page_units, valueOnly = TRUE)

  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================

   segments <- grid.segments(x0 = unit(new_x0, page_units), y0 = unit(page_height - new_y0, page_units), x1 = unit(new_x1, page_units),
                             y1 = unit(page_height - new_y1, page_units), gp = bb_segments$gp)

  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================

  bb_segments$grobs <- segments

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_segments[", segments$name, "]"))
  invisible(bb_segments)
}
