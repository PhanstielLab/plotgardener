#' Plot a polygon within a BentoBox layout
#'
#' @param x A numeric vector or unit object specifying polygon vertex x-locations.
#' @param y A numeric vector, unit object, or a character vector of values containing a "b" combined with a numeric value specifying polygon vertex y-locations.
#' The character vector will place polygon vertex y-locations relative to the bottom of the most recently plotted BentoBox plot according to the units of the BentoBox page.
#' @param default.units A string indicating the default units to use if \code{x} or \code{y} are only given as numeric vectors. Default value is \code{default.units = "inches"}.
#' @param linecolor A character value specifying polygon line color. Default value is \code{linecolor = "black"}.
#' @param lwd A numeric specifying polygon line width. Default value is \code{lwd = 1}.
#' @param lty A numeric specifying polygon line type. Default value is \code{lty = 1}.
#' @param fill A character value specifying polygon fill color. Default value is \code{fill = NA}.
#' @param alpha Numeric value specifying color transparency. Default value is \code{alpha = 1}.
#' @param id A numeric vector used to separate locations in \code{x} and \code{y} into multiple polygons. All locations with the same \code{id} belong to the same polygon.
#' @param id.lengths A numeric vector used to separate locations in \code{x} and \code{y} into multiple polygons. Specifies consecutive blocks of locations which make up separate polygons.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_polygon} object containing relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a BentoBox page
#' bb_pageCreate(width = 7.5, height = 6, default.units = "inches")
#'
#' ## Plot complex polygons one at a time
#' bb_plotPolygon(x = c(2.6, 4.65, 4.75, 6.05, 1.4, 1.3), y = c(2.5, 3.1, 3.5, 4, 3.15, 2.8),
#'                fill = "#4a168e", linecolor = NA)
#'
#' bb_plotPolygon(x = c(4.65, 4.75, 6.05, 5.05, 4.4), y = c(3.1, 3.5, 4, 1.45, 1.2),
#'                fill = "#9d28b0", linecolor = NA)
#'
#' ## Plot multiple triangles with different id's and colors
#' bb_plotPolygon(x = c(0.45, 6.05, 3, 3, 6.05, 5.25, 4.4, 5.05, 4.95,
#'                      1.3, 2.6, 1, 4.4, 4.95, 5, 4.95, 5, 6.25),
#'                y = c(2.85, 4, 5.55, 5.55, 4, 5.55, 1.2, 1.45, 1.1,
#'                      2.8, 2.5, 2.1, 1.2, 1.1, 0.45, 1.1, 0.45, 1.1),
#'                id = c(1,1,1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
#'                fill = c("#ce93d9", "#bb6ac9", "#4a168e", "#7b1fa0", "#bb6ac9", "#ce93d9"),
#'                linecolor = NA)
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @seealso \link[grid]{grid.polygon}
#'
#' @export
bb_plotPolygon <- function(x, y, default.units = "inches", linecolor = "black", lwd = 1, lty = 1, fill = NA, alpha = 1, id = NULL, id.lengths = NULL, params = NULL, ...){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(id)) id <- NULL
  if(missing(id.lengths)) id.lengths <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(fill)) fill <- NULL
  if(missing(lwd)) lwd <- NULL
  if(missing(lty)) lty <- NULL
  if(missing(alpha)) alpha <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if label/x/y arguments are missing (could be in object)
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL

  ## Compile all parameters into an internal object
  bb_polygonInternal <- structure(list(x = x, y = y, id = id, id.lengths = id.lengths, linecolor = linecolor, fill = fill,
                                    lwd = lwd, lty = lty, alpha = alpha, default.units = default.units), class = "bb_polygonInternal")

  bb_polygonInternal <- parseParams(bb_params = params, object_params = bb_polygonInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_polygonInternal$linecolor)) bb_polygonInternal$linecolor <- "black"
  if(is.null(bb_polygonInternal$fill)) bb_polygonInternal$fill <- NA
  if(is.null(bb_polygonInternal$lwd)) bb_polygonInternal$lwd <- 1
  if(is.null(bb_polygonInternal$lty)) bb_polygonInternal$lty <- 1
  if(is.null(bb_polygonInternal$alpha)) bb_polygonInternal$alpha <- 1
  if(is.null(bb_polygonInternal$default.units)) bb_polygonInternal$default.units <- "inches"

  ## Set gp
  bb_polygonInternal$gp <- gpar(col = bb_polygonInternal$linecolor, fill = bb_polygonInternal$fill, lwd = bb_polygonInternal$lwd,
                               lty = bb_polygonInternal$lty, alpha = bb_polygonInternal$alpha)
  bb_polygonInternal$gp <- setGP(gpList = bb_polygonInternal$gp, params = bb_polygonInternal, ...)

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_polygon <- structure(list(x = bb_polygonInternal$x, y = bb_polygonInternal$y, id = bb_polygonInternal$id, id.lengths = bb_polygonInternal$id.lengths,
                            grobs = NULL, gp = bb_polygonInternal$gp), class = "bb_polygon")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot plot polygon without a BentoBox page.")
  if(is.null(bb_polygon$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_polygon$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)


  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  if (!"unit" %in% class(bb_polygon$x)){

    if (!is.numeric(bb_polygon$x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot plot polygon.", call. = FALSE)

    }

    if (is.null(bb_polygonInternal$default.units)){

      stop("x-coordinate detected as numeric. \'default.units\' must be specified.", call. = FALSE)

    }

    bb_polygon$x <- unit(bb_polygon$x, bb_polygonInternal$default.units)

  }

  if (!"unit" %in% class(bb_polygon$y)){

    ## Check for "below" y-coords
    if (all(grepl("b", bb_polygon$y)) == TRUE){
      if (any(grepl("^[ac-zA-Z]+$", bb_polygon$y)) == TRUE){
        stop("\'below\' y-coordinate(s) detected with additional letters. Cannot parse y-coordinate(s).", call. = FALSE)
      }

      if(any(is.na(as.numeric(gsub("b","", bb_polygon$y))))){
        stop("\'below\' y-coordinate(s) does not have a numeric associated with it. Cannot parse y-coordinate(s).", call. = FALSE)
      }

      bb_polygon$y <- unit(unlist(lapply(bb_polygon$y, plot_belowY)), get("page_units", envir = bbEnv))

    } else {

      if (!is.numeric(bb_polygon$y)){

        stop("y-coordinate is neither a unit object or a numeric value. Cannot plot polygon.", call. = FALSE)

      }

      if (is.null(bb_polygonInternal$default.units)){

        stop("y-coordinate detected as numeric. \'default.units\' must be specified.", call. = FALSE)

      }

      bb_polygon$y <- unit(bb_polygon$y, bb_polygonInternal$default.units)


    }





  }

  ## Convert coordinates to page_units
  new_x <- convertX(bb_polygon$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(bb_polygon$y, unitTo = page_units, valueOnly = TRUE)

  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================

  polygon <- grid.polygon(x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), id = bb_polygon$id,
                          id.lengths = bb_polygon$id.lengths, gp = bb_polygon$gp)

  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================

  bb_polygon$grobs <- polygon

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_polygon[", polygon$name, "]"))
  invisible(bb_polygon)
}
