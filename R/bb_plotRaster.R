#' Plot a raster object within a BentoBox layout
#'
#' @param image Any R object that can be coerced to a raster object.
#' @param x A numeric vector or unit object specifying raster x-locations.
#' @param y A numeric vector, unit object, or a character vector of values
#' containing a "b" combined with a numeric value specifying
#' raster y-locations.
#' The character vector will place raster y relative to the bottom
#' of the most recently plotted BentoBox plot according to the units
#' of the BentoBox page.
#' @param width A numeric vector or unit object specifying raster widths.
#' @param height A numeric vector or unit object specifying raster heights.
#' @param just Justification of text relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and
#' \code{"top"}. Default value is \code{just = "center"}.
#' @param default.units A string indicating the default units
#' to use if \code{x}, \code{y}, \code{width}, or \code{height}
#' are only given as numerics or numeric vectors.
#' Default value is \code{default.units = "inches"}.
#' @param interpolate A logical value indicating whether to linearly
#' interpolate the image. Default value is \code{interpolate = TRUE}.
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_raster} object containing
#' relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' library(png)
#'
#' ## Load images
#' edamaman <- readPNG(system.file("images",
#'                                 "bento-edamaman.png",
#'                                  package = "BentoBox"))
#' logotype <- readPNG(system.file("images",
#'                                 "bento-logotype-singleline-black.png",
#'                                  package = "BentoBox"))
#' rlogo <- readPNG(system.file("images", "Rlogo.png", package = "BentoBox"))
#'
#' ## Create page
#' bb_pageCreate(width = 5, height = 6)
#'
#' ## Plot various images
#' bb_plotRaster(image = logotype,
#'               x = 2.5, y = 0.25, width = 3.25, height = 0.5, just = "top")
#'
#' bb_plotRaster(image = edamaman,
#'               x = 2.5, y = 5.5, width = 2, height = 4, just = "bottom")
#'
#' bb_plotRaster(image = rlogo,
#'               x = 2.5, y = 1, width = 0.5, height = 0.45,
#'               just = c("right", "top"))
#'
#' ## Hide page guies
#' bb_pageGuideHide()
#'
#' @seealso \link[grid]{grid.raster}
#'
#' @export
bb_plotRaster <- function(image, x, y, width, height, just = "center",
                          default.units = "inches", interpolate = TRUE,
                          params = NULL, ...){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(interpolate)) interpolate <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if image/x/y arguments are missing (could be in object)
  if(!hasArg(image)) image <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(width)) width <- NULL
  if(!hasArg(height)) height <- NULL

  ## Compile all parameters into an internal object
  bb_rastInternal <- structure(list(image = image, x = x, y = y,
                                    width = width, height = height,
                                    just = just, interpolate = interpolate,
                                    default.units = default.units,
                                    gp = gpar()), class = "bb_rastInternal")

  bb_rastInternal <- parseParams(bb_params = params,
                                 object_params = bb_rastInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_rastInternal$just)) bb_rastInternal$just <- "center"
  if(is.null(bb_rastInternal$interpolate)) bb_rastInternal$interpolate <- TRUE
  if(is.null(bb_rastInternal$default.units)) bb_rastInternal$default.units <- "inches"

  ## Set gp
  bb_rastInternal$gp <- setGP(gpList = bb_rastInternal$gp,
                              params = bb_rastInternal, ...)

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_rast <- structure(list(image = bb_rastInternal$image,
                            x = bb_rastInternal$x, y = bb_rastInternal$y,
                            width = bb_rastInternal$width,
                            height = bb_rastInternal$height,
                            just = bb_rastInternal$just,
                            interpolate = bb_rastInternal$interpolate,
                            grobs = NULL,
                            gp = bb_rastInternal$gp), class = "bb_raster")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot plot raster without a BentoBox page.")
  if(is.null(bb_rast$image)) stop ("arguement \"image\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rast$x)) stop("argument \"x\" is missing, with no default.",
                              call. = FALSE)
  if(is.null(bb_rast$y)) stop("argument \"y\" is missing, with no default.",
                              call. = FALSE)
  if(is.null(bb_rast$width)) stop("argument \"width\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rast$height)) stop("argument \"height\" is missing, with no default.", call. = FALSE)


  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  if (!"unit" %in% class(bb_rast$x)){

    if (!is.numeric(bb_rast$x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot plot raster.", call. = FALSE)

    }

    if (is.null(bb_rastInternal$default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_rast$x <- unit(bb_rast$x, bb_rastInternal$default.units)

  }

  if (!"unit" %in% class(bb_rast$y)){

    ## Check for "below" y-coords
    if (all(grepl("b", bb_rast$y)) == TRUE){
      if (any(grepl("^[ac-zA-Z]+$", bb_rast$y)) == TRUE){
        stop("\'below\' y-coordinate(s) detected with additional letters. Cannot parse y-coordinate(s).", call. = FALSE)
      }

      if(any(is.na(as.numeric(gsub("b","", bb_rast$y))))){
        stop("\'below\' y-coordinate(s) does not have a numeric associated with it. Cannot parse y-coordinate(s).", call. = FALSE)
      }

      bb_rast$y <- unit(unlist(lapply(bb_rast$y, plot_belowY)),
                        get("page_units", envir = bbEnv))

    } else {

      if (!is.numeric(bb_rast$y)){

        stop("y-coordinate is neither a unit object or a numeric value. Cannot plot raster.", call. = FALSE)

      }

      if (is.null(bb_rastInternal$default.units)){

        stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      bb_rast$y <- unit(bb_rast$y, bb_rastInternal$default.units)


    }

  }

  if (!"unit" %in% class(bb_rast$width)){

    if (!is.numeric(bb_rast$width)){

      stop("width is neither a unit object or a numeric value. Cannot plot raster.", call. = FALSE)

    }

    if (is.null(bb_rastInternal$default.units)){

      stop("width detected as numeric.\'default.units\' must be specified.",
           call. = FALSE)

    }

    bb_rast$width <- unit(bb_rast$width, bb_rastInternal$default.units)

  }

  if (!"unit" %in% class(bb_rast$height)){

    if (!is.numeric(bb_rast$height)){

      stop("height is neither a unit object or a numeric value. Cannot plot raster.",
           call. = FALSE)

    }

    if (is.null(bb_rastInternal$default.units)){

      stop("height detected as numeric.\'default.units\' must be specified.",
           call. = FALSE)

    }

    bb_rast$height <- unit(bb_rast$height, bb_rastInternal$default.units)

  }
  ## Convert coordinates to page_units
  new_x <- convertX(bb_rast$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(bb_rast$y, unitTo = page_units, valueOnly = TRUE)
  new_width <- convertWidth(bb_rast$width,
                            unitTo = page_units, valueOnly = TRUE)
  new_height <- convertHeight(bb_rast$height,
                              unitTo = page_units, valueOnly = TRUE)

  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================

  rast <- grid.raster(image = bb_rast$image, x = unit(new_x, page_units),
                      y = unit(page_height - new_y, page_units),
                      width = unit(new_width, page_units),
                      height = unit(new_height, page_units),
                      just = bb_rast$just,
                      interpolate = bb_rast$interpolate,
                      gp = bb_rast$gp)

  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================

  bb_rast$grobs <- rast

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_raster[", rast$name, "]"))
  invisible(bb_rast)
}
