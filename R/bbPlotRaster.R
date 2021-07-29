#' Plot a raster object within a BentoBox layout
#' 
#' @usage bbPlotRaster(
#'     image,
#'     x,
#'     y,
#'     width,
#'     height,
#'     just = "center",
#'     default.units = "inches",
#'     interpolate = TRUE,
#'     params = NULL,
#'     ...
#' )
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
#' @param params An optional \link[BentoBox]{bbParams} object
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
#'     "bento-edamaman.png",
#'     package = "BentoBox"
#' ))
#' logotype <- readPNG(system.file("images",
#'     "bento-logotype-singleline-black.png",
#'     package = "BentoBox"
#' ))
#' rlogo <- readPNG(system.file("images", "Rlogo.png", package = "BentoBox"))
#'
#' ## Create page
#' bbPageCreate(width = 5, height = 6)
#'
#' ## Plot various images
#' bbPlotRaster(
#'     image = logotype,
#'     x = 2.5, y = 0.25, width = 3.25, height = 0.5, just = "top"
#' )
#'
#' bbPlotRaster(
#'     image = edamaman,
#'     x = 2.5, y = 5.5, width = 2, height = 4, just = "bottom"
#' )
#'
#' bbPlotRaster(
#'     image = rlogo,
#'     x = 2.5, y = 1, width = 0.5, height = 0.45,
#'     just = c("right", "top")
#' )
#'
#' ## Hide page guies
#' bbPageGuideHide()
#' @seealso \link[grid]{grid.raster}
#'
#' @export
bbPlotRaster <- function(image, x, y, width, height, just = "center",
                        default.units = "inches", interpolate = TRUE,
                        params = NULL, ...) {


    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_rastInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_rastInternal"
    )

    ## Set gp
    bb_rastInternal$gp <- setGP(
        gpList = gpar(),
        params = bb_rastInternal, ...
    )
    
    ## Justification
    bb_rastInternal$just <- bb_justConversion(just = bb_rastInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    bb_rast <- structure(list(
        image = bb_rastInternal$image,
        x = bb_rastInternal$x, y = bb_rastInternal$y,
        width = bb_rastInternal$width,
        height = bb_rastInternal$height,
        just = bb_rastInternal$just,
        interpolate = bb_rastInternal$interpolate,
        grobs = NULL,
        gp = bb_rastInternal$gp
    ), class = "bb_raster")

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_bbpage(error = "Cannot plot raster without a BentoBox page.")
    if (is.null(bb_rast$image)) stop("arguement \"image\" is ",
                                    "missing, with no default.", call. = FALSE)
    if (is.null(bb_rast$x)) {
        stop("argument \"x\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(bb_rast$y)) {
        stop("argument \"y\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(bb_rast$width)) stop("argument \"width\" is missing, ",
                                    "with no default.", call. = FALSE)
    if (is.null(bb_rast$height)) stop("argument \"height\" is missing, ",
                                    "with no default.", call. = FALSE)


    # =========================================================================
    # DEFINE PARAMETERS
    # =========================================================================

    ## Get page_height and its units from bbEnv through bb_makepage
    page_height <- get("page_height", envir = bbEnv)
    page_units <- get("page_units", envir = bbEnv)
    
    bb_rast <- defaultUnits(
        object = bb_rast,
        default.units = bb_rastInternal$default.units
    )
    
    ## Convert coordinates to page_units
    new_x <- convertX(bb_rast$x, unitTo = page_units, valueOnly = TRUE)
    new_y <- convertY(bb_rast$y, unitTo = page_units, valueOnly = TRUE)
    new_width <- convertWidth(bb_rast$width,
        unitTo = page_units, valueOnly = TRUE
    )
    new_height <- convertHeight(bb_rast$height,
        unitTo = page_units, valueOnly = TRUE
    )

    # =========================================================================
    # MAKE GROB
    # =========================================================================
    name <- paste0(
        "bb_raster",
        length(grep(
            pattern = "bb_raster",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    rast <- grid.raster(
        image = bb_rast$image, x = unit(new_x, page_units),
        y = unit(page_height - new_y, page_units),
        width = unit(new_width, page_units),
        height = unit(new_height, page_units),
        just = bb_rast$just,
        interpolate = bb_rast$interpolate,
        gp = bb_rast$gp
    )

    # =========================================================================
    # ADD GROB TO OBJECT
    # =========================================================================

    bb_rast$grobs <- rast

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_raster[", name, "]")
    invisible(bb_rast)
}
