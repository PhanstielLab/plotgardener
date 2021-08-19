#' Plot a raster object within a plotgardener layout
#' 
#' @usage plotRaster(
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
#' of the most recently plotted plot according to the units
#' of the plotgardener page.
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
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{raster} object containing
#' relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' library(png)
#'
#' ## Load image
#' rlogo <- readPNG(system.file("images", "Rlogo.png", package = "plotgardener"))
#'
#' ## Create page
#' pageCreate(width = 5, height = 6)
#'
#' plotRaster(
#'     image = rlogo,
#'     x = 2.5, y = 1, width = 0.5, height = 0.45,
#'     just = c("right", "top")
#' )
#'
#' ## Hide page guies
#' pageGuideHide()
#' @seealso \link[grid]{grid.raster}
#'
#' @export
plotRaster <- function(image, x, y, width, height, just = "center",
                        default.units = "inches", interpolate = TRUE,
                        params = NULL, ...) {


    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    rastInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "rastInternal"
    )

    ## Set gp
    rastInternal$gp <- setGP(
        gpList = gpar(),
        params = rastInternal, ...
    )
    
    ## Justification
    rastInternal$just <- justConversion(just = rastInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    rast <- structure(list(
        image = rastInternal$image,
        x = rastInternal$x, y = rastInternal$y,
        width = rastInternal$width,
        height = rastInternal$height,
        just = rastInternal$just,
        interpolate = rastInternal$interpolate,
        grobs = NULL,
        gp = rastInternal$gp
    ), class = "raster")

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_page(error = "Cannot plot raster without a plotgardener page.")
    if (is.null(rast$image)) stop("argument \"image\" is ",
                                    "missing, with no default.", call. = FALSE)
    if (is.null(rast$x)) {
        stop("argument \"x\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(rast$y)) {
        stop("argument \"y\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(rast$width)) stop("argument \"width\" is missing, ",
                                    "with no default.", call. = FALSE)
    if (is.null(rast$height)) stop("argument \"height\" is missing, ",
                                    "with no default.", call. = FALSE)


    # =========================================================================
    # DEFINE PARAMETERS
    # =========================================================================

    ## Get page_height and its units from pgEnv
    page_height <- get("page_height", envir = pgEnv)
    page_units <- get("page_units", envir = pgEnv)
    
    rast <- defaultUnits(
        object = rast,
        default.units = rastInternal$default.units
    )
    
    ## Convert coordinates to page_units
    new_x <- convertX(rast$x, unitTo = page_units, valueOnly = TRUE)
    new_y <- convertY(rast$y, unitTo = page_units, valueOnly = TRUE)
    new_width <- convertWidth(rast$width,
        unitTo = page_units, valueOnly = TRUE
    )
    new_height <- convertHeight(rast$height,
        unitTo = page_units, valueOnly = TRUE
    )

    # =========================================================================
    # MAKE GROB
    # =========================================================================
    name <- paste0(
        "raster",
        length(grep(
            pattern = "raster",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    rast <- grid.raster(
        image = rast$image, x = unit(new_x, page_units),
        y = unit(page_height - new_y, page_units),
        width = unit(new_width, page_units),
        height = unit(new_height, page_units),
        just = rast$just,
        interpolate = rast$interpolate,
        gp = rast$gp
    )

    # =========================================================================
    # ADD GROB TO OBJECT
    # =========================================================================

    rast$grobs <- rast

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("raster[", name, "]")
    invisible(rast)
}
