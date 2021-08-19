#' Plot a rectangle within a plotgardener layout
#' 
#' @usage plotRect(
#'     x,
#'     y,
#'     width,
#'     height,
#'     just = "center",
#'     default.units = "inches",
#'     linecolor = "black",
#'     lwd = 1,
#'     lty = 1,
#'     fill = NA,
#'     alpha = 1,
#'     params = NULL,
#'     ...
#' )
#'
#' @param x A numeric vector or unit object specifying rectangle x-locations.
#' @param y A numeric vector, unit object, or a character vector of values
#' containing a "b" combined with a numeric value specifying
#' rectangle y-locations.
#' The character vector will place rectangle y-locations relative to
#' the bottom of the most recently plotted plot according to
#' the units of the plotgardener page.
#' @param width A numeric vector or unit object specifying rectangle widths.
#' @param height A numeric vector or unit object specifying rectangle heights.
#' @param just Justification of rectangle relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal justification
#' and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = "center"}.
#' @param default.units A string indicating the default units to use if
#' \code{x}, \code{y}, \code{width}, and \code{height} are only given as
#' numerics or numeric vectors.
#' Default value is \code{default.units = "inches"}.
#' @param linecolor A character value specifying rectangle line color.
#' Default value is \code{linecolor = "black"}.
#' @param lwd A numeric specifying rectangle line width.
#' Default value is \code{lwd = 1}.
#' @param lty A numeric specifying rectangle line type.
#' Default value is \code{lty = 1}.
#' @param fill A character value specifying rectangle fill color.
#' Default value is \code{fill = NA}.
#' @param alpha Numeric value specifying color transparency.
#' Default value is \code{alpha = 1}.
#' @param params An optional \link[plotgardener]{params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{rect} object containing
#' relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' pageCreate(width = 7.5, height = 6, default.units = "inches")
#'
#' ## Plot one rectangle with no fill
#' plotRect(
#'     x = 0.5, y = 0.5, width = 3, height = 3,
#'     just = c("left", "top"), default.units = "inches",
#'     lwd = 2, fill = NA
#' )
#'
#'
#' ## Plot two rectangles with same width and height at different locations
#' plotRect(
#'     x = 4, y = c(0.5, 2.25), width = 3, height = 1.25,
#'     just = c("left", "top"), default.units = "inches",
#'     fill = "#7ecdbb"
#' )
#'
#' ## Plot two rectangles with different widths, heights,
#' ## locations, and colors
#' plotRect(
#'     x = 3.75, y = c(4, 5.25), width = c(6.5, 4.5),
#'     height = c(1, 0.25),
#'     just = "top", default.units = "inches",
#'     fill = c("#7ecdbb", "#37a7db"), linecolor = NA, alpha = 0.4
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @seealso \link[grid]{grid.rect}
#'
#' @export
plotRect <- function(x, y, width, height, just = "center",
                        default.units = "inches", linecolor = "black",
                        lwd = 1, lty = 1, fill = NA, alpha = 1,
                        params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    rectInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "rectInternal"
    )

    ## Set gp
    rectInternal$gp <- gpar(
        col = rectInternal$linecolor,
        fill = rectInternal$fill,
        lwd = rectInternal$lwd,
        lty = rectInternal$lty,
        alpha = rectInternal$alpha
    )
    rectInternal$gp <- setGP(
        gpList = rectInternal$gp,
        params = rectInternal, ...
    )
    
    ## Justification
    rectInternal$just <- justConversion(just = rectInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    rect <- structure(list(
        x = rectInternal$x, y = rectInternal$y,
        width = rectInternal$width,
        height = rectInternal$height,
        just = rectInternal$just,
        grobs = NULL, gp = rectInternal$gp
    ),
    class = "rect"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_page(error = "Cannot plot rectangle without a `plotgardener` page.")
    if (is.null(rect$x)) {
        stop("argument \"x\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(rect$y)) {
        stop("argument \"y\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(rect$width)) {
        stop("argument \"width\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(rect$height)) {
        stop("argument \"height\" is missing, with no default.",
            call. = FALSE
        )
    }
    
    checkColorby(fill = rect$gp$fill,
                    colorby = FALSE)

    # =========================================================================
    # DEFINE PARAMETERS
    # =========================================================================

    ## Get page_height and its units from pgEnv
    page_height <- get("page_height", envir = pgEnv)
    page_units <- get("page_units", envir = pgEnv)
    
    rect$x <- misc_defaultUnits(
        value = rect$x,
        name = "x",
        default.units = rectInternal$default.units
    )
    rect$y <- misc_defaultUnits(
        value = rect$y,
        name = "y",
        default.units = rectInternal$default.units
    )
    rect$width <- misc_defaultUnits(
        value = rect$width,
        name = "width",
        default.units = rectInternal$default.units
    )
    rect$height <- misc_defaultUnits(
        value = rect$height,
        name = "height",
        default.units = rectInternal$default.units
    )

    ## Convert coordinates to page_units
    new_x <- convertX(rect$x, unitTo = page_units, valueOnly = TRUE)
    new_y <- convertY(rect$y, unitTo = page_units, valueOnly = TRUE)
    new_width <- convertWidth(rect$width,
        unitTo = page_units, valueOnly = TRUE
    )
    new_height <- convertHeight(rect$height,
        unitTo = page_units, valueOnly = TRUE
    )

    # =========================================================================
    # MAKE GROB
    # =========================================================================
    name <- paste0(
        "rect",
        length(grep(
            pattern = "rect",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    rect <- grid.rect(
        x = unit(new_x, page_units),
        y = unit(page_height - new_y, page_units),
        width = unit(new_width, page_units),
        height = unit(new_height, page_units),
        just = rect$just, gp = rect$gp,
        name = name
    )

    # =========================================================================
    # ADD GROB TO OBJECT
    # =========================================================================

    rect$grobs <- rect

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("rect[", rect$name, "]")
    invisible(rect)
}
