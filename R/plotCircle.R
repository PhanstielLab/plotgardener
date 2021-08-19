#' Plot a circle within a plotgardener layout
#' 
#' @usage plotCircle(
#'     x,
#'     y,
#'     r,
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
#' @param x A numeric vector or unit object specifying circle
#' x-locations relative to center.
#' @param y A numeric vector, unit object, or a character vector
#' of values containing a "b" combined with a numeric value
#' specifying circle y-locations relative to center.
#' The character vector will place circle y-locations relative to
#' the bottom of the most recently plotted plot according
#' to the units of the plotgardener page.
#' @param r A numeric vector or unit object specifying radii.
#' @param default.units A string indicating the default units to use
#' if \code{r}, \code{x}, or \code{y} are only given as numerics or
#' numeric vectors. Default value is \code{default.units = "inches"}.
#' @param linecolor A character value specifying circle line color.
#' Default value is \code{linecolor = "black"}.
#' @param lwd A numeric specifying circle line width.
#' Default value is \code{lwd = 1}.
#' @param lty A numeric specifying circle line type.
#' Default value is \code{lty = 1}.
#' @param fill A character value specifying circle fill color.
#' Default value is \code{fill = NA}.
#' @param alpha Numeric value specifying color transparency.
#' Default value is \code{alpha = 1}.
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{circle} object containing
#' relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' pageCreate(width = 2, height = 2, default.units = "inches")
#'
#' ## Plot two circles, one at a time
#' plotCircle(
#'     x = 0.6, y = 0.5, r = 0.1, fill = "black",
#'     default.units = "inches"
#' )
#' plotCircle(
#'     x = 1.4, y = 0.5, r = 0.1, fill = "black",
#'     default.units = "inches"
#' )
#'
#' ## Plot a vector of circles
#' xVals <- 1 + (0.5 * cos(seq(0, pi, pi / 8)))
#' yVals <- 1 + (0.5 * sin(seq(0, pi, pi / 8)))
#' plotCircle(x = xVals, y = yVals, r = 0.05, default.units = "inches")
#'
#' ## Hide page guides
#' pageGuideHide()
#' @seealso \link[grid]{grid.circle}
#'
#' @export
plotCircle <- function(x, y, r, default.units = "inches",
                        linecolor = "black", lwd = 1, lty = 1,
                        fill = NA, alpha = 1, params = NULL, ...) {


    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    circleInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "circleInternal"
    )

    ## Set gp
    circleInternal$gp <- gpar(
        col = circleInternal$linecolor,
        fill = circleInternal$fill,
        lwd = circleInternal$lwd,
        lty = circleInternal$lty,
        alpha = circleInternal$alpha
    )
    circleInternal$gp <- setGP(
        gpList = circleInternal$gp,
        params = circleInternal, ...
    )


    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    circle <- structure(list(
        x = circleInternal$x,
        y = circleInternal$y,
        r = circleInternal$r, grobs = NULL,
        gp = circleInternal$gp
    ), class = "circle")

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_page(error = "Cannot plot circle without a `plotgardener` page.")
    if (is.null(circle$x)) {
        stop("argument \"x\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(circle$y)) {
        stop("argument \"y\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(circle$r)) {
        stop("argument \"r\" is missing, with no default.",
            call. = FALSE
        )
    }
    
    checkColorby(fill = circle$gp$fill,
                colorby = FALSE)

    # =========================================================================
    # DEFINE PARAMETERS
    # =========================================================================

    ## Get page_height and its units from pgEnv
    page_height <- get("page_height", envir = pgEnv)
    page_units <- get("page_units", envir = pgEnv)
    
    circle$x <- misc_defaultUnits(value = circle$x,
                                    name = "x",
                                    default.units = 
                                        circleInternal$default.units)
    
    circle$y <- misc_defaultUnits(value = circle$y,
                                    name = "y",
                                    default.units = 
                                        circleInternal$default.units)
    
    circle$r <- misc_defaultUnits(value = circle$r,
                                    name = "r",
                                    default.units = 
                                        circleInternal$default.units)

    ## Convert coordinates to page_units
    new_x <- convertX(circle$x, unitTo = page_units, valueOnly = TRUE)
    new_y <- convertY(circle$y, unitTo = page_units, valueOnly = TRUE)
    new_r <- convertUnit(circle$r, unitTo = page_units, valueOnly = TRUE)

    # =========================================================================
    # MAKE GROB
    # =========================================================================
    name <- paste0(
        "circle",
        length(grep(
            pattern = "circle",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    circle <- grid.circle(
        x = unit(new_x, page_units),
        y = unit(page_height - new_y, page_units),
        r = unit(new_r, page_units),
        gp = circle$gp,
        name = name
    )

    # =========================================================================
    # ADD GROB TO OBJECT
    # =========================================================================

    circle$grobs <- circle

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("circle[", circle$name, "]")
    invisible(circle)
}
