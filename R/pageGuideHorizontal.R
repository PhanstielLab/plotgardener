#' Draw a horizontal guideline at a specified y-coordinate
#' on a plotgardener page
#' 
#' @usage pageGuideHorizontal(
#'     y,
#'     default.units = "inches",
#'     linecolor = "grey55",
#'     params = NULL,
#'     ...
#' )
#'
#' @param y A numeric or unit object specifying y-coordinate of guide.
#' @param default.units A string indicating the default units to use
#' if \code{y} is only given as a numeric.
#' Default value is \code{default.units = "inches"}.
#' @param linecolor Character value indicating color of guideline.
#' Default value is \code{linecolor = "grey55"}.
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @examples
#' ## Create a page
#' pageCreate(width = 6, height = 5, default.units = "inches")
#'
#' ## Add red horizontal guideline at y = 2.5 inches
#' pageGuideHorizontal(y = 2.5, linecolor = "red")
#' @return None.
#'
#' @export
pageGuideHorizontal <- function(y, default.units = "inches",
                                linecolor = "grey55", params = NULL, ...) {


    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    hguide <- parseParams(params = params, 
                            defaultArgs = formals(eval(match.call()[[1]])),
                            declaredArgs = lapply(match.call()[-1], 
                                                eval.parent, n = 2),
                            class = "hguide")

    ## Set gp
    hguide$gp <- gpar(col = hguide$linecolor)
    hguide$gp <- setGP(gpList = hguide$gp, params = hguide, ...)
    # =========================================================================
    # ERRORS
    # =========================================================================

    if (is.null(hguide$y)) {
        stop("argument \"y\" is missing, with no default.",
            call. = FALSE
        )
    }

    # =========================================================================
    # DEFAULT UNITS
    # =========================================================================

    y <- misc_defaultUnits(value = hguide$y, 
                        name = "y",
                        default.units = hguide$default.units,
                        funName = "pageGuideHorizontal",
                        yBelow = FALSE)

    # =========================================================================
    # MAKE GROB AND ASSIGN TO GTREE
    # =========================================================================

    y <- convertY(y,
        unitTo = get("page_units", envir = pgEnv),
        valueOnly = TRUE
    )

    guide <- grid.segments(
        x0 = unit(0, units = "npc"),
        x1 = unit(1, units = "npc"),
        y0 = get("page_height", envir = pgEnv) - y,
        y1 = get("page_height", envir = pgEnv) - y,
        default.units = get("page_units", envir = pgEnv),
        gp = hguide$gp
    )
    assign("guide_grobs",
        addGrob(
            gTree = get("guide_grobs", envir = pgEnv),
            child = guide
        ),
        envir = pgEnv
    )
}
