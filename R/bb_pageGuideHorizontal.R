#' Draw a horizontal guideline at a specified y-coordinate
#' on a BentoBox page
#' 
#' @usage bb_pageGuideHorizontal(
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
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @examples
#' ## Create a BentoBox page
#' bb_pageCreate(width = 6, height = 5, default.units = "inches")
#'
#' ## Add red horizontal guideline at y = 2.5 inches
#' bb_pageGuideHorizontal(y = 2.5, linecolor = "red")
#' @return None.
#'
#' @export
bb_pageGuideHorizontal <- function(y, default.units = "inches",
                                linecolor = "grey55", params = NULL, ...) {


    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_hguide <- parseParams(params = params, 
                            defaultArgs = formals(eval(match.call()[[1]])),
                            declaredArgs = lapply(match.call()[-1], 
                                                eval.parent, n = 2),
                            class = "bb_hguide")

    ## Set gp
    bb_hguide$gp <- gpar(col = bb_hguide$linecolor)
    bb_hguide$gp <- setGP(gpList = bb_hguide$gp, params = bb_hguide, ...)
    # =========================================================================
    # ERRORS
    # =========================================================================

    if (is.null(bb_hguide$y)) {
        stop("argument \"y\" is missing, with no default.",
            call. = FALSE
        )
    }

    # =========================================================================
    # DEFAULT UNITS
    # =========================================================================

    y <- misc_defaultUnits(value = bb_hguide$y, 
                        name = "y",
                        default.units = bb_hguide$default.units,
                        funName = "bb_pageGuideHorizontal",
                        yBelow = FALSE)

    # =========================================================================
    # MAKE GROB AND ASSIGN TO GTREE
    # =========================================================================

    y <- convertY(y,
        unitTo = get("page_units", envir = bbEnv),
        valueOnly = TRUE
    )

    guide <- grid.segments(
        x0 = unit(0, units = "npc"),
        x1 = unit(1, units = "npc"),
        y0 = get("page_height", envir = bbEnv) - y,
        y1 = get("page_height", envir = bbEnv) - y,
        default.units = get("page_units", envir = bbEnv),
        gp = bb_hguide$gp
    )
    assign("guide_grobs",
        addGrob(
            gTree = get("guide_grobs", envir = bbEnv),
            child = guide
        ),
        envir = bbEnv
    )
}
