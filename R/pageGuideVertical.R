#' Draw a vertical guideline at a specified x-coordinate on a 
#' plotgardener page
#' 
#' @usage pageGuideVertical(
#'     x,
#'     default.units = "inches",
#'     linecolor = "grey55",
#'     params = NULL,
#'     ...
#' )
#'
#' @param x A numeric or unit object specifying x-coordinate of guide.
#' @param default.units A string indicating the default units to use
#' if \code{x} is only given as a numeric.
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
#' ## Add blue vertical guideline at x = 1.7 inches
#' pageGuideVertical(x = 1.7, linecolor = "blue")
#' @return None.
#'
#' @export
pageGuideVertical <- function(x, default.units = "inches",
                                linecolor = "grey55", params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    vguide <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "vguide"
    )

    ## Set gp
    vguide$gp <- gpar(col = vguide$linecolor)
    vguide$gp <- setGP(
        gpList = vguide$gp,
        params = vguide, ...
    )

    # =========================================================================
    # ERRORS
    # =========================================================================

    if (is.null(vguide$x)) {
        stop("argument \"x\" is missing, with no default.",
            call. = FALSE
        )
    }
    
    # =========================================================================
    # DEFAULT UNITS
    # =========================================================================

    x <- misc_defaultUnits(value = vguide$x, 
                        name = "x", 
                        default.units = vguide$default.units)
    
    # =========================================================================
    # MAKE GROB AND ASSIGN TO GTREE
    # =========================================================================

    guide <- grid.segments(
        x0 = x, x1 = x,
        y0 = unit(0, units = "npc"),
        y1 = unit(1, units = "npc"),
        gp = vguide$gp
    )
    assign("guide_grobs",
        addGrob(
            gTree = get("guide_grobs", envir = pgEnv),
            child = guide
        ),
        envir = pgEnv
    )
}
