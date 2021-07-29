#' Draw a vertical guideline at a specified x-coordinate on a BentoBox page
#' 
#' @usage bbPageGuideVertical(
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
#' @param params An optional \link[BentoBox]{bbParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @examples
#' ## Create a BentoBox page
#' bbPageCreate(width = 6, height = 5, default.units = "inches")
#'
#' ## Add blue vertical guideline at x = 1.7 inches
#' bbPageGuideVertical(x = 1.7, linecolor = "blue")
#' @return None.
#'
#' @export
bbPageGuideVertical <- function(x, default.units = "inches",
                                linecolor = "grey55", params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_vguide <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_vguide"
    )

    ## Set gp
    bb_vguide$gp <- gpar(col = bb_vguide$linecolor)
    bb_vguide$gp <- setGP(
        gpList = bb_vguide$gp,
        params = bb_vguide, ...
    )

    # =========================================================================
    # ERRORS
    # =========================================================================

    if (is.null(bb_vguide$x)) {
        stop("argument \"x\" is missing, with no default.",
            call. = FALSE
        )
    }
    
    # =========================================================================
    # DEFAULT UNITS
    # =========================================================================

    x <- misc_defaultUnits(value = bb_vguide$x, 
                        name = "x", 
                        default.units = bb_vguide$default.units)
    
    # =========================================================================
    # MAKE GROB AND ASSIGN TO GTREE
    # =========================================================================

    guide <- grid.segments(
        x0 = x, x1 = x,
        y0 = unit(0, units = "npc"),
        y1 = unit(1, units = "npc"),
        gp = bb_vguide$gp
    )
    assign("guide_grobs",
        addGrob(
            gTree = get("guide_grobs", envir = bbEnv),
            child = guide
        ),
        envir = bbEnv
    )
}
