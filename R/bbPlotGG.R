#' Plot a ggplot2 plot, gtable, or grob object in a BentoBox layout
#' 
#' @usage bbPlotGG(
#'     plot,
#'     x,
#'     y,
#'     width,
#'     height,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     params = NULL
#' )
#'
#' @param plot ggplot, gtable, or grob object.
#' @param x A numeric or unit object specifying ggplot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying ggplot y-location.
#' The character value will
#' place the ggplot y relative to the bottom of the most recently
#' plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying ggplot width.
#' @param height A numeric or unit object specifying ggplot height.
#' @param just Justification of ggplot relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use
#' if \code{x}, \code{y}, \code{width}, or \code{height} are only given
#' as numerics. Default value is \code{default.units = "inches"}.
#' @param params An optional \link[BentoBox]{bbParams} object containing
#' relevant function parameters.
#'
#' @return Returns a \code{bb_gg} object containing
#' relevant placement and \link[grid]{grob} information.
#'
#' @seealso \link[ggplot2]{ggplot}
#'
#' @examples
#' ## Create a plot using ggplot2
#' library(ggplot2)
#' p <- ggplot(mtcars) +
#'     geom_point(aes(mpg, disp))
#'
#' ## Create a BentoBox page
#' bbPageCreate(width = 4, height = 4, default.units = "inches")
#'
#' ## Place ggplot in BentoBox page
#' bbPlotGG(
#'     plot = p, x = 0.5, y = 0.5, width = 3, height = 3,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Add title
#' bbPlotText(
#'     label = "mtcars", fontsize = 14, fontface = "bold",
#'     x = 1, y = 0.35
#' )
#'
#' ## Hide page guides
#' bbPageGuideHide()
#' @export
bbPlotGG <- function(plot, x, y, width, height, just = c("left", "top"),
                    default.units = "inches", params = NULL) {


    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_ggInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_ggInternal"
    )

    ## Justification
    bb_ggInternal$just <- bb_justConversion(just = bb_ggInternal$just)
    
    # =========================================================================
    # ERRORS
    # =========================================================================

    if (is.null(bb_ggInternal$plot)) stop("argument \"plot\" is missing, ",
                                        "with no default.", call. = FALSE)
    if (is.null(bb_ggInternal$x)) stop("argument \"x\" is missing, ",
                                        "with no default.", call. = FALSE)
    if (is.null(bb_ggInternal$y)) stop("argument \"y\" is missing, ",
                                        "with no default.", call. = FALSE)
    if (is.null(bb_ggInternal$width)) stop("argument \"width\" is missing, ",
                                        "with no default.", call. = FALSE)
    if (is.null(bb_ggInternal$height)) stop("argument \"height\" is missing, ",
                                            "with no default.", call. = FALSE)

    # =========================================================================
    # INITIALIZE PLOT OBJECT
    # =========================================================================

    ggPlot <- structure(list(
        width = bb_ggInternal$width,
        height = bb_ggInternal$height,
        x = bb_ggInternal$x, y = bb_ggInternal$y,
        just = bb_ggInternal$just, grobs = NULL
    ),
    class = "bb_gg"
    )

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    ggPlot <- defaultUnits(
        object = ggPlot,
        default.units = bb_ggInternal$default.units
    )

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = ggPlot)

    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_gg",
        length(grep(
            pattern = "bb_gg",
            x = currentViewports
        )) + 1
    )

    add_bbViewport(vp_name)

    ## Make viewport for ggplot
    vp <- viewport(
        height = page_coords$height, width = page_coords$width,
        x = page_coords$x, y = page_coords$y,
        just = bb_ggInternal$just, name = vp_name
    )

    # =========================================================================
    # PRINT GGPLOT
    # =========================================================================

    ggPlot$grobs <- as.grob(bb_ggInternal$plot)
    ggPlot$grobs$vp <- vp

    grid.draw(ggPlot$grobs)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_gg[", vp_name, "]")
    invisible(ggPlot)
}
