## Plot a base R plot in a plotgardener layout
#
# @param plot Plot formula of base R plotting functions.
# @param x A numeric or unit object specifying plot x-location.
# @param y A numeric, unit object, or character containing a "b"
# combined with a numeric value specifying plot y-location.
# The character value will
# place the plot y relative to the bottom of the most recently plotted
# plot according to the units of the plotgardener page.
# @param width A numeric or unit object specifying plot width.
# @param height A numeric or unit object specifying plot height.
# @param just Justification of base plot relative to its (x, y) location.
# If there are two values, the first value specifies horizontal
# justification and the second value specifies vertical justification.
# Possible string values are: \code{"left"}, \code{"right"},
# \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
# Default value is \code{just = c("left", "top")}.
# @param default.units A string indicating the default units to use
# if \code{x}, \code{y}, \code{width}, or \code{height} are only given
# as numerics. Default value is \code{default.units = "inches"}.
# @param bg Character value indicating background color.
# Default value is \code{bg = NA}.
# @param params An optional \link[plotgardener]{params} object
# containing relevant function parameters.
#
# @return Returns a \code{base} object containing
# relevant placement and \link[grid]{grob} information.
#
# @examples
# ## Define base R plot
# p <- ~ plot(1:10) + abline(v = 2)
#
# ## Create page
# pageCreate(width = 5, height = 4, default.units = "inches")
#
# ## Place base R plot in page
# plotBase(
#     plot = p,
#     x = 0.5, y = 0.5, width = 4, height = 3,
#     just = c("left", "top"), default.units = "inches"
# )
#
# ## Add title
# plotText(
#     label = "Base R Plot", fontsize = 14, fontface = "bold",
#     x = 2.75, y = 0.5
# )
#
# ## Remove page guides
# pageGuideHide()
plotBase <- function(plot, x, y, width, height, just = c("left", "top"),
                        default.units = "inches", bg = NA, params = NULL) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    baseInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "baseInternal"
    )
    
    ## Justification
    baseInternal$just <- justConversion(just = baseInternal$just)

    # =========================================================================
    # INITIALIZE PLOT OBJECT
    # =========================================================================

    base <- structure(list(
        x = baseInternal$x, y = baseInternal$y,
        width = baseInternal$width,
        height = baseInternal$height,
        just = baseInternal$just, grobs = NULL
    ),
    class = "base"
    )

    # =========================================================================
    # CALL ERRORS
    # =========================================================================

    if (is.null(baseInternal$plot)) {
        stop("argument \"plot\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(baseInternal$x)) {
        stop("argument \"x\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(baseInternal$y)) {
        stop("argument \"y\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(baseInternal$width)) {
        stop("argument \"width\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(baseInternal$height)) {
        stop("argument \"height\" is missing, with no default.",
            call. = FALSE
        )
    }

    check_page(error = "Must have a plotgardener page before
                adding a base R plot.")

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    base <- defaultUnits(
        object = base,
        default.units = baseInternal$default.units
    )

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = base)

    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "base",
        length(grep(
            pattern = "base",
            x = currentViewports
        )) + 1
    )
    addViewport(vp_name)

    ## Make viewport
    vp <- viewport(
        height = page_coords$height, width = page_coords$width,
        x = page_coords$x, y = page_coords$y,
        just = baseInternal$just, name = vp_name
    )

    # =========================================================================
    # CONVERT PLOT TO A GROB
    # =========================================================================

    gtree <- ggplotify::base2grob(baseInternal$plot)

    # =========================================================================
    # ASSIGN VIEWPORT TO GTREE
    # =========================================================================

    gtree$vp <- vp

    # =========================================================================
    # BACKGROUND AND BORDER
    # =========================================================================

    plotChildren <- gtree$childrenOrder
    backgroundGrob <- rectGrob(gp = gpar(
        fill = baseInternal$bg,
        col = NA
    ), name = "background")
    gtree <- addGrob(gTree = gtree, child = backgroundGrob)
    gtree$childrenOrder <- c("background", plotChildren)

    # =========================================================================
    # DRAW GTREE
    # =========================================================================

    grid.draw(gtree)

    # =========================================================================
    # ADD GTREE TO OBJECT
    # =========================================================================

    base$grobs <- gtree

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("base[", vp_name, "]")
    invisible(base)
}
