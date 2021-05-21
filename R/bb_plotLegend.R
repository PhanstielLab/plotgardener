#' Plot a legend
#'
#' @param legend A character or expression vector to appear in the legend.
#' @param fill If specified, this argument will produce boxes filled with
#' the specified colors to appear beside the legend text.
#' @param pch The plotting symbols appearing in the legend, as a
#' numeric vector.
#' @param lty The line types for lines appearing in the legend.
#' @param orientation A string specifying legend orientation.
#' Default value is \code{orientation = "v"}. Options are:
#' \itemize{
#' \item{\code{"v"}: }{Vertical legend orientation.}
#' \item{\code{"h"}: }{Horizontal legend orientation.}
#' }
#' @param title A character value giving a title to be placed at
#' the top of the legend.
#' @param fontsize A numeric specifying text fontsize in points.
#' Default value is \code{fontsize = 10}.
#' @param border Logical value indicating whether to add a border
#' around heatmap legend. Default value is \code{border = TRUE}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param x A numeric or unit object specifying legend x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying legend y-location.
#' The character value will
#' place the legend y relative to the bottom of the most recently
#' plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying legend width.
#' @param height A numeric or unit object specifying legend height.
#' @param just Justification of legend relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use
#' if \code{x}, \code{y}, \code{width}, or \code{height} are only given
#' as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should
#' be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_legend} object containing relevant
#' placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Load BED data
#' data("bb_bedData")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 7.5, height = 4, default.units = "inches")
#'
#' ## Plot a pileup plot, coloring elements by strand
#' pileupPlot <- bb_plotRanges(
#'     data = bb_bedData, chrom = "chr21",
#'     chromstart = 29072500, chromend = 29075000,
#'     fill = c("steel blue", "light salmon"),
#'     colorby = colorby("strand"),
#'     x = 0.5, y = 3.5, width = 6.5, height = 3.5,
#'     just = c("left", "bottom"),
#'     default.units = "inches"
#' )
#'
#' ## Add a legend depicting strand colors
#' legendPlot <- bb_plotLegend(
#'     legend = c("- strand", "+ strand"),
#'     fill = c("steel blue", "light salmon"),
#'     border = FALSE,
#'     x = 5, y = 0.5, width = 1.5, height = 0.7,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bb_annoGenomeLabel(
#'     plot = pileupPlot, x = 0.5, y = 3.5,
#'     just = c("left", "top")
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @export
bb_plotLegend <- function(legend, fill = NULL, pch = NULL, lty = NULL,
                        orientation = "v", title = NULL, fontsize = 10,
                        border = TRUE,
                        bg = NA, x = NULL, y = NULL, width = NULL,
                        height = NULL, just = c("left", "top"),
                        default.units = "inches", draw = TRUE,
                        params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    ## Check which defaults are not overwritten and set to NULL
    if (missing(orientation)) orientation <- NULL
    if (missing(fontsize)) fontsize <- NULL
    if (missing(border)) border <- NULL
    if (missing(bg)) bg <- NULL
    if (missing(just)) just <- NULL
    if (missing(default.units)) default.units <- NULL
    if (missing(draw)) draw <- NULL

    ## Check if legend argument is missing (could be in object)
    if (!hasArg(legend)) legend <- NULL

    ## Compile all parameters into an internal object
    bb_legInternal <- structure(list(
        legend = legend, orientation = orientation,
        fill = fill, pch = pch, lty = lty,
        title = title, fontsize = fontsize,
        border = border, bg = bg, x = x, y = y,
        width = width, height = height,
        just = just, default.units = default.units,
        draw = draw, gp = gpar()
    ),
    class = "bb_legInternal"
    )

    bb_legInternal <- parseParams(
        bb_params = params,
        object_params = bb_legInternal
    )

    ## For any defaults that are still NULL, set back to default
    if (is.null(bb_legInternal$orientation)) bb_legInternal$orientation <- "v"
    if (is.null(bb_legInternal$fontsize)) bb_legInternal$fontsize <- 10
    if (is.null(bb_legInternal$border)) bb_legInternal$border <- TRUE
    if (is.null(bb_legInternal$bg)) bb_legInternal$bg <- NA
    if (is.null(bb_legInternal$just)) bb_legInternal$just <- c("left", "top")
    if (is.null(bb_legInternal$default.units)) {
        bb_legInternal$default.units <- "inches"
    }
    if (is.null(bb_legInternal$draw)) bb_legInternal$draw <- TRUE

    ## Set gp
    bb_legInternal$gp <- setGP(
        gpList = bb_legInternal$gp,
        params = bb_legInternal, ...
    )

    ## Reset lty
    if (is.null(bb_legInternal$gp$lty)) {
        bb_legInternal$gp$lty <- NULL
    }


    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    legend_plot <- structure(list(
        x = bb_legInternal$x, y = bb_legInternal$y,
        width = bb_legInternal$width,
        height = bb_legInternal$height,
        just = bb_legInternal$just, grobs = NULL
    ),
    class = "bb_legend"
    )
    attr(x = legend_plot, which = "plotted") <- bb_legInternal$draw

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(bb_legInternal$legend)) stop("argument \"legend\" is missing,
                                            with no default.", call. = FALSE)

    check_placement(object = legend_plot)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    legend_plot <- defaultUnits(
        object = legend_plot,
        default.units = bb_legInternal$default.units
    )

    textHeight <- heightDetails(textGrob(
        label = "A",
        gp = gpar(fontsize = bb_legInternal$fontsize)
    ))
    textGrobs <- lapply(bb_legInternal$legend, textGrob,
        gp = gpar(fontsize = bb_legInternal$fontsize)
    )
    textWidths <- lapply(textGrobs, widthDetails)

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_legend",
        length(grep(
            pattern = "bb_legend",
            x = currentViewports
        )) + 1
    )

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(legend_plot$x) & is.null(legend_plot$y)) {
        vp <- viewport(
            height = unit(0.125, "snpc"), width = unit(0.20, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            just = "center",
            yscale = c(0, 0.125),
            xscale = c(0, 0.20),
            name = vp_name
        )
        pushViewport(vp)
        textHeight <- convertHeight(textHeight,
            unitTo = "native", valueOnly = TRUE
        )
        textWidths <- lapply(textWidths, convertWidth,
            unitTo = "native", valueOnly = TRUE
        )
        upViewport()

        height <- 0.125
        width <- 0.20

        if (bb_legInternal$draw == TRUE) {
            vp$name <- "bb_legend1"
            grid.newpage()
        }
    } else {

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = legend_plot)
        add_bbViewport(vp_name)

        height <- convertHeight(page_coords$height,
            unitTo = get("page_units", envir = bbEnv),
            valueOnly = TRUE
        )
        width <- convertWidth(page_coords$width,
            unitTo = get("page_units", envir = bbEnv),
            valueOnly = TRUE
        )
        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            just = bb_legInternal$just,
            yscale = c(0, height),
            xscale = c(0, width),
            name = vp_name
        )

        textHeight <- convertHeight(textHeight,
            unitTo = get("page_units", envir = bbEnv),
            valueOnly = TRUE
        )
        textWidths <- lapply(textWidths, convertWidth,
            unitTo = get("page_units", envir = bbEnv),
            valueOnly = TRUE
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS
    # =========================================================================

    assign("legend_grobs", gTree(vp = vp), envir = bbEnv)

    # =========================================================================
    # GROBS
    # =========================================================================

    ## Border
    if (bb_legInternal$border == TRUE) {
        bb_legInternal$gp$fill <- bb_legInternal$bg
        border <- rectGrob(gp = bb_legInternal$gp)
    } else {
        bb_legInternal$gp$fill <- bb_legInternal$bg
        bb_legInternal$gp$col <- NA
        border <- rectGrob(gp = bb_legInternal$gp)
    }

    assign("legend_grobs",
        addGrob(get("legend_grobs", envir = bbEnv),
            child = border
        ),
        envir = bbEnv
    )

    ## Title
    if (!is.null(bb_legInternal$title)) {
        if (bb_legInternal$orientation == "h") {
            remainingSpace <- height - textHeight * 2
            spaceNo <- 3
        } else {
            remainingSpace <- height - textHeight *
                (length(bb_legInternal$legend) + 1)
            spaceNo <- length(bb_legInternal$legend) + 2
        }

        ## Remove lty
        lty <- bb_legInternal$gp$lty
        bb_legInternal$gp$lty <- NULL
        ## Remove cex for text
        cex <- bb_legInternal$gp$cex
        bb_legInternal$gp$cex <- NULL

        bb_legInternal$gp$col <- NA

        titleGrob <- textGrob(
            label = bb_legInternal$title,
            x = unit(0.5, "npc"),
            y = height - remainingSpace / spaceNo,
            just = "top", gp = bb_legInternal$gp,
            default.units = "native"
        )
        assign("legend_grobs",
            addGrob(get("legend_grobs", envir = bbEnv),
                child = titleGrob
            ),
            envir = bbEnv
        )

        bb_legInternal$gp$lty <- lty
        bb_legInternal$gp$cex <- cex
    } else {
        if (bb_legInternal$orientation == "h") {
            remainingSpace <- height - textHeight
            spaceNo <- 2
        } else {
            remainingSpace <- height - textHeight *
                length(bb_legInternal$legend)
            spaceNo <- length(legend) + 1
        }
    }


    ## Colors
    ## Only take the first number of colors for the length of legend
    fillcolors <- bb_legInternal$fill[seq(1, length(bb_legInternal$legend))]
    ## Spacing and label coordinates
    spaceHeight <- remainingSpace / spaceNo

    if (bb_legInternal$orientation == "h") {
        ycoords <- spaceHeight
        remainingWidth <- width - (sum(unlist(textWidths)) +
            textHeight *
                length(bb_legInternal$legend))
        width_spaces <- 2 * length(bb_legInternal$legend) + 1
        widthSpace <- remainingWidth / width_spaces
        if (widthSpace < 0) {
            widthSpace <- 0.1 * textHeight
        }

        addTexts <- function(factor, textWidths) {
            totalText <- sum(textWidths[seq(0, factor)])
            return(totalText)
        }

        textWidth <- c(
            0,
            unlist(textWidths)[seq(
                1, (length(bb_legInternal$legend) - 1)
            )]
        )
        df <- as.data.frame(cbind(
            "factor" = seq(0, (length(bb_legInternal$legend) - 1)),
            "firstSpace" = rep(widthSpace, length(bb_legInternal$legend))
        ))
        df$totalText <- unlist(lapply(df$factor, addTexts,
            textWidths = unlist(textWidths)
        ))

        xcoords <- df$firstSpace + df$factor *
            (2 * widthSpace + textHeight) + df$totalText

        xpch <- xcoords + 0.5 * textHeight
        ypch <- rep(
            spaceHeight + 0.5 * textHeight,
            length(bb_legInternal$legend)
        )

        x0s <- xcoords
        y0s <- spaceHeight + 0.5 * textHeight
        x1s <- x0s + textHeight
        y1s <- spaceHeight + 0.5 * textHeight
    } else {
        ycoords <- seq(
            from = spaceHeight, to = height,
            by = (spaceHeight + textHeight)
        )[seq(1, length(bb_legInternal$legend))]
        xcoords <- textHeight * 2

        xpch <- rep(textHeight * 2, length(bb_legInternal$legend))
        ypch <- ycoords + 0.5 * textHeight

        x0s <- textHeight
        y0s <- ycoords + 0.5 * textHeight
        x1s <- 3 * textHeight
        y1s <- ycoords + 0.5 * textHeight

        fillcolors <- rev(fillcolors)
    }




    ## Symbols
    if (!is.null(bb_legInternal$pch)) {
        ## Only take the first number of symbols for the length of legend
        pchs <- bb_legInternal$pch[seq(1, length(bb_legInternal$legend))]
        if (bb_legInternal$orientation == "v") {
            pchs <- rev(pchs)
        }

        bb_legInternal$gp$col <- fillcolors

        labelSymbols <- pointsGrob(
            x = xpch, y = ypch,
            pch = pchs,
            gp = bb_legInternal$gp,
            default.units = "native"
        )
        assign("legend_grobs",
            addGrob(get("legend_grobs", envir = bbEnv),
                child = labelSymbols
            ),
            envir = bbEnv
        )

        ## Lines
    } else if (!is.null(lty)) {
        ## Only take the first number of symbols for the length of legend
        ltys <- lty[seq(1, length(legend))]
        lty <- bb_legInternal$gp$lty
        bb_legInternal$gp$lty <- ltys
        bb_legInternal$gp$col <- fillcolors
        labelLines <- segmentsGrob(
            x0 = x0s, y0 = y0s,
            x1 = x1s, y1 = y1s,
            gp = bb_legInternal$gp,
            default.units = "native"
        )
        assign("legend_grobs",
            addGrob(get("legend_grobs", envir = bbEnv),
                child = labelLines
            ),
            envir = bbEnv
        )

        bb_legInternal$gp$lty <- lty
    } else {
        labelColors <- rectGrob(
            x = xcoords, y = ycoords,
            width = textHeight, height = textHeight,
            just = c("left", "bottom"),
            gp = gpar(fill = fillcolors, col = NA),
            default.units = "native"
        )
        assign("legend_grobs",
            addGrob(get("legend_grobs", envir = bbEnv),
                child = labelColors
            ),
            envir = bbEnv
        )
    }

    bb_legInternal$gp$fontsize <- bb_legInternal$fontsize

    if (is.null(bb_legInternal$gp$fontcolor)) {
        bb_legInternal$gp$fontcolor <- "black"
        bb_legInternal$gp$col <- "black"
    } else {
        bb_legInternal$gp$col <- bb_legInternal$gp$fontcolor
    }
    bb_legInternal$gp$lty <- NULL
    bb_legInternal$gp$cex <- NULL

    ## Text labels
    if (bb_legInternal$orientation == "h") {
        labelText <- textGrob(
            label = bb_legInternal$legend,
            x = xcoords + textHeight + widthSpace,
            y = ycoords + 0.5 * textHeight,
            just = c("left", "center"),
            gp = bb_legInternal$gp, default.units = "native"
        )
    } else {
        labelText <- textGrob(
            label = rev(bb_legInternal$legend),
            x = xcoords + 2 * textHeight,
            y = ycoords + 0.5 * textHeight,
            just = c("left", "center"),
            gp = bb_legInternal$gp, default.units = "native"
        )
    }

    assign("legend_grobs",
        addGrob(get("legend_grobs", envir = bbEnv),
            child = labelText
        ),
        envir = bbEnv
    )


    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (bb_legInternal$draw == TRUE) {
        grid.draw(get("legend_grobs", envir = bbEnv))
    }
    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    legend_plot$grobs <- get("legend_grobs", envir = bbEnv)

    # ========================================================================
    # RETURN OBJECT
    # ========================================================================

    message("bb_legend[", vp$name, "]")
    invisible(legend_plot)
}
