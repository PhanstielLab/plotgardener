#' Plot a legend
#' 
#' @usage plotLegend(
#'     legend,
#'     fill = NULL,
#'     pch = NULL,
#'     lty = NULL,
#'     orientation = "v",
#'     title = NULL,
#'     fontsize = 10,
#'     border = TRUE,
#'     bg = NA,
#'     x = NULL,
#'     y = NULL,
#'     width = NULL,
#'     height = NULL,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     draw = TRUE,
#'     params = NULL,
#'     ...
#' )
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
#' plotted plot according to the units of the plotgardener page.
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
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{legend} object containing relevant
#' placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Load BED data
#' library(plotgardenerData)
#' data("IMR90_ChIP_CTCF_reads")
#'
#' ## Create page
#' pageCreate(width = 7.5, height = 4, default.units = "inches")
#'
#' ## Plot a pileup plot, coloring elements by strand
#' pileupPlot <- plotRanges(
#'     data = IMR90_ChIP_CTCF_reads, chrom = "chr21",
#'     chromstart = 29072500, chromend = 29075000,
#'     assembly = "hg19",
#'     fill = colorby("strand", palette = 
#'                 colorRampPalette(c("steel blue", "light salmon"))),
#'     x = 0.5, y = 3.5, width = 6.5, height = 3.5,
#'     just = c("left", "bottom"),
#'     default.units = "inches"
#' )
#'
#' ## Add a legend depicting strand colors
#' legendPlot <- plotLegend(
#'     legend = c("- strand", "+ strand"),
#'     fill = c("steel blue", "light salmon"),
#'     border = FALSE,
#'     x = 5, y = 0.5, width = 1.5, height = 0.7,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' annoGenomeLabel(
#'     plot = pileupPlot, x = 0.5, y = 3.5,
#'     just = c("left", "top")
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @export
plotLegend <- function(legend, fill = NULL, pch = NULL, lty = NULL,
                        orientation = "v", title = NULL, fontsize = 10,
                        border = TRUE,
                        bg = NA, x = NULL, y = NULL, width = NULL,
                        height = NULL, just = c("left", "top"),
                        default.units = "inches", draw = TRUE,
                        params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    legInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "legInternal"
    )

    ## Set gp
    legInternal$gp <- setGP(
        gpList = gpar(),
        params = legInternal, ...
    )

    
    
    ## Reset lty
    if (is.null(legInternal$gp$lty)) {
        legInternal$gp$lty <- NULL
    }

    ## Justification
    legInternal$just <- justConversion(just = legInternal$just)
    
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    legendPlot <- structure(list(
        x = legInternal$x, y = legInternal$y,
        width = legInternal$width,
        height = legInternal$height,
        just = legInternal$just, grobs = NULL
    ),
    class = "legend"
    )
    attr(x = legendPlot, which = "plotted") <- legInternal$draw

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(legInternal$legend)) stop("argument \"legend\" is missing, ",
                                            "with no default.", call. = FALSE)

    check_placement(object = legendPlot)
    
    checkColorby(fill = legInternal$fill,
                    colorby = FALSE)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    legendPlot <- defaultUnits(
        object = legendPlot,
        default.units = legInternal$default.units
    )

    textHeight <- heightDetails(textGrob(
        label = "A",
        gp = gpar(fontsize = legInternal$fontsize)
    ))
    textGrobs <- lapply(legInternal$legend, textGrob,
        gp = gpar(fontsize = legInternal$fontsize)
    )
    textWidths <- lapply(textGrobs, widthDetails)

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "legend",
        length(grep(
            pattern = "legend",
            x = currentViewports
        )) + 1
    )

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate
    ## Not translating into page_coordinates
    if (is.null(legendPlot$x) | is.null(legendPlot$y)) {
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

        if (legInternal$draw == TRUE) {
            vp$name <- "legend1"
            grid.newpage()
        }
    } else {

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = legendPlot)
        addViewport(vp_name)

        height <- convertHeight(page_coords$height,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
        width <- convertWidth(page_coords$width,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            just = legInternal$just,
            yscale = c(0, height),
            xscale = c(0, width),
            name = vp_name
        )

        textHeight <- convertHeight(textHeight,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
        textWidths <- lapply(textWidths, convertWidth,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS
    # =========================================================================

    assign("legend_grobs", gTree(vp = vp), envir = pgEnv)

    # =========================================================================
    # GROBS
    # =========================================================================

    ## Border
    if (legInternal$border == TRUE) {
        
        if (is.null(legInternal$gp$border.lty)){
            legInternal$gp$border.lty <- 1
        }
        if (is.null(legInternal$gp$border.linecolor)){
            legInternal$gp$border.linecolor <- "black"
        }
        
        border <- rectGrob(gp = gpar(fill = legInternal$bg,
                                     lty = legInternal$gp$border.lty,
                                     col = legInternal$gp$border.linecolor))
    } else {
        legInternal$gp$fill <- legInternal$bg
        legInternal$gp$col <- NA
        border <- rectGrob(gp = legInternal$gp)
    }

    assign("legend_grobs",
        addGrob(get("legend_grobs", envir = pgEnv),
            child = border
        ),
        envir = pgEnv
    )

    ## Title
    if (!is.null(legInternal$title)) {
        if (legInternal$orientation == "h") {
            remainingSpace <- height - textHeight * 2
            spaceNo <- 3
        } else {
            remainingSpace <- height - textHeight *
                (length(legInternal$legend) + 1)
            spaceNo <- length(legInternal$legend) + 2
        }

        ## Remove lty
        lty <- legInternal$gp$lty
        legInternal$gp$lty <- NULL
        ## Remove cex for text
        cex <- legInternal$gp$cex
        legInternal$gp$cex <- NULL

        legInternal$gp$col <- NA

        titleGrob <- textGrob(
            label = legInternal$title,
            x = unit(0.5, "npc"),
            y = height - remainingSpace / spaceNo,
            just = "top", gp = legInternal$gp,
            default.units = "native"
        )
        assign("legend_grobs",
            addGrob(get("legend_grobs", envir = pgEnv),
                child = titleGrob
            ),
            envir = pgEnv
        )

        legInternal$gp$lty <- lty
        legInternal$gp$cex <- cex
    } else {
        if (legInternal$orientation == "h") {
            remainingSpace <- height - textHeight
            spaceNo <- 2
        } else {
            remainingSpace <- height - textHeight *
                length(legInternal$legend)
            spaceNo <- length(legend) + 1
        }
    }


    ## Colors
    ## Only take the first number of colors for the length of legend
    fillcolors <- legInternal$fill[seq(1, length(legInternal$legend))]
    ## Spacing and label coordinates
    spaceHeight <- remainingSpace / spaceNo

    if (legInternal$orientation == "h") {
        ycoords <- spaceHeight
        remainingWidth <- width - (sum(unlist(textWidths)) +
            textHeight *
                length(legInternal$legend))
        width_spaces <- 2 * length(legInternal$legend) + 1
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
                1, (length(legInternal$legend) - 1)
            )]
        )
        df <- as.data.frame(cbind(
            "factor" = seq(0, (length(legInternal$legend) - 1)),
            "firstSpace" = rep(widthSpace, length(legInternal$legend))
        ))
        df$totalText <- unlist(lapply(df$factor, addTexts,
            textWidths = unlist(textWidths)
        ))

        xcoords <- df$firstSpace + df$factor *
            (2 * widthSpace + textHeight) + df$totalText

        xpch <- xcoords + 0.5 * textHeight
        ypch <- rep(
            spaceHeight + 0.5 * textHeight,
            length(legInternal$legend)
        )

        x0s <- xcoords
        y0s <- spaceHeight + 0.5 * textHeight
        x1s <- x0s + textHeight
        y1s <- spaceHeight + 0.5 * textHeight
    } else {
        ycoords <- seq(
            from = spaceHeight, to = height,
            by = (spaceHeight + textHeight)
        )[seq(1, length(legInternal$legend))]
        xcoords <- textHeight * 2

        xpch <- rep(textHeight * 2, length(legInternal$legend))
        ypch <- ycoords + 0.5 * textHeight

        x0s <- textHeight
        y0s <- ycoords + 0.5 * textHeight
        x1s <- 3 * textHeight
        y1s <- ycoords + 0.5 * textHeight

        fillcolors <- rev(fillcolors)
    }

    ## Symbols
    if (!is.null(legInternal$pch)) {
        ## Only take the first number of symbols for the length of legend
        pchs <- legInternal$pch[seq(1, length(legInternal$legend))]
        if (legInternal$orientation == "v") {
            pchs <- rev(pchs)
        }

        legInternal$gp$col <- fillcolors

        labelSymbols <- pointsGrob(
            x = xpch, y = ypch,
            pch = pchs,
            gp = legInternal$gp,
            default.units = "native"
        )
        assign("legend_grobs",
            addGrob(get("legend_grobs", envir = pgEnv),
                child = labelSymbols
            ),
            envir = pgEnv
        )

        ## Lines
    } else if (!is.null(lty)) {
        ## Only take the number of linetypes for the length of legend
        if (length(lty) < length(legend)){
            ltys <- rep(lty, length(legend))
        } else {
            ltys <- lty
        }
        ltys <- ltys[seq(1, length(legend))]
        lty <- legInternal$gp$lty
        legInternal$gp$lty <- ltys
        legInternal$gp$col <- fillcolors
        labelLines <- segmentsGrob(
            x0 = x0s, y0 = y0s,
            x1 = x1s, y1 = y1s,
            gp = legInternal$gp,
            default.units = "native"
        )
        assign("legend_grobs",
            addGrob(get("legend_grobs", envir = pgEnv),
                child = labelLines
            ),
            envir = pgEnv
        )

        legInternal$gp$lty <- lty
    } else {
        labelColors <- rectGrob(
            x = xcoords, y = ycoords,
            width = textHeight, height = textHeight,
            just = c("left", "bottom"),
            gp = gpar(fill = fillcolors, col = NA),
            default.units = "native"
        )
        assign("legend_grobs",
            addGrob(get("legend_grobs", envir = pgEnv),
                child = labelColors
            ),
            envir = pgEnv
        )
    }

    legInternal$gp$fontsize <- legInternal$fontsize

    if (is.null(legInternal$gp$fontcolor)) {
        legInternal$gp$fontcolor <- "black"
        legInternal$gp$col <- "black"
    } else {
        legInternal$gp$col <- legInternal$gp$fontcolor
    }
    legInternal$gp$lty <- NULL
    legInternal$gp$cex <- NULL

    ## Text labels
    if (legInternal$orientation == "h") {
        labelText <- textGrob(
            label = legInternal$legend,
            x = xcoords + textHeight + widthSpace,
            y = ycoords + 0.5 * textHeight,
            just = c("left", "center"),
            gp = legInternal$gp, default.units = "native"
        )
    } else {
        labelText <- textGrob(
            label = rev(legInternal$legend),
            x = xcoords + 2 * textHeight,
            y = ycoords + 0.5 * textHeight,
            just = c("left", "center"),
            gp = legInternal$gp, default.units = "native"
        )
    }

    assign("legend_grobs",
        addGrob(get("legend_grobs", envir = pgEnv),
            child = labelText
        ),
        envir = pgEnv
    )


    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (legInternal$draw == TRUE) {
        grid.draw(get("legend_grobs", envir = pgEnv))
    }
    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    legendPlot$grobs <- get("legend_grobs", envir = pgEnv)

    # ========================================================================
    # RETURN OBJECT
    # ========================================================================

    message("legend[", vp$name, "]")
    invisible(legendPlot)
}
