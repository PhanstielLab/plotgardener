#' Annotates a highlight box around a specified genomic region of a
#' BentoBox plot
#' 
#' @usage bbAnnoHighlight(
#'     plot,
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     fill = "grey",
#'     linecolor = NA,
#'     alpha = 0.4,
#'     y,
#'     height,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     params = NULL,
#'     ...
#' )
#'
#' @param plot Input BentoBox plot on which to annotate genomic region.
#' @param chrom Chromosome of region to be highlighted, as a string.
#' @param chromstart Integer start position on chromosome to be highlighted.
#' @param chromend Integer end position on chromosome to be highlighted.
#' @param fill A character value specifying highlight box fill color.
#' Default value is \code{fill = "grey"}.
#' @param linecolor A character value specifying highlight box line color.
#' Default value is \code{linecolor = NA}.
#' @param alpha Numeric value specifying color transparency.
#' Default value is \code{alpha = 0.4}.
#' @param y A numeric, unit object, or character containing a "b" combined
#' with a numeric value specifying square highlight box y-location.
#' The character value will place the highlight box y relative to the
#' bottom of the most recently plotted BentoBox plot according to the
#' units of the BentoBox page.
#' @param height A numeric or unit object specifying highlight box height.
#' @param just Justification of highlight box relative to its (x, y)
#' location. If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if
#' \code{y} or \code{height} are only given as numerics or numeric vectors.
#' Default value is \code{default.units = "inches"}.
#' @param params An optional \link[BentoBox]{bbParams} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_highlight} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' bbPageCreate(width = 7.5, height = 1.5, default.units = "inches")
#'
#' ## Plot and place a signal plot
#' library(BentoBoxData)
#' data("IMR90_ChIP_H3K27ac_signal")
#' region <- bbParams(
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     range = c(0, 45)
#' )
#' signalPlot <- bbPlotSignal(
#'     data = IMR90_ChIP_H3K27ac_signal, params = region,
#'     x = 0.5, y = 0.25, width = 6.5, height = 0.65,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Highlight genomic region on signal plot
#' bbAnnoHighlight(
#'     plot = signalPlot,
#'     chrom = "chr21",
#'     chromstart = 29000000, chromend = 29125000,
#'     y = 0.25, height = 1, just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Plot text label
#' bbPlotText(
#'     label = "region of interest", fontsize = 8, fontcolor = "black",
#'     x = 3.5, y = 0.2, just = "bottom", default.units = "inches"
#' )
#'
#' ## Plot genome label
#' bbPlotGenomeLabel(
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     x = 0.5, y = 1.3, length = 6.5, default.units = "inches"
#' )
#'
#' ## Hide page guides
#' bbPageGuideHide()
#' @export
bbAnnoHighlight <- function(plot, chrom, chromstart = NULL, chromend = NULL,
                            fill = "grey", linecolor = NA, alpha = 0.4,
                            y, height, just = c("left", "top"),
                            default.units = "inches", params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for bbAnnoHighlight
    bb_errorcheckAnnoHighlight <- function(object, fill) {
        if (!object$chrom %in% object$plot$chrom) {
            stop(object$chrom,
                "not found in input plot. Cannot annotate highlight.",
                call. = FALSE
            )
        }

        if (length(object$chrom) > 1) {
            stop("Cannot highlight multiple chromosome ",
            "regions in one function call.", call. = FALSE)
        }
        
        bb_checkColorby(fill = fill,
                        colorby = FALSE)
    }

    ## Define a function that resets y to a top-based justification
    y_Pagetop <- function(y, height, just) {
        page_height <- get("page_height", envir = bbEnv)
        page_units <- get("page_units", envir = bbEnv)
        y <- convertY(unit(page_height, units = page_units) - 
            convertY(y, unitTo = page_units), unitTo = page_units)
        height <- convertHeight(height, unitTo = page_units)
        if (identical(just, c(0, 0.5))) {
            ## left
            topy <- y + 0.5 * height
        } else if (identical(just, c(1, 0.5))) {
            ## right
            topy <- y + 0.5 * height
        } else if (identical(just, c(0.5, 0))) {
            ## bottom
            topy <- y + height
        } else if (identical(just, c(0.5, 1))) {
            ## top
            topy <- y
        } else if (identical(just, c(0, 1))) {
            ## left top
            topy <- y
        } else if (identical(just, c(1, 1))) {
            ## right top
            topy <- y
        } else if (identical(just, c(0, 0))) {
            ## left bottom
            topy <- y + height
        } else if (identical(just, c(1, 0))) {
            ## right bottom
            topy <- y + height
        } else {
            ## center
            topy <- y + 0.5 * height
        }
        
        return(topy)
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_highlightInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_highlightInternal"
    )

    ## Set gp
    bb_highlightInternal$gp <- gpar(
        fill = bb_highlightInternal$fill,
        col = bb_highlightInternal$linecolor,
        alpha = bb_highlightInternal$alpha
    )
    bb_highlightInternal$gp <- setGP(
        gpList = bb_highlightInternal$gp,
        params = bb_highlightInternal, ...
    )

    ## Justification
    bb_highlightInternal$just <- 
        bb_justConversion(just = bb_highlightInternal$just)
    
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    bb_highlight <- structure(list(
        chrom = bb_highlightInternal$chrom,
        chromstart = bb_highlightInternal$chromstart,
        chromend = bb_highlightInternal$chromend,
        assembly = bb_highlightInternal$plot$assembly,
        x = NULL, y = bb_highlightInternal$y,
        width = NULL,
        height = bb_highlightInternal$height,
        just = bb_highlightInternal$just,
        grobs = NULL
    ),
    class = "bb_highlight"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_bbpage(error = "Cannot add highlight annotation
                without a BentoBox page.")
    if (is.null(bb_highlightInternal$plot)) {
        stop("argument \"plot\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(bb_highlightInternal$chrom)) {
        stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(bb_highlight$y)) {
        stop("argument \"y\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(bb_highlight$height)) {
        stop("argument \"height\" is missing, with no default.",
            call. = FALSE
        )
    }
    bb_errorcheckAnnoHighlight(object = bb_highlightInternal,
                                fill = bb_highlightInternal$gp$fill)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    page_units <- get("page_units", envir = bbEnv)

    bb_highlight$y <- misc_defaultUnits(value = bb_highlight$y,
                                        name = "y",
                                        default.units = 
                                            bb_highlightInternal$default.units)
    bb_highlight$height <- misc_defaultUnits(value = bb_highlight$height,
                                        name = "height",
                                        default.units = 
                                            bb_highlightInternal$default.units)

    if (is(bb_highlightInternal$plot, "bb_genes")) {
        plotVP <- bb_highlightInternal$plot$grobs$children$background$vp
    } else if (is(bb_highlightInternal$plot, "bb_hicTriangle")  |
        is(bb_highlightInternal$plot, "bb_hicRectangle")) {
        plotVP <- bb_highlightInternal$plot$outsideVP
    } else {
        plotVP <- bb_highlightInternal$plot$grobs$vp
    }

    # =========================================================================
    # WHOLE CHROM GENOMIC SCALE
    # =========================================================================

    scaleChecks <- genomicScale(object = bb_highlight,
                                objectInternal = bb_highlightInternal,
                                plotType = "highlight")
    bb_highlight <- scaleChecks[[1]]
    
    # =========================================================================
    # DETERMINE X AND WIDTH BASED ON GENOMIC REGION
    # =========================================================================

    ## Get plot viewport
    if (is(bb_highlightInternal$plot, "bb_genes")) {
        plotVP <- bb_highlightInternal$plot$grobs$children$background$vp
    } else if (is(bb_highlightInternal$plot, "bb_hicTriangle")  |
        is(bb_highlightInternal$plot, "bb_hicRectangle")) {
        plotVP <- bb_highlightInternal$plot$outsideVP
    } else {
        plotVP <- bb_highlightInternal$plot$grobs$vp
    }

    ## Convert plot viewport to bottom left to get left position on entire page
    plotVP_bottomLeft <- vp_bottomLeft(viewport = plotVP)

    ## Convert plot genomic coordinates to position on page
    seekViewport(plotVP$name)

    if (!is.null(bb_highlight$chromstart) & !is.null(bb_highlight$chromend)) {
        if (is(bb_highlightInternal$plot, "bb_manhattan")) {

            ## Multiple chromosome manhattan plot
            if (length(bb_highlightInternal$plot$chrom) > 1) {
                convertedCoords <- convertManhattan(
                    object = bb_highlight,
                    manhattanPlot = bb_highlightInternal$plot
                )
                start <- convertedCoords[[1]]
                end <- convertedCoords[[2]]
            } else {
                start <- convertX(unit(bb_highlight$chromstart, "native"),
                    unitTo = page_units, valueOnly = TRUE
                )
                end <- convertX(unit(bb_highlight$chromend, "native"),
                    unitTo = page_units, valueOnly = TRUE
                )
            }
        } else {
            start <- convertX(unit(bb_highlight$chromstart, "native"),
                unitTo = page_units, valueOnly = TRUE
            )
            end <- convertX(unit(bb_highlight$chromend, "native"),
                unitTo = page_units, valueOnly = TRUE
            )
        }

        width <- end - start
        ## Add additional page units to start
        start <- as.numeric(plotVP_bottomLeft[[1]]) + start
    }

    seekViewport("bb_page")

    # =========================================================================
    # USE JUSTIFICATION TO DETERMINE PAGE-ADJUSTED TOP Y-COORD
    # =========================================================================

    top_y <- y_Pagetop(
        y = bb_highlight$y,
        height = bb_highlight$height,
        just = bb_highlight$just
    )
    print(bb_highlight$y)
    print(top_y)
    # =========================================================================
    # HIGHLIGHT GROB
    # =========================================================================

    name <- paste0(
        "bb_highlight",
        length(grep(
            pattern = "bb_highlight",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )

    if (!is.null(bb_highlight$chromstart) & !is.null(bb_highlight$chromend)) {
        
        if (width == 0){
            
            bb_highlightInternal$gp$col <- bb_highlightInternal$gp$fill
            
            highlightGrob <- grid.segments(
                x0 = unit(start, page_units), 
                x1 = unit(start, page_units),
                y0 = top_y,
                y1 = top_y - bb_highlight$height,
                gp = bb_highlightInternal$gp, name = name
            )
        } else {
            highlightGrob <- grid.rect(
                x = unit(start, page_units), y = top_y,
                width = unit(width, page_units),
                height = bb_highlight$height,
                just = c("left", "top"),
                gp = bb_highlightInternal$gp, name = name
            )
        }
        bb_highlight$grobs <- highlightGrob
    }

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_highlight[", name, "]")
    invisible(bb_highlight)
}
