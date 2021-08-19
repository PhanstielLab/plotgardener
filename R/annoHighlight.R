#' Annotates a highlight box around a specified genomic region of a
#' plot
#' 
#' @usage annoHighlight(
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
#' @param plot Input plot on which to annotate genomic region.
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
#' bottom of the most recently plotted plot according to the
#' units of the \code{plotgardener} page.
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
#' @param params An optional \link[plotgardener]{params} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{highlight} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' pageCreate(width = 7.5, height = 1.5, default.units = "inches")
#'
#' ## Plot and place a signal plot
#' library(plotgardenerData)
#' data("IMR90_ChIP_H3K27ac_signal")
#' region <- params(
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     range = c(0, 45)
#' )
#' signalPlot <- plotSignal(
#'     data = IMR90_ChIP_H3K27ac_signal, params = region,
#'     x = 0.5, y = 0.25, width = 6.5, height = 0.65,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Highlight genomic region on signal plot
#' annoHighlight(
#'     plot = signalPlot,
#'     chrom = "chr21",
#'     chromstart = 29000000, chromend = 29125000,
#'     y = 0.25, height = 1, just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Plot text label
#' plotText(
#'     label = "region of interest", fontsize = 8, fontcolor = "black",
#'     x = 3.5, y = 0.2, just = "bottom", default.units = "inches"
#' )
#'
#' ## Plot genome label
#' plotGenomeLabel(
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     x = 0.5, y = 1.3, length = 6.5, default.units = "inches"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @export
annoHighlight <- function(plot, chrom, chromstart = NULL, chromend = NULL,
                            fill = "grey", linecolor = NA, alpha = 0.4,
                            y, height, just = c("left", "top"),
                            default.units = "inches", params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for annoHighlight
    errorcheck_annoHighlight <- function(object, fill) {
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
        
        checkColorby(fill = fill,
                    colorby = FALSE)
    }

    ## Define a function that resets y to a top-based justification
    y_Pagetop <- function(y, height, just) {
        page_height <- get("page_height", envir = pgEnv)
        page_units <- get("page_units", envir = pgEnv)
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

    highlightInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "highlightInternal"
    )

    ## Set gp
    highlightInternal$gp <- gpar(
        fill = highlightInternal$fill,
        col = highlightInternal$linecolor,
        alpha = highlightInternal$alpha
    )
    highlightInternal$gp <- setGP(
        gpList = highlightInternal$gp,
        params = highlightInternal, ...
    )

    ## Justification
    highlightInternal$just <- 
        justConversion(just = highlightInternal$just)
    
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    highlight <- structure(list(
        chrom = highlightInternal$chrom,
        chromstart = highlightInternal$chromstart,
        chromend = highlightInternal$chromend,
        assembly = highlightInternal$plot$assembly,
        x = NULL, y = highlightInternal$y,
        width = NULL,
        height = highlightInternal$height,
        just = highlightInternal$just,
        grobs = NULL
    ),
    class = "highlight"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_page(error = "Cannot add highlight annotation
                without a `plotgardener` page.")
    if (is.null(highlightInternal$plot)) {
        stop("argument \"plot\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(highlightInternal$chrom)) {
        stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(highlight$y)) {
        stop("argument \"y\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(highlight$height)) {
        stop("argument \"height\" is missing, with no default.",
            call. = FALSE
        )
    }
    errorcheck_annoHighlight(object = highlightInternal,
                                fill = highlightInternal$gp$fill)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    page_units <- get("page_units", envir = pgEnv)

    highlight$y <- misc_defaultUnits(value = highlight$y,
                                        name = "y",
                                        default.units = 
                                            highlightInternal$default.units)
    highlight$height <- misc_defaultUnits(value = highlight$height,
                                        name = "height",
                                        default.units = 
                                            highlightInternal$default.units)

    if (is(highlightInternal$plot, "genes")) {
        plotVP <- highlightInternal$plot$grobs$children$background$vp
    } else if (is(highlightInternal$plot, "hicTriangle")  |
        is(highlightInternal$plot, "hicRectangle")) {
        plotVP <- highlightInternal$plot$outsideVP
    } else {
        plotVP <- highlightInternal$plot$grobs$vp
    }

    # =========================================================================
    # WHOLE CHROM GENOMIC SCALE
    # =========================================================================

    scaleChecks <- genomicScale(object = highlight,
                                objectInternal = highlightInternal,
                                plotType = "highlight")
    highlight <- scaleChecks[[1]]
    
    # =========================================================================
    # DETERMINE X AND WIDTH BASED ON GENOMIC REGION
    # =========================================================================

    ## Get plot viewport
    if (is(highlightInternal$plot, "genes")) {
        plotVP <- highlightInternal$plot$grobs$children$background$vp
    } else if (is(highlightInternal$plot, "hicTriangle")  |
        is(highlightInternal$plot, "hicRectangle")) {
        plotVP <- highlightInternal$plot$outsideVP
    } else {
        plotVP <- highlightInternal$plot$grobs$vp
    }

    ## Convert plot viewport to bottom left to get left position on entire page
    plotVP_bottomLeft <- vp_bottomLeft(viewport = plotVP)

    ## Convert plot genomic coordinates to position on page
    seekViewport(plotVP$name)

    if (!is.null(highlight$chromstart) & !is.null(highlight$chromend)) {
        if (is(highlightInternal$plot, "manhattan")) {

            ## Multiple chromosome manhattan plot
            if (length(highlightInternal$plot$chrom) > 1) {
                convertedCoords <- convertManhattan(
                    object = highlight,
                    manhattanPlot = highlightInternal$plot
                )
                start <- convertedCoords[[1]]
                end <- convertedCoords[[2]]
            } else {
                start <- convertX(unit(highlight$chromstart, "native"),
                    unitTo = page_units, valueOnly = TRUE
                )
                end <- convertX(unit(highlight$chromend, "native"),
                    unitTo = page_units, valueOnly = TRUE
                )
            }
        } else {
            start <- convertX(unit(highlight$chromstart, "native"),
                unitTo = page_units, valueOnly = TRUE
            )
            end <- convertX(unit(highlight$chromend, "native"),
                unitTo = page_units, valueOnly = TRUE
            )
        }

        width <- end - start
        ## Add additional page units to start
        start <- as.numeric(plotVP_bottomLeft[[1]]) + start
    }

    seekViewport("page")

    # =========================================================================
    # USE JUSTIFICATION TO DETERMINE PAGE-ADJUSTED TOP Y-COORD
    # =========================================================================

    top_y <- y_Pagetop(
        y = highlight$y,
        height = highlight$height,
        just = highlight$just
    )

    # =========================================================================
    # HIGHLIGHT GROB
    # =========================================================================

    name <- paste0(
        "highlight",
        length(grep(
            pattern = "highlight",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )

    if (!is.null(highlight$chromstart) & !is.null(highlight$chromend)) {
        
        if (width == 0){
            
            highlightInternal$gp$col <- highlightInternal$gp$fill
            
            highlightGrob <- grid.segments(
                x0 = unit(start, page_units), 
                x1 = unit(start, page_units),
                y0 = top_y,
                y1 = top_y - highlight$height,
                gp = highlightInternal$gp, name = name
            )
        } else {
            highlightGrob <- grid.rect(
                x = unit(start, page_units), y = top_y,
                width = unit(width, page_units),
                height = highlight$height,
                just = c("left", "top"),
                gp = highlightInternal$gp, name = name
            )
        }
        highlight$grobs <- highlightGrob
    }

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("highlight[", name, "]")
    invisible(highlight)
}
