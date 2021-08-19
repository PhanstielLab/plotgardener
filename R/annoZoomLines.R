#' Annotates zoom lines for a specified genomic region of a plot
#' 
#' @usage annoZoomLines(
#'     plot,
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     y0,
#'     x1 = NULL,
#'     y1,
#'     extend = 0,
#'     default.units = "inches",
#'     linecolor = "grey",
#'     lty = 2,
#'     params = NULL,
#'     ...
#' )
#'
#' @param plot Input plot to annotate genomic region zoom lines from.
#' @param chrom Chromosome of region to draw zoom lines from, as a string.
#' @param chromstart Integer start position on chromosome to draw
#' zoom lines from.
#' @param chromend Integer end position on chromosome to draw
#' zoom lines from.
#' @param y0 A numeric vector or unit object indicating the starting
#' y-values of the zoom line segments. If two values are given,
#' the first value will correspond to the left zoom line and the
#' second value will correspond to the right zoom line.
#' @param x1 A numeric vector or unit object indicating the stopping
#' x-values of the zoom line segments. If two values are given,
#' the first value will correspond to the left zoom line and the
#' second value will correspond to the right zoom line. If NULL,
#' straight lines from zoomed genomic region will be drawn.
#' @param y1 A numeric vector or unit object indicating the stopping
#' y-values of the zoom line segments. If two values are given,
#' the first value will correspond to the left zoom line and the
#' second value will correspond to the right zoom line.
#' @param extend A numeric vector or unit object indicating the length
#' to extend straight lines from each end
#' of the zoom line segments. If two values are given, the first value
#' will correspond to the top extension length
#' and the second value will correspond to the bottom extension length.
#' Default value is \code{extend = 0}.
#' @param default.units A string indicating the default units to use
#' if \code{y0}, \code{x1}, \code{y1}, or \code{extend} are only given
#' as numerics or numeric vectors.
#' Default value is \code{default.units = "inches"}.
#' @param linecolor A character value specifying zoom line color.
#' Default value is \code{linecolor = "grey"}.
#' @param lty A numeric specifying zoom line type.
#' Default value is \code{lty = 2}.
#' @param params An optional \link[plotgardener]{pgParams}
#' object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{zoom} object containing
#' relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' pageCreate(width = 7.5, height = 4.75, default.units = "inches")
#'
#' ## Plot and place a Manhattan plot
#' library(plotgardenerData)
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' data("hg19_insulin_GWAS")
#' manhattanPlot <- plotManhattan(
#'     data = hg19_insulin_GWAS, assembly = "hg19",
#'     fill = c("grey", "#37a7db"),
#'     sigLine = FALSE,
#'     col = "grey", lty = 2, range = c(0, 14),
#'     x = 0.5, y = 0, width = 6.5, height = 2,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#' annoYaxis(
#'     plot = manhattanPlot, at = c(0, 2, 4, 6, 8, 10, 12, 14),
#'     axisLine = TRUE, fontsize = 8
#' )
#'
#' ## Annotate zoom lines for a region on chromsome 21
#' zoomRegion <- pgParams(
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19"
#' )
#' annoZoomLines(
#'     plot = manhattanPlot, params = zoomRegion,
#'     y0 = 2, x1 = c(0.5, 7), y1 = 2.5, extend = c(0, 1.1),
#'     default.units = "inches",
#'     lty = 3
#' )
#'
#' ## Annotate highlight region for zoom region
#' annoHighlight(
#'     plot = manhattanPlot, params = zoomRegion,
#'     y = 2, height = 2, just = c("left", "bottom"),
#'     default.units = "inches",
#'     fill = "red", alpha = 0.8
#' )
#'
#' ## Plot Manhattan plot data and signal track under zoom lines
#' manhattanPlotZoom <- plotManhattan(
#'     data = hg19_insulin_GWAS,
#'     fill = "grey",
#'     sigLine = FALSE,
#'     baseline = TRUE,
#'     params = zoomRegion, range = c(0, 14),
#'     x = 0.5, y = 2.6,
#'     width = 6.5, height = 1
#' )
#' data("IMR90_ChIP_H3K27ac_signal")
#' signalPlot <- plotSignal(
#'     data = IMR90_ChIP_H3K27ac_signal, params = zoomRegion,
#'     range = c(0, 45),
#'     x = 0.5, y = "b0.1",
#'     width = 6.5, height = 0.65,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Plot genome label
#' plotGenomeLabel(
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     x = 0.5, y = 4.4, length = 6.5,
#'     default.units = "inches"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @export
annoZoomLines <- function(plot, chrom, chromstart = NULL, chromend = NULL,
                            y0, x1 = NULL, y1, extend = 0,
                            default.units = "inches", linecolor = "grey",
                            lty = 2, params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for annoZoomLines
    errorcheck_annoZoomLines <- function(object) {
        if (!object$chrom %in% object$plot$chrom) {
            stop(object$chrom, "not found in input plot. ",
                "Cannot annotate zoom lines.", call. = FALSE)
        }

        if (length(object$chrom) > 1) {
            stop("Cannot annotate zoom lines for multiple chromosome regions ",
                "in one function call.", call. = FALSE)
        }
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    zoomInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "zoomInternal"
    )

    ## Set gp
    zoomInternal$gp <- gpar(
        col = zoomInternal$linecolor,
        lty = zoomInternal$lty
    )
    zoomInternal$gp <- setGP(
        gpList = zoomInternal$gp,
        params = zoomInternal, ...
    )

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    zoom <- structure(list(
        chrom = zoomInternal$chrom,
        chromstart = zoomInternal$chromstart,
        chromend = zoomInternal$chromend,
        assembly = zoomInternal$plot$assembly,
        x0 = NULL, y0 = zoomInternal$y0,
        x1 = zoomInternal$x1, y1 = zoomInternal$y1,
        extend = zoomInternal$extend, grobs = NULL
    ),
    class = "zoom"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_page(error = "Cannot annotate zoom lines without a `plotgardener` page.")
    if (is.null(zoomInternal$plot)) stop("argument \"plot\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(zoom$chrom)) {
        stop("argument \"chrom\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(zoom$y0)) {
        stop("argument \"y0\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(zoom$y1)) {
        stop("argument \"y1\" is missing, with no default.",
            call. = FALSE
        )
    }
    errorcheck_annoZoomLines(object = zoomInternal)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    ## y0
    zoom$y0 <- misc_defaultUnits(value = zoom$y0,
                                    name = "y0",
                                    default.units = 
                                        zoomInternal$default.units,
                                    funName = "annoZoomLines",
                                    yBelow = FALSE)
    
    if (length(zoom$y0) == 1) {
        zoom$y0 <- rep(zoom$y0, 2)
    }

    ## y1
    zoom$y1 <- misc_defaultUnits(value = zoom$y1,
                                    name = "y1",
                                    default.units = 
                                        zoomInternal$default.units,
                                    funName = "annoZoomLines",
                                    yBelow = FALSE)

    if (length(zoom$y1) == 1) {
        zoom$y1 <- rep(zoom$y1, 2)
    }

    ## extend
    zoom$extend <- misc_defaultUnits(value = zoom$extend,
                                    name = "extend",
                                    default.units = 
                                        zoomInternal$default.units)

    if (length(zoom$extend) == 1) {
        zoom$extend <- rep(zoom$extend, 2)
    }

    # =========================================================================
    # WHOLE CHROM GENOMIC SCALE
    # =========================================================================
    
    scaleChecks <- genomicScale(object = zoom,
                                objectInternal = zoomInternal,
                                plotType = "zoom lines")
    zoom <- scaleChecks[[1]]

    # =========================================================================
    # PARSE GENOMIC REGION FOR X0
    # =========================================================================
    page_units <- get("page_units", envir = pgEnv)
    
    ## Get appropriate plot viewport
    plotVP <- getAnnoViewport(plot = zoomInternal$plot)
    
    ## Convert plot viewport to bottom left to get left position on entire page
    plotVP_bottomLeft <- vp_bottomLeft(viewport = plotVP)

    ## Convert plot genomic coordinates to position on page
    seekViewport(plotVP$name)

    if (!is.null(zoom$chromstart) & !is.null(zoom$chromend)) {
        if (is(zoomInternal$plot, "manhattan")) {

            ## Multiple chromosome manhattan plot
            if (length(zoomInternal$plot$chrom) > 1) {
                convertedCoords <- convertManhattan(
                    object = zoom,
                    manhattanPlot = zoomInternal$plot
                )
                x0_1 <- convertedCoords[[1]]
                x0_2 <- convertedCoords[[2]]
            } else {
                x0_1 <- convertX(unit(zoom$chromstart, "native"),
                    unitTo = page_units, valueOnly = TRUE
                )
                x0_2 <- convertX(unit(zoom$chromend, "native"),
                    unitTo = page_units, valueOnly = TRUE
                )
            }
        } else {
            x0_1 <- convertX(unit(zoom$chromstart, "native"),
                unitTo = page_units, valueOnly = TRUE
            )
            x0_2 <- convertX(unit(zoom$chromend, "native"),
                unitTo = page_units, valueOnly = TRUE
            )
        }

        ## Add additional page units to x0_1 and x0_2
        x0_1 <- as.numeric(plotVP_bottomLeft[[1]]) + x0_1
        x0_2 <- as.numeric(plotVP_bottomLeft[[1]]) + x0_2

        zoom$x0 <- unit(c(x0_1, x0_2), page_units)
    }

    seekViewport("page")

    # =========================================================================
    # PARSE X1
    # =========================================================================

    if (is.null(zoom$x1)) {
        zoom$x1 <- zoom$x0
    } else {
        if (!"unit" %in% class(zoom$x1)) {
            if (!is.numeric(zoom$x1)) {
                stop("x1-coordinate is neither a unit object or a numeric ",
                    "value. Cannot annotate zoom lines.", call. = FALSE)
            }

            if (is.null(zoomInternal$default.units)) {
                stop("x1-coordinate detected as numeric.\'default.units\' ",
                    "must be specified.", call. = FALSE)
            }

            zoom$x1 <- unit(zoom$x1, zoomInternal$default.units)
        }

        if (length(zoom$x1) == 1) {
            zoom$x1 <- rep(zoom$x1, 2)
        }
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS
    # =========================================================================
    name <- paste0(
        "zoom",
        length(grep(
            pattern = "zoom",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    assign("zoom_grobs", gTree(name = name), envir = pgEnv)

    if (!is.null(zoom$chromstart) & !is.null(zoom$chromend)) {
        # =====================================================================
        # ZOOM SEGMENTS
        # =====================================================================
        page_height <- get("page_height", envir = pgEnv)
        zoomSegment_left <- segmentsGrob(
            x0 = zoom$x0[1],
            y0 = unit(page_height, page_units) - zoom$y0[1],
            x1 = zoom$x1[1],
            y1 = unit(page_height, page_units) - zoom$y1[1],
            gp = zoomInternal$gp
        )
        zoomSegment_right <- segmentsGrob(
            x0 = zoom$x0[2],
            y0 = unit(page_height, page_units) - zoom$y0[2],
            x1 = zoom$x1[2],
            y1 = unit(page_height, page_units) - zoom$y1[2],
            gp = zoomInternal$gp
        )

        assign("zoom_grobs",
            addGrob(
                gTree = get("zoom_grobs", envir = pgEnv),
                child = zoomSegment_left
            ),
            envir = pgEnv
        )
        assign("zoom_grobs",
            addGrob(
                gTree = get("zoom_grobs", envir = pgEnv),
                child = zoomSegment_right
            ),
            envir = pgEnv
        )

        # =====================================================================
        # EXTEND SEGMENTS
        # =====================================================================

        topy1_left <- (unit(page_height, page_units) - zoom$y0[1]) +
            zoom$extend[1]
        topy1_right <- (unit(page_height, page_units) - zoom$y0[2]) +
            zoom$extend[1]
        bottomy1_left <- (unit(page_height, page_units) - zoom$y1[1]) -
            zoom$extend[2]
        bottomy1_right <- (unit(page_height, page_units) - zoom$y1[2]) -
            zoom$extend[2]


        extend_topLeft <- segmentsGrob(
            x0 = zoom$x0[1],
            y0 = unit(page_height, page_units) - zoom$y0[1],
            x1 = zoom$x0[1],
            y1 = topy1_left,
            gp = zoomInternal$gp
        )
        extend_topRight <- segmentsGrob(
            x0 = zoom$x0[2],
            y0 = unit(page_height, page_units) - zoom$y0[2],
            x1 = zoom$x0[2],
            y1 = topy1_right,
            gp = zoomInternal$gp
        )
        extend_bottomLeft <- segmentsGrob(
            x0 = zoom$x1[1],
            y0 = unit(page_height, page_units) - zoom$y1[1],
            x1 = zoom$x1[1],
            y1 = bottomy1_left,
            gp = zoomInternal$gp
        )
        extend_bottomRight <- segmentsGrob(
            x0 = zoom$x1[2],
            y0 = unit(page_height, page_units) - zoom$y1[2],
            x1 = zoom$x1[2],
            y1 = bottomy1_right,
            gp = zoomInternal$gp
        )

        assign("zoom_grobs",
            addGrob(
                gTree = get("zoom_grobs", envir = pgEnv),
                child = extend_topLeft
            ),
            envir = pgEnv
        )
        assign("zoom_grobs",
            addGrob(
                gTree = get("zoom_grobs", envir = pgEnv),
                child = extend_topRight
            ),
            envir = pgEnv
        )
        assign("zoom_grobs",
            addGrob(
                gTree = get("zoom_grobs", envir = pgEnv),
                child = extend_bottomLeft
            ),
            envir = pgEnv
        )
        assign("zoom_grobs",
            addGrob(
                gTree = get("zoom_grobs", envir = pgEnv),
                child = extend_bottomRight
            ),
            envir = pgEnv
        )
    }

    # =========================================================================
    # ADD GTREE TO OBJECT
    # =========================================================================

    zoom$grobs <- get("zoom_grobs", envir = pgEnv)
    grid.draw(zoom$grobs)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("zoom[", zoom$grobs$name, "]")
    invisible(zoom)
}
