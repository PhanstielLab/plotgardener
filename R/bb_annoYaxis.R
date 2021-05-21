#' Add a y-axis to a plot
#'
#' @param plot Plot object to annotate with y-axis.
#' @param at A numeric vector of y-value locations for tick marks.
#' @param label A logical value indicating whether to draw the labels
#' on the tick marks, or an expression or character vector which specify
#' the labels to use.
#' If not logical, must be the same length as the \code{at} argument.
#' @param main A logical value indicating whether to draw the y-axis at
#' the left of the plot. Default value is \code{main = TRUE}. Options are:
#' \itemize{
#' \item{\code{TRUE}: }{y-axis is drawn at the left of the plot.}
#' \item{\code{FALSE}: }{y-axis is drawn at the right of the plot.}
#' }
#' @param scipen An integer indicating the penalty to be applied when
#' deciding to print numeric values in fixed or exponential notation.
#' Default value is \code{scipen = 999}.
#' @param axisLine A logical value indicating whether to show the axis line.
#' Default value is \code{axisLine = FALSE}.
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_yaxis} object containing
#' relevant \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data
#' data("bb_imrHicData")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 4, height = 3.5, default.units = "inches")
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- bb_plotHicSquare(
#'     data = bb_imrHicData, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     x = 1, y = 0.5, width = 2.5, height = 2.5,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Add standard y-axis to Hi-C plot
#' bb_annoYaxis(
#'     plot = hicPlot, at = c(28000000, 29000000, 30300000),
#'     fontsize = 10
#' )
#'
#' ## Annotate genome label on x-axis
#' bb_annoGenomeLabel(plot = hicPlot, x = 1, y = 3.03)
#'
#' ## Annotate heatmap legend
#' bb_annoHeatmapLegend(
#'     plot = hicPlot,
#'     x = 3.6, y = 0.5, width = 0.12, height = 1.2
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @export
bb_annoYaxis <- function(plot, at = NULL, label = TRUE, main = TRUE,
                        scipen = 999, axisLine = FALSE, params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    ## Check which defaults are not overwritten and set to NULL
    if (missing(label)) label <- NULL
    if (missing(main)) main <- NULL
    if (missing(scipen)) scipen <- NULL
    if (missing(axisLine)) axisLine <- NULL

    ## Check if plot argument is missing (could be in object)
    if (!hasArg(plot)) plot <- NULL

    ## Compile all parameters into an internal object
    bb_yInternal <- structure(list(
        plot = plot, at = at, label = label,
        main = main, gp = gpar(), scipen = scipen,
        axisLine = axisLine
    ),
    class = "bb_yInternal"
    )

    bb_yInternal <- parseParams(
        bb_params = params,
        object_params = bb_yInternal
    )

    ## For any defaults that are still NULL, set back to default
    if (is.null(bb_yInternal$label)) bb_yInternal$label <- TRUE
    if (is.null(bb_yInternal$main)) bb_yInternal$main <- TRUE
    if (is.null(bb_yInternal$scipen)) bb_yInternal$scipen <- 999
    if (is.null(bb_yInternal$axisLine)) bb_yInternal$axisLine <- FALSE

    ## Set gp
    bb_yInternal$gp <- setGP(
        gpList = bb_yInternal$gp,
        params = bb_yInternal, ...
    )
    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(bb_yInternal$plot)) stop("argument \"plot\" is missing,
                                        with no default.", call. = FALSE)
    check_bbpage(error = "Cannot add an y-axis without a BentoBox page.")

    # =========================================================================
    # SAVE USER'S SCIPEN AND SET OURS TEMPORARILY
    # =========================================================================

    oo <- options(scipen = bb_yInternal$scipen)
    on.exit(options(oo))

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    yAxis <- structure(list(grobs = NULL), class = "bb_yaxis")

    # =========================================================================
    # CREATE GROB WITHOUT DRAWING
    # =========================================================================

    yGrob <- yaxisGrob(
        at = bb_yInternal$at, label = bb_yInternal$label,
        main = bb_yInternal$main, gp = bb_yInternal$gp,
        vp = bb_yInternal$plot$grobs$vp
    )

    # =========================================================================
    # GET CENTER OF INPUT PLOT VIEWPORT BASED ON INPUT PLOT TYPE AND JUST
    # =========================================================================

    if (class(bb_yInternal$plot) == "bb_genes") {
        plotVP <- bb_yInternal$plot$grobs$children$background$vp
    } else if (class(bb_yInternal$plot) == "bb_hicTriangle" |
        class(bb_yInternal$plot) == "bb_hicRectangle") {
        plotVP <- bb_yInternal$plot$outsideVP
    } else {
        plotVP <- bb_yInternal$plot$grobs$vp
    }

    adjusted_vp <- adjust_vpCoords(viewport = plotVP)


    # =========================================================================
    # SET AT
    # =========================================================================

    if (is.null(bb_yInternal$at)) {
        bb_yInternal$at <- c(plotVP$yscale[1], plotVP$yscale[2])
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Make viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_yaxis",
        length(grep(
            pattern = "bb_yaxis",
            x = currentViewports
        )) + 1
    )

    ## Define viewport
    if (bb_yInternal$main == TRUE) {
        vp <- viewport(
            width = widthDetails(yGrob), height = plotVP$height,
            x = adjusted_vp[[1]] - 0.5 * (plotVP$width),
            y = adjusted_vp[[2]] - 0.5 * (plotVP$height),
            just = c("right", "bottom"), yscale = plotVP$yscale,
            name = vp_name
        )
    } else if (bb_yInternal$main == FALSE) {
        vp <- viewport(
            width = widthDetails(yGrob), height = plotVP$height,
            x = adjusted_vp[[1]] + 0.5 * (plotVP$width),
            y = adjusted_vp[[2]] - 0.5 * (plotVP$height),
            just = c("left", "bottom"), yscale = plotVP$yscale,
            name = vp_name
        )
    }

    # =========================================================================
    # PLOT
    # =========================================================================

    yGrob <- grid.yaxis(
        at = bb_yInternal$at, label = bb_yInternal$label,
        main = bb_yInternal$main, gp = bb_yInternal$gp,
        vp = vp
    )
    if (!bb_yInternal$axisLine) grid.remove(paste0(yGrob$name, "::major"))
    yaxis_grobs <- gTree(vp = vp, children = gList(yGrob))
    yAxis$grobs <- yaxis_grobs

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_yaxis[", vp_name, "]")
    invisible(yAxis)
}
