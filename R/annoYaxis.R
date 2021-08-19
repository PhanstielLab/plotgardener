#' Add a y-axis to a plot
#' 
#' @usage annoYaxis(
#'     plot,
#'     at = NULL,
#'     label = TRUE,
#'     main = TRUE,
#'     scipen = 999,
#'     axisLine = FALSE,
#'     params = NULL,
#'     ...
#' )
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
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{yaxis} object containing
#' relevant \link[grid]{grob} information.
#'
#' @examples
#' ## Load Hi-C data
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create page
#' pageCreate(width = 4, height = 3.5, default.units = "inches")
#'
#' ## Plot and place a square Hi-C plot
#' hicPlot <- plotHicSquare(
#'     data = IMR90_HiC_10kb, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     assembly = "hg19",
#'     x = 1, y = 0.5, width = 2.5, height = 2.5,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Add standard y-axis to Hi-C plot
#' annoYaxis(
#'     plot = hicPlot, at = c(28000000, 29000000, 30300000),
#'     fontsize = 10
#' )
#'
#' ## Annotate genome label on x-axis
#' annoGenomeLabel(plot = hicPlot, x = 1, y = 3.03)
#'
#' ## Annotate heatmap legend
#' annoHeatmapLegend(
#'     plot = hicPlot,
#'     x = 3.6, y = 0.5, width = 0.12, height = 1.2
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @export
annoYaxis <- function(plot, at = NULL, label = TRUE, main = TRUE,
                        scipen = 999, axisLine = FALSE, params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    yInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "yInternal"
    )

    ## Set gp
    yInternal$gp <- setGP(
        gpList = gpar(),
        params = yInternal, ...
    )
    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(yInternal$plot)) stop("argument \"plot\" is missing, ",
                                        "with no default.", call. = FALSE)
    check_page(error = "Cannot add an y-axis without a `plotgardener` page.")

    # =========================================================================
    # SAVE USER'S SCIPEN AND SET OURS TEMPORARILY
    # =========================================================================

    oo <- options(scipen = yInternal$scipen)
    on.exit(options(oo))

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    yAxis <- structure(list(grobs = NULL), class = "yaxis")

    # =========================================================================
    # CREATE GROB WITHOUT DRAWING
    # =========================================================================

    yGrob <- yaxisGrob(
        at = yInternal$at, label = yInternal$label,
        main = yInternal$main, gp = yInternal$gp,
        vp = yInternal$plot$grobs$vp
    )

    # =========================================================================
    # GET CENTER OF INPUT PLOT VIEWPORT BASED ON INPUT PLOT TYPE AND JUST
    # =========================================================================

    # Get appropriate plot viewport
    plotVP <- getAnnoViewport(plot = yInternal$plot)
    
    adjusted_vp <- adjust_vpCoords(viewport = plotVP)

    # =========================================================================
    # SET AT
    # =========================================================================

    if (is.null(yInternal$at)) {
        yInternal$at <- c(plotVP$yscale[1], plotVP$yscale[2])
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Make viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "yaxis",
        length(grep(
            pattern = "yaxis",
            x = currentViewports
        )) + 1
    )

    ## Define viewport
    if (yInternal$main == TRUE) {
        vp <- viewport(
            width = widthDetails(yGrob), height = plotVP$height,
            x = adjusted_vp[[1]] - 0.5 * (plotVP$width),
            y = adjusted_vp[[2]] - 0.5 * (plotVP$height),
            just = c("right", "bottom"), yscale = plotVP$yscale,
            name = vp_name
        )
    } else if (yInternal$main == FALSE) {
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
        at = yInternal$at, label = yInternal$label,
        main = yInternal$main, gp = yInternal$gp,
        vp = vp
    )
    if (!yInternal$axisLine) grid.remove(paste0(yGrob$name, "::major"))
    yaxis_grobs <- gTree(vp = vp, children = gList(yGrob))
    yAxis$grobs <- yaxis_grobs

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("yaxis[", vp_name, "]")
    invisible(yAxis)
}
