#' Add an x-axis to a plot
#' 
#' @usage bbAnnoXaxis(
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
#' @param plot Plot object to annotate with x-axis.
#' @param at A numeric vector of x-value locations for tick marks.
#' @param label A logical value indicating whether to draw the labels on
#' the tick marks, or an expression or character vector which specify
#' the labels to use.
#' If not logical, must be the same length as the \code{at} argument.
#' @param main A logical value indicating whether to draw the x-axis at the
#' bottom of the plot. Default value is \code{main = TRUE}. Options are:
#' \itemize{
#' \item{\code{TRUE}: }{x-axis is drawn at the bottom of the plot.}
#' \item{\code{FALSE}: }{x-axis is drawn at the top of the plot.}
#' }
#' @param scipen An integer indicating the penalty to be applied when
#' deciding to print numeric values in fixed or exponential notation.
#' Default value is \code{scipen = 999}.
#' @param axisLine A logical value indicating whether to show the axis line.
#' Default value is \code{axisLine = FALSE}.
#' @param params An optional \link[BentoBox]{bbParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_xaxis} object containing
#' relevant \link[grid]{grob} information.
#'
#' @examples
#' ## Load transcript information
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Create BentoBox page
#' bbPageCreate(width = 7.5, height = 4.5, default.units = "inches")
#'
#' ## Plot gene transcripts
#' transcriptPlot <- bbPlotTranscripts(
#'     chrom = "chr1",
#'     chromstart = 1000000,
#'     chromend = 2000000,
#'     assembly = "hg19",
#'     x = 0.5, y = 0,
#'     width = 6.5, height = 4,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Add standard x-axis to transcript plot
#' bbAnnoXaxis(
#'     plot = transcriptPlot,
#'     at = c(1000000, 1250000, 1500000, 1750000, 2000000),
#'     fontsize = 8
#' )
#' bbPlotText(
#'     label = "Basepairs", fontsize = 10, fontface = "bold",
#'     x = 3.75, y = 4.3, just = "top"
#' )
#'
#' ## Hide page guides
#' bbPageGuideHide()
#' @export
bbAnnoXaxis <- function(plot, at = NULL, label = TRUE, main = TRUE,
                        scipen = 999, axisLine = FALSE, params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_xInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_xInternal"
    )

    ## Set gp
    bb_xInternal$gp <- setGP(
        gpList = gpar(),
        params = bb_xInternal, ...
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(bb_xInternal$plot)) stop("argument \"plot\" is missing, ",
                                        "with no default.", call. = FALSE)
    check_bbpage(error = "Cannot add an x-axis without a BentoBox page.")

    # =========================================================================
    # SAVE USER'S SCIPEN AND SET OURS TEMPORARILY
    # =========================================================================

    oo <- options(scipen = bb_xInternal$scipen)
    on.exit(options(oo))

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    xAxis <- structure(list(grobs = NULL), class = "bb_xaxis")

    # =========================================================================
    # CREATE GROB WITHOUT DRAWING
    # =========================================================================

    xGrob <- xaxisGrob(
        at = bb_xInternal$at, label = bb_xInternal$label,
        main = bb_xInternal$main, gp = bb_xInternal$gp,
        vp = bb_xInternal$plot$grobs$vp
    )

    # =========================================================================
    # GET CENTER OF INPUT PLOT VIEWPORT BASED ON INPUT PLOT TYPE AND JUST
    # =========================================================================

    # Get appropriate plot viewport
    plotVP <- getAnnoViewport(plot = bb_xInternal$plot)

    adjusted_vp <- adjust_vpCoords(viewport = plotVP)

    # =========================================================================
    # SET AT
    # =========================================================================

    if (is.null(bb_xInternal$at)) {
        bb_xInternal$at <- c(plotVP$xscale[1], plotVP$xscale[2])
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Make viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_xaxis",
        length(grep(
            pattern = "bb_xaxis",
            x = currentViewports
        )) + 1
    )

    ## Define viewport
    if (bb_xInternal$main == TRUE) {
        vp <- viewport(
            width = plotVP$width, height = heightDetails(xGrob),
            x = adjusted_vp[[1]] - 0.5 * (plotVP$width),
            y = adjusted_vp[[2]] - 0.5 * (plotVP$height),
            just = c("left", "top"), xscale = plotVP$xscale,
            name = vp_name
        )
    } else if (bb_xInternal$main == FALSE) {
        vp <- viewport(
            width = plotVP$width, height = heightDetails(xGrob),
            x = adjusted_vp[[1]] - 0.5 * (plotVP$width),
            y = adjusted_vp[[2]] + 0.5 * (plotVP$height),
            just = c("left", "bottom"), xscale = plotVP$xscale,
            name = vp_name
        )
    }


    # =========================================================================
    # PLOT
    # =========================================================================

    xGrob <- grid.xaxis(
        at = bb_xInternal$at, label = bb_xInternal$label,
        main = bb_xInternal$main, gp = bb_xInternal$gp, vp = vp
    )
    if (!bb_xInternal$axisLine) grid.remove(paste0(xGrob$name, "::major"))
    xaxis_grobs <- gTree(vp = vp, children = gList(xGrob))
    xAxis$grobs <- xaxis_grobs

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_xaxis[", vp_name, "]")
    invisible(xAxis)
}
