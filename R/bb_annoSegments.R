#' Annotates a line segment within a BentoBox plot
#' 
#' @usage bb_annoSegments(
#'     x0,
#'     y0,
#'     x1,
#'     y1,
#'     plot,
#'     default.units = "native",
#'     linecolor = "black",
#'     lwd = 1,
#'     lty = 1,
#'     lineend = "butt",
#'     linejoin = "mitre",
#'     arrow = NULL,
#'     params = NULL,
#'     ...
#' )
#'
#' @param x0 A numeric vector or unit object indicating the starting
#' x-values of the line segments.
#' @param y0 A numeric vector or unit object indicating the starting
#' y-values of the line segments.
#' @param x1 A numeric vector or unit object indicating the stopping
#' x-values of the line segments.
#' @param y1 A numeric vector or unit object indicating the stopping
#' y-values of the line segments.
#' @param plot Input BentoBox plot to internally plot line segments
#' relative to.
#' @param default.units A string indicating the default units to use
#' if \code{x0}, \code{y0}, \code{x1}, or \code{y1} are only given
#' as numeric vectors. Default value is \code{default.units = "native"}.
#' @param linecolor A character value specifying segment line color.
#' Default value is \code{linecolor = "black"}.
#' @param lwd A numeric specifying segment line width.
#' Default value is \code{lwd = 1}.
#' @param lty A numeric specifying segment line type.
#' Default value is \code{lty = 1}.
#' @param lineend A character value specifying line end style.
#' Default value is \code{lineend = "butt"}. Options are:
#' \itemize{
#' \item{\code{"round"}: Segment ends are rounded.}
#' \item{\code{"butt"}: Segment ends end exactly where ended.}
#' \item{\code{"square"}: Segment ends are squared.}
#' }
#' @param linejoin A character value specifying line join style.
#' Default value is \code{linejoin = "mitre"}. Options are:
#' \itemize{
#' \item{\code{"round"}: }{Line joins are rounded.}
#' \item{\code{"mitre"}: }{Line joins are sharp corners.}
#' \item{\code{"bevel"}: }{Line joins are flattened corners.}
#' }
#' @param arrow A list describing arrow heads to place at either
#' end of the line segments, as produced by the \link[grid]{arrow} function.
#' @param params An optional \link[BentoBox]{bb_params} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_segments} object containing relevant
#' placement and \link[grid]{grob} information.
#'
#' @examples
#' library(grid)
#' ## Create a BentoBox page
#' bb_pageCreate(width = 7.5, height = 2.5, default.units = "inches")
#'
#' ## Plot a Manhattan plot
#' library(BentoBoxData)
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' data("hg19_insulin_GWAS")
#' manhattanPlot <- bb_plotManhattan(
#'     data = hg19_insulin_GWAS, assembly = "hg19",
#'     fill = c("grey", "#37a7db"),
#'     sigLine = TRUE,
#'     col = "grey", lty = 2, range = c(0, 14),
#'     x = 0.5, y = 0, width = 6.5, height = 2,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bb_annoGenomeLabel(
#'     plot = manhattanPlot, x = 0.5, y = 2, fontsize = 8,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#' bb_plotText(
#'     label = "Chromosome", fontsize = 8,
#'     x = 3.75, y = 2.20, just = "center", default.units = "inches"
#' )
#'
#' ## Annotate y-axis
#' bb_annoYaxis(
#'     plot = manhattanPlot, at = c(0, 2, 4, 6, 8, 10, 12, 14),
#'     axisLine = TRUE, fontsize = 8
#' )
#'
#' ## Annotate a line segment for an additional significance line of
#' ## the Manhattan plot
#' bb_annoSegments(
#'     x0 = unit(0, "npc"), y0 = 10,
#'     x1 = unit(1, "npc"), y1 = 10,
#'     plot = manhattanPlot, default.units = "native",
#'     linecolor = "red", lty = 2
#' )
#'
#' ## Plot y-axis label
#' bb_plotText(
#'     label = "-log10(p-value)", x = 0.15, y = 1, rot = 90,
#'     fontsize = 8, fontface = "bold", just = "center",
#'     default.units = "inches"
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @seealso \link[grid]{grid.segments}, \link[grid]{arrow}
#' @export
bb_annoSegments <- function(x0, y0, x1, y1, plot, default.units = "native",
                            linecolor = "black", lwd = 1, lty = 1,
                            lineend = "butt", linejoin = "mitre",
                            arrow = NULL, params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_segmentsInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_segmentsInternal"
    )

    ## Set gp
    bb_segmentsInternal$gp <- gpar(
        col = bb_segmentsInternal$linecolor,
        lwd = bb_segmentsInternal$lwd,
        lty = bb_segmentsInternal$lty,
        lineend = bb_segmentsInternal$lineend,
        linejoin = bb_segmentsInternal$linejoin
    )
    bb_segmentsInternal$gp <- setGP(
        gpList = bb_segmentsInternal$gp,
        params = bb_segmentsInternal, ...
    )

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    bb_segments <- structure(list(
        x0 = bb_segmentsInternal$x0,
        y0 = bb_segmentsInternal$y0,
        x1 = bb_segmentsInternal$x1,
        y1 = bb_segmentsInternal$y1,
        arrow = bb_segmentsInternal$arrow,
        grobs = NULL,
        gp = bb_segmentsInternal$gp
    ),
    class = "bb_segmentsInternal"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_bbpage(error = "Cannot annotate segment without a BentoBox page.")
    if (is.null(bb_segmentsInternal$plot)) {
        stop("argument \"plot\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(bb_segments$x0)) stop("argument \"x0\" is missing, ",
                                    "with no default.", call. = FALSE)
    if (is.null(bb_segments$y0)) stop("argument \"y0\" is missing, ",
                                    "with no default.", call. = FALSE)
    if (is.null(bb_segments$x1)) stop("argument \"x1\" is missing, ",
                                    "with no default.", call. = FALSE)
    if (is.null(bb_segments$y1)) stop("argument \"y1\" is missing, ",
                                    "with no default.", call. = FALSE)

    # =========================================================================
    # DEFINE PARAMETERS
    # =========================================================================

    ## Get page_height and its units from bbEnv through bb_pageCreate
    page_height <- get("page_height", envir = bbEnv)
    page_units <- get("page_units", envir = bbEnv)
    
    bb_segments$x0 <- misc_defaultUnits(value = bb_segments$x0,
                                        name = "x0",
                                        default.units = 
                                            bb_segmentsInternal$default.units)
    bb_segments$y0 <- misc_defaultUnits(value = bb_segments$y0,
                                        name = "y0",
                                        default.units = 
                                            bb_segmentsInternal$default.units,
                                        funName = "bb_annoSegments",
                                        yBelow = FALSE)
    bb_segments$x1 <- misc_defaultUnits(value = bb_segments$x1,
                                        name = "x1",
                                        default.units = 
                                            bb_segmentsInternal$default.units)
    bb_segments$y1 <- misc_defaultUnits(value = bb_segments$y1,
                                        name = "y1",
                                        default.units = 
                                            bb_segmentsInternal$default.units,
                                        funName = "bb_annoSegments",
                                        yBelow = FALSE)
    
    ## Get appropriate plot viewport
    plotVP <- get_annoViewport(plot = bb_segmentsInternal$plot)

    ## Convert plot viewport to bottom left to get position on entire page
    plotVP_bottomLeft <- vp_bottomLeft(viewport = plotVP)

    ## Push plot viewport to convert x/y from plot native units to page units
    seekViewport(plotVP$name)
    new_x0 <- convertX(bb_segments$x0, unitTo = page_units, valueOnly = TRUE)
    new_x1 <- convertX(bb_segments$x1, unitTo = page_units, valueOnly = TRUE)
    new_y0 <- convertY(bb_segments$y0, unitTo = page_units, valueOnly = TRUE)
    new_y1 <- convertY(bb_segments$y1, unitTo = page_units, valueOnly = TRUE)

    seekViewport(name = "bb_page")

    ## Add additional page units to new_x0, new_x1, new_y0, and new_y1
    new_x0 <- as.numeric(plotVP_bottomLeft[[1]]) + new_x0
    new_x1 <- as.numeric(plotVP_bottomLeft[[1]]) + new_x1
    new_y0 <- as.numeric(plotVP_bottomLeft[[2]]) + new_y0
    new_y1 <- as.numeric(plotVP_bottomLeft[[2]]) + new_y1

    # =========================================================================
    # MAKE GROB
    # =========================================================================
    name <- paste0(
        "bb_segments",
        length(grep(
            pattern = "bb_segments",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    segments <- grid.segments(
        x0 = unit(new_x0, page_units),
        y0 = unit(new_y0, page_units),
        x1 = unit(new_x1, page_units),
        y1 = unit(new_y1, page_units),
        arrow = bb_segments$arrow, gp = bb_segments$gp,
        name = name
    )

    # =========================================================================
    # ADD GROB TO OBJECT
    # =========================================================================

    bb_segments$grobs <- segments

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_segments[", segments$name, "]")
    invisible(bb_segments)
}
