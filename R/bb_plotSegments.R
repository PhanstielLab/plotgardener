#' Draw a line segment within a BentoBox layout
#' 
#' @usage bb_plotSegments(
#'     x0,
#'     y0,
#'     x1,
#'     y1,
#'     default.units = "inches",
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
#' @param x0 A numeric vector or unit object indicating the
#' starting x-values of the line segments.
#' @param y0 A numeric vector, unit object, or a character vector
#' of values containing a "b" combined with a numeric value specifying
#' starting y-values of the line segments.
#' The character vector will place starting y-values relative to the
#' bottom of the most recently plotted BentoBox plot according to the
#' units of the BentoBox page.
#' @param x1 A numeric vector or unit object indicating the stopping
#' x-values of the line segments.
#' @param y1 A numeric vector, unit object, or a character vector of v
#' alues containing a "b" combined with a numeric value specifying
#' stopping y-values of the line segments.
#' The character vector will place stopping y-values relative to the
#' bottom of the most recently plotted BentoBox plot according to the
#' units of the BentoBox page.
#' @param default.units A string indicating the default units to use
#' if \code{x0}, \code{y0}, \code{x1}, or \code{y1} are only given as
#' numeric vectors. Default value is \code{default.units = "inches"}.
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
#' @param arrow A list describing arrow heads to place at either end of
#' the line segments, as produced by the \link[grid]{arrow} function.
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
#' bb_pageCreate(width = 7.5, height = 6, default.units = "inches")
#'
#' ## Plot one line segment
#' bb_plotSegments(
#'     x0 = 3.75, y0 = 0.25, x1 = 3.75, y1 = 5.75,
#'     default.units = "inches",
#'     lwd = 3, lty = 2
#' )
#'
#' ## Plot multiple line segments at different locations in different colors
#' bb_plotSegments(
#'     x0 = 0.5, y0 = c(1, 3, 5), x1 = 3.25, y1 = c(1, 3, 5),
#'     default.units = "inches",
#'     lwd = 2, linecolor = c("#7ecdbb", "#37a7db", "grey")
#' )
#'
#' ## Plot a line segment with an arrowhead
#' bb_plotSegments(
#'     x0 = 4.5, y0 = 0.5, x1 = 7, y1 = 3,
#'     default.units = "inches",
#'     arrow = arrow(type = "closed"), fill = "black"
#' )
#'
#' ## Plot lines with round lineends
#' bb_plotSegments(
#'     x0 = c(4, 7), y0 = 3.5, x1 = 5.5, y1 = 4.5,
#'     default.units = "inches",
#'     lwd = 5, lineend = "round"
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @seealso \link[grid]{grid.segments}, \link[grid]{arrow}
#'
#' @export
bb_plotSegments <- function(x0, y0, x1, y1, default.units = "inches",
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

    check_bbpage(error = "Cannot plot segment without a BentoBox page.")
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
    
    
    bb_segments$x0 <- misc_defaultUnits(
        value = bb_segments$x0,
        name = "x0",
        default.units = bb_segmentsInternal$default.units
    )
    
    bb_segments$y0 <- misc_defaultUnits(
        value = bb_segments$y0,
        name = "y0",
        default.units = bb_segmentsInternal$default.units
    )
    
    bb_segments$x1 <- misc_defaultUnits(
        value = bb_segments$x1,
        name = "x1",
        default.units = bb_segmentsInternal$default.units
    )

    bb_segments$y1 <- misc_defaultUnits(
        value = bb_segments$y1,
        name = "y1",
        default.units = bb_segmentsInternal$default.units
    )

    ## Convert coordinates to page_units
    new_x0 <- convertX(bb_segments$x0, unitTo = page_units, valueOnly = TRUE)
    new_y0 <- convertY(bb_segments$y0, unitTo = page_units, valueOnly = TRUE)
    new_x1 <- convertX(bb_segments$x1, unitTo = page_units, valueOnly = TRUE)
    new_y1 <- convertY(bb_segments$y1, unitTo = page_units, valueOnly = TRUE)

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
        y0 = unit(page_height - new_y0, page_units),
        x1 = unit(new_x1, page_units),
        y1 = unit(page_height - new_y1, page_units),
        arrow = bb_segments$arrow,
        gp = bb_segments$gp,
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
