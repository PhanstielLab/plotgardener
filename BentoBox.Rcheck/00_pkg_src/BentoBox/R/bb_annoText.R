#' Annotates text within a BentoBox plot
#'
#' @param label Character or expression of text to be plotted.
#' @param fontcolor A character value specifying text fontcolor.
#' Default value is \code{fontcolor = "black"}.
#' @param fontsize A numeric specifying text fontsize in points.
#' Default value is \code{fontsize = 12}.
#' @param rot A numeric specifying the angle to rotate the text.
#' Default value is \code{rot = 0}.
#' @param check.overlap A logical value to indicate whether to check for
#' and omit overlapping text. Default value is \code{check.overlap = FALSE}.
#' @param plot Input BentoBox plot to internally place text relative to.
#' @param x A numeric vector or unit object specifying text x-location.
#' @param y A numeric vector or unit object specifying text y-location.
#' @param just Justification of text relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = "center"}.
#' @param default.units A string indicating the default units to use if
#' \code{x} or \code{y} are only given as numerics or numeric vectors.
#' Default value is \code{default.units = "native"}.
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_text} object containing
#' relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a BentoBox page
#' bb_pageCreate(width = 4, height = 4, default.units = "inches")
#'
#' ## Plot text relative to a BentoBox plot
#' data("bb_imrHicData")
#' hicPlot <- bb_plotHicSquare(
#'     data = bb_imrHicData, chrom = "chr21",
#'     chromstart = 28000000, chromend = 29500000,
#'     zrange = c(0, 70),
#'     x = 0.5, y = 0.5, width = 3, height = 3,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#' bb_annoGenomeLabel(
#'     plot = hicPlot, x = 0.5, y = 3.55, scale = "Mb",
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' bb_annoText(
#'     label = "Loop", fontsize = 8, plot = hicPlot,
#'     x = 29075000, y = 28150000,
#'     just = "center", default.units = "native"
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @seealso \link[grid]{grid.text}
#'
#' @export
bb_annoText <- function(label, fontcolor = "black", fontsize = 12, rot = 0,
                        check.overlap = FALSE, plot, x, y, just = "center",
                        default.units = "native", params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    ## Check which defaults are not overwritten and set to NULL
    if (missing(just)) just <- NULL
    if (missing(fontcolor)) fontcolor <- NULL
    if (missing(fontsize)) fontsize <- NULL
    if (missing(rot)) rot <- NULL
    if (missing(check.overlap)) check.overlap <- NULL
    if (missing(default.units)) default.units <- NULL

    ## Check if label/x/y/plot arguments are missing (could be in object)
    if (!hasArg(label)) label <- NULL
    if (!hasArg(plot)) plot <- NULL
    if (!hasArg(x)) x <- NULL
    if (!hasArg(y)) y <- NULL

    ## Compile all parameters into an internal object
    bb_textInternal <- structure(list(
        label = label, x = x, y = y, plot = plot,
        just = just, fontcolor = fontcolor,
        fontsize = fontsize, rot = rot,
        check.overlap = check.overlap,
        default.units = default.units
    ),
    class = "bb_textInternal"
    )

    bb_textInternal <- parseParams(
        bb_params = params,
        object_params = bb_textInternal
    )

    ## For any defaults that are still NULL, set back to default
    if (is.null(bb_textInternal$just)) {
        bb_textInternal$just <- "center"
    }
    if (is.null(bb_textInternal$fontcolor)) {
        bb_textInternal$fontcolor <- "black"
    }
    if (is.null(bb_textInternal$fontsize)) {
        bb_textInternal$fontsize <- 12
    }
    if (is.null(bb_textInternal$rot)) {
        bb_textInternal$rot <- 0
    }
    if (is.null(bb_textInternal$check.overlap)) {
        bb_textInternal$check.overlap <- FALSE
    }
    if (is.null(bb_textInternal$default.units)) {
        bb_textInternal$default.units <- "native"
    }

    ## Set gp
    bb_textInternal$gp <- gpar(
        col = bb_textInternal$fontcolor,
        fontsize = bb_textInternal$fontsize
    )
    bb_textInternal$gp <- setGP(
        gpList = bb_textInternal$gp,
        params = bb_textInternal, ...
    )

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    bb_text <- structure(list(
        label = bb_textInternal$label,
        x = bb_textInternal$x, y = bb_textInternal$y,
        just = bb_textInternal$just, grobs = NULL
    ),
    class = "bb_text"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_bbpage(error = "Cannot annotate text without a BentoBox page.")
    if (is.null(bb_text$label)) stop("argument \"label\" is missing,
                                    with no default.", call. = FALSE)
    if (is.null(bb_textInternal$plot)) stop("argument \"plot\" is missing,
                                            with no default.", call. = FALSE)
    if (is.null(bb_text$x)) stop("argument \"x\" is missing,
                                with no default.", call. = FALSE)
    if (is.null(bb_text$y)) stop("argument \"y\" is missing,
                                with no default.", call. = FALSE)

    # =========================================================================
    # DEFINE PARAMETERS
    # =========================================================================

    ## Get page_height and its units from bbEnv through bb_pageCreate
    page_height <- get("page_height", envir = bbEnv)
    page_units <- get("page_units", envir = bbEnv)

    if (!"unit" %in% class(bb_text$x)) {
        if (!is.numeric(bb_text$x)) {
            stop("x-coordinate is neither a unit object or a numeric value.
                Cannot place object.", call. = FALSE)
        }

        if (is.null(bb_textInternal$default.units)) {
            stop("x-coordinate detected as numeric.\'default.units\'
                must be specified.", call. = FALSE)
        }

        bb_text$x <- unit(bb_text$x, bb_textInternal$default.units)
    }

    if (!"unit" %in% class(bb_text$y)) {

        ## Check for "below" y-coord
        if (grepl("b", bb_text$y) == TRUE) {
            stop("\'below\' y-coordinate detected. Cannot parse \'below\'
                y-coordinate for bb_annoText.", call. = FALSE)
        } else {
            if (!is.numeric(bb_text$y)) {
                stop("y-coordinate is neither a unit object or a
                    numeric value. Cannot place object.", call. = FALSE)
            }

            if (is.null(bb_textInternal$default.units)) {
                stop("y-coordinate detected as numeric.\'default.units\'
                    must be specified.", call. = FALSE)
            }

            bb_text$y <- unit(bb_text$y, bb_textInternal$default.units)
        }
    }

    if (class(bb_textInternal$plot) == "bb_genes") {
        plotVP <- bb_textInternal$plot$grobs$children$background$vp
    } else if (class(bb_textInternal$plot) == "bb_hicTriangle" |
        class(bb_textInternal$plot) == "bb_hicRectangle") {
        plotVP <- bb_textInternal$plot$outsideVP
    } else {
        plotVP <- bb_textInternal$plot$grobs$vp
    }

    ## Convert plot viewport to bottom left to get position on entire page
    plotVP_bottomLeft <- vp_bottomLeft(viewport = plotVP)

    ## Push plot viewport to convert x/y from plot native units to page units
    seekViewport(plotVP$name)
    new_x <- convertX(bb_text$x, unitTo = page_units, valueOnly = TRUE)
    new_y <- convertY(bb_text$y, unitTo = page_units, valueOnly = TRUE)

    seekViewport(name = "bb_page")

    ## Add additional page units to new_x and new_y
    new_x <- as.numeric(plotVP_bottomLeft[[1]]) + new_x
    new_y <- as.numeric(plotVP_bottomLeft[[2]]) + new_y

    # =========================================================================
    # MAKE GROB
    # =========================================================================
    name <- paste0(
        "bb_text",
        length(grep(
            pattern = "bb_text",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    text <- grid.text(
        label = bb_text$label, x = unit(new_x, page_units),
        y = unit(new_y, page_units), just = bb_text$just,
        gp = bb_textInternal$gp, rot = bb_textInternal$rot,
        check.overlap = bb_textInternal$check.overlap,
        name = name
    )

    # =========================================================================
    # ADD GROB TO OBJECT
    # =========================================================================

    bb_text$grobs <- text

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_text[", text$name, "]")
    invisible(bb_text)
}
