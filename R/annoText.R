#' Annotates text within a plot
#' 
#' @usage annoText(
#'     label,
#'     fontcolor = "black",
#'     fontsize = 12,
#'     rot = 0,
#'     check.overlap = FALSE,
#'     plot,
#'     x,
#'     y,
#'     just = "center",
#'     default.units = "native",
#'     params = NULL,
#'     ...
#' )
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
#' @param plot Input plotgardener plot to internally place text relative to.
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
#' @param params An optional \link[plotgardener]{params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{text} object containing
#' relevant placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' pageCreate(width = 4, height = 4, default.units = "inches")
#'
#' ## Plot text relative to a plotgardener plot
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
#' hicPlot <- plotHicSquare(
#'     data = IMR90_HiC_10kb, chrom = "chr21",
#'     chromstart = 28000000, chromend = 29500000,
#'     assembly = "hg19",
#'     zrange = c(0, 70),
#'     x = 0.5, y = 0.5, width = 3, height = 3,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#' annoGenomeLabel(
#'     plot = hicPlot, x = 0.5, y = 3.55, scale = "Mb",
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' annoText(
#'     label = "Loop", fontsize = 8, plot = hicPlot,
#'     x = 29075000, y = 28150000,
#'     just = "center", default.units = "native"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @seealso \link[grid]{grid.text}
#'
#' @export
annoText <- function(label, fontcolor = "black", fontsize = 12, rot = 0,
                        check.overlap = FALSE, plot, x, y, just = "center",
                        default.units = "native", params = NULL, ...) {

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    textInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "textInternal"
    )

    ## Set gp
    textInternal$gp <- gpar(
        col = textInternal$fontcolor,
        fontsize = textInternal$fontsize
    )
    textInternal$gp <- setGP(
        gpList = textInternal$gp,
        params = textInternal, ...
    )
    
    ## Justification
    textInternal$just <- justConversion(just = textInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    text <- structure(list(
        label = textInternal$label,
        x = textInternal$x, y = textInternal$y,
        just = textInternal$just, grobs = NULL
    ),
    class = "text"
    )

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    check_page(error = "Cannot annotate text without a `plotgardener` page.")
    if (is.null(text$label)) stop("argument \"label\" is missing, ",
                                    "with no default.", call. = FALSE)
    if (is.null(textInternal$plot)) stop("argument \"plot\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(text$x)) stop("argument \"x\" is missing, ",
                                "with no default.", call. = FALSE)
    if (is.null(text$y)) stop("argument \"y\" is missing, ",
                                "with no default.", call. = FALSE)

    # =========================================================================
    # DEFINE PARAMETERS
    # =========================================================================

    ## Get page_height and its units from pgEnv through bb_pageCreate
    page_height <- get("page_height", envir = pgEnv)
    page_units <- get("page_units", envir = pgEnv)
    
    text$x <- misc_defaultUnits(value = text$x,
                                name = "x",
                                default.units = 
                                    textInternal$default.units)
    text$y <- misc_defaultUnits(value = text$y,
                                name = "y",
                                default.units = 
                                    textInternal$default.units,
                                funName = "annoText",
                                yBelow = FALSE)

    # Get appropriate plot viewport
    plotVP <- getAnnoViewport(plot = textInternal$plot)
    
    ## Convert plot viewport to bottom left to get position on entire page
    plotVP_bottomLeft <- vp_bottomLeft(viewport = plotVP)

    ## Push plot viewport to convert x/y from plot native units to page units
    seekViewport(plotVP$name)
    new_x <- convertX(text$x, unitTo = page_units, valueOnly = TRUE)
    new_y <- convertY(text$y, unitTo = page_units, valueOnly = TRUE)

    seekViewport(name = "page")

    ## Add additional page units to new_x and new_y
    new_x <- as.numeric(plotVP_bottomLeft[[1]]) + new_x
    new_y <- as.numeric(plotVP_bottomLeft[[2]]) + new_y

    # =========================================================================
    # MAKE GROB
    # =========================================================================
    name <- paste0(
        "text",
        length(grep(
            pattern = "text",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    text <- grid.text(
        label = text$label, x = unit(new_x, page_units),
        y = unit(new_y, page_units), just = text$just,
        gp = textInternal$gp, rot = textInternal$rot,
        check.overlap = textInternal$check.overlap,
        name = name
    )

    # =========================================================================
    # ADD GROB TO OBJECT
    # =========================================================================

    text$grobs <- text

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("text[", text$name, "]")
    invisible(text)
}
