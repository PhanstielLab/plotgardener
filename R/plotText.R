#' Plot text within a plotgardener layout
#' 
#' @usage plotText(
#'     label,
#'     fontcolor = "black",
#'     fontsize = 12,
#'     rot = 0,
#'     check.overlap = FALSE,
#'     x,
#'     y,
#'     just = "center",
#'     default.units = "inches",
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
#' @param check.overlap A logical value to indicate whether to check
#' for and omit overlapping text.
#' Default value is \code{check.overlap = FALSE}.
#' @param x A numeric vector or unit object specifying text x-location.
#' @param y A numeric vector, unit object, or a character vector of
#' values containing a "b" combined with a numeric value
#' specifying text y-locations.
#' The character vector will place text y-locations relative to the
#' bottom of the most recently plotted plot according to the
#' units of the plotgardener page.
#' @param just Justification of text relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = "center"}.
#' @param default.units A string indicating the default units to use if
#' \code{x} or \code{y} are only given as numerics or numeric vectors.
#' Default value is \code{default.units = "inches"}.
#' @param params An optional \link[plotgardener]{params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{text} object containing relevant
#' placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' pageCreate(width = 4, height = 2, default.units = "inches")
#'
#' ## Plot text, adjusting fontsize and fontface
#' plotText(
#'     label = "plotgardener", fontsize = 14, fontface = "bold",
#'     x = 1, y = 1, just = "center", default.units = "inches"
#' )
#'
#' ## Plot text, adjusting color, rotation, and fontfamily
#' plotText(
#'     label = "coordinate-based", fontcolor = "#225EA8", rot = 90,
#'     fontfamily = "HersheyScript", x = 2, y = 1, just = "center",
#'     default.units = "inches"
#' )
#'
#' ## Plot a text label in multiple places at once
#' plotText(
#'     label = "R", x = c(0.5, 1, 1.5), y = 1.5, just = "center",
#'     default.units = "inches"
#' )
#'
#' ## Plot a vector of text labels
#' plotText(
#'     label = c("bb", "Bento", "Box"), x = 3, y = c(0.5, 1, 1.75),
#'     just = "center", default.units = "inches"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @seealso \link[grid]{grid.text}
#'
#' @export
plotText <- function(label, fontcolor = "black", fontsize = 12, rot = 0,
                        check.overlap = FALSE, x, y, just = "center",
                        default.units = "inches", params = NULL, ...) {


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

    check_page(error = "Cannot plot text without a `plotgardener` page.")
    if (is.null(text$label)) stop("argument \"label\" is missing, ",
                                    "with no default.", call. = FALSE)
    if (is.null(text$x)) {
        stop("argument \"x\" is missing, with no default.",
            call. = FALSE
        )
    }
    if (is.null(text$y)) {
        stop("argument \"y\" is missing, with no default.",
            call. = FALSE
        )
    }

    # =========================================================================
    # DEFINE PARAMETERS
    # =========================================================================

    ## Get page_height and its units from pgEnv
    page_height <- get("page_height", envir = pgEnv)
    page_units <- get("page_units", envir = pgEnv)
    
    text$x <- misc_defaultUnits(
        value = text$x,
        name = "x",
        default.units = textInternal$default.units
    )
    text$y <- misc_defaultUnits(
        value = text$y,
        name = "y",
        default.units = textInternal$default.units
    )

    ## Convert coordinates to page_units
    new_x <- convertX(text$x, unitTo = page_units, valueOnly = TRUE)
    new_y <- page_height - convertY(text$y,
        unitTo = page_units,
        valueOnly = TRUE
    )

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
