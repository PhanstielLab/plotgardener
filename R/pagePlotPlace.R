#' Place a plot that has been previously created but not drawn
#' 
#' @usage pagePlotPlace(
#'     plot,
#'     x = NULL,
#'     y = NULL,
#'     width = NULL,
#'     height = NULL,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     draw = TRUE,
#'     params = NULL
#' )
#'
#' @param plot Plot object to be placed, defined by the
#' output of a plotgardener plotting function.
#' @param x A numeric or unit object specifying plot x-location.
#' @param y A numeric, unit object, or character containing a "b" combined
#' with a numeric value specifying plot y-location.
#' The character value will place the plot y relative to the bottom
#' of the most recently plotted plot according to the units
#' of the plotgardener page.
#' @param width A numeric or unit object specifying plot width.
#' @param height A numeric or unit object specifying plot height.
#' @param just Justification of plot relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use
#' if \code{x}, \code{y}, \code{width}, or \code{height} are only
#' given as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output
#' should be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[plotgardener]{params} object
#' containing relevant function parameters.
#'
#' @return Function will update dimensions of an input plot and
#' return an updated plot object.
#'
#' @examples
#' ## Load Hi-C data
#' library(plotgardenerData)
#' data("IMR90_HiC_10kb")
#'
#' ## Create, but do not plot, square Hi-C plot
#' hicPlot <- plotHicSquare(
#'     data = IMR90_HiC_10kb, resolution = 10000,
#'     zrange = c(0, 70),
#'     chrom = "chr21",
#'     chromstart = 28000000, chromend = 30300000,
#'     draw = FALSE
#' )
#'
#' ## Create page
#' pageCreate(width = 3.75, height = 3.5, default.units = "inches")
#'
#' ## Place Hi-C plot on page
#' pagePlotPlace(
#'     plot = hicPlot,
#'     x = 0.25, y = 0.25, width = 3, height = 3,
#'     just = c("left", "top"),
#'     default.units = "inches", draw = TRUE
#' )
#'
#' ## Annotate heatmap legend
#' annoHeatmapLegend(
#'     plot = hicPlot,
#'     x = 3.4, y = 0.25, width = 0.12, height = 1.2,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @export
pagePlotPlace <- function(plot, x = NULL, y = NULL, width = NULL,
                            height = NULL, just = c("left", "top"),
                            default.units = "inches",
                            draw = TRUE, params = NULL) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    set_x <- function(object, val) {
        if (is.null(val)) {
            object["x"] <- list(NULL)
        } else {
            object$x <- val
        }

        return(object)
    }

    set_y <- function(object, val) {
        if (is.null(val)) {
            object["y"] <- list(NULL)
        } else {
            object$y <- val
        }

        return(object)
    }

    set_width <- function(object, val) {
        if (is.null(val)) {
            object["width"] <- list(NULL)
        } else {
            object$width <- val
        }

        return(object)
    }

    set_height <- function(object, val) {
        if (is.null(val)) {
            object["height"] <- list(NULL)
        } else {
            object$height <- val
        }

        return(object)
    }

    set_values <- function(object, x, y, width, height) {
        object <- set_x(object = object, val = x)
        object <- set_y(object = object, val = y)
        object <- set_width(object = object, val = width)
        object <- set_height(object = object, val = height)

        return(object)
    }

    replace_value <- function(val, new) {
        return(new[val])
    }

    ## Define a function to parse coordinates
    parse_coordinates <- function(inputPlot, outputPlot) {

        ## Make sublists of the dimensions/coords of the input/output plots
        inputCoords <- list(
            x = inputPlot$x, y = inputPlot$y,
            width = inputPlot$width, height = inputPlot$height,
            just = inputPlot$just
        )
        outputCoords <- list(
            x = outputPlot$x, y = outputPlot$y,
            width = outputPlot$width,
            height = outputPlot$height, just = outputPlot$just
        )

        ## Determine which values in the output plot are NULL
        to_replace <- names(outputCoords[vapply(
            outputCoords,
            is.null, logical(1)
        )])
        not_replace <- outputCoords[!vapply(
            outputCoords,
            is.null, logical(1)
        )]
        ## Get corresponding values for those that are NULL from the input plot
        replaced <- unlist(lapply(to_replace, replace_value,
            new = inputCoords
        ),
        recursive = FALSE
        )

        ## Recombine values that weren't replaced and those that were
        new_coords <- c(not_replace, replaced)

        outputPlot <- set_values(
            object = outputPlot, x = new_coords$x,
            y = new_coords$y, width = new_coords$width,
            height = new_coords$height
        )
        outputPlot$just <- new_coords$just

        ## Return object
        return(outputPlot)
    }

    ## Define a function that renames grobs and copies to a new gtree
    copy_grobs <- function(grob) {
        new_name <- grobName(grob)

        grob$name <- new_name
        assign("new_gtree",
            addGrob(
                gTree = get("new_gtree", envir = pgEnv),
                child = grob
            ),
            envir = pgEnv
        )
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    place <- parseParams(params = params, 
                            defaultArgs = formals(eval(match.call()[[1]])),
                            declaredArgs = lapply(match.call()[-1], eval.parent,
                                                n = 2),
                            class = "place")
    ## Justification
    place$just <- justConversion(place$just)
    
    # =========================================================================
    # ERRORS
    # =========================================================================

    if (is.null(place$plot)) stop("argument \"plot\" is missing, ",
                                    "with no default.", call. = FALSE)

    # =========================================================================
    # INITIALIZE PLOT OBJECT COPY
    # =========================================================================

    object <- place$plot

    # =========================================================================
    # PARSE UNITS FOR INPUTS
    # =========================================================================

    if (!is.null(place$x)) {
        place$x <- misc_defaultUnits(value = place$x, 
                                        name = "x",
                                        default.units = place$default.units)
    }
    if (!is.null(place$y)) {
        place$y <- misc_defaultUnits(value = place$y, 
                                        name = "y",
                                        default.units = place$default.units)
    }
    if (!is.null(place$width)) {
        place$width <- misc_defaultUnits(value = place$width, 
                                        name = "width",
                                        default.units = place$default.units)
    }
    if (!is.null(place$height)) {
        place$height <- misc_defaultUnits(value = place$height, 
                                            name = "height",
                                            default.units = 
                                                place$default.units)
    }

    # =========================================================================
    # UPDATE DIMENSIONS AND COORDINATES OF PLOT OBJECT BASED ON INPUTS
    # =========================================================================

    object <- set_values(
        object = object, x = place$x, y = place$y,
        width = place$width, height = place$height
    )
    object$just <- place$just
    attr(x = object, which = "plotted") <- place$draw

    # =========================================================================
    # INHERIT DIMENSIONS/COOORDINATES WHERE NULL
    # =========================================================================

    object <- parse_coordinates(
        inputPlot = place$plot,
        outputPlot = object
    )

    # =========================================================================
    # CALL ERRORS
    # =========================================================================

    check_placement(object = object)

    # =========================================================================
    # DEFINE A NEW VIEWPORT
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    num <-
        length(grep(
            pattern = gsub(
                pattern = "[0-9]",
                replacement = "",
                x = object$grobs$vp$name
            ),
            x = currentViewports
        )) + 1

    vp_name <- gsub("[0-9]", replacement = num, x = object$grobs$vp$name)

    ## If full placing information isn't provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(object$x) | is.null(object$y) |
        is.null(object$width) | is.null(object$height)) {
        new_vp <- viewport(
            height = unit(1, "snpc"), width = unit(1, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            clip = "on",
            xscale = object$grobs$vp$xscale,
            yscale = object$grobs$vp$yscale,
            just = "center",
            name = vp_name
        )

        if (place$draw == TRUE) {
            grid.newpage()
            warning("Plot placement will only fill up the graphical device.",
                call. = FALSE
            )
        }
    } else {
        addViewport(vp_name)
        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = object)

        ## Make viewport
        new_vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = object$grobs$vp$xscale,
            yscale = object$grobs$vp$yscale,
            just = place$just,
            name = vp_name
        )
    }

    # =========================================================================
    # RENAME GROBS
    # =========================================================================

    assign("new_gtree", gTree(vp = new_vp), envir = pgEnv)
    invisible(lapply(object$grobs$children, copy_grobs))

    # =========================================================================
    # ASSIGN NEW GROBS TO OBJECT
    # =========================================================================

    object$grobs <- get("new_gtree", envir = pgEnv)

    # =========================================================================
    # IF DRAW == TRUE, DRAW GROBS
    # =========================================================================

    if (place$draw == TRUE) {
        grid.draw(object$grobs)
    }

    # =========================================================================
    # RETURN UPDATED OBJECT
    # =========================================================================

    message(gsub(
        pattern = "[0-9]",
        replacement = "",
        x = vp_name
    ), "[", vp_name, "]")
    invisible(object)
}
