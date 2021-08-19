#' Annotate genomic coordinates along the x or y-axis of a plot
#' 
#' @usage annoGenomeLabel(
#'     plot,
#'     fontsize = 10,
#'     fontcolor = "black",
#'     linecolor = "black",
#'     margin = unit(1, "mm"),
#'     scale = "bp",
#'     commas = TRUE,
#'     sequence = TRUE,
#'     boxWidth = 0.5,
#'     axis = "x",
#'     at = NULL,
#'     tcl = 0.5,
#'     x,
#'     y,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     params = NULL,
#'     ...
#' )
#'
#' @param plot Input plot to annotate genomic coordinates.
#' Genomic coordinates and assembly will be inherited from \code{plot}.
#' @param fontsize A numeric specifying text fontsize in points.
#' Default value is \code{fontsize = 10}.
#' @param fontcolor A character value indicating the color for text.
#' Default value is \code{fontcolor = "black"}.
#' @param linecolor A character value indicating the color of
#' the genome label axis. Default value is \code{linecolor = "black"}.
#' @param margin A numeric or unit vector specifying space between axis
#' and coordinate labels. Default value is \code{margin = unit(1, "mm")}.
#' @param scale A character value indicating the scale of the coordinates
#' along the genome label. Default value is \code{scale = "bp"}. Options are:
#' \itemize{
#' \item{\code{"bp"}: }{base pairs.}
#' \item{\code{"Kb"}: }{kilobase pairs. 1 kilobase pair is equal to
#' 1000 base pairs.}
#' \item{\code{"Mb"}: }{megabase pairs. 1 megabase pair is equal to
#' 1000000 base pairs.}
#' }
#' @param commas A logical value indicating whether to include commas in
#' start and stop labels. Default value is \code{commas = TRUE}.
#' @param sequence A logical value indicating whether to include sequence
#' information above the label of an x-axis (only at appropriate resolutions).
#' @param boxWidth A numeric value indicating the width of the boxes
#' representing sequence information at appropriate resolutions.
#' Default value is \code{boxWidth = 0.5}.
#' @param axis A character value indicating along which axis to
#' add genome label. Sequence information will not be displayed along a y-axis.
#' Default value is \code{axis = "x"}.
#' Options are:
#' \itemize{
#' \item{\code{"x"}: }{Genome label will be plotted along the x-axis.}
#' \item{\code{"y"}: }{Genome label will be plotted along the y-axis.
#' This is typically used for a square Hi-C plot made with
#' \code{plotHicSquare}.}
#' }
#' @param at A numeric vector of x-value locations for tick marks.
#' @param tcl A numeric specifying the length of tickmarks as a fraction of
#' text height. Default value is \code{tcl = 0.5}.
#' @param x A numeric or unit object specifying genome label x-location.
#' @param y A numeric, unit object, or character containing a "b" combined
#' with a numeric value specifying genome label y-location.
#' The character value will place the genome label y relative to the bottom
#' of the most recently plotted plot according to the units of the
#' plotgardener page.
#' @param just Justification of genome label relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal justification
#' and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"},
#' \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use
#' if \code{x} or \code{y} are only given as numerics.
#' Default value is \code{default.units = "inches"}.
#' @param params An optional \link[plotgardener]{params} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters or digit specifications.
#' See \link[grid]{gpar} and \link[base]{formatC}.
#'
#' @return Returns a \code{genomeLabel} object containing
#' relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load hg19 genomic annotation packages
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Create page
#' pageCreate(width = 5, height = 2, default.units = "inches")
#'
#' ## Plot and place gene track on page
#' genesPlot <- plotGenes(
#'     chrom = "chr8",
#'     chromstart = 1000000, chromend = 2000000,
#'     assembly = "hg19", fill = c("grey", "grey"),
#'     fontcolor = c("grey", "grey"),
#'     x = 0.5, y = 0.25, width = 4, height = 1,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate x-axis genome labels at different scales
#' annoGenomeLabel(
#'     plot = genesPlot, scale = "Mb",
#'     x = 0.5, y = 1.25, just = c("left", "top"),
#'     default.units = "inches"
#' )
#' annoGenomeLabel(
#'     plot = genesPlot, scale = "Kb",
#'     x = 0.5, y = 1.5, just = c("left", "top"),
#'     default.units = "inches"
#' )
#' annoGenomeLabel(
#'     plot = genesPlot, scale = "bp",
#'     x = 0.5, y = 1.75, just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @export
annoGenomeLabel <- function(plot, fontsize = 10, fontcolor = "black",
                            linecolor = "black", margin = unit(1, "mm"),
                            scale = "bp", commas = TRUE, sequence = TRUE,
                            boxWidth = 0.5, axis = "x", at = NULL,
                            tcl = 0.5, x, y, just = c("left", "top"),
                            default.units = "inches", params = NULL, ...) {


    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================
    
    genomeLabelInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "genomeLabelInternal"
    )
    
    # =========================================================================
    # CATCH ARGUMENT/PLOT INPUT ERRORS
    # =========================================================================
    if (is.null(genomeLabelInternal$plot)) {
        stop("argument \"plot\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(genomeLabelInternal$x)) {
        stop("argument \"x\" is missing, with no default.", call. = FALSE)
    }
    if (is.null(genomeLabelInternal$y)) {
        stop("argument \"y\" is missing, with no default.", call. = FALSE)
    }
    

    ## Check that input plot is a valid type of plot to be annotated
    if (!is(genomeLabelInternal$plot, "manhattan")) {

        ## Manhattan plots can do whole genome assembly but other plots can't
        inputNames <- attributes(genomeLabelInternal$plot)$names
        if (!("chrom" %in% inputNames) |
            !("chromstart" %in% inputNames) |
            !("chromend" %in% inputNames)) {
            stop("Invalid input plot. Please input a plot that has genomic ",
            "coordinates associated with it.", call. = FALSE)
        }
    }

    # =========================================================================
    # ASSIGN PARAMETERS BASED ON PLOT INPUT
    # =========================================================================

    chrom <- genomeLabelInternal$plot$chrom
    chromstart <- genomeLabelInternal$plot$chromstart
    chromend <- genomeLabelInternal$plot$chromend
    assembly <- genomeLabelInternal$plot$assembly
    x <- genomeLabelInternal$x
    y <- genomeLabelInternal$y
    length <- genomeLabelInternal$plot$width
    if (genomeLabelInternal$axis == "y") {
        length <- genomeLabelInternal$plot$height
    }
    
    just <- justConversion(just = genomeLabelInternal$just)
    
    ## Whole genome Manhattan plot spacing
    space <- genomeLabelInternal$plot$space

    # =========================================================================
    # CALL PLOTGENOMELABEL
    # =========================================================================

    genomeLabel <- bbPlotGenomeLabel(
        chrom = chrom, chromstart = chromstart,
        chromend = chromend, assembly = assembly,
        fontsize = genomeLabelInternal$fontsize,
        fontcolor = genomeLabelInternal$fontcolor,
        linecolor = genomeLabelInternal$linecolor,
        margin = genomeLabelInternal$margin,
        scale = genomeLabelInternal$scale,
        commas = genomeLabelInternal$commas,
        sequence = genomeLabelInternal$sequence,
        boxWidth = genomeLabelInternal$boxWidth,
        axis = genomeLabelInternal$axis,
        at = genomeLabelInternal$at,
        tcl = genomeLabelInternal$tcl,
        x = x, y = y, length = length,
        just = just,
        default.units = genomeLabelInternal$default.units,
        space = space, params = params, ...
    )

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    invisible(genomeLabel)
}
