#' Plot genomic coordinates along the x or y-axis of a BentoBox plot
#' 
#' @usage bb_plotGenomeLabel(
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
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
#'     length,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     params = NULL,
#'     ...
#' )
#'
#' @param chrom Chromosome of genome label, as a string,
#' or a character vector of chromosomes for a whole genome Manhattan plot.
#' @param chromstart Integer start of genome label.
#' @param chromend Integer end of genome label.
#' @param assembly Default genome assembly as a string or a
#' \link[BentoBox]{bb_assembly} object.
#' @param fontsize A numeric specifying text fontsize in points.
#' Default value is \code{fontsize = 10}.
#' @param fontcolor A character value indicating the color for text.
#' Default value is \code{fontcolor = "black"}.
#' @param linecolor A character value indicating the color of
#' the genome label axis. Default value is \code{linecolor = "black"}.
#' @param margin A numeric or unit vector specifying space between axis
#' and coordinate labels. Default value is \code{margin = unit(1, "mm")},
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
#' \code{bb_plotHicSquare}.}
#' }
#' @param at A numeric vector of x-value locations for tick marks.
#' @param tcl A numeric specifying the length of tickmarks as a
#' fraction of text height. Default value is \code{tcl = 0.5}.
#' @param x A numeric or unit object specifying genome label x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying genome label y-location.
#' The character value will
#' place the genome label y relative to the bottom of the most recently
#' plotted BentoBox plot according to the units of the BentoBox page.
#' @param length A numeric or unit object specifying length of
#' genome label axis.
#' @param just Justification of genome label relative to its (x, y)
#' location. If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"},
#' and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to
#' use if \code{x}, \code{y}, or \code{length} are only given as numerics.
#' Default value is \code{default.units = "inches"}.
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters or digit specifications.
#' See \link[grid]{gpar} and \link[base]{formatC}.
#'
#' @return Returns a \code{bb_genomeLabel} object containing
#' relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load hg19 genomic annotation packages
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#' library("BSgenome.Hsapiens.UCSC.hg19")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 5, height = 3, default.units = "inches")
#'
#' ## Plot and place gene track on a BentoBox page
#' genesPlot <- bb_plotGenes(
#'     chrom = "chr8",
#'     chromstart = 1000000, chromend = 2000000,
#'     assembly = "hg19", fill = c("grey", "grey"),
#'     fontcolor = c("grey", "grey"),
#'     x = 0.5, y = 0.25, width = 4, height = 1,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Plot x-axis genome labels at different scales
#' bb_plotGenomeLabel(
#'     chrom = "chr8",
#'     chromstart = 1000000, chromend = 2000000,
#'     assembly = "hg19",
#'     scale = "Mb",
#'     x = 0.5, y = 1.25, length = 4, just = c("left", "top"),
#'     default.units = "inches"
#' )
#' bb_plotGenomeLabel(
#'     chrom = "chr8",
#'     chromstart = 1000000, chromend = 2000000,
#'     assembly = "hg19",
#'     scale = "Kb",
#'     x = 0.5, y = 1.5, length = 4, just = c("left", "top"),
#'     default.units = "inches"
#' )
#' bb_plotGenomeLabel(
#'     chrom = "chr8",
#'     chromstart = 1000000, chromend = 2000000,
#'     assembly = "hg19",
#'     scale = "bp",
#'     x = 0.5, y = 1.75, length = 4, just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Plot a different genomic label region, zooming in enough
#' ## to see base pairs
#' bb_plotGenomeLabel(
#'     chrom = "chr8",
#'     chromstart = 1000000, chromend = 1000050,
#'     assembly = "hg19",
#'     x = 0.25, y = 2.2, length = 4.5
#' )
#' bb_plotGenomeLabel(
#'     chrom = "chr8",
#'     chromstart = 1000000, chromend = 1000020,
#'     assembly = "hg19",
#'     x = 0, y = 2.6, length = 5
#' )
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#' @export
bb_plotGenomeLabel <- function(chrom, chromstart = NULL, chromend = NULL,
                            assembly = "hg38", fontsize = 10,
                            fontcolor = "black", linecolor = "black",
                            margin = unit(1, "mm"),
                            scale = "bp", commas = TRUE, sequence = TRUE,
                            boxWidth = 0.5, axis = "x", at = NULL,
                            tcl = 0.5, x, y, length,
                            just = c("left", "top"),
                            default.units = "inches", params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors for bb_plotGenomeLabel
    errorcheck_bb_genomeLabel <- function(scale, ticks, object, axis) {

        ## Check that scale is an appropriate value
        if (!scale %in% c("bp", "Kb", "Mb")) {
            stop("Invalid \'scale\'. Options are \'bp\', \'Kb\', or \'Mb\'.",
                call. = FALSE
            )
        }

        if (!is.null(ticks)) {

            ## Can't have ticks if label is genome assembly for manhattan plot
            if (is.null(object$chrom)) {
                stop("Cannot add tick marks to a genome label of entire ",
                    "genome assembly.", call. = FALSE)
            }

            ## Make sure ticks fall within the chromstart to chromend range
            if (range(ticks)[1] < object$chromstart |
                range(ticks)[2] > object$chromend) {
                stop("Given tick locations do not fall within ",
                    "the genomic range.",
                    call. = FALSE
                )
            }
        }

        if (!axis %in% c("x", "y")) {
            stop("Invalid \'axis\'. Options are \'x\' or \'y\'.", call. = FALSE)
        }
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_genomeLabelInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_genomeLabelInternal"
    )

    ## Parsing for "space" from input Manhattan plot from bb_annoGenomeLabel
    additionalParams <- list(...)
    if ("space" %in% names(additionalParams)) {
        bb_genomeLabelInternal$space <- additionalParams$space
    }

    ## Assign "gp"
    bb_genomeLabelInternal$gp <- gpar(
        fontsize = bb_genomeLabelInternal$fontsize,
        col = bb_genomeLabelInternal$linecolor,
        fontcolor = bb_genomeLabelInternal$fontcolor,
        lineend = "butt"
    )
    bb_genomeLabelInternal$gp <- setGP(
        gpList = bb_genomeLabelInternal$gp,
        params = bb_genomeLabelInternal, ...
    )
    
    ## Justification
    bb_genomeLabelInternal$just <- 
        bb_justConversion(bb_genomeLabelInternal$just)
    
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    bb_genomeLabel <- structure(list(
        chrom = bb_genomeLabelInternal$chrom,
        chromstart = bb_genomeLabelInternal$chromstart,
        chromend = bb_genomeLabelInternal$chromend,
        assembly = bb_genomeLabelInternal$assembly,
        x = bb_genomeLabelInternal$x,
        y = bb_genomeLabelInternal$y,
        width = NULL, height = NULL,
        just = bb_genomeLabelInternal$just,
        grobs = NULL
    ), class = "bb_genomeLabel")

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================
    
    if (is.null(bb_genomeLabel$x)) stop("argument \"x\" is missing, ",
                                        "with no default.", call. = FALSE)
    
    if (is.null(bb_genomeLabel$y)) stop("argument \"y\" is missing, ",
                                        "with no default.", call. = FALSE)
    
    if (is.null(bb_genomeLabel$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)
    
    if (is.null(bb_genomeLabelInternal$length)) {
        stop("argument \"length\" is missing, ",
            "with no default.",
            call. = FALSE
        )
    }
    
    if (base::length(bb_genomeLabel$chrom) == 1) {
        if (is.null(bb_genomeLabel$chromstart)) stop("argument \"chromstart\" ",
                                                    "is missing, with no ",
                                                    "default.", call. = FALSE)

        if (is.null(bb_genomeLabel$chromend)) stop("argument \"chromend\" ",
                                                "is missing, with no ",
                                                "default.", call. = FALSE)
    } else {
        
        if (is.null(bb_genomeLabelInternal$space)) {
            bb_genomeLabelInternal$space <- 0.01
        }
    }

    check_bbpage(error = "Cannot plot a genome label without a BentoBox page.")
    errorcheck_bb_genomeLabel(
        scale = bb_genomeLabelInternal$scale,
        ticks = bb_genomeLabelInternal$at,
        object = bb_genomeLabel,
        axis = bb_genomeLabelInternal$axis
    )
    
    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    bb_genomeLabel$assembly <-
        parse_bbAssembly(assembly = bb_genomeLabel$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================
    
    bb_genomeLabel$x <- misc_defaultUnits(value = bb_genomeLabel$x,
                                        name = "x",
                                        default.units = 
                                        bb_genomeLabelInternal$default.units)
    bb_genomeLabel$y <- misc_defaultUnits(
        value = bb_genomeLabel$y,
        name = "y",
        default.units = bb_genomeLabelInternal$default.units)
    bb_genomeLabelInternal$length <- misc_defaultUnits(
        value = bb_genomeLabelInternal$length,
        name = "length",
        default.units = bb_genomeLabelInternal$default.units)
    bb_genomeLabelInternal$margin <- misc_defaultUnits(
        value = bb_genomeLabelInternal$margin,
        name = "margin",
        default.units = bb_genomeLabelInternal$default.units)
    
    # =========================================================================
    # LABEL DEPTH
    # =========================================================================

    bb_genomeLabelInternal$margin <- 
        convertHeight(bb_genomeLabelInternal$margin,
                    unitTo = get("page_units", envir = bbEnv)
        )
    bb_genomeLabelInternal$tgH <- convertHeight(heightDetails(textGrob(
        label = bb_genomeLabelInternal$scale,
        x = 0.5, y = 0.5,
        default.units = "npc",
        gp = bb_genomeLabelInternal$gp
    )),
    unitTo = get("page_units", envir = bbEnv)
    )
    
    if (!is.null(bb_genomeLabelInternal$at)) {
        bb_genomeLabelInternal$tick_height <- 
            bb_genomeLabelInternal$tgH * (bb_genomeLabelInternal$tcl)
        bb_genomeLabelInternal$depth <- 
            convertHeight(bb_genomeLabelInternal$tgH + 
                        bb_genomeLabelInternal$tick_height + 
                        0.5 * bb_genomeLabelInternal$tgH + 
                        bb_genomeLabelInternal$margin,
            unitTo = get("page_units", envir = bbEnv),
            valueOnly = TRUE
        )
    } else {
        bb_genomeLabelInternal$tick_height <- NULL
        bb_genomeLabelInternal$depth <- 
            convertHeight(bb_genomeLabelInternal$tgH + 
                        bb_genomeLabelInternal$margin,
            unitTo = get("page_units", envir = bbEnv),
            valueOnly = TRUE
        )
    }

    # =========================================================================
    # GENOME LABEL FOR SINGLE CHROMOSOME
    # =========================================================================
    
    if (base::length(bb_genomeLabel$chrom) == 1){
        
        vp_name <- bb_plotChromGenomeLabel(bb_genomeLabel = bb_genomeLabel,
                                bb_genomeLabelInternal = bb_genomeLabelInternal,
                                ...)
        
    } else {
    # =========================================================================
    # GENOME LABEL FOR MANHATTAN PLOT
    # =========================================================================
        vp_name <- bb_plotManhattanGenomeLabel(bb_genomeLabel = bb_genomeLabel,
                                bb_genomeLabelInternal = bb_genomeLabelInternal)
        
    }

    # =========================================================================
    # ASSIGN GROBS TO OBJECT
    # =========================================================================

    bb_genomeLabel$grobs <- get("genomeLabel_grobs", envir = bbEnv)
    grid.draw(bb_genomeLabel$grobs)

    # =========================================================================
    # ASSIGN DIMENSIONS BASED ON AXIS
    # =========================================================================

    bb_genomeLabel$width <- bb_genomeLabelInternal$length
    bb_genomeLabel$height <- bb_genomeLabelInternal$depth
    if (bb_genomeLabelInternal$axis == "y") {
        bb_genomeLabel$width <- bb_genomeLabelInternal$depth
        bb_genomeLabel$height <- bb_genomeLabelInternal$length
    }

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_genomeLabel[", vp_name, "]")
    invisible(bb_genomeLabel)
}

# Plots a standard genome label for a single chromosomal region
# @param bb_genomeLabel bb_genomeLabel object from bb_plotGenomeLabel
# @param bb_genomeLabelInternal bb_genomeLabelInternal object 
# from bb_plotGenomeLabel
bb_plotChromGenomeLabel <- function(bb_genomeLabel, 
                                    bb_genomeLabelInternal, ...){
    
    # =========================================================================
    # FUNCTIONS
    # =========================================================================
    
    ## Define a function that adds commas to chromstart/chromend labels
    comma_labels <- function(object, commas, fact, ...) {
        digits <- list(...)$digits
        if (is.null(digits)){
            if (fact == 1){
                digits <- rep(0, 2)
            } else {
                digits <- rep(1, 2)
            }
        }
        if (length(digits) == 1){
            digits <- rep(digits, 2)
        }
        
        roundedStart <- round(object$chromstart / fact, digits[1])
        if ((roundedStart * fact) != object$chromstart) {
            roundedStart <- paste0("~", roundedStart)
            warning("Start label is rounded.", call. = FALSE)
        }
        
        roundedEnd <- round(object$chromend / fact, digits[2])
        if ((roundedEnd * fact) != object$chromend) {
            roundedEnd <- paste0("~", roundedEnd)
            warning("End label is rounded.", call. = FALSE)
        }

        if (commas == TRUE) {
            chromstartlabel <- formatC(roundedStart,
                                    format = "f",
                                    big.mark = ",", digits = digits[1]
            )
            chromendlabel <- formatC(roundedEnd,
                                    format = "f",
                                    big.mark = ",", digits = digits[2]
            )
        } else {
            chromstartlabel <- roundedStart
            chromendlabel <- roundedEnd
        }
        
        return(list(chromstartlabel, chromendlabel))
    }
    
    ## Define a function that makes the label viewport
    chrom_viewport <- function(object, length, depth, seqType,
                            vp_name, seqHeight, just, axis){
        
        ## No matter the orientation, convert to page units
        convertedPageCoords <- convert_page(object = structure(list(
            width = length,
            height = unit(
                depth,
                get("page_units", envir = bbEnv)
            ),
            x = object$x,
            y = object$y
        ),
        class = "bb_genomeLabelInternal"
        ))
        ## Add "length" and "depth" into converted dimensions
        convertedPageCoords$length <- convertedPageCoords$width
        convertedPageCoords$depth <- convertedPageCoords$height
        
        ## Compile new dimensions into a new dummy viewport,
        ## where the default is along the x-axis
        convertedViewport <- viewport(
            width = convertedPageCoords$length,
            height = convertedPageCoords$depth,
            x = convertedPageCoords$x,
            y = convertedPageCoords$y, just = just
        )
        
        if (!is.null(seqType)) {
            
            ## Get x and y coordinates of top left of what
            ## would be the entire viewport
            topLeftViewport <- vp_topLeft(viewport = convertedViewport)
            seq_height <- unit(seqHeight, get("page_units", envir = bbEnv))
            ## One vp for genome
            vp1 <- viewport(
                width = convertedPageCoords$width,
                height = unit(depth, get("page_units", envir = bbEnv)),
                x = topLeftViewport[[1]],
                y = topLeftViewport[[2]] - seq_height,
                just = c("left", "top"),
                name = paste0(vp_name, "_01"),
                xscale = c(object$chromstart, object$chromend),
                yscale = c(0, depth)
            )
            ## One vp for sequence
            vp2 <- viewport(
                width = convertedPageCoords$width,
                height = seq_height,
                x = topLeftViewport[[1]],
                y = topLeftViewport[[2]],
                just = c("left", "top"),
                name = paste0(vp_name, "_02"),
                clip = "on",
                xscale = c(object$chromstart, object$chromend)
            )
            
            ## Combine viewports into one
            vp <- vpList(vp1, vp2)
        } else {
            if (axis == "y") {
                
                
                ## Update converted viewport for y-axis
                convertedViewport <- viewport(
                    width = convertedPageCoords$depth,
                    height = convertedPageCoords$length,
                    x = convertedPageCoords$x,
                    y = convertedPageCoords$y,
                    just = just
                )
                ## Get x and y coords of bottom right to rotate vp
                bottomRightViewport <-
                    vp_bottomRight(viewport = convertedViewport)
                
                ## Make x-axis equivalent viewport and rotate into
                ## dimensions of given y-axis viewport
                vp <- viewport(
                    width = convertedPageCoords$length,
                    height = convertedPageCoords$depth,
                    x = bottomRightViewport[[1]] -
                        convertedPageCoords$depth,
                    y = bottomRightViewport[[2]],
                    just = c("left", "top"),
                    name = vp_name,
                    xscale = c(object$chromstart, object$chromend),
                    yscale = c(0, depth),
                    angle = 90
                )
            } else {
                vp <- viewport(
                    width = convertedPageCoords$width,
                    height = convertedPageCoords$height,
                    x = convertedPageCoords$x,
                    y = convertedPageCoords$y,
                    just = just,
                    name = vp_name,
                    xscale = c(object$chromstart, object$chromend),
                    yscale = c(0, depth)
                )
            }
        }
        
        return(vp)
    }
    
    ## Define a function that makes tick, line, and text grobs for
    ## chrom/chromstart/chromend labels
    chrom_grobs <- function(ticks, seqType, scale, chromLabel,
                            startLabel, endLabel, object, vp,
                            yaxis) {
        margin <- convertHeight(object$margin,
                                unitTo = get("page_units", envir = bbEnv),
                                valueOnly = TRUE
        )
        
        height <- convertHeight(object$depth,
                                unitTo = get("page_units", envir = bbEnv),
                                valueOnly = TRUE
        )
        
        if (!is.null(seqType)) {
            assign("genomeLabel_grobs", gTree(), envir = bbEnv)
            chrom_vp <- vp[[1]]
            
            if (!is.null(ticks)) {
                tgH <- convertHeight(object$tgH,
                                    unitTo = get("page_units", envir = bbEnv),
                                    valueOnly = TRUE
                )
                tick_height <- convertHeight(object$tick_height,
                                            unitTo = get("page_units", 
                                                        envir = bbEnv),
                                            valueOnly = TRUE
                )
                x_coords <- ticks
                y0_coord <- height
                y1_coords <- rep(height - tick_height, length(ticks))
                yLabel <- unit(height - (tick_height + margin), "native")
                
                tickGrobs <- segmentsGrob(
                    x0 = x_coords,
                    y0 = rep(y0_coord, length(ticks)),
                    x1 = x_coords,
                    y1 = y1_coords,
                    vp = chrom_vp,
                    gp = object$gp,
                    default.units = "native"
                )
                line <- segmentsGrob(
                    x0 = unit(0, "npc"),
                    x1 = unit(1, "npc"),
                    y0 = y0_coord,
                    y1 = y0_coord,
                    vp = chrom_vp,
                    gp = object$gp,
                    default.units = "native"
                )
                object$gp$col <- object$gp$fontcolor
                startLab <- textGrob(
                    label = paste(startLabel, scale),
                    x = unit(0, "npc"), y = yLabel,
                    just = c("left", "top"),
                    vp = chrom_vp,
                    gp = object$gp
                )
                endLab <- textGrob(
                    label = paste(endLabel, scale),
                    x = unit(1, "npc"), y = yLabel,
                    just = c("right", "top"),
                    vp = chrom_vp,
                    gp = object$gp
                )
                object$gp$fontface <- "bold"
                chromLab <- textGrob(
                    label = chromLabel,
                    x = unit(0.5, "npc"), y = yLabel,
                    vp = chrom_vp,
                    gp = object$gp, just = c("center", "top")
                )
                
                assign("genomeLabel_grobs",
                    setChildren(get("genomeLabel_grobs", envir = bbEnv),
                                    children = gList(
                                                line, chromLab, startLab,
                                                endLab, tickGrobs
                                        )
                        ),
                        envir = bbEnv
                )
            } else {
                yLabel <- unit(height - margin, "native")
                line <- segmentsGrob(
                    x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                    y0 = height, y1 = height,
                    vp = chrom_vp,
                    gp = object$gp, default.units = "native"
                )
                object$gp$col <- object$gp$fontcolor
                startLab <- textGrob(
                    label = paste(startLabel, scale, sep = " "),
                    x = unit(0, "npc"), y = yLabel,
                    vp = chrom_vp,
                    just = c("left", "top"),
                    gp = object$gp
                )
                endLab <- textGrob(
                    label = paste(endLabel, scale, sep = " "),
                    x = unit(1, "npc"), y = yLabel,
                    vp = chrom_vp,
                    just = c("right", "top"),
                    gp = object$gp
                )
                object$gp$fontface <- "bold"
                chromLab <- textGrob(
                    label = chromLabel, x = unit(0.5, "npc"),
                    y = yLabel,
                    vp = chrom_vp,
                    gp = object$gp,
                    just = c("center", "top")
                )
                
                assign("genomeLabel_grobs",
                    setChildren(get("genomeLabel_grobs", envir = bbEnv),
                            children = gList(
                                        line, chromLab,
                                        startLab, endLab
                                        )
                    ),
                    envir = bbEnv
                )
            }
        } else {
            assign("genomeLabel_grobs", gTree(vp = vp), envir = bbEnv)
            
            if (!is.null(ticks)) {
                tgH <- convertHeight(object$tgH,
                                    unitTo = get("page_units", envir = bbEnv),
                                    valueOnly = TRUE
                )
                tick_height <- convertHeight(object$tick_height,
                                            unitTo = get("page_units", 
                                                        envir = bbEnv),
                                            valueOnly = TRUE
                )
                x_coords <- ticks
                
                if (yaxis == TRUE) {
                    y0_coord <- 0
                    y1_coords <- rep(tick_height, length(ticks))
                    yLabel <- unit(tick_height + margin, "native")
                    
                    tickGrobs <- segmentsGrob(
                        x0 = x_coords,
                        y0 = rep(y0_coord, length(ticks)),
                        x1 = x_coords, y1 = y1_coords,
                        gp = object$gp, default.units = "native"
                    )
                    line <- segmentsGrob(
                        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                        y0 = y0_coord, y1 = y0_coord,
                        gp = object$gp, default.units = "native"
                    )
                    object$gp$col <- object$gp$fontcolor
                    startLab <- textGrob(
                        label = paste(startLabel, scale),
                        x = unit(0, "npc"), y = yLabel,
                        just = c("left", "bottom"),
                        gp = object$gp
                    )
                    endLab <- textGrob(
                        label = paste(endLabel, scale),
                        x = unit(1, "npc"), y = yLabel,
                        just = c("right", "bottom"),
                        gp = object$gp
                    )
                    object$gp$fontface <- "bold"
                    chromLab <- textGrob(
                        label = chromLabel,
                        x = unit(0.5, "npc"), y = yLabel,
                        gp = object$gp, just = c("center", "bottom")
                    )
                } else {
                    y0_coord <- height
                    y1_coords <- rep(height - tick_height, length(ticks))
                    yLabel <- unit(height - (tick_height + margin), "native")
                    tickGrobs <- segmentsGrob(
                        x0 = x_coords,
                        y0 = rep(y0_coord, length(ticks)),
                        x1 = x_coords, y1 = y1_coords,
                        gp = object$gp, default.units = "native"
                    )
                    line <- segmentsGrob(
                        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                        y0 = y0_coord, y1 = y0_coord,
                        gp = object$gp, default.units = "native"
                    )
                    object$gp$col <- object$gp$fontcolor
                    startLab <- textGrob(
                        label = paste(startLabel, scale),
                        x = unit(0, "npc"),
                        y = unit(tgH + 0.25 * tgH, "native"),
                        just = c("left", "top"),
                        gp = object$gp
                    )
                    endLab <- textGrob(
                        label = paste(endLabel, scale),
                        x = unit(1, "npc"),
                        y = unit(tgH + 0.25 * tgH, "native"),
                        just = c("right", "top"),
                        gp = object$gp
                    )
                    object$gp$fontface <- "bold"
                    chromLab <- textGrob(
                        label = chromLabel, x = unit(0.5, "npc"),
                        y = unit(tgH + 0.25 * tgH, "native"),
                        gp = object$gp, just = c("center", "top")
                    )
                }
                
                assign("genomeLabel_grobs",
                    setChildren(get("genomeLabel_grobs", envir = bbEnv),
                                children = gList(
                                            line, chromLab, startLab,
                                            endLab, tickGrobs
                                        )
                        ),
                    envir = bbEnv
                )
            } else {
                if (yaxis == TRUE) {
                    yLabel <- unit(margin, "native")
                    line <- segmentsGrob(
                        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                        y0 = unit(0, "npc"), y1 = unit(0, "npc"),
                        gp = object$gp
                    )
                    object$gp$col <- object$gp$fontcolor
                    startLab <- textGrob(
                        label = paste(startLabel, scale, sep = " "),
                        x = unit(0, "npc"), y = yLabel,
                        just = c("left", "bottom"),
                        gp = object$gp
                    )
                    endLab <- textGrob(
                        label = paste(endLabel, scale, sep = " "),
                        x = unit(1, "npc"), y = yLabel,
                        just = c("right", "bottom"),
                        gp = object$gp
                    )
                    object$gp$fontface <- "bold"
                    chromLab <- textGrob(
                        label = chromLabel, x = unit(0.5, "npc"),
                        y = yLabel,
                        gp = object$gp,
                        just = c("center", "bottom")
                    )
                } else {
                    yLabel <- unit(height - margin, "native")
                    line <- segmentsGrob(
                        x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                        y0 = height, y1 = height,
                        gp = object$gp, default.units = "native"
                    )
                    object$gp$col <- object$gp$fontcolor
                    startLab <- textGrob(
                        label = paste(startLabel, scale, sep = " "),
                        x = unit(0, "npc"), y = yLabel,
                        just = c("left", "top"),
                        gp = object$gp
                    )
                    endLab <- textGrob(
                        label = paste(endLabel, scale, sep = " "),
                        x = unit(1, "npc"), y = yLabel,
                        just = c("right", "top"),
                        gp = object$gp
                    )
                    object$gp$fontface <- "bold"
                    chromLab <- textGrob(
                        label = chromLabel, x = unit(0.5, "npc"),
                        y = yLabel,
                        gp = object$gp,
                        just = c("center", "top")
                    )
                }
                
                
                
                assign("genomeLabel_grobs",
                    setChildren(get("genomeLabel_grobs", envir = bbEnv),
                                children = gList(
                                            line, chromLab,
                                            startLab, endLab
                                )
                        ),
                    envir = bbEnv
                )
            }
        }
    }
    
    ## Define a function that makes sequence grobs (boxes or letters)
    seq_grobs <- function(object, seqHeight, seqType, assembly, chromLabel, vp,
                        boxWidth, gparParams) {
        
        bsgenome <- eval(parse(text = paste0(as.name(object$assembly$BSgenome),
                                            "::",
                                            as.name(object$assembly$BSgenome))))
        ## Get sequence in that region
        sequence <- strsplit(as.character(BSgenome::getSeq(
            bsgenome,
            GenomicRanges::GRanges(
                seqnames = chromLabel,
                ranges = IRanges::IRanges(
                    start = object$chromstart,
                    end = object$chromend
                )
            )
        )),
        split = ""
        )
        ## Make dataframe of sequence letter, position, and color
        dfSequence <- data.frame(
            "nucleotide" = unlist(sequence),
            "pos" = seq(object$chromstart, object$chromend),
            "col" = "grey"
        )
        
        ## Make colors A = green, T = red, G = orange, C = blue
        invisible(tryCatch(dfSequence[which(
            dfSequence$nucleotide == "A"
        ), ]$col <- "#7CD95B",
        error = function(e) {}
        ))
        invisible(tryCatch(dfSequence[which(
            dfSequence$nucleotide == "T"
        ), ]$col <- "#F1686C",
        error = function(e) {}
        ))
        invisible(tryCatch(dfSequence[which(
            dfSequence$nucleotide == "G"
        ), ]$col <- "#FFDD12",
        error = function(e) {}
        ))
        invisible(tryCatch(dfSequence[which(
            dfSequence$nucleotide == "C"
        ), ]$col <- "#5566FF",
        error = function(e) {}
        ))
        
        seq_vp <- vp[[2]]
        
        ## Make grobs based on seqType
        if (seqType == "letters") {
            seqGrobs <- textGrob(
                label = dfSequence$nucleotide, x = dfSequence$pos,
                y = unit(0.5, "npc"), just = "center",
                vp = seq_vp,
                default.units = "native",
                gp = gpar(
                    col = dfSequence$col,
                    fontsize = gparParams$gp$fontsize - 2
                )
            )
        } else if (seqType == "boxes") {
            seqGrobs <- rectGrob(
                x = dfSequence$pos, y = unit(1, "npc"),
                width = boxWidth,
                height = unit(
                    seqHeight - 0.05 * seqHeight,
                    get("page_units", envir = bbEnv)
                ),
                just = c("center", "top"),
                vp = seq_vp,
                default.units = "native",
                gp = gpar(col = NA, fill = dfSequence$col)
            )
        }
        
        assign("genomeLabel_grobs",
                addGrob(
                    gTree = get("genomeLabel_grobs", envir = bbEnv),
                    child = seqGrobs
                ),
                envir = bbEnv
        )
    }
    
    # =========================================================================
    # COMMAS
    # =========================================================================
    
    ## Determine scale of labels
    if (bb_genomeLabelInternal$scale == "bp") {
        fact <- 1
    }
    if (bb_genomeLabelInternal$scale == "Mb") {
        fact <- 1000000
    }
    if (bb_genomeLabelInternal$scale == "Kb") {
        fact <- 1000
    }
    
    commaLabels <- comma_labels(
            object = bb_genomeLabel,
            commas = bb_genomeLabelInternal$commas,
            fact = fact, ...
        )
    chromstartlabel <- commaLabels[[1]]
    chromendlabel <- commaLabels[[2]]
    
    # =========================================================================
    # NUCLEOTIDE SEQUENCE INFORMATION
    # =========================================================================
    
    seq_height <- heightDetails(textGrob(
        label = "A",
        x = 0.5, y = 0.5,
        default.units = "npc",
        gp = gpar(fontsize = bb_genomeLabelInternal$gp$fontsize - 2)
    ))
    seq_height <- convertHeight(seq_height + 0.05 * seq_height,
                                unitTo = get("page_units", envir = bbEnv)
    )
    
    seqType <- NULL
    if (bb_genomeLabelInternal$sequence == TRUE) {
        if (bb_genomeLabelInternal$axis == "x") {
            labelWidth <- convertWidth(bb_genomeLabelInternal$length,
                                    unitTo = "inches",
                                    valueOnly = TRUE
            )
            bpWidth <- convertWidth(widthDetails(textGrob(
                label = "A",
                x = 0.5, y = 0.5,
                default.units = "npc",
                gp = gpar(fontsize = bb_genomeLabelInternal$gp$fontsize - 2)
            )),
            unitTo = "inches",
            valueOnly = TRUE
            )
            seqRange <- bb_genomeLabel$chromend - bb_genomeLabel$chromstart
            seqWidth <- bpWidth * seqRange
            
            if (seqWidth <= labelWidth) {
                seqType <- "letters"
            } else if (seqWidth / labelWidth <= 9) {
                seqType <- "boxes"
            }
        }
    }
    
    ## Check for BSgenome packages and reset seqType if necessary
    if (!is.null(seqType)) {
        if (!is.null(bb_genomeLabel$assembly$BSgenome)) {
            
            if (!requireNamespace(bb_genomeLabel$assembly$BSgenome, 
                                quietly = TRUE)){
                bsChecks <- FALSE
                warning("`", bb_genomeLabel$assembly$BSgenome, 
                        "` not available. ",
                        "Sequence information will not be displayed.", 
                        call. = FALSE)
            } else {
                bsChecks <- TRUE
            }
            
            if (bsChecks == FALSE) {
                seqType <- NULL
            }
        } else {
            warning("No `BSgenome` package found for the input assembly. ",
                    "Sequence information cannot be displayed.", call. = FALSE)
            seqType <- NULL
        }
    }
    
    if (!is.null(seqType)) {
        seq_height <- convertHeight(seq_height,
                                    unitTo = get("page_units", envir = bbEnv),
                                    valueOnly = TRUE
        )
        bb_genomeLabelInternal$depth <- unit(
            bb_genomeLabelInternal$depth + seq_height,
            get("page_units", envir = bbEnv)
        )
    } else {
        bb_genomeLabelInternal$depth <- unit(
            bb_genomeLabelInternal$depth,
            get("page_units", envir = bbEnv)
        )
    }
    
    # =========================================================================
    # VIEWPORTS
    # =========================================================================
    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_genomeLabel",
        base::length(grep(
            pattern = "bb_genomeLabel",
            x = currentViewports
        )) + 1
    )
    
    ## Make viewport
    vp <- chrom_viewport(
        object = bb_genomeLabel,
        length = bb_genomeLabelInternal$length,
        depth = bb_genomeLabelInternal$depth,
        seqType = seqType,
        seqHeight = seq_height,
        vp_name = vp_name, just = bb_genomeLabel$just,
        axis = bb_genomeLabelInternal$axis
    )
    
    # =========================================================================
    # GROBS AND GTREE
    # =========================================================================
    
    chrom_grobs(
        ticks = bb_genomeLabelInternal$at,
        seqType = seqType,
        scale = bb_genomeLabelInternal$scale,
        chromLabel = bb_genomeLabel$chrom, 
        startLabel = chromstartlabel, endLabel = chromendlabel,
        object = bb_genomeLabelInternal, vp = vp,
        yaxis = (bb_genomeLabelInternal$axis == "y")
    )
    
    ## Sequence grobs if applicable
    if (!is.null(seqType)) {
        seq_grobs(
            object = bb_genomeLabel, seqHeight = seq_height,
            seqType = seqType, assembly = bb_genomeLabel$assembly,
            chromLabel = bb_genomeLabel$chrom, vp = vp,
            boxWidth = bb_genomeLabelInternal$boxWidth,
            gparParams = bb_genomeLabelInternal
        )
    }
    
    return(vp_name)
}

# Plots a genome label for a multi-chromosomal Manhattan plot
# @param bb_genomeLabel bb_genomeLabel object from bb_plotGenomeLabel
# @param bb_genomeLabelInternal bb_genomeLabelInternal object 
# from bb_plotGenomeLabel
bb_plotManhattanGenomeLabel <- function(bb_genomeLabel, bb_genomeLabelInternal){
    
    # =========================================================================
    # FUNCTIONS
    # =========================================================================
    
    ## Define a function that makes the label viewport
    manhattan_viewport <- function(object, length, depth,
                                vp_name, just, axis, space){
        
        ## Convert to page units
        convertedPageCoords <- convert_page(object = structure(list(
            width = length,
            height = unit(
                depth,
                get("page_units", envir = bbEnv)
            ),
            x = object$x,
            y = object$y
        ),
        class = "bb_genomeLabelInternal"
        ))
        ## Add "length" and "depth" into converted dimensions
        convertedPageCoords$length <- convertedPageCoords$width
        convertedPageCoords$depth <- convertedPageCoords$height
        
        ## Compile new dimensions into a new dummy viewport,
        ## where the default is along the x-axis
        convertedViewport <- viewport(
            width = convertedPageCoords$length,
            height = convertedPageCoords$depth,
            x = convertedPageCoords$x,
            y = convertedPageCoords$y, just = just
        )
        
        
        ## Get assembly data
        if (is(object$assembly$TxDb, "TxDb")) {
            txdbChecks <- TRUE
        } else {
            
            if (!requireNamespace(object$assembly$TxDb, quietly = TRUE)){
                txdbChecks <- FALSE
                warning("`", object$assembly$TxDb, "` not available. Please ",
                "install to label genome.", call. = FALSE)
            } else {
                txdbChecks <- TRUE
            }
            
        }
        
        if (txdbChecks == TRUE) {
            if (is(object$assembly$TxDb, "TxDb")) {
                tx_db <- object$assembly$TxDb
            } else {
                tx_db <- eval(parse(text = paste0(as.name(object$assembly$TxDb),
                                        "::",
                                        as.name(object$assembly$TxDb))))
            }
            
            assembly_data <- as.data.frame(setDT(as.data.frame(
                GenomeInfoDb::seqlengths(tx_db)
            ),
            keep.rownames = TRUE
            ))
            colnames(assembly_data) <- c("chrom", "length")
            assembly_data <- assembly_data[which(assembly_data[, "chrom"] %in%
                                                object$chrom), ]
            ## get the offsets based on spacer for the assembly
            offsetAssembly <- spaceChroms(
                assemblyData = assembly_data,
                space = space
            )
            cumsums <- cumsum(as.numeric(assembly_data[, "length"]))
            spacer <- cumsums[length(cumsum(
                as.numeric(assembly_data[, "length"])
            ))] * space
            xscale <- c(0, max(offsetAssembly[, "end"]) + spacer)
        } else {
            xscale <- c(0, 1)
        }
        
        
        if (axis == "y") {
            ## Update converted viewport for y-axis
            convertedViewport <- viewport(
                width = convertedPageCoords$depth,
                height = convertedPageCoords$length,
                x = convertedPageCoords$x,
                y = convertedPageCoords$y, just = just
            )
            ## Get x and y coordinates of bottom right to rotate
            ## x-axis viewport
            bottomRightViewport <-
                vp_bottomRight(viewport = convertedViewport)
            ## Make x-axis equivalent viewport and rotate into
            ## dimensions of given y-axis viewport
            vp <- viewport(
                width = convertedPageCoords$length,
                height = convertedPageCoords$depth,
                x = bottomRightViewport[[1]] - convertedPageCoords$depth,
                y = bottomRightViewport[[2]],
                just = c("left", "top"),
                name = vp_name,
                xscale = c(object$chromstart, object$chromend),
                yscale = c(0, depth),
                angle = 90
            )
        } else {
            vp <- viewport(
                width = convertedPageCoords$width,
                height = convertedPageCoords$height,
                x = convertedPageCoords$x,
                y = convertedPageCoords$y,
                just = just,
                name = vp_name,
                xscale = xscale,
                yscale = c(0, depth)
            )
        }
        
        return(vp)
        
    }
    
    ## Define a function that makes line and text grobs for whole assembly
    ## labels in Manhattan plots
    genome_grobs <- function(object, vp, gp, space) {
        
        ## Initialize gTree
        assign("genomeLabel_grobs", gTree(vp = vp), envir = bbEnv)
        
        ## Get assembly data
        if (is(object$assembly$TxDb, "TxDb")) {
            txdbChecks <- TRUE
        } else {
            
            if (!requireNamespace(object$assembly$TxDb, quietly = TRUE)){
                txdbChecks <- FALSE
            } else {
                txdbChecks <- TRUE
            }
        }
        
        if (txdbChecks == TRUE) {
            if (is(object$assembly$TxDb, "TxDb")) {
                tx_db <- object$assembly$TxDb
            } else {
                tx_db <- eval(parse(text = 
                                        paste0(as.name(object$assembly$TxDb),
                                            "::",
                                            as.name(object$assembly$TxDb))))
            }
            
            assembly_data <- as.data.frame(setDT(as.data.frame(
                GenomeInfoDb::seqlengths(tx_db)
            ),
            keep.rownames = TRUE
            ))
            colnames(assembly_data) <- c("chrom", "length")
            assembly_data <- assembly_data[which(assembly_data[, "chrom"] %in%
                                                    object$chrom), ]
            ## Get the offsets based on spacer for the assembly
            offsetAssembly <- spaceChroms(
                assemblyData = assembly_data,
                space = space
            )
            
            ## Get the centers of each chrom
            chromCenters <- (offsetAssembly[, "start"] + 
                                offsetAssembly[, "end"]) / 2
            
            margin <- convertHeight(object$margin,
                                    unitTo = get("page_units", envir = bbEnv),
                                    valueOnly = TRUE
            )
            
            line <- segmentsGrob(
                x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                y0 = unit(1, "npc"), y1 = unit(1, "npc"), gp = gp
            )
            gp$col <- gp$fontcolor
            labels <- textGrob(
                label = gsub("chr", "", offsetAssembly[, 1]),
                x = chromCenters,
                y = unit(1, "npc") - unit(margin, "native"),
                just = c("center", "top"),
                gp = gp,
                default.units = "native"
            )
            assign("genomeLabel_grobs",
                setChildren(get("genomeLabel_grobs", envir = bbEnv),
                        children = gList(line, labels)
                    ),
                envir = bbEnv
            )
        }
    }
    
    # =========================================================================
    # VIEWPORTS
    # =========================================================================
    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_genomeLabel",
        base::length(grep(
            pattern = "bb_genomeLabel",
            x = currentViewports
        )) + 1
    )
    
    ## Make viewport
    vp <- manhattan_viewport(
        object = bb_genomeLabel,
        length = bb_genomeLabelInternal$length,
        depth = bb_genomeLabelInternal$depth,
        vp_name = vp_name, just = bb_genomeLabel$just,
        axis = bb_genomeLabelInternal$axis,
        space = bb_genomeLabelInternal$space
    )
    
    # =========================================================================
    # GROBS AND GTREE
    # =========================================================================
    
    genome_grobs(
        object = bb_genomeLabelInternal, vp = vp,
        gp = bb_genomeLabelInternal$gp,
        space = bb_genomeLabelInternal$space
    )
    
    return(vp_name)
}
