#' Plot paired-end genomic range elements
#' 
#' @usage bbPlotPairs(
#'     data,
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
#'     fill = "#1f4297",
#'     linecolor = NA,
#'     bg = NA,
#'     boxHeight = unit(2, "mm"),
#'     spaceWidth = 0.02,
#'     spaceHeight = 0.3,
#'     limitLabel = TRUE,
#'     baseline = FALSE,
#'     baseline.color = "grey",
#'     baseline.lwd = 1,
#'     x = NULL,
#'     y = NULL,
#'     width = NULL,
#'     height = NULL,
#'     just = c("left", "top"),
#'     default.units = "inches",
#'     draw = TRUE,
#'     params = NULL,
#'     ...
#' )
#'
#' @param data A string specifying the BEDPE file path, a dataframe
#' in BEDPE format specifying data to be plotted, or a
#' \link[InteractionSet]{GInteractions} object.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[BentoBox]{bbAssembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param fill A single character value, a vector, or 
#' a \link[BentoBox]{colorby} object specifying fill colors of
#' paired range elements. Default value is \code{fill = "#1f4297"}.
#' @param linecolor A single character value, a vector, or a
#' \link[BentoBox]{colorby} object specifying the color of the lines
#' outlining paired range elements. Default value is \code{linecolor = NA}.
#' Special options include:
#' \itemize{
#' \item{\code{NA}: }{No line color.}
#' \item{\code{"fill"}: }{Same color as \code{fill}.}
#' }
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param boxHeight A numeric or unit object specifying height of boxes
#' at either end of paired range elements.
#' Default value is \code{boxHeight = unit(2, "mm")}.
#' @param spaceWidth A numeric specifying the width of spacing between
#' paired range elements, as a fraction of the plot's genomic range.
#' Default value is \code{spaceWidth = 0.02}.
#' @param spaceHeight A numeric specifying the height of space between
#' boxes of paired range elements on different rows.
#' Default value is \code{spaceHeight = 0.3}.
#' @param limitLabel A logical value indicating whether to draw a "+"
#' when not all elements can be plotted in the plotting space. Default 
#' value is \code{limitLabel = TRUE}.
#' @param baseline Logical value indicating whether to include a baseline
#' along the x-axis. Default value is \code{baseline = FALSE}.
#' @param baseline.color Baseline color.
#' Default value is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width.
#' Default value is \code{baseline.lwd = 1}.
#' @param x A numeric or unit object specifying paired range plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying paired range plot y-location.
#' The character value will
#' place the paired range plot y relative to the bottom of the most recently
#' plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying paired range plot width.
#' @param height A numeric or unit object specifying paired range plot height.
#' @param just Justification of paired range plot relative
#' to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use
#' if \code{x}, \code{y}, \code{width}, or \code{height} are only given
#' as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics
#' output should be produced.
#' @param params An optional \link[BentoBox]{bbParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_pairs} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load paired ranges data in BEDPE format
#' library(BentoBoxData)
#' data("IMR90_DNAloops_pairs")
#'
#' ## Set the coordinates
#' params <- bbParams(
#'     chrom = "chr21",
#'     chromstart = 27900000, chromend = 30700000,
#'     assembly = "hg19",
#'     width = 7
#' )
#'
#' ## Create a page
#' bbPageCreate(width = 7.5, height = 2.1, default.units = "inches")
#'
#' ## Add a length column
#' IMR90_DNAloops_pairs$length <- 
#'         (IMR90_DNAloops_pairs$start2 - IMR90_DNAloops_pairs$start1) / 1000
#'
#' ## Plot the data
#' bedpePlot <- bbPlotPairs(
#'     data = IMR90_DNAloops_pairs, params = params,
#'     fill = colorby("length", palette = 
#'                 colorRampPalette(c("dodgerblue2", "firebrick2"))),
#'     lwd = 2, spaceHeight = .7,
#'     x = 0.25, y = 0.25, height = 1.5,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bbAnnoGenomeLabel(plot = bedpePlot, x = 0.25, y = 1.78, scale = "Mb")
#'
#' ## Add heatmap legend
#' bbAnnoHeatmapLegend(
#'     plot = bedpePlot, fontcolor = "black",
#'     x = 7.0, y = 0.25,
#'     width = 0.10, height = 1, fontsize = 10
#' )
#'
#' ## Add heatmap legend label
#' bbPlotText(
#'     label = "Kb", rot = 90, x = 6.9, y = 0.75,
#'     just = c("center", "center"), fontsize = 10
#' )
#'
#' ## Hide page guides
#' bbPageGuideHide()
#' @details
#' #' A paired ranges plot can be placed on a BentoBox coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' bbPlotPairs(data, chrom,
#'             chromstart = NULL, chromend = NULL,
#'             x, y, width, height, just = c("left", "top"),
#'             default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated paired
#' ranges plot by ignoring plot placement parameters:
#' \preformatted{
#' bbPlotPairs(data, chrom,
#'             chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
bbPlotPairs <- function(data, chrom, chromstart = NULL, chromend = NULL,
                        assembly = "hg38", fill = "#1f4297",
                        linecolor = NA, bg = NA, boxHeight = unit(2, "mm"),
                        spaceWidth = 0.02, spaceHeight = 0.3,
                        limitLabel = TRUE,
                        baseline = FALSE, baseline.color = "grey",
                        baseline.lwd = 1,
                        x = NULL, y = NULL, width = NULL, height = NULL,
                        just = c("left", "top"), default.units = "inches",
                        draw = TRUE, params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors
    errorcheck_bbPlotPairs <- function(bedpe, bedpePlot, fill) {

        ## Genomic region
        bb_regionErrors(chromstart = bedpePlot$chromstart,
                    chromend = bedpePlot$chromend)
        
        ## Fill colorby checks
        bb_checkColorby(fill = fill,
                        colorby = TRUE,
                        data = bedpe)
            
        }
        
    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_bedpeInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_bedpeInternal"
    )

    ## Parse gp
    bb_bedpeInternal$gp <- setGP(
        gpList = gpar(),
        params = bb_bedpeInternal, ...
    )
    
    ## Justification
    bb_bedpeInternal$just <- bb_justConversion(just = bb_bedpeInternal$just)

    # =========================================================================
    # CHECK ARGUMENT ERRORS
    # =========================================================================
    if (is.null(bb_bedpeInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(bb_bedpeInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    bb_bedpe <- structure(list(
        bedpe = NULL, chrom = bb_bedpeInternal$chrom,
        chromstart = bb_bedpeInternal$chromstart,
        chromend = bb_bedpeInternal$chromend,
        assembly = bb_bedpeInternal$assembly,
        color_palette = NULL,
        zrange = NULL,
        x = bb_bedpeInternal$x, y = bb_bedpeInternal$y,
        width = bb_bedpeInternal$width,
        height = bb_bedpeInternal$height,
        just = bb_bedpeInternal$just, grobs = NULL
    ),
    class = "bb_pairs"
    )
    attr(x = bb_bedpe, which = "plotted") <- bb_bedpeInternal$draw

    # =========================================================================
    # CHECK PLACEMENT
    # =========================================================================

    check_placement(object = bb_bedpe)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    bb_bedpe$assembly <- parse_bbAssembly(assembly = bb_bedpe$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    bb_bedpe <- defaultUnits(
        object = bb_bedpe,
        default.units = bb_bedpeInternal$default.units
    )
    
    bb_bedpeInternal$boxHeight <- misc_defaultUnits(
        value = bb_bedpeInternal$boxHeight,
        name = "boxHeight",
        default.units = bb_bedpeInternal$default.units
    )
    
    # =========================================================================
    # READ IN FILE OR DATAFRAME
    # =========================================================================

    bedpe <- read_pairedData(data = bb_bedpeInternal$data,
                            assembly = bb_bedpe$assembly)
    
    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    errorcheck_bbPlotPairs(
        bedpe = bedpe, bedpePlot = bb_bedpe,
        fill = bb_bedpeInternal$fill
    )
    
    ## chrom format and data chrom format
    chromDataAgreement(data = bedpe, chrom = bb_bedpe$chrom,
                        type = "pairs")

    # =========================================================================
    # ORGANIZE DATA
    # =========================================================================

    ## Get appropriate starts/stops
    start1 <- apply(bedpe[, c("start1", "end1")], 1, min)
    stop1 <- apply(bedpe[, c("start1", "end1")], 1, max)
    start2 <- apply(bedpe[, c("start2", "end2")], 1, min)
    stop2 <- apply(bedpe[, c("start2", "end2")], 1, max)
    bedpe$start1 <- start1
    bedpe$end1 <- stop1
    bedpe$start2 <- start2
    bedpe$end2 <- stop2

    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================

    scaleChecks <- genomicScale(object = bb_bedpe,
                                objectInternal = bb_bedpeInternal,
                                plotType = "paired data plot")
    bb_bedpe <- scaleChecks[[1]]
    bb_bedpeInternal <- scaleChecks[[2]]
    
    # =========================================================================
    # COLORS
    # =========================================================================
    
    pairColors <- bb_parseColors(data = bedpe,
                                fill = bb_bedpeInternal$fill,
                                object = bb_bedpe,
                                subset = "pairs")
    if (length(pairColors[[1]]) > 0){
        bedpe$color <- pairColors[[1]]
    } else {
        bedpe$color <- rep("#1f4297", nrow(bedpe))
    }
    
    bb_bedpe <- pairColors[[2]]
    bedpe$linecolor <- bb_lineColors(linecolor = bb_bedpeInternal$linecolor,
                                    fillcolors = bedpe$color,
                                    data = bedpe,
                                    object = bb_bedpe,
                                    subset = "pairs")

    # =========================================================================
    # SUBSET DATA FOR CHROMOSOME AND ANY OVERLAPPING REGIONS
    # =========================================================================

    if (!is.null(bb_bedpe$chromstart) & !is.null(bb_bedpe$chromend)) {
        bedpe <- bedpe[which(bedpe[, "chrom1"] == bb_bedpe$chrom &
                            bedpe[, "chrom2"] == bb_bedpe$chrom &
                            ((bedpe[, "end1"] >= bb_bedpe$chromstart &
                            bedpe[, "end1"] <= bb_bedpe$chromend) |
                            (bedpe[, "start2"] <= bb_bedpe$chromstart &
                            bedpe[, "start2"] >= bb_bedpe$chromend))), ]
    } else {
        bedpe <- data.frame(matrix(nrow = 0, ncol = 6))
    }

    # =========================================================================
    # GET BOX WIDTHS AND TOTAL DISTANCES
    # =========================================================================

    bedpe$width1 <- bedpe[, "end1"] - bedpe[, "start1"]
    bedpe$width2 <- bedpe[, "end2"] - bedpe[, "start2"]
    bedpe$pos1 <- rowMeans(bedpe[, c("start1", "end1")])
    bedpe$pos2 <- rowMeans(bedpe[, c("start2", "end2")])
    bedpe$distance <- abs(bedpe$pos2 - bedpe$pos1)

    # =========================================================================
    # SORT BY DISTANCE FOR PRETTIER PLOTTING
    # =========================================================================

    bedpe <- bedpe[order(bedpe$distance, decreasing = TRUE), ]

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_pairs",
        length(grep(
            pattern = "bb_pairs",
            x = currentViewports
        )) + 1
    )

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(bb_bedpe$x) | is.null(bb_bedpe$y)) {
        vp <- viewport(
            height = unit(0.5, "snpc"), width = unit(1, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            clip = "on",
            xscale = bb_bedpeInternal$xscale,
            yscale = c(0, 1),
            just = "center",
            name = vp_name
        )

        if (bb_bedpeInternal$draw == TRUE) {
            vp$name <- "bb_pairs1"
            grid.newpage()
        }
    } else {
        add_bbViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = bb_bedpe)

        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = bb_bedpeInternal$xscale,
            yscale = c(0, convertHeight(page_coords$height,
                unitTo = get("page_units",
                    envir = bbEnv
                ),
                valueOnly = TRUE
            )),
            just = bb_bedpeInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================

    backgroundGrob <- rectGrob(gp = gpar(
        fill = bb_bedpeInternal$bg,
        col = NA
    ), name = "background")
    assign("bedpe_grobs", gTree(
        vp = vp,
        children = gList(backgroundGrob)
    ),
    envir = bbEnv
    )

    # =========================================================================
    # DETERMINE ROWS FOR EACH ELEMENT
    # =========================================================================
    if (nrow(bedpe) > 0) {
        
        if (is.null(bb_bedpe$x) & is.null(bb_bedpe$y)) {
            pushViewport(vp)
            boxHeight <- convertHeight(bb_bedpeInternal$boxHeight,
                                    unitTo = "npc", valueOnly = TRUE
            )
            spaceHeight <- boxHeight * (bb_bedpeInternal$spaceHeight)
            upViewport()
        } else {
            boxHeight <- convertHeight(bb_bedpeInternal$boxHeight,
                                    unitTo = get("page_units", 
                                                    envir = bbEnv),
                                    valueOnly = TRUE
            )
            spaceHeight <- boxHeight * (bb_bedpeInternal$spaceHeight)
        }

        ## Determine how many pair elements are going to fit
        maxRows <- floor((as.numeric(vp$height) + spaceHeight) /
                            (boxHeight + spaceHeight))
        wiggle <- abs(bb_bedpe$chromend - bb_bedpe$chromstart) *
            bb_bedpeInternal$spaceWidth
        
        ## Assign rows
        rowData <- assignRows(data = bedpe[,c("start1","end2","start2")],
                            maxRows = maxRows,
                            wiggle = wiggle,
                            rowCol = 3,
                            limitLabel = bb_bedpeInternal$limitLabel,
                            gTree = "bedpe_grobs",
                            extraData = bedpe[,c("color", "linecolor", 
                                                "width1",
                                                "width2", "pos1",
                                                "pos2", "distance")],
                            colNames = c("color", "linecolor", 
                                        "width1", "width2",
                                        "pos1", "pos2", "distance"))
        
        ## Calculate y-coordinates
        rowData$y <- rowData$row * (boxHeight + spaceHeight)
        
        # =====================================================================
        # MAKE GROBS
        # =====================================================================

        if (bb_bedpeInternal$baseline == TRUE) {
            baselineGrob <- segmentsGrob(
                x0 = unit(0, "npc"), y0 = 0,
                x1 = unit(1, "npc"), y1 = 0,
                default.units = "native",
                gp = gpar(
                    col = bb_bedpeInternal$baseline.color,
                    lwd = bb_bedpeInternal$baseline.lwd
                )
            )
            assign("bedpe_grobs",
                addGrob(
                    gTree = get("bedpe_grobs", envir = bbEnv),
                    child = baselineGrob
                ),
                envir = bbEnv
            )
        }

        bb_bedpeInternal$gp$fill <- rowData$color
        bb_bedpeInternal$gp$col <- rowData$linecolor

        bedpeRect1 <- rectGrob(
            x = rowData[,"start1"],
            y = rowData$y,
            width = rowData$width1,
            height = boxHeight,
            just = c("left", "bottom"),
            default.units = "native",
            gp = bb_bedpeInternal$gp
        )

        bedpeRect2 <- rectGrob(
            x = rowData[,"start2"],
            y = rowData$y,
            width = rowData$width2,
            height = boxHeight,
            just = c("left", "bottom"),
            default.units = "native",
            gp = bb_bedpeInternal$gp
        )

        bb_bedpeInternal$gp$col <- rowData$color
        bb_bedpeInternal$gp$lineend <- "butt"

        bedpeLine <- segmentsGrob(
            x0 = rowData$pos1,
            y0 = rowData$y + 0.5 * boxHeight,
            x1 = rowData$pos2,
            y1 = rowData$y + 0.5 * boxHeight,
            default.units = "native",
            gp = bb_bedpeInternal$gp
        )

        assign("bedpe_grobs",
            addGrob(
                gTree = get("bedpe_grobs", envir = bbEnv),
                child = bedpeLine
            ),
            envir = bbEnv
        )
        assign("bedpe_grobs",
            addGrob(
                gTree = get("bedpe_grobs", envir = bbEnv),
                child = bedpeRect1
            ),
            envir = bbEnv
        )
        assign("bedpe_grobs",
            addGrob(
                gTree = get("bedpe_grobs", envir = bbEnv),
                child = bedpeRect2
            ),
            envir = bbEnv
        )
    } else {
        if (bb_bedpeInternal$txdbChecks == TRUE) {
            warning("Data contains no values.", call. = FALSE)
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (bb_bedpeInternal$draw == TRUE) {
        grid.draw(get("bedpe_grobs", envir = bbEnv))
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    bb_bedpe$grobs <- get("bedpe_grobs", envir = bbEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================
    message("bb_pairs[", vp$name, "]")
    invisible(bb_bedpe)
}
