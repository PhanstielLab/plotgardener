#' Plot paired-end genomic range elements
#' 
#' @usage plotPairs(
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
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param fill A single character value, a vector, or 
#' a \link[plotgardener]{colorby} object specifying fill colors of
#' paired range elements. Default value is \code{fill = "#1f4297"}.
#' @param linecolor A single character value, a vector, or a
#' \link[plotgardener]{colorby} object specifying the color of the lines
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
#' plotted plot according to the units of the plotgardener page.
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
#' @param params An optional \link[plotgardener]{params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{pairs} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load paired ranges data in BEDPE format
#' library(plotgardenerData)
#' data("IMR90_DNAloops_pairs")
#'
#' ## Set the coordinates
#' params <- params(
#'     chrom = "chr21",
#'     chromstart = 27900000, chromend = 30700000,
#'     assembly = "hg19",
#'     width = 7
#' )
#'
#' ## Create a page
#' pageCreate(width = 7.5, height = 2.1, default.units = "inches")
#'
#' ## Add a length column
#' IMR90_DNAloops_pairs$length <- 
#'         (IMR90_DNAloops_pairs$start2 - IMR90_DNAloops_pairs$start1) / 1000
#'
#' ## Plot the data
#' bedpePlot <- plotPairs(
#'     data = IMR90_DNAloops_pairs, params = params,
#'     fill = colorby("length", palette = 
#'                 colorRampPalette(c("dodgerblue2", "firebrick2"))),
#'     lwd = 2, spaceHeight = .7,
#'     x = 0.25, y = 0.25, height = 1.5,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' annoGenomeLabel(plot = bedpePlot, x = 0.25, y = 1.78, scale = "Mb")
#'
#' ## Add heatmap legend
#' annoHeatmapLegend(
#'     plot = bedpePlot, fontcolor = "black",
#'     x = 7.0, y = 0.25,
#'     width = 0.10, height = 1, fontsize = 10
#' )
#'
#' ## Add heatmap legend label
#' plotText(
#'     label = "Kb", rot = 90, x = 6.9, y = 0.75,
#'     just = c("center", "center"), fontsize = 10
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' #' A paired ranges plot can be placed on a plotgardener coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' plotPairs(data, chrom,
#'             chromstart = NULL, chromend = NULL,
#'             x, y, width, height, just = c("left", "top"),
#'             default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated paired
#' ranges plot by ignoring plot placement parameters:
#' \preformatted{
#' plotPairs(data, chrom,
#'             chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
plotPairs <- function(data, chrom, chromstart = NULL, chromend = NULL,
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
    errorcheck_plotPairs <- function(bedpeData, bedpePlot, fill) {

        ## Genomic region
        regionErrors(chromstart = bedpePlot$chromstart,
                    chromend = bedpePlot$chromend)
        
        ## Fill colorby checks
        checkColorby(fill = fill,
                        colorby = TRUE,
                        data = bedpeData)
            
        }
        
    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bedpeInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bedpeInternal"
    )

    ## Parse gp
    bedpeInternal$gp <- setGP(
        gpList = gpar(),
        params = bedpeInternal, ...
    )
    
    ## Justification
    bedpeInternal$just <- justConversion(just = bedpeInternal$just)

    # =========================================================================
    # CHECK ARGUMENT ERRORS
    # =========================================================================
    if (is.null(bedpeInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(bedpeInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    bedpe <- structure(list(
        bedpe = NULL, chrom = bedpeInternal$chrom,
        chromstart = bedpeInternal$chromstart,
        chromend = bedpeInternal$chromend,
        assembly = bedpeInternal$assembly,
        color_palette = NULL,
        zrange = NULL,
        x = bedpeInternal$x, y = bedpeInternal$y,
        width = bedpeInternal$width,
        height = bedpeInternal$height,
        just = bedpeInternal$just, grobs = NULL
    ),
    class = "pairs"
    )
    attr(x = bedpe, which = "plotted") <- bedpeInternal$draw

    # =========================================================================
    # CHECK PLACEMENT
    # =========================================================================

    check_placement(object = bedpe)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    bedpe$assembly <- parseAssembly(assembly = bedpe$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    bedpe <- defaultUnits(
        object = bedpe,
        default.units = bedpeInternal$default.units
    )
    
    bedpeInternal$boxHeight <- misc_defaultUnits(
        value = bedpeInternal$boxHeight,
        name = "boxHeight",
        default.units = bedpeInternal$default.units
    )
    
    # =========================================================================
    # READ IN FILE OR DATAFRAME
    # =========================================================================

    bedpeData <- read_pairedData(data = bedpeInternal$data,
                            assembly = bedpe$assembly)
    
    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    errorcheck_plotPairs(
        bedpeData = bedpeData, bedpePlot = bedpe,
        fill = bedpeInternal$fill
    )
    
    ## chrom format and data chrom format
    chromDataAgreement(data = bedpeData, chrom = bedpe$chrom,
                        type = "pairs")

    # =========================================================================
    # ORGANIZE DATA
    # =========================================================================

    ## Get appropriate starts/stops
    start1 <- apply(bedpeData[, c("start1", "end1")], 1, min)
    stop1 <- apply(bedpeData[, c("start1", "end1")], 1, max)
    start2 <- apply(bedpeData[, c("start2", "end2")], 1, min)
    stop2 <- apply(bedpeData[, c("start2", "end2")], 1, max)
    bedpeData$start1 <- start1
    bedpeData$end1 <- stop1
    bedpeData$start2 <- start2
    bedpeData$end2 <- stop2

    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================

    scaleChecks <- genomicScale(object = bedpe,
                                objectInternal = bedpeInternal,
                                plotType = "paired data plot")
    bedpe <- scaleChecks[[1]]
    bedpeInternal <- scaleChecks[[2]]
    
    # =========================================================================
    # COLORS
    # =========================================================================
    
    pairColors <- parseColors(data = bedpeData,
                                fill = bedpeInternal$fill,
                                object = bedpe,
                                subset = "pairs")
    if (length(pairColors[[1]]) > 0){
        bedpe$color <- pairColors[[1]]
    } else {
        bedpe$color <- rep("#1f4297", nrow(bedpe))
    }
    
    bedpe <- pairColors[[2]]
    bedpe$linecolor <- bb_lineColors(linecolor = bedpeInternal$linecolor,
                                    fillcolors = bedpe$color,
                                    data = bedpeData,
                                    object = bedpe,
                                    subset = "pairs")

    # =========================================================================
    # SUBSET DATA FOR CHROMOSOME AND ANY OVERLAPPING REGIONS
    # =========================================================================

    if (!is.null(bedpe$chromstart) & !is.null(bedpe$chromend)) {
        bedpeData <- bedpeData[which(bedpeData[, "chrom1"] == bedpeData$chrom &
                            bedpeData[, "chrom2"] == bedpeData$chrom &
                            ((bedpeData[, "end1"] >= bedpeData$chromstart &
                            bedpeData[, "end1"] <= bedpeData$chromend) |
                            (bedpeData[, "start2"] <= bedpeData$chromstart &
                            bedpeData[, "start2"] >= bedpeData$chromend))), ]
    } else {
        bedpeData <- data.frame(matrix(nrow = 0, ncol = 6))
    }

    # =========================================================================
    # GET BOX WIDTHS AND TOTAL DISTANCES
    # =========================================================================

    bedpeData$width1 <- bedpeData[, "end1"] - bedpeData[, "start1"]
    bedpeData$width2 <- bedpeData[, "end2"] - bedpeData[, "start2"]
    bedpeData$pos1 <- rowMeans(bedpeData[, c("start1", "end1")])
    bedpeData$pos2 <- rowMeans(bedpeData[, c("start2", "end2")])
    bedpeData$distance <- abs(bedpeData$pos2 - bedpeData$pos1)

    # =========================================================================
    # SORT BY DISTANCE FOR PRETTIER PLOTTING
    # =========================================================================

    bedpeData <- bedpeData[order(bedpeData$distance, decreasing = TRUE), ]

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
    if (is.null(bedpe$x) | is.null(bedpe$y)) {
        vp <- viewport(
            height = unit(0.5, "snpc"), width = unit(1, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            clip = "on",
            xscale = bedpeInternal$xscale,
            yscale = c(0, 1),
            just = "center",
            name = vp_name
        )

        if (bedpeInternal$draw == TRUE) {
            vp$name <- "bb_pairs1"
            grid.newpage()
        }
    } else {
        addViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = bedpe)

        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = bedpeInternal$xscale,
            yscale = c(0, convertHeight(page_coords$height,
                unitTo = get("page_units",
                    envir = pgEnv
                ),
                valueOnly = TRUE
            )),
            just = bedpeInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================

    backgroundGrob <- rectGrob(gp = gpar(
        fill = bedpeInternal$bg,
        col = NA
    ), name = "background")
    assign("bedpe_grobs", gTree(
        vp = vp,
        children = gList(backgroundGrob)
    ),
    envir = pgEnv
    )

    # =========================================================================
    # DETERMINE ROWS FOR EACH ELEMENT
    # =========================================================================
    if (nrow(bedpeData) > 0) {
        
        if (is.null(bedpe$x) & is.null(bedpe$y)) {
            pushViewport(vp)
            boxHeight <- convertHeight(bedpeInternal$boxHeight,
                                    unitTo = "npc", valueOnly = TRUE
            )
            spaceHeight <- boxHeight * (bedpeInternal$spaceHeight)
            upViewport()
        } else {
            boxHeight <- convertHeight(bedpeInternal$boxHeight,
                                    unitTo = get("page_units", 
                                                    envir = pgEnv),
                                    valueOnly = TRUE
            )
            spaceHeight <- boxHeight * (bedpeInternal$spaceHeight)
        }

        ## Determine how many pair elements are going to fit
        maxRows <- floor((as.numeric(vp$height) + spaceHeight) /
                            (boxHeight + spaceHeight))
        wiggle <- abs(bedpe$chromend - bedpe$chromstart) *
            bedpeInternal$spaceWidth
        
        ## Assign rows
        rowData <- assignRows(data = bedpeData[,c("start1","end2","start2")],
                            maxRows = maxRows,
                            wiggle = wiggle,
                            rowCol = 3,
                            limitLabel = bedpeInternal$limitLabel,
                            gTree = "bedpe_grobs",
                            extraData = bedpeData[,c("color", "linecolor", 
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

        if (bedpeInternal$baseline == TRUE) {
            baselineGrob <- segmentsGrob(
                x0 = unit(0, "npc"), y0 = 0,
                x1 = unit(1, "npc"), y1 = 0,
                default.units = "native",
                gp = gpar(
                    col = bedpeInternal$baseline.color,
                    lwd = bedpeInternal$baseline.lwd
                )
            )
            assign("bedpe_grobs",
                addGrob(
                    gTree = get("bedpe_grobs", envir = pgEnv),
                    child = baselineGrob
                ),
                envir = pgEnv
            )
        }

        bedpeInternal$gp$fill <- rowData$color
        bedpeInternal$gp$col <- rowData$linecolor

        bedpeRect1 <- rectGrob(
            x = rowData[,"start1"],
            y = rowData$y,
            width = rowData$width1,
            height = boxHeight,
            just = c("left", "bottom"),
            default.units = "native",
            gp = bedpeInternal$gp
        )

        bedpeRect2 <- rectGrob(
            x = rowData[,"start2"],
            y = rowData$y,
            width = rowData$width2,
            height = boxHeight,
            just = c("left", "bottom"),
            default.units = "native",
            gp = bedpeInternal$gp
        )

        bedpeInternal$gp$col <- rowData$color
        bedpeInternal$gp$lineend <- "butt"

        bedpeLine <- segmentsGrob(
            x0 = rowData$pos1,
            y0 = rowData$y + 0.5 * boxHeight,
            x1 = rowData$pos2,
            y1 = rowData$y + 0.5 * boxHeight,
            default.units = "native",
            gp = bedpeInternal$gp
        )

        assign("bedpe_grobs",
            addGrob(
                gTree = get("bedpe_grobs", envir = pgEnv),
                child = bedpeLine
            ),
            envir = pgEnv
        )
        assign("bedpe_grobs",
            addGrob(
                gTree = get("bedpe_grobs", envir = pgEnv),
                child = bedpeRect1
            ),
            envir = pgEnv
        )
        assign("bedpe_grobs",
            addGrob(
                gTree = get("bedpe_grobs", envir = pgEnv),
                child = bedpeRect2
            ),
            envir = pgEnv
        )
    } else {
        if (bedpeInternal$txdbChecks == TRUE) {
            warning("Data contains no values.", call. = FALSE)
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (bedpeInternal$draw == TRUE) {
        grid.draw(get("bedpe_grobs", envir = pgEnv))
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    bedpe$grobs <- get("bedpe_grobs", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================
    message("pairs[", vp$name, "]")
    invisible(bedpe)
}
