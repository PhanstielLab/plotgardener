#' Plot genomic range elements in a pileup or collapsed format
#' 
#' @usage bbPlotRanges(
#'     data,
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
#'     fill = "#7ecdbb",
#'     linecolor = NA,
#'     collapse = FALSE,
#'     boxHeight = unit(2, "mm"),
#'     spaceWidth = 0.02,
#'     spaceHeight = 0.3,
#'     limitLabel = TRUE,
#'     strandSplit = FALSE,
#'     bg = NA,
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
#' @param data Data to be plotted; as a character value specifying
#' a BED file path, a data frame in BED format, a character value
#' specifying a .bam file path where a bam index file (.bam.bai)
#' is in the same directory, or a \link[GenomicRanges]{GRanges} object.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[BentoBox]{bbAssembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param fill A single character value, a vector, or a 
#' \link[BentoBox]{colorby} object specifying fill colors of range elements.
#' Default value is \code{fill = "#7ecdbb"}.
#' @param linecolor A single character value, a vector, or a
#' \link[BentoBox]{colorby} object specifying the color of the lines
#' outlining range elements. Default value is \code{linecolor = NA}.
#' Special options include:
#' \itemize{
#' \item{\code{NA}: }{No line color.}
#' \item{\code{"fill"}: }{Same color as \code{fill}.}
#' } .
#' @param collapse A logical value indicating whether to collapse
#' range elements into a single row, or into
#' two rows if \code{strandSplit = TRUE}.
#' If \code{collapse = TRUE}, \code{boxHeight} will be ignored and elements
#' will be the height of the entire plot if \code{strandSplit = FALSE} or
#' be the height of half of the entire plot if \code{strandSplit = TRUE}.
#' Default value is \code{collapse = FALSE}.
#' @param boxHeight A numeric or unit object specifying height of range element
#' boxes. Default value is \code{boxHeight = unit(2, "mm")}.
#' @param spaceWidth A numeric value specifying the width of minimum spacing
#' between range element boxes, as a fraction of the plot's genomic range.
#' Default value is \code{spaceWidth = 0.02}.
#' @param spaceHeight A numeric value specifying the height of spacing between
#' range element boxes on different rows, as a fraction of boxHeight.
#' Default value is \code{spaceHeight = 0.3}.
#' @param limitLabel A logical value indicating whether to draw a "+"
#' when not all elements can be plotted in the plotting space. Default 
#' value is \code{limitLabel = TRUE}.
#' @param strandSplit A logical value indicating whether plus and
#' minus-stranded elements should be separated. Elements can only be
#' split by strand if a \code{strand} column is found in \code{data}.
#' Default value is \code{strandSplit = FALSE}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a
#' baseline along the x-axis. Default value is \code{baseline = FALSE}.
#' @param baseline.color Baseline color.
#' Default value is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width.
#' Default value is \code{baseline.lwd = 1}.
#' @param x A numeric or unit object specifying ranges plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying ranges plot y-location.
#' The character value will
#' place the ranges plot y relative to the bottom of the most recently
#' plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying ranges plot width.
#' @param height A numeric or unit object specifying ranges plot height.
#' @param just Justification of ranges plot relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use
#' if \code{x}, \code{y}, \code{width}, or \code{height} are only given
#' as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should
#' be produced. Default value \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bbParams} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_ranges} object containing relevant
#' genomic region, coloring data, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load ranges data in BED format
#' library(BentoBoxData)
#' data("IMR90_ChIP_CTCF_reads")
#'
#' ## Create page
#' bbPageCreate(width = 7.5, height = 5, default.units = "inches")
#'
#' ## Plot and place a pileup ranges plot
#' pileupPlot <- bbPlotRanges(
#'     data = IMR90_ChIP_CTCF_reads, chrom = "chr21",
#'     chromstart = 29073000, chromend = 29074000,
#'     assembly = "hg19",
#'     fill = colorby("strand", palette = 
#'                 colorRampPalette(c("#7ecdbb", "#37a7db"))),
#'     strandSplit = TRUE, 
#'     x = 0.5, y = 0.25, width = 6.5, height = 4.25,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bbAnnoGenomeLabel(
#'     plot = pileupPlot, x = 0.5, y = 4.5,
#'     just = c("left", "top")
#' )
#'
#' ## Add text labels
#' bbPlotText(
#'     label = "+ strand", fontcolor = "#37a7db", fontsize = 12,
#'     x = 0.5, y = 1.25, just = "left"
#' )
#' bbPlotText(
#'     label = "- strand", fontcolor = "#7ecdbb", fontsize = 12,
#'     x = 0.5, y = 3.5, just = "left"
#' )
#'
#' ## Hide page guides
#' bbPageGuideHide()
#' @details
#' A ranges plot can be placed on a BentoBox coordinate page by providing
#' plot placement parameters:
#' \preformatted{
#' bbPlotRanges(data, chrom,
#'             chromstart = NULL, chromend = NULL,
#'             x, y, width, height, just = c("left", "top"),
#'             default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated BED plot
#' by ignoring plot placement parameters:
#' \preformatted{
#' bbPlotRanges(data, chrom,
#'             chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
bbPlotRanges <- function(data, chrom, chromstart = NULL, chromend = NULL,
                        assembly = "hg38", fill = "#7ecdbb",
                        linecolor = NA, collapse = FALSE, 
                        boxHeight = unit(2, "mm"), spaceWidth = 0.02,
                        spaceHeight = 0.3, limitLabel = TRUE,
                        strandSplit = FALSE, bg = NA,
                        baseline = FALSE, baseline.color = "grey",
                        baseline.lwd = 1,
                        x = NULL, y = NULL, width = NULL, height = NULL,
                        just = c("left", "top"), default.units = "inches",
                        draw = TRUE, params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that catches errors
    errorcheck_bbPlotpileup <- function(pileupPlot, fill, bed) {

        ## Genomic region
        bb_regionErrors(chromstart = pileupPlot$chromstart,
                    chromend = pileupPlot$chromend)
        
        ## Fill/colorby checks
        bb_checkColorby(fill = fill,
                        colorby = TRUE,
                        data = bed)
    }

    ## Define a function that parses the yscale based on split strands
    strand_scale <- function(strandSplit, height) {
        if (strandSplit == TRUE) {
            yscale <- c(-height / 2, height / 2)
        } else {
            yscale <- c(0, height)
        }

        return(yscale)
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_pileInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_pileInternal"
    )

    ## Parse gp
    bb_pileInternal$gp <- setGP(
        gpList = gpar(),
        params = bb_pileInternal, ...
    )
    
    ## Justification
    bb_pileInternal$just <- bb_justConversion(just = bb_pileInternal$just)

    # =========================================================================
    # CHECK ARGUMENT ERROS
    # =========================================================================

    if (is.null(bb_pileInternal$data)) stop("argument \"data\" is missing, ",
                                            "with no default.", call. = FALSE)
    if (is.null(bb_pileInternal$chrom)) stop("argument \"chrom\" is missing, ",
                                            "with no default.", call. = FALSE)
    
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    pileupPlot <- structure(list(
        chrom = bb_pileInternal$chrom,
        chromstart = bb_pileInternal$chromstart,
        chromend = bb_pileInternal$chromend,
        assembly = bb_pileInternal$assembly,
        color_palette = NULL,
        zrange = NULL,
        x = bb_pileInternal$x, y = bb_pileInternal$y,
        width = bb_pileInternal$width,
        height = bb_pileInternal$height,
        just = bb_pileInternal$just, grobs = NULL
    ),
    class = "bb_ranges"
    )
    attr(x = pileupPlot, which = "plotted") <- bb_pileInternal$draw

    # =========================================================================
    # CHECK PLACEMENT ERRORS
    # =========================================================================

    check_placement(object = pileupPlot)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    pileupPlot$assembly <- parse_bbAssembly(assembly = pileupPlot$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    pileupPlot <- defaultUnits(
        object = pileupPlot,
        default.units = bb_pileInternal$default.units
    )
    
    bb_pileInternal$boxHeight <- misc_defaultUnits(
        value = bb_pileInternal$boxHeight,
        name = "boxHeight",
        default.units = bb_pileInternal$default.units
    )

    # =========================================================================
    # READ IN FILE OR DATAFRAME
    # =========================================================================

    bed <- read_rangeData(data = bb_pileInternal$data,
                        assembly = pileupPlot$assembly,
                        chrom = pileupPlot$chrom,
                        start = pileupPlot$chromstart,
                        end = pileupPlot$chromend)

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    errorcheck_bbPlotpileup(
        pileupPlot = pileupPlot,
        fill = bb_pileInternal$fill,
        bed = bed
    )

    ## chrom format and data chrom format
    chromDataAgreement(data = bed, chrom = pileupPlot$chrom,
                    type = "ranges")
    
    # =========================================================================
    # GENOMIC SCALE
    # =========================================================================
    
    scaleChecks <- genomicScale(object = pileupPlot,
                                objectInternal = bb_pileInternal,
                                plotType = "ranges plot")
    pileupPlot <- scaleChecks[[1]]
    bb_pileInternal <- scaleChecks[[2]]
    
    # =========================================================================
    # COLORS
    # =========================================================================
    
    rangeColors <- bb_parseColors(data = bed,
                                fill = bb_pileInternal$fill,
                                object = pileupPlot,
                                subset = "ranges")
    if (length(rangeColors[[1]]) > 0){
        bed$color <- rangeColors[[1]]
    } else {
        bed$color <- rep("#7ecdbb", nrow(bed))
    }
    
    pileupPlot <- rangeColors[[2]]
    bed$linecolor <- bb_lineColors(linecolor = bb_pileInternal$linecolor,
                                fillcolors = bed$color,
                                data = bed,
                                object = pileupPlot,
                                subset = "ranges")
    
    # =========================================================================
    # SUBSET DATA FOR CHROMOSOME AND ANY OVERLAPPING REGIONS
    # =========================================================================

    if (!is.null(pileupPlot$chromstart) & !is.null(pileupPlot$chromend)) {
        bed <- bed[which(bed[, "chrom"] == pileupPlot$chrom 
                    & bed[, "start"] <=
            pileupPlot$chromend & bed[, "end"] >=
            pileupPlot$chromstart), ]
    } else {
        bed <- data.frame(matrix(nrow = 0, ncol = 3))
    }

    # =========================================================================
    # SEPARATE DATA INTO STRANDS
    # =========================================================================

    if (bb_pileInternal$strandSplit == TRUE) {

        ## Look for column named 'strand'
        if (!any(colnames(bed) == "strand")) {
            stop("No `strand` column found in data. ",
                "Cannot split data based on strand.", call. = FALSE)
        } else {
            posStrand <- bed[which(bed[, "strand"] == "+"), ]
            minStrand <- bed[which(bed[, "strand"] == "-"), ]
        }
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_ranges",
        length(grep(
            pattern = "bb_ranges",
            x = currentViewports
        )) + 1
    )


    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(pileupPlot$x) | is.null(pileupPlot$y)) {
        yscale <- strand_scale(
            strandSplit = bb_pileInternal$strandSplit,
            height = 0.5
        )

        vp <- viewport(
            height = unit(0.5, "snpc"), width = unit(1, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            clip = "on",
            xscale = bb_pileInternal$xscale,
            yscale = yscale,
            just = "center",
            name = vp_name
        )

        if (bb_pileInternal$draw == TRUE) {
            vp$name <- "bb_ranges1"
            grid.newpage()
        }
    } else {
        add_bbViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = pileupPlot)

        yscale <- strand_scale(
            strandSplit = bb_pileInternal$strandSplit,
            height = convertHeight(page_coords$height,
                unitTo = get("page_units",
                    envir = bbEnv
                ),
                valueOnly = TRUE
            )
        )

        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = bb_pileInternal$xscale,
            yscale = yscale,
            just = bb_pileInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================

    backgroundGrob <- rectGrob(
        gp = gpar(
            fill = bb_pileInternal$bg,
            col = NA
        ),
        name = "background"
    )
    assign("pileup_grobs", gTree(
        vp = vp,
        children = gList(backgroundGrob)
    ),
    envir = bbEnv
    )

    # =========================================================================
    # DETERMINE ROWS FOR EACH ELEMENT
    # =========================================================================

    if (is.null(pileupPlot$x) & is.null(pileupPlot$y)) {
        pushViewport(vp)
        boxHeight <- convertHeight(bb_pileInternal$boxHeight,
            unitTo = "npc", valueOnly = TRUE
        )
        spaceHeight <- boxHeight * (bb_pileInternal$spaceHeight)
        upViewport()
    } else {
        boxHeight <- convertHeight(bb_pileInternal$boxHeight,
            unitTo = get("page_units",
                envir = bbEnv
            ),
            valueOnly = TRUE
        )
        spaceHeight <- boxHeight * (bb_pileInternal$spaceHeight)
    }


    if (bb_pileInternal$collapse == FALSE) {
        
        ## Calculate number of element rows that will fit
        maxRows <- floor((as.numeric(vp$height) + spaceHeight) /
            (boxHeight + spaceHeight))
        wiggle <- abs(pileupPlot$chromend - pileupPlot$chromstart) *
            bb_pileInternal$spaceWidth

        if (bb_pileInternal$strandSplit == FALSE) {
            
            ## Randomize order of data
            bed <- bed[sample(nrow(bed)), ]
            
            ## Assign rows
            rowData <- assignRows(data = bed[,c("start", "end")], 
                                maxRows = maxRows,
                                wiggle = wiggle, rowCol = 2,
                                limitLabel = bb_pileInternal$limitLabel,
                                gTree = "pileup_grobs",
                                extraData = bed[,c("color", "linecolor")],
                                colNames = c("color", "linecolor"))
            
            ## Calculate y-coordinates
            rowData$y <- rowData$row * (boxHeight + spaceHeight)

        } else {
            
            ## Randomize order of data
            posStrand <- posStrand[sample(nrow(posStrand)), ]
            minStrand <- minStrand[sample(nrow(minStrand)), ]
            
            ## Assign rows
            posData <- assignRows(data = posStrand[,c("start", "end")],
                                maxRows = maxRows * 0.5,
                                wiggle = wiggle, rowCol = 2,
                                limitLabel = bb_pileInternal$limitLabel,
                                gTree = "pileup_grobs",
                                extraData = posStrand[,c("color", "linecolor")],
                                colNames = c("color", "linecolor"))
            minData <- assignRows(data = minStrand[,c("start", "end")],
                                maxRows = maxRows * 0.5,
                                wiggle = wiggle, rowCol = 2,
                                limitLabel = bb_pileInternal$limitLabel,
                                side = "bottom",
                                gTree = "pileup_grobs",
                                extraData = minStrand[,c("color", "linecolor")],
                                colNames = c("color", "linecolor"))
            
            ## Calculate y-coordinates
            posData$y <- (0.5 * spaceHeight) + posData$row *
                (boxHeight + spaceHeight)
            minData$y <- ((0.5 * spaceHeight + boxHeight) + minData$row *
                                (boxHeight + spaceHeight)) * -1
            
            ## Combine plus and minus strand data
            rowData <- rbind(posData, minData)
        }

    } else {
        if (bb_pileInternal$strandSplit == FALSE) {
            boxHeight <- as.numeric(vp$height)
            rowData <- bed
            ## y-coordinates all the same
            rowData$y <- rep(0, nrow(rowData))
        } else {
            boxHeight <- (1 - spaceHeight) * as.numeric(vp$height) * 0.5
            if (nrow(posStrand) > 0) {
                posStrand$y <- vp$yscale[2] - boxHeight
            } else {
                posStrand <- data.frame()
            }
            if (nrow(minStrand) > 0) {
                minStrand$y <- vp$yscale[1]
            } else {
                minStrand <- data.frame()
            }

            ## Combine plus and minus strand data
            rowData <- rbind(posStrand, minStrand)
            
        }
    }
    
    
    # =====================================================================
    # MAKE GROBS
    # =====================================================================
    
    if (nrow(rowData) > 0) {
        
        rowData$width <- rowData[,"end"] - rowData[,"start"]

        alpha <- 1
        if ("alpha" %in% names(bb_pileInternal$gp)) {
            alpha <- bb_pileInternal$gp$alpha
        }

        bedRects <- rectGrob(
            x = rowData$start,
            y = rowData$y,
            width = rowData$width,
            height = boxHeight,
            just = c("left", "bottom"),
            default.units = "native",
            gp = gpar(
                fill = rowData$color,
                col = rowData$linecolor,
                alpha = alpha
            )
        )
        assign("pileup_grobs",
            addGrob(
                gTree = get("pileup_grobs", envir = bbEnv),
                child = bedRects
            ),
            envir = bbEnv
        )

        if ((bb_pileInternal$strandSplit == TRUE |
            bb_pileInternal$baseline == TRUE) &
            bb_pileInternal$collapse == FALSE) {
            lineGrob <- segmentsGrob(
                x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                y0 = unit(0, "native"), y1 = unit(0, "native"),
                gp = gpar(
                    col = bb_pileInternal$baseline.color,
                    lwd = bb_pileInternal$baseline.lwd
                )
            )
            assign("pileup_grobs",
                addGrob(
                    gTree = get("pileup_grobs", envir = bbEnv),
                    child = lineGrob
                ),
                envir = bbEnv
            )
        }
    } else {
        if (bb_pileInternal$txdbChecks == TRUE) {
            warning("No range data to plot.", call. = FALSE)
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (bb_pileInternal$draw == TRUE) {
        grid.draw(get("pileup_grobs", envir = bbEnv))
    }
    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    pileupPlot$grobs <- get("pileup_grobs", envir = bbEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_ranges[", vp$name, "]")
    invisible(pileupPlot)
}
