#' Plot gene transcripts in a pileup style for a single chromosome
#' 
#' @usage plotTranscripts(
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
#'     fill = c("#669fd9", "#abcc8e"),
#'     colorbyStrand = TRUE,
#'     strandSplit = FALSE,
#'     boxHeight = unit(2, "mm"),
#'     spaceWidth = 0.02,
#'     spaceHeight = 0.3,
#'     limitLabel = TRUE,
#'     transcriptHighlights = NULL,
#'     fontsize = 8,
#'     labels = "transcript",
#'     stroke = 0.1,
#'     bg = NA,
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
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param fill Character value(s) as a single value or vector
#' specifying fill colors of transcripts.
#' Default value is \code{fill = c("#669fd9", "#abcc8e")}.
#' @param colorbyStrand A logical value indicating whether to
#' color plus and minus strands by the first two colors in
#' a \code{fill} vector, where plus strand transcripts will be
#' colored by the first \code{fill} color and
#' minus strand transcripts will be colored by the second \code{fill}
#' color. Default value is \code{colorbyStrand = TRUE}.
#' @param strandSplit A logical value indicating whether plus and
#' minus-stranded transcripts should be separated, with plus strand
#' transcripts plotted above the x-axis and minus strand transcripts
#' plotted below the x-axis. Default value is \code{strandSplit = FALSE}.
#' @param boxHeight A numeric or unit object specifying height of transcripts.
#' Default value is \code{boxHeight = unit(2, "mm")}.
#' @param spaceWidth A numeric value specifying the width of minimum spacing
#' between transcripts, as a fraction of the plot's genomic range.
#' Default value is \code{spaceWidth = 0.02}.
#' @param spaceHeight A numeric value specifying the height of spacing
#' between transcripts on different rows, as a fraction of \code{boxHeight}.
#' Default value is \code{spaceHeight = 0.3}.
#' @param limitLabel A logical value indicating whether to draw a "+"
#' when not all elements can be plotted in the plotting space. Default 
#' value is \code{limitLabel = TRUE}.
#' @param transcriptHighlights A two-column dataframe with a column named
#' "transcript" or "gene" containing transcript names or their associated gene 
#' names as strings to highlight and a column named "color" containing 
#' corresponding highlight colors.
#' @param fontsize A numeric specifying text fontsize in points.
#' Default value is \code{fontsize = 8}.
#' @param labels A character value describing the format of
#' transcript text labels. Default value is \code{labels = "trancript"}.
#' Options are:
#' \itemize{
#' \item{\code{NULL}: }{No labels.}
#' \item{\code{"transcript"}: }{Transcript name labels.}
#' \item{\code{"gene"}: }{Gene name labels.}
#' \item{\code{"both"}: }{Combined transcript and gene name labels with
#' the format "gene name:transcript name".}
#' }
#' @param stroke A numeric value indicating the stroke width for
#' transcript body outlines. Default value is \code{stroke = 0.1}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param x A numeric or unit object specifying transcript plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying transcript plot y-location.
#' The character value will
#' place the transcript plot y relative to the bottom of the most recently
#' plotted plot according to the units of the plotgardener page.
#' @param width A numeric or unit object specifying transcript plot width.
#' @param height A numeric or unit object specifying transcript plot height.
#' @param just Justification of transcript plot relative to
#' its (x, y) location. If there are two values, the first value specifies
#' horizontal justification and the second value specifies vertical
#' justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if
#' \code{x}, \code{y}, \code{width}, or \code{height} are only given as
#' numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be
#' produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[plotgardener]{pgParams} object containing
#' relevant function parameters.
#'
#' @return Returns a \code{transcripts} object containing relevant
#' genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load hg19 genomic annotation packages
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Create page
#' pageCreate(width = 7.5, height = 3.5, default.units = "inches")
#'
#' ## Plot and place transcripts
#' plotTranscripts(
#'     chrom = "chr8", chromstart = 1000000, chromend = 2000000,
#'     assembly = "hg19", labels = "gene",
#'     x = 0.5, y = 0.5, width = 6.5, height = 2.5,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Plot genome label
#' plotGenomeLabel(
#'     chrom = "chr8", chromstart = 1000000, chromend = 2000000,
#'     assembly = "hg19",
#'     x = 0.5, y = 3.03, length = 6.5, default.units = "inches"
#' )
#'
#' ## Plot a legend
#' plotLegend(
#'     legend = c("+ strand", "- strand"),
#'     fill = c("#669fd9", "#abcc8e"), border = FALSE,
#'     x = 0.5, y = 1, width = 1, height = 0.5,
#'     just = c("left", "top")
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' A transcripts plot can be placed on a plotgardener coordinate page
#' by providing plot placement parameters:
#' \preformatted{
#' plotTranscripts(chrom, chromstart = NULL, chromend = NULL,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated
#' transcripts plot by ignoring plot placement parameters:
#' \preformatted{
#' plotTranscripts(chrom, chromstart = NULL, chromend = NULL)
#' }
#'
#' Genomic annotation information is acquired through
#' \link[GenomicFeatures]{TxDb} and \link[AnnotationDbi]{OrgDb-class} packages,
#' as determined through the \code{assembly} parameter.
#'
#' @seealso \link[plotgardener]{assembly},
#' \link[plotgardener]{genomes}, \link[plotgardener]{defaultPackages}
#'
#' @export
plotTranscripts <- function(chrom, chromstart = NULL, chromend = NULL,
                            assembly = "hg38",
                            fill = c("#669fd9", "#abcc8e"),
                            colorbyStrand = TRUE, strandSplit = FALSE,
                            boxHeight = unit(2, "mm"), spaceWidth = 0.02,
                            spaceHeight = 0.3, limitLabel = TRUE,
                            transcriptHighlights = NULL,
                            fontsize = 8,
                            labels = "transcript", stroke = 0.1, bg = NA,
                            x = NULL, y = NULL, width = NULL,
                            height = NULL, just = c("left", "top"),
                            default.units = "inches", draw = TRUE,
                            params = NULL) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that checks errors for plotTranscripts
    errorcheck_plotTranscripts <- function(transcriptPlot, labels, fill) {

        ## Genomic region
        regionErrors(chromstart = transcriptPlot$chromstart,
                    chromend = transcriptPlot$chromend)


        if (!labels %in% c(NULL, "transcript", "gene", "both")) {
            stop("Invalid \'labels\' input. Options are \'NULL\', ",
                "\'transcript\', \'gene\', or \'both\'.", call. = FALSE)
        }
        
        checkColorby(fill = fill,
                        colorby = FALSE)
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

    transcriptsInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "transcriptsInternal"
    )
    
    ## Justification
    transcriptsInternal$just <- 
        justConversion(just = transcriptsInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    transcripts <- structure(list(
        chrom = transcriptsInternal$chrom,
        chromstart = transcriptsInternal$chromstart,
        chromend = transcriptsInternal$chromend,
        assembly = transcriptsInternal$assembly,
        x = transcriptsInternal$x,
        y = transcriptsInternal$y,
        width = transcriptsInternal$width,
        height = transcriptsInternal$height,
        just = transcriptsInternal$just,
        grobs = NULL
    ), class = "transcripts")
    attr(x = transcripts, which = "plotted") <- transcriptsInternal$draw

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(transcripts$chrom)) {
        stop("argument \"chrom\" is missing, with no default.",
            call. = FALSE
        )
    }
    check_placement(object = transcripts)
    errorcheck_plotTranscripts(
        transcriptPlot = transcripts,
        labels = transcriptsInternal$labels,
        fill = transcriptsInternal$fill
    )

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    transcripts$assembly <-
        parseAssembly(assembly = transcripts$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    transcripts <- defaultUnits(
        object = transcripts,
        default.units = transcriptsInternal$default.units
    )
    
    transcriptsInternal$boxHeight <- misc_defaultUnits(
        value = transcriptsInternal$boxHeight,
        name = "boxHeight",
        default.units = transcriptsInternal$default.units
    )
    
    # =========================================================================
    # GET APPROPRIATE BUILD DATA
    # =========================================================================

    buildData <- geneData(object = transcripts,
                        objectInternal = transcriptsInternal)
    transcripts <- buildData[[1]]
    transcriptsInternal <- buildData[[2]]
    data <- transcriptsInternal$data
    ## Get transcript lengths
    data$length <- data$TXEND - data$TXSTART
    
    # =========================================================================
    # COLORS AND HIGHLIGHTS
    # =========================================================================

    if (transcriptsInternal$colorbyStrand == TRUE) {
        if (length(transcriptsInternal$fill) == 1) {
            posCol <- transcriptsInternal$fill
            negCol <- transcriptsInternal$fill
        } else {
            posCol <- transcriptsInternal$fill[1]
            negCol <- transcriptsInternal$fill[2]
        }

        pos <- data[which(data$TXSTRAND == "+"), ]
        pos$color <- rep(posCol, nrow(pos))
        neg <- data[which(data$TXSTRAND == "-"), ]
        neg$color <- rep(negCol, nrow(neg))
        data <- rbind(pos, neg)
    } else {
        data$color <- rep(transcriptsInternal$fill[1], nrow(data))
    }
    
    ## Find transcripts to be highlighted
    if (!is.null(transcriptsInternal$transcriptHighlights)){
        
        # check column name of transcriptHighlights
        if (colnames(transcriptsInternal$transcriptHighlights)[1] == "transcript"){
            transcriptHighlights <- data[which(data[,"TXNAME"] %in% 
                transcriptsInternal$transcriptHighlights[,"transcript"]),]
            transcriptHighlights$color <- NULL
            transcriptBackground <- data[which(!data[,"TXNAME"] %in% 
                transcriptsInternal$transcriptHighlights[,"transcript"]),]
            # Change highlight transcripts to highlight color
            transcriptHighlights <- 
                merge(transcriptHighlights, 
                    transcriptsInternal$transcriptHighlights,
                    by.x = "TXNAME", by.y = "transcript")
            
        } else if 
        (colnames(transcriptsInternal$transcriptHighlights)[1] == "gene"){
            transcriptHighlights <- 
                data[which(data[[transcripts$assembly$display.column]] %in% 
                    transcriptsInternal$transcriptHighlights[,"gene"]),]
            transcriptHighlights$color <- NULL
            transcriptBackground <- 
                data[which(!data[[transcripts$assembly$display.column]] %in% 
                    transcriptsInternal$transcriptHighlights[,"gene"]),]
            # Change highlight transcripts to highlight color
            transcriptHighlights <- 
                merge(transcriptHighlights, 
                    transcriptsInternal$transcriptHighlights,
                    by.x = transcripts$assembly$display.column, by.y = "gene")
            
        }
        
        ## Put highlight and other transcripts back together
        data <- rbind(transcriptHighlights, transcriptBackground)
    }
    

    # =========================================================================
    # SEPARATE DATA INTO STRANDS
    # =========================================================================

    if (transcriptsInternal$strandSplit == TRUE) {
        posStrand <- data[which(data$TXSTRAND == "+"), ]
        minStrand <- data[which(data$TXSTRAND == "-"), ]
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    if (is.null(transcripts$x) | is.null(transcripts$y)) {
        height <- 0.5
        yscale <- strand_scale(
            strandSplit = transcriptsInternal$strandSplit,
            height = height
        )

        vp <- viewport(
            height = unit(0.5, "snpc"), width = unit(1, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            clip = "on",
            xscale = transcriptsInternal$xscale,
            yscale = yscale,
            just = "center",
            name = "transcripts1"
        )

        if (transcriptsInternal$draw == TRUE) {
            grid.newpage()
        }
    } else {
        ## Name viewport
        currentViewports <- current_viewports()
        vp_name <- paste0(
            "transcripts",
            length(grep(
                pattern = "transcripts",
                x = currentViewports
            )) + 1
        )
        
        addViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = transcripts)

        height <- convertHeight(page_coords$height,
            unitTo = get("page_units",
                envir = pgEnv
            ),
            valueOnly = TRUE
        )
        yscale <- strand_scale(
            strandSplit = transcriptsInternal$strandSplit,
            height = height
        )

        ## Make viewport for gene track
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = transcriptsInternal$xscale,
            yscale = yscale,
            just = transcriptsInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================

    backgroundGrob <- rectGrob(gp = gpar(
        fill = transcriptsInternal$bg,
        col = NA
    ), name = "background")
    assign("transcript_grobs", gTree(vp = vp, children = gList(backgroundGrob)),
        envir = pgEnv
    )

    # =========================================================================
    # DETERMINE ROWS FOR EACH ELEMENT
    # =========================================================================

    ## Determine how many rows are going to fit based on
    ## boxHeight, spaceHeight, and fontsize
    if (is.null(transcriptsInternal$labels)) {
        textHeight <- unit(0, "npc")
    } else {
        textHeight <- heightDetails(textGrob(
            label = "A",
            gp = gpar(fontsize = transcriptsInternal$fontsize)
        ))
    }


    if (is.null(transcripts$x) & is.null(transcripts$y)) {
        pushViewport(vp)
        boxHeight <- convertHeight(transcriptsInternal$boxHeight,
            unitTo = "npc", valueOnly = TRUE
        )
        spaceHeight <- boxHeight * (transcriptsInternal$spaceHeight)
        textHeight <- convertHeight(textHeight,
            unitTo = "npc",
            valueOnly = TRUE
        )
        upViewport()
        unit <- "npc"
    } else {
        boxHeight <- convertHeight(transcriptsInternal$boxHeight,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
        spaceHeight <- boxHeight * (transcriptsInternal$spaceHeight)
        textHeight <- convertHeight(textHeight,
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
        unit <- get("page_units", envir = pgEnv)
    }

    maxRows <- floor(height / (boxHeight + spaceHeight +
        textHeight + 0.25 * textHeight))
    wiggle <- abs(transcripts$chromend - transcripts$chromstart) *
        transcriptsInternal$spaceWidth

    if (transcriptsInternal$strandSplit == FALSE) {
        
        ## Get one representative row per transcript for ordering
        repData <- data[duplicated(data$TXNAME) == FALSE, ]
        
        ## Access default transcript prioritization
        repData <- defaultGenePriorities(
            data = repData,
            assembly = transcripts$assembly,
            transcript = TRUE
        )
        
        ## Assign rows
        rowData <- assignRows(data = repData[,c("TXSTART", "TXEND", "TXID")],
                            maxRows = maxRows,
                            wiggle = wiggle, 
                            rowCol = 3,
                            limitLabel = transcriptsInternal$limitLabel,
                            gTree = "transcript_grobs")
        
        ## Recombine row data with original data with all transcript details
        rowData <- suppressMessages(dplyr::left_join(
            x = data,
            y = rowData[,c("row", "TXID")],
            by = "TXID"
        ))
        
        ## Calculate y-coordinates
        rowData$y <- rowData$row * (boxHeight + spaceHeight + textHeight +
                                        0.25 * textHeight)
    } else {
        
        ## Get one rep set of data per unique transcript for ordering
        repPos <- posStrand[duplicated(posStrand$TXNAME) == FALSE, ]
        repMin <- minStrand[duplicated(minStrand$TXNAME) == FALSE, ]
        
        ## Order transcripts based on gene priorities
        repPos <- defaultGenePriorities(
            data = repPos,
            assembly = transcripts$assembly,
            geneHighlights = 
                transcriptsInternal$transcriptHighlights[,"transcript"],
            transcript = TRUE
        )
        repMin <- defaultGenePriorities(
            data = repMin,
            assembly = transcripts$assembly,
            geneHighlights = 
                transcriptsInternal$transcriptHighlights[,"transcript"],
            transcript = TRUE
        )
        
        ## Assign rows
        posData <- assignRows(data = repPos[,c("TXSTART", "TXEND", "TXID")],
                            maxRows = floor(maxRows / 2),
                            wiggle = wiggle, rowCol = 3,
                            limitLabel = transcriptsInternal$limitLabel,
                            gTree = "transcript_grobs")
        minData <- assignRows(data = repMin[,c("TXSTART", "TXEND", "TXID")],
                            maxRows = floor(maxRows / 2),
                            wiggle = wiggle, rowCol = 3,
                            limitLabel = transcriptsInternal$limitLabel,
                            gTree = "transcript_grobs",
                            side = "bottom")
        
        ## Recombine row data with original data
        posData <- suppressMessages(dplyr::left_join(
            x = posStrand,
            y = posData[,c("row", "TXID")],
            by = "TXID"
        ))
        minData <- suppressMessages(dplyr::left_join(
            x = minStrand,
            y = minData[,c("row", "TXID")],
            by = "TXID"
        ))
        
        ## Calculate y-coordinates
        posData$y <- (0.5 * spaceHeight) + posData$row *
            (boxHeight + spaceHeight + textHeight + 0.25 * textHeight)
        minData$y <- ((0.5 * spaceHeight + boxHeight + textHeight +
                        0.25 * textHeight) + minData$row *
                            (boxHeight + spaceHeight +
                                textHeight + 0.25 * textHeight)) * -1
        
        ## Combine pos and neg strand data into one 
        rowData <- rbind(posData, minData)
        
    }
    
    # =========================================================================
    # MAKE GROBS
    # =========================================================================

    if (nrow(rowData) > 0) {

        ##########################################################
        ## TRANSCRIPT LINES ONLY
        ##########################################################

        if ((transcripts$chromend - transcripts$chromstart) >= 25000000) {
            transcriptLine <- rectGrob(
                x = rowData$TXSTART,
                y = rowData$y,
                width = rowData$TXEND - rowData$TXSTART,
                height = boxHeight,
                just = c("left", "bottom"),
                gp = gpar(
                    fill = rowData$color,
                    col = rowData$color,
                    lwd = transcriptsInternal$stroke
                ),
                default.units = "native"
            )
        } else {

            ##########################################################
            ## TRANSCRIPT LINES
            ##########################################################

            transcriptLine <- rectGrob(
                x = rowData$TXSTART,
                y = rowData$y + 0.5 * boxHeight,
                width = rowData$TXEND - rowData$TXSTART,
                height = boxHeight * 0.2,
                just = "left",
                gp = gpar(
                    fill = rowData$color,
                    col = rowData$color,
                    lwd = transcriptsInternal$stroke
                ),
                default.units = "native"
            )

            ##########################################################
            ## TRANSCRIPT EXONS
            ##########################################################

            transcriptExons <- rectGrob(
                x = rowData$EXONSTART,
                y = rowData$y + 0.5 * boxHeight,
                width = rowData$EXONEND - rowData$EXONSTART,
                height = boxHeight * 0.65,
                just = "left",
                gp = gpar(fill = rowData$color, col = NA),
                default.units = "native"
            )
            assign("transcript_grobs",
                addGrob(get("transcript_grobs", envir = pgEnv),
                    child = transcriptExons
                ),
                envir = pgEnv
            )

            ##########################################################
            ## TRANSCRIPT CDS
            ##########################################################

            ## Get CDS regions that aren't NA
            cdsData <- rowData[which(!is.na(rowData$CDSSTART)), ]
            if (nrow(cdsData) > 0) {
                transcriptCds <- rectGrob(
                    x = cdsData$CDSSTART,
                    y = cdsData$y + 0.5 * boxHeight,
                    width = cdsData$CDSEND - cdsData$CDSSTART,
                    height = boxHeight,
                    just = "left",
                    gp = gpar(
                        fill = cdsData$color,
                        col = cdsData$color,
                        lwd = 1.25
                    ),
                    default.units = "native"
                )
                assign("transcript_grobs",
                    addGrob(get("transcript_grobs", envir = pgEnv),
                        child = transcriptCds
                    ),
                    envir = pgEnv
                )
            }
        }

        assign("transcript_grobs",
            addGrob(get("transcript_grobs", envir = pgEnv),
                child = transcriptLine
            ),
            envir = pgEnv
        )


        ##########################################################
        ## TRANSCRIPT NAME LABELS
        ##########################################################

        if (!is.null(transcriptsInternal$labels)) {

            ## Get representative transcript names
            transcriptNames <- rowData[duplicated(rowData$TXNAME) == FALSE, ]

            ## Add column with center location of each transcript label
            transcriptNames$labelLoc <- rowMeans(
                transcriptNames[c("TXSTART", "TXEND")]
            )

            if (transcriptsInternal$labels == "transcript") {
                label <- transcriptNames$TXNAME
            } else if (transcriptsInternal$labels == "gene") {
                label <- transcriptNames[[
                transcripts$assembly$display.column]]
            } else {
                label <- paste0(
                    transcriptNames[[
                    transcripts$assembly$display.column]],
                    ":", transcriptNames$TXNAME
                )
            }

            ## Get which names aren't cutoff
            checkedLabels <- apply(data.frame(
                "label" = label,
                "labelLoc" =
                    transcriptNames$labelLoc
            ),
            1, cutoffLabel,
            fontsize = transcriptsInternal$fontsize,
            xscale = transcriptsInternal$xscale,
            vp = vp, unit = unit
            )
            transcriptNames$label <- checkedLabels
            transcriptNames <- transcriptNames[!is.na(transcriptNames$label), ]

            if (nrow(transcriptNames) > 0) {
                transcript_names <- textGrob(
                    label = transcriptNames$label,
                    x = transcriptNames$labelLoc,
                    y = transcriptNames$y + boxHeight + textHeight * 0.25,
                    just = "bottom",
                    gp = gpar(
                        col = transcriptNames$color,
                        fontsize = transcriptsInternal$fontsize
                    ),
                    default.units = "native",
                    check.overlap = TRUE
                )

                assign("transcript_grobs",
                    addGrob(get("transcript_grobs", envir = pgEnv),
                        child = transcript_names
                    ),
                    envir = pgEnv
                )
            }
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (transcriptsInternal$draw == TRUE) {
        grid.draw(get("transcript_grobs", envir = pgEnv))
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    transcripts$grobs <- get("transcript_grobs", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("transcripts[", vp$name, "]")
    invisible(transcripts)
}
