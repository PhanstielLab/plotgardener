#' Plot a gene track for a specified genomic region
#' 
#' @usage bbPlotGenes(
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
#'     fontsize = 8,
#'     fontcolor = c("#669fd9", "#abcc8e"),
#'     fill = c("#669fd9", "#abcc8e"),
#'     geneOrder = NULL,
#'     geneHighlights = NULL,
#'     geneBackground = "grey",
#'     strandLabels = TRUE,
#'     stroke = 0.1,
#'     bg = NA,
#'     x = NULL,
#'     y = NULL,
#'     width = NULL,
#'     height = unit(0.6, "inches"),
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
#' \link[BentoBox]{bbAssembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param fontsize A numeric specifying text fontsize in points.
#' Default value is \code{fontsize = 8}.
#' @param fontcolor A character value or vector of length 2 indicating
#' the fontcolors for the plus strand and minus strand gene labels.
#' The first value will color the plus strand gene labels and
#' the second value will color the minus strand gene labels.
#' Default value is \code{fontcolor = c("#669fd9", "#abcc8e")}.
#' @param fill A character value or vector of length 2 indicating the
#' strand fill colors for the plus strand and minus strand plot elements.
#' The first value will color the plus strand plot elements and
#' the second label will color the minus strand plot elements.
#' Default value is \code{fill = c("#669fd9", "#abcc8e")}.
#' @param geneOrder An ordered character vector of gene names to
#' prioritize when labeling genes.
#' @param geneHighlights A two-column dataframe with a column named "gene"
#' containing gene names as strings to highlight and a named column "color"
#' containing corresponding highlight colors.
#' @param geneBackground If \code{geneHighlights} is given, a character
#' value indicating the color for genes that are not highlighted.
#' @param strandLabels A logical value indicating whether to include
#' + and - strand labels to the left of the gene track.
#' @param stroke A numeric value indicating the stroke width for gene
#' body outlines. Default value is \code{stroke = 0.1}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param x A numeric or unit object specifying genes plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying genes plot y-location.
#' The character value will
#' place the genes plot y relative to the bottom of the most recently
#' plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying genes plot width.
#' @param height A numeric or unit object specifying genes plot height.
#' @param just Justification of genes plot relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if
#' \code{x}, \code{y}, \code{width}, or \code{height} are only given
#' as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output
#' should be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bbParams} object
#' containing relevant function parameters.
#'
#' @return Returns a \code{bb_genes} object containing
#' relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load hg19 genomic annotation packages
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Set genomic coordinates
#' paramssmall <- bbParams(
#'     chrom = "chr8",
#'     chromstart = 1, chromend = 3000000,
#'     assembly = "hg19", width = 7
#' )
#' paramsbig <- bbParams(
#'     chrom = "chr8",
#'     chromstart = 1, chromend = 146364022,
#'     assembly = "hg19", width = 7
#' )
#' ## Set colors
#' cols <- c("#41B6C4", "#225EA8")
#'
#' ## Create page
#' bbPageCreate(width = 7.5, height = 3.5, default.units = "inches")
#'
#' ## Plot genes big
#' genesPlot <- bbPlotGenes(
#'     params = paramsbig, fill = cols,
#'     fontcolor = cols,
#'     x = 0.25, y = 0.25, height = 0.75,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bbAnnoGenomeLabel(
#'     plot = genesPlot, x = 0.25, y = 1.0,
#'     scale = "Mb", just = c("left", "top")
#' )
#'
#' ## Plot genes small
#' genesPlot <- bbPlotGenes(
#'     params = paramssmall,
#'     geneHighlights = data.frame(
#'         "gene" = c("DLGAP2"),
#'         "color" = c("#225EA8")
#'     ),
#'     geneBackground = "grey",
#'     x = 0.25, y = 2.25, height = 0.75,
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' bbAnnoGenomeLabel(
#'     plot = genesPlot, x = 0.25, y = 3.0, scale = "Mb",
#'     just = c("left", "top")
#' )
#'
#' ## Hide page guides
#' bbPageGuideHide()
#' @details
#' A gene track can be placed on a BentoBox coordinate page by providing
#' plot placement parameters:
#' \preformatted{
#' bbPlotGenes(chrom, chromstart = NULL, chromend = NULL,
#'             x, y, width, height, just = c("left", "top"),
#'             default.units = "inches")
#' }
#' This function can be used to quickly plot an unnannotated gene track
#' by ignoring plot placement parameters:
#' \preformatted{
#' bbPlotGenes(chrom, chromstart = NULL, chromend = NULL)
#' }
#'
#' Genomic annotation information is acquired through
#' \link[GenomicFeatures]{TxDb} and \link[AnnotationDbi]{OrgDb-class}
#' packages, as determined
#' through the \code{assembly} parameter. To avoid overcrowding of gene name
#' labels, plotted gene labels are by default prioritized according to
#' citation counts.
#'
#' @seealso \link[BentoBox]{bbAssembly},
#' \link[BentoBox]{bbGenomes}, \link[BentoBox]{bbDefaultPackages}
#'
#' @export
bbPlotGenes <- function(chrom, chromstart = NULL, chromend = NULL,
                        assembly = "hg38", fontsize = 8,
                        fontcolor = c("#669fd9", "#abcc8e"),
                        fill = c("#669fd9", "#abcc8e"),
                        geneOrder = NULL, geneHighlights = NULL,
                        geneBackground = "grey", strandLabels = TRUE,
                        stroke = 0.1, bg = NA, x = NULL, y = NULL,
                        width = NULL, height = unit(0.6, "inches"),
                        just = c("left", "top"), default.units = "inches",
                        draw = TRUE, params = NULL) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that checks errors for bbPlotGenes
    errorcheck_bbPlotGenes <- function(chromstart, chromend, fill) {
        
        ## Genomic region errors
        bb_regionErrors(chromstart = chromstart, chromend = chromend)
        
        ## Make sure fill isn't colorby
        bb_checkColorby(fill = fill,
                        colorby = FALSE)
    }

    ## Define a function that gets total lengths and centers of visible genes
    get_geneCenters <- function(geneData, object) {
        minStart <- min(geneData$TXSTART)
        if (minStart < object$chromstart){
            minStart <- object$chromstart
        }
        maxEnd <- max(geneData$TXEND)
        if (maxEnd > object$chromend){
            maxEnd <- object$chromend
        }
        
        geneLength <- maxEnd - minStart
        geneCenter <- mean(c(minStart, maxEnd))
        return(list(geneLength, geneCenter))
    }

    ## Define a function that prioritizes genes based on internal
    ## priorities and/or user-defined geneOrder
    gene_priorities <- function(genes, geneOrder, displayCol) {

        ## Split list into the ones in geneOrder and the ones not
        subset <- genes[which(genes[[displayCol]] %in% geneOrder), ]
        remaining <- genes[which(!genes[[displayCol]] %in% geneOrder), ]

        ## Put the geneOrder subset into the same order
        subset <- subset[match(geneOrder, subset[[displayCol]]), ]
        subset <- subset[which(!is.na(subset[[displayCol]])), ]

        ## Recombine lists, putting geneOrder section at the top
        combined <- rbind(subset, remaining)

        return(combined)
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    bb_genesInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "bb_genesInternal"
    )
    
    ## Justification
    bb_genesInternal$just <- bb_justConversion(just = bb_genesInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    bb_genes <- structure(list(
        chrom = bb_genesInternal$chrom,
        chromstart = bb_genesInternal$chromstart,
        chromend = bb_genesInternal$chromend,
        assembly = bb_genesInternal$assembly,
        x = bb_genesInternal$x, y = bb_genesInternal$y,
        width = bb_genesInternal$width,
        height = bb_genesInternal$height,
        just = bb_genesInternal$just, grobs = NULL
    ),
    class = "bb_genes"
    )
    attr(x = bb_genes, which = "plotted") <- bb_genesInternal$draw

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(bb_genes$chrom)) stop("argument \"chrom\" is missing, ",
                                    "with no default.", call. = FALSE)
    errorcheck_bbPlotGenes(
        chromstart = bb_genes$chromstart,
        chromend = bb_genes$chromend,
        fill = bb_genesInternal$fill
    )
    check_placement(object = bb_genes)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    bb_genes$assembly <- parse_bbAssembly(assembly = bb_genes$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    bb_genes <- defaultUnits(
        object = bb_genes,
        default.units = bb_genesInternal$default.units
    )

    # =========================================================================
    # GET APPROPRIATE BUILD DATA
    # =========================================================================

    buildData <- geneData(object = bb_genes,
                        objectInternal = bb_genesInternal)
    bb_genes <- buildData[[1]]
    bb_genesInternal <- buildData[[2]]
    displayCol <- bb_genesInternal$displayCol

    # =========================================================================
    # SUBSET DATA BY STRAND
    # =========================================================================
    data <- bb_genesInternal$data
    ## Genes on plus strand and genes on minus strand
    plus_genes <- data[which(data$TXSTRAND == "+"), ]
    minus_genes <- data[which(data$TXSTRAND == "-"), ]

    # =========================================================================
    # COLORS AND HIGHLIGHTS
    # =========================================================================

    ## If fontcolor and fill are length 1, change to length 2
    if (length(bb_genesInternal$fontcolor) == 1) {
        bb_genesInternal$fontcolor <- rep(bb_genesInternal$fontcolor, 2)
    }
    if (length(bb_genesInternal$fill) == 1) {
        bb_genesInternal$fill <- rep(bb_genesInternal$fill, 2)
    }


    ## Add strandColor and fontColors
    if (nrow(plus_genes) > 0) {
        plus_genes$strandColor <- bb_genesInternal$fill[1]
        plus_genes$fontColor <- bb_genesInternal$fontcolor[1]
    }

    if (nrow(minus_genes) > 0) {
        minus_genes$strandColor <- bb_genesInternal$fill[2]
        minus_genes$fontColor <- bb_genesInternal$fontcolor[2]
    }


    ## Find genes to be highlighted
    if (!is.null(bb_genesInternal$geneHighlights)) {
        colnames(bb_genesInternal$geneHighlights) <- c(
            displayCol,
            "strandColor"
        )
        plus_genes$strandColor <- NULL
        minus_genes$strandColor <- NULL

        plusHighlight <- plus_genes[which(plus_genes[[displayCol]] %in%
            bb_genesInternal$geneHighlights
            [, displayCol]), ]
        minusHighlight <- minus_genes[which(minus_genes[[displayCol]] %in%
            bb_genesInternal$geneHighlights
            [, displayCol]), ]
        plusBackground <- plus_genes[which(!plus_genes[[displayCol]] %in%
            bb_genesInternal$geneHighlights
            [, displayCol]), ]
        minusBackground <- minus_genes[which(!minus_genes[[displayCol]] %in%
            bb_genesInternal$geneHighlights
            [, displayCol]), ]

        # Change highlight genes to highlight color
        plusHighlight <- merge(plusHighlight, bb_genesInternal$geneHighlights,
            by = displayCol
        )
        plusHighlight$fontColor <- plusHighlight$strandColor

        minusHighlight <- merge(minusHighlight, bb_genesInternal$geneHighlights,
            by = displayCol
        )
        minusHighlight$fontColor <- minusHighlight$strandColor

        # Change background genes to background color
        plusBackground$strandColor <- rep(
            bb_genesInternal$geneBackground[1],
            nrow(plusBackground)
        )
        plusBackground$fontColor <- rep(
            bb_genesInternal$geneBackground[1],
            nrow(plusBackground)
        )
        minusBackground$strandColor <- rep(
            bb_genesInternal$geneBackground[1],
            nrow(minusBackground)
        )
        minusBackground$fontColor <- rep(
            bb_genesInternal$geneBackground[1],
            nrow(minusBackground)
        )

        ## Automatically change fontcolor for +/- strand label
        bb_genesInternal$fontcolor[1] <- bb_genesInternal$geneBackground[1]
        bb_genesInternal$fontcolor[2] <- bb_genesInternal$geneBackground[1]

        ## Put highlight and background genes back together
        plus_genes <- rbind(plusHighlight, plusBackground)
        minus_genes <- rbind(minusHighlight, minusBackground)
    }

    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Define text grob for "+" and "-" viewport scaling
    tG <- textGrob(label = "+", gp = gpar(
        fontsize =
            bb_genesInternal$fontsize
    ))

    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "bb_genes",
        length(intersect(
            grep(
                pattern = "bb_genes",
                x = currentViewports
            ),
            grep(
                pattern = "label",
                x = currentViewports,
                invert = TRUE
            )
        )) + 1
    )

    if (is.null(bb_genes$x) | is.null(bb_genes$y)) {
        if (bb_genesInternal$strandLabels == TRUE) {

            ## Make viewport for "+" and "-" labels to the left of genetrack
            vp_labelW <- convertWidth(widthDetails(tG) * 2, unitTo = "npc")

            vp_label <- viewport(
                height = unit(.12, "npc"),
                width = vp_labelW,
                x = unit(0, "npc"), y = unit(0.5, "npc"),
                just = "left",
                name = paste0(vp_name, "_label")
            )

            vp_gene <- viewport(
                height = unit(.12, "npc"),
                width = unit(1, "npc") - vp_labelW,
                x = vp_labelW, y = unit(0.5, "npc"),
                clip = "on",
                xscale = bb_genesInternal$xscale,
                just = "left",
                name = vp_name
            )
        } else {
            vp_gene <- viewport(
                height = unit(.12, "npc"),
                width = unit(1, "npc"),
                x = unit(0, "npc"), y = unit(0.5, "npc"),
                clip = "on",
                xscale = bb_genesInternal$xscale,
                just = "left",
                name = vp_name
            )
        }


        if (bb_genesInternal$draw == TRUE) {
            vp_gene$name <- "bb_genes1"
            if (bb_genesInternal$strandLabels == TRUE) {
                vp_label$name <- "bb_genes1_label"
            }
            grid.newpage()
        }
    } else {
        add_bbViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = bb_genes)

        ## Make viewport for gene track
        vp_gene <- viewport(
            height = page_coords$height,
            width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = bb_genesInternal$xscale,
            just = bb_genesInternal$just,
            name = vp_name
        )

        if (bb_genesInternal$strandLabels == TRUE) {

            ## Make viewport for "+" and "-" labels based on above viewport
            topLeft_vp <- vp_topLeft(viewport = vp_gene)

            vp_label <- viewport(
                height = page_coords$height,
                width = convertWidth(widthDetails(tG) * 2,
                    unitTo = get("page_units",
                        envir = bbEnv
                    )
                ),
                x = topLeft_vp[[1]], y = topLeft_vp[[2]],
                just = c("right", "top"),
                name = paste0(vp_name, "_label")
            )
        }
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================

    backgroundGrob <- rectGrob(
        gp = gpar(
            fill = bb_genesInternal$bg,
            col = NA
        ),
        name = "background", vp = vp_gene
    )
    assign("gene_grobs", gTree(children = gList(backgroundGrob)),
        envir = bbEnv
    )

    # =========================================================================
    # MAKE GROBS
    # =========================================================================

    if (nrow(data) > 0) {

        ##########################################################
        ## GENE LINES ONLY
        ##########################################################

        if ((bb_genes$chromend - bb_genes$chromstart) >= 25000000) {
            if (nrow(plus_genes) > 0) {
                plus_geneGrobs <- rectGrob(
                    x = plus_genes$TXSTART,
                    y = unit(0.63, "npc"),
                    width = plus_genes$TXEND - plus_genes$TXSTART,
                    height = unit(0.18, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = plus_genes$strandColor,
                        col = plus_genes$strandColor,
                        lwd = bb_genesInternal$stroke
                    ),
                    vp = vp_gene, default.units = "native"
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = bbEnv),
                        child = plus_geneGrobs
                    ),
                    envir = bbEnv
                )
            }

            if (nrow(minus_genes) > 0) {
                minus_geneGrobs <- rectGrob(
                    x = minus_genes$TXSTART,
                    y = unit(0.37, "npc"),
                    width = minus_genes$TXEND - minus_genes$TXSTART,
                    height = unit(0.18, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = minus_genes$strandColor,
                        col = minus_genes$strandColor,
                        lwd = bb_genesInternal$stroke
                    ),
                    vp = vp_gene, default.units = "native"
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = bbEnv),
                        child = minus_geneGrobs
                    ),
                    envir = bbEnv
                )
            }
        } else {
            if (nrow(plus_genes) > 0) {
                ##########################################################
                ## GENE LINES
                ##########################################################

                plus_geneGrobs <- rectGrob(
                    x = plus_genes$TXSTART,
                    y = unit(0.63, "npc"),
                    width = plus_genes$TXEND - plus_genes$TXSTART,
                    height = unit(0.05, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = plus_genes$strandColor,
                        col = plus_genes$strandColor,
                        lwd = bb_genesInternal$stroke
                    ),
                    vp = vp_gene, default.units = "native"
                )

                ##########################################################
                ## GENE EXONS
                ###########################################################

                plus_exonsGrobs <- rectGrob(
                    x = plus_genes$EXONSTART,
                    y = unit(0.63, "npc"),
                    width = plus_genes$EXONEND - plus_genes$EXONSTART,
                    height = unit(0.1, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = plus_genes$strandColor,
                        col = NA
                    ),
                    vp = vp_gene, default.units = "native"
                )


                ##########################################################
                ## GENE CDS
                ##########################################################

                ## Get CDS regions that aren't NA
                poscdsData <- plus_genes[which(!is.na(plus_genes$CDSSTART)), ]
                if (nrow(poscdsData) > 0) {
                    plus_cdsGrobs <- rectGrob(
                        x = poscdsData$CDSSTART,
                        y = unit(0.63, "npc"),
                        width = poscdsData$CDSEND - poscdsData$CDSSTART,
                        height = unit(0.18, "npc"),
                        just = "left",
                        gp = gpar(
                            fill = poscdsData$strandColor,
                            col = poscdsData$strandColor,
                            lwd = 1.25
                        ),
                        vp = vp_gene, default.units = "native"
                    )
                    assign("gene_grobs",
                        addGrob(get("gene_grobs", envir = bbEnv),
                            child = plus_cdsGrobs
                        ),
                        envir = bbEnv
                    )
                }

                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = bbEnv),
                        child = plus_geneGrobs
                    ),
                    envir = bbEnv
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = bbEnv),
                        child = plus_exonsGrobs
                    ),
                    envir = bbEnv
                )
            }

            if (nrow(minus_genes) > 0) {

                ##########################################################
                ## GENE LINES
                ##########################################################

                minus_geneGrobs <- rectGrob(
                    x = minus_genes$TXSTART,
                    y = unit(0.37, "npc"),
                    width = minus_genes$TXEND - minus_genes$TXSTART,
                    height = unit(0.05, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = minus_genes$strandColor,
                        col = minus_genes$strandColor,
                        lwd = bb_genesInternal$stroke
                    ),
                    vp = vp_gene, default.units = "native"
                )

                ##########################################################
                ## GENE EXONS
                ###########################################################

                minus_exonsGrobs <- rectGrob(
                    x = minus_genes$EXONSTART,
                    y = unit(0.37, "npc"),
                    width = minus_genes$EXONEND - minus_genes$EXONSTART,
                    height = unit(0.1, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = minus_genes$strandColor,
                        col = NA
                    ),
                    vp = vp_gene, default.units = "native"
                )

                ##########################################################
                ## GENE CDS
                ##########################################################

                ## Get CDS regions that aren't NA
                mincdsData <- minus_genes[which(!is.na(minus_genes$CDSSTART)), ]
                if (nrow(mincdsData) > 0) {
                    minus_cdsGrobs <- rectGrob(
                        x = mincdsData$CDSSTART,
                        y = unit(0.37, "npc"),
                        width = mincdsData$CDSEND - mincdsData$CDSSTART,
                        height = unit(0.18, "npc"),
                        just = "left",
                        gp = gpar(
                            fill = mincdsData$strandColor,
                            col = mincdsData$strandColor,
                            lwd = 1.25
                        ),
                        vp = vp_gene, default.units = "native"
                    )
                    assign("gene_grobs",
                        addGrob(get("gene_grobs", envir = bbEnv),
                            child = minus_cdsGrobs
                        ),
                        envir = bbEnv
                    )
                }

                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = bbEnv),
                        child = minus_geneGrobs
                    ),
                    envir = bbEnv
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = bbEnv),
                        child = minus_exonsGrobs
                    ),
                    envir = bbEnv
                )
            }
        }

        ##########################################################
        ## GENE NAME LABELS
        ##########################################################

        ## Split into mini dataframes based on GENEID
        separatedplusGenes <- split(plus_genes, plus_genes$GENEID)
        separatedminusGenes <- split(minus_genes, minus_genes$GENEID)

        ## For every GENEID dataframe, get length and center
        ## location of each total gene
        plusgeneCenters <- lapply(separatedplusGenes, get_geneCenters,
                                object = bb_genes)
        plusgeneCenters <- data.frame(
            GENEID = names(plusgeneCenters),
            length = unlist(purrr::map(plusgeneCenters, 1)),
            labelLoc = unlist(purrr::map(plusgeneCenters, 2))
        )
        minusgeneCenters <- lapply(separatedminusGenes, get_geneCenters, 
                                object = bb_genes)
        minusgeneCenters <- data.frame(
            GENEID = names(minusgeneCenters),
            length = unlist(purrr::map(minusgeneCenters, 1)),
            labelLoc = unlist(purrr::map(minusgeneCenters, 2))
        )

        ## Get representative value of each shown gene name
        plusgeneNames <- plus_genes[duplicated(plus_genes[[displayCol]]) ==
            FALSE, ]
        minusgeneNames <- minus_genes[duplicated(minus_genes[[displayCol]]) ==
            FALSE, ]

        ## Add column with center location of each gene label
        plusgeneNames <- suppressMessages(dplyr::left_join(
            x = plusgeneNames,
            y = plusgeneCenters,
            by = "GENEID"
        ))
        minusgeneNames <- suppressMessages(dplyr::left_join(
            x = minusgeneNames,
            y = minusgeneCenters,
            by = "GENEID"
        ))

        ## Add all the gene names in the region to object
        geneNames <- c(
            plusgeneNames[[displayCol]],
            minusgeneNames[[displayCol]]
        )
        bb_genes$genes <- geneNames

        ## DECLUTTER LABELS

        ## Put genes in order of default gene prioritization
        ## based on citation/gene length
        plusgeneNames <- defaultGenePriorities(
            data = plusgeneNames,
            assembly = bb_genes$assembly
        )
        minusgeneNames <- defaultGenePriorities(
            data = minusgeneNames,
            assembly = bb_genes$assembly
        )

        if (!is.null(bb_genesInternal$geneOrder)) {

            ## Integrate geneOrder and default prioritization
            plusgeneNames <- gene_priorities(
                genes = plusgeneNames,
                geneOrder = bb_genesInternal$geneOrder,
                displayCol = displayCol
            )
            minusgeneNames <- gene_priorities(
                genes = minusgeneNames,
                geneOrder = bb_genesInternal$geneOrder,
                displayCol = displayCol
            )
        }

        ## Get which gene names aren't cutoff
        if (nrow(plusgeneNames) > 0) {
            checkedplusLabels <- apply(data.frame(
                "label" = plusgeneNames[[displayCol]],
                "labelLoc" = plusgeneNames$labelLoc
            ),
            1, cutoffLabel,
            fontsize = bb_genesInternal$fontsize,
            xscale = bb_genesInternal$xscale,
            vp = vp_gene,
            unit = "inches"
            )
            plusgeneNames[[displayCol]] <- checkedplusLabels
            plusgeneNames <- plusgeneNames[!is.na(
                plusgeneNames[[displayCol]]
            ), ]

            if (nrow(plusgeneNames) > 0) {
                plus_names <- textGrob(
                    label = plusgeneNames[[displayCol]],
                    x = plusgeneNames$labelLoc,
                    y = unit(0.85, "npc"),
                    gp = gpar(
                        col = plusgeneNames$fontColor,
                        fontsize = bb_genesInternal$fontsize
                    ),
                    vp = vp_gene,
                    default.units = "native",
                    check.overlap = TRUE
                )

                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = bbEnv),
                        child = plus_names
                    ),
                    envir = bbEnv
                )
            }
        }

        if (nrow(minusgeneNames) > 0) {
            checkedminusLabels <- apply(data.frame(
                "label" = minusgeneNames[[displayCol]],
                "labelLoc" = minusgeneNames$labelLoc
            ),
            1, cutoffLabel,
            fontsize = bb_genesInternal$fontsize,
            xscale = bb_genesInternal$xscale,
            vp = vp_gene,
            unit = "inches"
            )
            minusgeneNames[[displayCol]] <- checkedminusLabels
            minusgeneNames <- minusgeneNames[!is.na(
                minusgeneNames[[displayCol]]
            ), ]

            if (nrow(minusgeneNames) > 0) {
                minus_names <- textGrob(
                    label = minusgeneNames[[displayCol]],
                    x = minusgeneNames$labelLoc,
                    y = unit(0.15, "npc"),
                    gp = gpar(
                        col = minusgeneNames$fontColor,
                        fontsize = bb_genesInternal$fontsize
                    ),
                    vp = vp_gene,
                    default.units = "native",
                    check.overlap = TRUE
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = bbEnv),
                        child = minus_names
                    ),
                    envir = bbEnv
                )
            }
        }

        ##########################################################
        ## + AND - LABELS
        ##########################################################

        if (bb_genesInternal$strandLabels == TRUE) {
            plus_label <- textGrob(
                label = "+", y = unit(0.63, "npc"),
                gp = gpar(
                    fontsize = bb_genesInternal$fontsize + 2,
                    col = bb_genesInternal$fontcolor[1],
                    fontface = "bold"
                ), vp = vp_label
            )
            minus_label <- textGrob(
                label = "-", y = unit(0.37, "npc"),
                gp = gpar(
                    fontsize = bb_genesInternal$fontsize + 2,
                    col = bb_genesInternal$fontcolor[2],
                    fontface = "bold"
                ), vp = vp_label
            )

            assign("gene_grobs",
                addGrob(get("gene_grobs", envir = bbEnv),
                    child = plus_label
                ),
                envir = bbEnv
            )
            assign("gene_grobs",
                addGrob(get("gene_grobs", envir = bbEnv),
                    child = minus_label
                ),
                envir = bbEnv
            )
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (bb_genesInternal$draw == TRUE) {
        grid.draw(get("gene_grobs", envir = bbEnv))
    }


    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    bb_genes$grobs <- get("gene_grobs", envir = bbEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("bb_genes[", vp_gene$name, "]")
    invisible(bb_genes)
}

## Define a function that removes gene and transcript
## name labels that will be cutoff
## @param df data.frame of genes in region
## @param fontsize fontsize of labels
## @param xscale vector of genomic region of viewport
## @param vp associated viewport where genes/transcripts are plotted
## @param unit for bbPlotTranscripts, unit indicator
cutoffLabel <- function(df, fontsize, xscale, vp, unit) {
    label <- df[1]
    location <- utils::type.convert(df[2], as.is = TRUE)
    
    ## Update viewport fontsize for proper text size calculation
    vp$gp <- gpar(fontsize = fontsize)
    
    if (unit == "npc"){
        downViewport(name = vp$name)
        labelWidth <- convertWidth(widthDetails(textGrob(
            label = label,
            gp = gpar(fontsize = fontsize)
        )),
        unitTo = "native", valueOnly = TRUE
        )
        upViewport()
    } else {
        pushViewport(vp)
        labelWidth <- convertWidth(widthDetails(textGrob(
            label = label,
            gp = gpar(fontsize = fontsize)
        )),
        unitTo = "native", valueOnly = TRUE
        )
        upViewport()
    }
    
    leftBound <- location - 0.5 * labelWidth
    rightBound <- location + 0.5 * labelWidth
    
    if (leftBound < xscale[1] | rightBound > xscale[2]) {
        return(NA)
    } else {
        return(label)
    }
}
