#' Plot a gene track for a specified genomic region
#' 
#' @usage plotGenes(
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
#' \link[plotgardener]{assembly} object.
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
#' plotted plot according to the units of the plotgardener page.
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
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#'
#' @return Returns a \code{genes} object containing
#' relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load hg19 genomic annotation packages
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Set genomic coordinates
#' paramssmall <- pgParams(
#'     chrom = "chr8",
#'     chromstart = 1, chromend = 3000000,
#'     assembly = "hg19", width = 7
#' )
#' paramsbig <- pgParams(
#'     chrom = "chr8",
#'     chromstart = 1, chromend = 146364022,
#'     assembly = "hg19", width = 7
#' )
#' ## Set colors
#' cols <- c("#41B6C4", "#225EA8")
#'
#' ## Create page
#' pageCreate(width = 7.5, height = 3.5, default.units = "inches")
#'
#' ## Plot genes big
#' genesPlot <- plotGenes(
#'     params = paramsbig, fill = cols,
#'     fontcolor = cols,
#'     x = 0.25, y = 0.25, height = 0.75,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Annotate genome label
#' annoGenomeLabel(
#'     plot = genesPlot, x = 0.25, y = 1.0,
#'     scale = "Mb", just = c("left", "top")
#' )
#'
#' ## Plot genes small
#' genesPlot <- plotGenes(
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
#' annoGenomeLabel(
#'     plot = genesPlot, x = 0.25, y = 3.0, scale = "Mb",
#'     just = c("left", "top")
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' A gene track can be placed on a page by providing
#' plot placement parameters:
#' \preformatted{
#' plotGenes(chrom, chromstart = NULL, chromend = NULL,
#'             x, y, width, height, just = c("left", "top"),
#'             default.units = "inches")
#' }
#' This function can be used to quickly plot an unnannotated gene track
#' by ignoring plot placement parameters:
#' \preformatted{
#' plotGenes(chrom, chromstart = NULL, chromend = NULL)
#' }
#'
#' Genomic annotation information is acquired through
#' \link[GenomicFeatures]{TxDb} and \link[AnnotationDbi]{OrgDb-class}
#' packages, as determined
#' through the \code{assembly} parameter. To avoid overcrowding of gene name
#' labels, plotted gene labels are by default prioritized according to
#' citation counts.
#'
#' @seealso \link[plotgardener]{assembly},
#' \link[plotgardener]{genomes}, \link[plotgardener]{defaultPackages}
#'
#' @export
plotGenes <- function(chrom, chromstart = NULL, chromend = NULL,
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

    ## Define a function that checks errors for plotGenes
    errorcheck_plotGenes <- function(chromstart, chromend, fill) {
        
        ## Genomic region errors
        regionErrors(chromstart = chromstart, chromend = chromend)
        
        ## Make sure fill isn't colorby
        checkColorby(fill = fill,
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

    genesInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "genesInternal"
    )
    
    ## Justification
    genesInternal$just <- justConversion(just = genesInternal$just)

    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    genes <- structure(list(
        chrom = genesInternal$chrom,
        chromstart = genesInternal$chromstart,
        chromend = genesInternal$chromend,
        assembly = genesInternal$assembly,
        x = genesInternal$x, y = genesInternal$y,
        width = genesInternal$width,
        height = genesInternal$height,
        just = genesInternal$just, grobs = NULL
    ),
    class = "genes"
    )
    attr(x = genes, which = "plotted") <- genesInternal$draw

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    if (is.null(genes$chrom)) stop("argument \"chrom\" is missing, ",
                                    "with no default.", call. = FALSE)
    errorcheck_plotGenes(
        chromstart = genes$chromstart,
        chromend = genes$chromend,
        fill = genesInternal$fill
    )
    check_placement(object = genes)

    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================

    genes$assembly <- parseAssembly(assembly = genes$assembly)

    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    genes <- defaultUnits(
        object = genes,
        default.units = genesInternal$default.units
    )

    # =========================================================================
    # GET APPROPRIATE BUILD DATA
    # =========================================================================

    buildData <- geneData(object = genes,
                        objectInternal = genesInternal)
    genes <- buildData[[1]]
    genesInternal <- buildData[[2]]
    displayCol <- genesInternal$displayCol

    # =========================================================================
    # SUBSET DATA BY STRAND
    # =========================================================================
    data <- genesInternal$data
    ## Genes on plus strand and genes on minus strand
    plus_genes <- data[which(data$TXSTRAND == "+"), ]
    minus_genes <- data[which(data$TXSTRAND == "-"), ]

    # =========================================================================
    # COLORS AND HIGHLIGHTS
    # =========================================================================

    ## If fontcolor and fill are length 1, change to length 2
    if (length(genesInternal$fontcolor) == 1) {
        genesInternal$fontcolor <- rep(genesInternal$fontcolor, 2)
    }
    if (length(genesInternal$fill) == 1) {
        genesInternal$fill <- rep(genesInternal$fill, 2)
    }


    ## Add strandColor and fontColors
    if (nrow(plus_genes) > 0) {
        plus_genes$strandColor <- genesInternal$fill[1]
        plus_genes$fontColor <- genesInternal$fontcolor[1]
    }

    if (nrow(minus_genes) > 0) {
        minus_genes$strandColor <- genesInternal$fill[2]
        minus_genes$fontColor <- genesInternal$fontcolor[2]
    }


    ## Find genes to be highlighted
    if (!is.null(genesInternal$geneHighlights)) {
        colnames(genesInternal$geneHighlights) <- c(
            displayCol,
            "strandColor"
        )
        plus_genes$strandColor <- NULL
        minus_genes$strandColor <- NULL

        plusHighlight <- plus_genes[which(plus_genes[[displayCol]] %in%
            genesInternal$geneHighlights
            [, displayCol]), ]
        minusHighlight <- minus_genes[which(minus_genes[[displayCol]] %in%
            genesInternal$geneHighlights
            [, displayCol]), ]
        plusBackground <- plus_genes[which(!plus_genes[[displayCol]] %in%
            genesInternal$geneHighlights
            [, displayCol]), ]
        minusBackground <- minus_genes[which(!minus_genes[[displayCol]] %in%
            genesInternal$geneHighlights
            [, displayCol]), ]

        # Change highlight genes to highlight color
        plusHighlight <- merge(plusHighlight, genesInternal$geneHighlights,
            by = displayCol
        )
        plusHighlight$fontColor <- plusHighlight$strandColor

        minusHighlight <- merge(minusHighlight, genesInternal$geneHighlights,
            by = displayCol
        )
        minusHighlight$fontColor <- minusHighlight$strandColor

        # Change background genes to background color
        plusBackground$strandColor <- rep(
            genesInternal$geneBackground[1],
            nrow(plusBackground)
        )
        plusBackground$fontColor <- rep(
            genesInternal$geneBackground[1],
            nrow(plusBackground)
        )
        minusBackground$strandColor <- rep(
            genesInternal$geneBackground[1],
            nrow(minusBackground)
        )
        minusBackground$fontColor <- rep(
            genesInternal$geneBackground[1],
            nrow(minusBackground)
        )

        ## Automatically change fontcolor for +/- strand label
        genesInternal$fontcolor[1] <- genesInternal$geneBackground[1]
        genesInternal$fontcolor[2] <- genesInternal$geneBackground[1]

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
            genesInternal$fontsize
    ))

    ## Name viewport
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "genes",
        length(intersect(
            grep(
                pattern = "genes",
                x = currentViewports
            ),
            grep(
                pattern = "label",
                x = currentViewports,
                invert = TRUE
            )
        )) + 1
    )

    if (is.null(genes$x) | is.null(genes$y)) {
        if (genesInternal$strandLabels == TRUE) {

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
                xscale = genesInternal$xscale,
                just = "left",
                name = vp_name
            )
        } else {
            vp_gene <- viewport(
                height = unit(.12, "npc"),
                width = unit(1, "npc"),
                x = unit(0, "npc"), y = unit(0.5, "npc"),
                clip = "on",
                xscale = genesInternal$xscale,
                just = "left",
                name = vp_name
            )
        }


        if (genesInternal$draw == TRUE) {
            vp_gene$name <- "genes1"
            if (genesInternal$strandLabels == TRUE) {
                vp_label$name <- "genes1_label"
            }
            grid.newpage()
        }
    } else {
        addViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = genes)

        ## Make viewport for gene track
        vp_gene <- viewport(
            height = page_coords$height,
            width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = genesInternal$xscale,
            just = genesInternal$just,
            name = vp_name
        )

        if (genesInternal$strandLabels == TRUE) {

            ## Make viewport for "+" and "-" labels based on above viewport
            topLeft_vp <- vp_topLeft(viewport = vp_gene)

            vp_label <- viewport(
                height = page_coords$height,
                width = convertWidth(widthDetails(tG) * 2,
                    unitTo = get("page_units",
                        envir = pgEnv
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
            fill = genesInternal$bg,
            col = NA
        ),
        name = "background", vp = vp_gene
    )
    assign("gene_grobs", gTree(children = gList(backgroundGrob)),
        envir = pgEnv
    )

    # =========================================================================
    # MAKE GROBS
    # =========================================================================

    if (nrow(data) > 0) {

        ##########################################################
        ## GENE LINES ONLY
        ##########################################################

        if ((genes$chromend - genes$chromstart) >= 25000000) {
            if (nrow(plus_genes) > 0) {
                plus_lines <- unique(plus_genes[,c("TXSTART", "TXEND", 
                                                "fontColor", "strandColor")])
                plus_geneGrobs <- rectGrob(
                    x = plus_lines$TXSTART,
                    y = unit(0.63, "npc"),
                    width = plus_lines$TXEND - plus_lines$TXSTART,
                    height = unit(0.18, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = plus_lines$strandColor,
                        col = plus_lines$strandColor,
                        lwd = genesInternal$stroke
                    ),
                    vp = vp_gene, default.units = "native"
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = pgEnv),
                        child = plus_geneGrobs
                    ),
                    envir = pgEnv
                )
            }

            if (nrow(minus_genes) > 0) {
                minus_lines <- unique(minus_genes[,c("TXSTART", "TXEND", 
                                                "fontColor", "strandColor")])
                minus_geneGrobs <- rectGrob(
                    x = minus_lines$TXSTART,
                    y = unit(0.37, "npc"),
                    width = minus_lines$TXEND - minus_lines$TXSTART,
                    height = unit(0.18, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = minus_lines$strandColor,
                        col = minus_lines$strandColor,
                        lwd = genesInternal$stroke
                    ),
                    vp = vp_gene, default.units = "native"
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = pgEnv),
                        child = minus_geneGrobs
                    ),
                    envir = pgEnv
                )
            }
        } else {
            if (nrow(plus_genes) > 0) {
                ##########################################################
                ## GENE LINES
                ##########################################################
                plus_lines <- unique(plus_genes[,c("TXSTART", "TXEND", 
                                                "fontColor", "strandColor")])
                plus_geneGrobs <- rectGrob(
                    x = plus_lines$TXSTART,
                    y = unit(0.63, "npc"),
                    width = plus_lines$TXEND - plus_lines$TXSTART,
                    height = unit(0.05, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = plus_lines$strandColor,
                        col = plus_lines$strandColor,
                        lwd = genesInternal$stroke
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
                        addGrob(get("gene_grobs", envir = pgEnv),
                            child = plus_cdsGrobs
                        ),
                        envir = pgEnv
                    )
                }

                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = pgEnv),
                        child = plus_geneGrobs
                    ),
                    envir = pgEnv
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = pgEnv),
                        child = plus_exonsGrobs
                    ),
                    envir = pgEnv
                )
            }

            if (nrow(minus_genes) > 0) {

                ##########################################################
                ## GENE LINES
                ##########################################################
                minus_lines <- unique(minus_genes[,c("TXSTART", "TXEND", 
                                                "fontColor", "strandColor")])
                minus_geneGrobs <- rectGrob(
                    x = minus_lines$TXSTART,
                    y = unit(0.37, "npc"),
                    width = minus_lines$TXEND - minus_lines$TXSTART,
                    height = unit(0.05, "npc"),
                    just = "left",
                    gp = gpar(
                        fill = minus_lines$strandColor,
                        col = minus_lines$strandColor,
                        lwd = genesInternal$stroke
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
                        addGrob(get("gene_grobs", envir = pgEnv),
                            child = minus_cdsGrobs
                        ),
                        envir = pgEnv
                    )
                }

                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = pgEnv),
                        child = minus_geneGrobs
                    ),
                    envir = pgEnv
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = pgEnv),
                        child = minus_exonsGrobs
                    ),
                    envir = pgEnv
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
                                object = genes)
        plusgeneCenters <- data.frame(
            GENEID = names(plusgeneCenters),
            length = unlist(purrr::map(plusgeneCenters, 1)),
            labelLoc = unlist(purrr::map(plusgeneCenters, 2))
        )
        minusgeneCenters <- lapply(separatedminusGenes, get_geneCenters, 
                                object = genes)
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
        genes$genes <- geneNames

        ## DECLUTTER LABELS

        ## Put genes in order of default gene prioritization
        ## based on citation/gene length
        plusgeneNames <- defaultGenePriorities(
            data = plusgeneNames,
            assembly = genes$assembly,
            geneHighlights = genesInternal$geneHighlights[[displayCol]],
            displayCol = displayCol
        )
        minusgeneNames <- defaultGenePriorities(
            data = minusgeneNames,
            assembly = genes$assembly,
            geneHighlights = genesInternal$geneHighlights[[displayCol]],
            displayCol = displayCol
        )

        if (!is.null(genesInternal$geneOrder)) {

            ## Integrate geneOrder and default prioritization
            plusgeneNames <- gene_priorities(
                genes = plusgeneNames,
                geneOrder = genesInternal$geneOrder,
                displayCol = displayCol
            )
            minusgeneNames <- gene_priorities(
                genes = minusgeneNames,
                geneOrder = genesInternal$geneOrder,
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
            fontsize = genesInternal$fontsize,
            xscale = genesInternal$xscale,
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
                        fontsize = genesInternal$fontsize
                    ),
                    vp = vp_gene,
                    default.units = "native",
                    check.overlap = TRUE
                )

                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = pgEnv),
                        child = plus_names
                    ),
                    envir = pgEnv
                )
            }
        }

        if (nrow(minusgeneNames) > 0) {
            checkedminusLabels <- apply(data.frame(
                "label" = minusgeneNames[[displayCol]],
                "labelLoc" = minusgeneNames$labelLoc
            ),
            1, cutoffLabel,
            fontsize = genesInternal$fontsize,
            xscale = genesInternal$xscale,
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
                        fontsize = genesInternal$fontsize
                    ),
                    vp = vp_gene,
                    default.units = "native",
                    check.overlap = TRUE
                )
                assign("gene_grobs",
                    addGrob(get("gene_grobs", envir = pgEnv),
                        child = minus_names
                    ),
                    envir = pgEnv
                )
            }
        }

        ##########################################################
        ## + AND - LABELS
        ##########################################################

        if (genesInternal$strandLabels == TRUE) {
            plus_label <- textGrob(
                label = "+", y = unit(0.63, "npc"),
                gp = gpar(
                    fontsize = genesInternal$fontsize + 2,
                    col = genesInternal$fontcolor[1],
                    fontface = "bold"
                ), vp = vp_label
            )
            minus_label <- textGrob(
                label = "-", y = unit(0.37, "npc"),
                gp = gpar(
                    fontsize = genesInternal$fontsize + 2,
                    col = genesInternal$fontcolor[2],
                    fontface = "bold"
                ), vp = vp_label
            )

            assign("gene_grobs",
                addGrob(get("gene_grobs", envir = pgEnv),
                    child = plus_label
                ),
                envir = pgEnv
            )
            assign("gene_grobs",
                addGrob(get("gene_grobs", envir = pgEnv),
                    child = minus_label
                ),
                envir = pgEnv
            )
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (genesInternal$draw == TRUE) {
        grid.draw(get("gene_grobs", envir = pgEnv))
    }


    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    genes$grobs <- get("gene_grobs", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("genes[", vp_gene$name, "]")
    invisible(genes)
}

## Define a function that removes gene and transcript
## name labels that will be cutoff
## @param df data.frame of genes in region
## @param fontsize fontsize of labels
## @param xscale vector of genomic region of viewport
## @param vp associated viewport where genes/transcripts are plotted
## @param unit for plotTranscripts, unit indicator
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
