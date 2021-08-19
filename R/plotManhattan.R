#' Plot a Manhattan plot
#' 
#' @usage plotManhattan(
#'     data,
#'     sigVal = 5e-08,
#'     chrom = NULL,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     assembly = "hg38",
#'     fill = "black",
#'     pch = 19,
#'     cex = 0.25,
#'     leadSNP = NULL,
#'     sigLine = FALSE,
#'     sigCol = NULL,
#'     ymax = 1,
#'     range = NULL,
#'     space = 0.01,
#'     bg = NA,
#'     baseline = FALSE,
#'     baseline.color = "grey",
#'     baseline.lwd = 1,
#'     x = NULL,
#'     y = NULL,
#'     width = NULL,
#'     height = NULL,
#'     just = c("left", "top"),
#'     flip = FALSE,
#'     default.units = "inches",
#'     draw = TRUE,
#'     params = NULL,
#'     ...
#' )
#'
#' @param data Data to be plotted, as a character value specifying a
#' file path of GWAS data, a dataframe, or a \link[GenomicRanges]{GRanges}
#' object. Each of these data types must have the following columns:
#' \itemize{
#' \item{\code{"chrom"}: }{Chromosome names. This column must be a character.}
#' \item{\code{"pos"}: }{Chromosomal position. This column must be
#' an integer or numeric.}
#' \item{\code{"p"}: }{p-value. This column must be numeric.
#' p-values will be converted to -log(10) space.}
#' \item{\code{"snp"}(optional): }{SNP name or rsid.
#' This column should be a character.}
#' }
#' @param sigVal A numeric specifying the significance level of p-values.
#' Along with data p-values, this value will be converted to -log(10) space.
#' Default value is \code{sigVal = 5e-08}.
#' @param chrom Chromosome of region to be plotted, as a string.
#' If left \code{NULL}, all chromosomes found in data will be plotted.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param fill A single character value, a vector, or a 
#' \link[plotgardener]{colorby} object specifying fill colors of data points.
#' For a Manhattan plot with multiple chromosomes, a vector of colors 
#' will be used to color points of different chromosomes.
#' Default value is \code{fill = "black"}.
#' @param pch A numeric value or numeric vector specifying point symbols.
#' If \link[plotgardener]{colorby} object is supplied for \code{fill},
#' point symbols will be mapped to
#' \code{colorby} values. Default value is \code{pch = 19}.
#' @param cex A numeric indicating the amount by which points should be
#' scaled relative to the default. Default value is \code{cex = 0.25}.
#' @param leadSNP A list specifying the lead SNP in the desired region and
#' any associated aesthetic features of the lead SNP data point and text label.
#' The lead SNP should be specified as a character with the name slot
#' \code{"snp"} in the list. Accepted lead SNP aesthetic
#' features in the list include
#' \code{fill}, \code{pch}, \code{cex}, \code{fontcolor}, and \code{fontsize}.
#' @param sigLine Logical value indicating whether to draw a line at the
#' significance level indicated with \code{sigVal}.
#' Default value is \code{sigLine = FALSE}.
#' @param sigCol Single character value specifying the color of
#' significant data points. If \code{scaleLD} is supplied,
#' \code{sigCol} will be ignored.
#' @param ymax A numeric specifying the fraction of the max y-value to
#' set as the height of the plot. Default value is \code{ymax = 1}.
#' @param range A numeric vector of length 2 specifying the y-range
#' of p-values to plot (c(min, max)).
#' @param space A numeric value indicating the space between each
#' chromosome as a fraction of the width of the plot, if plotting multiple
#' chromosomes. Default value is \code{space = 0.01}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a
#' baseline along the x-axis. Default value is \code{baseline = FALSE}.
#' @param baseline.color Baseline color. Default value
#' is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width. Default value
#' is \code{baseline.lwd = 1}.
#' @param x A numeric or unit object specifying Manhattan plot x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying Manhattan plot y-location.
#' The character value will
#' place the Manhattan plot y relative to the bottom of the most
#' recently plotted plot according to the units of the plotgardener page.
#' @param width A numeric or unit object specifying Manhattan plot width.
#' @param height A numeric or unit object specifying Manhattan plot height.
#' @param just Justification of Manhattan plot relative to its (x, y)
#' location. If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param flip Logical value indicating whether to reflect Manhattan plot
#' over the x-axis. Default value is \code{flip = FALSE}.
#' @param default.units A string indicating the default units to use
#' if \code{x}, \code{y}, \code{width}, or \code{height} are only given
#' as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should
#' be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[plotgardener]{pgParams} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{manhattan} object containing
#' relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load genomic assembly information
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' ## Load GWAS data
#' library(plotgardenerData)
#' data("hg19_insulin_GWAS")
#'
#' ## Create a page
#' pageCreate(width = 7.5, height = 4.5, default.units = "inches")
#'
#' ## Plot all GWAS data
#' manhattanPlot <- plotManhattan(
#'     data = hg19_insulin_GWAS, assembly = "hg19",
#'     fill = c("grey", "#37a7db"),
#'     sigLine = TRUE,
#'     col = "grey", lty = 2, range = c(0, 14),
#'     x = 0.5, y = 0, width = 6.5, height = 2,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#' ## Annotate genome label
#' annoGenomeLabel(
#'     plot = manhattanPlot, x = 0.5, y = 2, fontsize = 8,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#' plotText(
#'     label = "Chromosome", fontsize = 8,
#'     x = 3.75, y = 2.20, just = "center", default.units = "inches"
#' )
#'
#' ## Annotate y-axis
#' annoYaxis(
#'     plot = manhattanPlot, at = c(0, 2, 4, 6, 8, 10, 12, 14),
#'     axisLine = TRUE, fontsize = 8
#' )
#'
#' ## Plot y-axis label
#' plotText(
#'     label = "-log10(p-value)", x = 0.15, y = 1, rot = 90,
#'     fontsize = 8, fontface = "bold", just = "center",
#'     default.units = "inches"
#' )
#'
#'
#' ## Plot GWAS data zooming in on chromosome 11
#' ## highlighting a lead SNP, and coloring by LD score
#' hg19_insulin_GWAS$LD <- as.numeric(hg19_insulin_GWAS$LD)
#' ## Group LD column into LD ranges
#' hg19_insulin_GWAS <- as.data.frame(dplyr::group_by(hg19_insulin_GWAS,
#'                                         LDgrp = cut(
#'                                         hg19_insulin_GWAS$LD,
#'                                         c(0, 0.2, 0.4, 0.6, 0.8, 1))))
#' hg19_insulin_GWAS$LDgrp <- addNA(hg19_insulin_GWAS$LDgrp)
#' leadSNP_p <- min(hg19_insulin_GWAS[
#'     which(hg19_insulin_GWAS$chrom == "chr11"), ]$p)
#' leadSNP <- hg19_insulin_GWAS[which(hg19_insulin_GWAS$p == leadSNP_p), ]$snp
#' chr11_manhattanPlot <- plotManhattan(
#'     data = hg19_insulin_GWAS, chrom = "chr11",
#'     chromstart = 60000000,
#'     chromend = 130000000,
#'     assembly = "hg19",
#'     fill = colorby("LDgrp",
#'     palette = colorRampPalette(c(
#'         "#1f4297",
#'         "#37a7db", "green",
#'         "orange", "red", "grey"
#'     ))),
#'     sigLine = TRUE, col = "grey",
#'     lty = 2, range = c(0, 16),
#'     leadSNP = list(
#'         snp = leadSNP,
#'         pch = 18,
#'         cex = 0.75,
#'         fill = "#7ecdbb",
#'         fontsize = 8
#'     ),
#'     scaleLD = "LD",
#'     x = 0.5, y = 2.5, width = 6.5,
#'     height = 1.5,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )
#'
#' ## Plot legend for LD scores
#' plotLegend(
#'     legend = c(
#'         "LD Ref Var",
#'         paste("0.4", ">", "r^2", 
#'         "", ">=", "0.2"),
#'         paste("0.2", ">", "r^2", 
#'         "", ">=", "0"),
#'         "no LD data"
#'     ),
#'     fill = c("#7ecdbb", "#37a7db", "#1f4297", "grey"), cex = 0.75,
#'     pch = c(18, 19, 19, 19), border = FALSE, x = 7, y = 2.5,
#'     width = 1.5, height = 0.6, just = c("right", "top"),
#'     default.units = "inches"
#' )
#'
#'
#' ## Annotate genome label
#' annoGenomeLabel(
#'     plot = chr11_manhattanPlot, x = 0.5, y = 4.01,
#'     fontsize = 8, scale = "Mb",
#'     just = c("left", "top"), default.units = "inches"
#' )
#'
#' ## Annotate y-axis
#' annoYaxis(
#'     plot = chr11_manhattanPlot,
#'     at = c(0, 2, 4, 6, 8, 10, 12, 14, 16),
#'     axisLine = TRUE, fontsize = 8
#' )
#'
#' ## Plot y-axis label
#' plotText(
#'     label = "-log10(p-value)", x = 0.15, y = 3.25, rot = 90,
#'     fontsize = 8, fontface = "bold", just = "center",
#'     default.units = "inches"
#' )
#'
#' ## Hide page guides
#' pageGuideHide()
#' @details
#' A Manhattan plot can be placed on a plotgardener coordinate page by
#' providing plot placement parameters:
#' \preformatted{
#' plotManhattan(data,
#'                 chrom = NULL,
#'                 chromstart = NULL, chromend = NULL,
#'                 x, y, width, height, just = c("left", "top"),
#'                 default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated
#' Manhattan plot by ignoring plot placement parameters:
#' \preformatted{
#' plotManhattan(data,
#'                 chrom = NULL,
#'                 chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
plotManhattan <- function(data, sigVal = 5e-08, chrom = NULL,
                            chromstart = NULL, chromend = NULL,
                            assembly = "hg38", fill = "black", pch = 19,
                            cex = 0.25, leadSNP = NULL,
                            sigLine = FALSE, sigCol = NULL, ymax = 1,
                            range = NULL, space = 0.01, bg = NA,
                            baseline = FALSE, baseline.color = "grey",
                            baseline.lwd = 1, x = NULL, y = NULL,
                            width = NULL, height = NULL,
                            just = c("left", "top"),
                            flip = FALSE, default.units = "inches",
                            draw = TRUE, params = NULL, ...) {

    # =========================================================================
    # FUNCTIONS
    # =========================================================================

    ## Define a function that checks for errors in plotManhattan
    errorcheck_plotManhattan <- function(bedfile, chrom, chromstart,
                                            chromend, object,
                                            leadSNP, fill) {

        ## check bedfile columns
        if (!"chrom" %in% colnames(bedfile)) {
            stop("\'chrom\' column not found in data.", call. = FALSE)
        } else {
            if (!is(bedfile$chrom, "character")) {
                stop("\'chrom\' column must be a character.", call. = FALSE)
            }
        }
        
        if (!"pos" %in% colnames(bedfile)) {
            stop("\'pos\' column not found in data.", call. = FALSE)
        } else {
            if (!is(bedfile$pos, "numeric") &
                !is(bedfile$pos, "integer")) {
                stop("\'pos\' column must be an integer or numeric.",
                    call. = FALSE
                )
            }
        }
        if (!"p" %in% colnames(bedfile)) {
            stop("\'p\' column not found in data.", call. = FALSE)
        } else {
            if (!is(bedfile$p, "numeric")) {
                stop("\'p\' column must be numeric.", call. = FALSE)
            }
        }

        ## Genomic region
        if (!is.null(chrom)) {
            
            regionErrors(chromstart = chromstart,
                        chromend = chromend)
            
        } else {
            if (!is.null(chromstart) | !is.null(chromend)) {
                warning("Plotting multiple chromosomes. \'chromstart\' ",
                        "and \'chromend\' inputs will be ignored.",
                    call. = FALSE
                )
            }
        }

        ## range
        rangeErrors(range = object$range)

        ## lead SNP
        if (!is.null(leadSNP)) {
            if (!is(leadSNP, "list")) {
                stop("\'leadSNP\' must be a list with a \'snp\' name ",
                "slot and any other aesthetic options for that SNP, ",
                "like \'fill\', \'pch\', \'cex\', \'fontcolor\', ",
                "and \'fontsize\'.", call. = FALSE)
            }

            if (!"snp" %in% names(leadSNP)) {
                stop("\'leadSNP\' must be a list with a \'snp\' name ",
                "slot and any other aesthetic options for that SNP, ",
                "like \'fill\', \'pch\', \'cex\', \'fontcolor\', ",
                "and \'fontsize\'.", call. = FALSE)
            }
        }
        
        ## Colorby
        checkColorby(fill = fill,
                        colorby = TRUE,
                        data = bedfile)
    }

    ## Define a function that parses the data into an internal format
    parse_data <- function(bedfile, fill) {

        ## Rearrange/reformat columns in case not in order
        data <- data.frame(
            "chrom" = bedfile[, "chrom"],
            "pos" = bedfile[, "pos"],
            "p" = bedfile[, "p"]
        )

        ## Parse optional 'snp' column
        if ("snp" %in% colnames(bedfile)) {
            data$snp <- bedfile[, "snp"]
        }

        ## Get colorby column if applicable
        if (is(fill, "colorby")){
            colorbyCol <- bedfile[, fill$column]
            data[[fill$column]] <- colorbyCol
        }
        
        return(data)
    }

    ## Define a function that adds offsets to the beddata based on the
    ## offsets of the genome assembly
    bed_offset <- function(bedData, offsetAssembly) {
        offsetChrom <- function(offsetChrom, bedData) {
            chromMatch <- offsetChrom["chrom"]
            offset <- as.numeric(offsetChrom["start"])
            bedMatches <- bedData[which(bedData[, "chrom"] == chromMatch), ]
            bedMatches[, "pos"] <- bedMatches[, "pos"] + offset

            return(bedMatches)
        }

        updatedBed <- apply(offsetAssembly, 1, offsetChrom, bedData = bedData)
        updatedBed <- dplyr::bind_rows(updatedBed)

        return(updatedBed)
    }

    ## Define a function that adjusts the yrange of the plot
    manhattan_range <- function(bedData, object) {
        if (is.null(object$range)) {
            object$range <- c(0, object$ymax * max(-log10(bedData[, "p"])))
        }

        return(object)
    }
    
    ## Define a function that maps colors to a multi-chrom Manhattan plot
    genome_color <- function(fill, offsetAssembly, bedData){
        
        newCol <- rep(
            fill,
            ceiling(nrow(offsetAssembly) / length(fill))
        )
        
        ## Assign associated color number
        colNum <- length(newCol)
        colNum_vector <- rep_len(seq(1, colNum),
                            length.out = nrow(offsetAssembly)
        )
        
        ## Get color based on number
        colVec <- newCol[colNum_vector]
        offsetAssembly <- cbind(offsetAssembly, colVec)
        
        ## Get associated chroms in bedfile and assign the color
        chromColor <- function(offsetChrom, bedData) {
            chromMatch <- offsetChrom["chrom"]
            chromCol <- as.character(offsetChrom["colVec"])
            bedMatches <- bedData[which(bedData[, "chrom"] == chromMatch), ]
            bedCols <- colnames(bedMatches)
            bedMatches <- cbind(bedMatches, rep(chromCol, nrow(bedMatches)))
            colnames(bedMatches) <- c(bedCols, "color")
            return(bedMatches)
        }
        
        colorBed <- apply(offsetAssembly, 1, chromColor, bedData = bedData)
        colorBed <- dplyr::bind_rows(colorBed)
        return(colorBed)
        
    }
    
    ## Define a function that maps a vector of pch values to a colorby column
    mapPch <- function(bedData, fill, pch){
        ## Map pch to unique colorby column levels
        pchCol <- bedData[,which(colnames(bedData) == fill$column)]
        if (!is(pchCol, "factor")){
            pchCol <- as.factor(pchCol)
        }
        pch <- rep(pch,
                    ceiling(length(levels(pchCol)) / length(pch))
        )[seq(1, length(pchCol))]
        names(pch) <- levels(pchCol)
        bedData$pch <- pch[pchCol]
        
        return(bedData)
    }

    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================

    manInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "manInternal"
    )

    ## Set gp
    manInternal$gp <- gpar(cex = manInternal$cex)
    manInternal$gp <- setGP(
        gpList = manInternal$gp,
        params = manInternal, ...
    )
    manInternal$gp$fill <- NA
    
    ## Justification
    manInternal$just <- justConversion(just = manInternal$just)
    
    # =========================================================================
    # INITIALIZE OBJECT
    # =========================================================================

    man_plot <- structure(list(
        chrom = manInternal$chrom,
        chromstart = manInternal$chromstart,
        chromend = manInternal$chromend,
        assembly = NULL,
        range = manInternal$range,
        ymax = manInternal$ymax,
        space = manInternal$space,
        color_palette = NULL,
        zrange = NULL,
        x = manInternal$x, y = manInternal$y,
        width = manInternal$width,
        height = manInternal$height,
        just = manInternal$just, grobs = NULL
    ),
    class = "manhattan"
    )
    attr(x = man_plot, which = "plotted") <- manInternal$draw

    # =========================================================================
    # CATCH MISSING ARGUMENT AND PLACEMENT ERRORS
    # =========================================================================

    if (is.null(manInternal$data)) stop("argument \"data\" is missing, ",
                                        "with no default.", call. = FALSE)

    check_placement(object = man_plot)
    
    # ========================================================================
    # PARSE ASSEMBLY
    # ========================================================================
    
    man_plot$assembly <- parseAssembly(assembly = manInternal$assembly)

    # =========================================================================
    # READ IN DATA
    # =========================================================================

    bedfile <- read_rangeData(data = manInternal$data,
                            assembly = man_plot$assembly,
                            type = "GWAS")
    
    # =========================================================================
    # CATCH MORE ERRORS
    # =========================================================================

    errorcheck_plotManhattan(
        bedfile = bedfile, chrom = man_plot$chrom,
        chromstart = man_plot$chromstart,
        chromend = man_plot$chromend,
        object = man_plot,
        leadSNP = manInternal$leadSNP,
        fill = manInternal$fill
    )

    # =========================================================================
    # FORMAT DATA
    # =========================================================================
    
    bed_data <- parse_data(
        bedfile = bedfile,
        fill = manInternal$fill
    )
    
    # =========================================================================
    # PARSE UNITS
    # =========================================================================

    man_plot <- defaultUnits(
        object = man_plot,
        default.units = manInternal$default.units
    )
    
    # =====================================================================
    # GENOMIC SCALE
    # =====================================================================
    
    if (nrow(bed_data) > 0) {

        # =====================================================================
        # MULTIPLE CHROMOSOMES
        # =====================================================================

        if (is.null(man_plot$chrom)) {
            man_plot$chromstart <- NULL
            man_plot$chromend <- NULL
            chroms <- as.character(unique(bed_data$chr))

            ## Get chrom sizes based on assembly data

            if (is(man_plot$assembly$TxDb, "TxDb")) {
                txdbChecks <- TRUE
            } else {
                
                if (!requireNamespace(man_plot$assembly$TxDb, quietly = TRUE)){
                    txdbChecks <- FALSE
                    warning("`", man_plot$assembly$TxDb, "` not available. ",
                    "Please load to generate Manhattan plot.", call. = FALSE)
                } else {
                    txdbChecks <- TRUE
                }
                
            }

            if (txdbChecks == TRUE) {
                if (is(man_plot$assembly$TxDb, "TxDb")) {
                    tx_db <- man_plot$assembly$TxDb
                } else {
                    tx_db <- eval(parse(text = 
                                paste0(as.name(man_plot$assembly$TxDb),
                                "::",
                                as.name(man_plot$assembly$TxDb))))
                }

                assembly_data <- as.data.frame(setDT(as.data.frame(
                    GenomeInfoDb::seqlengths(tx_db)
                ),
                keep.rownames = TRUE
                ))
                colnames(assembly_data) <- c("chrom", "length")
                assembly_data <- assembly_data[which(
                    assembly_data[, "chrom"] %in% chroms
                ), ]
                man_plot$chrom <- assembly_data[, "chrom"]

                if (any(!chroms %in% assembly_data[, "chrom"])) {
                    non_txdb <- chroms[which(!chroms 
                                        %in% assembly_data[, "chrom"])]
                    warning("Chromosome(s)",
                        "'", non_txdb, "'",
                        collapse = ", ",
                        "not found in",
                        "`", man_plot$assembly$TxDb$packageName, "`",
                        "and will be ignored.",
                        call. = FALSE
                    )
                }

                ## get the offsets based on spacer for the assembly
                offsetAssembly <- spaceChroms(
                    assemblyData = assembly_data,
                    space = manInternal$space
                )

                ## remove bed_data data that aren't in the genome assembly
                bed_data <- bed_data[bed_data$chrom %in% 
                                                offsetAssembly[, "chrom"], ]
                ## Add chromosome offsets to bed_data
                bed_data <- bed_offset(
                    bedData = bed_data,
                    offsetAssembly = offsetAssembly
                )
                
                ## Set viewport xscale
                cumsums <- cumsum(as.numeric(assembly_data[, "length"]))
                spacer <- cumsums[length(cumsum(as.numeric(
                    assembly_data[, "length"]
                )))] * manInternal$space

                xscale <- c(0, max(offsetAssembly[, "end"]) + spacer)
                
                # =============================================================
                # COLORS
                # =============================================================
                
                ## Vector of colors behaves differently for multi-chrom 
                ## Manhattan plot
                if (!is(manInternal$fill, "colorby")){
                    
                    if (length(manInternal$fill) > 1){
                        bed_data <- genome_color(fill = manInternal$fill,
                                            offsetAssembly = offsetAssembly,
                                            bedData = bed_data)
                    }
                    
                } else {
                    multichromColors <- parseColors(data = bed_data,
                                                fill = manInternal$fill,
                                                object = man_plot)
                    bed_data$color <- multichromColors[[1]]
                    man_plot <- multichromColors[[2]]
                }
                
            } else {
                xscale <- c(0, 1)
            }
            
            manInternal$xscale <- xscale
            manInternal$txdbChecks <- txdbChecks
        } else {

            # ==================================================================
            # SINGLE CHROMOSOME
            # ==================================================================
            
            scaleChecks <- genomicScale(object = man_plot,
                                        objectInternal = manInternal,
                                        plotType = "Manhattan plot")
            man_plot <- scaleChecks[[1]]
            manInternal <- scaleChecks[[2]]
            
            # ==================================================================
            # COLORS
            # ==================================================================
            ## Apply fill/colorby for a single chromosome
            chromColors <- parseColors(data = bed_data,
                                        fill = manInternal$fill,
                                        object = man_plot,
                                        subset = "manhattan")
            bed_data$color <- chromColors[[1]]
            man_plot <- chromColors[[2]]
            
            # ==================================================================
            # SUBSET DATA
            # ==================================================================
            
            bed_data <- bed_data[which(bed_data$chrom == man_plot$chrom), ]
            if (!is.null(man_plot$chromstart) & !is.null(man_plot$chromend)) {
                bed_data <- bed_data[which(
                    bed_data$pos >= man_plot$chromstart &
                        bed_data$pos <= man_plot$chromend), ]
            }
        }
        
        if (manInternal$txdbChecks != FALSE) {
            
            # =================================================================
            # PCH
            # =================================================================
            
            if (is(manInternal$fill, "colorby")){
                ## Map pch to colorby column values
                bed_data <- mapPch(bedData = bed_data,
                                fill = manInternal$fill,
                                pch = manInternal$pch)
            } else {
                bed_data$pch <- rep(manInternal$pch[1], nrow(bed_data))
            }
            
            # =================================================================
            # SIGNIFICANCE COLORING
            # =================================================================
            
            ## Significance coloring with sigCol
            if (!is.null(manInternal$sigCol)){
                if (nrow(bed_data[bed_data$p <= manInternal$sigVal, ]) > 0){
                    bed_data[bed_data$p <= manInternal$sigVal, ]$color <- 
                        manInternal$sigCol[1]
                }
            }

            # =================================================================
            # Y-LIMITS
            # =================================================================

            man_plot <- manhattan_range(bedData = bed_data, object = man_plot)
            
        } else {
            manInternal$xscale <- c(0, 1)
            man_plot$range <- c(0, 1)
        }
    } else {
        manInternal$txdbChecks <- TRUE
        manInternal$xscale <- c(0, 1)
        man_plot$range <- c(0, 1)
        warning("No data found in region.", call. = FALSE)
    }


    # =========================================================================
    # VIEWPORTS
    # =========================================================================

    ## Get viewport name
    currentViewports <- current_viewports()
    vp_name <- paste0(
        "manhattan",
        length(grep(
            pattern = "manhattan",
            x = currentViewports
        )) + 1
    )

    ## y-scale
    yscale <- c(man_plot$range[1], man_plot$range[2])
    if (manInternal$flip == TRUE) {
        yscale <- rev(yscale)
    }

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from bb_makepage
    ## Not translating into page_coordinates
    if (is.null(man_plot$x) | is.null(man_plot[["y"]])) {
        vp <- viewport(
            height = unit(0.25, "snpc"), width = unit(1, "snpc"),
            x = unit(0.5, "npc"), y = unit(0.5, "npc"),
            clip = "on",
            xscale = manInternal$xscale, yscale = yscale,
            just = "center",
            name = vp_name
        )

        if (manInternal$draw == TRUE) {
            vp$name <- "manhattan1"
            grid.newpage()
        }
    } else {
        addViewport(vp_name)

        ## Convert coordinates into same units as page
        page_coords <- convert_page(object = man_plot)

        ## Make viewport
        vp <- viewport(
            height = page_coords$height, width = page_coords$width,
            x = page_coords$x, y = page_coords$y,
            clip = "on",
            xscale = manInternal$xscale, yscale = yscale,
            just = manInternal$just,
            name = vp_name
        )
    }

    # =========================================================================
    # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
    # =========================================================================
    backgroundGrob <- rectGrob(gp = gpar(
        fill = manInternal$bg,
        col = NA
    ), name = "background")
    assign("manhattan_grobs", gTree(vp = vp, children = gList(backgroundGrob)),
        envir = pgEnv
    )

    if (nrow(bed_data) > 0 & manInternal$txdbChecks == TRUE) {

        # =====================================================================
        # SUBSET DATA FOR LEAD SNP
        # =====================================================================
        leadSNP_row <- data.frame(matrix(nrow = 0, ncol = 1))
        if (!is.null(manInternal$leadSNP)) {
            if ("snp" %in% colnames(bed_data)) {
                # Find index SNP in data
                leadSNP_row <- bed_data[which(
                    bed_data$snp == manInternal$leadSNP$snp
                ), ]
                if (nrow(leadSNP_row) > 0) {

                    ## Remove from data to plot separately
                    bed_data <- suppressMessages(dplyr::anti_join(
                        bed_data, leadSNP_row
                    ))
                } else {
                    warning("Specified lead SNP not found in data.",
                        call. = FALSE
                    )
                }
            }
        }
        # =====================================================================
        # POINTS
        # =====================================================================
        
        if ("col" %in% names(manInternal$gp)) {
            manInternal$gp$linecolor <- manInternal$gp$col
        }
        
        ## Assign color column to gp
        manInternal$gp$col <- bed_data$color
        
        ## Remove any additional lty information for point plotting
        manInternal$gp$linetype <- manInternal$gp$lty
        manInternal$gp$lty <- NULL
        
        points <- pointsGrob(
            x = bed_data$pos, y = -log10(bed_data$p),
            pch = bed_data$pch,
            gp = manInternal$gp,
            default.units = "native"
        )
        assign("manhattan_grobs",
            addGrob(
                gTree = get("manhattan_grobs", envir = pgEnv),
                child = points
            ),
            envir = pgEnv
        )

        # =====================================================================
        # LEAD SNP
        # =====================================================================

        if (nrow(leadSNP_row) > 0) {
            manInternal$gp$col <- manInternal$leadSNP$fill
            manInternal$gp$cex <- manInternal$leadSNP$cex
            if (is.null(manInternal$leadSNP$pch)) {
                manInternal$leadSNP$pch <- manInternal$pch[1]
            }
            if (is.null(manInternal$leadSNP$cex)) {
                manInternal$gp$cex <- manInternal$cex
            }

            point <- pointsGrob(
                x = leadSNP_row$pos, y = -log10(leadSNP_row$p),
                pch = manInternal$leadSNP$pch,
                gp = manInternal$gp,
                default.units = "native"
            )
            snp <- textGrob(
                label = leadSNP_row$snp, x = leadSNP_row$pos,
                y = unit(-log10(leadSNP_row$p), "native") + unit(1.5, "mm"),
                just = "bottom",
                gp = gpar(
                    fontsize = manInternal$leadSNP$fontsize,
                    col = manInternal$leadSNP$fontcolor
                ),
                default.units = "native"
            )
            assign("manhattan_grobs",
                addGrob(
                    gTree = get("manhattan_grobs", envir = pgEnv),
                    child = point
                ),
                envir = pgEnv
            )
            assign("manhattan_grobs",
                addGrob(
                    gTree = get("manhattan_grobs", envir = pgEnv),
                    child = snp
                ),
                envir = pgEnv
            )
        }

        # =====================================================================
        # SIGLINE
        # =====================================================================

        if (manInternal$sigLine == TRUE) {
            manInternal$gp$col <- manInternal$gp$linecolor
            manInternal$gp$lty <- manInternal$gp$linetype
            sigGrob <- segmentsGrob(
                x0 = unit(0, "npc"),
                y0 = unit(
                    -log10(manInternal$sigVal),
                    "native"
                ),
                x1 = unit(1, "npc"),
                y1 = unit(
                    -log10(manInternal$sigVal),
                    "native"
                ),
                gp = manInternal$gp
            )
            assign("manhattan_grobs",
                addGrob(
                    gTree = get("manhattan_grobs", envir = pgEnv),
                    child = sigGrob
                ),
                envir = pgEnv
            )
        }

        if (manInternal$baseline == TRUE) {
            baselineGrob <- segmentsGrob(
                x0 = unit(0, "npc"),
                y0 = 0,
                x1 = unit(1, "npc"),
                y1 = 0,
                gp = gpar(
                    col = manInternal$baseline.color,
                    lwd = manInternal$baseline.lwd
                ),
                default.units = "native"
            )
            assign("manhattan_grobs",
                addGrob(
                    gTree = get("manhattan_grobs", envir = pgEnv),
                    child = baselineGrob
                ),
                envir = pgEnv
            )
        }
    }

    # =========================================================================
    # IF PLOT == TRUE, DRAW GROBS
    # =========================================================================

    if (manInternal$draw == TRUE) {
        grid.draw(get("manhattan_grobs", envir = pgEnv))
    }

    # =========================================================================
    # ADD GROBS TO OBJECT
    # =========================================================================

    man_plot$grobs <- get("manhattan_grobs", envir = pgEnv)

    # =========================================================================
    # RETURN OBJECT
    # =========================================================================

    message("manhattan[", vp$name, "]")
    invisible(man_plot)
}
