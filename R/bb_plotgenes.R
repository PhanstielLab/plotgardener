#' plots a gene track for a specified region
#'
#' @param chrom chromsome of region to be plotted
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param assembly desired genome assembly
#' @param chromstart start position
#' @param chromend end position
#' @param fontcolors a vector indicating the font colors for the plus strand and minus strand gene labels
#' @param strandcolors a vector indicating the strand colors for the plus strand and minus strand
#' @param geneOrder an ordered vector of gene names to prioritize labeling
#' @param geneHighlights a two-column dataframe with gene names to highlight and their corresponding highlight color
#' @param geneBackground if geneHighlights is given, background color for genes that are not highlighted
#' @param stroke numerical value indicating the stroke width for gene body outlines
#' @param fontsize the size of gene label text (in points)
#' @param strandLabels A logical value indicating whether to include +/- strand labels
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param just string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @return Function will plot a gene track and return a bb_genes object
#'
#' @export
bb_plotGenes <- function(chrom, params = NULL, assembly = "hg19", chromstart = NULL, chromend = NULL, fontcolors = c("#2929ff", "#ff3434"),
                         strandcolors = c("#8a8aff", "#ff7e7e"), geneOrder = NULL, geneHighlights = NULL, geneBackground = "grey",
                         stroke = 0.1, fontsize = 8, strandLabels = T, x = NULL, y = NULL, width = NULL, height = unit(0.6, "inches"),
                         just = c("left", "top"), default.units = "inches", draw = T){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  errorcheck_genes <- function(chromstart, chromend){

    if (!is.null(chromstart) & is.null(chromend)){

      stop("If specifying \'chromstart\', need to provide \'chromend\'.", call. = FALSE)

    }

    if (!is.null(chromend) & is.null(chromstart)){

      stop("If specifying \'chromend\', need to provide \'chromstart\'.", call. = FALSE)

    }

    if (!is.null(chromstart) & !is.null(chromend)){

      ## chromstart cannot be larger than chromend

      if (chromstart > chromend){

        stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)
      }

    }

  }

  parse_starts <- function(range){
    ## Separate character range into two numeric coords
    range <- strsplit(range, "-")[[1]]
    start <- as.numeric(range[1])

    return(start)
  }

  parse_widths <- function(range){
    ## Separate character range into two numeric coords
    range <- strsplit(range, "-")[[1]]
    start <- as.numeric(range[1])
    end <- as.numeric(range[2])
    width <- end - start

    return(width)
  }

  exon_grobs <- function(df){

    exon_ranges <- as.list(strsplit(as.character(df[8]), ",")[[1]])

    if (length(exon_ranges) > 0){

      starts <- lapply(exon_ranges, parse_starts)
      widths <- lapply(exon_ranges, parse_widths)
      exons_dataframe <- cbind(unlist(starts), unlist(widths))

      if (df[4] == "+"){

        exons <- rectGrob(x = exons_dataframe[,1],
                          y = unit(0.63, "npc"),
                          just = "left",
                          width = exons_dataframe[,2],
                          height = unit(0.18, "npc"),
                          gp = gpar(fill = df[11],
                                    col = df[11],
                                    lwd = 1.25, alpha = 0.5),
                          vp = vp_gene,
                          default.units = "native")

      } else if (df[4] == "-"){

        exons <- rectGrob(x = exons_dataframe[,1],
                          y = unit(0.37, "npc"),
                          just = "left",
                          width = exons_dataframe[,2],
                          height = unit(0.18, "npc"),
                          gp = gpar(fill = df[11],
                                    col = df[11],
                                    lwd = 1.25, alpha = 0.5),
                          vp = vp_gene,
                          default.units = "native")

      }

      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = exons), envir = bbEnv)

    }

  }

  utr_grobs <- function(df){

    utr_ranges <- as.list(strsplit(as.character(df[9]), ",")[[1]])

    if (length(utr_ranges) > 0){

      starts <- lapply(utr_ranges, parse_starts)
      widths <- lapply(utr_ranges, parse_widths)
      utrs_dataframe <- cbind(unlist(starts), unlist(widths))

      if (df[4] == "+"){

        # invisible(lapply(utr_ranges, utr_grobs, yCoord = unit(0.63, "npc"), strandcolor = strandcolors[1]))
        utrs <- rectGrob(x = utrs_dataframe[,1],
                         y = unit(0.63, "npc"),
                         just = "left",
                         width = utrs_dataframe[,2],
                         height = unit(0.1, "npc"),
                         gp = gpar(fill = df[11], col = NA, alpha = 0.5),
                         vp = vp_gene,
                         default.units = "native")

      } else if (df[4] == "-"){

        #invisible(lapply(utr_ranges, utr_grobs, yCoord = unit(0.37, "npc"), strandcolor = strandcolors[2]))
        utrs <- rectGrob(x = utrs_dataframe[,1],
                         y = unit(0.37, "npc"),
                         just = "left",
                         width = utrs_dataframe[,2],
                         height = unit(0.1, "npc"),
                         gp = gpar(fill = df[11], col = NA, alpha = 0.5),
                         vp = vp_gene,
                         default.units = "native")

      }

      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = utrs), envir = bbEnv)
    }


  }

  cutoffLabel <- function(df, fontsize, xscale, vp){

    label <- df[1]
    location <- df[2]

    pushViewport(vp)
    labelWidth <- convertWidth(widthDetails(textGrob(label = label, gp = gpar(fontsize = fontsize))), unitTo = "native", valueOnly = T)
    upViewport()
    leftBound <- as.numeric(location) - 0.5*labelWidth
    rightBound <- as.numeric(location) + 0.5*labelWidth

    if (leftBound < xscale[1] | rightBound > xscale[2]){
      return(NA)
    } else {
      return(label)
    }

  }

  gene_priorities <- function(genes, geneOrder){

    ## Split list into the ones in geneOrder and the ones not
    subset <- genes[which(genes$Gene %in% geneOrder),]
    remaining <- genes[which(!genes$Gene %in% geneOrder),]

    ## Put the geneOrder subset into the same order
    subset <- subset[match(geneOrder, subset$Gene),]
    subset <- subset[which(!is.na(subset$Gene)),]

    ## Put the remaining genes in order of citation
    remaining <- remaining[order(remaining$Citations, decreasing = TRUE),]

    ## Recombine lists, putting geneOrder section at the top
    combined <- rbind(subset, remaining)

    return(combined)
  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(fontcolors)) fontcolors <- NULL
  if(missing(strandcolors)) strandcolors <- NULL
  if(missing(geneBackground)) geneBackground <- NULL
  if(missing(stroke)) stroke <- NULL
  if(missing(fontsize)) fontsize <- NULL
  if(missing(strandLabels)) strandLabels <- NULL
  if(missing(height)) height <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if chrom argument is missing (could be in object)
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_geneInternal <- structure(list(assembly = assembly, chrom = chrom, chromstart = chromstart, chromend = chromend, fontcolors = fontcolors,
                                    strandcolors = strandcolors, geneOrder = geneOrder, geneHighlights = geneHighlights, geneBackground = geneBackground,
                                    stroke = stroke, fontsize = fontsize, strandLabels = strandLabels, x = x, y = y, width = width, height = height,
                                    just = just, default.units = default.units, draw = draw), class = "bb_geneInternal")

  bb_geneInternal <- parseParams(bb_params = params, object_params = bb_geneInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_geneInternal$assembly)) bb_geneInternal$assembly <- "hg19"
  if(is.null(bb_geneInternal$fontcolors)) bb_geneInternal$fontcolors <- c("#2929ff", "#ff3434")
  if(is.null(bb_geneInternal$strandcolors)) bb_geneInternal$strandcolors <- c("#8a8aff", "#ff7e7e")
  if(is.null(bb_geneInternal$geneBackground)) bb_geneInternal$geneBackground <- "grey"
  if(is.null(bb_geneInternal$stroke)) bb_geneInternal$stroke <- 0.1
  if(is.null(bb_geneInternal$fontsize)) bb_geneInternal$fontsize <- 8
  if(is.null(bb_geneInternal$strandLabels)) bb_geneInternal$strandLabels <- TRUE
  if(is.null(bb_geneInternal$height)) bb_geneInternal$height <- unit(0.6, "inches")
  if(is.null(bb_geneInternal$just)) bb_geneInternal$just <- c("left", "top")
  if(is.null(bb_geneInternal$default.units)) bb_geneInternal$default.units <- "inches"
  if(is.null(bb_geneInternal$draw)) bb_geneInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  genes_plot <- structure(list(chrom = bb_geneInternal$chrom, chromstart = bb_geneInternal$chromstart, chromend = bb_geneInternal$chromend, width = bb_geneInternal$width,
                               height = bb_geneInternal$height, x = bb_geneInternal$x, y = bb_geneInternal$y, justification = bb_geneInternal$just, grobs = NULL,
                               assembly = bb_geneInternal$assembly), class = "bb_genes")
  attr(x = genes_plot, which = "plotted") <- bb_geneInternal$draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(genes_plot$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
  errorcheck_genes(chromstart = genes_plot$chromstart, chromend = genes_plot$chromend)
  check_placement(object = genes_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  genes_plot <- defaultUnits(object = genes_plot, default.units = bb_geneInternal$default.units)

  # ======================================================================================================================================================================================
  # GET APPROPRIATE BUILD DATA
  # ======================================================================================================================================================================================

  if (genes_plot$assembly == "hg19"){

    data <- bb_hg19gtf
    genome <- bb_hg19

  }

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  if (is.null(genes_plot$chromstart) & is.null(genes_plot$chromend)){

    ## Just chromosome
    data <- data[which(data$Chromosome == genes_plot$chrom),]
    genes_plot$chromstart <- 1
    genes_plot$chromend <- genome[which(genome$chrom == genes_plot$chrom),]$length

  } else {

    ## Chromosome and any overlapping regions of chromstart/chromend
    data <- data[which(data$Chromosome == genes_plot$chrom & data$Start <= genes_plot$chromend & data$Stop >= genes_plot$chromstart),]

  }

  ## Genes on plus strand and genes on minus strand
  plus_genes <- data[which(data$Strand == "+"),]
  minus_genes <- data[which(data$Strand == "-"),]


  ## Add strandColor and fontColors
  if (nrow(plus_genes) > 0){
    plus_genes$strandColor <- bb_geneInternal$strandcolors[1]
    plus_genes$fontColor <- bb_geneInternal$fontcolors[1]
  }

  if (nrow(minus_genes) > 0){
    minus_genes$strandColor <- bb_geneInternal$strandcolors[2]
    minus_genes$fontColor <- bb_geneInternal$fontcolors[2]
  }

  ## Get width of each gene
  plus_genes$width <- plus_genes$Stop - plus_genes$Start
  minus_genes$width <- minus_genes$Stop - minus_genes$Start

  ## Find genes to be highlighted
  if (!is.null(bb_geneInternal$geneHighlights)){

    colnames(bb_geneInternal$geneHighlights) = c("Gene", "strandColor")
    plus_genes$strandColor <- NULL
    minus_genes$strandColor <- NULL

    plusHighlight <- plus_genes[which(plus_genes$Gene %in% bb_geneInternal$geneHighlights[,1]),]
    minusHighlight <- minus_genes[which(minus_genes$Gene %in% bb_geneInternal$geneHighlights[,1]),]
    plusBackground <- plus_genes[which(!plus_genes$Gene %in% bb_geneInternal$geneHighlights[,1]),]
    minusBackground <- minus_genes[which(!minus_genes$Gene %in% bb_geneInternal$geneHighlights[,1]),]

    # Change highlight genes to highlight color
    plusHighlight <- merge(plusHighlight, bb_geneInternal$geneHighlights, by = "Gene")
    plusHighlight$fontColor <- plusHighlight$strandColor

    minusHighlight <- merge(minusHighlight, bb_geneInternal$geneHighlights, by = "Gene")
    minusHighlight$fontColor <- minusHighlight$strandColor

    # Change background genes to background color
    plusBackground$strandColor <- bb_geneInternal$geneBackground[1]
    plusBackground$fontColor <- bb_geneInternal$geneBackground[1]
    minusBackground$strandColor <- bb_geneInternal$geneBackground[1]
    minusBackground$fontColor <- bb_geneInternal$geneBackground[1]

    ## Put highlight and background genes back together
    plus_genes <- rbind(plusHighlight, plusBackground)
    minus_genes <- rbind(minusHighlight, minusBackground)
  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Define text grob for "+" and "-" viewport scaling
  tG <- textGrob(label = "+", gp = gpar(fontsize = bb_geneInternal$fontsize))

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_genes", length(grep(pattern = "bb_genes", x = currentViewports)) + 1)

  if (is.null(genes_plot$x) & is.null(genes_plot$y)){

    if (bb_geneInternal$strandLabels == TRUE){

      ## Make viewport for "+" and "-" labels to the left of the gene track
      vp_labelW <- convertWidth(widthDetails(tG) * 2, unitTo = "npc")

      vp_label <- viewport(height = unit(.12, "npc"),
                           width = vp_labelW,
                           x = unit(0, "npc"), y = unit(0.5, "npc"), just = "left",
                           name = paste0(vp_name, "_label"))

      vp_gene <- viewport(height = unit(.12, "npc"), width = unit(1, "npc") - vp_labelW,
                          x = vp_labelW, y = unit(0.5, "npc"),
                          clip = "on",
                          xscale = c(genes_plot$chromstart, genes_plot$chromend),
                          just = "left",
                          name = vp_name)
    } else {

      vp_gene <- viewport(height = unit(.12, "npc"), width = unit(1, "npc"),
                          x = unit(0, "npc"), y = unit(0.5, "npc"),
                          clip = "on",
                          xscale = c(genes_plot$chromstart, genes_plot$chromend),
                          just = "left",
                          name = vp_name)
    }


    if (bb_geneInternal$draw == TRUE){

      vp_gene$name <- "bb_genes1"
      if (bb_geneInternal$strandLabels == TRUE){
        vp_label$name <- "bb_genes1_label"
      }
      grid.newpage()
    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = genes_plot)

    ## Make viewport for gene track
    vp_gene <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(genes_plot$chromstart, genes_plot$chromend),
                   just = bb_geneInternal$just,
                   name = vp_name)

    if (bb_geneInternal$strandLabels == TRUE){

      ## Make viewport for "+" and "-" labels based on above viewport
      topLeft_vp <- vp_topLeft(viewport = vp_gene)

      vp_label <- viewport(height = page_coords$height,
                           width = convertWidth(widthDetails(tG) * 2, unitTo = get("page_units", envir = bbEnv)),
                           x = topLeft_vp[[1]], y = topLeft_vp[[2]], just = c("right", "top"),
                           name = paste0(vp_name, "_label"))
    }

  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("gene_grobs", gTree(), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  ##########################################################
  ## GENE LINES
  ##########################################################

  if ((genes_plot$chromend - genes_plot$chromstart) >= 25000000){

    if (nrow(plus_genes) > 0){

      plus_geneGrobs <- rectGrob(x = plus_genes$Start, y = unit(0.63, "npc"),
                                 width = plus_genes$width, height = unit(0.18, "npc"),
                                 just = "left", gp = gpar(fill = plus_genes$strandColor, col = makeTransparent(plus_genes$strandColor, alpha = 0.5), lwd = bb_geneInternal$stroke, alpha = 0.5),
                                 vp = vp_gene, default.units = "native")
      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = plus_geneGrobs), envir = bbEnv)


    }

    if (nrow(minus_genes) > 0 ){

      minus_geneGrobs <- rectGrob(x = minus_genes$Start, y = unit(0.37, "npc"),
                                  width = minus_genes$width, height = unit(0.18, "npc"),
                                  just = "left", gp = gpar(fill = minus_genes$strandColor, col = makeTransparent(minus_genes$strandColor, alpha = 0.5), lwd = bb_geneInternal$stroke, alpha = 0.5),
                                  vp = vp_gene, default.units = "native")
      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = minus_geneGrobs), envir = bbEnv)

    }


  } else {

    if (nrow(plus_genes) > 0){

      plus_geneGrobs <- rectGrob(x = plus_genes$Start, y = unit(0.63, "npc"),
                                 width = plus_genes$width, height = unit(0.05, "npc"),
                                 just = "left", gp = gpar(fill = plus_genes$strandColor, col = makeTransparent(plus_genes$strandColor, alpha = 0.5), lwd = bb_geneInternal$stroke, alpha = 0.5),
                                 vp = vp_gene, default.units = "native")
      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = plus_geneGrobs), envir = bbEnv)


    }

    if (nrow(minus_genes) > 0 ){

      minus_geneGrobs <- rectGrob(x = minus_genes$Start, y = unit(0.37, "npc"),
                                  width = minus_genes$width, height = unit(0.05, "npc"),
                                  just = "left", gp = gpar(fill = minus_genes$strandColor, col = makeTransparent(minus_genes$strandColor, alpha = 0.5), lwd = bb_geneInternal$stroke, alpha = 0.5),
                                  vp = vp_gene, default.units = "native")
      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = minus_geneGrobs), envir = bbEnv)

    }


    ##########################################################
    ## GENE EXONS AND UTRS
    ##########################################################

    if (nrow(plus_genes) > 0){
      invisible(apply(plus_genes, 1, exon_grobs))
      invisible(apply(plus_genes, 1, utr_grobs))
    }

    if (nrow(minus_genes) > 0){
      invisible(apply(minus_genes, 1, exon_grobs))
      invisible(apply(minus_genes, 1, utr_grobs))
    }

  }

  ##########################################################
  ## GENE NAME LABELS
  ##########################################################

  ## Add column with center location of each gene label
  plus_genes$label <- rowMeans(plus_genes[c("Start", "Stop")])
  minus_genes$label <- rowMeans(minus_genes[c("Start", "Stop")])

  ## Add all the gene names in the region to object
  geneNames <- c(plus_genes$Gene, minus_genes$Gene)
  genes_plot$genes <- geneNames

  ## Declutter labels
  if (is.null(bb_geneInternal$geneOrder)){

    ## No given gene order, just sort genes according to their citation number
    plus_genes <- plus_genes[order(plus_genes$Citations, decreasing = TRUE),]
    minus_genes <- minus_genes[order(minus_genes$Citations, decreasing = TRUE),]

  } else {

    ## Integrate geneOrder and citation number prioritization
    plus_genes <- gene_priorities(genes = plus_genes, geneOrder = bb_geneInternal$geneOrder)
    minus_genes <- gene_priorities(genes = minus_genes, geneOrder = bb_geneInternal$geneOrder)

  }

  ## Grobs
  if (nrow(plus_genes) > 0){

    checkedplusLabels <- apply(data.frame("label" = plus_genes$Gene, "labelLoc" = plus_genes$label), 1, cutoffLabel, fontsize = bb_geneInternal$fontsize,
                           xscale = c(genes_plot$chromstart, genes_plot$chromend),
                           vp = vp_gene)
    plus_genes$Gene <- checkedplusLabels
    plus_genes <- plus_genes[!is.na(plus_genes$Gene), ]

    if (nrow(plus_genes) > 0){

      plus_names <- textGrob(label = plus_genes$Gene,
                             x = plus_genes$label, y = unit(0.85, "npc"),
                             gp = gpar(col = plus_genes$fontColor, fontsize = bb_geneInternal$fontsize),
                             vp = vp_gene,
                             default.units = "native",
                             check.overlap = TRUE)

      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = plus_names), envir = bbEnv)

    }


  }

  if (nrow(minus_genes) > 0){

    checkedminusLabels <- apply(data.frame("label" = minus_genes$Gene, "labelLoc" = minus_genes$label), 1, cutoffLabel, fontsize = bb_geneInternal$fontsize,
                               xscale = c(genes_plot$chromstart, genes_plot$chromend),
                               vp = vp_gene)
    minus_genes$Gene <- checkedminusLabels
    minus_genes <- minus_genes[!is.na(minus_genes$Gene), ]

    if(nrow(minus_genes) > 0){

      minus_names <- textGrob(label = minus_genes$Gene,
                              x = minus_genes$label, y = unit(0.15, "npc"),
                              gp = gpar(col = minus_genes$fontColor, fontsize = bb_geneInternal$fontsize),
                              vp = vp_gene,
                              default.units = "native",
                              check.overlap = TRUE)
      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = minus_names), envir = bbEnv)

    }

  }

  ##########################################################
  ## + AND - LABELS
  ##########################################################

  if (bb_geneInternal$strandLabels == TRUE){

    plus_label <- textGrob(label = "+", y = unit(0.63, "npc"), gp = gpar(fontsize = bb_geneInternal$fontsize + 2, col = bb_geneInternal$fontcolors[1], fontface = "bold"), vp = vp_label)
    minus_label <- textGrob(label = "-", y = unit(0.37, "npc"), gp = gpar(fontsize = bb_geneInternal$fontsize + 2, col = bb_geneInternal$fontcolors[2], fontface = "bold"), vp = vp_label)

    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = plus_label), envir = bbEnv)
    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = minus_label), envir = bbEnv)

  }
  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_geneInternal$draw == TRUE){

    grid.draw(get("gene_grobs", envir = bbEnv))

  }


  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  genes_plot$grobs <- get("gene_grobs", envir = bbEnv)

  return(genes_plot)
}
