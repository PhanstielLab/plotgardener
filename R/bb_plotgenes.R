

#' @export
bb_plotGenes <- function(assembly = "hg19", chrom, chromstart, chromend, fontcolors = c("#2929ff", "#ff3434"),
                            strandcolors = c("#8a8aff", "#ff7e7e"), width = NULL, height = unit(0.6, "inches"), x = NULL, y = NULL,
                            just = c("left", "top"), fontsize = 8, exclude = "", strandLabels = T, draw = T){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

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

  exon_grobs <- function(df, strandcolors){

    exon_ranges <- as.list(strsplit(as.character(df[7]), ",")[[1]])

    if (length(exon_ranges) > 0){

      starts <- lapply(exon_ranges, parse_starts)
      widths <- lapply(exon_ranges, parse_widths)
      exons_dataframe <- cbind(unlist(starts), unlist(widths))


      if (df[3] == "+"){

        exons <- rectGrob(x = exons_dataframe[,1],
                          y = unit(0.63, "npc"),
                          just = "left",
                          width = exons_dataframe[,2],
                          height = unit(0.18, "npc"),
                          gp = gpar(fill = strandcolors[1],
                                    col = strandcolors[1],
                                    lwd = 1.25, alpha = 0.5),
                          vp = vp_gene,
                          default.units = "native")

      } else if (df[3] == "-"){

        exons <- rectGrob(x = exons_dataframe[,1],
                          y = unit(0.37, "npc"),
                          just = "left",
                          width = exons_dataframe[,2],
                          height = unit(0.18, "npc"),
                          gp = gpar(fill = strandcolors[2],
                                    col = strandcolors[2],
                                    lwd = 1.25, alpha = 0.5),
                          vp = vp_gene,
                          default.units = "native")

      }

      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = exons), envir = bbEnv)

    }

  }

  utr_grobs <- function(df, strandcolors){

    utr_ranges <- as.list(strsplit(as.character(df[8]), ",")[[1]])

    if (length(utr_ranges) > 0){

      starts <- lapply(utr_ranges, parse_starts)
      widths <- lapply(utr_ranges, parse_widths)
      utrs_dataframe <- cbind(unlist(starts), unlist(widths))

      if (df[3] == "+"){

        # invisible(lapply(utr_ranges, utr_grobs, yCoord = unit(0.63, "npc"), strandcolor = strandcolors[1]))
        utrs <- rectGrob(x = utrs_dataframe[,1],
                         y = unit(0.63, "npc"),
                         just = "left",
                         width = utrs_dataframe[,2],
                         height = unit(0.1, "npc"),
                         gp = gpar(fill = strandcolors[1], col = NA, alpha = 0.5),
                         vp = vp_gene,
                         default.units = "native")


      } else if (df[3] == "-"){

        #invisible(lapply(utr_ranges, utr_grobs, yCoord = unit(0.37, "npc"), strandcolor = strandcolors[2]))
        utrs <- rectGrob(x = utrs_dataframe[,1],
                         y = unit(0.37, "npc"),
                         just = "left",
                         width = utrs_dataframe[,2],
                         height = unit(0.1, "npc"),
                         gp = gpar(fill = strandcolors[2], col = NA, alpha = 0.5),
                         vp = vp_gene,
                         default.units = "native")

      }

      assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = utrs), envir = bbEnv)
    }


  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  genes_plot <- structure(list(chrom = gsub(pattern = "chr", replacement = "", x = chrom), chromstart = chromstart, chromend = chromend, width = width, height = height,
                               x = x, y = y, justification = just, grobs = NULL), class = "bb_genes")
  attr(x = genes_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = genes_plot)

  # ======================================================================================================================================================================================
  # GET APPROPRIATE BUILD DATA
  # ======================================================================================================================================================================================

  if (assembly == "hg19"){

    data <- bb_gene_data
    #data <- get("GENE_DATA", envir = globalenv())

  }

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  ## Chromosome and any overlapping regions
  data <- data[which(data$Chromosome == chrom & data$Start <= chromend & data$Stop >= chromstart),]

  ## Genes on plus strand and genes on minus strand
  plus_genes <- data[which(data$Strand == "+"),]
  minus_genes <- data[which(data$Strand == "-"),]

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Define text grob for "+" and "-" viewport scaling
  tG <- textGrob(label = "+", gp = gpar(fontsize = fontsize))

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_genes", length(grep(pattern = "bb_genes", x = currentViewports)) + 1)

  if (is.null(x) & is.null(y)){

    if (strandLabels == TRUE){

      ## Make viewport for "+" and "-" labels to the left of the gene track
      vp_labelW <- convertWidth(widthDetails(tG) * 2, unitTo = "npc")

      vp_label <- viewport(height = unit(.12, "npc"),
                           width = vp_labelW,
                           x = unit(0, "npc"), y = unit(0.5, "npc"), just = "left",
                           name = paste0(vp_name, "_label"))

      vp_gene <- viewport(height = unit(.12, "npc"), width = unit(1, "npc") - vp_labelW,
                          x = vp_labelW, y = unit(0.5, "npc"),
                          clip = "on",
                          xscale = c(chromstart, chromend),
                          just = "left",
                          name = vp_name)
    } else {

      vp_gene <- viewport(height = unit(.12, "npc"), width = unit(1, "npc"),
                          x = unit(0, "npc"), y = unit(0.5, "npc"),
                          clip = "on",
                          xscale = c(chromstart, chromend),
                          just = "left",
                          name = vp_name)


    }


    if (draw == TRUE){

      vp_gene$name <- "bb_genes1"
      if (strandLabels == TRUE){
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
                   xscale = c(chromstart, chromend),
                   just = just,
                   name = vp_name)

    if (strandLabels == TRUE){

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

  plus_genes$width <- plus_genes$Stop - plus_genes$Start

  minus_genes$width <- minus_genes$Stop - minus_genes$Start

  if (nrow(plus_genes) > 0){

    plus_geneGrobs <- rectGrob(x = plus_genes$Start, y = unit(0.63, "npc"),
                               width = plus_genes$width, height = unit(0.05, "npc"),
                               just = "left", gp = gpar(fill = strandcolors[1], col = NA, alpha = 0.5),
                               vp = vp_gene, default.units = "native")
    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = plus_geneGrobs), envir = bbEnv)

  }

  if (nrow(minus_genes) > 0 ){

    minus_geneGrobs <- rectGrob(x = minus_genes$Start, y = unit(0.37, "npc"),
                                width = minus_genes$width, height = unit(0.05, "npc"),
                                just = "left", gp = gpar(fill = strandcolors[2], col = NA, alpha = 0.5),
                                vp = vp_gene, default.units = "native")
    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = minus_geneGrobs), envir = bbEnv)

  }

  ##########################################################
  ## GENE EXONS
  ##########################################################

  invisible(apply(data, 1, exon_grobs, strandcolors = strandcolors))

  ##########################################################
  ## GENE UTRS
  ##########################################################

  invisible(apply(data, 1, utr_grobs, strandcolors = strandcolors))

  ##########################################################
  ## GENE NAME LABELS
  ##########################################################
  assign("plus_genes", plus_genes, envir = globalenv())
  assign("minus_genes", minus_genes, envir = globalenv())

  ## Add column with center location of each gene label
  plus_genes$label <- rowMeans(plus_genes[c("Start", "Stop")])
  minus_genes$label <- rowMeans(minus_genes[c("Start", "Stop")])


  ## Declutter labels
  ## Sort data according to their gene length so longer ones will be labeled first
  plus_genes <- plus_genes[order(plus_genes$width, decreasing = TRUE),]
  minus_genes <- minus_genes[order(minus_genes$width, decreasing = TRUE),]

  ## Additional prioritization can go here


  ## Grobs
  if (nrow(plus_genes) > 0){
    plus_names <- textGrob(label = plus_genes$Gene,
                           x = plus_genes$label, y = unit(0.85, "npc"),
                           gp = gpar(col = fontcolors[1], fontsize = fontsize),
                           vp = vp_gene,
                           default.units = "native",
                           check.overlap = TRUE)

    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = plus_names), envir = bbEnv)
  }

  if (nrow(minus_genes) > 0){

    minus_names <- textGrob(label = minus_genes$Gene,
                            x = minus_genes$label, y = unit(0.15, "npc"),
                            gp = gpar(col = fontcolors[2], fontsize = fontsize),
                            vp = vp_gene,
                            default.units = "native",
                            check.overlap = TRUE)
    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = minus_names), envir = bbEnv)
  }

  ##########################################################
  ## + AND - LABELS
  ##########################################################

  if (strandLabels == TRUE){

    plus_label <- textGrob(label = "+", y = unit(0.63, "npc"), gp = gpar(fontsize = fontsize + 2, col = fontcolors[1], fontface = "bold"), vp = vp_label)
    minus_label <- textGrob(label = "-", y = unit(0.37, "npc"), gp = gpar(fontsize = fontsize + 2, col = fontcolors[2], fontface = "bold"), vp = vp_label)

    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = plus_label), envir = bbEnv)
    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = minus_label), envir = bbEnv)

  }
  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("gene_grobs", envir = bbEnv))

  }


  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  genes_plot$grobs <- get("gene_grobs", envir = bbEnv)

  return(genes_plot)
}
