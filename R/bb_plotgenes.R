#' plots a track of genes
#' @param gtf gtf file or dataframe
#' @param chrom chromsome of gene region as a string
#' @param chromstart chromstart of gene region
#' @param chromend chromend of gene region
#' @param fontcolors a vector of 2 indicating the font colors for the plus and minus strands
#' @param strandcolors a vector of 2 indicating the strang color for the plus and minus strands
#' @param width width of plot
#' @param height height of plot
#' @param x x-coordinate of plot
#' @param y y-coordinate of plot
#' @param fontsize fontsize of gene names
#' @param exclude names of genes to exclude from gene labeling
#' @export
#'
bb_plotGenes <- function(gtf, chrom = "chr8", chromstart = 133600000, chromend = 134800000, fontcolors = c("#ff3434", "#2929ff"),
                         strandcolors = c("#ff7e7e", "#8a8aff"), width = unit(3, "inches"), height = unit(1, "inches"), x = unit(1, "inches"),
                         y = unit(1, "inches"), just = c("left", "top"),
                         fontsize = 12, exclude = ""){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================
  ## Define a function to draw boxes where exons are
  draw_exon <- function(df, chromstart, chromend, strandcolors){

    start <- as.numeric(df[2])
    stop <- as.numeric(df[3])
    strand <- df[4]

    ## Separate plus and minus strand
    if (strand == "+"){
      yBottom <- 0.5
      yTop <- 0.7
      fill <- strandcolors[1]

    } else if (strand == "-"){

      yBottom <- 0.5
      yTop <- 0.3
      fill <- strandcolors[2]

    }

    ## Draw exon boxes
    exon <- grid.polygon(x = c(start, start, stop, stop), y = c(yBottom, yTop, yTop, yBottom),
                 gp = gpar(fill = fill, col = NA), default.units = "native")
    assign("genes_grobs", addGrob(get("genes_grobs", envir = bbEnv), child = exon), envir = bbEnv)
  }

  ## Define a function to draw lines where genes are
  draw_gene <- function(df, chromstart, chromend, fontsize, fontcolors, strandcolors, exclude){

    start <- as.numeric(df[2])
    stop <- as.numeric(df[3])
    strand <- df[4]
    geneName <- df[6]

    ## Get center of gene
    xText <- mean(c(start, stop))

    ## Filter out unwanted gene names
    if(geneName %in% exclude){
      geneName <- ""
    }

    ## Define y positions for lines and text for + and - strands
    if (strand == "+"){

      yBottom <- 0.575
      yTop <- 0.625
      yText <- 0.85
      fill <- strandcolors[1]
      name <- grid.text(geneName, x = xText, y = yText, just = c("center", "bottom"),
                gp = gpar(fontsize = fontsize, col = fontcolors[1]), default.units = "native")

    } else if (strand == "-"){

      yBottom <- 0.375
      yTop <- 0.425
      yText <- 0.15
      fill <- strandcolors[2]
      name <- grid.text(geneName, x = xText, y = yText, just = c("center", "top"),
                gp = gpar(fontsize = fontsize, col = fontcolors[2]), default.units = "native")
    }

    gene <- grid.polygon(x = c(start, start, stop, stop),
                 y = c(yBottom, yTop, yTop, yBottom),
                 gp = gpar(fill = fill, col = NA), default.units = "native")

    assign("genes_grobs", setChildren(get("genes_grobs", envir = bbEnv), children = gList(name, gene)), envir = bbEnv)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  genes_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, width = width, height = height,
                               x = x, y = y, just = just, grobs = NULL), class = "genes_plot")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  #check_bbpage()

  # ======================================================================================================================================================================================
  # READ IN GTF
  # ======================================================================================================================================================================================

  if (class(gtf) %in% "data.frame"){
    gtf_df = gtf

  } else {
    gtf <- rtracklayer::import(gtf)
    gtf_df <- as.data.frame(gtf)
  }

  gtf_df <- data.frame(gtf_df$seqnames, gtf_df$start, gtf_df$end, gtf_df$strand, gtf_df$type, gtf_df$gene_name)
  colnames(gtf_df) <- c("seqnames", "start", "end", "strand", "type", "gene_name")

  # ======================================================================================================================================================================================
  # SUBSET GTF
  # ======================================================================================================================================================================================

  ## Subset for region
  gtf_subset <- gtf_df[which(gtf_df$seqnames == chrom & gtf_df$start >= chromstart & gtf_df$end <= chromend ), ]

  ## Split into exons and genes
  gtf_exons <- gtf_subset[which(gtf_subset$type == "exon"), ]
  gtf_genes <- gtf_subset[which(gtf_subset$type == "gene"), ]

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = genes_plot)

  ## Name viewport
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_genes", length(grep(pattern = "bb_genes", x = current_viewports)) + 1)

  ## Make viewport
  vp <- viewport(height = page_coords$height, width = page_coords$width,
                 x = page_coords$x, y = page_coords$y,
                 xscale = c(chromstart, chromend), just = just, name = vp_name)

  pushViewport(vp)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("genes_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  ## Plot exons
  if (nrow(gtf_exons) > 0){

    invisible(apply(gtf_exons, 1, draw_exon, chromstart = chromstart, chromend = chromend, strandcolors = strandcolors))

  }

  if(nrow(gtf_genes) > 0){

  ## Plot and label genes
  invisible(apply(gtf_genes, 1, draw_gene, chromstart = chromstart, chromend = chromend, fontsize = fontsize,
                  fontcolors = fontcolors, strandcolors = strandcolors, exclude = exclude))
  }



  ## Go back up viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  ## Add grobs to scale object
  genes_plot$grobs <- get("genes_grobs", envir = bbEnv)
  #grid.draw(genes_plot$grobs)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(genes_plot)

}
