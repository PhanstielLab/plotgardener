#' @export
#'
bb_plotgenes <- function(gtf, chrom = "chr8", chromstart = 133600000, chromend = 134800000, fontcolors = c("#ff3434", "#2929ff"),
                         strandcolors = c("#ff7e7e", "#8a8aff"), width = 3, height = 1, x = 1, y = 1, just = c("left", "top"),
                         units = "inches", fontsize = 12, exclude = ""){

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
    grid.polygon(x = c(start, start, stop, stop), y = c(yBottom, yTop, yTop, yBottom),
                 gp = gpar(fill = fill, col = NA), default.units = "native")
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
      grid.text(geneName, x = xText, y = yText, just = c("center", "bottom"),
                gp = gpar(fontsize = fontsize, col = fontcolors[1]), default.units = "native")

    } else if (strand == "-"){

      yBottom <- 0.375
      yTop <- 0.425
      yText <- 0.15
      fill <- strandcolors[2]
      grid.text(geneName, x = xText, y = yText, just = c("center", "top"),
                gp = gpar(fontsize = fontsize, col = fontcolors[2]), default.units = "native")
    }

    grid.polygon(x = c(start, start, stop, stop),
                 y = c(yBottom, yTop, yTop, yBottom),
                 gp = gpar(fill = fill, col = NA), default.units = "native")

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  genes_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, width = width, height = height,
                               x = x, y = y, units = units, just = just), class = "genes_plot")

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


  ## Make viewport
  vp <- viewport(height = unit(page_coords[[1]]$height, page_coords[[3]]), width = unit(page_coords[[1]]$width, page_coords[[3]]),
                 x = unit(page_coords[[1]]$x, page_coords[[3]]), y = unit((page_coords[[2]]-page_coords[[1]]$y), page_coords[[3]]),
                 xscale = c(chromstart, chromend), just = just)
  pushViewport(vp)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  ## Plot exons
  invisible(apply(gtf_exons, 1, draw_exon, chromstart = chromstart, chromend = chromend, strandcolors = strandcolors))

  ## Plot and label genes
  invisible(apply(gtf_genes, 1, draw_gene, chromstart = chromstart, chromend = chromend, fontsize = fontsize,
                  fontcolors = fontcolors, strandcolors = strandcolors, exclude = exclude))

  ## Go back up viewport
  upViewport()

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(genes_plot)

}
