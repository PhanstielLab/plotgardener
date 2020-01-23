

#' @export
bb_plotGenesNEW <- function(assembly = "hg19", chrom, chromstart, chromend, fontcolors = c("#2929ff", "#ff3434"),
                            strandcolors = c("#8a8aff", "#ff7e7e"), width, height = unit(0.6, "inches"), x, y, just = c("left", "top"),
                            fontsize = 12, exclude = "", draw = T){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================


  exon_grobs <- function(range, yCoord, strandcolor){

    ## Separate character range into two numeric coords
    range <- strsplit(range, "-")[[1]]
    start <- as.numeric(range[1])
    end <- as.numeric(range[2])

    exonGrob <- rectGrob(x = unit(start, "native"), y = yCoord,
                         width = unit(end-start, "native"), height = unit(0.15, "npc"),
                         just = "left", gp = gpar(fill = strandcolor, col = NA))
    assign("gene_grobs", addGrob(get("gene_grobs", envir = bbEnv), child = exonGrob), envir = bbEnv)
  }


  plus_exon_grobs <- function(ranges, strandcolor){

    exon_ranges <- as.list(strsplit(as.character(ranges[9]), ",")[[1]])

    invisible(lapply(exon_ranges, exon_grobs, yCoord = unit(0.65, "npc"), strandcolor = strandcolor))

  }

  minus_exon_grobs <- function(ranges, strandcolor){

    exon_ranges <- as.list(strsplit(as.character(ranges[9]), ",")[[1]])
    invisible(lapply(exon_ranges, exon_grobs, yCoord = unit(0.35, "npc"), strandcolor = strandcolor))

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  genes_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, width = width, height = height,
                               x = x, y = y, justification = just, grobs = NULL), class = "bb_genes")
  attr(x = genes_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # GET APPROPRIATE BUILD DATA
  # ======================================================================================================================================================================================

  if (assembly == "hg19"){

    data <- bb_gene_data

  }


  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  data <- data[which(data$chromosome == chrom & data$start >= chromstart & data$stop <= chromend),]

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



  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("gene_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  plus_strand <- data[which(data$strand == "+"),]
  minus_strand <- data[which(data$strand == "-"),]

  #assign("plus_strand", plus_strand, envir = globalenv())
  #assign("minus_strand", minus_strand, envir = globalenv())

  invisible(apply(plus_strand,1, plus_exon_grobs, strandcolor = strandcolors[1]))
  invisible(apply(minus_strand, 1, minus_exon_grobs, strandcolor = strandcolors[2]))


  grid.draw(get("gene_grobs", envir = bbEnv))



}
