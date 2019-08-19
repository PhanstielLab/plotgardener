#' @export
#'
bb_plotgenes <- function(gtf, chrom = "chr8", chromstart = 133600000, chromend = 134800000, strand = "+", col = "black",
                         height = 1, width = 3, x, y, fill = NA, gene_height = 0.5){

  gtf <- rtracklayer::import(gtf)
  gtf_df <- as.data.frame(gtf)

  ## Subset for exons in region
  gtf_subset <- gtf_df[which(gtf_df$seqnames == chrom & gtf_df$start >= chromstart & gtf_df$end <= chromend & gtf_df$type == "exon" & gtf_df$strand == strand),]

  ## Get data in correct format
  gtf_subset <- c(gtf_subset$start, gtf_subset$stop, gtf_subset$gene_name)

  ## Viewport navigation
  if(is.null(current.vpPath()) == FALSE){
    upViewport()
  }

  ## Get page_height and margins from bbEnv
  page_height <- get("page_height", envir = bbEnv)

  ## Convert coordinated for viewport
  converted_coords = convert_coordinates(height = height, width = width, x = x, y = y, pageheight = page_height)

  vp <- viewport(height = unit(height, "in"), width = unit(width, "in"), x = unit(converted_coords[1], "in"), y = unit(converted_coords[2], "in"))
  pushViewport(vp)

  ## Plot long line through entire strand
  grid.segments(x0 = unit(0, "npc"), y0 = unit(0.5, "npc"), x1 = unit(1, "npc"), y1 = unit(0.5, "npc"))

  ## Function to plot boxes where exons are
  draw_exon <- function(df, chromstart, chromend, col, fill, gene_height){

    gene_name <- df[3]
    start <- as.numeric(df[1])
    start.normalized <- normalize(start, chromstart, chromend)
    stop <- as.numeric(df[2])
    stop.normalized <- normalize(stop, chromstart, chromend)
    ybottom <- 0.5 - (0.5 * gene_height)
    ytop <- 0.5 + (0.5 * gene_height)
    name_location <- mean(c(start.normalized, stop.normalized))

    grid.polygon(x = c(start.normalized, start.normalized, stop.normalized, stop.normalized), y = c(ybottom, ytop, ytop, ybottom),
                 gp = gpar(col = col, fill = fill))

    ## Gene label
    grid.text(gene_name, x = unit(name_location, "npc"), y = unit(0.5, "npc"))

  }

  invisible(apply(gtf_subset, 1, draw_exon, chromstart = chromstart, chromend = chromend, col = col, fill = fill, gene_height = gene_height))


}
