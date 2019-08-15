#' @export
#'
bb_plotgenes <- function(gtf, chrom = "chr8", chromstart = 133600000, chromend = 134800000, strand = "+", col = "black", height, width, x, y){

  gtf <- rtracklayer::import(gtf)
  gtf_df <- as.data.frame(gtf)

  ## Subset for exons in region
  gtf_subset <- gtf_df[which(gtf_df$seqnames == chrom & gtf_df$start >= chromstart & gtf_df$end <= chromend & gtf_df$type == "exon" & gtf_df$strand == strand),]

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

  ## Plot boxes where exons are
  draw_exon <- function(df, chromstart, chromend){
    start <- df$start
    start.normalized <- normalize(start, chromstart, chromend)
    stop <- df$stop
    stop.normalized <- normalize(stop, chromstart, chromend)
    ybottom <- 0.25
    ytop <- 0.75

    grid.polygon(x = c(start.normalized, start.normalized, stop.normalized, stop.normalized), y = c())


  }





}
