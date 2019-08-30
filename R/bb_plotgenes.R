#' @export
#'
bb_plotgenes <- function(gtf, chrom = "chr8", chromstart = 133600000, chromend = 134800000, color = "black",
                        width = 3, height = 1, x = 1, y = 1, units = "inches", fontsize = 12){

  ## Function to plot boxes where exons are
  draw_exon <- function(df, chromstart, chromend, color){

    start <- as.numeric(df[2])
    start.normalized <- normalize(start, chromstart, chromend)
    stop <- as.numeric(df[3])
    stop.normalized <- normalize(stop, chromstart, chromend)
    strand <- df[5]

    ## Separate plus and minus strand
    if (strand == "+"){
      ybottom <- 0.5
      ytop <- 0.7

    } else if (strand == "-"){
      ybottom <- 0
      ytop <- 0.2
    }

    grid.polygon(x = c(start.normalized, start.normalized, stop.normalized, stop.normalized), y = c(ybottom, ytop, ytop, ybottom),
                 gp = gpar(fill = color))
  }

  ## Function to plot lines where genes are
  draw_gene <- function(df, chromstart, chromend, color, fontsize){

    start <- as.numeric(df[2])
    start.normalized <- normalize(start, chromstart, chromend)
    stop <- as.numeric(df[3])
    stop.normalized <- normalize(stop, chromstart, chromend)
    strand <- df[5]
    gene_name <- df[14]
    text_x <- mean(c(start.normalized, stop.normalized))

    ## Separate plus and minus strand
    if (strand == "+"){
      y_coord <- 0.6
      text_y <- 0.85

    } else if (strand == "-"){
      y_coord <- 0.1
      text_y <- 0.35
    }
    grid.segments(x0 = start.normalized, y0 = y_coord, x1 = stop.normalized, y1 = y_coord, gp = gpar(col = color))
    grid.text(gene_name, x = text_x, y = text_y, gp = gpar(fontsize = fontsize))

  }

  gtf <- rtracklayer::import(gtf)
  gtf_df <- as.data.frame(gtf)

  ## Subset for region
  gtf_subset <- gtf_df[which(gtf_df$seqnames == chrom & gtf_df$start >= chromstart & gtf_df$end <= chromend ), ]

  ## Split into exons and genes
  gtf_exons <- gtf_subset[which(gtf_subset$type == "exon"), ]
  gtf_genes <- gtf_subset[which(gtf_subset$type == "gene"), ]

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  ## Convert x and y coordinates and height and width to same page_units
  old_x <- unit(x, units = units)
  old_y <- unit(y, units = units)
  old_height <- unit(height, units = units)
  old_width <- unit(width, units = units)
  new_x <- convertX(old_x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(old_y, unitTo = page_units, valueOnly = TRUE)
  new_height <- convertHeight(old_height, unitTo = page_units, valueOnly = TRUE)
  new_width <- convertWidth(old_width, unitTo = page_units, valueOnly = TRUE)

  ## Convert coordinates for viewport
  converted_coords = convert_coordinates(height = new_height, width = new_width, x = new_x, y = new_y, pageheight = page_height)
  vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                 y = unit(converted_coords[2], units = page_units))
  pushViewport(vp)

  ## Plot exons
  invisible(apply(gtf_exons, 1, draw_exon, chromstart = chromstart, chromend = chromend, col = color))

  ## Plot and label genes
  invisible(apply(gtf_genes, 1, draw_gene, chromstart = chromstart, chromend = chromend, col = color, fontsize = fontsize))

  ## Go back up viewport
  upViewport()

}
