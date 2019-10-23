bb_convertCoordinates <- function (x, y, width, height, units, pageHeight = NULL, pageUnits = NULL)
{
  ## Get pageHeight and pageUnits from environment if NULL
  if(is.null(pageHeight)) pageHeight <- get("page_height", envir = bbEnv)
  if(is.null(pageUnits)) pageUnits <- get("page_units", envir = bbEnv)

  ## Convert inputs to unit objects
  x       <- unit(x, units = units)
  y       <- unit(y, units = units)
  width   <- unit(width, units = units)
  height  <- unit(height, units = units)

  ## Convert unit objects to pageUnits
  x       <- convertX(x, unitTo = pageUnits, valueOnly = TRUE)
  y       <- convertY(y, unitTo = pageUnits, valueOnly = TRUE)
  width   <- convertWidth(width, unitTo = pageUnits, valueOnly = TRUE)
  height  <- convertHeight(height, unitTo = pageUnits, valueOnly = TRUE)

  ## Adjust positions to top left corner
  ybottom <- pageHeight - y
  x      <- x + (0.5 * width)
  y      <- ybottom - (0.5 * height)

  ## Return new x and y
  return(list(x = x, y = y, width = width, height = height))
}

bb_plotGenes <- function (gtf, chrom = "chr8", chromstart = 133600000, chromend = 134800000,
                          color = "black", width = 3, height = 1, x = 1, y = 1, units = "inches",
                          fontsize = 12, exclude="")
{
  ## Function for drawing exon boxes
  draw_exon <- function(df, chromstart, chromend, color) {

    ## Extract and convert df columns to appropriate data type
    start     <- as.numeric(df[2])
    stop      <- as.numeric(df[3])
    strand    <- df[5]

    ## Normalize start and stop coordinates
    startNorm <- BentoBox::normalize(start, chromstart, chromend)
    stopNorm  <- BentoBox::normalize(stop, chromstart, chromend)

    ## Define y positions for boxes for + and - strands
    if (strand == "+") {
      yBottom <- 0.5
      yTop    <- 0.7
      col   <- "firebrick"
    }
    else if (strand == "-") {
      yBottom <- 0.5
      yTop    <- 0.3
      col  <- "steelblue"
    }

    ## Draw grid exon boxes
    grid.polygon(x = c(startNorm, startNorm, stopNorm, stopNorm),
                 y = c(yBottom, yTop, yTop, yBottom), gp = gpar(fill = col, col = NA))
  }

  ## Function for drawing gene lines
  draw_gene <- function(df, chromstart, chromend, color, fontsize, height) {

    ## Extract and convert df columns to appropriate data type
    start     <- as.numeric(df[2])
    stop      <- as.numeric(df[3])
    strand    <- df[5]
    geneName  <- df[14]

    ## Normalize start and stop coordinates
    startNorm <- BentoBox::normalize(start, chromstart, chromend)
    stopNorm  <- BentoBox::normalize(stop, chromstart, chromend)

    ## Get center of gene
    xText     <- mean(c(startNorm, stopNorm))

    print(geneName)
    ## Filter out unwanted gene names
    if(geneName %in% exclude){
      geneName <- ""
    }

    ## Define y positions for lines and text for + and - strands
    if (strand == "+") {
      yBottom   <- 0.575
      yTop      <- 0.625
      yText     <- 0.85
      col   <- "firebrick"
      grid.text(geneName, x = xText, y = yText, just = c("center", "bottom"), gp = gpar(fontsize = fontsize, col = col))
    }
    else if (strand == "-") {
      yBottom   <- 0.375
      yTop      <- 0.425
      yText     <- 0.3
      col  <- "steelblue"
      grid.text(geneName, x = xText, y = yText, just = c("center", "top"), gp = gpar(fontsize = fontsize, col = col))
    }

    grid.polygon(x = c(startNorm, startNorm, stopNorm, stopNorm),
                 y = c(yBottom, yTop, yTop, yBottom),
                 gp = gpar(fill = col, col = NA))


    # print((0.2)/convertHeight(x = stringHeight(string = "PBX"), unitTo = "npc", valueOnly = T))

  }

  ## Read in gtf as either data.frame or import from path
  if (class(gtf) %in% "data.frame") {
    gtf_df = gtf
  }
  else {
    gtf <- rtracklayer::import(gtf)
    gtf_df <- as.data.frame(gtf)
  }

  ## Subset gtf for region of interest; separate exons and genes
  gtf_subset  <- gtf_df[which(gtf_df$seqnames == chrom & gtf_df$start >= chromstart & gtf_df$end <= chromend), ]
  gtf_exons   <- gtf_subset[which(gtf_subset$type == "exon"), ]
  gtf_genes   <- gtf_subset[which(gtf_subset$type == "gene"), ]

  ## Get page height and unit parameters from the BentoBox environment
  page_height <- get("page_height", envir = bbEnv)
  page_units  <- get("page_units", envir = bbEnv)

  ## Convert coordinate system
  converted_coords <- bb_convertCoordinates(x = x, y = y, width = width, height = height, units = units,
                                            pageHeight = page_height, pageUnits = page_units)

  ## Define & push viewport with converted coordinates
  vp <- viewport(height = unit(converted_coords$height, page_units),
                 width = unit(converted_coords$width, page_units),
                 x = unit(converted_coords$x, units = page_units),
                 y = unit(converted_coords$y, units = page_units))
  pushViewport(vp)

  # grid.xaxis()
  # grid.yaxis()

  ## Draw exons and genes on subsetted gtf file
  invisible(apply(gtf_genes, 1, draw_gene, chromstart = chromstart,
                  chromend = chromend, col = color, fontsize = fontsize, height = height))
  invisible(apply(gtf_exons, 1, draw_exon, chromstart = chromstart,
                  chromend = chromend, col = color))

  ## Exit viewport
  upViewport()
}
