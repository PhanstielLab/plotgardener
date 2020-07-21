




bb_plotPileup <- function(bed, chrom, chromstart, chromend, assembly = "hg19", fillcolor = "black", colorby = NULL, strandSplit = FALSE,
                          boxHeight =  unit(0.025, "inches"), spaceHeight = unit(.025, "inches"), x = NULL,
                          y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================




  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  pileup_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, width = width,
                                height = height, x = x, y = y, justification = just, grobs = NULL, assembly = assembly), class = "bb_pileup")
  attr(x = pileup_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = pileup_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  pileup_plot <- defaultUnits(object = pileup_plot, default.units = default.units)

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if (!"data.frame" %in% class(bed)){

    bed <- as.data.frame(data.table::fread(bed))

  }

  # ======================================================================================================================================================================================
  # SUBSET DATA FOR CHROMOSOME AND ANY OVERLAPPING REGIONS
  # ======================================================================================================================================================================================

   if (is.null(chromstart) & is.null(chromend)){

    if (assembly == "hg19"){
      genome <- bb_hg19
    }

    pileup_plot$chromstart <- 1
    pileup_plot$chromend <- genome[which(genome$chrom == chrom),]$length

  }

  bed <- bed[which(bed[,1] == pileup_plot$chrom & bed[,2] <= pileup_plot$chromend & bed[,3] >= pileup_plot$chromstart),]

  # ======================================================================================================================================================================================
  # COLORS
  # ======================================================================================================================================================================================

  if (is.null(colorby)){

    bed$color <- fillcolor

  } else {

    ## Find associated vector to colorby
    colorbyCol <- which(colnames(bed) == colorby)
    colorbyCol <- bed[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      colorbyCol <- as.numeric(colorbyCol)
    }

    if (class(fillcolor) == "function"){
      colorVec <- bb_maptocolors(colorbyCol, fillcolor)
    } else {
      colorbyCol <- factor(colorbyCol)
      mappedColors <- rep(fillcolor, ceiling(length(levels(colorbyCol))/length(fillcolor)))
      colorVec <- mappedColors[colorbyCol]
    }

    bed$color <- colorVec

  }

  # ======================================================================================================================================================================================
  # SEPARATE DATA INTO STRANDS
  # ======================================================================================================================================================================================

  if (strandSplit == TRUE){

    ## assuming strand is in the 6th column
    posStrand <- bed[which(bed[,6] == "+"),]
    minStrand <- bed[which(bed[,6] == "-"),]

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_pileup", length(grep(pattern = "bb_pileup", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(1, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = c(pileup_plot$chromstart, pileup_plot$chromend),
                   yscale = c(0, 1),
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      vp$name <- "bb_pileup1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = pileup_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(pileup_plot$chromstart, pileup_plot$chromend),
                   yscale = c(0, convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)),
                   just = just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("pileup_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # DETERMINE ROWS FOR EACH ELEMENT
  # ======================================================================================================================================================================================

  ## Determine how many rows are going to fit based on boxHeight and spaceHeight
  if (is.null(x) & is.null(y)){

    pushViewport(vp)
    boxHeight <- convertHeight(boxHeight, unitTo = "npc", valueOnly = T)
    spaceHeight <- convertHeight(spaceHeight, unitTo = "npc", valueOnly = T)
    upViewport()

  } else {

    boxHeight <- convertHeight(boxHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    spaceHeight <- convertHeight(spaceHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

  }

  maxRows <- floor((vp$yscale[2] + spaceHeight)/(boxHeight + spaceHeight))


  if (strandSplit == FALSE){

    bed$row <- 0

    ## Randomize order of data
    bed <- bed[sample(nrow(bed)),]


  } else {
    posStrand <- posStrand[sample(nrow(posStrand)),]
    posStrand$row <- 0
    minStrand <- minStrand[sample(nrow(minStrand)),]
    minStrand$row <- 0


  }



  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================


  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================



}
