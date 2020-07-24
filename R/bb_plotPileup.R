#' plots data stored in BED format in a pileup style
#'
#' @param bed genomic data in BED format
#' @param chrom chromosome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param assembly desired genome assembly
#' @param fillcolor single value, vector, or palette specifying colors of pileup elements
#' @param linecolor linecolor of line outlining pileup elements
#' @param colorby name of column in bed data to scale colors by
#' @param colorbyrange the range of values to apply a colorby palette scale to, if colorby values are numeric
#' @param strandSplit logical indicating whether plus and minus-stranded elements should be separated
#' @param boxHeight height of pileup element boxes, as a numeric value with default units or a unit value
#' @param spaceHeight height of spacing between pileup element boxes, as a fraction of boxHeight
#' @param spaceWidth width of spacing between pileup element boxes, as a fraction of the plot's genomic range
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param just A string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @return Function will plot a pileup style plot and return a bb_pileup object
#'
#' @export

bb_plotPileup <- function(bed, chrom, chromstart = NULL, chromend = NULL, assembly = "hg19", fillcolor = "black", linecolor = NA, colorby = NULL, colorbyrange = NULL, strandSplit = FALSE,
                          boxHeight =  unit(2, "mm"), spaceHeight = 0.3, spaceWidth = 0.02, x = NULL,
                          y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotpileup <- function(pileup_plot){


    ## Can't have only one NULL chromstart or chromend
    if ((is.null(pileup_plot$chromstart) & !is.null(pileup_plot$chromend)) | (is.null(pileup_plot$chromend) & !is.null(pileup_plot$chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.", call. = FALSE)

    }


    if (!is.null(pileup_plot$chromstart) & !is.null(pileup_plot$chromend)){

      ## chromend > chromstart
      if (pileup_plot$chromend < pileup_plot$chromstart){

        stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)


      }

    }





  }

  strand_scale <- function(strandSplit, height){

    if (strandSplit == TRUE){

      yscale <- c(-height/2, height/2)

    } else {
      yscale <- c(0, height)
    }

    return(yscale)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  pileup_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, color_palette = NULL, zrange = colorbyrange, width = width,
                                height = height, x = x, y = y, justification = just, grobs = NULL, assembly = assembly), class = "bb_pileup")
  attr(x = pileup_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = pileup_plot)
  errorcheck_bb_plotpileup(pileup_plot = pileup_plot)

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
  # SGET COLORBY DATA
  # ======================================================================================================================================================================================

  if (!is.null(colorby)){
    colorbyCol <- which(colnames(bed) == colorby)
    colorbyCol <- bed[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      bed$colorby <- as.numeric(colorbyCol)
    } else {

      if (is.null(colorbyrange)){
        colorbyrange <- c(min(bed$colorbyvalue), max(bed$colorbyvalue))
        bed_plot$zrange <- colorbyrange
      }

      bed$colorby <- colorbyCol
    }

  } else {
    bed$colorby <- NA
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

    yscale <- strand_scale(strandSplit = strandSplit, height = 0.5)

    vp <- viewport(height = unit(0.5, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = c(pileup_plot$chromstart, pileup_plot$chromend),
                   yscale = yscale,
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      vp$name <- "bb_pileup1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = pileup_plot)

    yscale <- strand_scale(strandSplit = strandSplit, height = convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE))

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(pileup_plot$chromstart, pileup_plot$chromend),
                   yscale = yscale,
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
    spaceHeight <- boxHeight*spaceHeight
    upViewport()

  } else {

    boxHeight <- convertHeight(boxHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    spaceHeight <- boxHeight*spaceHeight

  }

  maxRows <- floor((vp$yscale[2] + spaceHeight)/(boxHeight + spaceHeight))
  wiggle <- abs(pileup_plot$chromend - pileup_plot$chromstart) * spaceWidth


  if (strandSplit == FALSE){

    if (nrow(bed) > 0){

      bed$row <- 0

      ## Randomize order of data
      bed <- bed[sample(nrow(bed)),]

      ## Convert to numeric matrix for Rcpp function parsing
      bedMatrix <- as.matrix(bed[,c(2,3,ncol(bed)-1, ncol(bed))])

      ## Assign a row for each element
      rowDF <- checkRow(bedMatrix, maxRows, 3, wiggle)

      rowDF <- as.data.frame(rowDF)
      colnames(rowDF) <- c("start", "stop", "colorby", "row")


      if (any(rowDF$row == 0)){
        rowDF <- rowDF[which(rowDF$row != 0),]
        warning("Not enough plotting space for all provided pileup elements.", call. = FALSE)

        limitGrob <- textGrob(label = "+", x = unit(1, "npc"), y = unit(1, "npc"),
                              just = c("right", "top"), gp = gpar(col = "black"))
        assign("pileup_grobs", addGrob(gTree = get("pileup_grobs", envir = bbEnv), child = limitGrob), envir = bbEnv)

      }

      ## Change row index to 0
      rowDF$row <- rowDF$row - 1
      rowDF$width <- rowDF$stop - rowDF$start
      rowDF$y <- rowDF$row*(boxHeight + spaceHeight)

    } else {
      rowDF <- data.frame()
    }


  } else {

    if (nrow(posStrand) > 0){

      posStrand <- posStrand[sample(nrow(posStrand)),]
      posStrand$row <- 0
      ## Convert to numeric matrix for Rcpp function parsing
      posMatrix <- as.matrix(posStrand[,c(2,3,ncol(posStrand)-1, ncol(posStrand))])
      posDF <- checkRow(posMatrix, floor(maxRows/2), 3, wiggle)
      posDF <- as.data.frame(posDF)
      colnames(posDF) <- c("start", "stop", "colorby", "row")
      if (any(posDF$row == 0)){
        posDF <- posDF[which(posDF$row != 0),]
        warning("Not enough plotting space for all provided plus strand pileup elements.", call. = FALSE)
        limitGrob1 <- textGrob(label = "+", x = unit(1, "npc"), y = unit(1, "npc"),
                              just = c("right", "top"), gp = gpar(col = "black"))
        assign("pileup_grobs", addGrob(gTree = get("pileup_grobs", envir = bbEnv), child = limitGrob1), envir = bbEnv)
      }


      posDF$row <- posDF$row - 1
      posDF$width <- posDF$stop - posDF$start
      posDF$y <- (0.5*spaceHeight) + posDF$row*(boxHeight + spaceHeight)
      posDF$row <- posDF$row + 1
      posDF$row <- posDF$row + floor(maxRows/2)

    } else {
      posDF <- data.frame()
    }


    if (nrow(minStrand) > 0){

      minStrand <- minStrand[sample(nrow(minStrand)),]
      minStrand$row <- 0
      minMatrix <- as.matrix(minStrand[,c(2,3,ncol(minStrand)-1,ncol(minStrand))])
      minDF <- checkRow(minMatrix, floor(maxRows/2), 3, wiggle)
      minDF <- as.data.frame(minDF)
      colnames(minDF) <- c("start", "stop", "colorby", "row")
      if (any(minDF$row == 0)){
        minDF <- minDF[which(minDF$row != 0),]
        warning("Not enough plotting space for all provided minus strand pileup elements.", call. = FALSE)
        limitGrob2 <- textGrob(label = "+", x = unit(1, "npc"), y = unit(0, "npc"),
                               just = c("right", "bottom"), gp = gpar(col = "black"))
        assign("pileup_grobs", addGrob(gTree = get("pileup_grobs", envir = bbEnv), child = limitGrob2), envir = bbEnv)

      }

      minDF$row <- minDF$row - 1
      minDF$width <- minDF$stop - minDF$start
      minDF$y <- ((0.5*spaceHeight + boxHeight) + minDF$row*(boxHeight + spaceHeight))*-1
      rowIndex <- minDF$row + 1
      rowRange <- floor(maxRows/2):1
      minDF$row <- rowRange[rowIndex]

    } else {
      minDF <- data.frame()
    }

    rowDF <- rbind(posDF, minDF)
  }

  # ======================================================================================================================================================================================
  # COLORS
  # ======================================================================================================================================================================================

  if (is.null(colorby)){

    if (class(fillcolor) == "function"){
      colors <- fillcolor(maxRows)
      indeces <- rowDF$row
      rowDF$color <- colors[indeces]

    } else {

      if (length(fillcolor) == 1){
        rowDF$color <- fillcolor
      } else {

        colors <- rep(fillcolor, ceiling(maxRows/length(fillcolor)))[1:maxRows]
        indeces <- rowDF$row
        rowDF$color <- colors[indeces]

      }

    }

  } else {

    if (class(fillcolor) == "function"){

      rowDF$color <- bb_maptocolors(rowDF$colorby, fillcolor, range = colorbyrange)
      sorted_colors <- unique(rowDF[order(rowDF$colorby),]$color)

    } else {

      colorbyCol <- factor(rowDF$colorby)
      mappedColors <- rep(fillcolor, ceiling(length(levels(colorbyCol))/length(fillcolor)))
      rowDF$color <- mappedColors[rowDF$colorby]
    }


  }


  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (nrow(rowDF) > 0){

    bedRects <- rectGrob(x = rowDF$start,
                         y = rowDF$y,
                         width = rowDF$width,
                         height = boxHeight,
                         just = c("left", "bottom"),
                         default.units = "native",
                         gp = gpar(fill = rowDF$color, col = linecolor))
    assign("pileup_grobs", addGrob(gTree = get("pileup_grobs", envir = bbEnv), child = bedRects), envir = bbEnv)

  } else {

    warning("No pileup data to plot.", call. = FALSE)
  }

  if (strandSplit == TRUE){

    lineGrob <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = unit(0, "native"), y1 = unit(0, "native"),
                             gp = gpar(col = "black"))
    assign("pileup_grobs", addGrob(gTree = get("pileup_grobs", envir = bbEnv), child = lineGrob), envir = bbEnv)

  }


  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("pileup_grobs", envir = bbEnv))

  }
  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  pileup_plot$grobs <-  get("pileup_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(pileup_plot)
}
