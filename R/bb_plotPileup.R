#' plots data stored in BED format in a pileup style
#'
#' @param bed genomic data in BED format
#' @param chrom chromosome of region to be plotted
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param chromstart start position
#' @param chromend end position
#' @param assembly desired genome assembly
#' @param fillcolor single value, vector, or palette specifying colors of pileup elements
#' @param linecolor linecolor of line outlining pileup elements
#' @param colorby name of column in bed data to scale colors by
#' @param colorbyrange the range of values to apply a colorby palette scale to, if colorby values are numeric
#' @param strandSplit logical indicating whether plus and minus-stranded elements should be separated
#' @param boxHeight A numeric or unit object specifying height of pileup element boxes
#' @param spaceHeight height of spacing between pileup element boxes, as a fraction of boxHeight
#' @param spaceWidth width of minimum spacing between pileup element boxes, as a fraction of the plot's genomic range
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

bb_plotPileup <- function(bed, chrom, params = NULL, chromstart = NULL, chromend = NULL, assembly = "hg19", fillcolor = "black", linecolor = NA,
                          colorby = NULL, colorbyrange = NULL, strandSplit = FALSE, boxHeight =  unit(2, "mm"), spaceHeight = 0.3,
                          spaceWidth = 0.02, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"),
                          default.units = "inches", draw = TRUE, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors
  errorcheck_bb_plotpileup <- function(bed, pileup_plot, colorby){

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

    if (!is.null(colorby)){
      if (!any(colnames(bed) == colorby)){
        stop("Colorby column not found in data. Check colorby column name.", call. = FALSE)
      }

      if (length(which(colnames(bed) == colorby)) > 1){
        stop("Multiple matching colorby columns found in data. Please provide colorby column name with only one occurrence.", call. = FALSE)
      }
    }




  }

  ## Define a function that parses the yscale based on split strands
  strand_scale <- function(strandSplit, height){

    if (strandSplit == TRUE){

      yscale <- c(-height/2, height/2)

    } else {
      yscale <- c(0, height)
    }

    return(yscale)
  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(fillcolor)) fillcolor <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(strandSplit)) strandSplit <- NULL
  if(missing(boxHeight)) boxHeight <- NULL
  if(missing(spaceHeight)) spaceHeight <- NULL
  if(missing(spaceWidth)) spaceWidth <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if bed/chrom arguments are missing (could be in object)
  if(!hasArg(bed)) bed <- NULL
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_pileInternal <- structure(list(bed = bed, chrom = chrom, chromstart = chromstart, chromend = chromend, assembly = assembly, fillcolor = fillcolor, linecolor = linecolor,
                                    colorby = colorby, colorbyrange = colorbyrange, strandSplit = strandSplit, boxHeight = boxHeight, spaceHeight = spaceHeight, spaceWidth = spaceWidth,
                                    x = x, y = y, width = width, height = height, just = just, default.units = default.units, draw = draw), class = "bb_pileInternal")

  bb_pileInternal <- parseParams(bb_params = params, object_params = bb_pileInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_pileInternal$assembly)) bb_pileInternal$assembly <- "hg19"
  if(is.null(bb_pileInternal$fillcolor)) bb_pileInternal$fillcolor <- "black"
  if(is.null(bb_pileInternal$linecolor)) bb_pileInternal$linecolor <- NA
  if(is.null(bb_pileInternal$strandSplit)) bb_pileInternal$strandSplit <- FALSE
  if(is.null(bb_pileInternal$boxHeight)) bb_pileInternal$boxHeight <- unit(2, "mm")
  if(is.null(bb_pileInternal$spaceHeight)) bb_pileInternal$spaceHeight <- 0.3
  if(is.null(bb_pileInternal$spaceWidth)) bb_pileInternal$spaceWidth <- 0.02
  if(is.null(bb_pileInternal$just)) bb_pileInternal$just <- c("left", "top")
  if(is.null(bb_pileInternal$default.units)) bb_pileInternal$default.units <- "inches"
  if(is.null(bb_pileInternal$draw)) bb_pileInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  pileup_plot <- structure(list(chrom = bb_pileInternal$chrom, chromstart = bb_pileInternal$chromstart, chromend = bb_pileInternal$chromend, color_palette = NULL,
                                zrange = bb_pileInternal$colorbyrange, width = bb_pileInternal$width, height = bb_pileInternal$height, x = bb_pileInternal$x, y = bb_pileInternal$y,
                                justification = bb_pileInternal$just, grobs = NULL, assembly = bb_pileInternal$assembly), class = "bb_pileup")
  attr(x = pileup_plot, which = "plotted") <- bb_pileInternal$draw

  # ======================================================================================================================================================================================
  # CHECK PLACEMENT/ARGUMENT ERROS
  # ======================================================================================================================================================================================

  if(is.null(bb_pileInternal$bed)) stop("argument \"bed\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_pileInternal$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
  check_placement(object = pileup_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  pileup_plot <- defaultUnits(object = pileup_plot, default.units = bb_pileInternal$default.units)
  if (!"unit" %in% class(bb_pileInternal$boxHeight)){

    if (!is.numeric(bb_pileInternal$boxHeight)){

      stop("\'boxHeight\' is neither a unit object or a numeric value. Cannot make pileup plot.", call. = FALSE)

    }

    if (is.null(bb_pileInternal$default.units)){

      stop("\'boxHeight\' detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_pileInternal$boxHeight <- unit(bb_pileInternal$boxHeight, bb_pileInternal$default.units)

  }

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if (!"data.frame" %in% class(bb_pileInternal$bed)){
    bed <- as.data.frame(data.table::fread(bb_pileInternal$bed))
  } else {
    bed <- as.data.frame(bb_pileInternal$bed)
  }

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotpileup(bed = bed, pileup_plot = pileup_plot, colorby = bb_pileInternal$colorby)

  # ======================================================================================================================================================================================
  # SUBSET DATA FOR CHROMOSOME AND ANY OVERLAPPING REGIONS
  # ======================================================================================================================================================================================

  ## EDIT HERE
   if (is.null(pileup_plot$chromstart) & is.null(pileup_plot$chromend)){

    if (bb_pileInternal$assembly == "hg19"){
      genome <- bb_hg19
    }

    pileup_plot$chromstart <- 1
    pileup_plot$chromend <- genome[which(genome$chrom == pieleup_plot$chrom),]$length

  }

  bed <- bed[which(bed[,1] == pileup_plot$chrom & bed[,2] <= pileup_plot$chromend & bed[,3] >= pileup_plot$chromstart),]

  # ======================================================================================================================================================================================
  # SET COLORBY DATA
  # ======================================================================================================================================================================================

  if (!is.null(bb_pileInternal$colorby)){
    colorbyCol <- which(colnames(bed) == bb_pileInternal$colorby)
    colorbyCol <- bed[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      bed$colorby <- as.numeric(colorbyCol)
    } else {

      bed$colorby <- colorbyCol
    }

    if (is.null(bb_pileInternal$colorbyrange)){
      colorbyrange <- c(min(bed$colorby), max(bed$colorby))
      pileup_plot$zrange <- colorbyrange
    }

  } else {
    bed$colorby <- rep(NA, nrow(bed))
  }

  # ======================================================================================================================================================================================
  # SEPARATE DATA INTO STRANDS
  # ======================================================================================================================================================================================

  if (bb_pileInternal$strandSplit == TRUE){

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
  if (is.null(pileup_plot$x) & is.null(pileup_plot$y)){

    yscale <- strand_scale(strandSplit = bb_pileInternal$strandSplit, height = 0.5)

    vp <- viewport(height = unit(0.5, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = c(pileup_plot$chromstart, pileup_plot$chromend),
                   yscale = yscale,
                   just = "center",
                   name = vp_name)

    if (bb_pileInternal$draw == TRUE){

      vp$name <- "bb_pileup1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = pileup_plot)

    yscale <- strand_scale(strandSplit = bb_pileInternal$strandSplit, height = convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE))

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(pileup_plot$chromstart, pileup_plot$chromend),
                   yscale = yscale,
                   just = bb_pileInternal$just,
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
  if (is.null(pileup_plot$x) & is.null(pileup_plot$y)){

    pushViewport(vp)
    boxHeight <- convertHeight(bb_pileInternal$boxHeight, unitTo = "npc", valueOnly = T)
    spaceHeight <- boxHeight*(bb_pileInternal$spaceHeight)
    upViewport()

  } else {

    boxHeight <- convertHeight(bb_pileInternal$boxHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    spaceHeight <- boxHeight*(bb_pileInternal$spaceHeight)

  }

  maxRows <- floor((vp$yscale[2] + spaceHeight)/(boxHeight + spaceHeight))
  wiggle <- abs(pileup_plot$chromend - pileup_plot$chromstart) * bb_pileInternal$spaceWidth


  if (bb_pileInternal$strandSplit == FALSE){

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

      ## Reset row for colors
      rowDF$row <- rowDF$row + 1

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

  if (is.null(bb_pileInternal$colorby)){

    if (class(bb_pileInternal$fillcolor) == "function"){
      colors <- bb_pileInternal$fillcolor(maxRows)
      indeces <- rowDF$row
      rowDF$color <- colors[indeces]

    } else {

      if (length(bb_pileInternal$fillcolor) == 1){
        rowDF$color <- rep(bb_pileInternal$fillcolor, nrow(rowDF))
      } else {

        colors <- rep(bb_pileInternal$fillcolor, ceiling(maxRows/length(bb_pileInternal$fillcolor)))[1:maxRows]
        indeces <- rowDF$row
        rowDF$color <- colors[indeces]

      }

    }

  } else {

    if (class(bb_pileInternal$fillcolor) == "function"){

      rowDF$color <- bb_maptocolors(rowDF$colorby, bb_pileInternal$fillcolor, range = pileup_plot$zrange)
      sorted_colors <- unique(rowDF[order(rowDF$colorby),]$color)
      pileup_plot$color_palette <- sorted_colors

    } else {

      colorbyCol <- factor(rowDF$colorby)
      mappedColors <- rep(bb_pileInternal$fillcolor, ceiling(length(levels(colorbyCol))/length(bb_pileInternal$fillcolor)))
      rowDF$color <- mappedColors[rowDF$colorby]
    }


  }


  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (nrow(rowDF) > 0){

    alpha <- 1
    if (!length(list(...)) == 0){
      if ("alpha" %in% names(list(...))){
        alpha <- list(...)$alpha
      }
    }

    bedRects <- rectGrob(x = rowDF$start,
                         y = rowDF$y,
                         width = rowDF$width,
                         height = boxHeight,
                         just = c("left", "bottom"),
                         default.units = "native",
                         gp = gpar(fill = rowDF$color, col = bb_pileInternal$linecolor, alpha = alpha))
    assign("pileup_grobs", addGrob(gTree = get("pileup_grobs", envir = bbEnv), child = bedRects), envir = bbEnv)

  } else {

    warning("No pileup data to plot.", call. = FALSE)
  }

  if (bb_pileInternal$strandSplit == TRUE){

    lineGrob <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = unit(0, "native"), y1 = unit(0, "native"),
                             gp = gpar(...))
    assign("pileup_grobs", addGrob(gTree = get("pileup_grobs", envir = bbEnv), child = lineGrob), envir = bbEnv)

  }


  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_pileInternal$draw == TRUE){

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
