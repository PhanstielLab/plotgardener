#' plots paired-end data for a single chromosome
#'
#' @param bedpe bed paired end data to be plotted
#' @param chrom chromsome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param fillcolor single value, vector, or palette specifying colors of bedpe elements
#' @param colorby name of column in bedpe data to scale fillcolors by
#' @param colorbyrange the range of values to apply a colorby palette scale to, if colorby values are numeric
#' @param linecolor border color
#' @param assembly desired genome assembly
#' @param boxHeight height of boxes at either end of bedpe element
#' @param spaceHeight height of space between boxes of different bedpe elements
#' @param spaceWidth width of spacing between bedpe elements, as a fraction of the plot's genomic range
#' @param x A numeric vector or unit object specifying x-location
#' @param y A numeric vector or unit object specifying y-location
#' @param width A numeric vector or unit object specifying width
#' @param height A numeric vector or unit object specifying height
#' @param just string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @export

bb_plotBedpe <- function(bedpe, chrom, chromstart = NULL, chromend = NULL, fillcolor = "black", colorby = NULL, colorbyrange = NULL, linecolor = NA, assembly = "hg19",
                         boxHeight = unit(2, "mm"), spaceHeight = 0.3, spaceWidth = 0.02, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a funciton that catches errors
  errorcheck_bb_plotBedpe <- function(bedpe, bedpe_plot, colorby){

    ## Can't have only one NULL chromstart or chromend
    if ((is.null(bedpe_plot$chromstart) & !is.null(bedpe_plot$chromend)) | (is.null(bedpe_plot$chromend) & !is.null(bedpe_plot$chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.", call. = FALSE)

    }




    if (!is.null(bedpe_plot$chromstart) & !is.null(bedpe_plot$chromend)){

      ## chromend > chromstart
      if (bedpe_plot$chromend < bedpe_plot$chromstart){

        stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)


      }

    }


    if (!is.null(colorby)){
      if (!any(colnames(bedpe) == colorby)){
        stop("Colorby column not found in data. Check colorby column name.", call. = FALSE)
      }

      if (length(which(colnames(bedpe) == colorby)) > 1){
        stop("Multiple matching colorby columns found in data. Please provide colorby column name with only one occurrence.", call. = FALSE)
      }
    }


  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bedpe_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, color_palette = NULL, zrange = colorbyrange,
                                 width = width, height = height, x = x, y = y, justification = just, grobs = NULL, assembly = assembly), class = "bb_bedpe")
  attr(x = bedpe_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CHECK PLACEMENT
  # ======================================================================================================================================================================================

  check_placement(object = bedpe_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  bedpe_plot <- defaultUnits(object = bedpe_plot, default.units = default.units)

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if (!"data.frame" %in% class(bedpe)){
    bedpe <- as.data.frame(data.table::fread(bedpe))
  } else {
    bedpe <- as.data.frame(bedpe)
  }

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotBedpe(bedpe = bedpe, bedpe_plot = bedpe_plot, colorby = colorby)

  # ======================================================================================================================================================================================
  # ORGANIZE DATA
  # ======================================================================================================================================================================================

  ## Get appropriate starts/stops
  start1 <- apply(bedpe[,c(2,3)], 1, min)
  stop1 <-  apply(bedpe[,c(2,3)], 1, max)
  start2 <- apply(bedpe[,c(5,6)], 1, min)
  stop2 <-  apply(bedpe[,c(5,6)], 1, max)
  bedpe[,2] <- start1
  bedpe[,3] <- stop1
  bedpe[,5] <- start2
  bedpe[,6] <- stop2


  # ======================================================================================================================================================================================
  # SUBSET DATA FOR CHROMOSOME AND ANY OVERLAPPING REGIONS
  # ======================================================================================================================================================================================
  if (is.null(chromstart) & is.null(chromend)){

    if (assembly == "hg19"){
      genome <- bb_hg19
    }

    bedpe_plot$chromstart <- 1
    bedpe_plot$chromend <- genome[which(genome$chrom == chrom),]$length

  }

  bedpe <- bedpe[which(bedpe[,1] == bedpe_plot$chrom & bedpe[,4] == bedpe_plot$chrom & bedpe[,2] <= bedpe_plot$chromend & bedpe[,6] >= bedpe_plot$chromstart),]


  # ======================================================================================================================================================================================
  # GET BOX WIDTHS AND TOTAL DISTANCES
  # ======================================================================================================================================================================================

  bedpe$width1 <- bedpe[,3] - bedpe[,2]
  bedpe$width2 <- bedpe[,6] - bedpe[,5]
  bedpe$pos1 <- rowMeans(bedpe[,2:3])
  bedpe$pos2 <- rowMeans(bedpe[,5:6])
  bedpe$distance <- abs(bedpe$pos2- bedpe$pos1)

  # ======================================================================================================================================================================================
  # SORT BY DISTANCE FOR PRETTIER PLOTTING
  # ======================================================================================================================================================================================

  bedpe <- bedpe[order(bedpe$distance, decreasing = TRUE),]

  # ======================================================================================================================================================================================
  # SET COLORBY DATA
  # ======================================================================================================================================================================================
  if (!is.null(colorby)){

    colorbyCol <- which(colnames(bedpe) == colorby)
    colorbyCol <- bedpe[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      bedpe$colorby <- as.numeric(colorbyCol)
    } else {

      bedpe$colorby <- colorbyCol
      if (is.null(colorbyrange)){
        colorbyrange <- c(min(bedpe$colorby), max(bedpe$colorby))
        bedpe_plot$zrange <- colorbyrange
      }

    }

  } else {

    bedpe$colorby <- NA

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_bedpe", length(grep(pattern = "bb_bedpe", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(0.5, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = c(bedpe_plot$chromstart, bedpe_plot$chromend),
                   yscale = c(0, 1),
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      vp$name <- "bb_bedpe1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = bedpe_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(bedpe_plot$chromstart, bedpe_plot$chromend),
                   yscale = c(0, convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)),
                   just = just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("bedpe_grobs", gTree(vp = vp), envir = bbEnv)


  # ======================================================================================================================================================================================
  # DETERMINE ROWS FOR EACH ELEMENT
  # ======================================================================================================================================================================================

  ## Determine how many bepe elements are going to fit based on boxHeight and space
  if (is.null(x) & is.null(y)){

    pushViewport(vp)
    boxHeight <- convertHeight(boxHeight, unitTo = "npc", valueOnly = T)
    spaceHeight <- boxHeight*spaceHeight
    upViewport()

  } else {

    boxHeight <- convertHeight(boxHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    spaceHeight <- boxHeight*spaceHeight

  }


  limit <- floor((vp$yscale[2] + spaceHeight)/(boxHeight + spaceHeight))
  wiggle <- abs(bedpe_plot$chromend - bedpe_plot$chromstart) * spaceWidth


  if (nrow(bedpe) > 0){

    bedpe$row <- 0

    ## Convert to numeric matrix for Rcpp function parsing
    bedpeMatrix <- as.matrix(bedpe[,c(2,6,5, (ncol(bedpe)-6):ncol(bedpe))])

    ## Assign a row for each element
    rowBedpe <- checkRow(bedpeMatrix, limit, 9, wiggle)

    rowBedpe <- as.data.frame(rowBedpe)
    colnames(rowBedpe) <- c("start1", "stop2", "start2", "width1", "width2", "pos1", "pos2", "distance", "colorby", "row")

    if (any(rowBedpe$row == 0)){
      rowBedpe <- rowBedpe[which(rowBedpe$row != 0),]
      warning("Not enough plotting space for all provided bedpe elements.", call. = FALSE)

      limitGrob <- textGrob(label = "+", x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
      assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = limitGrob), envir = bbEnv)

    }

    ## Change row index to 0
    rowBedpe$row <- rowBedpe$row - 1
    rowBedpe$y <- rowBedpe$row*(boxHeight + spaceHeight)

  } else {
    rowBedpe <- data.frame()
  }

  # ======================================================================================================================================================================================
  # COLORS
  # ======================================================================================================================================================================================

  if (is.null(colorby)){

    if (class(fillcolor) == "function"){
      colors <- fillcolor(limit)
      indeces <- rowBedpe$row + 1
      rowBedpe$color <- colors[indeces]

    } else {

      if (length(fillcolor) == 1){
        rowBedpe$color <- rep(fillcolor, nrow(rowBedpe))
      } else {

        colors <- rep(fillcolor, ceiling(limit/length(fillcolor)))[1:limit]
        indeces <- rowBedpe$row + 1
        rowBedpe$color <- colors[indeces]

      }

    }

  } else {

    if (class(fillcolor) == "function"){

      rowBedpe$color <- bb_maptocolors(rowBedpe$colorby, fillcolor, range = colorbyrange)
      sorted_colors <- unique(rowBedpe[order(rowBedpe$colorby),]$color)
      bedpe_plot$color_palette <- sorted_colors

    } else {

      colorbyCol <- factor(rowBedpe$colorby)
      mappedColors <- rep(fillcolor, ceiling(length(levels(colorbyCol))/length(fillcolor)))
      rowBedpe$color <- mappedColors[rowBedpe$colorby]
    }

  }


  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (nrow(rowBedpe) > 0){


    bedpeRect1 <- rectGrob(x = rowBedpe$start1,
                           y = rowBedpe$y,
                           width = rowBedpe$width1,
                           height = boxHeight,
                           just = c("left", "bottom"),
                           default.units = "native",
                           gp = gpar(fill = rowBedpe$color, col = linecolor))

    bedpeRect2 <- rectGrob(x = rowBedpe$start2,
                           y = rowBedpe$y,
                           width = rowBedpe$width2,
                           height = boxHeight,
                           just = c("left", "bottom"),
                           default.units = "native",
                           gp = gpar(fill = rowBedpe$color, col = linecolor))

    bedpeLine <- segmentsGrob(x0 = rowBedpe$pos1,
                              y0 = rowBedpe$y + 0.5*boxHeight,
                              x1 = rowBedpe$pos2,
                              y1 = rowBedpe$y + 0.5*boxHeight,
                              default.units = "native",
                              gp = gpar(col = rowBedpe$color))


    assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = bedpeRect1), envir = bbEnv)
    assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = bedpeRect2), envir = bbEnv)
    assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = bedpeLine), envir = bbEnv)

  } else {

    warning("Bedpe contains no values.", call. = FALSE)

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("bedpe_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  bedpe_plot$grobs <-  get("bedpe_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bedpe_plot)

}
