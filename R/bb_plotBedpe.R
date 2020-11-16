#' plots paired-end data for a single chromosome
#'
#' @param bedpe bed paired end data to be plotted
#' @param chrom chromsome of region to be plotted
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param chromstart start position
#' @param chromend end position
#' @param fillcolor single value, vector, or palette specifying colors of bedpe elements
#' @param colorby name of column in bedpe data to scale fillcolors by
#' @param colorbyrange the range of values to apply a colorby palette scale to, if colorby values are numeric
#' @param linecolor border color
#' @param assembly default genome assembly as a string or a bb_assembly object
#' @param boxHeight A numeric or unit object specifying height of boxes at either end of bedpe element
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

bb_plotBedpe <- function(bedpe, chrom, params = NULL, chromstart = NULL, chromend = NULL, fillcolor = "black", colorby = NULL, colorbyrange = NULL, linecolor = NA, assembly = "hg19",
                         boxHeight = unit(2, "mm"), spaceHeight = 0.3, spaceWidth = 0.02, x = NULL, y = NULL, width = NULL, height = NULL,
                         just = c("left", "top"), default.units = "inches", draw = TRUE, ...){

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
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(fillcolor)) fillcolor <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(boxHeight)) boxHeight <- NULL
  if(missing(spaceHeight)) spaceHeight <- NULL
  if(missing(spaceWidth)) spaceWidth <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if bedpe/chrom arguments are missing (could be in object)
  if(!hasArg(bedpe)) bedpe <- NULL
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_bedpeInternal <- structure(list(bedpe = bedpe, chrom = chrom, chromstart = chromstart, chromend = chromend,
                                     fillcolor = fillcolor, colorby = colorby, colorbyrange = colorbyrange, linecolor = linecolor,
                                     assembly = assembly, boxHeight = boxHeight, spaceHeight = spaceHeight, spaceWidth = spaceWidth,
                                     x = x, y = y, width = width, height = height, just = just, default.units = default.units, draw = draw), class = "bb_bedpeInternal")

  bb_bedpeInternal <- parseParams(bb_params = params, object_params = bb_bedpeInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_bedpeInternal$fillcolor)) bb_bedpeInternal$fillcolor <- "black"
  if(is.null(bb_bedpeInternal$linecolor)) bb_bedpeInternal$linecolor <- NA
  if(is.null(bb_bedpeInternal$assembly)) bb_bedpeInternal$assembly <- "hg19"
  if(is.null(bb_bedpeInternal$boxHeight)) bb_bedpeInternal$boxHeight <- unit(2, "mm")
  if(is.null(bb_bedpeInternal$spaceHeight)) bb_bedpeInternal$spaceHeight <- 0.3
  if(is.null(bb_bedpeInternal$spaceWidth)) bb_bedpeInternal$spaceWidth <- 0.02
  if(is.null(bb_bedpeInternal$just)) bb_bedpeInternal$just <- c("left", "top")
  if(is.null(bb_bedpeInternal$default.units)) bb_bedpeInternal$default.units <- "inches"
  if(is.null(bb_bedpeInternal$draw)) bb_bedpeInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bedpe_plot <- structure(list(chrom = bb_bedpeInternal$chrom, chromstart = bb_bedpeInternal$chromstart, chromend = bb_bedpeInternal$chromend, color_palette = NULL,
                               zrange = bb_bedpeInternal$colorbyrange, width = bb_bedpeInternal$width, height = bb_bedpeInternal$height, x = bb_bedpeInternal$x,
                               y = bb_bedpeInternal$y, justification = bb_bedpeInternal$just,
                               grobs = NULL, assembly = bb_bedpeInternal$assembly), class = "bb_bedpe")
  attr(x = bedpe_plot, which = "plotted") <- bb_bedpeInternal$draw

  # ======================================================================================================================================================================================
  # CHECK INPUTS/PLACEMENT
  # ======================================================================================================================================================================================

  if(is.null(bb_bedpeInternal$bedpe)) stop("argument \"bedpe\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_bedpeInternal$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)

  check_placement(object = bedpe_plot)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  bedpe_plot$assembly <- parse_bbAssembly(assembly = bedpe_plot$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  bedpe_plot <- defaultUnits(object = bedpe_plot, default.units = bb_bedpeInternal$default.units)
  if(!"unit" %in% class(bb_bedpeInternal$boxHeight)){

    if (!is.numeric(bb_bedpeInternal$boxHeight)){

      stop("\'boxHeight\' is neither a unit object or a numeric value. Cannot make bedpe plot.", call. = FALSE)

    }

    if (is.null(bb_bedpeInternal$default.units)){

      stop("\'boxHeight\' detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_bedpeInternal$boxHeight <- unit(bb_bedpeInternal$boxHeight, bb_bedpeInternal$default.units)
  }

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if (!"data.frame" %in% class(bb_bedpeInternal$bedpe)){
    bedpe <- as.data.frame(data.table::fread(bb_bedpeInternal$bedpe))
  } else {
    bedpe <- as.data.frame(bb_bedpeInternal$bedpe)
  }

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotBedpe(bedpe = bedpe, bedpe_plot = bedpe_plot, colorby = bb_bedpeInternal$colorby)

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
  # WHOLE CHROMOSOME DATA AND XSCALE
  # ======================================================================================================================================================================================

  if (is.null(bedpe_plot$chromstart) & is.null(bedpe_plot$chromend)){

    txdbChecks <- check_loadedPackage(package = bedpe_plot$assembly$TxDb, message = paste(paste0("`", bedpe_plot$assembly$TxDb,"`"),
                                                                                            "not loaded. Please install and load to plot full chromosome paired-end data."))
    xscale <- c(0, 1)
    if (txdbChecks == TRUE){
      tx_db <- eval(parse(text = bedpe_plot$assembly$TxDb))
      assembly_data <- seqlengths(tx_db)

      if (!bedpe_plot$chrom %in% names(assembly_data)){
        txdbChecks <- FALSE
        warning(paste("Chromosome", paste0("'", bedpe_plot$chrom, "'"), "not found in", paste0("`", bedpe_plot$assembly$TxDb, "`"), "and data for entire chromosome cannot be plotted."), call. = FALSE)
      } else {
        bedpe_plot$chromstart <- 1
        bedpe_plot$chromend <- assembly_data[[bedpe_plot$chrom]]
        xscale <- c(bedpe_plot$chromstart, bedpe_plot$chromend)

      }

    }

  } else {
    txdbChecks <- TRUE
    xscale <- c(bedpe_plot$chromstart, bedpe_plot$chromend)
  }

  # ======================================================================================================================================================================================
  # SUBSET DATA FOR CHROMOSOME AND ANY OVERLAPPING REGIONS
  # ======================================================================================================================================================================================

  if (!is.null(bedpe_plot$chromstart) & !is.null(bedpe_plot$chromend)){
    bedpe <- bedpe[which(bedpe[,1] == bedpe_plot$chrom & bedpe[,4] == bedpe_plot$chrom & bedpe[,2] <= bedpe_plot$chromend & bedpe[,6] >= bedpe_plot$chromstart),]
  } else {
    bedpe <- data.frame(matrix(nrow = 0, ncol = 6))
  }

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
  if (!is.null(bb_bedpeInternal$colorby) & nrow(bedpe) > 0){

    colorbyCol <- which(colnames(bedpe) == bb_bedpeInternal$colorby)
    colorbyCol <- bedpe[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      bedpe$colorby <- as.numeric(colorbyCol)
    } else {

      bedpe$colorby <- colorbyCol
    }

    if (is.null(bb_bedpeInternal$colorbyrange)){
      colorbyrange <- c(min(bedpe$colorby), max(bedpe$colorby))
      bedpe_plot$zrange <- colorbyrange
    }

  } else {

    bedpe$colorby <- rep(NA, nrow(bedpe))

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_bedpe", length(grep(pattern = "bb_bedpe", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(bedpe_plot$x) & is.null(bedpe_plot$y)){

    vp <- viewport(height = unit(0.5, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = xscale,
                   yscale = c(0, 1),
                   just = "center",
                   name = vp_name)

    if (bb_bedpeInternal$draw == TRUE){

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
                   xscale = xscale,
                   yscale = c(0, convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)),
                   just = bb_bedpeInternal$just,
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
  if (is.null(bedpe_plot$x) & is.null(bedpe_plot$y)){

    pushViewport(vp)
    boxHeight <- convertHeight(bb_bedpeInternal$boxHeight, unitTo = "npc", valueOnly = T)
    spaceHeight <- boxHeight*(bb_bedpeInternal$spaceHeight)
    upViewport()

  } else {

    boxHeight <- convertHeight(bb_bedpeInternal$boxHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    spaceHeight <- boxHeight*(bb_bedpeInternal$spaceHeight)

  }


  limit <- floor((as.numeric(vp$height) + spaceHeight)/(boxHeight + spaceHeight))
  wiggle <- abs(bedpe_plot$chromend - bedpe_plot$chromstart) * bb_bedpeInternal$spaceWidth


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

    # ======================================================================================================================================================================================
    # COLORS
    # ======================================================================================================================================================================================

    if (is.null(bb_bedpeInternal$colorby)){

      if (class(bb_bedpeInternal$fillcolor) == "function"){
        colors <- bb_bedpeInternal$fillcolor(limit)
        indeces <- rowBedpe$row + 1
        rowBedpe$color <- colors[indeces]

      } else {

        if (length(bb_bedpeInternal$fillcolor) == 1){
          rowBedpe$color <- rep(bb_bedpeInternal$fillcolor, nrow(rowBedpe))
        } else {

          colors <- rep(bb_bedpeInternal$fillcolor, ceiling(limit/length(bb_bedpeInternal$fillcolor)))[1:limit]
          indeces <- rowBedpe$row + 1
          rowBedpe$color <- colors[indeces]

        }

      }

    } else {

      if (class(bb_bedpeInternal$fillcolor) == "function"){

        rowBedpe$color <- bb_maptocolors(rowBedpe$colorby, bb_bedpeInternal$fillcolor, range = bedpe_plot$zrange)
        sorted_colors <- unique(rowBedpe[order(rowBedpe$colorby),]$color)
        bedpe_plot$color_palette <- sorted_colors

      } else {

        colorbyCol <- factor(rowBedpe$colorby)
        mappedColors <- rep(bb_bedpeInternal$fillcolor, ceiling(length(levels(colorbyCol))/length(bb_bedpeInternal$fillcolor)))
        rowBedpe$color <- mappedColors[rowBedpe$colorby]
      }

    }


    # ======================================================================================================================================================================================
    # MAKE GROBS
    # ======================================================================================================================================================================================
    bedpeRect1 <- rectGrob(x = rowBedpe$start1,
                           y = rowBedpe$y,
                           width = rowBedpe$width1,
                           height = boxHeight,
                           just = c("left", "bottom"),
                           default.units = "native",
                           gp = gpar(fill = rowBedpe$color, col = bb_bedpeInternal$linecolor, ...))

    bedpeRect2 <- rectGrob(x = rowBedpe$start2,
                           y = rowBedpe$y,
                           width = rowBedpe$width2,
                           height = boxHeight,
                           just = c("left", "bottom"),
                           default.units = "native",
                           gp = gpar(fill = rowBedpe$color, col = bb_bedpeInternal$linecolor, ...))

    bedpeLine <- segmentsGrob(x0 = rowBedpe$pos1,
                              y0 = rowBedpe$y + 0.5*boxHeight,
                              x1 = rowBedpe$pos2,
                              y1 = rowBedpe$y + 0.5*boxHeight,
                              default.units = "native",
                              gp = gpar(col = rowBedpe$color, ...))

    assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = bedpeRect1), envir = bbEnv)
    assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = bedpeRect2), envir = bbEnv)
    assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = bedpeLine), envir = bbEnv)


  } else {
    if (txdbChecks == TRUE){
      warning("Bedpe contains no values.", call. = FALSE)
    }

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_bedpeInternal$draw == TRUE){

    grid.draw(get("bedpe_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  bedpe_plot$grobs <- get("bedpe_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bedpe_plot)

}
