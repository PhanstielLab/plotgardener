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
                         boxHeight = unit(2, "mm"), spaceHeight = 0.3, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotBedpe <- function(bedpe, bedpe_plot){

    ## bedpe file
    ## if it's a file path, it needs to exist
    if (!"data.frame" %in% class(bedpe)){

      # ## File extension
      if (file_ext(bedpe) != "bedpe"){

        stop("Invalid input. File must have a \".bedpe\" extension")

      }

      ## File existence
      if (!file.exists(bedpe)){

        stop(paste("File", bedpe, "does not exist."), call. = FALSE)

      }

    }



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


  }


  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bedpe_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, color_palette = NULL, zrange = colorbyrange,
                                 width = width, height = height, x = x, y = y, justification = just, grobs = NULL, assembly = assembly), class = "bb_bedpe")
  attr(x = bedpe_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = bedpe_plot)
  errorcheck_bb_plotBedpe(bedpe = bedpe, bedpe_plot = bedpe_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  bedpe_plot <- defaultUnits(object = bedpe_plot, default.units = default.units)

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if (!"data.frame" %in% class(bedpe)){

    bedpe <- data.table::fread(bedpe)

  }

  # ======================================================================================================================================================================================
  # ORGANIZE DATA
  # ======================================================================================================================================================================================

  ## Assuming data is in first six columns only

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
  # COLORS
  # ======================================================================================================================================================================================
  if (is.null(colorby)){

    if (class(fillcolor) == "function"){
      bedpe$color <- fillcolor(nrow(bedpe))
    } else {

      if (length(fillcolor) == 1){
        bedpe$color <- fillcolor
      } else {
        bedpe$color <- rep(fillcolor, ceiling(nrow(bedpe)/length(fillcolor)))[1:nrow(bedpe)]
      }

    }

  } else {

    ## Find associated vector to colorby
    colorbyCol <- which(colnames(bedpe) == colorby)
    colorbyCol <- bedpe[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      bedpe$colorbyvalue <- as.numeric(colorbyCol)
    } else {

      if (is.null(colorbyrange)){
        colorbyrange <- c(min(bedpe$colorbyvalue), max(bedpe$colorbyvalue))
        bedpe_plot$zrange <- colorbyrange
      }
      bedpe$colorbyvalue <- colorbyCol
    }

    if (class(fillcolor) == "function"){

      bedpe$color <- bb_maptocolors(bedpe$colorbyvalue, fillcolor, range = colorbyrange)
      sorted_colors <- unique(bedpe[order(bedpe$colorbyvalue),]$color)
      bedpe_plot$color_palette <- sorted_colors

    } else {
      colorbyCol <- factor(colorbyCol)
      mappedColors <- rep(fillcolor, ceiling(length(levels(colorbyCol))/length(fillcolor)))
      bedpe$color <- mappedColors[colorbyCol]
    }

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
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (nrow(bedpe) > 0){

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

    if (nrow(bedpe) > limit){
      bedpe <- bedpe[1:limit,]
      warning("Not enough plotting space for all provided bedpe elements.", call. = FALSE)

      limitGrob <- textGrob(label = "+", x = unit(1, "npc"), y = unit(1, "npc"), just = c("right", "top"))
      assign("bedpe_grobs", addGrob(gTree = get("bedpe_grobs", envir = bbEnv), child = limitGrob), envir = bbEnv)


    } else if (nrow(bedpe) < limit){
      warning("Excess plotting space for the provided bedpe elements.", call. = FALSE)
    }

    bedpe$rowNum <- seq(0, nrow(bedpe) - 1, 1)
    bedpe$y <- bedpe$rowNum*(boxHeight + spaceHeight)

    bedpeRect1 <- rectGrob(x = bedpe$start1,
                           y = bedpe$y,
                           width = bedpe$width1,
                           height = boxHeight,
                           just = c("left", "bottom"),
                           default.units = "native",
                           gp = gpar(fill = bedpe$color, col = linecolor))

    bedpeRect2 <- rectGrob(x = bedpe$start2,
                           y = bedpe$y,
                           width = bedpe$width2,
                           height = boxHeight,
                           just = c("left", "bottom"),
                           default.units = "native",
                           gp = gpar(fill = bedpe$color, col = linecolor))

    bedpeLine <- segmentsGrob(x0 = bedpe$pos1,
                              y0 = bedpe$y + 0.5*boxHeight,
                              x1 = bedpe$pos2,
                              y1 = bedpe$y + 0.5*boxHeight,
                              default.units = "native",
                              gp = gpar(col = bedpe$color))


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
