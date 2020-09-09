#' plots a Manhattan plot
#'
#' @param bed bedfile for Manhattan plot, either .bed file or dataframe in bed format
#' @param pVals name of column in bedfile of corresponding p-values (will be converted to -log(10) space)
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param chrom chromosome of region to be plotted, if specific region desired
#' @param chromstart start of region to be plotted, if specific region desired
#' @param chromend end of region to be plotted, if specific region desired
#' @param assembly genome assembly, for entire genome plotting
#' @param colors single color, vector of colors, or color palette
#' @param space the space between each chromsome as a fraction of the width of the plot
#' @param cex number indiciating the amount by which points should be scaled relative to the default
#' @param transp transparency of colors
#' @param ymax fraction of max y value to set as height of plot
#' @param range y-range of p-values to plot (c(min, max))
#' @param sigVal numeric indicating the significance level
#' @param sigLine logical indicating whether to draw a line at the significance level
#' @param sigCol single color for coloring significant values
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param just A string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @export
bb_plotManhattan <- function(bed, pVals, params = NULL, chrom = NULL, chromstart = NULL, chromend = NULL, assembly = "hg19", fillcolor = "black", pch = 19, space = 0.01,
                             cex = 0.25, ymax = 1, range = NULL, sigVal = 5e-08, sigLine = F, sigCol = NULL, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"),
                             default.units = "inches", draw = T, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that checks for errors in bb_plotManhattan
  errorcheck_bb_plotmanhattan <- function(bedfile, chrom, assembly, chromstart, chromend, pVals, object, fillcolor){

    # errorcheck valid assembly
    if (is.null(chrom)){

      if (assembly != "hg19"){

        stop("Invalid genome assembly.", call. = FALSE)
      }

    }

    if (!is.null(chrom)){

      ## Need both chromstart and chromend if trying to do a region within a chrom

      if (!is.null(chromstart) & is.null(chromend)){

        stop("If specifying \'chromstart\', need to provide \'chromend\'.", call. = FALSE)

      }

      if (!is.null(chromend) & is.null(chromstart)){

        stop("If specifying \'chromend\', need to provide \'chromstart\'.", call. = FALSE)

      }

      if (!is.null(chromstart) & !is.null(chromend)){

        ## chromstart cannot be larger than chromend

        if (chromstart > chromend){

          stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)
        }

      }

    }

    ## pVals input
    if (class(pVals) == "character"){

      ## make sure can find name of pvalue column in the bedfile
      if (!pVals %in% colnames(bedfile)){

        stop("Name of p-value column not found in bedfile.", call. = FALSE)
      }

    } else {

      stop("P-value input must the name of a column in the bedfile.", call. = FALSE)

    }

    ## colors and corresponding types of plotted regions
    if (!is.null(chrom)){

      if (class(fillcolor) == "function"){

        stop("\'fillcolor\' cannot be a palette when plotting a specified chromsome range.", call. = FALSE)

      }

      if (length(fillcolor) > 1){

        stop("\'fillcolor\' cannot be a vector when plotting a specified chromsome range.", call. = FALSE)
      }

    }

    ## range
    if (!is.null(object$range)){

      ## zrange needs to be a vector
      if (!is.vector(object$range)){

        stop("\'range\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be length 2
      if (length(man_plot$range) != 2){

        stop("\'range\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(man_plot$range)){

        stop("\'range\' must be a vector of two numbers.", call. = FALSE)

      }

      ## second value should be larger than the first value
      if (man_plot$range[1] >= man_plot$range[2]){

        stop("\'range\' must be a vector of two numbers in which the 2nd value is larger than the 1st.", call. = FALSE)

      }

    }
  }

  ## Define a function that parses the data into an internal format
  parse_data <- function(bedfile, chrom, chromstart, chromend, pVals){

    ## Read in data if it's not a dataframe
    if (!"data.frame" %in% class(bedfile)){
      bedfile <- fread(bedfile)
    }

    ## Subset data
    if (!is.null(chrom)){

      if (is.null(chromstart) & is.null(chromend)){

        bedfile <- bedfile[which(bedfile[,1] == chrom),]

      } else {

        bedfile <- bedfile[which(bedfile[,1] == chrom & bedfile[,2] >= chromstart & bedfile[,2] <= chromend),]

      }

    }

    ## Parse pVals
    pVal_col <- which(colnames(bedfile) == pVals)
    pvals <- bedfile[,pVal_col]

    ## Put data into internal format
    bed_data <- data.frame("chr" = bedfile[,1], "pos" = bedfile[,2], "pval" = pvals)
    bed_data$chr <- as.character(bed_data$chr)

    return(bed_data)
  }

  ## Define a function that parses the genome assembly and gets its internal data
  internal_assembly <- function(assembly){

    if (assembly == "hg19"){

      return(bb_hg19)
    }

  }

  ## Define a function that adds offsets to the beddata based on the offsets of the genome assembly
  bed_offset <- function(bedData, offsetAssembly){

    offsetChrom <- function(offsetChrom, bedData){

      chromMatch <- as.character(offsetChrom[1])
      offset <- as.numeric(offsetChrom[3])
      bedMatches <- bedData[which(bedData[,1] == chromMatch),]

      bedMatches[,2] <- bedMatches[,2] + offset

      return(bedMatches)
    }

    updatedBed <- apply(offsetAssembly, 1, offsetChrom, bedData = bedData)
    updatedBed <- dplyr::bind_rows(updatedBed)

    return(updatedBed)

  }

  ## Define a function that adjusts the yrange of the plot
  manhattan_range <- function(bedData, object){

    if (is.null(object$range)){

      object$range <- c(0, object$ymax * max(-log10(bedData[,3])))

    }

    return(object)
  }

  ## Define a function that parses colors
  parse_color <- function(fillcolor, offsetAssembly, bedData){

    if (!is.null(offsetAssembly)){

      ## parse type of color input
      if (class(fillcolor) == "function"){

        newCol <- fillcolor(nrow(offsetAssembly))

      } else {

        newCol <- rep(fillcolor, ceiling(nrow(offsetAssembly)/length(fillcolor)))

      }

      ## Assign associated color number
      colNum <- length(newCol)
      colNum_vector <- rep_len(1:colNum, length.out = nrow(offsetAssembly))

      ## Get color based on number
      colVec <- newCol[colNum_vector]
      offsetAssembly <- cbind(offsetAssembly, colVec)

      ## Get associated chroms in bedfile and assign the color
      chromColor <- function(offsetChrom, bedData){

        chromMatch <- as.character(offsetChrom[1])
        chromCol <- as.character(offsetChrom[5])
        bedMatches <- bedData[which(bedData[,1] == chromMatch),]

        bedMatches <- cbind(bedMatches, rep(chromCol, nrow(bedMatches)))
        bedMatches[,4] <- as.character(bedMatches[,4])

        return(bedMatches)
      }


      colorBed <- apply(offsetAssembly, 1, chromColor, bedData = bedData)
      colorBed <- dplyr::bind_rows(colorBed)

    } else {

      colorBed <- cbind(bedData, rep(fillcolor, nrow(bedData)))
      colorBed[,4] <- as.character(colorBed[,4])

    }

    colnames(colorBed) <- c("chr", "pos", "pval", "color")

    return(colorBed)


  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(fillcolor)) fillcolor <- NULL
  if(missing(pch)) pch <- NULL
  if(missing(space)) space <- NULL
  if(missing(cex)) cex <- NULL
  if(missing(ymax)) ymax <- NULL
  if(missing(sigVal)) sigVal <- NULL
  if(missing(sigLine)) sigLine <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if bed/pVals arguments are missing (could be in object)
  if(!hasArg(bed)) bed <- NULL
  if(!hasArg(pVals)) pVals <- NULL

  ## Compile all parameters into an internal object
  bb_manInternal <- structure(list(bed = bed, pVals = pVals, chrom = chrom, chromstart = chromstart, chromend = chromend, assembly = assembly,
                                   fillcolor = fillcolor, pch = pch, space = space, cex = cex, ymax = ymax, range = range, sigVal = sigVal,
                                   sigLine = sigLine, sigCol = sigCol, x = x, y = y, width = width, height = height, just = just,
                                   default.units = default.units, draw = draw), class = "bb_manInternal")

  bb_manInternal <- parseParams(bb_params = params, object_params = bb_manInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_manInternal$assembly)) bb_manInternal$assembly <- "hg19"
  if(is.null(bb_manInternal$fillcolor)) bb_manInternal$fillcolor <- "black"
  if(is.null(bb_manInternal$pch)) bb_manInternal$pch <- 19
  if(is.null(bb_manInternal$space)) bb_manInternal$space <- 0.01
  if(is.null(bb_manInternal$cex)) bb_manInternal$cex <- 0.25
  if(is.null(bb_manInternal$ymax)) bb_manInternal$ymax <- 1
  if(is.null(bb_manInternal$sigVal)) bb_manInternal$sigVal <- 5e-08
  if(is.null(bb_manInternal$sigLine)) bb_manInternal$sigLine <- FALSE
  if(is.null(bb_manInternal$just)) bb_manInternal$just <- c("left", "top")
  if(is.null(bb_manInternal$default.units)) bb_manInternal$default.units <- "inches"
  if(is.null(bb_manInternal$draw)) bb_manInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  man_plot <- structure(list(range = bb_manInternal$range, ymax = bb_manInternal$ymax, space = bb_manInternal$space,
                             width = bb_manInternal$width, height = bb_manInternal$height, x = bb_manInternal$x, y = bb_manInternal$y,
                             justification = bb_manInternal$just, grobs = NULL, assembly = bb_manInternal$assembly), class = "bb_manhattan")
  attr(x = man_plot, which = "plotted") <- bb_manInternal$draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_manInternal$bed)) stop("argument \"bed\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_manInternal$pVals)) stop("argument \"pVals\" is missing, with no default.", call. = FALSE)

  check_placement(object = man_plot)
  errorcheck_bb_plotmanhattan(bedfile = bb_manInternal$bed, chrom = man_plot$chrom, assembly = man_plot$assembly, chromstart = man_plot$chromstart, chromend = man_plot$chromend,
                              pVals = bb_manInternal$pVals, object = man_plot,
                              fillcolor = bb_manInternal$fillcolor)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  man_plot <- defaultUnits(object = man_plot, default.units = bb_manInternal$default.units)

  # ======================================================================================================================================================================================
  # READ AND SUBSET DATA
  # ======================================================================================================================================================================================

  bed_data <- parse_data(bedfile = bb_manInternal$bed, chrom = man_plot$chrom, chromstart = man_plot$chromstart, chromend = man_plot$chromend, pVals = bb_manInternal$pVals)

  # ======================================================================================================================================================================================
  # WHOLE GENOME
  # ======================================================================================================================================================================================

  if (is.null(man_plot$chrom)){

    ## Add assembly to object to denote as a whole genome
    man_plot$chromstart <- NULL
    man_plot$chromend <- NULL

    ## Access internal assembly
    assembly_data <- internal_assembly(assembly = man_plot$assembly)

    ## get the offsets based on spacer for the assembly
    offsetAssembly <- parse_assembly(assemblyData = assembly_data, space = bb_manInternal$space)

    ## remove bed_data data that aren't in the genome assembly
    bed_data <- bed_data[bed_data[,1] %in% offsetAssembly[,1],]

    ## Add chromosome offsets to bed_data
    bed_data <- bed_offset(bedData = bed_data, offsetAssembly = offsetAssembly)

    ## Set viewport xscale
    cumsums <- cumsum(as.numeric(assembly_data[,2]))
    spacer <- cumsums[length(cumsum(as.numeric(assembly_data[,2])))] * bb_manInternal$space

    xscale <- c(0, max(offsetAssembly[,4]) + spacer)

  } else {

  # ======================================================================================================================================================================================
  # SINGLE CHROMOSOME
  # ======================================================================================================================================================================================
    man_plot$assembly <- NULL
    offsetAssembly <- NULL

    ## Whole single chromosome
    if (is.null(man_plot$chromstart) & is.null(man_plot$chromend)){

      man_plot$chromstart <- min(bed_data[,2])
      man_plot$chromend <- max(bed_data[,2])

      ## Set viewport xscale
      xscale <- c(min(bed_data[,2]), max(bed_data[,2]))


    } else {

      ## Set viewport xscale
      xscale <- c(man_plot$chromstart, man_plot$chromend)

    }

  }

  # ======================================================================================================================================================================================
  # Y-LIMITS
  # ======================================================================================================================================================================================

  man_plot <- manhattan_range(bedData = bed_data, object = man_plot)

  # ======================================================================================================================================================================================
  # SPLIT DATA INTO SIG AND NON-SIG VALUES FOR ADDITIONAL COLORING
  # ======================================================================================================================================================================================

  sigBed <- bed_data[bed_data$pval <= bb_manInternal$sigVal,]
  nonsigBed <- bed_data[bed_data$pval > bb_manInternal$sigVal,]

  # ======================================================================================================================================================================================
  # COLORS
  # ======================================================================================================================================================================================

  color_nonsig <- parse_color(fillcolor = bb_manInternal$fillcolor, offsetAssembly = offsetAssembly, bedData = nonsigBed)

  if (!is.null(bb_manInternal$sigCol)){

    color_sig <- parse_color(fillcolor = bb_manInternal$sigCol, offsetAssembly = offsetAssembly, bedData = sigBed)

  } else {

    color_sig <- parse_color(fillcolor = bb_manInternal$fillcolor, offsetAssembly = offsetAssembly, bedData = sigBed)
  }


  colorBed <- rbind(color_nonsig, color_sig)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_manhattan", length(grep(pattern = "bb_manhattan", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(man_plot$x) & is.null(man_plot$y)){

    vp <- viewport(height = unit(0.25, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = xscale, yscale = c(man_plot$range[1], man_plot$range[2]),
                   just = "center",
                   name = vp_name)

    if (bb_manInternal$draw == TRUE){

      vp$name <- "bb_manhattan1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = man_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = xscale, yscale = c(man_plot$range[1], man_plot$range[2]),
                   just = bb_manInternal$just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("manhattan_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================
  gp = gpar(...)
  if (length(gp) != 0){
    if ("col" %in% names(gp)){

      gp$linecolor <- gp$col

    }
  }
  gp$col <- colorBed$color
  gp$cex <- cex

  points <- pointsGrob(x = colorBed$pos, y = -log10(colorBed$pval), pch = bb_manInternal$pch,
                       gp = gp,
                       default.units = "native")
  assign("manhattan_grobs", addGrob(gTree = get("manhattan_grobs", envir = bbEnv), child = points), envir = bbEnv)

  # ======================================================================================================================================================================================
  # SIGLINE
  # ======================================================================================================================================================================================

  if (bb_manInternal$sigLine == TRUE){

    gp$col <- gp$linecolor
    sigGrob <- segmentsGrob(x0 = unit(0, "npc"), y0 = unit(-log10(bb_manInternal$sigVal), "native"), x1 = unit(1, "npc"), y1 = unit(-log10(bb_manInternal$sigVal), "native"),
                            gp = gp)
    assign("manhattan_grobs", addGrob(gTree = get("manhattan_grobs", envir = bbEnv), child = sigGrob), envir = bbEnv)
  }


  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_manInternal$draw == TRUE){

    grid.draw(get("manhattan_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  man_plot$grobs <- get("manhattan_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(man_plot)

}
