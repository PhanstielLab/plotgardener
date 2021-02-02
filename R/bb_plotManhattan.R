#' Plot a Manhattan plot
#'
#' @param data Data to be plotted; as a character value specifying a BED file path, a dataframe in BED format, or a \link[GenomicRanges]{GRanges} object.
#' @param pVals Character value specifying the name of the \code{data} column of corresponding p-values (will be converted to -log(10) space).
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a \link[BentoBox]{bb_assembly} object. Default value is \code{assembly = "hg19"}.
#' @param sigLine Logical value indicating whether to draw a line at the significance level indicated with \code{sigVal}. Default value is \code{sigLine = FALSE}.
#' @param sigCol Single character value specifying the color of significant data points.
#' @param fill Character value(s) as a single value, vector, or palette specifying fill colors of data points. Default value is \code{fill = "black"}.
#' @param pch A numeric specifying point symbols. Default value is \code{pch = 19}.
#' @param cex A numeric indiciating the amount by which points should be scaled relative to the default. Default value is \code{cex = 0.25}.
#' @param ymax A numeric specifying the fraction of the max y-value to set as the height of the plot. Default value is \code{ymax = 1}.
#' @param range A numeric vector of length 2 specifying the y-range of p-values to plot (c(min, max)).
#' @param space A numeric value indicating the space between each chromsome as a fraction of the width of the plot, if plotting multiple chromosomes. Default value is \code{space = 0.01}.
#' @param bg Character value indicating background color. Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a baseline along the x-axis. Default value is \code{baseline = FALSE}.
#' @param x A numeric or unit object specifying Manhattan plot x-location.
#' @param y A numeric or unit object specifying Manhattan plot y-location.
#' @param width A numeric or unit object specifying Manhattan plot width.
#' @param height A numeric or unit object specifying Manhattan plot height.
#' @param just Justification of Manhattan plot relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{x}, \code{y}, \code{width}, or \code{height} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_manhattan} object containing relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load GWAS data
#' data("bb_gwasData")
#'
#' ## Plot Manhattan plot filling up entire graphic device
#' bb_plotManhattan(data = bb_gwasData, pVals = "pVal", chrom = "chr21",
#'                  chromstart = 28000000, chromend = 30300000, ymax = 1.1, cex = 0.20)
#'
#' ## Plot and place Manhattan plot on a BentoBox page
#' bb_pageCreate(width = 5, height = 2, default.units = "inches", xgrid = 0, ygrid = 0)
#' bb_plotManhattan(data = bb_gwasData, pVals = "pVal", chrom = "chr21",
#'                  chromstart = 28000000, chromend = 30300000, ymax = 1.1, cex = 0.20,
#'                  x = 0.5, y = 0.5, width = 4, height = 1.5,
#'                  just = c("left", "top"), default.units = "inches")
#'
#' @details
#' This function can be used to quickly plot a Manhattan plot by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotManhattan(data, pVals,
#'                  chrom = NULL,
#'                  chromstart = NULL, chromend = NULL)
#' }
#' A Manhattan plot can be placed on a BentoBox coordinate page by providing plot placement parameters:
#' \preformatted{
#' bb_plotManhattan(data, pVals,
#'                  chrom = NULL,
#'                  chromstart = NULL, chromend = NULL,
#'                  x, y, width, height, just = c("left", "top"),
#'                  default.units = "inches")
#' }
#'
#' @export
bb_plotManhattan <- function(data, pVals, sigVal = 5e-08, chrom = NULL, chromstart = NULL, chromend = NULL, assembly = "hg19", sigLine = FALSE, sigCol = NULL, fill = "black", pch = 19,
                             cex = 0.25, ymax = 1, range = NULL, space = 0.01, bg = NA, baseline = FALSE, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"),
                             default.units = "inches", draw = TRUE, params = NULL, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that checks for errors in bb_plotManhattan
  errorcheck_bb_plotmanhattan <- function(bedfile, chrom, chromstart, chromend, pVals, object, fillcolor){

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

    } else {

      if (!is.null(chromstart) | !is.null(chromend)){
        warning("Plotting multiple chromosomes. \'chromstart\' and \'chromend\' inputs will be ignored.", call. = FALSE)
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

        stop("\'fillcolor\' cannot be a palette when plotting a single, specified chromsome range.", call. = FALSE)

      }

      if (length(fillcolor) > 1){

        stop("\'fillcolor\' cannot be a vector when plotting a single, specified chromsome range.", call. = FALSE)
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

    ## Subset data
    if (!is.null(chrom)){
      bedfile <- bedfile[which(bedfile[,1] == chrom),]
      if (!is.null(chromstart) & !is.null(chromend)){
        bedfile <- bedfile[which(bedfile[,2] >= chromstart & bedfile[,2] <= chromend),]
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
  if(missing(fill)) fill <- NULL
  if(missing(pch)) pch <- NULL
  if(missing(space)) space <- NULL
  if(missing(cex)) cex <- NULL
  if(missing(ymax)) ymax <- NULL
  if(missing(sigVal)) sigVal <- NULL
  if(missing(sigLine)) sigLine <- NULL
  if(missing(bg)) bg <- NULL
  if(missing(baseline)) baseline <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if bed/pVals arguments are missing (could be in object)
  if(!hasArg(data)) data <- NULL
  if(!hasArg(pVals)) pVals <- NULL

  ## Compile all parameters into an internal object
  bb_manInternal <- structure(list(data = data, pVals = pVals, chrom = chrom, chromstart = chromstart, chromend = chromend, assembly = assembly,
                                   fill = fill, pch = pch, space = space, cex = cex, ymax = ymax, range = range, sigVal = sigVal,
                                   sigLine = sigLine, sigCol = sigCol, bg = bg, baseline = baseline, x = x, y = y, width = width, height = height, just = just,
                                   default.units = default.units, draw = draw), class = "bb_manInternal")

  bb_manInternal <- parseParams(bb_params = params, object_params = bb_manInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_manInternal$assembly)) bb_manInternal$assembly <- "hg19"
  if(is.null(bb_manInternal$fill)) bb_manInternal$fill <- "black"
  if(is.null(bb_manInternal$pch)) bb_manInternal$pch <- 19
  if(is.null(bb_manInternal$space)) bb_manInternal$space <- 0.01
  if(is.null(bb_manInternal$cex)) bb_manInternal$cex <- 0.25
  if(is.null(bb_manInternal$ymax)) bb_manInternal$ymax <- 1
  if(is.null(bb_manInternal$sigVal)) bb_manInternal$sigVal <- 5e-08
  if(is.null(bb_manInternal$sigLine)) bb_manInternal$sigLine <- FALSE
  if(is.null(bb_manInternal$bg)) bb_manInternal$bg <- NA
  if(is.null(bb_manInternal$baseline)) bb_manInternal$baseline <- FALSE
  if(is.null(bb_manInternal$just)) bb_manInternal$just <- c("left", "top")
  if(is.null(bb_manInternal$default.units)) bb_manInternal$default.units <- "inches"
  if(is.null(bb_manInternal$draw)) bb_manInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  man_plot <- structure(list(chrom = bb_manInternal$chrom, chromstart = bb_manInternal$chromstart, chromend = bb_manInternal$chromend, assembly = NULL,
                             range = bb_manInternal$range, ymax = bb_manInternal$ymax, space = bb_manInternal$space, x = bb_manInternal$x, y = bb_manInternal$y,
                             width = bb_manInternal$width, height = bb_manInternal$height, just = bb_manInternal$just, grobs = NULL), class = "bb_manhattan")
  attr(x = man_plot, which = "plotted") <- bb_manInternal$draw

  # ======================================================================================================================================================================================
  # CATCH MISSING ARGUMENT AND PLACEMENT ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_manInternal$data)) stop("argument \"data\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_manInternal$pVals)) stop("argument \"pVals\" is missing, with no default.", call. = FALSE)

  check_placement(object = man_plot)

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  bedfile <- bb_manInternal$data
  ## Read in data if it's not a dataframe or data.table
  if (!"data.frame" %in% class(bedfile)){
    if (!"GRanges" %in% class(bedfile)){
      bedfile <- fread(bedfile)
    }

  }

  bedfile <- as.data.frame(bedfile)

  # ======================================================================================================================================================================================
  # CATCH MORE ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotmanhattan(bedfile = bedfile, chrom = man_plot$chrom, chromstart = man_plot$chromstart, chromend = man_plot$chromend,
                              pVals = bb_manInternal$pVals, object = man_plot,
                              fillcolor = bb_manInternal$fill)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  man_plot$assembly <- parse_bbAssembly(assembly = bb_manInternal$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  man_plot <- defaultUnits(object = man_plot, default.units = bb_manInternal$default.units)

  # ======================================================================================================================================================================================
  # READ AND SUBSET DATA
  # ======================================================================================================================================================================================

  bed_data <- parse_data(bedfile = bedfile, chrom = man_plot$chrom, chromstart = man_plot$chromstart, chromend = man_plot$chromend, pVals = bb_manInternal$pVals)

  if (nrow(bed_data) > 0){

    # ======================================================================================================================================================================================
    # MULTIPLE CHROMOSOMES
    # ======================================================================================================================================================================================

    if (is.null(man_plot$chrom)){

      man_plot$chromstart <- NULL
      man_plot$chromend <- NULL
      chroms <- as.character(unique(bed_data[,1]))

      ## Get chrom sizes based on assembly data
      txdbChecks <- check_loadedPackage(package = man_plot$assembly$TxDb, message = paste(paste0("`", man_plot$assembly$TxDb,"`"), "not loaded. Please install and load to generate full genome assembly Manhattan plot."))
      if (txdbChecks == TRUE){

        tx_db <- eval(parse(text = man_plot$assembly$TxDb))
        assembly_data <- as.data.frame(setDT(as.data.frame(seqlengths(tx_db)), keep.rownames = TRUE))
        assembly_data <- assembly_data[which(assembly_data[,1] %in% chroms),]
        man_plot$chrom <- assembly_data[,1]

        if (any(!chroms %in% assembly_data[,1])){
          non_txdb <- chroms[which(!chroms %in% assembly_data[,1])]
          warning(paste("Chromosome(s)", paste0("'", non_txdb, "'", collapse = ", "), "not found in", paste0("`", man_plot$assembly$TxDb, "`"), "and will be ignored."), call. = FALSE)
        }

        ## get the offsets based on spacer for the assembly
        offsetAssembly <- spaceChroms(assemblyData = assembly_data, space = bb_manInternal$space)

        ## remove bed_data data that aren't in the genome assembly
        bed_data <- bed_data[bed_data[,1] %in% offsetAssembly[,1],]

        ## Add chromosome offsets to bed_data
        bed_data <- bed_offset(bedData = bed_data, offsetAssembly = offsetAssembly)

        ## Set viewport xscale
        cumsums <- cumsum(as.numeric(assembly_data[,2]))
        spacer <- cumsums[length(cumsum(as.numeric(assembly_data[,2])))] * bb_manInternal$space

        xscale <- c(0, max(offsetAssembly[,4]) + spacer)

      } else {
        xscale <- c(0, 1)
      }



    } else {

      # ======================================================================================================================================================================================
      # SINGLE CHROMOSOME
      # ======================================================================================================================================================================================
      offsetAssembly <- NULL

      ## Whole single chromosome needs chromosome length information
      if (is.null(man_plot$chromstart) & is.null(man_plot$chromend)){

        txdbChecks <- check_loadedPackage(package = man_plot$assembly$TxDb, message = paste(paste0("`", man_plot$assembly$TxDb,"`"), "not loaded. Please install and load to generate full chromosome Manhattan plot."))
        if (txdbChecks == TRUE){

          tx_db <- eval(parse(text = man_plot$assembly$TxDb))
          assembly_data <- seqlengths(tx_db)

          if (!man_plot$chrom %in% names(assembly_data)){
            warning(paste("Chromosome", paste0("'", man_plot$chrom, "'"), "not found in", paste0("`", man_plot$assembly$TxDb, "`"), "and data for entire chromosome cannot be plotted."), call. = FALSE)
            xscale <- c(0, 1)
          } else {
            man_plot$chromstart <- 1
            man_plot$chromend <- assembly_data[[man_plot$chrom]]
            xscale <- c(man_plot$chromstart, man_plot$chromend)

          }


        }



      } else {
        xscale <- c(man_plot$chromstart, man_plot$chromend)
        txdbChecks <- TRUE
      }


    }

    if (txdbChecks != FALSE){

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

      color_nonsig <- parse_color(fillcolor = bb_manInternal$fill, offsetAssembly = offsetAssembly, bedData = nonsigBed)

      if (!is.null(bb_manInternal$sigCol)){

        color_sig <- parse_color(fillcolor = bb_manInternal$sigCol, offsetAssembly = offsetAssembly, bedData = sigBed)

      } else {

        color_sig <- parse_color(fillcolor = bb_manInternal$fill, offsetAssembly = offsetAssembly, bedData = sigBed)
      }


      colorBed <- rbind(color_nonsig, color_sig)




    } else {
      xscale <- c(0, 1)
      man_plot$range <- c(0, 1)
    }




  } else {
    txdbChecks <- TRUE
    xscale <- c(0, 1)
    man_plot$range <- c(0, 1)
    warning("No data found in region.", call. = FALSE)

  }


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
  # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
  # ======================================================================================================================================================================================
  backgroundGrob <- rectGrob(gp = gpar(fill = bb_manInternal$bg, col = NA), name = "background")
  assign("manhattan_grobs", gTree(vp = vp, children = gList(backgroundGrob)), envir = bbEnv)

  if (nrow(bed_data) > 0 & txdbChecks == TRUE){

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

    if (bb_manInternal$baseline == TRUE){
      gp$col <- gp$linecolor
      gp$lwd <- 1.5
      baselineGrob <- segmentsGrob(x0 = unit(0, "npc"), y0 = 0, x1 = unit(1, "npc"), y1 = 0,
                              gp = gp, default.units = "native")
      assign("manhattan_grobs", addGrob(gTree = get("manhattan_grobs", envir = bbEnv), child = baselineGrob), envir = bbEnv)

    }
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

  message(paste0("bb_manhattan[", vp$name, "]"))
  invisible(man_plot)

}
