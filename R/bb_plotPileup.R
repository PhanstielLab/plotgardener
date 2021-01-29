#' Plot data elements in a pileup style for a single chromosome
#'
#' @param data Data to be plotted; as a character value specifying a BED file path, a data frame in BED format, a character value specifying a .bam file path where a bam index file (.bam.bai) is in the same directory, or a \link[GenomicRanges]{GRanges} object.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a \link[BentoBox]{bb_assembly} object. Default value is \code{assembly = "hg19"}.
#' @param fill Character value(s) as a single value, vector, or palette specifying fill colors of pileup elements. Default value is \code{fill = "black"}.
#' @param colorby A "\link[BentoBox]{colorby}" object specifying information for scaling colors in data.
#' @param linecolor A character value specifying the color of the lines outlining pileup elements. Default value is \code{linecolor = NA}.
#' @param boxHeight A numeric or unit object specifying height of pileup element boxes. Default value is \code{boxHeight = unit(2, "mm")}.
#' @param spaceWidth A numeric value specifying the width of minimum spacing between pileup element boxes, as a fraction of the plot's genomic range. Default value is \code{spaceWidth = 0.02}.
#' @param spaceHeight A numeric value specifying the height of spacing between pileup element boxes on different rows, as a fraction of boxHeight. Default value is \code{spaceHeight = 0.3}.
#' @param strandSplit A logical value indicating whether plus and minus-stranded elements should be separated. Elements can only be split by strand if a \code{strand} column is found in \code{data}. Default value is \code{strandSplit = FALSE}.
#' @param bg Character value indicating background color. Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a baseline along the x-axis. Default value is \code{baseline = FALSE}.
#' @param x A numeric or unit object specifying pileup plot x-location.
#' @param y A numeric or unit object specifying pileup plot y-location.
#' @param width A numeric or unit object specifying pileup plot width.
#' @param height A numeric or unit object specifying pileup plot height.
#' @param just Justification of pileup plot relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{x}, \code{y}, \code{width}, or \code{height} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be produced. Default value \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_pileup} object containing relevant genomic region, coloring data, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load BED data
#' data("bb_bedData")
#'
#' ## Plot pileup plot filling up entire graphic device
#' bb_plotPileup(data = bb_bedData, chrom = "chr21", chromstart = 29073000, chromend = 29074000,
#'               fill = "grey")
#'
#' ## Plot and place pileup plot on a BentoBox page
#' bb_pageCreate(width = 5, height = 5, default.units = "inches", xgrid = 0, ygrid = 0)
#' bb_plotPileup(data = bb_bedData, chrom = "chr21", chromstart = 29073000, chromend = 29074000,
#'               fill = "grey", baseline = TRUE,
#'               x = 0.5, y = 0.5, width = 4, height = 4,
#'               just = c("left", "top"), default.units = "inches")
#'
#' @details
#' This function can be used to quickly plot a pileup plot by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotPileup(data, chrom,
#'               chromstart = NULL, chromend = NULL)
#' }
#' A pileup plot can be placed on a BentoBox coordinate page by providing plot placement parameters:
#' \preformatted{
#' bb_plotPileup(data, chrom,
#'               chromstart = NULL, chromend = NULL,
#'               x, y, width, height, just = c("left", "top"),
#'               default.units = "inches")
#' }
#'
#' @export
bb_plotPileup <- function(data, chrom, chromstart = NULL, chromend = NULL, assembly = "hg19", fill = "black", colorby = NULL, linecolor = NA,
                          boxHeight =  unit(2, "mm"), spaceWidth = 0.02, spaceHeight = 0.3, strandSplit = FALSE, bg = NA, baseline = FALSE, x = NULL, y = NULL, width = NULL,
                          height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, params = NULL, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors
  errorcheck_bb_plotpileup <- function(pileup_plot, colorby){


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

      if (!any(colnames(bed) == colorby$column)){

        stop(paste("Colorby column", paste0('`', colorby$column, '`'), "not found in data. Check colorby column name."), call. = FALSE)
      }

      if (length(which(colnames(bed) == colorby$column)) > 1){
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
  if(missing(fill)) fill <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(strandSplit)) strandSplit <- NULL
  if(missing(boxHeight)) boxHeight <- NULL
  if(missing(spaceHeight)) spaceHeight <- NULL
  if(missing(spaceWidth)) spaceWidth <- NULL
  if(missing(bg)) bg <- NULL
  if(missing(baseline)) baseline <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if data/chrom arguments are missing (could be in object)
  if(!hasArg(data)) data <- NULL
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_pileInternal <- structure(list(data = data, chrom = chrom, chromstart = chromstart, chromend = chromend, assembly = assembly, fill = fill, linecolor = linecolor,
                                    colorby = colorby, strandSplit = strandSplit, boxHeight = boxHeight, spaceHeight = spaceHeight, spaceWidth = spaceWidth, bg = bg,
                                    baseline = baseline, x = x, y = y, width = width, height = height, just = just, default.units = default.units, draw = draw), class = "bb_pileInternal")

  bb_pileInternal <- parseParams(bb_params = params, object_params = bb_pileInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_pileInternal$assembly)) bb_pileInternal$assembly <- "hg19"
  if(is.null(bb_pileInternal$fill)) bb_pileInternal$fill <- "black"
  if(is.null(bb_pileInternal$linecolor)) bb_pileInternal$linecolor <- NA
  if(is.null(bb_pileInternal$strandSplit)) bb_pileInternal$strandSplit <- FALSE
  if(is.null(bb_pileInternal$boxHeight)) bb_pileInternal$boxHeight <- unit(2, "mm")
  if(is.null(bb_pileInternal$spaceHeight)) bb_pileInternal$spaceHeight <- 0.3
  if(is.null(bb_pileInternal$spaceWidth)) bb_pileInternal$spaceWidth <- 0.02
  if(is.null(bb_pileInternal$bg)) bb_pileInternal$bg <- NA
  if(is.null(bb_pileInternal$baseline)) bb_pileInternal$baseline <- FALSE
  if(is.null(bb_pileInternal$just)) bb_pileInternal$just <- c("left", "top")
  if(is.null(bb_pileInternal$default.units)) bb_pileInternal$default.units <- "inches"
  if(is.null(bb_pileInternal$draw)) bb_pileInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # CHECK ARGUMENT ERROS
  # ======================================================================================================================================================================================

  if(is.null(bb_pileInternal$data)) stop("argument \"data\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_pileInternal$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)

  if (!is.null(bb_pileInternal$colorby)){
    if(class(bb_pileInternal$colorby) != "bb_colorby"){
      stop("\"colorby\" not of class \"bb_colorby\". Input colorby information with \"colorby()\".", call. = FALSE)
    }
  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  pileup_plot <- structure(list(chrom = bb_pileInternal$chrom, chromstart = bb_pileInternal$chromstart, chromend = bb_pileInternal$chromend, assembly = bb_pileInternal$assembly,
                                color_palette = NULL, zrange = bb_pileInternal$colorby$range, x = bb_pileInternal$x, y = bb_pileInternal$y, width = bb_pileInternal$width, height = bb_pileInternal$height,
                                just = bb_pileInternal$just, grobs = NULL), class = "bb_pileup")
  attr(x = pileup_plot, which = "plotted") <- bb_pileInternal$draw

  # ======================================================================================================================================================================================
  # CHECK PLACEMENT ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = pileup_plot)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  pileup_plot$assembly <- parse_bbAssembly(assembly = pileup_plot$assembly)

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

  bed <- bb_pileInternal$data
  if (!"data.frame" %in% class(bed)){
    if (!"GRanges" %in% class(bed)){

      if (file_ext(bed) == "bam"){
        indexFile <- paste0(bed, ".bai")
        if (!file.exists(indexFile)){
          stop("Cannot read in bam file without a corresponding bam index file (.bai) in the same directory.", call. = FALSE)
        }
        bed <- read_bam(bed) %>% filter_by_overlaps(GRanges(seqnames = pileup_plot$chrom, ranges = IRanges(start = pileup_plot$chromstart, end = pileup_plot$chromend))) %>% mutate()
      } else {
        bed <- fread(bed)
      }

    }

  }

  bed <- as.data.frame(bed)


  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotpileup(pileup_plot = pileup_plot, colorby = bb_pileInternal$colorby)

  # ======================================================================================================================================================================================
  # WHOLE CHROMOSOME DATA AND XSCALE
  # ======================================================================================================================================================================================

  if (is.null(pileup_plot$chromstart) & is.null(pileup_plot$chromend)){

    txdbChecks <- check_loadedPackage(package = pileup_plot$assembly$TxDb, message = paste(paste0("`", pileup_plot$assembly$TxDb,"`"),
                                                                                          "not loaded. Please install and load to plot full chromosome pileup plot."))
    xscale <- c(0, 1)
    if (txdbChecks == TRUE){
      tx_db <- eval(parse(text = pileup_plot$assembly$TxDb))
      assembly_data <- seqlengths(tx_db)

      if (!pileup_plot$chrom %in% names(assembly_data)){
        txdbChecks <- FALSE
        warning(paste("Chromosome", paste0("'", pileup_plot$chrom, "'"), "not found in", paste0("`", pileup_plot$assembly$TxDb, "`"), "and data for entire chromosome cannot be plotted."), call. = FALSE)
      } else {
        pileup_plot$chromstart <- 1
        pileup_plot$chromend <- assembly_data[[pileup_plot$chrom]]
        xscale <- c(pileup_plot$chromstart, pileup_plot$chromend)

      }

    }

  } else {
    txdbChecks <- TRUE
    xscale <- c(pileup_plot$chromstart, pileup_plot$chromend)
  }

  # ======================================================================================================================================================================================
  # SUBSET DATA FOR CHROMOSOME AND ANY OVERLAPPING REGIONS
  # ======================================================================================================================================================================================

  if (!is.null(pileup_plot$chromstart) & !is.null(pileup_plot$chromend)){
    bed <- bed[which(bed[,1] == pileup_plot$chrom & bed[,2] <= pileup_plot$chromend & bed[,3] >= pileup_plot$chromstart),]
  } else {
    bed <- data.frame(matrix(nrow = 0, ncol = 3))
  }

  # ======================================================================================================================================================================================
  # SET COLORBY DATA
  # ======================================================================================================================================================================================

  if (!is.null(bb_pileInternal$colorby) & nrow(bed) > 0){
    colorbyCol <- which(colnames(bed) == bb_pileInternal$colorby$column)
    colorbyCol <- bed[,colorbyCol]

    ## if the associated column isn't numbers, convert unique values to a set of numbers
    if (class(colorbyCol) != "numeric" | class(colorbyCol) != "integer"){
      colorbyCol <- factor(colorbyCol)
      bed$colorby <- as.numeric(colorbyCol)
    } else {

      bed$colorby <- colorbyCol
    }

    if (is.null(bb_pileInternal$colorby$range)){
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

    ## Look for column named 'strand'
    if (!any(colnames(bed) == "strand")){
      stop("No `strand` column found in data. Cannot split data based on strand.", call. = FALSE)
    } else {
      strand_col <- which(colnames(bed) == "strand")
      posStrand <- bed[which(bed[,strand_col] == "+"),]
      minStrand <- bed[which(bed[,strand_col] == "-"),]
    }

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
                   xscale = xscale,
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
                   xscale = xscale,
                   yscale = yscale,
                   just = bb_pileInternal$just,
                   name = vp_name)
  }


  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
  # ======================================================================================================================================================================================

  backgroundGrob <- rectGrob(gp = gpar(fill = bb_pileInternal$bg, col = NA), name = "background")
  assign("pileup_grobs", gTree(vp = vp, children = gList(backgroundGrob)), envir = bbEnv)

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

  maxRows <- floor((as.numeric(vp$height) + spaceHeight)/(boxHeight + spaceHeight))
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
      posDF <- checkRow(posMatrix, maxRows*0.5, 3, wiggle)
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
      minDF <- checkRow(minMatrix, maxRows*0.5, 3, wiggle)
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

  if (nrow(rowDF) > 0){

    # ======================================================================================================================================================================================
    # COLORS
    # ======================================================================================================================================================================================

    if (is.null(bb_pileInternal$colorby)){

      if (class(bb_pileInternal$fill) == "function"){
        colors <- bb_pileInternal$fill(maxRows)
        indeces <- rowDF$row
        rowDF$color <- colors[indeces]

      } else {

        if (length(bb_pileInternal$fill) == 1){
          rowDF$color <- rep(bb_pileInternal$fill, nrow(rowDF))
        } else {

          colors <- rep(bb_pileInternal$fill, ceiling(maxRows/length(bb_pileInternal$fill)))[1:maxRows]
          indeces <- rowDF$row
          rowDF$color <- colors[indeces]

        }

      }

    } else {

      if (class(bb_pileInternal$fill) == "function"){

        rowDF$color <- bb_maptocolors(rowDF$colorby, bb_pileInternal$fill, range = pileup_plot$zrange)
        #sorted_colors <- unique(rowDF[order(rowDF$colorby),]$color)
        pileup_plot$color_palette <- bb_pileInternal$fill

      } else {

        colorbyCol <- factor(rowDF$colorby)
        mappedColors <- rep(bb_pileInternal$fill, ceiling(length(levels(colorbyCol))/length(bb_pileInternal$fill)))
        rowDF$color <- mappedColors[rowDF$colorby]
      }


    }


    # ======================================================================================================================================================================================
    # MAKE GROBS
    # ======================================================================================================================================================================================


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

    if (bb_pileInternal$strandSplit == TRUE | bb_pileInternal$baseline == TRUE){

      lineGrob <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                               y0 = unit(0, "native"), y1 = unit(0, "native"),
                               gp = gpar(...))
      assign("pileup_grobs", addGrob(gTree = get("pileup_grobs", envir = bbEnv), child = lineGrob), envir = bbEnv)

    }


  } else {

    if (txdbChecks == TRUE){
      warning("No pileup data to plot.", call. = FALSE)
    }

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
