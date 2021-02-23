#' Plot gene transcripts in a pileup style for a single chromosome
#'
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a \link[BentoBox]{bb_assembly} object. Default value is \code{assembly = "hg19"}.
#' @param fill Character value(s) as a single value or vector specifying fill colors of transcripts. Default value is \code{fill = c("#669fd9", "#abcc8e")}.
#' @param colorbyStrand A logical value indicating whether to color plus and minus strands by the first two colors in a \code{fill} vector, where plus strand transcripts will be colored by the first \code{fill} color and
#' minus strand transcripts will be colored by the second \code{fill} color. Default value is \code{colorbyStrand = TRUE}.
#' @param strandSplit A logical value indicating whether plus and minus-stranded transcripts should be separated, with plus strand transcripts plotted above the x-axis and minus strand transcripts plotted below the x-axis. Default value is \code{strandSplit = FALSE}.
#' @param boxHeight A numeric or unit object specifying height of transcripts. Default value is \code{boxHeight = unit(2, "mm")}.
#' @param spaceWidth A numeric value specifying the width of minimum spacing between transcripts, as a fraction of the plot's genomic range. Default value is \code{spaceWidth = 0.02}.
#' @param spaceHeight A numeric value specifying the height of spacing between transcripts on different rows, as a fraction of \code{boxHeight}. Default value is \code{spaceHeight = 0.3}.
#' @param fontsize A numeric specifying text fontsize in points. Default value is \code{fontsize = 8}.
#' @param labels A character value describing the format of transcript text labels. Default value is \code{labels = "trancript"}. Options are:
#' \itemize{
#' \item{\code{NULL}: }{No labels.}
#' \item{\code{"transcript"}: }{Transcript name labels.}
#' \item{\code{"gene"}: }{Gene name labels.}
#' \item{\code{"both"}: }{Combined transcript and gene name labels with the format "gene name:transcript name".}
#' }
#' @param stroke A numeric value indicating the stroke width for transcript body outlines. Default value is \code{stroke = 0.1}.
#' @param bg Character value indicating background color. Default value is \code{bg = NA}.
#' @param x A numeric or unit object specifying transcript plot x-location.
#' @param y A numeric, unit object, or character containing a "b" combined with a numeric value specifying transcript plot y-location. The character value will
#' place the transcript plot y relative to the bottom of the most recently plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying transcript plot width.
#' @param height A numeric or unit object specifying transcript plot height.
#' @param just Justification of transcript plot relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{x}, \code{y}, \code{width}, or \code{height} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#'
#' @return Returns a \code{bb_transcripts} object containing relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load hg19 genomic annotation packages
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Create page
#' bb_pageCreate(width = 7.5, height = 3.5, default.units = "inches")
#'
#' ## Plot and place transcripts
#' bb_plotTranscripts(chrom = "chr8", chromstart = 1000000, chromend = 2000000,
#'                    assembly = "hg19", labels = "gene",
#'                    x = 0.5, y = 0.5, width = 6.5, height = 2.5,
#'                    just = c("left", "top"), default.units = "inches")
#'
#' ## Plot genome label
#' bb_plotGenomeLabel(chrom = "chr8", chromstart = 1000000, chromend = 2000000,
#'                    x = 0.5, y = 3.03, length = 6.5, default.units = "inches")
#'
#' ## Plot a legend
#' bb_plotLegend(legend = c("+ strand", "- strand"), fill = c("#669fd9", "#abcc8e"), border = FALSE,
#'               x = 0.5, y = 1, width = 1, height = 0.5, just = c("left", "top"))
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @details
#' A transcripts plot can be placed on a BentoBox coordinate page by providing plot placement parameters:
#' \preformatted{
#' bb_plotTranscripts(chrom, chromstart = NULL, chromend = NULL,
#'                    x, y, width, height, just = c("left", "top"),
#'                    default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated transcripts plot by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotTranscripts(chrom, chromstart = NULL, chromend = NULL)
#' }
#'
#' Genomic annotation information is acquired through \link[GenomicFeatures]{TxDb} and \link[AnnotationDbi]{OrgDb-class} packages, as determined
#' through the \code{assembly} parameter.
#'
#' @seealso \link[BentoBox]{bb_assembly}, \link[BentoBox]{bb_genomes}, \link[BentoBox]{bb_defaultPackages}
#'
#' @export
bb_plotTranscripts <- function(chrom, chromstart = NULL, chromend = NULL, assembly = "hg19", fill = c("#669fd9", "#abcc8e"), colorbyStrand = TRUE, strandSplit = FALSE,
                               boxHeight = unit(2, "mm"), spaceWidth = 0.02, spaceHeight = 0.3, fontsize = 8, labels = "transcript", stroke = 0.1, bg = NA,
                               x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, params = NULL){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that checks errors for bb_plotTranscripts
  errorcheck_bb_plotTranscripts <- function(transcript_plot, labels){

    ## Can't have only one NULL chromstart or chromend
    if ((is.null(transcript_plot$chromstart) & !is.null(transcript_plot$chromend)) | (is.null(transcript_plot$chromend) & !is.null(transcript_plot$chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.", call. = FALSE)

    }


    if (!is.null(transcript_plot$chromstart) & !is.null(transcript_plot$chromend)){

      ## chromend > chromstart
      if (transcript_plot$chromend < transcript_plot$chromstart){

        stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)


      }

    }


    if (!labels %in% c(NULL, "transcript", "gene", "both")){
      stop("Invalid \'labels\' input. Options are \'NULL\', \'transcript\', \'gene\', or \'both\'.", call. = FALSE)
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

  ## Define a function that makes utr grobs
  utr_grobs <- function(df, boxHeight){

    utr_ranges <- as.list(strsplit(as.character(df[6]), ",")[[1]])

    if (length(utr_ranges) > 0){

      starts <- lapply(utr_ranges, parse_starts)
      widths <- lapply(utr_ranges, parse_widths)
      utrs_dataframe <- cbind(unlist(starts), unlist(widths))

      utrs <- rectGrob(x = utrs_dataframe[,1],
                       y = as.numeric(df[10]) + 0.5*boxHeight,
                       width = utrs_dataframe[,2],
                       height = boxHeight*0.65,
                       just = "left",
                       gp = gpar(fill = df[8], col = NA),
                       default.units = "native")

      assign("transcript_grobs", addGrob(get("transcript_grobs", envir = bbEnv), child = utrs), envir = bbEnv)
    }


  }

  ## Define a function that determines if a label will be cut off
  cutoffLabel <- function(df, fontsize, xscale, vp, unit){

    label <- df[1]
    location <- df[2]

    if (unit == "npc"){
      downViewport(name = vp$name)
      labelWidth <- convertWidth(widthDetails(textGrob(label = label, gp = gpar(fontsize = fontsize))), unitTo = "native", valueOnly = T)
      upViewport()
    } else {
      pushViewport(vp)
      labelWidth <- convertWidth(widthDetails(textGrob(label = label, gp = gpar(fontsize = fontsize))), unitTo = "native", valueOnly = T)
      upViewport()
    }


    leftBound <- as.numeric(location) - 0.5*labelWidth
    rightBound <- as.numeric(location) + 0.5*labelWidth

    if (leftBound < xscale[1] | rightBound > xscale[2]){
      return(NA)
    } else {
      return(label)
    }

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(boxHeight)) boxHeight <- NULL
  if(missing(spaceHeight)) spaceHeight <- NULL
  if(missing(spaceWidth)) spaceWidth <- NULL
  if(missing(fill)) fill <- NULL
  if(missing(colorbyStrand)) colorbyStrand <- NULL
  if(missing(labels)) labels <- NULL
  if(missing(fontsize)) fontsize <- NULL
  if(missing(strandSplit)) strandSplit <- NULL
  if(missing(stroke)) stroke <- NULL
  if(missing(bg)) bg <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if chrom argument is missing (could be in object)
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_transcriptsInternal <- structure(list(assembly = assembly, chrom = chrom, chromstart = chromstart, chromend = chromend, boxHeight = boxHeight, spaceHeight = spaceHeight,
                                     spaceWidth = spaceWidth, fill = fill, colorbyStrand = colorbyStrand, labels = labels, fontsize = fontsize,
                                     strandSplit = strandSplit, stroke = stroke, bg = bg, x = x, y = y, width = width, height = height, just = just, default.units = default.units,
                                     draw = draw), class = "bb_transcriptsInternal")

  bb_transcriptsInternal <- parseParams(bb_params = params, object_params = bb_transcriptsInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_transcriptsInternal$assembly)) bb_transcriptsInternal$assembly <- "hg19"
  if(is.null(bb_transcriptsInternal$boxHeight)) bb_transcriptsInternal$boxHeight <- unit(2, "mm")
  if(is.null(bb_transcriptsInternal$spaceHeight)) bb_transcriptsInternal$spaceHeight <- 0.3
  if(is.null(bb_transcriptsInternal$spaceWidth)) bb_transcriptsInternal$spaceWidth <- 0.02
  if(is.null(bb_transcriptsInternal$fill)) bb_transcriptsInternal$fill <- c("#669fd9", "#abcc8e")
  if(is.null(bb_transcriptsInternal$colorbyStrand)) bb_transcriptsInternal$colorbyStrand <- TRUE
  if(is.null(bb_transcriptsInternal$labels)) bb_transcriptsInternal$labels <- "transcript"
  if(is.null(bb_transcriptsInternal$fontsize)) bb_transcriptsInternal$fontsize <- 8
  if(is.null(bb_transcriptsInternal$strandSplit)) bb_transcriptsInternal$strandSplit <- FALSE
  if(is.null(bb_transcriptsInternal$stroke)) bb_transcriptsInternal$stroke <- 0.1
  if(is.null(bb_transcriptsInternal$bg)) bb_transcriptsInternal$bg <- NA
  if(is.null(bb_transcriptsInternal$just)) bb_transcriptsInternal$just <- c("left", "top")
  if(is.null(bb_transcriptsInternal$default.units)) bb_transcriptsInternal$default.units <- "inches"
  if(is.null(bb_transcriptsInternal$draw)) bb_transcriptsInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_transcripts <- structure(list(chrom = bb_transcriptsInternal$chrom, chromstart = bb_transcriptsInternal$chromstart, chromend = bb_transcriptsInternal$chromend, assembly = bb_transcriptsInternal$assembly,
                                   x = bb_transcriptsInternal$x, y = bb_transcriptsInternal$y, width = bb_transcriptsInternal$width, height = bb_transcriptsInternal$height,
                                   just = bb_transcriptsInternal$just, grobs = NULL), class = "bb_transcripts")
  attr(x = bb_transcripts, which = "plotted") <- bb_transcriptsInternal$draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_transcripts$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)
  check_placement(object = bb_transcripts)
  errorcheck_bb_plotTranscripts(transcript_plot = bb_transcripts, labels = bb_transcriptsInternal$labels)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  bb_transcripts$assembly <- parse_bbAssembly(assembly = bb_transcripts$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  bb_transcripts <- defaultUnits(object = bb_transcripts, default.units = bb_transcriptsInternal$default.units)
  if(!"unit" %in% class(bb_transcriptsInternal$boxHeight)){

    if (!is.numeric(bb_transcriptsInternal$boxHeight)){

      stop("\'boxHeight\' is neither a unit object or a numeric value. Cannot make transcript plot.", call. = FALSE)

    }

    if (is.null(bb_transcriptsInternal$default.units)){

      stop("\'boxHeight\' detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_transcriptsInternal$boxHeight <- unit(bb_transcriptsInternal$boxHeight, bb_transcriptsInternal$default.units)
  }


  # ======================================================================================================================================================================================
  # GET APPROPRIATE BUILD DATA
  # ======================================================================================================================================================================================

  txdbChecks <- check_loadedPackage(package = bb_transcripts$assembly$TxDb, message = paste(paste0("`", bb_transcripts$assembly$TxDb,"`"),
                                                                                      "not loaded. Please install and load to plot gene transcripts."))
  orgdbChecks <- check_loadedPackage(package = bb_transcripts$assembly$OrgDb, message = paste(paste0("`", bb_transcripts$assembly$OrgDb,"`"),
                                                                                        "not loaded. Please install and load to plot gene transcripts."))
  data <- data.frame(matrix(ncol = 22, nrow = 0))
  xscale <- c(0, 1)
  if (txdbChecks == TRUE & orgdbChecks == TRUE){

    tx_db <- eval(parse(text = bb_transcripts$assembly$TxDb))
    genome <- seqlengths(tx_db)
    displayCol <- bb_transcripts$assembly$display.column

    if (!bb_transcripts$chrom %in% names(genome)){
      warning(paste("Chromosome", paste0("'", bb_transcripts$chrom, "'"), "not found in", paste0("`", bb_transcripts$assembly$TxDb, "`"), "and gene transcripts cannot be plotted."), call. = FALSE)
    } else {
      if (is.null(bb_transcripts$chromstart) & is.null(bb_transcripts$chromend)){
        bb_transcripts$chromstart <- 1
        bb_transcripts$chromend <- genome[[bb_transcripts$chrom]]
      }

      data <- bb_getExons(assembly = bb_transcripts$assembly, chromosome = bb_transcripts$chrom, start = bb_transcripts$chromstart, stop = bb_transcripts$chromend)
      xscale <- c(bb_transcripts$chromstart, bb_transcripts$chromend)

    }

  }

  # ======================================================================================================================================================================================
  # COLORS
  # ======================================================================================================================================================================================

  if (bb_transcriptsInternal$colorbyStrand == TRUE){
    if (length(bb_transcriptsInternal$fill) == 1){
      posCol <- bb_transcriptsInternal$fill
      negCol <- bb_transcriptsInternal$fill
    } else {
      posCol <- bb_transcriptsInternal$fill[1]
      negCol <- bb_transcriptsInternal$fill[2]
    }

    pos <- data[which(data$TXSTRAND == "+"),]
    pos$color <- rep(posCol, nrow(pos))
    neg <- data[which(data$TXSTRAND == "-"),]
    neg$color <- rep(negCol, nrow(neg))
    data <- rbind(pos, neg)

  } else {

    data$color <- rep(bb_transcriptsInternal$fill[1], nrow(data))
  }


  # ======================================================================================================================================================================================
  # SEPARATE DATA INTO STRANDS
  # ======================================================================================================================================================================================

  if (bb_transcriptsInternal$strandSplit == TRUE){

    posStrand <- data[which(data$TXSTRAND == "+"),]
    minStrand <- data[which(data$TXSTRAND == "-"),]

  }


  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_transcripts", length(grep(pattern = "bb_transcripts", x = currentViewports)) + 1)

  if (is.null(bb_transcripts$x) & is.null(bb_transcripts$y)){

    height <- 0.5
    yscale <- strand_scale(strandSplit = bb_transcriptsInternal$strandSplit, height = height)

    vp <- viewport(height = unit(0.5, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = xscale,
                   yscale = yscale,
                   just = "center",
                   name = vp_name)


    if (bb_transcriptsInternal$draw == TRUE){

      vp$name <- "bb_transcripts1"
      grid.newpage()
    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = bb_transcripts)

    height <- convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)
    yscale <- strand_scale(strandSplit = bb_transcriptsInternal$strandSplit, height = height)

    ## Make viewport for gene track
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = xscale,
                   yscale = yscale,
                   just = bb_transcriptsInternal$just,
                   name = vp_name)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
  # ======================================================================================================================================================================================

  backgroundGrob <- rectGrob(gp = gpar(fill = bb_transcriptsInternal$bg, col = NA), name = "background")
  assign("transcript_grobs", gTree(vp = vp, children = gList(backgroundGrob)), envir = bbEnv)

  # ======================================================================================================================================================================================
  # DETERMINE ROWS FOR EACH ELEMENT
  # ======================================================================================================================================================================================

  ## Determine how many rows are going to fit based on boxHeight, spaceHeight, and fontsize
  if (is.null(bb_transcriptsInternal$labels)){
    textHeight <- unit(0, "npc")
  } else {
    textHeight <- heightDetails(textGrob(label = "A", gp = gpar(fontsize = bb_transcriptsInternal$fontsize)))
  }


  if (is.null(bb_transcripts$x) & is.null(bb_transcripts$y)){

    pushViewport(vp)
    boxHeight <- convertHeight(bb_transcriptsInternal$boxHeight, unitTo = "npc", valueOnly = T)
    spaceHeight <- boxHeight*(bb_transcriptsInternal$spaceHeight)
    textHeight <- convertHeight(textHeight, unitTo = "npc", valueOnly = T)
    upViewport()
    unit <- "npc"

  } else {

    boxHeight <- convertHeight(bb_transcriptsInternal$boxHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    spaceHeight <- boxHeight*(bb_transcriptsInternal$spaceHeight)
    textHeight <- convertHeight(textHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    unit <- get("page_units", envir = bbEnv)

  }

  maxRows <- floor(height/(boxHeight + spaceHeight + textHeight + 0.25*textHeight))
  wiggle <- abs(bb_transcripts$chromend - bb_transcripts$chromstart) * bb_transcriptsInternal$spaceWidth

  if (bb_transcriptsInternal$strandSplit == FALSE){

    if (nrow(data) > 0){

      ## Get one representative row per transcript
      repData <- data[duplicated(data$TXNAME) == F,]

      repData$row <- 0

      repData$length <- repData$TXEND - repData$TXSTART

      ## Access default transcript prioritization based on citation/transcript length
      repData <- defaultGenePriorities(data = repData, assembly = bb_transcripts$assembly, transcript = TRUE)

      ## Convert to numeric matrix with txstart, txend, and row for Rcpp function parsing
      dataMatrix <- as.matrix(repData[c("TXSTART", "TXEND", "row")])

      ## Assign a row for each representative element
      rowData <- checkRow(dataMatrix, maxRows, 2, wiggle)
      rowData <- as.data.frame(rowData)

      ## Recombine assigned rows with original data
      rowData <- cbind(rowData[c("row")], repData$TXNAME)
      colnames(rowData) <- c("row", "TXNAME")
      rowData <- suppressMessages(dplyr::left_join(x = data, y = rowData, by = "TXNAME"))


      if (any(rowData$row == 0)){
        rowData <- rowData[which(rowData$row != 0),]
        warning("Not all transcripts shown.", call. = FALSE)

        limitGrob <- textGrob(label = "+", x = unit(1, "npc"), y = unit(1, "npc"),
                              just = c("right", "top"), gp = gpar(col = "black"))
        assign("transcript_grobs", addGrob(gTree = get("transcript_grobs", envir = bbEnv), child = limitGrob), envir = bbEnv)

      }

      ## Change row index to 0 for y-coordinate setting
      rowData$row <- rowData$row - 1
      rowData$y <- rowData$row*(boxHeight + spaceHeight + textHeight + 0.25*textHeight)

      ## Reset rows for colors
      rowData$row <- rowData$row + 1

    } else {
      rowData <- data.frame()
    }

  } else {

    if (nrow(posStrand) > 0){

      ## Get one representative row per transcript
      repPosData <- posStrand[duplicated(posStrand$TXNAME) == F,]

      repPosData$row <- 0

      repPosData$length <- repPosData$TXEND - repPosData$TXSTART
      ## Access default transcript prioritization based on citation/transcript length
      repPosData <- defaultGenePriorities(data = repPosData, assembly = bb_transcripts$assembly, transcript = TRUE)

      ## Convert to numeric matrix with txstart, txend, and row for Rcpp function parsing
      posMatrix <- as.matrix(repPosData[c("TXSTART", "TXEND", "row")])

      ## Assign a row for each representative element
      posRowData <- checkRow(posMatrix, floor(maxRows/2), 2, wiggle)
      posRowData <- as.data.frame(posRowData)

      ## Recombine assigned rows with original data
      posRowData <- cbind(posRowData[c("row")], repPosData$TXNAME)
      colnames(posRowData) <- c("row", "TXNAME")
      posRowData <- suppressMessages(dplyr::left_join(x = posStrand, y = posRowData, by = "TXNAME"))

      if (any(posRowData$row == 0)){
        posRowData <- posRowData[which(posRowData$row != 0),]
        warning("Not all plus strand transcripts shown.", call. = FALSE)
        limitGrob1 <- textGrob(label = "+", x = unit(1, "npc"), y = unit(1, "npc"),
                               just = c("right", "top"), gp = gpar(col = "black"))
        assign("transcript_grobs", addGrob(gTree = get("transcript_grobs", envir = bbEnv), child = limitGrob1), envir = bbEnv)
      }

      ## Set row index to 0 for y-coordinate setting
      posRowData$row <- posRowData$row - 1
      posRowData$y <- (0.5*spaceHeight) + posRowData$row*(boxHeight + spaceHeight + textHeight + 0.25*textHeight)

      ## Reset rows for colors
      posRowData$row <- posRowData$row + 1
      posRowData$row <- posRowData$row + floor(maxRows/2)

    } else {
      posRowData <- data.frame()
    }

    if (nrow(minStrand) > 0){

      ## Get one representative row per transcript
      repMinData <- minStrand[duplicated(minStrand$TXNAME) == F,]

      repMinData$row <- 0

      repMinData$length <- repMinData$TXEND - repMinData$TXSTART
      ## Access default transcript prioritization based on citation/transcript length
      repMinData <- defaultGenePriorities(data = repMinData, assembly = bb_transcripts$assembly, transcript = TRUE)
      ## Convert to numeric matrix with txstart, txend, and row for Rcpp function parsing
      minMatrix <- as.matrix(repMinData[c("TXSTART", "TXEND", "row")])

      ## Assign a row for each representative element
      minRowData <- checkRow(minMatrix, floor(maxRows/2), 2, wiggle)
      minRowData <- as.data.frame(minRowData)

      ## Recombine assigned rows with original data
      minRowData <- cbind(minRowData[c("row")], repMinData$TXNAME)
      colnames(minRowData) <- c("row", "TXNAME")
      minRowData <- suppressMessages(dplyr::left_join(x = minStrand, y = minRowData, by = "TXNAME"))

      if (any(minRowData$row == 0)){
        minRowData <- minRowData[which(minRowData$row != 0),]
        warning("Not all minus strand transcripts shown.", call. = FALSE)
        limitGrob2 <- textGrob(label = "+", x = unit(1, "npc"), y = unit(0, "npc"),
                               just = c("right", "bottom"), gp = gpar(col = "black"))
        assign("transcript_grobs", addGrob(gTree = get("transcript_grobs", envir = bbEnv), child = limitGrob2), envir = bbEnv)

      }

      ## Set row index to 0 for y-coordinate setting
      minRowData$row <- minRowData$row - 1
      minRowData$y <- ((0.5*spaceHeight + boxHeight + textHeight + 0.25*textHeight) + minRowData$row*(boxHeight + spaceHeight + textHeight + 0.25*textHeight))*-1

      ## Reset rows for colors
      rowIndex <- minRowData$row + 1
      rowRange <- floor(maxRows/2):1
      minRowData$row <- rowRange[rowIndex]

    } else {
      minRowData <- data.frame()
    }

    rowData <- rbind(posRowData, minRowData)

  }

  # ======================================================================================================================================================================================
  # UPDATE COLORS IF NECESSARY
  # ======================================================================================================================================================================================

  if (bb_transcriptsInternal$colorbyStrand == FALSE & length(bb_transcriptsInternal$fill) > 1){

    colors <- rep(bb_transcriptsInternal$fill, ceiling(maxRows/length(bb_transcriptsInternal$fill)))[1:maxRows]
    indeces <- rowData$row
    rowData$color <- colors[indeces]

  }


  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (nrow(rowData) > 0){

    ##########################################################
    ## TRANSCRIPT LINES ONLY
    ##########################################################

    if ((bb_transcripts$chromend - bb_transcripts$chromstart) >= 25000000){

      transcriptLine <- rectGrob(x = rowData$TXSTART,
                                 y = rowData$y,
                                 width = rowData$TXEND - rowData$TXSTART,
                                 height = boxHeight,
                                 just = c("left", "bottom"),
                                 gp = gpar(fill = rowData$color, col = rowData$color, lwd = bb_transcriptsInternal$stroke),
                                 default.units = "native")

    } else {

      ##########################################################
      ## TRANSCRIPT LINES
      ##########################################################

      transcriptLine <- rectGrob(x = rowData$TXSTART,
                                 y = rowData$y + 0.5*boxHeight,
                                 width = rowData$TXEND - rowData$TXSTART,
                                 height = boxHeight*0.2,
                                 just = "left",
                                 gp = gpar(fill = rowData$color, col = rowData$color, lwd = bb_transcriptsInternal$stroke),
                                 default.units = "native")

      ##########################################################
      ## TRANSCRIPT EXONS
      ##########################################################

      transcriptExons <- rectGrob(x = rowData$EXONSTART,
                                 y = rowData$y + 0.5*boxHeight,
                                 width = rowData$EXONEND - rowData$EXONSTART,
                                 height = boxHeight*0.65,
                                 just = "left",
                                 gp = gpar(fill = rowData$color, col = NA),
                                 default.units = "native")
      assign("transcript_grobs", addGrob(get("transcript_grobs", envir = bbEnv), child = transcriptExons), envir = bbEnv)

      ##########################################################
      ## TRANSCRIPT CDS
      ##########################################################

      ## Get CDS regions that aren't NA
      cdsData <- rowData[which(!is.na(rowData$CDSSTART)),]
      if (nrow(cdsData) > 0){

        transcriptCds <- rectGrob(x = cdsData$CDSSTART,
                                  y = cdsData$y + 0.5*boxHeight,
                                  width = cdsData$CDSEND - cdsData$CDSSTART,
                                  height = boxHeight,
                                  just = "left",
                                  gp = gpar(fill = cdsData$color, col = cdsData$color, lwd = 1.25),
                                  default.units = "native")
        assign("transcript_grobs", addGrob(get("transcript_grobs", envir = bbEnv), child = transcriptCds), envir = bbEnv)

      }

    }

    assign("transcript_grobs", addGrob(get("transcript_grobs", envir = bbEnv), child = transcriptLine), envir = bbEnv)


    ##########################################################
    ## TRANSCRIPT NAME LABELS
    ##########################################################

    if (!is.null(bb_transcriptsInternal$labels)){

      ## Get representative transcript names
      transcriptNames <- rowData[duplicated(rowData$TXNAME) == F,]

      ## Add column with center location of each transcript label
      transcriptNames$labelLoc <- rowMeans(transcriptNames[c("TXSTART", "TXEND")])

      if (bb_transcriptsInternal$labels == "transcript"){
        label <- transcriptNames$TXNAME
      } else if (bb_transcriptsInternal$labels == "gene"){
        label <- transcriptNames[[bb_transcripts$assembly$display.column]]
      } else {
        label <- paste0(transcriptNames[[bb_transcripts$assembly$display.column]], ":", transcriptNames$TXNAME)
      }

      ## Get which names aren't cutoff
      checkedLabels <- apply(data.frame("label" = label, "labelLoc" = transcriptNames$labelLoc), 1, cutoffLabel, fontsize = bb_transcriptsInternal$fontsize,
                             xscale = xscale,
                             vp = vp, unit = unit)
      transcriptNames$label <- checkedLabels
      transcriptNames <- transcriptNames[!is.na(transcriptNames$label), ]

      if (nrow(transcriptNames) > 0){

        transcript_names <- textGrob(label = transcriptNames$label,
                                    x = transcriptNames$labelLoc,
                                    y = transcriptNames$y + boxHeight + textHeight*0.25,
                                    just = "bottom",
                                    gp = gpar(col = transcriptNames$color, fontsize = bb_transcriptsInternal$fontsize),
                                    default.units = "native",
                                    check.overlap = TRUE)

        assign("transcript_grobs", addGrob(get("transcript_grobs", envir = bbEnv), child = transcript_names), envir = bbEnv)

      }

    }

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_transcriptsInternal$draw == TRUE){

    grid.draw(get("transcript_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  bb_transcripts$grobs <- get("transcript_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_transcripts[", vp$name, "]"))
  invisible(bb_transcripts)

}
