#' plots gene transcripts in a pileup style
#'
#' @param assembly desired genome assembly
#' @param chrom chromsome of region to be plotted
#' @param chromstart start position
#' @param chromend end position
#' @param boxHeight total height of transcripts, as a numeric value with default units or a unit value
#' @param spaceHeight height of spacing between transcripts, as a fraction of boxHeight
#' @param spaceWidth width of minimum spacing between transcripts, as a fraction of the plot's genomic range
#' @param fillcolor single value or vector specifying colors of transcripts
#' @param colorbyStrand a logical value indicating whether to color strands by the first two colors in a fillcolor vector
#' @param labels how labels should be determined; options are NULL for no labels, "transcript" for transcript names, "gene" for gene names,
#' or "both" for combined transcript and gene names
#' @param fontsize the size of gene label text (in points)
#' @param strandSplit logical indicating whether plus and minus-stranded elements should be separated
#' @param stroke numerical value indicating the stroke width for transcript body outlines
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param just string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @return Function will plot a pileup of transcripts and return a bb_transcripts object

#' @export

bb_plotTranscripts <- function(assembly = "hg19", chrom, chromstart = NULL, chromend = NULL, boxHeight = unit(2, "mm"), spaceHeight = 0.3,
                               spaceWidth = 0.02, fillcolor = c("#8a8aff", "#ff7e7e"), colorbyStrand = TRUE, labels = "transcript",
                               fontsize = 8, strandSplit = FALSE, stroke = 0.1, x = NULL, y = NULL, width = NULL, height = NULL,
                               just = c("left", "top"), default.units = "inches", draw = TRUE){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that checks errors for bb_plotTranscripts
  errorcheck_bbTranscripts <- function(transcript_plot, labels){

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

  parse_starts <- function(range){
    ## Separate character range into two numeric coords
    range <- strsplit(range, "-")[[1]]
    start <- as.numeric(range[1])

    return(start)
  }

  parse_widths <- function(range){
    ## Separate character range into two numeric coords
    range <- strsplit(range, "-")[[1]]
    start <- as.numeric(range[1])
    end <- as.numeric(range[2])
    width <- end - start

    return(width)
  }

  ## Define a function that makes exon grobs
  exon_grobs <- function(df, boxHeight){

    exon_ranges <- as.list(strsplit(as.character(df[5]), ",")[[1]])

    if (length(exon_ranges) > 0){

      starts <- lapply(exon_ranges, parse_starts)
      widths <- lapply(exon_ranges, parse_widths)
      exons_dataframe <- cbind(unlist(starts), unlist(widths))
      exons <- rectGrob(x = exons_dataframe[,1],
                        y = as.numeric(df[10]),
                        width = exons_dataframe[,2],
                        height = boxHeight,
                        just = c("left", "bottom"),
                        gp = gpar(fill = df[8],
                                  col = df[8],
                                  lwd = 1.25, alpha = 0.5),
                        default.units = "native")


      assign("transcript_grobs", addGrob(get("transcript_grobs", envir = bbEnv), child = exons), envir = bbEnv)

    }

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
                       gp = gpar(fill = df[8], col = NA, alpha = 0.5),
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
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  transcript_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, width = width,
                                height = height, x = x, y = y, justification = just, grobs = NULL, assembly = assembly), class = "bb_transcripts")
  attr(x = transcript_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = transcript_plot)
  errorcheck_bbTranscripts(transcript_plot = transcript_plot, labels = labels)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  transcript_plot <- defaultUnits(object = transcript_plot, default.units = default.units)

  # ======================================================================================================================================================================================
  # GET APPROPRIATE BUILD DATA
  # ======================================================================================================================================================================================

  if (assembly == "hg19"){

    data <- bb_hg19transcripts
    genome <- bb_hg19

  }

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  if (is.null(chromstart) & is.null(chromend)){

    ## Just chromosome
    data <- data[which(data$Chromosome == chrom),]
    transcript_plot$chromstart <- 1
    transcript_plot$chromend <- genome[which(genome$chrom == chrom),]$length

  } else {

    ## Chromosome and any overlapping regions of chromstart/chromend
    data <- data[which(data$Chromosome == chrom & data$Start <= chromend & data$Stop >= chromstart),]

  }

  ## Get width of each transcript
  data$width <- data$Stop - data$Start

  # ======================================================================================================================================================================================
  # COLORS
  # ======================================================================================================================================================================================

  if (colorbyStrand == TRUE){
    if (length(fillcolor) == 1){
      posCol <- fillcolor
      negCol <- fillcolor
    } else {
      posCol <- fillcolor[1]
      negCol <- fillcolor[2]
    }

    pos <- data[which(data$Strand == "+"),]
    pos$color <- posCol
    neg <- data[which(data$Strand == "-"),]
    neg$color <- negCol
    data <- rbind(pos, neg)

  } else {

    data$color <- rep(fillcolor[1], nrow(data))
  }


  # ======================================================================================================================================================================================
  # SEPARATE DATA INTO STRANDS
  # ======================================================================================================================================================================================

  if (strandSplit == TRUE){

    ## assuming strand is in the 6th column
    posStrand <- data[which(data$Strand == "+"),]
    minStrand <- data[which(data$Strand == "-"),]

  }


  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_transcripts", length(grep(pattern = "bb_transcripts", x = currentViewports)) + 1)

  if (is.null(x) & is.null(y)){

    height <- 0.5
    yscale <- strand_scale(strandSplit = strandSplit, height = height)

    vp <- viewport(height = unit(0.5, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = c(transcript_plot$chromstart, transcript_plot$chromend),
                   yscale = yscale,
                   just = "center",
                   name = vp_name)


    if (draw == TRUE){

      vp$name <- "bb_transcripts1"
      grid.newpage()
    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = transcript_plot)

    height <- convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = TRUE)
    yscale <- strand_scale(strandSplit = strandSplit, height = height)

    ## Make viewport for gene track
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(transcript_plot$chromstart, transcript_plot$chromend),
                   yscale = yscale,
                   just = just,
                   name = vp_name)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("transcript_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # DETERMINE ROWS FOR EACH ELEMENT
  # ======================================================================================================================================================================================

  ## Determine how many rows are going to fit based on boxHeight, spaceHeight, and fontsize
  if (is.null(labels)){
    textHeight <- unit(0, "npc")
  } else {
    textHeight <- heightDetails(textGrob(label = "A", gp = gpar(fontsize = fontsize)))
  }


  if (is.null(x) & is.null(y)){

    pushViewport(vp)
    boxHeight <- convertHeight(boxHeight, unitTo = "npc", valueOnly = T)
    spaceHeight <- boxHeight*spaceHeight
    textHeight <- convertHeight(textHeight, unitTo = "npc", valueOnly = T)
    upViewport()
    unit <- "npc"

  } else {

    boxHeight <- convertHeight(boxHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    spaceHeight <- boxHeight*spaceHeight
    textHeight <- convertHeight(textHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    unit <- get("page_units", envir = bbEnv)

  }

  maxRows <- floor(height/(boxHeight + spaceHeight + textHeight + 0.25*textHeight))
  wiggle <- abs(transcript_plot$chromend - transcript_plot$chromstart) * spaceWidth

  if (strandSplit == FALSE){

    if (nrow(data) > 0){

      data$row <- 0

      ## Put in order of citation number
      data <- data[order(data$Citation, decreasing = TRUE),]

      ## Convert to numeric matrix for Rcpp function parsing
      dataMatrix <- as.matrix(data[,c(5,6,13)])

      ## Assign a row for each element
      rowData <- checkRow(dataMatrix, maxRows, 2, wiggle)

      rowData <- as.data.frame(rowData)
      rowData <- cbind(rowData, data$Transcript, data$Exons, data$UTRs, data$width, data$color, data$Gene)
      colnames(rowData) <- c("Start", "Stop", "row", "Transcript", "Exons", "UTRs", "width", "color", "Gene")

      if (any(rowData$row == 0)){
        rowData <- rowData[which(rowData$row != 0),]
        warning("Not all transcripts shown.", call. = FALSE)

        limitGrob <- textGrob(label = "+", x = unit(1, "npc"), y = unit(1, "npc"),
                              just = c("right", "top"), gp = gpar(col = "black"))
        assign("transcript_grobs", addGrob(gTree = get("transcript_grobs", envir = bbEnv), child = limitGrob), envir = bbEnv)

      }

      ## Change row index to 0 fo y-coordinate setting
      rowData$row <- rowData$row - 1
      rowData$y <- rowData$row*(boxHeight + spaceHeight + textHeight + 0.25*textHeight)

      ## Reset rows for colors
      rowData$row <- rowData$row + 1



    } else {
      rowData <- data.frame()
    }

  } else {

    if (nrow(posStrand) > 0){

      ## Put in order of citation number
      posStrand <- posStrand[order(posStrand$Citation, decreasing = TRUE),]
      posStrand$row <- 0
      ## Convert to numeric matrix for Rcpp function parsing
      posMatrix <- as.matrix(posStrand[,c(5,6,13)])
      posData <- checkRow(posMatrix, floor(maxRows/2), 2, wiggle)
      posData <- as.data.frame(posData)
      posData <- cbind(posData, posStrand$Transcript, posStrand$Exons, posStrand$UTRs, posStrand$width, posStrand$color, posStrand$Gene)
      colnames(posData) <- c("Start", "Stop", "row", "Transcript", "Exons", "UTRs", "width", "color", "Gene")
      if (any(posData$row == 0)){
        posData <- posData[which(posData$row != 0),]
        warning("Not all plus strand transcripts shown.", call. = FALSE)
        limitGrob1 <- textGrob(label = "+", x = unit(1, "npc"), y = unit(1, "npc"),
                               just = c("right", "top"), gp = gpar(col = "black"))
        assign("transcript_grobs", addGrob(gTree = get("transcript_grobs", envir = bbEnv), child = limitGrob1), envir = bbEnv)
      }

      ## Set row index to 0 for y-coordinate setting
      posData$row <- posData$row - 1
      posData$y <- (0.5*spaceHeight) + posData$row*(boxHeight + spaceHeight + textHeight + 0.25*textHeight)

      ## Reset rows for colors
      posDF$row <- posDF$row + 1
      posDF$row <- posDF$row + floor(maxRows/2)

    } else {
      posData <- data.frame()
    }

    if (nrow(minStrand) > 0){

      ## Put in order of citation number
      minStrand <- minStrand[order(minStrand$Citation, decreasing = TRUE),]
      minStrand$row <- 0
      minMatrix <- as.matrix(minStrand[,c(5,6,13)])
      minData <- checkRow(minMatrix, floor(maxRows/2), 2, wiggle)
      minData <- as.data.frame(minData)
      minData <- cbind(minData, minStrand$Transcript, minStrand$Exons, minStrand$UTRs, minStrand$width, minStrand$color, minStrand$Gene)
      colnames(minData) <- c("Start", "Stop", "row", "Transcript", "Exons", "UTRs", "width", "color", "Gene")
      if (any(minData$row == 0)){
        minData <- minData[which(minData$row != 0),]
        warning("Not all minus strand transcripts shown.", call. = FALSE)
        limitGrob2 <- textGrob(label = "+", x = unit(1, "npc"), y = unit(0, "npc"),
                               just = c("right", "bottom"), gp = gpar(col = "black"))
        assign("transcript_grobs", addGrob(gTree = get("transcript_grobs", envir = bbEnv), child = limitGrob2), envir = bbEnv)

      }

      ## Set row index to 0 for y-coordinate setting
      minData$row <- minData$row - 1
      minData$y <- ((0.5*spaceHeight + boxHeight + textHeight + 0.25*textHeight) + minData$row*(boxHeight + spaceHeight + textHeight + 0.25*textHeight))*-1

      ## Reset rows for colors
      rowIndex <- minDF$row + 1
      rowRange <- floor(maxRows/2):1
      minDF$row <- rowRange[rowIndex]

    } else {
      minData <- data.frame()
    }

    rowData <- rbind(posData, minData)

  }

  # ======================================================================================================================================================================================
  # UPDATE COLORS IF NECESSARY
  # ======================================================================================================================================================================================

  if (colorbyStrand == FALSE & length(fillcolor) > 1){

    colors <- rep(fillcolor, ceiling(maxRows/length(fillcolor)))[1:maxRows]
    indeces <- rowData$row
    rowData$color <- colors[indeces]

  }


  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (nrow(rowData) > 0){


    if ((transcript_plot$chromend - transcript_plot$chromstart) >= 25000000){

      ##########################################################
      ## JUST TRANSCRIPT LINES
      ##########################################################

      transcriptLine <- rectGrob(x = rowData$Start,
                                 y = rowData$y,
                                 width = rowData$width,
                                 height = boxHeight,
                                 just = c("left", "bottom"),
                                 gp = gpar(fill = rowData$color, col = rowData$color, lwd = stroke, alpha = 0.5),
                                 default.units = "native")

    } else {

      ##########################################################
      ## TRANSCRIPT LINES, EXONS, AND UTRS
      ##########################################################

      transcriptLine <- rectGrob(x = rowData$Start,
                                 y = rowData$y + 0.5*boxHeight,
                                 width = rowData$width,
                                 height = boxHeight*0.2,
                                 just = "left",
                                 gp = gpar(fill = rowData$color, col = makeTransparent(rowData$color, alpha = 0.5), lwd = stroke, alpha = 0.5),
                                 default.units = "native")


      invisible(apply(rowData, 1, exon_grobs, boxHeight = boxHeight))
      invisible(apply(rowData, 1, utr_grobs, boxHeight = boxHeight))


    }

    assign("transcript_grobs", addGrob(get("transcript_grobs", envir = bbEnv), child = transcriptLine), envir = bbEnv)


    ##########################################################
    ## TRANSCRIPT NAME LABELS
    ##########################################################

    if (!is.null(labels)){

      ## Add column with center location of each gene label
      rowData$labelLoc <- rowMeans(rowData[c("Start", "Stop")])

      if (labels == "transcript"){
        label <- rowData$Transcript
      } else if (labels == "gene"){
        label <- rowData$Gene
      } else {
        label <- paste0(rowData$Gene, ":", rowData$Transcript)
      }

      checkedLabels <- apply(data.frame("label" = label, "labelLoc" = rowData$labelLoc), 1, cutoffLabel, fontsize = fontsize,
                             xscale = c(transcript_plot$chromstart, transcript_plot$chromend),
                             vp = vp, unit = unit)
      rowData$label <- checkedLabels
      rowData <- rowData[!is.na(rowData$label), ]

      if (nrow(rowData) > 0){

        transcriptNames <- textGrob(label = rowData$label,
                                    x = rowData$labelLoc,
                                    y = rowData$y + boxHeight + textHeight*0.25,
                                    just = "bottom",
                                    gp = gpar(col = rowData$color, fontsize = fontsize),
                                    default.units = "native",
                                    check.overlap = TRUE)

        assign("transcript_grobs", addGrob(get("transcript_grobs", envir = bbEnv), child = transcriptNames), envir = bbEnv)

      }

    }

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("transcript_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  transcript_plot$grobs <- get("transcript_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(transcript_plot)

}
