#' Adds genome coordinates to the axis of a BentoBox plot
#'
#' @param plot BentoBox plot to annotate.
#' @param x A numeric or unit object specifying x-location.
#' @param y A numeric or unit object specifying y-location.
#' @param just The justification of the label relative to its (x, y) location.If there are two values, the first specifies horizontal justification and
#' the second value specifies vertical justification. Possible string values are: "left", "right", "centre", "center", "bottom", and "top".
#' @param scale scale of the genome of the label. Options are "bp", "Kb", or "Mb"
#' @param fontsize Fontsize for labels (in points).
#' @param fontcolor Color for all text and lines (with the exception of sequence information).
#' @param linecolor Axis linecolor.
#' @param lwd Axis linewidth.
#' @param fontfamily Text fontfamily.
#' @param commas A logical value indicating whether to include commas in start and stop labels.
#' @param sequence A logical value indicating whether to include sequence information above the label axis (only at appropriate resolutions).
#' @param assembly A character value indicating the genome assembly of the label. Possible values are: "hg19".
#' @param ticks A numeric vector of x-value locations for tick marks.
#' @param tcl Length of tickmark as fraction of text height.
#' @param default.units A string indicating the default units to use if x or y are only given as numerics.
#' @export

bb_labelGenome <- function(plot, x, y, just = c("left", "top"),
                           scale = "bp", fontsize = 10, fontcolor = "black", linecolor = "black", lwd = 1,
                           fontfamily = "", commas = TRUE, sequence = TRUE, boxWidth = 0.5, assembly = "hg19", ticks = NULL, tcl = 0.5, default.units = "inches"){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_labelgenome
  errorcheck_bb_labelgenome <- function(plot, scale, ticks, object){

    ## Check that input plot is a valid type of plot to be annotated
    if (class(plot) != "bb_manhattan"){

      ## Manhattan plots can do whole genome assembly but other plots can't
      inputNames <- attributes(plot)$names
      if (!("chrom" %in% inputNames) | !("chromstart" %in% inputNames) | !("chromend" %in% inputNames)){

        stop("Invalid input plot. Please input a plot that has genomic coordinates associated with it.")

      }

    }

    ## Check that scale is an appropriate value
    if (!scale %in% c("bp", "Kb", "Mb")){

      stop("Invalid \'scale\'. Options are \'bp\', \'Kb\', or \'Mb\'.", call. = FALSE)

    }

    if (!is.null(ticks)){

      ## Can't have ticks if label is genome assembly for manhattan plot
      if (is.null(object$chrom)){

        stop("Cannot add tick marks to a genome label of entire genome assembly.", call. = FALSE)

      }

      ## Make sure ticks fall within the chromstart to chromend range
      if (range(ticks)[1] < object$chromstart | range(ticks)[2] > object$chromend){

        stop("Given tick locations do not fall within the genomic range.", call. = FALSE)

      }

    }

  }

  ## Define a function that parses genome assembly vs. chrom/chromstart/chromend label
  parse_plot <- function(plot, object){

    ## Manhattan plots will require different types of genome labels
    if (class(plot) == "bb_manhattan"){

      ## whole genome option
      if (is.null(plot$chrom)){

        object$assembly <- plot$assembly

      } else {

        chromNumber <- gsub(pattern = "chr", replacement = "", x = plot$chrom)
        chrom <- paste0("chr", chromNumber)
        object$chrom <- chrom
        object$chromstart <- plot$chromstart
        object$chromend <- plot$chromend

      }

    } else {

      chromNumber <- gsub(pattern = "chr", replacement = "", x = plot$chrom)
      chrom <- paste0("chr", chromNumber)
      object$chrom <- chrom
      object$chromstart <- plot$chromstart
      object$chromend <- plot$chromend

    }

    return(object)

  }

  ## Define a function that adds commas to chromstart/chromend labels
  comma_labels <- function(object, commas, format, fact){

    if (commas == TRUE){

      chromstartlabel <- formatC(round(object$chromstart/fact, 1), format = format, big.mark = ",")
      chromendlabel <- formatC(round(object$chromend/fact, 1), format = format, big.mark = ",")

    } else {

      chromstartlabel <- round(object$chromstart/fact, 1)
      chromendlabel <- round(object$chromend/fact, 1)

    }

    return(list(chromstartlabel, chromendlabel))

  }

  ## Define a function that parses the viewport for genome assembly vs. chrom/chromstart/chromend label w/ or w/o sequence viewport
  parse_viewport <- function(plot, object, height, sequence, seqType, seqHeight, vp_name, just){

    convertedPageCoords <- convert_page(object = object)
    convertedViewport <- viewport(width = convertedPageCoords$width, height = convertedPageCoords$height,
                                  x = convertedPageCoords$x, y = convertedPageCoords$y, just = just)

    ## Get x and y coordinates of top left of what would be the entire viewport
    topLeftViewport <- vp_topLeft(viewport = convertedViewport)

    seq_height <- unit(seqHeight, get("page_units", envir = bbEnv))

    if (!is.null(object$chrom)){

      if (sequence == TRUE & !is.null(seqType)){

        ## One vp for genome
        vp1 <- viewport(width = convertedPageCoords$width, height = unit(height, get("page_units", envir = bbEnv)),
                        x = topLeftViewport[[1]], y = topLeftViewport[[2]] - seq_height,
                        just = c("left", "top"),
                        name = paste0(vp_name, "_01"),
                        xscale = c(object$chromstart, object$chromend),
                        yscale = c(0, height))
        ## One vp for sequence
        vp2 <- viewport(width = convertedPageCoords$width, height = seq_height,
                        x = topLeftViewport[[1]], y = topLeftViewport[[2]],
                        just = c("left", "top"),
                        name = paste0(vp_name, "_02"),
                        clip = "on",
                        xscale = c(object$chromstart, object$chromend))

        ## Combine viewports into one
        vp <- vpList(vp1, vp2)

      } else {

        vp <- viewport(width = convertedPageCoords$width, height = convertedPageCoords$height,
                        x = convertedPageCoords$x, y = convertedPageCoords$y,
                        just = just,
                        name = vp_name,
                        xscale = c(object$chromstart, object$chromend),
                        yscale = c(0, height))
      }



    } else {

      ## get the offsets based on spacer for the assembly
      if (object$assembly == "hg19"){
        assembly_data <- bb_hg19
      }

      ## get the offsets based on spacer for the assembly
      offsetAssembly <- parse_assembly(assemblyData = assembly_data, space = plot$space)
      cumsums <- cumsum(as.numeric(assembly_data[,2]))
      spacer <- cumsums[length(cumsum(as.numeric(assembly_data[,2])))] * plot$space

      vp <- viewport(width = convertedPageCoords$width, height = convertedPageCoords$height,
                     x = convertedPageCoords$x, y = convertedPageCoords$y,
                     just = just,
                     name = vp_name,
                     xscale = c(0, max(offsetAssembly[,4]) + spacer),
                     yscale = c(0, height))
    }

    return(vp)
  }

  ## Define a function that makes tick, line, and text grobs for chrom/chromstart/chromend labels
  chrom_grobs <- function(tgH, ticks, tickHeight, sequence, seqType, scale, chromLabel, startLabel, endLabel, height, object, vp){

    if (sequence == TRUE & !is.null(seqType)){
      assign("label_grobs", gTree(), envir = bbEnv)
      chrom_vp <- vp[[1]]

      if (!is.null(ticks)){
        tgH <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
        tick_height <- convertHeight(tickHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
        x_coords <- ticks
        y0_coord <- tgH + tick_height + 0.5*tgH
        y1_coords <- rep(tgH + 0.5*tgH, length(ticks))

        tickGrobs <- segmentsGrob(x0 = x_coords, y0 = rep(y0_coord, length(ticks)),
                                  x1 = x_coords, y1 = y1_coords,
                                  vp = chrom_vp,
                                  gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd), default.units = "native")
        line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = y0_coord, y1 = y0_coord,
                             vp = chrom_vp,
                             gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd), default.units = "native")
        chromLab <- textGrob(label = chromLabel, x = unit(0.5, "npc"), y = unit(tgH + 0.25*tgH, "native"),
                             vp = chrom_vp,
                             gp = gpar(fontface = "bold", fontsize = object$gpar$fontsize, col = object$gpar$fontcolor,
                                       fontfamily = object$gpar$fontfamily), just = c("center", "top"))
        startLab <- textGrob(label = paste(startLabel, scale), x = unit(0, "npc"), y =  unit(tgH + 0.25*tgH, "native"),
                             just = c("left", "top"),
                             vp = chrom_vp,
                             gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))
        endLab <- textGrob(label = paste(endLabel, scale), x = unit(1, "npc"), y = unit(tgH + 0.25*tgH, "native"),
                           just = c("right", "top"),
                           vp = chrom_vp,
                           gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))
        assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children =  gList(line, chromLab, startLab, endLab, tickGrobs)), envir = bbEnv)



      } else {

        line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = height, y1 = height,
                             vp = chrom_vp,
                             gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd), default.units = "native")
        chromLab <- textGrob(label = chromLabel, x = unit(0.5, "npc"), y = unit(0, "npc"),
                             vp = chrom_vp,
                             gp = gpar(fontface = "bold", fontsize = object$gpar$fontsize, col = object$gpar$fontcolor, fontfamily = object$gpar$fontfamily),
                             just = c("center", "bottom"))
        startLab <- textGrob(label = paste(startLabel, scale, sep = " "), x = unit(0, "npc"), y = unit(0.85, "npc"),
                             vp = chrom_vp,
                             just = c("left", "top"),
                             gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))
        endLab <- textGrob(label = paste(endLabel, scale, sep = " "), x = unit(1, "npc"), y = unit(0.85, "npc"),
                           vp = chrom_vp,
                           just = c("right", "top"),
                           gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))

        assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)


      }

    } else {
      assign("label_grobs", gTree(vp = vp), envir = bbEnv)

      if (!is.null(ticks)){
        tgH <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
        tick_height <- convertHeight(tickHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
        x_coords <- ticks
        y0_coord <- tgH + tick_height + 0.5*tgH
        y1_coords <- rep(tgH + 0.5*tgH, length(ticks))

        tickGrobs <- segmentsGrob(x0 = x_coords, y0 = rep(y0_coord, length(ticks)),
                                  x1 = x_coords, y1 = y1_coords,
                                  gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd), default.units = "native")
        line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = y0_coord, y1 = y0_coord,
                             gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd), default.units = "native")
        chromLab <- textGrob(label = chromLabel, x = unit(0.5, "npc"), y = unit(tgH + 0.25*tgH, "native"),
                             gp = gpar(fontface = "bold", fontsize = object$gpar$fontsize, col = object$gpar$fontcolor,
                                       fontfamily = object$gpar$fontfamily), just = c("center", "top"))
        startLab <- textGrob(label = paste(startLabel, scale), x = unit(0, "npc"), y =  unit(tgH + 0.25*tgH, "native"),
                             just = c("left", "top"),
                             gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))
        endLab <- textGrob(label = paste(endLabel, scale), x = unit(1, "npc"), y = unit(tgH + 0.25*tgH, "native"),
                           just = c("right", "top"),
                           gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))
        assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children =  gList(line, chromLab, startLab, endLab, tickGrobs)), envir = bbEnv)

      } else {

        line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = height, y1 = height,
                             gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd), default.units = "native")
        chromLab <- textGrob(label = chromLabel, x = unit(0.5, "npc"), y = unit(0.85, "npc"),
                             gp = gpar(fontface = "bold", fontsize = object$gpar$fontsize, col = object$gpar$fontcolor, fontfamily = object$gpar$fontfamily),
                             just = c("center", "top"))
        startLab <- textGrob(label = paste(startLabel, scale, sep = " "), x = unit(0, "npc"), y = unit(0.85, "npc"),
                             just = c("left", "top"),
                             gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))
        endLab <- textGrob(label = paste(endLabel, scale, sep = " "), x = unit(1, "npc"), y = unit(0.85, "npc"),
                           just = c("right", "top"),
                           gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))

        assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)

      }

    }


  }

  ## Define a function that makes line and text grobs for whole assembly labels
  genome_grobs <- function(plot, object, tgH, vp){

    ## Get internal assembly data
    if (object$assembly == "hg19"){
      assembly_data <- bb_hg19
    }

    assign("label_grobs", gTree(vp = vp), envir = bbEnv)

    ## Get the offsets based on spacer for the assembly
    offsetAssembly <- parse_assembly(assemblyData = assembly_data, space = plot$space)

    ## Get the centers of each chrom
    chromCenters <- (offsetAssembly[,3] + offsetAssembly[,4]) / 2

    tgH <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

    line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = tgH, y1 = tgH, gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd),
                         default.units = "native")
    labels <- textGrob(label = gsub("chr", "", offsetAssembly[,1]), x = chromCenters, y = unit(0, "npc"), just = c("center", "bottom"),
                       gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$fontcolor, fontfamily = object$gpar$fontfamily),
                       default.units = "native")
    assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, labels)), envir = bbEnv)

  }

  ## Define a function that makes sequence grobs (boxes or letters)
  seq_grobs <- function(object, seqHeight, seqType, assembly, chromLabel, vp, boxWidth){

    if (assembly == "hg19"){
      seqAssembly <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    }

    ## Get sequence in that region
    sequence <- strsplit(as.character(BSgenome::getSeq(seqAssembly,
                                             GenomicRanges::GRanges(seqnames = chromLabel, ranges = IRanges::IRanges(start = object$chromstart, end = object$chromend)))),
                         split = "")
    ## Make dataframe of sequence letter, position, and color
    dfSequence <- data.frame("nucleotide" = unlist(sequence), "pos" = seq(object$chromstart, object$chromend), "col" = "black")

    ## Make colors A = green, T = red, G = orange, C = blue
    dfSequence[which(dfSequence$nucleotide == "A"),]$col <- "#009600"
    dfSequence[which(dfSequence$nucleotide == "T"),]$col <- "#ff0000"
    dfSequence[which(dfSequence$nucleotide == "G"),]$col <- "#d17105"
    dfSequence[which(dfSequence$nucleotide == "C"),]$col <- "#0000ff"

    seq_vp <- vp[[2]]

    ## Make grobs based on seqType
    if (seqType == "letters"){
      seqGrobs <- textGrob(label = dfSequence$nucleotide, x = dfSequence$pos, y = unit(0.5, "npc"), just = "center",
                           vp = seq_vp,
                           default.units = "native",
                           gp = gpar(col = dfSequence$col, fontsize = object$gpar$fontsize - 2, fontfamily = object$gpar$fontfamily))

    } else if (seqType == "boxes"){
      #seq_height <- convertHeight(seqHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
      seqGrobs <- rectGrob(x = dfSequence$pos, y = unit(1, "npc"), width = boxWidth, height = unit(seqHeight - 0.05*seqHeight, get("page_units", envir = bbEnv)),
                           just = c("center", "top"),
                           vp = seq_vp,
                           default.units = "native",
                           gp = gpar(col = NA, fill = dfSequence$col))

    }

    assign("label_grobs", addGrob(gTree = get("label_grobs", envir = bbEnv), child = seqGrobs), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_genome_label <- structure(list(assembly = NULL, chrom = NULL, chromstart = NULL, chromend = NULL, x = x, y = y, width = NULL, height = NULL, scale = scale, grobs = NULL,
                                    gpar = list(fontsize = fontsize, fontcolor = fontcolor, linecolor = linecolor, lwd = lwd, fontfamily = fontfamily)), class = "genome_label")

  # ======================================================================================================================================================================================
  # PARSE TYPE OF INPUT PLOT
  # ======================================================================================================================================================================================

  ## Determine whole genome assembly labeling vs chrom/chromstart/chromend
  bb_genome_label <- parse_plot(plot = plot, object = bb_genome_label)

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot add a genome label without a BentoBox page.")
  errorcheck_bb_labelgenome(plot = plot, scale = scale, ticks = ticks, object = bb_genome_label)

  # ======================================================================================================================================================================================
  # SET UP PAGE/SCALE
  # ======================================================================================================================================================================================

  ## Determine scale of labels
  if (scale == "bp"){
    fact = 1
    format = "d"
  }
  if (scale == "Mb"){
    fact = 1000000
    format = NULL
    warning("Chromosome start and stop will be rounded.", call. = FALSE)
  }
  if (scale == "Kb"){
    fact = 1000
    format = "d"
    warning("Chromosome start and stop will be rounded.", call. = FALSE)
  }

  tgH <- convertHeight(heightDetails(textGrob(label = scale, x = 0.5, y = 0.5, default.units = "npc", gp = gpar(fontsize = fontsize, fontfamily = fontfamily))),
                       unitTo = get("page_units", envir = bbEnv))
  seq_height <- heightDetails(textGrob(label = "A", x = 0.5, y = 0.5, default.units = "npc", gp = gpar(fontsize = fontsize - 2, fontfamily = fontfamily)))
  seq_height <- convertHeight(seq_height + 0.05*seq_height, unitTo = get("page_units", envir = bbEnv))

  # ======================================================================================================================================================================================
  # SET PARAMETERS
  # ======================================================================================================================================================================================

  ## If chrom/chromstart/chromend label - comma parsing
  if (!is.null(bb_genome_label$chrom)){

    commaLabels <- comma_labels(object = bb_genome_label, commas = commas, format = format, fact = fact)
    chromstartlabel <- commaLabels[[1]]
    chromendlabel <- commaLabels[[2]]

  }

  ## Total label dimensions, taking into account tick and sequence height
  bb_genome_label$width <- plot$width
  if (!is.null(ticks)){
    tick_height <- tgH*tcl
    height <- convertHeight(tgH + tick_height + 0.5*tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

  } else {
    tick_height <- NULL
    height <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
  }

  if (sequence == TRUE){

    seq_height <- convertHeight(seq_height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    bb_genome_label$height <- unit(height + seq_height, get("page_units", envir = bbEnv))

  } else {

    bb_genome_label$height <- unit(height, get("page_units", envir = bbEnv))
  }


  ## Determine appropriate scaling of nucleotides
  if (!is.null(bb_genome_label$chrom)){

    labelWidth <- convertWidth(bb_genome_label$width, unitTo = "inches", valueOnly = T)
    bpWidth <- convertWidth(widthDetails(textGrob(label = "A", x = 0.5, y = 0.5, default.units = "npc",
                                                   gp = gpar(fontsize = fontsize - 2, fonfamily = fontfamily))),
                             unitTo = "inches", valueOnly = T)
    seqRange <- bb_genome_label$chromend - bb_genome_label$chromstart
    seqWidth <- bpWidth*seqRange

    if (seqWidth <= labelWidth){
      seqType <- "letters"
    } else if (seqWidth/labelWidth <= 9){
      seqType <- "boxes"
    } else {
      seqType <- NULL
    }
  }

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  if (class(x) != "unit"){

    if (!is.numeric(x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_genome_label$x <- unit(x, default.units)

  }

  if (class(y) != "unit"){

    if (!is.numeric(y)){

      stop("y-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(default.units)){

      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_genome_label$y <- unit(y, default.units)

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_genomeLabel", length(grep(pattern = "bb_genomeLabel", x = currentViewports)) + 1)

  ## Make viewport
  vp <- parse_viewport(plot = plot, object = bb_genome_label, height = height, sequence = sequence, seqType = seqType, seqHeight = seq_height,
                       vp_name = vp_name, just = just)

  # ======================================================================================================================================================================================
  # GROBS AND GTREE
  # ======================================================================================================================================================================================

  ## Chrom/chromstart/chromend grobs
  if (!is.null(bb_genome_label$chrom)){

    chrom_grobs(tgH = tgH, ticks = ticks, tickHeight = tick_height, sequence = sequence, seqType = seqType, scale = scale, chromLabel = bb_genome_label$chrom,
                startLabel = chromstartlabel, endLabel = chromendlabel, height = height, object = bb_genome_label, vp = vp)

    ## Sequence grobs
    if (sequence == TRUE & !is.null(seqType)){

      seq_grobs(object = bb_genome_label, seqHeight = seq_height, seqType = seqType, assembly = assembly, chromLabel = bb_genome_label$chrom, vp = vp, boxWidth = boxWidth)

    }

  } else {
    ## Whole genome grobs
    genome_grobs(plot = plot, object = bb_genome_label, vp = vp)

  }



  # ======================================================================================================================================================================================
  # ASSIGN GROBS TO SCALE OBJECT
  # ======================================================================================================================================================================================

  bb_genome_label$grobs <- get("label_grobs", envir = bbEnv)
  grid.draw(bb_genome_label$grobs)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_genome_label)
}
