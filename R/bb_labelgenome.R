#' Adds genome coordinates to the axis of a BentoBox plot
#'
#' @param plot BentoBox plot to annotate.
#' @param x A numeric or unit object specifying x-location.
#' @param y A numeric or unit object specifying y-location.
#' @param params an optional "bb_params" object space containing relevant function parameters
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
#' @param ticks A numeric vector of x-value locations for tick marks.
#' @param tcl Length of tickmark as fraction of text height.
#' @param default.units A string indicating the default units to use if x or y are only given as numerics.
#' @export

bb_labelGenome <- function(plot, x, y, params = NULL, just = c("left", "top"), scale = "bp", fontsize = 10, fontcolor = "black",
                           linecolor = "black", commas = TRUE, sequence = TRUE, boxWidth = 0.5, ticks = NULL, tcl = 0.5, default.units = "inches", ...){

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

    if (length(object$chrom) == 1){

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


      ## Get assembly data
      txdbChecks <- check_loadedPackage(package = object$assembly$TxDb, message = paste(paste0("`", object$assembly$TxDb,"`"), "not loaded. Please install and load to label genome."))
      if (txdbChecks == TRUE){
        tx_db <- eval(parse(text = object$assembly$TxDb))
        assembly_data <- as.data.frame(setDT(as.data.frame(seqlengths(tx_db)), keep.rownames = TRUE))
        assembly_data <- assembly_data[which(assembly_data[,1] %in% object$chrom),]
        ## get the offsets based on spacer for the assembly
        offsetAssembly <- spaceChroms(assemblyData = assembly_data, space = plot$space)
        cumsums <- cumsum(as.numeric(assembly_data[,2]))
        spacer <- cumsums[length(cumsum(as.numeric(assembly_data[,2])))] * plot$space
        xscale <- c(0, max(offsetAssembly[,4]) + spacer)
      } else {
        xscale <- c(0, 1)
      }



      vp <- viewport(width = convertedPageCoords$width, height = convertedPageCoords$height,
                     x = convertedPageCoords$x, y = convertedPageCoords$y,
                     just = just,
                     name = vp_name,
                     xscale = xscale,
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
                                  gp = object$gp, default.units = "native")
        line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = y0_coord, y1 = y0_coord,
                             vp = chrom_vp,
                             gp = object$gp, default.units = "native")
        object$gp$col <- object$gp$fontcolor
        startLab <- textGrob(label = paste(startLabel, scale), x = unit(0, "npc"), y =  unit(tgH + 0.25*tgH, "native"),
                             just = c("left", "top"),
                             vp = chrom_vp,
                             gp = object$gp)
        endLab <- textGrob(label = paste(endLabel, scale), x = unit(1, "npc"), y = unit(tgH + 0.25*tgH, "native"),
                           just = c("right", "top"),
                           vp = chrom_vp,
                           gp = object$gp)
        object$gp$fontface <- "bold"
        chromLab <- textGrob(label = chromLabel, x = unit(0.5, "npc"), y = unit(tgH + 0.25*tgH, "native"),
                             vp = chrom_vp,
                             gp = object$gp, just = c("center", "top"))
        assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children =  gList(line, chromLab, startLab, endLab, tickGrobs)), envir = bbEnv)



      } else {

        line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = height, y1 = height,
                             vp = chrom_vp,
                             gp = object$gp, default.units = "native")
        object$gp$col <- object$gp$fontcolor
        startLab <- textGrob(label = paste(startLabel, scale, sep = " "), x = unit(0, "npc"), y = unit(0.85, "npc"),
                             vp = chrom_vp,
                             just = c("left", "top"),
                             gp = object$gp)
        endLab <- textGrob(label = paste(endLabel, scale, sep = " "), x = unit(1, "npc"), y = unit(0.85, "npc"),
                           vp = chrom_vp,
                           just = c("right", "top"),
                           gp = object$gp)
        object$gp$fontface <- "bold"
        chromLab <- textGrob(label = chromLabel, x = unit(0.5, "npc"), y = unit(0, "npc"),
                             vp = chrom_vp,
                             gp = object$gp,
                             just = c("center", "bottom"))

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
                                  gp = object$gp, default.units = "native")
        line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = y0_coord, y1 = y0_coord,
                             gp = object$gp, default.units = "native")
        object$gp$col <- object$gp$fontcolor
        startLab <- textGrob(label = paste(startLabel, scale), x = unit(0, "npc"), y =  unit(tgH + 0.25*tgH, "native"),
                             just = c("left", "top"),
                             gp = object$gp)
        endLab <- textGrob(label = paste(endLabel, scale), x = unit(1, "npc"), y = unit(tgH + 0.25*tgH, "native"),
                           just = c("right", "top"),
                           gp = object$gp)
        object$gp$fontface <- "bold"
        chromLab <- textGrob(label = chromLabel, x = unit(0.5, "npc"), y = unit(tgH + 0.25*tgH, "native"),
                             gp = object$gp$fontface, just = c("center", "top"))
        assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children =  gList(line, chromLab, startLab, endLab, tickGrobs)), envir = bbEnv)

      } else {

        line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                             y0 = height, y1 = height,
                             gp = object$gp, default.units = "native")
        object$gp$col <- object$gp$fontcolor
        startLab <- textGrob(label = paste(startLabel, scale, sep = " "), x = unit(0, "npc"), y = unit(0.85, "npc"),
                             just = c("left", "top"),
                             gp = object$gp)
        endLab <- textGrob(label = paste(endLabel, scale, sep = " "), x = unit(1, "npc"), y = unit(0.85, "npc"),
                           just = c("right", "top"),
                           gp = object$gp)
        object$gp$fontface <- "bold"
        chromLab <- textGrob(label = chromLabel, x = unit(0.5, "npc"), y = unit(0.85, "npc"),
                             gp = object$gp,
                             just = c("center", "top"))

        assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)

      }

    }


  }

  ## Define a function that makes line and text grobs for whole assembly labels
  genome_grobs <- function(plot, object, tgH, vp){

    ## Initialize gTree
    assign("label_grobs", gTree(vp = vp), envir = bbEnv)

    ## Get assembly data
    txdbChecks <- suppressWarnings(check_loadedPackage(package = object$assembly$TxDb, message = NULL))
    if (txdbChecks == TRUE){
      tx_db <- eval(parse(text = object$assembly$TxDb))
      assembly_data <- as.data.frame(setDT(as.data.frame(seqlengths(tx_db)), keep.rownames = TRUE))
      assembly_data <- assembly_data[which(assembly_data[,1] %in% object$chrom),]
      ## Get the offsets based on spacer for the assembly
      offsetAssembly <- spaceChroms(assemblyData = assembly_data, space = plot$space)

      ## Get the centers of each chrom
      chromCenters <- (offsetAssembly[,3] + offsetAssembly[,4]) / 2

      tgH <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

      line <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = tgH, y1 = tgH, gp = object$gp,
                           default.units = "native")
      object$gp$col <- object$gp$fontcolor
      labels <- textGrob(label = gsub("chr", "", offsetAssembly[,1]), x = chromCenters, y = unit(0, "npc"), just = c("center", "bottom"),
                         gp = object$gp,
                         default.units = "native")
      assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, labels)), envir = bbEnv)

    }


  }

  ## Define a function that makes sequence grobs (boxes or letters)
  seq_grobs <- function(object, seqHeight, seqType, assembly, chromLabel, vp, boxWidth){

    bsgenome <- eval(parse(text = object$assembly$BSgenome))
    ## Get sequence in that region
    sequence <- strsplit(as.character(BSgenome::getSeq(bsgenome,
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
                           gp = gpar(col = dfSequence$col, fontsize = object$gp$fontsize - 2))

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
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(scale)) scale <- NULL
  if(missing(fontsize)) fontsize <- NULL
  if(missing(fontcolor)) fontcolor <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(commas)) commas <- NULL
  if(missing(sequence)) sequence <- NULL
  if(missing(boxWidth)) boxWidth <- NULL
  if(missing(tcl)) tcl <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if plot/x/y arguments are missing (could be in object)
  if(!hasArg(plot)) plot <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL

  ## Compile all parameters into an internal object
  bb_glabelInternal <- structure(list(plot = plot, x = x, y = y, just = just, scale = scale, fontsize = fontsize, fontcolor = fontcolor,
                                      linecolor = linecolor, commas = commas, sequence = sequence, boxWidth = boxWidth,
                                      ticks = ticks, tcl = tcl, default.units = default.units), class = "bb_glabelInternal")
  bb_glabelInternal <- parseParams(bb_params = params, object_params = bb_glabelInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_glabelInternal$just)) bb_glabelInternal$just <- c("left", "top")
  if(is.null(bb_glabelInternal$scale)) bb_glabelInternal$scale <- "bp"
  if(is.null(bb_glabelInternal$fontsize)) bb_glabelInternal$fontsize <- 10
  if(is.null(bb_glabelInternal$fontcolor)) bb_glabelInternal$fontcolor <- "black"
  if(is.null(bb_glabelInternal$linecolor)) bb_glabelInternal$linecolor <- "black"
  if(is.null(bb_glabelInternal$commas)) bb_glabelInternal$commas <- TRUE
  if(is.null(bb_glabelInternal$sequence)) bb_glabelInternal$sequence <- TRUE
  if(is.null(bb_glabelInternal$boxWidth)) bb_glabelInternal$boxWidth <- 0.5
  if(is.null(bb_glabelInternal$tcl)) bb_glabelInternal$tcl <- 0.5
  if(is.null(bb_glabelInternal$default.units)) bb_glabelInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_genome_label <- structure(list(assembly = plot$assembly, chrom = plot$chrom, chromstart = plot$chromstart, chromend = plot$chromend, x = bb_glabelInternal$x, y = bb_glabelInternal$y, width = NULL, height = NULL,
                                    just = bb_glabelInternal$just, scale = bb_glabelInternal$scale, grobs = NULL,
                                    gp = gpar(fontsize = bb_glabelInternal$fontsize, col = bb_glabelInternal$linecolor, fontcolor = bb_glabelInternal$fontcolor, ...)), class = "bb_genomeLabel")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_glabelInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_genome_label$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_genome_label$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)

  check_bbpage(error = "Cannot add a genome label without a BentoBox page.")
  errorcheck_bb_labelgenome(plot = bb_glabelInternal$plot, scale = bb_genome_label$scale, ticks = bb_glabelInternal$ticks, object = bb_genome_label)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  bb_genome_label$assembly <- parse_bbAssembly(assembly = bb_genome_label$assembly)

  # ======================================================================================================================================================================================
  # SET UP PAGE/SCALE
  # ======================================================================================================================================================================================

  ## Determine scale of labels
  if (bb_genome_label$scale == "bp"){
    fact = 1
    format = "d"
  }
  if (bb_genome_label$scale == "Mb"){
    fact = 1000000
    format = NULL
    warning("Chromosome start and stop will be rounded.", call. = FALSE)
  }
  if (bb_genome_label$scale == "Kb"){
    fact = 1000
    format = "d"
    warning("Chromosome start and stop will be rounded.", call. = FALSE)
  }

  tgH <- convertHeight(heightDetails(textGrob(label = bb_genome_label$scale, x = 0.5, y = 0.5, default.units = "npc", gp = bb_genome_label$gp)),
                       unitTo = get("page_units", envir = bbEnv))
  seq_height <- heightDetails(textGrob(label = "A", x = 0.5, y = 0.5, default.units = "npc", gp = gpar(fontsize = bb_genome_label$gp$fontsize - 2)))
  seq_height <- convertHeight(seq_height + 0.05*seq_height, unitTo = get("page_units", envir = bbEnv))

  # ======================================================================================================================================================================================
  # SET PARAMETERS
  # ======================================================================================================================================================================================
  ## If single chrom/chromstart/chromend label - comma parsing
  if (length(bb_genome_label$chrom) == 1){

    commaLabels <- comma_labels(object = bb_genome_label, commas = bb_glabelInternal$commas, format = format, fact = fact)
    chromstartlabel <- commaLabels[[1]]
    chromendlabel <- commaLabels[[2]]

  }

  ## Total label dimensions, taking into account tick and sequence height
  bb_genome_label$width <- bb_glabelInternal$plot$width
  if (!is.null(bb_glabelInternal$ticks)){
    tick_height <- tgH*(bb_glabelInternal$tcl)
    height <- convertHeight(tgH + tick_height + 0.5*tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

  } else {
    tick_height <- NULL
    height <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
  }

  if (bb_glabelInternal$sequence == TRUE){

    seq_height <- convertHeight(seq_height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    bb_genome_label$height <- unit(height + seq_height, get("page_units", envir = bbEnv))

  } else {

    bb_genome_label$height <- unit(height, get("page_units", envir = bbEnv))
  }

  ## Determine appropriate scaling of nucleotides and check for BSgenome packages
  if (length(bb_genome_label$chrom) == 1){

    seqType <- NULL

    if (!is.null(bb_genome_label$assembly$BSgenome)){

      bsChecks <- check_loadedPackage(package = bb_genome_label$assembly$BSgenome, message = paste(paste0("`", bb_genome_label$assembly$BSgenome,"`"), "not loaded. Sequence information will not be displayed."))
      if (bsChecks == TRUE){

        labelWidth <- convertWidth(bb_genome_label$width, unitTo = "inches", valueOnly = T)
        bpWidth <- convertWidth(widthDetails(textGrob(label = "A", x = 0.5, y = 0.5, default.units = "npc",
                                                      gp = gpar(fontsize = bb_genome_label$gp$fontsize - 2))),
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

    } else {

      warning("No `BSgenome` package found for the input assembly. Sequence information cannot be displayed.", call. = FALSE)

    }


  }

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================
  if (!"unit" %in% class(bb_genome_label$x)){

    if (!is.numeric(bb_genome_label$x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_glabelInternal$default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_genome_label$x <- unit(bb_genome_label$x, bb_glabelInternal$default.units)

  }

  if (!"unit" %in% class(bb_genome_label$y)){

    if (!is.numeric(bb_genome_label$y)){

      stop("y-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_glabelInternal$default.units)){

      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_genome_label$y <- unit(bb_genome_label$y, bb_glabelInternal$default.units)

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================
  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_genomeLabel", length(grep(pattern = "bb_genomeLabel", x = currentViewports)) + 1)

  ## Make viewport
  vp <- parse_viewport(plot = bb_glabelInternal$plot, object = bb_genome_label, height = height, sequence = bb_glabelInternal$sequence, seqType = seqType, seqHeight = seq_height,
                       vp_name = vp_name, just = bb_genome_label$just)

  # ======================================================================================================================================================================================
  # GROBS AND GTREE
  # ======================================================================================================================================================================================

  ## Chrom/chromstart/chromend grobs
  if (length(bb_genome_label$chrom) == 1){

    chrom_grobs(tgH = tgH, ticks = bb_glabelInternal$ticks, tickHeight = tick_height, sequence = bb_glabelInternal$sequence, seqType = seqType, scale = bb_genome_label$scale,
                chromLabel = bb_genome_label$chrom,
                startLabel = chromstartlabel, endLabel = chromendlabel, height = height, object = bb_genome_label, vp = vp)

    ## Sequence grobs
    if (bb_glabelInternal$sequence == TRUE & !is.null(seqType)){

      seq_grobs(object = bb_genome_label, seqHeight = seq_height, seqType = seqType, assembly = bb_genome_label$assembly, chromLabel = bb_genome_label$chrom, vp = vp,
                boxWidth = bb_glabelInternal$boxWidth)

    }

  } else {
    ## Whole genome grobs
    genome_grobs(plot = bb_glabelInternal$plot, object = bb_genome_label, tgH = tgH, vp = vp)

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
