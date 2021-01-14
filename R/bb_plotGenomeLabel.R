#' Plots genomic coordinates for the x or y axis of a BentoBox plot
#'
#' @param x A numeric or unit object specifying x-location.
#' @param y A numeric or unit object specifying y-location.
#' @param plot BentoBox plot to add genome label
#' @param length if not specifying plot, length of genome label axis, as a numeric or unit object
#' @param chrom if not specifying plot, chromosome of genome label
#' @param chromstart if not specifying plot, chromstart of genome label
#' @param chromend if not specifying plot, chromend of genome label
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param just The justification of the label relative to its (x, y) location.If there are two values, the first specifies horizontal justification and
#' the second value specifies vertical justification. Possible string values are: "left", "right", "centre", "center", "bottom", and "top".
#' @param scale scale of the genome of the label. Options are "bp", "Kb", or "Mb"
#' @param assembly default genome assembly as a string or a bb_assembly object (will be inherited from plot object if provided)
#' @param fontsize Fontsize for labels (in points).
#' @param fontcolor Color for all text and lines (with the exception of sequence information).
#' @param linecolor Axis linecolor.
#' @param commas A logical value indicating whether to include commas in start and stop labels.
#' @param sequence A logical value indicating whether to include sequence information above the label of an x-axis (only at appropriate resolutions).
#' @param axis "x" (x-axis) or "y" (y-axis). Sequence information will not be displayed along a y-axis.
#' @param ticks A numeric vector of x-value locations for tick marks.
#' @param tcl Length of tickmark as fraction of text height.
#' @param default.units A string indicating the default units to use if x or y are only given as numerics.
#' @export

bb_plotGenomeLabel <- function(x, y, plot = NULL, length = NULL, chrom = NULL, chromstart = NULL, chromend = NULL, params = NULL, just = c("left", "top"),
                               scale = "bp", assembly = "hg19", fontsize = 10, fontcolor = "black", linecolor = "black", commas = TRUE, sequence = TRUE, axis = "x",
                               boxWidth = 0.5, ticks = NULL, tcl = 0.5, default.units = "inches", ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_labelgenome
  errorcheck_bb_genomeLabel <- function(scale, ticks, object, axis){

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

    if(!axis %in% c("x", "y")){
      stop("Invalid \'axis\'. Options are \'x\' or \'y\'.", call. = FALSE)
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
  parse_viewport <- function(plot, object, length, depth, seqType, seqHeight, vp_name, just, axis){

    ## No matter the orientation, convert length, depth, x, and y to page units
    convertedPageCoords <- convert_page(object = structure(list(width = length, height = unit(depth, get("page_units", envir = bbEnv)), x = object$x, y = object$y), class = "bb_genomeLabelInternal"))
    ## Add "length" and "depth" into converted dimensions for better understanding
    convertedPageCoords$length <- convertedPageCoords$width
    convertedPageCoords$depth <- convertedPageCoords$height

    ## Compile new dimensions into a new dummy viewport, where the default is along the x-axis
    convertedViewport <- viewport(width = convertedPageCoords$length, height = convertedPageCoords$depth,
                                  x = convertedPageCoords$x, y = convertedPageCoords$y, just = just)
    if (length(object$chrom) == 1){

      if (!is.null(seqType)){

        ## Get x and y coordinates of top left of what would be the entire viewport
        topLeftViewport <- vp_topLeft(viewport = convertedViewport)
        seq_height <- unit(seqHeight, get("page_units", envir = bbEnv))
        ## One vp for genome
        vp1 <- viewport(width = convertedPageCoords$width, height = unit(depth, get("page_units", envir = bbEnv)),
                        x = topLeftViewport[[1]], y = topLeftViewport[[2]] - seq_height,
                        just = c("left", "top"),
                        name = paste0(vp_name, "_01"),
                        xscale = c(object$chromstart, object$chromend),
                        yscale = c(0, depth))
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

        if (axis == "y"){


          ## Update converted viewport for y-axis
          convertedViewport <- viewport(width = convertedPageCoords$depth, height = convertedPageCoords$length,
                                        x = convertedPageCoords$x, y = convertedPageCoords$y, just = just)
          ## Get x and y coordinates of bottom right to rotate x-axis viewport
          bottomRightViewport <- vp_bottomRight(viewport = convertedViewport)

          ## Make x-axis equivalent viewport and rotate into dimensions of given y-axis viewport
          vp <- viewport(width = convertedPageCoords$length, height = convertedPageCoords$depth,
                         x = bottomRightViewport[[1]], y = bottomRightViewport[[2]] + convertedPageCoords$length,
                         just = c("left", "top"),
                         name = vp_name,
                         xscale = c(object$chromstart, object$chromend),
                         yscale = c(0, depth),
                         angle = -90)

        } else {
          vp <- viewport(width = convertedPageCoords$width, height = convertedPageCoords$height,
                         x = convertedPageCoords$x, y = convertedPageCoords$y,
                         just = just,
                         name = vp_name,
                         xscale = c(object$chromstart, object$chromend),
                         yscale = c(0, depth))
        }


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


      if (axis == "y"){
        ## Update converted viewport for y-axis
        convertedViewport <- viewport(width = convertedPageCoords$depth, height = convertedPageCoords$length,
                                      x = convertedPageCoords$x, y = convertedPageCoords$y, just = just)
        ## Get x and y coordinates of bottom right to rotate x-axis viewport
        bottomRightViewport <- vp_bottomRight(viewport = convertedViewport)
        ## Make x-axis equivalent viewport and rotate into dimensions of given y-axis viewport
        vp <- viewport(width = convertedPageCoords$length, height = convertedPageCoords$depth,
                       x = bottomRightViewport[[1]], y = bottomRightViewport[[2]] + convertedPageCoords$length,
                       just = c("left", "top"),
                       name = vp_name,
                       xscale = c(object$chromstart, object$chromend),
                       yscale = c(0, depth),
                       angle = -90)

      } else {
        vp <- viewport(width = convertedPageCoords$width, height = convertedPageCoords$height,
                       x = convertedPageCoords$x, y = convertedPageCoords$y,
                       just = just,
                       name = vp_name,
                       xscale = xscale,
                       yscale = c(0, depth))
      }


    }

    return(vp)
  }

  ## Define a function that makes tick, line, and text grobs for chrom/chromstart/chromend labels
  chrom_grobs <- function(tgH, ticks, tickHeight, seqType, scale, chromLabel, startLabel, endLabel, height, object, vp){

    if (!is.null(seqType)){
      assign("genomeLabel_grobs", gTree(), envir = bbEnv)
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
        assign("genomeLabel_grobs", setChildren(get("genomeLabel_grobs", envir = bbEnv), children =  gList(line, chromLab, startLab, endLab, tickGrobs)), envir = bbEnv)



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

        assign("genomeLabel_grobs", setChildren(get("genomeLabel_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)


      }

    } else {
      assign("genomeLabel_grobs", gTree(vp = vp), envir = bbEnv)

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
        assign("genomeLabel_grobs", setChildren(get("genomeLabel_grobs", envir = bbEnv), children =  gList(line, chromLab, startLab, endLab, tickGrobs)), envir = bbEnv)

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

        assign("genomeLabel_grobs", setChildren(get("genomeLabel_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)

      }

    }


  }

  ## Define a function that makes line and text grobs for whole assembly labels
  genome_grobs <- function(plot, object, tgH, vp){

    ## Initialize gTree
    assign("genomeLabel_grobs", gTree(vp = vp), envir = bbEnv)

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
      assign("genomeLabel_grobs", setChildren(get("genomeLabel_grobs", envir = bbEnv), children = gList(line, labels)), envir = bbEnv)

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

    assign("genomeLabel_grobs", addGrob(gTree = get("genomeLabel_grobs", envir = bbEnv), child = seqGrobs), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(scale)) scale <- NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(fontsize)) fontsize <- NULL
  if(missing(fontcolor)) fontcolor <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(commas)) commas <- NULL
  if(missing(sequence)) sequence <- NULL
  if(missing(axis)) axis <- NULL
  if(missing(boxWidth)) boxWidth <- NULL
  if(missing(tcl)) tcl <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if x/y arguments are missing (could be in object)
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL

  ## Compile all parameters into an internal object
  bb_genomeLabelInternal <- structure(list(x = x, y = y, length = length, plot = plot, chrom = chrom, chromstart = chromstart, chromend = chromend, just = just, scale = scale,
                                           assembly = assembly, fontsize = fontsize, fontcolor = fontcolor, linecolor = linecolor, commas = commas,
                                           sequence = sequence, axis = axis, boxWidth = boxWidth, ticks = ticks, tcl = tcl, default.units = default.units), class = "bb_genomeLabelInternal")
  bb_genomeLabelInternal <- parseParams(bb_params = params, object_params = bb_genomeLabelInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_genomeLabelInternal$just)) bb_genomeLabelInternal$just <- c("left", "top")
  if(is.null(bb_genomeLabelInternal$scale)) bb_genomeLabelInternal$scale <- "bp"
  if(is.null(bb_genomeLabelInternal$assembly)) bb_genomeLabelInternal$assembly <- "hg19"
  if(is.null(bb_genomeLabelInternal$fontsize)) bb_genomeLabelInternal$fontsize <- 10
  if(is.null(bb_genomeLabelInternal$fontcolor)) bb_genomeLabelInternal$fontcolor <- "black"
  if(is.null(bb_genomeLabelInternal$linecolor)) bb_genomeLabelInternal$linecolor <- "black"
  if(is.null(bb_genomeLabelInternal$commas)) bb_genomeLabelInternal$commas <- TRUE
  if(is.null(bb_genomeLabelInternal$sequence)) bb_genomeLabelInternal$sequence <- TRUE
  if(is.null(bb_genomeLabelInternal$axis)) bb_genomeLabelInternal$axis <- "x"
  if(is.null(bb_genomeLabelInternal$boxWidth)) bb_genomeLabelInternal$boxWidth <- 0.5
  if(is.null(bb_genomeLabelInternal$tcl)) bb_genomeLabelInternal$tcl <- 0.5
  if(is.null(bb_genomeLabelInternal$default.units)) bb_genomeLabelInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================


  bb_genomeLabel <- structure(list(chrom = NULL, chromstart = NULL, chromend = NULL,
                                   assembly = NULL, x = bb_genomeLabelInternal$x, y = bb_genomeLabelInternal$y,
                                   width = NULL, height = NULL, just = bb_genomeLabelInternal$just, scale = bb_genomeLabelInternal$scale, grobs = NULL,
                                   gp = gpar(fontsize = bb_genomeLabelInternal$fontsize, col = bb_genomeLabelInternal$linecolor, fontcolor = bb_genomeLabelInternal$fontcolor, ...)), class = "bb_genomeLabel")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_genomeLabel$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_genomeLabel$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)


  if(is.null(bb_genomeLabelInternal$plot)){
    if (is.null(bb_genomeLabelInternal$chrom)){
      stop("Neither `plot` nor `chrom` arguments specified.", call. = FALSE)
    }

    if(is.null(bb_genomeLabelInternal$chromstart)){
      stop("Neither `plot` nor `chromstart` arguments specified.", call. = FALSE)
    }

    if(is.null(bb_genomeLabelInternal$chromend)){
      stop("Neither `plot` nor `chromend` arguments specified.", call. = FALSE)
    }

    if (is.null(bb_genomeLabelInternal$length)){
      stop("If not specifying `plot` input, please provide `length`.", call. = FALSE)
    }

    bb_genomeLabel$chrom <- bb_genomeLabelInternal$chrom
    bb_genomeLabel$chromstart <- bb_genomeLabelInternal$chromstart
    bb_genomeLabel$chromend <- bb_genomeLabelInternal$chromend
    bb_genomeLabel$assembly <- bb_genomeLabelInternal$assembly

  } else {

    ## Check that input plot is a valid type of plot to be annotated
    if (class(bb_genomeLabelInternal$plot) != "bb_manhattan"){

      ## Manhattan plots can do whole genome assembly but other plots can't
      inputNames <- attributes(bb_genomeLabelInternal$plot)$names
      if (!("chrom" %in% inputNames) | !("chromstart" %in% inputNames) | !("chromend" %in% inputNames)){

        stop("Invalid input plot. Please input a plot that has genomic coordinates associated with it.", call. = FALSE)

      }

    }

    bb_genomeLabel$chrom <- bb_genomeLabelInternal$plot$chrom
    bb_genomeLabel$chromstart <- bb_genomeLabelInternal$plot$chromstart
    bb_genomeLabel$chromend <- bb_genomeLabelInternal$plot$chromend
    bb_genomeLabel$assembly <- bb_genomeLabelInternal$plot$assembly
    bb_genomeLabelInternal$length <- bb_genomeLabelInternal$plot$width
    if (bb_genomeLabelInternal$axis == "y"){
      bb_genomeLabelInternal$length <- bb_genomeLabelInternal$plot$height
    }

  }

  check_bbpage(error = "Cannot plot a genome label without a BentoBox page.")
  errorcheck_bb_genomeLabel(scale = bb_genomeLabel$scale, ticks = bb_genomeLabelInternal$ticks, object = bb_genomeLabel, axis = bb_genomeLabelInternal$axis)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  bb_genomeLabel$assembly <- parse_bbAssembly(assembly = bb_genomeLabel$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================
  if (!"unit" %in% class(bb_genomeLabel$x)){

    if (!is.numeric(bb_genomeLabel$x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_genomeLabelInternal$default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_genomeLabel$x <- unit(bb_genomeLabel$x, bb_genomeLabelInternal$default.units)

  }

  if (!"unit" %in% class(bb_genomeLabel$y)){

    if (!is.numeric(bb_genomeLabel$y)){

      stop("y-coordinate is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_genomeLabelInternal$default.units)){

      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_genomeLabel$y <- unit(bb_genomeLabel$y, bb_genomeLabelInternal$default.units)

  }


  if (!"unit" %in% class(bb_genomeLabelInternal$length)){

    if (!is.numeric(bb_genomeLabelInternal$length)){

      stop("Length is neither a unit object or a numeric value. Cannot place object.", call. = FALSE)

    }

    if (is.null(bb_genomeLabelInternal$default.units)){

      stop("Length detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_genomeLabelInternal$length <- unit(bb_genomeLabelInternal$length, bb_genomeLabelInternal$default.units)

  }


  # ======================================================================================================================================================================================
  # SET UP PAGE/SCALE
  # ======================================================================================================================================================================================

  ## Determine scale of labels
  if (bb_genomeLabel$scale == "bp"){
    fact = 1
    format = "d"
  }
  if (bb_genomeLabel$scale == "Mb"){
    fact = 1000000
    format = NULL
    warning("Chromosome start and stop will be rounded.", call. = FALSE)
  }
  if (bb_genomeLabel$scale == "Kb"){
    fact = 1000
    format = "d"
    warning("Chromosome start and stop will be rounded.", call. = FALSE)
  }

  tgH <- convertHeight(heightDetails(textGrob(label = bb_genomeLabel$scale, x = 0.5, y = 0.5, default.units = "npc", gp = bb_genomeLabel$gp)),
                       unitTo = get("page_units", envir = bbEnv))
  seq_height <- heightDetails(textGrob(label = "A", x = 0.5, y = 0.5, default.units = "npc", gp = gpar(fontsize = bb_genomeLabel$gp$fontsize - 2)))
  seq_height <- convertHeight(seq_height + 0.05*seq_height, unitTo = get("page_units", envir = bbEnv))

  # ======================================================================================================================================================================================
  # SET PARAMETERS
  # ======================================================================================================================================================================================
  ########## If single chrom/chromstart/chromend label - comma parsing
  if (length(bb_genomeLabel$chrom) == 1){

    commaLabels <- comma_labels(object = bb_genomeLabel, commas = bb_genomeLabelInternal$commas, format = format, fact = fact)
    chromstartlabel <- commaLabels[[1]]
    chromendlabel <- commaLabels[[2]]

  }
  ########## END comma parsing

  ########## Determine appropriate scaling of nucleotides
  seqType <- NULL
  if (length(bb_genomeLabel$chrom) == 1 & bb_genomeLabelInternal$sequence == TRUE){

    if (bb_genomeLabelInternal$axis == "x"){
      labelWidth <- convertWidth(bb_genomeLabelInternal$length, unitTo = "inches", valueOnly = T)
      bpWidth <- convertWidth(widthDetails(textGrob(label = "A", x = 0.5, y = 0.5, default.units = "npc",
                                                    gp = gpar(fontsize = bb_genomeLabel$gp$fontsize - 2))),
                              unitTo = "inches", valueOnly = T)
      seqRange <- bb_genomeLabel$chromend - bb_genomeLabel$chromstart
      seqWidth <- bpWidth*seqRange


      if (seqWidth <= labelWidth){
        seqType <- "letters"
      } else if (seqWidth/labelWidth <= 9){
        seqType <- "boxes"
      }
    }
  }
  ########## END nucleotide scaling

  ########## Check for BSgenome packages and reset seqType if necessary
  if (!is.null(seqType)){
    if (!is.null(bb_genomeLabel$assembly$BSgenome)){
      bsChecks <- check_loadedPackage(package = bb_genomeLabel$assembly$BSgenome, message = paste(paste0("`", bb_genomeLabel$assembly$BSgenome,"`"), "not loaded. Sequence information will not be displayed."))
      if (bsChecks == FALSE){
        seqType <- NULL
      }
    } else {
      warning("No `BSgenome` package found for the input assembly. Sequence information cannot be displayed.", call. = FALSE)
      seqType <- NULL
    }
  }
  ########## END check for BSgenome packages

  ##########  Total label dimensions, taking into account tick and sequence height
  if (!is.null(bb_genomeLabelInternal$ticks)){
    tick_height <- tgH*(bb_genomeLabelInternal$tcl)
    depth <- convertHeight(tgH + tick_height + 0.5*tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

  } else {
    tick_height <- NULL
    depth <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
  }

  if (!is.null(seqType)){

    seq_height <- convertHeight(seq_height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    bb_genomeLabelInternal$depth <- unit(depth + seq_height, get("page_units", envir = bbEnv))

  } else {

    bb_genomeLabelInternal$depth <- unit(depth, get("page_units", envir = bbEnv))
  }
  ##########  END label length and depth


  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================
  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_genomeLabel", length(grep(pattern = "bb_genomeLabel", x = currentViewports)) + 1)

  ## Make viewport
  vp <- parse_viewport(plot = bb_genomeLabelInternal$plot, object = bb_genomeLabel, length = bb_genomeLabelInternal$length, depth = depth, seqType = seqType,
                       seqHeight = seq_height, vp_name = vp_name, just = bb_genomeLabel$just, axis = bb_genomeLabelInternal$axis)

  # ======================================================================================================================================================================================
  # GROBS AND GTREE
  # ======================================================================================================================================================================================

  ## Chrom/chromstart/chromend grobs
  if (length(bb_genomeLabel$chrom) == 1){

    chrom_grobs(tgH = tgH, ticks = bb_genomeLabelInternal$ticks, tickHeight = tick_height, seqType = seqType, scale = bb_genomeLabel$scale,
                chromLabel = bb_genomeLabel$chrom,
                startLabel = chromstartlabel, endLabel = chromendlabel, height = depth, object = bb_genomeLabel, vp = vp)

    ## Sequence grobs if applicable
    if (!is.null(seqType)){

      seq_grobs(object = bb_genomeLabel, seqHeight = seq_height, seqType = seqType, assembly = bb_genomeLabel$assembly, chromLabel = bb_genomeLabel$chrom, vp = vp,
                boxWidth = bb_genomeLabelInternal$boxWidth)

    }

  } else {
    ## Whole genome grobs for input Manhattan plot
    genome_grobs(plot = bb_genomeLabelInternal$plot, object = bb_genomeLabel, tgH = tgH, vp = vp)

  }



  # ======================================================================================================================================================================================
  # ASSIGN GROBS TO SCALE OBJECT
  # ======================================================================================================================================================================================

  bb_genomeLabel$grobs <- get("genomeLabel_grobs", envir = bbEnv)
  grid.draw(bb_genomeLabel$grobs)

  # ======================================================================================================================================================================================
  # ASSIGN DIMENSIONS BASED ON AXIS
  # ======================================================================================================================================================================================

  bb_genomeLabel$width <- bb_genomeLabelInternal$length
  bb_genomeLabel$height <- bb_genomeLabelInternal$depth
  if (bb_genomeLabelInternal$axis == "y"){
    bb_genomeLabel$width <- bb_genomeLabelInternal$depth
    bb_genomeLabel$height <- bb_genomeLabelInternal$length
  }

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_genomeLabel)
}
