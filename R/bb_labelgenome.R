#' Adds genome coordinates to the axis of a plot
#'
#' @param plot plot to annotate
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param just justification
#' @param rotation angle of rotation of label
#' @param scale scale of the plot; options are "bp", "Kb", or "Mb"
#' @param fontsize fontsize for labels (in points)
#' @param fontcolor color for all text/lines
#' @param linecolor linecolor
#' @param lwd linewidth
#' @param fontfamily fontfamily for text
#' @param commas A logical value indicating whether to include commas in start and stop labels
#' @param ticks Specified locations of ticks
#' @param tcl Length of tickmark as fraction of text height
#' @param default.units A string indicating the default units to use if x or y are only given as numeric vectors
#' @export

bb_labelGenome <- function(plot, x, y, just = c("left", "top"), rotation = 0,
                           scale = "bp", fontsize = 10, fontcolor = "black", linecolor = "black", lwd = 1, fontfamily = "", commas = TRUE, ticks = NULL, tcl = 0.5, default.units = "inches"){

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

  ## Define a function that parses the viewport for genome assembly vs. chrom/chromstart/chromend label
  parse_viewport <- function(plot, object, height, page_coords, vp_name, just, rotation){

    if (!is.null(object$chrom)){

      vp <- viewport(width = page_coords$width, height = page_coords$height,
                     x = page_coords$x, y = page_coords$y,
                     just = just,
                     name = vp_name,
                     xscale = c(object$chromstart, object$chromend),
                     yscale = c(0, height),
                     angle = rotation)
    } else {

      ## get the offsets based on spacer for the assembly
      if (object$assembly == "hg19"){
        assembly_data <- bb_hg19
      }

      ## get the offsets based on spacer for the assembly
      offsetAssembly <- parse_assembly(assemblyData = assembly_data, space = plot$space)
      cumsums <- cumsum(as.numeric(assembly_data[,2]))
      spacer <- cumsums[length(cumsum(as.numeric(assembly_data[,2])))] * plot$space

      vp <- viewport(width = page_coords$width, height = page_coords$height,
                     x = page_coords$x, y = page_coords$y,
                     just = just,
                     name = vp_name,
                     xscale = c(0, max(offsetAssembly[,4]) + spacer),
                     yscale = c(0, height),
                     angle = rotation)
    }

    return(vp)
  }

  ## Define a function that makes tick, line, and text grobs for chrom/chromstart/chromend labels
  chrom_grobs <- function(tgH, ticks, tickHeight, scale, chromLabel, startLabel, endLabel, height, object){

    if (!is.null(ticks)){

      ## check to make sure the ticks actually fall within the region
      tgH <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
      tick_height <- convertHeight(tickHeight, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
      x_coords <- ticks
      y_coords <- rep(height-tick_height, length(ticks))
      tickGrobs <- segmentsGrob(x0 = x_coords, y0 = rep(height, length(ticks)),
                                x1 = x_coords, y1 = y_coords,
                                gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd), default.units = "native")
      line <- segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd))
      chromLab <- textGrob(label = chromLabel, x = 0.5, y = unit(height-(tick_height + 0.5*tgH), "native"),
                           gp = gpar(fontface = "bold", fontsize = object$gpar$fontsize, col = object$gpar$fontcolor,
                                     fontfamily = object$gpar$fontfamily), just = c("center", "top"))
      startLab <- textGrob(label = paste(startLabel, scale, sep = " "), x = 0, y = unit(height-(tick_height + 0.5*tgH), "native"), just = c("left", "top"),
                           gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))
      endLab <- textGrob(label = paste(endLabel, scale, sep = " "), x = 1, y = unit(height-(tick_height + 0.5*tgH), "native"), just = c("right","top"),
                         gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))

      assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children =  gList(line, chromLab, startLab, endLab, tickGrobs)), envir = bbEnv)

    } else {

      line <- segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd))
      chromLab <- textGrob(label = chromLabel, x = 0.5, y = 0.85,
                           gp = gpar(fontface = "bold", fontsize = object$gpar$fontsize, col = object$gpar$fontcolor, fontfamily = object$gpar$fontfamily), just = c("center", "top"))
      startLab <- textGrob(label = paste(startLabel, scale, sep = " "), x = 0, y = 0.85, just = c("left", "top"),
                           gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))
      endLab <- textGrob(label = paste(endLabel, scale, sep = " "), x = 1, y = 0.85, just = c("right","top"),
                         gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$linecolor, fontfamily = object$gpar$fontfamily))

      assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)

    }

  }

  ## define a function that makes line and text grobs for whole assembly labels
  genome_grobs <- function(plot, object){

    ## Get internal assembly data
    if (object$assembly == "hg19"){
      assembly_data <- bb_hg19
    }

    ## Get the offsets based on spacer for the assembly
    offsetAssembly <- parse_assembly(assemblyData = assembly_data, space = plot$space)

    ## Get the centers of each chrom
    chromCenters <- (offsetAssembly[,3] + offsetAssembly[,4]) / 2

    line <- segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = object$gpar$linecolor, lwd = object$gpar$lwd))
    labels <- textGrob(label = gsub("chr", "", offsetAssembly[,1]), x = chromCenters, y = unit(0.85, "npc"), just = c("center", "top"),
                       gp = gpar(fontsize = object$gpar$fontsize, col = object$gpar$fontcolor, fontfamily = object$gpar$fontfamily),
                       default.units = "native")
    assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, labels)), envir = bbEnv)

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

  # ======================================================================================================================================================================================
  # SET PARAMETERS
  # ======================================================================================================================================================================================

  ## If chrom/chromstart/chromend label - comma parsing
  if (!is.null(bb_genome_label$chrom)){

    if (commas == TRUE){
      chromstartlabel <- formatC(round(plot$chromstart/fact, 1), format = format, big.mark = ",")
      chromendlabel <- formatC(round(plot$chromend/fact, 1), format = format, big.mark = ",")
    } else {
      chromstartlabel <- round(plot$chromstart/fact, 1)
      chromendlabel <- round(plot$chromend/fact, 1)
    }

  }

  ## Label dimensions
  bb_genome_label$width <- plot$width
  if (!is.null(ticks)){
    tick_height <- tgH*tcl
    height <- convertHeight(tgH + tick_height + 0.5*tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    bb_genome_label$height <- unit(height, get("page_units", envir = bbEnv))
  } else {
    tick_height <- NULL
    height <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    bb_genome_label$height <- unit(tgH, get("page_units", envir = bbEnv))
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

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = bb_genome_label)

  ## Make viewport
  vp <- parse_viewport(plot = plot, object = bb_genome_label, height = height,
                       page_coords = page_coords, vp_name = vp_name, just = just, rotation = rotation)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("label_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # GROBS
  # ======================================================================================================================================================================================

  ## Chrom/chromstart/chromend grobs
  if (!is.null(bb_genome_label$chrom)){

    chrom_grobs(tgH = tgH, ticks = ticks, tickHeight = tick_height, scale = scale, chromLabel = bb_genome_label$chrom,
                startLabel = chromstartlabel, endLabel = chromendlabel, height = height, object = bb_genome_label)

  } else {
    ## Whole genome grobs
    genome_grobs(plot = plot, object = bb_genome_label)

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
