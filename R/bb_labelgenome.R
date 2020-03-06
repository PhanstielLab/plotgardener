#' Adds genome coordinates to the axis of a plot
#'
#' @param plot plot to annotate
#' @param x A unit object specifying x-location
#' @param y A unit object specifying y-location
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
#' @export

bb_labelGenome <- function(plot, x, y, just = c("left", "top"), rotation = 0,
                           scale = "bp", fontsize = 10, fontcolor = "black", linecolor = "black", lwd = 1, fontfamily = "", commas = TRUE, ticks = NULL, tcl = 0.5){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_labelgenome
  errorcheck_bb_labelgenome <- function(plot, scale, ticks){

    ## Check that input plot is actually plotted
    currentViewports <- current_viewports()
    if (!plot$grobs$vp$name %in% currentViewports){

      stop("Input plot is not plotted.", call. = FALSE)

    }

    ## Check that scale is an appropriate value
    if (!scale %in% c("bp", "Kb", "Mb")){

      stop("Invalid \'scale\'. Options are \'bp\', \'Kb\', or \'Mb\'.", call. = FALSE)

    }

    ## If giving tick values, make sure they fall within the chromstart to chromend range
    if (!is.null(ticks)){

      if (range(ticks)[1] < plot$chromstart | range(ticks[2]) > plot$chromend){

        stop("Given tick locations do not fall within the genomic range.", call. = FALSE)

      }

    }

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_genome_label <- structure(list(chrom = NULL, chromstart = plot$chromstart, chromend = plot$chromend, x = x, y = y, width = NULL, height = NULL, scale = scale, grobs = NULL,
                                    gpar = list(fontsize = fontsize, fontcolor = fontcolor, linecolor = linecolor, lwd = lwd, fontfamily = fontfamily)), class = "genome_label")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot add a genome label without a BentoBox page.")
  errorcheck_bb_labelgenome(plot = plot, scale = scale, ticks = ticks)

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
  chromNumber <- gsub(pattern = "chr", replacement = "", x = plot$chrom)
  chrom <- paste0("chr", chromNumber)
  bb_genome_label$chrom <- chrom

  if (commas == TRUE){
    chromstartlabel <- formatC(round(plot$chromstart/fact, 1), format = format, big.mark = ",")
    chromendlabel <- formatC(round(plot$chromend/fact, 1), format = format, big.mark = ",")
  } else {
    chromstartlabel <- round(plot$chromstart/fact, 1)
    chromendlabel <- round(plot$chromend/fact, 1)
  }

  bb_genome_label$width <- plot$width
  if (!is.null(ticks)){
    tick_height <- tgH*tcl
    height <- convertHeight(tgH + tick_height + 0.5*tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    bb_genome_label$height <- unit(height, get("page_units", envir = bbEnv))
  } else {
    height <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    bb_genome_label$height <- unit(tgH, get("page_units", envir = bbEnv))
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
  vp <- viewport(width = page_coords$width, height = page_coords$height,
                 x = page_coords$x, y = page_coords$y,
                 just = just,
                 name = vp_name,
                 xscale = c(plot$chromstart, plot$chromend),
                 yscale = c(0, height),
                 angle = rotation)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("label_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # GROBS
  # ======================================================================================================================================================================================

  if (!is.null(ticks)){

    ## check to make sure the ticks actually fall within the region
    tgH <- convertHeight(tgH, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    tick_height <- convertHeight(tick_height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    x_coords <- ticks
    y_coords <- rep(height-tick_height, length(ticks))
    tickGrobs <- segmentsGrob(x0 = x_coords, y0 = rep(height, length(ticks)), x1 = x_coords, y1 = y_coords, gp = gpar(col = linecolor, lwd = lwd), default.units = "native")

    line <- segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = linecolor, lwd = lwd))
    chromLab <- textGrob(label = chrom, x = 0.5, y = unit(height-(tick_height + 0.5*tgH), "native"),
                         gp = gpar(fontface = "bold", fontsize = fontsize, col = fontcolor, fontfamily = fontfamily), just = c("center", "top"))
    startLab <- textGrob(label = paste(chromstartlabel, scale, sep = " "), x = 0, y = unit(height-(tick_height + 0.5*tgH), "native"), just = c("left", "top"),
                         gp = gpar(fontsize = fontsize, col = linecolor, fontfamily = fontfamily))
    endLab <- textGrob(label = paste(chromendlabel, scale, sep = " "), x = 1, y = unit(height-(tick_height + 0.5*tgH), "native"), just = c("right","top"),
                       gp = gpar(fontsize = fontsize, col = linecolor, fontfamily = fontfamily))

    assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children =  gList(line, chromLab, startLab, endLab, tickGrobs)), envir = bbEnv)

  } else {

    line <- segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = linecolor, lwd = lwd))
    chromLab <- textGrob(label = chrom, x = 0.5, y = 0.85, gp = gpar(fontface = "bold", fontsize = fontsize, col = fontcolor, fontfamily = fontfamily), just = c("center", "top"))
    startLab <- textGrob(label = paste(chromstartlabel, scale, sep = " "), x = 0, y = 0.85, just = c("left", "top"),
                         gp = gpar(fontsize = fontsize, col = linecolor, fontfamily = fontfamily))
    endLab <- textGrob(label = paste(chromendlabel, scale, sep = " "), x = 1, y = 0.85, just = c("right","top"),
                       gp = gpar(fontsize = fontsize, col = linecolor, fontfamily = fontfamily))

    assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)

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
