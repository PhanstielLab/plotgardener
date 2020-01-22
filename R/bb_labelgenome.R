#' Adds genome coordinates to the axis of a plot
#'
#' @param plot plot to annotate
#' @param x A unit object specifying x-location
#' @param y A unit object specifying y-location
#' @param just justification
#' @param units units of inputted x and y
#' @param scale scale of the plot; options are "bp", "Kb", or "Mb"
#' @param side side of the plot to annotate, options are "bottom", "top", "right", or "left"
#' @param fontsize fontsize for labels (in points)
#' @param fontcolor color for all text/lines
#' @param lineAbove TRUE/FALSE indicating if line should be above or below text
#' @param linecolor linecolor
#' @param lwd linewidth


#' @export

bb_labelGenome <- function(plot, x, y, just = c("left", "top"), rotation = 0,
                           scale = "bp", fontsize = 10, fontcolor = "black", lineAbove = T, linecolor = "black", lwd = 1){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_labelgenome
  # errorcheck_bb_labelgenome <- function()

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_genome_label <- structure(list(x = x, y = y, width = NULL, height = NULL, scale = scale, grobs = NULL,
                                    gpar = list(fontsize = fontsize, fontcolor = fontcolor, linecolor = linecolor, lwd = lwd)), class = "genome_label")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot add a genome label without a BentoBox page.")

  # ======================================================================================================================================================================================
  # SET UP PAGE/SCALE
  # ======================================================================================================================================================================================

  ## Determine scale of labels
  if (scale == "bp"){
    fact = 1
  }
  if (scale == "Mb"){
    fact = 1000000
  }
  if (scale == "Kb"){
    fact = 1000
  }

  tg <- textGrob(label = scale, x = 0.5, y = 0.5, default.units = "npc", gp = gpar(fontsize = fontsize))

  # ======================================================================================================================================================================================
  # SET PARAMETERS
  # ======================================================================================================================================================================================

  chrom <- paste0("chr", plot$chrom)
  chromstartlabel <- plot$chromstart/fact
  chromendlabel <- plot$chromend/fact
  bb_genome_label$width <- plot$width
  bb_genome_label$height <- convertHeight(heightDetails(tg), unitTo = get("page_units", envir = bbEnv))

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Name viewport
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_genomeLabel", length(grep(pattern = "bb_genomeLabel", x = current_viewports)) + 1)

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = bb_genome_label)

  ## Make viewport
  vp <- viewport(width = page_coords$width, height = page_coords$height,
                 x = page_coords$x, y = page_coords$y,
                 just = just,
                 name = vp_name)
  #pushViewport(vp)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  assign("label_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE
  # ======================================================================================================================================================================================

  if (lineAbove == T){

    line <- segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1, gp = gpar(col = linecolor, lwd = lwd))
    chromLab <- textGrob(label = chrom, x = 0.5, y = 0.85, gp = gpar(fontface = "bold", fontsize = fontsize, col = fontcolor), just = c("center", "top"))
    startLab <- textGrob(label = paste(round(chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.85, just = c("left", "top"), gp = gpar(fontsize = fontsize, col = linecolor))
    endLab <- textGrob(label = paste(round(chromendlabel, 1), scale, sep = " "), x = 1, y = 0.85, just = c("right","top"), gp = gpar(fontsize = fontsize, col = linecolor))

  } else if (lineAbove == F){

    line <- segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0, gp = gpar(col = linecolor, lwd = lwd))
    chromLab <- textGrob(label = chrom, x = 0.5, y = 0.15, gp = gpar(fontface = "bold", fontsize = fontsize, col = fontcolor), just = c("center", "bottom"))
    startLab <- textGrob(label = paste(round(chromstartlabel, 1), scale, sep = " "), x = 0, y = 0.15, just = c("left", "bottom"), gp = gpar(fontsize = fontsize, col = linecolor))
    endLab <- textGrob(label = paste(round(chromendlabel, 1), scale, sep = " "), x = 1, y = 0.15, just = c("right", "bottom"), gp = gpar(fontsize = fontsize, col = linecolor))

  }

  assign("label_grobs", setChildren(get("label_grobs", envir = bbEnv), children = gList(line, chromLab, startLab, endLab)), envir = bbEnv)

  ## Go back to root viewport
  #upViewport()

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
