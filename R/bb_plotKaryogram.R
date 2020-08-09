#' plots and highlights a chromosome with its cytobands
#'
#' @param assembly desired genome assembly
#' @param chrom chromsome to plot
#' @param start highlight start
#' @param end highlight end
#' @param highlightCol fillcolor and linecolor for highlight box
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param just A string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @return Function will return a bb_karyogram object
#' @export
bb_plotKaryogram <- function(assembly = "hg19", chrom, start = NULL, end = NULL, highlightCol = "red", x = NULL, y = NULL, width = NULL, height = NULL,
                             just = c("left", "top"), default.units = "inches", draw = TRUE,...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  errorcheck_bbKaryogram <- function(start, end){

    if (!is.null(start)){
      if (is.null(end)){
        stop('\'start\' provided without \'end\'.', call. = FALSE)
      }
    }

    if (!is.null(end)){
      if (is.null(start)){
        stop('\'end\' provided without \'start\'.', call. = FALSE)
      }
    }

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  karyogram_plot <- structure(list(chrom = chrom, width = width, height = height,
                               x = x, y = y, justification = just, grobs = NULL, assembly = assembly), class = "bb_karyogram")
  attr(x = karyogram_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = karyogram_plot)
  errorcheck_bbKaryogram(start = start, end = end)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  karyogram_plot <- defaultUnits(object = karyogram_plot, default.units = default.units)

  # ======================================================================================================================================================================================
  # GET APPROPRIATE BUILD DATA
  # ======================================================================================================================================================================================

  if (assembly == "hg19"){

    data <- bb_hg19cyto
    genome <- bb_hg19
    leftBezier <- -750000
    rightBezier <- 500000

  }

  # ======================================================================================================================================================================================
  # SUBSET FOR CHROMOSOME
  # ======================================================================================================================================================================================

  minFirstBand <- min(data[which(data$chromStart == 0),]$chromEnd)
  minLastBand <- data[which(data$chromEnd %in% genome$length),]
  minLastBand <- min(minLastBand$chromEnd - minLastBand$chromStart)

  data <- data[which(data[,1] == chrom),]
  chromLength <- genome[which(genome$chrom == chrom),]$length

  # ======================================================================================================================================================================================
  # ASSIGN COLORS
  # ======================================================================================================================================================================================

  data$color <- "black"
  data[which(data$gieStain == "acen"),]$color <- "#802c28"
  data[which(data$gieStain == "stalk"),]$color <- "#6a7ea1"
  data[which(data$gieStain == "gneg"),]$color <- "#FFFFFF"
  data[which(data$gieStain == "gpos25"),]$color <- "#CFCFCF"
  data[which(data$gieStain == "gpos50"),]$color <- "#9F9F9F"
  data[which(data$gieStain == "gpos75"),]$color <- "#6F6F6F"
  data[which(data$gieStain == "gpos100"),]$color <- "#404040"

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_karyogram", length(grep(pattern = "bb_karyogram", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(0.15, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   xscale = c(0, chromLength),
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      vp$name <- "bb_karyogram1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = karyogram_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   xscale = c(0, chromLength),
                   just = just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("karyogram_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # CHROMOSOME GROBS
  # ======================================================================================================================================================================================

  ## Each end of chrom is curved
  firstBand <- data[which(data$chromStart == 0),]
  lastBand <- data[which(data$chromEnd == chromLength),]

  band1a <- bezierGrob(x = c(minFirstBand, leftBezier, leftBezier, minFirstBand),
                              y = unit(c(0, 0, 1, 1), "npc"),
                              default.units = "native",
                              gp = gpar(col = "black"))
  band1a_points <- bezierPoints(band1a)
  band1aGrob <- polygonGrob(x = convertX(band1a_points$x, "native"), y = convertY(band1a_points$y, "npc"),
                               gp = gpar(fill = firstBand$color, col = NA))
  band1b <- rectGrob(x = minFirstBand,
                     y = unit(0.5, "npc"),
                     width = firstBand$chromEnd-minFirstBand,
                     height = unit(1, "npc"),
                     just = "left",
                     gp = gpar(fill = firstBand$color, col = NA),
                     default.units = "native")

  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = band1aGrob), envir = bbEnv)
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = band1b), envir = bbEnv)
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = band1a), envir = bbEnv)

  band2a <- bezierGrob(x = c(lastBand$chromEnd - minLastBand, lastBand$chromEnd + rightBezier, lastBand$chromEnd + rightBezier, lastBand$chromEnd - minLastBand),
                       y = unit(c(0, 0, 1, 1), "npc"),
                       default.units = "native",
                       gp = gpar(col = "black"))
  band2a_points <- bezierPoints(band2a)
  band2aGrob <- polygonGrob(x = convertX(band2a_points$x, "native"), y = convertY(band2a_points$y, "npc"),
                               gp = gpar(fill = lastBand$color, col = NA))
  band2b <- rectGrob(x = lastBand$chromEnd - minLastBand,
                     y = unit(0.5, "npc"),
                     width = lastBand$chromEnd - minLastBand - lastBand$chromStart,
                     height = unit(1, "npc"),
                     just = "right",
                     gp = gpar(fill = lastBand$color, col = NA),
                     default.units = "native")

  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = band2aGrob), envir = bbEnv)
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = band2b), envir = bbEnv)
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = band2a), envir = bbEnv)

  ## Centromere regions are triangles
  leftCent <- data[which(data$gieStain == "acen"),][1]
  rightCent <- data[which(data$gieStain == "acen"),][2]

  leftCentTri <- polygonGrob(x = c(leftCent$chromStart, leftCent$chromStart, leftCent$chromEnd),
                             y = unit(c(0, 1, 0.5), "npc"),
                             gp = gpar(fill = leftCent$color, col = NA),
                             default.units = "native")
  rightCentTri <- polygonGrob(x = c(rightCent$chromStart, rightCent$chromEnd, rightCent$chromEnd),
                             y = unit(c(0.5, 1, 0), "npc"),
                             gp = gpar(fill = rightCent$color, col = NA),
                             default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = leftCentTri), envir = bbEnv)
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = rightCentTri), envir = bbEnv)

  ## Other chrom regions are rectangles
  innerChrom <- data[which(data$gieStain != "acen"),]
  innerChrom <- innerChrom[2:(nrow(innerChrom)-1),]
  innerChrom$width <- innerChrom$chromEnd - innerChrom$chromStart

  innerCyto <- rectGrob(x = innerChrom$chromStart, y = unit(0.5, "npc"),
                        width = innerChrom$width, height = unit(1, "npc"),
                        just = "left", gp = gpar(fill = innerChrom$color, col = NA),
                        default.units = "native")

  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = innerCyto), envir = bbEnv)

  ## Black border around chromosome
  topL_border <- segmentsGrob(x0 = minFirstBand, y0 = unit(1, "npc"),
                              x1 = leftCent$chromStart, y1 = unit(1, "npc"),
                              gp = gpar(col = "black"),
                              default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = topL_border), envir = bbEnv)

  topR_border <- segmentsGrob(x0 = rightCent$chromEnd, y0 = unit(1, "npc"),
                              x1 = lastBand$chromEnd - minLastBand, y1 = unit(1, "npc"),
                              gp = gpar(col = "black"),
                              default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = topR_border), envir = bbEnv)

  bottomL_border <- segmentsGrob(x0 = minFirstBand, y0 = unit(0, "npc"),
                              x1 = leftCent$chromStart, y1 = unit(0, "npc"),
                              gp = gpar(col = "black"),
                              default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = bottomL_border), envir = bbEnv)

  bottomR_border <- segmentsGrob(x0 = rightCent$chromEnd, y0 = unit(0, "npc"),
                              x1 = lastBand$chromEnd - minLastBand, y1 = unit(0, "npc"),
                              gp = gpar(col = "black"),
                              default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = bottomR_border), envir = bbEnv)

  topLeftCent <- segmentsGrob(x0 = leftCent$chromStart, y0 = unit(1, "npc"),
                               x1 = leftCent$chromEnd, y1 = unit(0.5, "npc"),
                               gp = gpar(col = "black"),
                               default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = topLeftCent), envir = bbEnv)

  topRightCent <- segmentsGrob(x0 = rightCent$chromStart, y0 = unit(0.5, "npc"),
                              x1 = rightCent$chromEnd, y1 = unit(1, "npc"),
                              gp = gpar(col = "black"),
                              default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = topRightCent), envir = bbEnv)

  bottomLeftCent <- segmentsGrob(x0 = leftCent$chromStart, y0 = unit(0, "npc"),
                              x1 = leftCent$chromEnd, y1 = unit(0.5, "npc"),
                              gp = gpar(col = "black"),
                              default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = bottomLeftCent), envir = bbEnv)

  bottomRightCent <- segmentsGrob(x0 = rightCent$chromStart, y0 = unit(0.5, "npc"),
                               x1 = rightCent$chromEnd, y1 = unit(0, "npc"),
                               gp = gpar(col = "black"),
                               default.units = "native")
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = bottomRightCent), envir = bbEnv)

  # ======================================================================================================================================================================================
  # HIGHLIGHT BOX
  # ======================================================================================================================================================================================

  if (!is.null(start) & !is.null(end)){


    highlightGrob <- rectGrob(x = start,
                              y = unit(0.5, "npc"),
                              width = end-start,
                              height = unit(1, "npc"),
                              just = "left",
                              gp = gpar(fill = highlightCol, col = highlightCol, alpha = 0.5,...),
                              default.units = "native")
    assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = highlightGrob), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("karyogram_grobs", envir = bbEnv))

  }
  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  karyogram_plot$grobs <-  get("karyogram_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(karyogram_plot)

}
