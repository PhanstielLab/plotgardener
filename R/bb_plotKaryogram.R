#' plots and highlights a chromosome with its cytobands
#'
#' @param chrom chromsome to plot
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param assembly desired genome assembly
#' @param orientation "v" (vertical) or "h" (horizontal) orientation
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
bb_plotKaryogram <- function(chrom, params = NULL, assembly = "hg19", orientation = "h", start = NULL, end = NULL, highlightCol = "red", x = NULL, y = NULL, width = NULL, height = NULL,
                             just = c("left", "top"), default.units = "inches", draw = TRUE,...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  errorcheck_bbKaryogram <- function(start, end, orientation){

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


    if(!orientation %in% c("v", "h")){
      stop("Invalid /'orientation/' parameter. Options are 'v' or 'h'.", call. = FALSE)

    }

  }

  ## Define a function that draws bands that fall within left curved regions
  curvedBands_left <- function(df, xCurve, yCurve, ymax){
    start <- as.numeric(df[2])
    end <- as.numeric(df[3])
    col <- df[6]
    if(end > max(xCurve)){
      xpoints <- c(xCurve[which(xCurve >= start)], end, end)
      ypoints <- c(yCurve[which(xCurve >= start)], 0, ymax)
    } else {
      xpoints <- xCurve[which(xCurve >= start & xCurve <= end)]
      ypoints <- yCurve[which(xCurve >= start & xCurve <= end)]
    }

    curvedGrob <- polygonGrob(x = xpoints, y = ypoints,
                              default.units = "native",
                              gp = gpar(fill = col, col = NA))

    assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = curvedGrob), envir = bbEnv)

  }

  ## Define a function that draws bands that fall within right curved regions
  curvedBands_right <- function(df, xCurve, yCurve, ymax){
    start <- as.numeric(df[2])
    end <- as.numeric(df[3])
    col <- df[6]
    if(start < min(xCurve)){
      xpoints <- c(xCurve[which(xCurve <= end)], start, start)
      ypoints <- c(yCurve[which(xCurve <= end)], ymax, 0)
    } else {
      xpoints <- xCurve[which(xCurve >= start & xCurve <= end)]
      ypoints <- yCurve[which(xCurve >= start & xCurve <= end)]
    }

    curvedGrob <- polygonGrob(x = xpoints, y = ypoints,
                              default.units = "native",
                              gp = gpar(fill = col, col = NA))

    assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = curvedGrob), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(orientation)) orientation <- NULL
  if(missing(highlightCol)) highlightCol <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if chrom argument is missing (could be in object)
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_karyInternal <- structure(list(chrom = chrom, assembly = assembly, orientation = orientation, start = start, end = end, highlightCol = highlightCol,
                                    x = x, y = y, width = width, height = height, just = just, default.units = default.units,
                                    draw = draw), class = "bb_karyInternal")

  bb_karyInternal <- parseParams(bb_params = params, object_params = bb_karyInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_karyInternal$assembly)) bb_karyInternal$assembly <- "hg19"
  if(is.null(bb_karyInternal$orientation)) bb_karyInternal$orientation <- "h"
  if(is.null(bb_karyInternal$highlightCol)) bb_karyInternal$highlightCol <- "red"
  if(is.null(bb_karyInternal$just)) bb_karyInternal$just <- c("left", "top")
  if(is.null(bb_karyInternal$default.units)) bb_karyInternal$default.units <- "inches"
  if(is.null(bb_karyInternal$draw)) bb_karyInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  karyogram_plot <- structure(list(chrom = bb_karyInternal$chrom, width = bb_karyInternal$width, height = bb_karyInternal$height,
                               x = bb_karyInternal$x, y = bb_karyInternal$y, justification = bb_karyInternal$just, grobs = NULL, assembly = bb_karyInternal$assembly), class = "bb_karyogram")
  attr(x = karyogram_plot, which = "plotted") <- bb_karyInternal$draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(karyogram_plot$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)

  check_placement(object = karyogram_plot)
  errorcheck_bbKaryogram(start = bb_karyInternal$start, end = bb_karyInternal$end, orientation = bb_karyInternal$orientation)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  karyogram_plot <- defaultUnits(object = karyogram_plot, default.units = bb_karyInternal$default.units)

  # ======================================================================================================================================================================================
  # GET APPROPRIATE BUILD DATA
  # ======================================================================================================================================================================================

  ## EDIT HERE
  if (karyogram_plot$assembly == "hg19"){

    data <- bb_hg19cyto
    genome <- bb_hg19

  }

  # ======================================================================================================================================================================================
  # SUBSET FOR CHROMOSOME
  # ======================================================================================================================================================================================

  minFirstBand <- min(data[which(data$chromStart == 0),]$chromEnd)
  minLastBand <- data[which(data$chromEnd %in% genome$length),]
  minLastBand <- min(minLastBand$chromEnd - minLastBand$chromStart)

  data <- data[which(data[,1] == karyogram_plot$chrom),]
  chromLength <- genome[which(genome$chrom == karyogram_plot$chrom),]$length

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
  if (is.null(karyogram_plot$x) & is.null(karyogram_plot$y)){

    height <- 0.10
    width <- 0.9

    scaleRatio <- width/height
    yscale <- chromLength/scaleRatio

    if (bb_karyInternal$orientation == "h"){
      vp <- viewport(height = unit(height, "snpc"), width = unit(width, "snpc"),
                     x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                     xscale = c(0, chromLength),
                     yscale = c(0, yscale),
                     just = "center",
                     name = vp_name)
    } else {
      height <- 0.9
      width <- 0.10
      vp <- viewport(height = unit(width, "snpc"), width = unit(height, "snpc"),
                     x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                     xscale = c(0, chromLength),
                     yscale = c(0, yscale),
                     angle = -90,
                     just = "center",
                     name = vp_name)

    }


    if (bb_karyInternal$draw == TRUE){

      vp$name <- "bb_karyogram1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = karyogram_plot)

    height <- convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    width <- convertWidth(page_coords$width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    scaleRatio <- max(c(width, height))/min(c(width, height))
    yscale <- chromLength/scaleRatio

    if (bb_karyInternal$orientation == "h"){

      vp <- viewport(height = page_coords$height, width = page_coords$width,
                     x = page_coords$x, y = page_coords$y,
                     xscale = c(0, chromLength),
                     yscale = c(0, yscale),
                     just = bb_karyInternal$just,
                     name = vp_name)

    } else {
      ## Make viewport based on user inputs
      vpOG <- viewport(height = page_coords$height, width = page_coords$width,
                     x = page_coords$x, y = page_coords$y,
                     just = bb_karyInternal$just)

      ## Convert viewport to bottom left (bottom right of horizontal)
      vp_bottom <- vp_bottomLeft(viewport = vpOG)

      ## Make new rotated viewport
      vp <- viewport(height = page_coords$width, width = page_coords$height,
                     x = vp_bottom[[1]], y = vp_bottom[[2]],
                     xscale = c(0, chromLength),
                     yscale = c(0, yscale),
                     angle = -90,
                     just = c("right", "bottom"),
                     name = vp_name)

    }

  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("karyogram_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # CHROMOSOME GROBS
  # ======================================================================================================================================================================================

  ## Generate points along curves for the ends
  r = vp$yscale[2] * 0.5
  leftAngles <- seq(pi/2, 3*pi/2, pi/500)
  rightAngles <- seq(3*pi/2, 5*pi/2, pi/500)

  leftXpoints <- r + r*cos(leftAngles)
  leftYpoints <- r + r*sin(leftAngles)

  rightXpoints <- (chromLength - r) + r*cos(rightAngles)
  rightYpoints <- r + r*sin(rightAngles)

  ## FIRST BAND ##
  firstBand <- data[which(data$chromStart == 0),]
  data <- subset(data, data$chromStart != 0)

  if (firstBand$chromEnd > max(leftXpoints)){
    firstBand_Xpoints <- c(leftXpoints, firstBand$chromEnd, firstBand$chromEnd)
    firstBand_Ypoints <- c(leftYpoints, 0, vp$yscale[2])
  } else {
    firstBand_Xpoints <- leftXpoints[which(leftXpoints <= firstBand$chromEnd)]
    firstBand_Ypoints <- leftYpoints[which(leftXpoints <= firstBand$chromEnd)]
  }

  firstBand_grob <- polygonGrob(x = firstBand_Xpoints, y = firstBand_Ypoints,
                                default.units = "native",
                                gp = gpar(fill = firstBand$color, col = NA))

  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = firstBand_grob), envir = bbEnv)

  ## CENTER BANDS ##
  leftCent <- data[which(data$gieStain == "acen"),][1]
  leftCent_length <- leftCent$chromEnd - leftCent$chromStart
  rightCent <- data[which(data$gieStain == "acen"),][2]
  rightCent_length <- rightCent$chromEnd - rightCent$chromStart
  centerX <- leftCent$chromEnd
  data <- subset(data, data$gieStain != "acen")

  ## Generate points along curves for the centers
  centerleftXpoints <- (centerX - r*0.75) + r*cos(rightAngles)
  centerleftYpoints <- r + r*sin(rightAngles)
  centerleftYpoints <- centerleftYpoints[which(centerleftXpoints <= (rightCent$chromEnd - 0.5*rightCent_length))]
  centerleftXpoints <- centerleftXpoints[which(centerleftXpoints <= (rightCent$chromEnd - 0.5*rightCent_length))]

  centerrightXpoints <- (centerX + r*0.75) + r*cos(leftAngles)
  centerrightYpoints <- r + r*sin(leftAngles)
  centerrightYpoints <- centerrightYpoints[which(centerrightXpoints >= (leftCent$chromStart + 0.5*leftCent_length))]
  centerrightXpoints <- centerrightXpoints[which(centerrightXpoints >= (leftCent$chromStart + 0.5*leftCent_length))]

  ## CENTER LEFT BAND ##
  if (leftCent$chromStart < min(centerleftXpoints)){
    leftCent_Xpoints <- c(centerleftXpoints, leftCent$chromStart, leftCent$chromStart)
    leftCent_Ypoints <- c(centerleftYpoints, vp$yscale[2], 0)
  } else {
    leftCent_Xpoints <- centerleftXpoints[which(centerleftXpoints >= leftCent$chromStart)]
    leftCent_Ypoints <- centerleftYpoints[which(centerleftXpoints >= leftCent$chromStart)]
  }

  leftCent_grob <- polygonGrob(x = leftCent_Xpoints, y = leftCent_Ypoints,
                               default.units = "native",
                               gp = gpar(fill = leftCent$color, col = NA))

  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = leftCent_grob), envir = bbEnv)

  ## CENTER RIGHT BAND ##

  if (rightCent$chromEnd > max(centerrightXpoints)){
    rightCent_Xpoints <- c(centerrightXpoints, rightCent$chromEnd, rightCent$chromEnd)
    rightCent_Ypoints <- c(centerrightYpoints, 0, vp$yscale[2])
  } else {
    rightCent_Xpoints <- centerrightXpoints[which(centerrightXpoints <= rightCent$chromEnd)]
    rightCent_Ypoints <- centerrightYpoints[which(centerrightXpoints <= rightCent$chromEnd)]
  }

  rightCent_grob <- polygonGrob(x = rightCent_Xpoints, y = rightCent_Ypoints,
                                default.units = "native",
                                gp = gpar(fill = rightCent$color, col = NA))

  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = rightCent_grob), envir = bbEnv)

  ## LAST BAND ##
  lastBand <- data[which(data$chromEnd == chromLength),]
  data <- subset(data, data$chromEnd != chromLength)

  if (lastBand$chromStart < min(rightXpoints)){
    lastBand_Xpoints <- c(rightXpoints, lastBand$chromStart, lastBand$chromStart)
    lastBand_Ypoints <- c(rightYpoints, vp$yscale[2], 0)
  } else {
    lastBand_Xpoints <- rightXpoints[which(rightXpoints >= lastBand$chromStart)]
    lastBand_Ypoints <- rightYpoints[which(rightXpoints >= lastBand$chromStart)]
  }

  lastBand_grob <- polygonGrob(x = lastBand_Xpoints, y = lastBand_Ypoints,
                               default.units = "native",
                               gp = gpar(fill = lastBand$color, col = NA))

  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = lastBand_grob), envir = bbEnv)

  ## GET ANY BANDS THAT FALL WITHIN CURVED REGIONS ##
  leftcurvedBands <- data[which(data$chromStart < max(leftXpoints)),]
  inleftcurvedBands <- data[which(data$chromEnd > min(centerleftXpoints) & data$chromEnd <= centerX),]
  inrightcurvedBands <- data[which(data$chromStart < max(centerrightXpoints) & data$chromStart >= centerX),]
  rightcurvedBands <- data[which(data$chromStart >= min(rightXpoints)),]

  if(nrow(leftcurvedBands > 0)) invisible(apply(leftcurvedBands, 1, curvedBands_left, xCurve = leftXpoints, yCurve = leftYpoints, ymax = vp$yscale[2]))
  if(nrow(inleftcurvedBands > 0)) invisible(apply(inleftcurvedBands, 1, curvedBands_right, xCurve = centerleftXpoints, yCurve = centerleftYpoints, ymax = vp$yscale[2]))
  if(nrow(inrightcurvedBands > 0)) invisible(apply(inrightcurvedBands, 1, curvedBands_left, xCurve = centerrightXpoints, yCurve = centerrightYpoints, ymax = vp$yscale[2]))
  if(nrow(rightcurvedBands > 0)) invisible(apply(rightcurvedBands, 1, curvedBands_right, xCurve = rightXpoints, yCurve = rightYpoints, ymax = vp$yscale[2]))

  ## REMAINING BANDS ##
  data <- suppressMessages(dplyr::anti_join(data, leftcurvedBands))
  data <- suppressMessages(dplyr::anti_join(data, inleftcurvedBands))
  data <- suppressMessages(dplyr::anti_join(data, inrightcurvedBands))
  data <- suppressMessages(dplyr::anti_join(data, rightcurvedBands))

  data$width <- data$chromEnd - data$chromStart

  rectBands <- rectGrob(x = data$chromStart, y = unit(0.5, "npc"),
                           width = data$width, height = unit(1, "npc"),
                           just = "left", default.units = "native",
                           gp = gpar(fill = data$color, col = NA))
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = rectBands), envir = bbEnv)


  ## OUTLINE ##

  topIntersectY <- sqrt(r^2 -((r*0.75)^2)) + r
  bottomIntersectY <- -1*sqrt(r^2 -((r*0.75)^2)) + r

  lbottomX <- centerleftXpoints[which(centerleftXpoints <= centerX & centerleftYpoints <= bottomIntersectY)]
  lbottomY <- centerleftYpoints[which(centerleftXpoints <= centerX & centerleftYpoints <= bottomIntersectY)]

  rbottomX <- centerrightXpoints[which(centerrightXpoints >= centerX & centerrightYpoints <= bottomIntersectY)]
  rbottomY <- centerrightYpoints[which(centerrightXpoints >= centerX & centerrightYpoints <= bottomIntersectY)]

  rtopX <- centerrightXpoints[which(centerrightXpoints >= centerX & centerrightYpoints >= topIntersectY)]
  rtopY <- centerrightYpoints[which(centerrightXpoints >= centerX & centerrightYpoints >= topIntersectY)]

  ltopX <- centerleftXpoints[which(centerleftXpoints <= centerX & centerleftYpoints >= topIntersectY)]
  ltopY <- centerleftYpoints[which(centerleftXpoints <= centerX & centerleftYpoints >= topIntersectY)]

  Xoutline <- c(leftXpoints, lbottomX, centerX, rbottomX, rightXpoints, rtopX, centerX, ltopX)
  Youtline <- c(leftYpoints, lbottomY, bottomIntersectY, rbottomY, rightYpoints, rtopY, topIntersectY, ltopY)

  outlineGrob <- polygonGrob(x = Xoutline, y = Youtline,
                             default.units = "native",
                             gp = gpar(fill = NA, col = "black"))
  assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = outlineGrob), envir = bbEnv)



  # ======================================================================================================================================================================================
  # HIGHLIGHT BOX
  # ======================================================================================================================================================================================

  if (!is.null(bb_karyInternal$start) & !is.null(bb_karyInternal$end)){


    highlightGrob <- rectGrob(x = bb_karyInternal$start,
                              y = unit(0.5, "npc"),
                              width = bb_karyInternal$end-bb_karyInternal$start,
                              height = unit(1, "npc"),
                              just = "left",
                              gp = gpar(fill = bb_karyInternal$highlightCol, col = bb_karyInternal$highlightCol, alpha = 0.5,...),
                              default.units = "native")
    assign("karyogram_grobs", addGrob(get("karyogram_grobs", envir = bbEnv), child = highlightGrob), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_karyInternal$draw == TRUE){

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
