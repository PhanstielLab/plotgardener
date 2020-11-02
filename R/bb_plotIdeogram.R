#' plots and highlights a chromosome with its cytobands
#'
#' @param chrom chromsome to plot
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param assembly default genome assembly as a string or a bb_assembly object
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
bb_plotIdeogram <- function(chrom, params = NULL, assembly = "hg19", orientation = "h", start = NULL, end = NULL, highlightCol = "red", x = NULL, y = NULL, width = NULL, height = NULL,
                             just = c("left", "top"), default.units = "inches", draw = TRUE,...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that checks errors for bb_plotIdeogram
  errorcheck_bbIdeogram <- function(start, end, orientation){

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

  ## Define a function to get cytoBand data for a genome assembly
  cytoAssembly <- function(assembly){
    availCytos <- list(hg18 = "cytoBand.Hsapiens.UCSC.hg18", hg19 = "cytoBand.Hsapiens.UCSC.hg19", hg38 = "cytoBand.Hsapiens.UCSC.hg38",
                       mm9 = "cytoBand.Mmusculus.UCSC.mm9", mm10 = "cytoBand.Mmusculus.UCSC.mm10", dm6 = "cytoBand.Dmelanogaster.UCSC.dm6",
                       rn5 = "cytoBand.Rnorvegicus.UCSC.rn5", rn6 = "cytoBand.Rnorvegicus.UCSC.rn6", danRer10 = "cytoBand.Drerio.UCSC.danRer10")

    ## Split TxDb assembly name
    assemblyName <- unlist(strsplit(assembly$TxDb, split = "[.]"))

    if (!any(names(availCytos) %in% assemblyName)){
      warning(paste("CytoBand data not available for the given genome assembly. Ideograms can only be plotted for the following assemblies:", cat(names(availCytos), sep = ", ")), call. = FALSE)
      cytoData <- NULL
      genomeData <- NULL
    } else {

      ## Get name of associated cytoband data (included in package)
      cytoData <- availCytos[[which(names(availCytos) %in% assemblyName)]]

      ## Check that included cytoband data is loaded
      currentLoaded <- ls(envir = globalenv())
      if (!cytoData %in% currentLoaded){
        warning(paste("Assembly cytoBand data not loaded. Run", paste0('`data("',cytoData,'")`'), "to load data."), call. = FALSE)
        cytoData <- NULL
      }

      ## Check that TxDb package is loaded
      txdbChecks <- check_loadedPackage(package = assembly$TxDb, message = paste(paste0("`", assembly$TxDb,"`"), "not loaded. Please install and load to plot Ideogram."))
      if(txdbChecks == FALSE){
        genomeData <- NULL
      } else {
        genomeData <- assembly$TxDb
      }

    }

    return(list(cytoData, genomeData))
  }

  ## Define a function to check that a chromosome name is in an associated TxDb
  checkChroms <- function(chrom, txdb){

    tx_db <- eval(parse(text = txdb))
    txdbChroms <- seqlevels(tx_db)
    if (chrom %in% txdbChroms){
      return(TRUE)
    } else {
      warning(paste(paste0("'", chrom, "'"), "not found in", paste0(txdb, ".")), call. = FALSE)
      return(FALSE)
    }

  }

  ## Define a function to give colors to different gieStains for the various assemblies
  assignColors <- function(data, assembly){

    human <- c("hg18", "hg19", "hg38")
    mouse <- c("mm9", "mm10")
    rat <- c("rn5", "rn6")
    fly <- c("dm6")
    zebrafish <- c("danRer10")

    data$color <- "black"

    if (assembly %in% human){
      if(any(data$gieStain == "acen")) data[which(data$gieStain == "acen"),]$color <- "#802c28"
      if(any(data$gieStain == "stalk")) data[which(data$gieStain == "stalk"),]$color <- "#6a7ea1"
      if(any(data$gieStain == "gneg")) data[which(data$gieStain == "gneg"),]$color <- "#FFFFFF"
      if(any(data$gieStain == "gpos25")) data[which(data$gieStain == "gpos25"),]$color <- "#CFCFCF"
      if(any(data$gieStain == "gpos50")) data[which(data$gieStain == "gpos50"),]$color <- "#9F9F9F"
      if(any(data$gieStain == "gpos75")) data[which(data$gieStain == "gpos75"),]$color <- "#6F6F6F"
      if(any(data$gieStain == "gpos100")) data[which(data$gieStain == "gpos100"),]$color <- "#404040"

    } else if (assembly %in% mouse){
      if(any(data$gieStain == "gneg")) data[which(data$gieStain == "gneg"),]$color <- "#FFFFFF"
      if(any(data$gieStain == "gpos33")) data[which(data$gieStain == "gpos33"),]$color <- "#CFCFCF"
      if(any(data$gieStain == "gpos66")) data[which(data$gieStain == "gpos66"),]$color <- "#9F9F9F"
      if(any(data$gieStain == "gpos75")) data[which(data$gieStain == "gpos75"),]$color <- "#6F6F6F"
      if(any(data$gieStain == "gpos100")) data[which(data$gieStain == "gpos100"),]$color <- "#404040"

    } else if (assembly %in% rat){
      if (assembly == "rn5"){
        if(any(data$gieStain == "gneg")) data[which(data$gieStain == "gneg"),]$color <- "#FFFFFF"
        if(any(data$gieStain == "gpos")) data[which(data$gieStain == "gpos"),]$color <- "#9F9F9F"

      } else {
        data$color <- rep(c("#FFFFFF", "#9F9F9F"), ceiling(nrow(data)/2))[1:nrow(data)]
      }
    } else if (assembly %in% fly){

      data$color <- rep(c("#FFFFFF", "#9F9F9F"), ceiling(nrow(data)/2))[1:nrow(data)]

    } else {
      if(any(data$gieStain == "gneg")) data[which(data$gieStain == "gneg"),]$color <- "#FFFFFF"
    }

    return(data)

  }

  ## Define a function that draws bands that fall within left curved regions
  curvedBands_left <- function(df, xCurve, yCurve, ymax){
    start <- as.numeric(df[2])
    end <- as.numeric(df[3])
    col <- df[8]

    if(end > max(xCurve)){
      xpoints <- c(xCurve[which(xCurve >= start)], end, end)
      ypoints <- c(yCurve[which(xCurve >= start)], 0, ymax)
    } else {
      xpoints <- xCurve[which(xCurve >= start & xCurve <= end)]
      ypoints <- yCurve[which(xCurve >= start & xCurve <= end)]
    }

    if (length(xpoints) > 0 & length(ypoints) > 0){
      curvedGrob <- polygonGrob(x = xpoints, y = ypoints,
                                default.units = "native",
                                gp = gpar(fill = col, col = NA))

      assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = curvedGrob), envir = bbEnv)
    }


  }

  ## Define a function that draws bands that fall within right curved regions
  curvedBands_right <- function(df, xCurve, yCurve, ymax){
    start <- as.numeric(df[2])
    end <- as.numeric(df[3])
    col <- df[8]
    if(start < min(xCurve)){
      xpoints <- c(xCurve[which(xCurve <= end)], start, start)
      ypoints <- c(yCurve[which(xCurve <= end)], ymax, 0)
    } else {
      xpoints <- xCurve[which(xCurve >= start & xCurve <= end)]
      ypoints <- yCurve[which(xCurve >= start & xCurve <= end)]
    }

    if (length(xpoints) > 0 & length(ypoints) > 0){
      curvedGrob <- polygonGrob(x = xpoints, y = ypoints,
                                default.units = "native",
                                gp = gpar(fill = col, col = NA))

      assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = curvedGrob), envir = bbEnv)

    }


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
  bb_ideoInternal <- structure(list(chrom = chrom, assembly = assembly, orientation = orientation, start = start, end = end, highlightCol = highlightCol,
                                    x = x, y = y, width = width, height = height, just = just, default.units = default.units,
                                    draw = draw), class = "bb_ideoInternal")

  bb_ideoInternal <- parseParams(bb_params = params, object_params = bb_ideoInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_ideoInternal$assembly)) bb_ideoInternal$assembly <- "hg19"
  if(is.null(bb_ideoInternal$orientation)) bb_ideoInternal$orientation <- "h"
  if(is.null(bb_ideoInternal$highlightCol)) bb_ideoInternal$highlightCol <- "red"
  if(is.null(bb_ideoInternal$just)) bb_ideoInternal$just <- c("left", "top")
  if(is.null(bb_ideoInternal$default.units)) bb_ideoInternal$default.units <- "inches"
  if(is.null(bb_ideoInternal$draw)) bb_ideoInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  ideogram_plot <- structure(list(chrom = bb_ideoInternal$chrom, width = bb_ideoInternal$width, height = bb_ideoInternal$height,
                               x = bb_ideoInternal$x, y = bb_ideoInternal$y, justification = bb_ideoInternal$just, grobs = NULL, assembly = bb_ideoInternal$assembly), class = "bb_ideogram")
  attr(x = ideogram_plot, which = "plotted") <- bb_ideoInternal$draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(ideogram_plot$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)

  check_placement(object = ideogram_plot)
  errorcheck_bbIdeogram(start = bb_ideoInternal$start, end = bb_ideoInternal$end, orientation = bb_ideoInternal$orientation)


  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  ideogram_plot$assembly <- parse_bbAssembly(assembly = ideogram_plot$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  ideogram_plot <- defaultUnits(object = ideogram_plot, default.units = bb_ideoInternal$default.units)

  # ======================================================================================================================================================================================
  # GET APPROPRIATE BUILD DATA
  # ======================================================================================================================================================================================

  cytoData <- cytoAssembly(assembly = ideogram_plot$assembly)
  data <- cytoData[[1]]
  if(!is.null(data)) data <- get(data)
  genome <- cytoData[[2]]
  if (!is.null(genome)){
    chromCheck <- checkChroms(chrom = bb_ideoInternal$chrom, txdb = genome)
    genome <- eval(parse(text = genome))

  }

  chromLength <- 1
  if (!is.null(data) & !is.null(genome)){

    if (chromCheck == TRUE){

      # ======================================================================================================================================================================================
      # SUBSET FOR CHROMOSOME
      # ======================================================================================================================================================================================

      data <- data[which(data[,1] == ideogram_plot$chrom),]
      data$seqnames <- as.character(data$seqnames)
      data$strand <- as.character(data$strand)
      data$name <- as.character(data$name)
      data$gieStain <- as.character(data$gieStain)
      chromLength <- seqlengths(genome)[[ideogram_plot$chrom]]

      # ======================================================================================================================================================================================
      # ASSIGN COLORS
      # ======================================================================================================================================================================================

      data <- assignColors(data = data, assembly = bb_ideoInternal$assembly)


    }

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_ideogram", length(grep(pattern = "bb_ideogram", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(ideogram_plot$x) & is.null(ideogram_plot$y)){

    height <- 0.10
    width <- 0.9

    scaleRatio <- width/height
    yscale <- chromLength/scaleRatio

    if (bb_ideoInternal$orientation == "h"){
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


    if (bb_ideoInternal$draw == TRUE){

      vp$name <- "bb_ideogram1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = ideogram_plot)

    height <- convertHeight(page_coords$height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    width <- convertWidth(page_coords$width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    scaleRatio <- max(c(width, height))/min(c(width, height))
    yscale <- chromLength/scaleRatio

    if (bb_ideoInternal$orientation == "h"){

      vp <- viewport(height = page_coords$height, width = page_coords$width,
                     x = page_coords$x, y = page_coords$y,
                     xscale = c(0, chromLength),
                     yscale = c(0, yscale),
                     just = bb_ideoInternal$just,
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

  assign("ideogram_grobs", gTree(vp = vp), envir = bbEnv)

  if (!is.null(data) & !is.null(genome)){

    if (chromCheck == TRUE){
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


      if (nrow(data) > 1){

        ## FIRST BAND ##
        firstBand <- data[which(data$start == 1),]
        data <- subset(data, data$start != 1)

        if (firstBand$end > max(leftXpoints)){
          firstBand_Xpoints <- c(leftXpoints, firstBand$end, firstBand$end)
          firstBand_Ypoints <- c(leftYpoints, 0, vp$yscale[2])
        } else {
          firstBand_Xpoints <- leftXpoints[which(leftXpoints <= firstBand$end)]
          firstBand_Ypoints <- leftYpoints[which(leftXpoints <= firstBand$end)]
        }

        firstBand_grob <- polygonGrob(x = firstBand_Xpoints, y = firstBand_Ypoints,
                                      default.units = "native",
                                      gp = gpar(fill = firstBand$color, col = NA))

        assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = firstBand_grob), envir = bbEnv)

        ## LAST BAND ##
        lastBand <- data[which(data$end == chromLength),]
        data <- subset(data, data$end != chromLength)

        if (lastBand$start < min(rightXpoints)){
          lastBand_Xpoints <- c(rightXpoints, lastBand$start, lastBand$start)
          lastBand_Ypoints <- c(rightYpoints, vp$yscale[2], 0)
        } else {
          lastBand_Xpoints <- rightXpoints[which(rightXpoints >= lastBand$start)]
          lastBand_Ypoints <- rightYpoints[which(rightXpoints >= lastBand$start)]
        }

        lastBand_grob <- polygonGrob(x = lastBand_Xpoints, y = lastBand_Ypoints,
                                     default.units = "native",
                                     gp = gpar(fill = lastBand$color, col = NA))

        assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = lastBand_grob), envir = bbEnv)


        if (bb_ideoInternal$assembly %in% c("hg18", "hg19", "hg38")){
          ## CENTER BANDS ##
          leftCent <- data[which(data$gieStain == "acen"),][1,]
          leftCent_length <- leftCent$end - leftCent$start
          rightCent <- data[which(data$gieStain == "acen"),][2,]
          rightCent_length <- rightCent$end - rightCent$start
          centerX <- leftCent$end
          data <- subset(data, data$gieStain != "acen")

          ## Generate points along curves for the centers
          centerleftXpoints <- (centerX - r*0.75) + r*cos(rightAngles)
          centerleftYpoints <- r + r*sin(rightAngles)
          centerleftYpoints <- centerleftYpoints[which(centerleftXpoints <= (rightCent$end - 0.5*rightCent_length))]
          centerleftXpoints <- centerleftXpoints[which(centerleftXpoints <= (rightCent$end - 0.5*rightCent_length))]

          centerrightXpoints <- (centerX + r*0.75) + r*cos(leftAngles)
          centerrightYpoints <- r + r*sin(leftAngles)
          centerrightYpoints <- centerrightYpoints[which(centerrightXpoints >= (leftCent$start + 0.5*leftCent_length))]
          centerrightXpoints <- centerrightXpoints[which(centerrightXpoints >= (leftCent$start + 0.5*leftCent_length))]

          ## CENTER LEFT BAND ##
          if (leftCent$start < min(centerleftXpoints)){
            leftCent_Xpoints <- c(centerleftXpoints, leftCent$start, leftCent$start)
            leftCent_Ypoints <- c(centerleftYpoints, vp$yscale[2], 0)
          } else {
            leftCent_Xpoints <- centerleftXpoints[which(centerleftXpoints >= leftCent$start)]
            leftCent_Ypoints <- centerleftYpoints[which(centerleftXpoints >= leftCent$start)]
          }

          leftCent_grob <- polygonGrob(x = leftCent_Xpoints, y = leftCent_Ypoints,
                                       default.units = "native",
                                       gp = gpar(fill = leftCent$color, col = NA))

          assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = leftCent_grob), envir = bbEnv)

          ## CENTER RIGHT BAND ##

          if (rightCent$end > max(centerrightXpoints)){
            rightCent_Xpoints <- c(centerrightXpoints, rightCent$end, rightCent$end)
            rightCent_Ypoints <- c(centerrightYpoints, 0, vp$yscale[2])
          } else {
            rightCent_Xpoints <- centerrightXpoints[which(centerrightXpoints <= rightCent$end)]
            rightCent_Ypoints <- centerrightYpoints[which(centerrightXpoints <= rightCent$end)]
          }

          rightCent_grob <- polygonGrob(x = rightCent_Xpoints, y = rightCent_Ypoints,
                                        default.units = "native",
                                        gp = gpar(fill = rightCent$color, col = NA))

          assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = rightCent_grob), envir = bbEnv)



          ## GET ANY BANDS THAT FALL WITHIN CENTER CURVED REGIONS ##
          inleftcurvedBands <- data[which(data$end > min(centerleftXpoints) & data$end <= centerX),]
          inrightcurvedBands <- data[which(data$start < max(centerrightXpoints) & data$start >= centerX),]
          if(nrow(inleftcurvedBands > 0)) invisible(apply(inleftcurvedBands, 1, curvedBands_right, xCurve = centerleftXpoints, yCurve = centerleftYpoints, ymax = vp$yscale[2]))
          if(nrow(inrightcurvedBands > 0)) invisible(apply(inrightcurvedBands, 1, curvedBands_left, xCurve = centerrightXpoints, yCurve = centerrightYpoints, ymax = vp$yscale[2]))

          ## REMAINING BANDS ##
          data <- suppressMessages(dplyr::anti_join(data, inleftcurvedBands))
          data <- suppressMessages(dplyr::anti_join(data, inrightcurvedBands))

        }


        ## GET ANY BANDS THAT FALL WITHIN OUTSIDE CURVED REGIONS ##
        leftcurvedBands <- data[which(data$start < max(leftXpoints)),]
        rightcurvedBands <- data[which(data$start >= min(rightXpoints)),]

        if(nrow(leftcurvedBands > 0)) invisible(apply(leftcurvedBands, 1, curvedBands_left, xCurve = leftXpoints, yCurve = leftYpoints, ymax = vp$yscale[2]))
        if(nrow(rightcurvedBands > 0)) invisible(apply(rightcurvedBands, 1, curvedBands_right, xCurve = rightXpoints, yCurve = rightYpoints, ymax = vp$yscale[2]))

        ## REMAINING BANDS ##
        data <- suppressMessages(dplyr::anti_join(data, leftcurvedBands))
        data <- suppressMessages(dplyr::anti_join(data, rightcurvedBands))

        rectBands <- rectGrob(x = data$start, y = unit(0.5, "npc"),
                              width = data$width, height = unit(1, "npc"),
                              just = "left", default.units = "native",
                              gp = gpar(fill = data$color, col = NA))
        assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = rectBands), envir = bbEnv)

      }


      # ======================================================================================================================================================================================
      # OUTLINE GROBS
      # ======================================================================================================================================================================================

      if (bb_ideoInternal$assembly %in% c("hg18", "hg19", "hg38")){

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

      } else {

        Xoutline <- c(leftXpoints, rightXpoints)
        Youtline <- c(leftYpoints, rightYpoints)

      }

      outlineGrob <- polygonGrob(x = Xoutline, y = Youtline,
                                 default.units = "native",
                                 gp = gpar(fill = NA, col = "black"))
      assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = outlineGrob), envir = bbEnv)


    }

    # ======================================================================================================================================================================================
    # HIGHLIGHT BOX
    # ======================================================================================================================================================================================

    if (!is.null(bb_ideoInternal$start) & !is.null(bb_ideoInternal$end)){


      highlightGrob <- rectGrob(x = bb_ideoInternal$start,
                                y = unit(0.5, "npc"),
                                width = bb_ideoInternal$end-bb_ideoInternal$start,
                                height = unit(1, "npc"),
                                just = "left",
                                gp = gpar(fill = bb_ideoInternal$highlightCol, col = bb_ideoInternal$highlightCol, alpha = 0.5,...),
                                default.units = "native")
      assign("ideogram_grobs", addGrob(get("ideogram_grobs", envir = bbEnv), child = highlightGrob), envir = bbEnv)

    }

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_ideoInternal$draw == TRUE){

    grid.draw(get("ideogram_grobs", envir = bbEnv))

  }
  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  ideogram_plot$grobs <-  get("ideogram_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(ideogram_plot)

}
