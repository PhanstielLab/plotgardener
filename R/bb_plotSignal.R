#' plots signal track data
#'
#' @param signal signal track data to be plotted (bigwig file, bedgraph format dataframe, or bedgraph); can be one or a list of two
#' @param chrom chromosome of region to be plotted as a string (i.e. "chr3")
#' @param chromstart start position
#' @param chromend end position
#' @param range y-range to plot (c(min, max))
#' @param linecolor color(s) of line outlining signal track(s); if using fillcolor and no outline is desired, set to NA
#' @param lwd linewidth of line outlining signal track
#' @param binSize the length of each bin in bp
#' @param binCap TRUE/FALSE whether the function will limit the number of bins to 8,000
#' @param fillcolor color(s) filling in signal track(s)
#' @param assembly desired genome assembly
#' @param transparency numerical value for fillcolor transparency, if fillcolor is specified
#' @param ymax fraction of max y-value to set as height of plot
#' @param scale A logical value indicating whether to include a data scale label in the top left corner
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param just a string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @export
bb_plotSignal <- function(signal, chrom, chromstart = NULL, chromend = NULL, range = NULL, linecolor = "grey", lwd = 1,
                          binSize = NA, binCap = TRUE, fillcolor = NULL, assembly = "hg19", transparency = 1, ymax = 1, scale = FALSE, width = NULL,
                           height = NULL, x = NULL, y = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, ...  ){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_plotSignal
  errorcheck_bb_signaltrack <- function(signal, signaltrack){

    dfChecks <- function(signal){

      if ("data.frame" %in% class(signal) && ncol(signal) != 4){

        stop("Invalid dataframe format.  Input a dataframe with 4 columns in bedgraph format: chrom, chromstart, chromend, counts.", call. = FALSE)

      }

      if (!"data.frame" %in% class(signal)){

        if (!file_ext(signal) %in% c("bw", "bigWig", "bigwig", "bedgraph")){

          stop("Invalid input. File must have a valid bigwig or bedgraph extension.", call. = FALSE)

        }

      }

    }

    if (length(class(signal)) == 1 && class(signal) == "list"){

      if (length(signal) > 2){

        stop("Invalid signal input. More than 2 signals provided.", call. = FALSE)

      }

      invisible(lapply(signal, dfChecks))

      } else {

        dfChecks(signal = signal)

      }

    ## Can't have only one NULL chromstart or chromend
    if ((is.null(signal_track$chromstart) & !is.null(signal_track$chromend)) | (is.null(signal_track$chromend) & !is.null(signal_track$chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.", call. = FALSE)

    }

    if (!is.null(signal_track$chromstart) & !is.null(signal_track$chromend)){

      if (signal_track$chromstart > signal_track$chromend){

        stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)

      }
    }

    if (!is.null(signaltrack$range)){

      ## range needs to be a vector
      if (!is.vector(signaltrack$range)){

        stop("\'range\' must be a vector of length 2.", call. = FALSE)

      }

      ## range vector needs to be length 2
      if (length(signaltrack$range) != 2){

        stop("\'range\' must be a vector of length 2.", call. = FALSE)

      }

      ## range vector needs to be numbers
      if (!is.numeric(signaltrack$range)){

        stop("\'range\' must be a vector of two numbers.", call. = FALSE)

      }

      ## second value should be larger than the first value
      if (signaltrack$range[1] >= signaltrack$range[2]){

        stop("\'range\' must be a vector of two numbers in which the 2nd value is larger than the 1st.", call. = FALSE)

      }

    }

  }

## Define a function that reads in signal data for bb_plotSignal
  read_signal <- function(signal, signaltrack){

    ## if .bw file, read in with bb_readBigwig
    if (!"data.frame" %in% class(signal)){

      signal <- bb_readBigwig(filename = signal, chrom = signaltrack$chrom, chromstart = signaltrack$chromstart, chromend = signaltrack$chromend)
      signal <- signal[,c(1:3,6)]

    } else {

      message("Reading in signal dataframe.  Assuming dataframe in bedgraph format with \'chrom\' in column1, \'chromstart\' in column 2,
       \'chromend\' in column 3, and \'counts\' in column 4.")

      ## Make dataframe
      signal <- as.data.frame(signal)

    }

    return(signal)

  }

  ## Define a function that formats/filters signal data
  format_data <- function(signal, signaltrack){

    ## Ensure the chromosome is a character
    signal[,1] <- as.character(signal[,1])

    ## Filter for desired region
    signal <- signal[which(signal[,1] == signaltrack$chrom & ((signal[,2] > signaltrack$chromstart & signal[,2] < signaltrack$chromend | signal[,3] > signaltrack$chromstart &
                                                                 signal[,3] < signaltrack$chromend ))), (2:4)]
    ## Remove any duplicate rows
    signal <- signal[!duplicated(signal),]

    return(signal)

  }

  ## Define a function that checks and adjust the number of bins
  check_binNum <- function(signaltrack, binCap){

    binNum = (signaltrack$chromend - signaltrack$chromstart)/signaltrack$binSize
    signaltrack$binNum <- binNum

    ## Scale back binNum and print warning if binNum is greater than 8000
    if (binNum > 8000 && signaltrack$binCap == TRUE){
      updated_binNum <- 8000
      updated_binSize <- (signaltrack$chromend - signaltrack$chromstart)/binNum
      signaltrack$binNum <- updated_binNum
      signaltrack$binSize <- updated_binSize
      warning(paste0("Too many bins: adjusting to 8000 bins of size ", binSize, ". To override try binCap = FALSE."), call. = FALSE)
    }

    ## Scale bin size to 1 if binNum is larger than span
    if (binNum > (signaltrack$chromend - signaltrack$chromstart)){
      updated_binNum <- (signaltrack$chromend - signaltrack$chromstart)
      updated_binSize <- 1
      signaltrack$binNum <- updated_binNum
      signaltrack$binSize <- updated_Size
      warning(paste0("Number of bins larger than plot length: adjusting to ", binNum, " bins of size 1."), call. = FALSE)
    }

    return(signaltrack)
  }

  ## Define a function that bins, links, sorts, and combines data
  parseData <- function(signal, signaltrack){

    ## Define a function that finds the max signal values for each bin
    bin_signal <- function(line, signal){

      line <- as.integer(line)
      list <- c(0)
      list <- append(list, signal[,3][which((signal[,1] >= line[1] & signal[,1] < line[2]) |
                                              signal[,2] > line[1] & signal[,2] <= line[2] |
                                              signal[,1] < line[1] & signal[,2] > line[2])])

      return(max(list))
    }

    ## Define a function that adds slightly negative values for polygon plotting if using fill
    add_neg_vals <- function(signal, signaltrack){

      if (!is.null(signaltrack$fillcolor)){

        signal <- rbind(c(min(signal[,1]), -0.00001), signal)
        signal <- rbind(signal, c(max(signal[,1]), -0.00001))

      }

      return(signal)
    }

    # ===============================================================================================================================================
    # BIN DATA
    # ===============================================================================================================================================
    binned_signal <- data.frame("chromstart" = seq(signaltrack$chromstart, signaltrack$chromend - signaltrack$binSize, signaltrack$binSize),
                                "chromend" = seq(signaltrack$chromstart + signaltrack$binSize, signaltrack$chromend, signaltrack$binSize),
                                "counts" = 0)
    binned_signal[,3] = apply(binned_signal, 1, bin_signal, signal = signal)

    ## Use binned data as signal track
    newSignal <- binned_signal

    # ================================================================================================================================================
    # LINKING REGIONS
    # ================================================================================================================================================

    linking_regions <- cbind(newSignal[1:(nrow(newSignal) - 1), 2], newSignal[2:nrow(newSignal), 1])

    linking_regions <- matrix(linking_regions[which(linking_regions[,1] != linking_regions[,2]),], ncol = 2)

    if (nrow(linking_regions) > 0){

      linking_regions <- cbind(linking_regions, 0)
      ## Make column names the same
      colnames(linking_regions)[(1:3)] <- c("chromstart", "chromend", "counts")

      ## Add linking regions to signaltrack
      newSignal <- rbind(newSignal, linking_regions)
    }

    # ================================================================================================================================================
    # SORT AND COMBINE DATA
    # ================================================================================================================================================

    ## Sort data
    newSignal <- newSignal[order(newSignal[,1]),]

    ## Convert two columns to one
    newSignal <- cbind(as.vector(t(newSignal[,c(1, 2)])), as.vector(t(newSignal[,c(3, 3)])))

    # ================================================================================================================================================
    # NEG VALUES FOR PROPER POLYGON PLOTTING
    # ================================================================================================================================================

    newSignal <- add_neg_vals(signal = newSignal, signaltrack = signaltrack)

    return(newSignal)
  }

  ## Define a function that parses out one vs. two fillcolors/linecolors
  parseColors <- function(color){

    if (length(color) >= 2){
      posCol <- color[1]
      negCol <- color[2]
    } else {
      posCol <- color
      negCol <- color
    }

    return(list(posCol, negCol))

  }

  ## Define a function that makes a grob for a signal (pos/neg)
  sigGrob <- function(signal, fillCol, lineCol, lwd, transparency){

    if (!is.null(fillCol)){

      finalcolor <- makeTransparent(color = fillCol, alpha = transparency)
      sigGrob <- polygonGrob(x = signal[,1], y = signal[,2], gp = gpar(fill = finalcolor, lwd = lwd, col = lineCol), default.units = "native")

    } else {

      sigGrob <- segmentsGrob(x0 = signal[c(1:1-length(signal[,1])), 1], y0 = signal[c(1:1-length(signal[,2])), 2],
                              x1 = signal[c(2:length(signal[,1])), 1], y1 = signal[c(2:length(signal[,2])), 2], gp = gpar(col = lineCol, lwd = lwd), default.units = "native")

    }


    ## Add grob to gtree
    assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = sigGrob), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  signal_track <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, range = range,
                                  linecolor = linecolor, lwd = lwd, fillcolor = fillcolor,
                                  binSize = binSize, binNum = NULL, ymax = ymax,
                                  width = width, height = height, x = x, y = y, justification = just, grobs = NULL, assembly = "hg19"), class = "bb_signal")
  attr(x = signal_track, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = signal_track)
  errorcheck_bb_signaltrack(signal = signal, signaltrack = signal_track)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  signal_track <- defaultUnits(object = signal_track, default.units = default.units)


  # ======================================================================================================================================================================================
  # WHOLE CHROM
  # ======================================================================================================================================================================================

  if (is.null(chromstart) & is.null(chromend)){
    if (assembly == "hg19"){
      genome <- bb_hg19
    }

    signal_track$chromstart <- 1
    signal_track$chromend <- genome[which(genome$chrom == chrom),]$length

  }

  # ======================================================================================================================================================================================
  # SET BINSIZE
  # ======================================================================================================================================================================================

  if (is.na(binSize) == TRUE){

    binSize <- (signal_track$chromend - signal_track$chromstart)/2000
    signal_track$binSize <- binSize

  }

  # ======================================================================================================================================================================================
  # CHECK AND ADJUST BIN NUMBER
  # ======================================================================================================================================================================================

  signal_track <- check_binNum(signaltrack = signal_track, binCap = binCap)

  # ======================================================================================================================================================================================
  # READ IN, FORMAT, FILTER, BIN, LINK AND SORT DATA
  # ======================================================================================================================================================================================

  if (length(class(signal)) == 1 && class(signal) == "list"){

    signal <- lapply(signal, read_signal, signaltrack = signal_track)
    signal <- lapply(signal, format_data, signaltrack = signal_track)
    posSignal <- signal[[1]]
    negSignal <- signal[[2]]
    split <- TRUE

  } else {

    signal <- read_signal(signal = signal, signaltrack = signal_track)
    signal <- format_data(signal = signal, signaltrack = signal_track)

    if (any(signal[,3] < 0)){

      posSignal <- signal[which(signal[,3] >= 0),]
      negSignal <- signal[which(signal[,3] < 0),]
      negSignal[,3] <- negSignal[,3]*-1
      split <- TRUE

    } else {

      posSignal <- signal
      split <- FALSE
    }

  }

  # ======================================================================================================================================================================================
  # BIN, LINK, AND SORT DATA AND FIX Y-LIMITS
  # ======================================================================================================================================================================================

  if (split == TRUE){

    if (nrow(posSignal) >= 2){
      posSignal2 <- parseData(signal = posSignal, signaltrack = signal_track)
      if (is.null(signal_track$range)){
        signal_track$range[2] <- signal_track$ymax * max(posSignal2[,2])
      }

    } else {
      signal_track$range[2] <- 1
    }

    if (nrow(negSignal) >= 2){
      negSignal2 <- parseData(signal = negSignal, signaltrack = signal_track)
      negSignal2[,2] <- negSignal2[,2]*-1
      if (is.na(signal_track$range[1])){
        signal_track$range[1] <- signal_track$ymax * min(negSignal2[,2])
      }

    } else {
      signal_track$range[1] <- -1
    }


  } else {

    if (nrow(posSignal) >= 2){
      posSignal2 <- parseData(signal = posSignal, signaltrack = signal_track)

      if (is.null(signal_track$range)){
        signal_track$range <- c(0, ymax*max(posSignal2[,2]))
      }

    } else {
      signal_track$range <- c(0, 1)
    }

  }


  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_signal", length(grep(pattern = "bb_signal", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(0.25, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = c(signal_track$chromstart, signal_track$chromend),
                   yscale = c(signal_track$range[1], signal_track$range[2]),
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      vp$name <- "bb_signal1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = signal_track)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(signal_track$chromstart, signal_track$chromend),
                   yscale = c(signal_track$range[1], signal_track$range[2]),
                   just = just,
                   name = vp_name)
  }


  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("signal_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (split == TRUE){

    fills <- parseColors(signal_track$fillcolor)
    lines <- parseColors(signal_track$linecolor)

    if (nrow(posSignal) >= 2){
      sigGrob(signal = posSignal2, fillCol = fills[[1]], lineCol = lines[[1]], lwd = signal_track$lwd, transparency = transparency)
    } else {
      posGrob <- segmentsGrob(x0 = 0, y0 = 0, x1 = 1, y1 = 0, gp = gpar(col = lines[[1]], lwd = signal_track$lwd))
      assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = posGrob), envir = bbEnv)
      warning("Not enough top signal data to plot.", call. = FALSE)
    }


    if (nrow(negSignal) >= 2){
      sigGrob(signal = negSignal2, fillCol = fills[[2]], lineCol = lines[[2]], lwd = signal_track$lwd, transparency = transparency)
    } else {
      negGrob <- segmentsGrob(x0 = 0, y0 = 0, x1 = 1, y1 = 0, gp = gpar(col = lines[[2]], lwd = signal_track$lwd))
      assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = negGrob), envir = bbEnv)
      warning("Not enough bottom signal data to plot.", call. = FALSE)
    }

    lineGrob <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = 0, y1 = 0,
                             gp = gpar(col = "black", lwd = signal_track$lwd + 0.5), default.units = "native")
    assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = lineGrob), envir = bbEnv)

  } else {

    if (nrow(posSignal) >= 2){
      sigGrob(signal = posSignal2, fillCol = signal_track$fillcolor, lineCol = signal_track$linecolor,
              lwd = signal_track$lwd, transparency = transparency)
    } else {
      signalGrob <- segmentsGrob(x0 = 0, y0 = 0, x1 = 1, y1 = 0, gp = gpar(col = signal_track$linecolor, lwd = signal_track$lwd))
      assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = signalGrob), envir = bbEnv)
      warning("Not enough data within range to plot.", call. = FALSE)
    }

  }
  # ======================================================================================================================================================================================
  # SCALE
  # ======================================================================================================================================================================================

  ## Add scale of the range of data in the top left corner
  if (scale == TRUE){
    upperLim <- round(signal_track$range[2], digits = 4)
    lowerLim <- round(signal_track$range[1], digits = 4)

    scaleGrob <- textGrob(label = paste0("[", lowerLim, " - ", upperLim, "]"), just = c("left", "top"), x = 0, y = 1, gp = gpar(col = "black", fontsize = 8))

    ## Add grob to gtree
    assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = scaleGrob), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("signal_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  signal_track$grobs <-  get("signal_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(signal_track)

}

