#' plots signal track data
#'
#' @param signal signal track data to be plotted (bigwig file, bedgraph format dataframe, or bedgraph)
#' @param chrom chromosome of region to be plotted as a string (i.e. "chr3")
#' @param chromstart start position
#' @param chromend end position
#' @param range y-range to plot (c(min, max))
#' @param linecolor color of line outlining signal track
#' @param lwd linewidth of line outlining signal track
#' @param binSize the length of each bin in bp
#' @param binCap TRUE/FALSE whether the function will limit the number of bins to 8,000
#' @param fillcolor if want plot to be filled, specify fillcolor
#' @param transparency if want to specify a transparency for the fill, specify transparency
#' @param ymax fraction of max y value to set as height of plot
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param just a string or numeric vector specifying the justification of the viewport relative to its (x, y) location: "left", "right", "centre", "center", "bottom", "top"
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param draw A logical value indicating whether graphics output should be produced
#'
#' @export
bb_plotSignal <- function(signal, chrom, chromstart, chromend, range = NULL, linecolor = "grey", lwd = 1, yaxis = FALSE,
                           binSize = NA, binCap = TRUE, fillcolor = NULL, transparency = NULL, ymax = 1, width = NULL,
                           height = NULL, x = NULL, y = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, ...  ){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_plotSignal
  errorcheck_bb_signaltrack <- function(signal, signaltrack){

    if (class(signal) %in% "data.frame" && ncol(signal) != 4){

      stop("Invalid dataframe format.  Input a dataframe with 4 columns in bedgraph format: chrom, chromstart, chromend, counts.")

    }

    if (!class(signal) %in% "data.frame"){

      if (!file_ext(signal) %in% c("bw", "bigWig", "bigwig", "bedgraph")){

        stop("Invalid input. File must have a valid bigwig or bedgraph extension")

      }

    }

    if (chromstart > chromend){

      stop("\'chromstart\' should not be larger than \'chromend\'.")

    }

    if (!is.null(range)){

      ## range needs to be a vector
      if (!is.vector(signaltrack$range)){

        stop("\'range\' must be a vector of length 2.")

      }

      ## range vector needs to be length 2
      if (length(signaltrack$range) != 2){

        stop("\'range\' must be a vector of length 2.")

      }

      ## range vector needs to be numbers
      if (!is.numeric(signaltrack$range)){

        stop("\'range\' must be a vector of two numbers.")

      }

      ## second value should be larger than the first value
      if (signaltrack$range[1] >= signaltrack$range[2]){

        stop("\'range\' must be a vector of two numbers in which the 2nd value is larger than the 1st.")

      }



    }




  }

  ## Define a function to check range of data in dataframe
  check_signal_dataframe <- function(signal, signaltrack){

    if (min(signal[,2]) > signaltrack$chromstart | max(signal[,2]) < signaltrack$chromend | min(signal[,3]) > signaltrack$chromstart | max(signal[,3]) < signaltrack$chromend){

        warning("Data is incomplete for the specified range.")

    }

  }

  ## Define a function that reads in signal data for bb_plotSignal
  read_signal <- function(signal, signaltrack){

    ## if .hic file, read in with bb_rhic
    if (!(class(signal) %in% "data.frame")){

      signal <- bb_readBigwig(filename = signal, chrom = chrom, chromstart = chromstart, chromend = chromend)
      signal <- signal[,c(1,2,3,6)]

    } else {

      message("Reading in dataframe.  Assuming dataframe in bedgraph format with \'chrom\' in column1, \'chromstart\' in column 2,
       \'chromend\' in column 3, and \'counts\' in column 4.")

      ## check range of data in dataframe
      #check_signal_dataframe(signal = signal, signaltrack = signaltrack)

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
      warning(paste0("Too many bins: adjusting to 8000 bins of size ", binSize, ". To override try binCap = FALSE."))
    }

    ## Scale bin size to 1 if binNum is larger than span
    if (binNum > (signaltrack$chromend - signaltrack$chromstart)){
      updated_binNum <- (signaltrack$chromend - signaltrack$chromstart)
      updated_binSize <- 1
      signaltrack$binNum <- updated_binNum
      signaltrack$binSize <- updated_Size
      warning(paste0("Number of bins larger than plot length: adjusting to ", binNum, " bins of size 1."))
    }

    return(signaltrack)
  }

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

  ## Define a function that adjusts the range
  adjust_range <- function(signal, signaltrack){

    if (is.null(signaltrack$range)){

      signaltrack$range <- c(0, signaltrack$ymax * max(signal[,2]))

    }

    return(signaltrack)
  }

  ## Define a function that makes the signal grobs
  signal_grobs <- function(signal, signaltrack, transparency){

    if (!is.null(signaltrack$fillcolor)){

      rgbcol = col2rgb(signaltrack$fillcolor)

      if (is.null(transparency)){

        transparency = 1
      }

      finalcolor = rgb(rgbcol[1], rgbcol[2], rgbcol[3], alpha = transparency * 255, maxColorValue = 255)

      signalGrob <- polygonGrob(x = signal[,1], y = signal[,2], gp = gpar(fill = finalcolor, lwd = signaltrack$lwd, col = signaltrack$linecolor), default.units = "native")

    } else {

      #signalGrob <- polygonGrob(x = signal[,1], y = signal[,2], gp = gpar(fill = NA, lwd = signaltrack$lwd, col = signaltrack$linecolor), default.units = "native")

      signalGrob <- segmentsGrob(x0 = signal[c(1:1-length(signal[,1])), 1], y0 = signal[c(1:1-length(signal[,2])), 2],
                    x1 = signal[c(2:length(signal[,1])), 1], y1 = signal[c(2:length(signal[,2])), 2], gp = gpar(col = signaltrack$linecolor, lwd = signaltrack$lwd), default.units = "native")

    }


    ## Add grob to gtree
    assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = signalGrob), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  signal_track <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, range = range,
                                  linecolor = linecolor, lwd = lwd, fillcolor = fillcolor,
                                  binSize = binSize, binNum = NULL, ymax = ymax,
                                  width = width, height = height, x = x, y = y, justification = just, grobs = NULL), class = "bb_signal")
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
  # SET BINSIZE
  # ======================================================================================================================================================================================

  if (is.na(binSize) == TRUE){

    binSize <- (chromend - chromstart)/2000
    signal_track$binSize <- binSize

  }

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  signal <- read_signal(signal = signal, signaltrack = signal_track)

  # ======================================================================================================================================================================================
  # FORMAT AND FILTER DATA
  # ======================================================================================================================================================================================

  signal <- format_data(signal = signal, signaltrack = signal_track)

  if (nrow(signal) >= 2){

    # ======================================================================================================================================================================================
    # CHECK AND ADJUST BIN NUMBER
    # ======================================================================================================================================================================================

    signal_track <- check_binNum(signaltrack = signal_track, binCap = binCap)

    # ======================================================================================================================================================================================
    # BIN DATA
    # ======================================================================================================================================================================================

    binned_signal <- data.frame(seq(signal_track$chromstart, signal_track$chromend - signal_track$binSize, signal_track$binSize),
                                seq(signal_track$chromstart + signal_track$binSize, signal_track$chromend, signal_track$binSize),
                                rep(0, times = signal_track$binNum))

    ## Add column names
    colnames(binned_signal) = c("chromstart", "chromend", "counts")

    binned_signal[,3] = apply(binned_signal, 1, bin_signal, signal = signal)

    ## Use binned data as signal track
    signal <- binned_signal

    # ======================================================================================================================================================================================
    # LINKING REGIONS
    # ======================================================================================================================================================================================

    linking_regions <- cbind(signal[1:(nrow(signal) - 1), 2], signal[2:nrow(signal), 1])

    linking_regions <- matrix(linking_regions[which(linking_regions[,1] != linking_regions[,2]),], ncol = 2)

    if (nrow(linking_regions) > 0){

      linking_regions <- cbind(linking_regions, 0)
      ## Make column names the same
      colnames(linking_regions)[(1:3)] <- c("chromstart", "chromend", "counts")

      ## Add linking regions to signaltrack
      signal <- rbind(signal, linking_regions)
    }

    # ======================================================================================================================================================================================
    # SORT AND COMBINE DATA
    # ======================================================================================================================================================================================

    ## Sort data
    signal <- signal[order(signal[,1]),]

    ## Convert two columns to one
    signal <- cbind(as.vector(t(signal[,c(1, 2)])), as.vector(t(signal[,c(3, 3)])))

    # ======================================================================================================================================================================================
    # Y-LIMITS
    # ======================================================================================================================================================================================

    ## Determine add slightly negative value to both ends to ensure proper polygon plotting
    signal <- add_neg_vals(signal = signal, signaltrack = signal_track)

    ## Determine the y-limits
    signal_track <- adjust_range(signal = signal, signaltrack = signal_track)


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
                   xscale = c(chromstart, chromend), yscale = c(signal_track$range[1], signal_track$range[2]),
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = signal_track)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = c(chromstart, chromend), yscale = c(signal_track$range[1], signal_track$range[2]),
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

  if (nrow(signal) >= 2){

    signal_grobs(signal = signal, signaltrack = signal_track, transparency = transparency)

  } else {

    ## just making a flat line
    signalGrob <- segmentsGrob(x0 = 0, y0 = 0, x1 = 1, y1 = 0, gp = gpar(col = signal_track$linecolor, lwd = signal_track$lwd))
    ## Add grob to gtree
    assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = signalGrob), envir = bbEnv)
    warning("Not enough data within range to plot.")

  }


  # ======================================================================================================================================================================================
  # SCALE
  # ======================================================================================================================================================================================

  ## Add y-axis scale
  if (yaxis == TRUE){

    # scaleGrob <- textGrob(label = paste(min(signal[,1]), max(signal[,1]), sep = "-"), x = 1, y = 1,
    #                       just = c("right", "top"), gp = gpar(col = "grey"))
    scaleGrob <- segmentsGrob(x0 = 0, x1 = 0, y0 = 0, y1 = 1, gp = gpar(col = "grey"))
    scalelabelGrob <- textGrob(label = signal_track$range[2], just = c("left", "top"), x = 0, y = 1, gp = gpar(col = "grey"))

    ## Add grob to gtree
    assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = scaleGrob), envir = bbEnv)
    assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = scalelabelGrob), envir = bbEnv)

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

