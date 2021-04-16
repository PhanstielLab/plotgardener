#' Plot any kind of signal track data for a single chromosome
#'
#' @param data Data to be plotted as a character value specifying a bigwig file path, a dataframe in BED format, or a \link[GenomicRanges]{GRanges} object with metadata column \code{counts}.
#' Either one \code{data} argument or a list of two can be provided, where the second \code{data} will be plotted below the x-axis.
#' @param binSize A numeric specifying the length of each data bin in basepairs. Default value is \code{binSize = NA}.
#' @param binCap A logical value indicating whether the function will limit the number of data bins to 8,000. Default value is \code{binCap = TRUE}.
#' @param negData A logical value indicating whether the data has both positive and negative scores and the y-axis should be split. Default value is \code{negData = FALSE}.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a \link[BentoBox]{bb_assembly} object. Default value is \code{assembly = "hg19"}.
#' @param linecolor A character value or vector of length 2 specifying the line color(s) outlining the signal track(s). Default value is \code{linecolor = "#37a7db"}.
#' @param fill A character value or vector of length 2 specifying the fill color(s) of the signal track(s). Default value is \code{fill = NA}.
#' @param ymax A numeric specifying the fraction of the max y-value to set as the height of the plot. Default value is \code{ymax = 1}.
#' @param range A numeric vector of length 2 specifying the y-range of data to plot (c(min, max)).
#' @param scale A logical value indicating whether to include a data scale label in the top left corner of the plot. Default value is \code{scale = FALSE}.
#' @param bg Character value indicating background color. Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a baseline along the x-axis. Default value is \code{baseline = TRUE}.
#' @param baseline.color Baseline color. Default value is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width. Default value is \code{baseline.lwd = 1}.
#' @param x A numeric or unit object specifying signal plot x-location.
#' @param y A numeric, unit object, or character containing a "b" combined with a numeric value specifying signal plot y-location. The character value will
#' place the signal plot y relative to the bottom of the most recently plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying signal plot width.
#' @param height A numeric or unit object specifying signal plot height.
#' @param just Justification of signal plot relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{x}, \code{y}, \code{width}, or \code{height} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be produced. Default value \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_signal} object containing relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load signal data
#' data("bb_imrH3K27acData")
#' data("bb_gmH3K27acData")
#'
#' ## Create a page
#' bb_pageCreate(width = 7.5, height = 2.1, default.units = "inches")
#'
#' ## Define region
#' region <- bb_params(chrom = "chr21", chromstart = 28000000, chromend = 30300000,
#'                     range = c(0, 45))
#'
#' ## Plot and place signal plots
#' signal1 <- bb_plotSignal(data = bb_imrH3K27acData, params = region,
#'                          x = 0.5, y = 0.25, width = 6.5, height = 0.65,
#'                          just = c("left", "top"), default.units = "inches")
#'
#' signal2 <- bb_plotSignal(data = bb_gmH3K27acData, params = region,
#'                          linecolor = "#7ecdbb",
#'                          x = 0.5, y = 1, width = 6.5, height = 0.65,
#'                          just = c("left", "top"), default.units = "inches")
#'
#' ## Plot genome label
#' bb_plotGenomeLabel(chrom = "chr21", chromstart = 28000000, chromend = 30300000,
#'                    x = 0.5, y = 1.68, length = 6.5, default.units = "inches")
#'
#' ## Add text labels
#' bb_plotText(label = "IMR90", fonsize = 10, fontcolor = "#37a7db",
#'             x = 0.5, y = 0.25, just = c("left", "top"), default.units = "inches")
#' bb_plotText(label = "GM12878", fonsize = 10, fontcolor = "#7ecdbb",
#'             x = 0.5, y = 1, just = c("left", "top"), default.units = "inches")
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @details
#' #A signal track can be placed on a BentoBox coordinate page by providing plot placement parameters:
#' \preformatted{
#' bb_plotSignal(data, chrom,
#'               chromstart = NULL, chromend = NULL,
#'               x, y, width, height, just = c("left", "top"),
#'               default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated signal track by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotSignal(data, chrom,
#'               chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
bb_plotSignal <- function(data, binSize = NA, binCap = TRUE, negData = FALSE, chrom, chromstart = NULL, chromend = NULL, assembly = "hg19", linecolor = "#37a7db", fill = NA,
                          ymax = 1, range = NULL, scale = FALSE, bg = NA, baseline = TRUE, baseline.color = "grey", baseline.lwd = 1, x = NULL, y = NULL, width = NULL, height = NULL,
                          just = c("left", "top"), default.units = "inches", draw = TRUE, params = NULL, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that catches errors for bb_plotSignal
  errorcheck_bb_signaltrack <- function(signal, signaltrack){

    dfChecks <- function(signal){

      if (!"data.frame" %in% class(signal)){

        if (!"GRanges" %in% class(signal)){

          if (!file_ext(signal) %in% c("bw", "bigWig", "bigwig", "bedgraph")){

            stop("Invalid input. File must have a valid bigwig or bedgraph extension.", call. = FALSE)

          }
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

      if (signal_track$chromstart == signal_track$chromend){
        stop("Genomic region is 0 bp long.", call. = FALSE)
      }

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

      if (!"GRanges" %in% class(signal)){
        signal <- bb_readBigwig(file = signal, chrom = signaltrack$chrom, chromstart = signaltrack$chromstart, chromend = signaltrack$chromend)
      }

    }

    signal <- as.data.frame(signal)

    ## Check for overlapping data ranges
    if (any(IRanges::overlapsAny(GenomicRanges::makeGRangesFromDataFrame(signal), drop.self = TRUE) == TRUE)){
      stop("Data ranges cannot overlap. Please check `start` and `end` column ranges.", call. = FALSE)
    }

    return(signal)

  }

  ## Define a function that formats/filters signal data
  format_data <- function(signal, signaltrack){

    if (!any(colnames(signal) == "score")){
      stop("Cannot find associated `score` column in data.", call. = FALSE)
    } else {
      ## Grab chrom, start, end and a "score" column
      scoreCol <- which(colnames(signal) == "score")
      signal <- signal[c(1:3, scoreCol)]

    }

    ## Ensure the chromosome is a character
    signal[,1] <- as.character(signal[,1])

    ## Filter for desired region
    signal <- signal[which(signal[,1] == signaltrack$chrom & ((signal[,2] > signaltrack$chromstart & signal[,2] < signaltrack$chromend | signal[,3] > signaltrack$chromstart &
                                                                 signal[,3] < signaltrack$chromend ))), (2:4)]
    ## Remove any duplicate rows
    signal <- signal[!duplicated(signal),]

    ## Remove any NaN score values
    signal <- na.omit(signal)

    return(signal)

  }

  ## Define a function that checks and adjust the number/sizes of bins
  check_binNum <- function(signaltrack, binCap){

    if (!is.na(signaltrack$binSize)){

      if (signaltrack$binSize %% 0.25 != 0){
        updated_binSize <- round(signaltrack$binSize/0.25) * 0.25
        signaltrack$binSize <- updated_binSize
      }


      binNum = (signaltrack$chromend - signaltrack$chromstart)/signaltrack$binSize
      signaltrack$binNum <- binNum

      if (!is.nan(binNum)){

        ## Scale back binNum and print warning if binNum is greater than 8000
        if (binNum > 8000 && binCap == TRUE){
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
          signaltrack$binSize <- updated_binSize
          warning(paste0("Number of bins larger than plot length: adjusting to ", binNum, " bins of size 1."), call. = FALSE)
        }


      } else {
        signaltrack$binSize <- NA
      }


    }

    return(signaltrack)
  }

  ## Define a function that bins, links, sorts, and combines data
  parseData <- function(signal, signaltrack){

    if (!is.na(signaltrack$binSize)){

      # ===============================================================================================================================================
      # BIN DATA
      # ===============================================================================================================================================
      ## Find the max signal value for each bin

      binChromend <- signaltrack$binSize*signaltrack$binNum + signaltrack$chromstart + signaltrack$binSize

      binDF <- data.frame("start" = seq(signaltrack$chromstart, binChromend - signaltrack$binSize, signaltrack$binSize),
                          "end" = seq(signaltrack$chromstart + signaltrack$binSize, binChromend, signaltrack$binSize))
      binDF$maxScore <- rebinBigwig(signal, binDF)

      ## Use binned data as new signal data
      newSignal <- binDF
      colnames(newSignal) <- c("chromstart", "chromend", "counts")

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

      signal <- newSignal

    }


    return(signal)
  }

  ## Define a function that adjusts the range
  set_range <- function(signal1, signal2, signaltrack, split, pos = TRUE){

    if (split == TRUE){

      ## posSignal
      if (pos == TRUE){

        if (is.null(signaltrack$range)){

          if (nrow(signal1) >= 2){

            signaltrack$range[2] <- signaltrack$ymax * max(signal2[,2])

          } else {

            signaltrack$range[2] <- 1

          }

        }

      }

      ## negSignal
      if (pos == FALSE){

        if (is.na(signaltrack$range[1])){

          if (nrow(signal1) >= 2){

            signaltrack$range[1] <- signaltrack$ymax * min(signal2[,2])

          } else {

            signaltrack$range[1] <- -1

          }

        }

      }


    } else {

      if (is.null(signaltrack$range)){

        if (nrow(signal1) >= 2){

          signaltrack$range <- c(0, signaltrack$ymax*max(signal2[,2]))

        } else {

          signaltrack$range <- c(0, 1)

        }

      }

    }

    return(signaltrack)

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
  sigGrob <- function(signal, fillCol, lineCol, gp){

    gp$col <- lineCol
    if (!is.null(fillCol) & !is.na(fillCol)){

      if ("alpha" %in% names(gp)){
        fillCol <- makeTransparent(color = fillCol, alpha = gp$alpha)
      }

      gp$fill <- fillCol

      sigGrob <- polygonGrob(x = c(signal[1,1], signal[,1], signal[nrow(signal),1]),
                             y = c(0, signal[,2], 0), gp = gp, default.units = "native")

    } else {

      sigGrob <- segmentsGrob(x0 = signal[c(1:1-length(signal[,1])), 1], y0 = signal[c(1:1-length(signal[,2])), 2],
                              x1 = signal[c(2:length(signal[,1])), 1], y1 = signal[c(2:length(signal[,2])), 2], gp = gp, default.units = "native")

    }


    ## Add grob to gtree
    assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = sigGrob), envir = bbEnv)

  }

  ## Define a function that finds data that falls out of the plot's range and draws a tiny black line to indicate it
  cutoffGrobs <- function(signal, signaltrack, side){

    grobCutoffs <- function(df, side){

      x0 <- df[1]
      x1 <- df[2]

      if (side == "top"){
        y <- 1
      } else {
        y <- 0
      }

      cutoffGrob <- segmentsGrob(x0 = x0, x1 = x1, y0 = unit(y, "npc"), y1 = unit(y, "npc"),
                                 gp = gpar(lwd = 1, col = "grey"),
                                 default.units = "native")
      assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = cutoffGrob), envir = bbEnv)

    }

    if (side == "top"){
      outsideData <- which(signal[,2] > signaltrack$range[2])
    } else {
      outsideData <- which(signal[,2] < signaltrack$range[1])
    }

    ## Get index pairs of outside data
    outsidePairs <- as.integer(outsideData + 1)
    ## Combine and order
    Outsidei <- c(outsideData, outsidePairs)
    Outsidei <- Outsidei[order(Outsidei)]
    ## Get xcoords
    signal <- as.data.frame(signal)
    Outside <- signal[Outsidei,1]
    ## x0s are odd indeces and x1s are even
    x0s <- Outside[c(TRUE, FALSE)]
    x1s <- Outside[c(FALSE, TRUE)]
    pairs <- data.frame("x0" = x0s, "x1" = x1s)

    invisible(apply(pairs, 1, grobCutoffs, side = side))

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(fill)) fill <- NULL
  if(missing(binSize)) binSize <- NULL
  if(missing(binCap)) binCap <- NULL
  if(missing(negData)) negData <- NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(ymax)) ymax <- NULL
  if(missing(scale)) scale <- NULL
  if(missing(bg)) bg <- NULL
  if(missing(baseline)) baseline <- NULL
  if(missing(baseline.color)) baseline.color <- NULL
  if(missing(baseline.lwd)) baseline.lwd <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if data/chrom arguments are missing (could be in object)
  if(!hasArg(data)) data <- NULL
  if(!hasArg(chrom)) chrom <- NULL

  ## Compile all parameters into an internal object
  bb_sigInternal <- structure(list(data = data, chrom = chrom, chromstart = chromstart, chromend = chromend, negData = negData, range = range, linecolor = linecolor, binSize = binSize,
                                   binCap = binCap, fill = fill, assembly = assembly, ymax = ymax, scale = scale, bg = bg, baseline = baseline,
                                   width = width, height = height, baseline.color = baseline.color, baseline.lwd = baseline.lwd,
                                   x = x, y = y, just = just, default.units = default.units, draw = draw, gp = gpar()), class = "bb_sigInternal")

  bb_sigInternal <- parseParams(bb_params = params, object_params = bb_sigInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_sigInternal$linecolor)) bb_sigInternal$linecolor <-"#37a7db"
  if(is.null(bb_sigInternal$fill)) bb_sigInternal$fill <- NA
  if(is.null(bb_sigInternal$binSize)) bb_sigInternal$binSize <- NA
  if(is.null(bb_sigInternal$binCap)) bb_sigInternal$binCap <- TRUE
  if(is.null(bb_sigInternal$negData)) bb_sigInternal$negData <- FALSE
  if(is.null(bb_sigInternal$assembly)) bb_sigInternal$assembly <- "hg19"
  if(is.null(bb_sigInternal$ymax)) bb_sigInternal$ymax <- 1
  if(is.null(bb_sigInternal$scale)) bb_sigInternal$scale <- FALSE
  if(is.null(bb_sigInternal$bg)) bb_sigInternal$bg <- NA
  if(is.null(bb_sigInternal$baseline)) bb_sigInternal$baseline <- TRUE
  if(is.null(bb_sigInternal$baseline.color)) bb_sigInternal$baseline.color <- "grey"
  if(is.null(bb_sigInternal$baseline.lwd)) bb_sigInternal$baseline.lwd <- 1
  if(is.null(bb_sigInternal$just)) bb_sigInternal$just <- c("left", "top")
  if(is.null(bb_sigInternal$default.units)) bb_sigInternal$default.units <- "inches"
  if(is.null(bb_sigInternal$draw)) bb_sigInternal$draw <- TRUE

  ## Set gp
  bb_sigInternal$gp <- setGP(gpList = bb_sigInternal$gp, params = bb_sigInternal, ...)

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  signal_track <- structure(list(chrom = bb_sigInternal$chrom, chromstart = bb_sigInternal$chromstart, chromend = bb_sigInternal$chromend, assembly = bb_sigInternal$assembly,
                                 binSize = bb_sigInternal$binSize, binNum = NULL, range = bb_sigInternal$range, ymax = bb_sigInternal$ymax,
                                 x = bb_sigInternal$x, y = bb_sigInternal$y, width = bb_sigInternal$width, height = bb_sigInternal$height,
                                 just = bb_sigInternal$just, grobs = NULL), class = "bb_signal")
  attr(x = signal_track, which = "plotted") <- bb_sigInternal$draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_sigInternal$data)) stop("argument \"data\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_sigInternal$chrom)) stop("argument \"chrom\" is missing, with no default.", call. = FALSE)

  check_placement(object = signal_track)
  errorcheck_bb_signaltrack(signal = bb_sigInternal$data, signaltrack = signal_track)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  signal_track$assembly <- parse_bbAssembly(assembly = signal_track$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  signal_track <- defaultUnits(object = signal_track, default.units = bb_sigInternal$default.units)

  # ======================================================================================================================================================================================
  # WHOLE CHROM
  # ======================================================================================================================================================================================

  if (is.null(signal_track$chromstart) & is.null(signal_track$chromend)){

    if (class(signal_track$assembly$TxDb) == "TxDb"){
      txdbChecks <- TRUE
    } else {
      txdbChecks <- check_loadedPackage(package = signal_track$assembly$TxDb, message = paste(paste0("`", signal_track$assembly$TxDb,"`"),
                                                                                              "not loaded. Please install and load to generate full chromosome signal track."))
    }

    xscale <- c(0, 1)
    if (txdbChecks == TRUE){

      if (class(signal_track$assembly$TxDb) == "TxDb"){
        tx_db <- signal_track$assembly$TxDb
      } else {
        tx_db <- eval(parse(text = signal_track$assembly$TxDb))
      }

      assembly_data <- seqlengths(tx_db)
      if (!signal_track$chrom %in% names(assembly_data)){
        warning(paste("Chromosome", paste0("'", signal_track$chrom, "'"), "not found in", paste0("`", signal_track$assembly$TxDb, "`"), "and data for entire chromosome cannot be plotted."), call. = FALSE)
      } else {
        signal_track$chromstart <- 1
        signal_track$chromend <- assembly_data[[signal_track$chrom]]
        xscale <- c(signal_track$chromstart, signal_track$chromend)

      }


    }

  } else {
    xscale <- c(signal_track$chromstart, signal_track$chromend)
  }

  # ======================================================================================================================================================================================
  # SET BINSIZE
  # ======================================================================================================================================================================================

  if (is.na(signal_track$binSize) == TRUE){

    if (!is.null(signal_track$chromstart)){
      binSize <- (signal_track$chromend - signal_track$chromstart)/2000
      signal_track$binSize <- binSize
    }

  }

  # ======================================================================================================================================================================================
  # CHECK AND ADJUST BIN NUMBER/BIN SIZE
  # ======================================================================================================================================================================================

  signal_track <- check_binNum(signaltrack = signal_track, binCap = bb_sigInternal$binCap)

  # ======================================================================================================================================================================================
  # READ IN, FORMAT, FILTER, BIN, LINK AND SORT DATA
  # ======================================================================================================================================================================================

  if (length(class(bb_sigInternal$data)) == 1 && class(bb_sigInternal$data) == "list"){

    signal <- lapply(bb_sigInternal$data, read_signal, signaltrack = signal_track)
    signal <- lapply(signal, format_data, signaltrack = signal_track)
    posSignal <- signal[[1]]
    negSignal <- signal[[2]]
    split <- TRUE

  } else {

    signal <- read_signal(signal = bb_sigInternal$data, signaltrack = signal_track)
    signal <- format_data(signal = signal, signaltrack = signal_track)

    if (any(signal[,3] < 0)){

      if (bb_sigInternal$negData == FALSE){
        warning("Negative scores detected in signal data. To make an entirely positive signal track,
                please remove negative scores from data.", call. = FALSE)
      }

      posSignal <- signal[which(signal[,3] >= 0),]
      negSignal <- signal[which(signal[,3] < 0),]
      negSignal[,3] <- negSignal[,3]*-1
      split <- TRUE

    } else {

      posSignal <- signal
      split <- FALSE
      if (bb_sigInternal$negData == TRUE){
        negSignal <- data.frame()
        split <- TRUE
      }

    }

  }


  # ======================================================================================================================================================================================
  # BIN, LINK, AND SORT DATA AND FIX Y-LIMITS
  # ======================================================================================================================================================================================

  if (split == TRUE){

    if (nrow(posSignal) >= 2){
      posSignal2 <- parseData(signal = posSignal, signaltrack = signal_track)
    }

    signal_track <- set_range(signal1 = posSignal, signal2 = posSignal2, signaltrack = signal_track, split = TRUE)

    if (nrow(negSignal) >= 2){
      negSignal2 <- parseData(signal = negSignal, signaltrack = signal_track)
      negSignal2[,2] <- negSignal2[,2]*-1
    }


    signal_track <- set_range(signal1 = negSignal, signal2 = negSignal2, signaltrack = signal_track, split = TRUE, pos = FALSE)


    if (signal_track$range[1] == 0 & signal_track$range[2] == 0){
      signal_track$range <- c(-1, 1)
    }



  } else {

    if (nrow(posSignal) >= 2){
      posSignal2 <- parseData(signal = posSignal, signaltrack = signal_track)
    }

    signal_track <- set_range(signal1 = posSignal, signal2 = posSignal2, signaltrack = signal_track, split = FALSE)


  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_signal", length(grep(pattern = "bb_signal", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(signal_track$x) & is.null(signal_track$y)){

    vp <- viewport(height = unit(0.25, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = xscale,
                   yscale = c(signal_track$range[1], signal_track$range[2]),
                   just = "center",
                   name = vp_name)

    if (bb_sigInternal$draw == TRUE){

      vp$name <- "bb_signal1"
      grid.newpage()

    }

  } else {

    add_bbViewport(vp_name)

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = signal_track)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = xscale,
                   yscale = c(signal_track$range[1], signal_track$range[2]),
                   just = bb_sigInternal$just,
                   name = vp_name)
  }


  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
  # ======================================================================================================================================================================================

  backgroundGrob <- rectGrob(gp = gpar(fill = bb_sigInternal$bg, col = NA), name = "background")
  assign("signal_grobs", gTree(vp = vp, children = gList(backgroundGrob)), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  if (!is.na(signal_track$binSize)){

    if (split == TRUE){

      fills <- parseColors(bb_sigInternal$fill)
      lines <- parseColors(bb_sigInternal$linecolor)

      if (nrow(posSignal) >= 2){
        sigGrob(signal = posSignal2, fillCol = fills[[1]], lineCol = lines[[1]], gp = bb_sigInternal$gp)
        ## Find and make cutoff lines
        cutoffGrobs(signal = posSignal2, signaltrack = signal_track, side = "top")

      } else {
        bb_sigInternal$gp$col <- lines[[1]]
        posGrob <- segmentsGrob(x0 = 0, y0 = unit(0, "native"), x1 = 1, y1 = unit(0, "native"), gp = bb_sigInternal$gp)
        assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = posGrob), envir = bbEnv)
        warning("Not enough top signal data to plot.", call. = FALSE)
      }


      if (nrow(negSignal) >= 2){
        sigGrob(signal = negSignal2, fillCol = fills[[2]], lineCol = lines[[2]], gp = bb_sigInternal$gp)
        ## Find and make cutoff lines
        cutoffGrobs(signal = negSignal2, signaltrack = signal_track, side = "bottom")
      } else {
        bb_sigInternal$gp$col <- lines[[2]]
        negGrob <- segmentsGrob(x0 = 0, y0 = unit(0, "native"), x1 = 1, y1 = unit(0, "native"), gp = bb_sigInternal$gp)
        assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = negGrob), envir = bbEnv)
        warning("Not enough bottom signal data to plot.", call. = FALSE)
      }


      lineGrob <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = 0, y1 = 0,
                               gp = gpar(col = bb_sigInternal$baseline.color, lwd = bb_sigInternal$baseline.lwd),
                               default.units = "native")
      assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = lineGrob), envir = bbEnv)

    } else {

      if (nrow(posSignal) >= 2){

        if(bb_sigInternal$baseline == TRUE){
          baselineGrob <- segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"), y0 = 0, y1 = 0,
                                   gp = gpar(col = bb_sigInternal$baseline.color, lwd = bb_sigInternal$baseline.lwd),
                                   default.units = "native")
          assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = baselineGrob), envir = bbEnv)
        }

        sigGrob(signal = posSignal2, fillCol = bb_sigInternal$fill[1], lineCol = bb_sigInternal$linecolor[1], gp = bb_sigInternal$gp)
        ## Find and make cutoff lines
        cutoffGrobs(signal = posSignal2, signaltrack = signal_track, side = "top")
      } else {
        bb_sigInternal$gp$col <- bb_sigInternal$linecolor
        signalGrob <- segmentsGrob(x0 = 0, y0 = unit(0, "native"), x1 = 1, y1 = unit(0, "native"), gp = bb_sigInternal$gp)
        assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = signalGrob), envir = bbEnv)
        warning("Not enough data within range to plot.", call. = FALSE)
      }

    }
    # ======================================================================================================================================================================================
    # SCALE
    # ======================================================================================================================================================================================

    ## Add scale of the range of data in the top left corner
    if (bb_sigInternal$scale == TRUE){
      upperLim <- round(signal_track$range[2], digits = 4)
      lowerLim <- round(signal_track$range[1], digits = 4)
      scaleGrob <- textGrob(label = paste0("[", lowerLim, " - ", upperLim, "]"), just = c("left", "top"), x = 0, y = 1,
                            gp = bb_sigInternal$gp)

      ## Add grob to gtree
      assign("signal_grobs", addGrob(gTree = get("signal_grobs", envir = bbEnv), child = scaleGrob), envir = bbEnv)

    }


  }


  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_sigInternal$draw == TRUE){

    grid.draw(get("signal_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  signal_track$grobs <-  get("signal_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_signal[", vp$name, "]"))
  invisible(signal_track)

}

