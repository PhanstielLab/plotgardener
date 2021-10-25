#' Plot multiple signal tracks inline with each other
#' 
#' @usage plotMultiSignal(
#'    data, 
#'    x = NULL, 
#'    y = NULL,
#'    orientation = "h", 
#'    width = NULL,
#'    height = NULL, 
#'    gapdistance = .2,
#'    just = c("left", "top"),
#'    default.units = "inches", 
#'    binSize = NA, 
#'    binCap = TRUE, 
#'    negData = FALSE,
#'    chrom, 
#'    chromstart = NULL, 
#'    chromend = NULL,
#'    assembly = "hg38", 
#'    linecolor,
#'    fill = NA, 
#'    ymax = 1, 
#'    range = NULL, 
#'    scale = FALSE,
#'    bg = NA, 
#'    baseline = TRUE, 
#'    baseline.color = "grey",
#'    baseline.lwd = 1, 
#'    draw = TRUE,
#'    params = NULL, ...
#' )
#'
#' @param data List of data to be plotted as a character value specifying a
#' bigwig file path, a dataframe in BED format, or a
#' \link[GenomicRanges]{GRanges} object with metadata column \code{score}.
#' Either one \code{data} argument or a list of two can be provided, where
#' the second \code{data} will be plotted below the x-axis.
#' @param x A numeric vector or unit object specifying the starting
#' vertex x-locations.
#' @param y A numeric vector, unit object, or a character vector
#' of values containing a "b" combined with a numeric value specifying
#' polygon vertex y-locations.
#' The character vector will place polygon vertex y-locations relative
#' to the bottom of the most recently plotted plot according
#' to the units of the plotgardener page.
#' @param orientation A string specifying signal track orientation.
#' Default value is \code{orientation = "h"}. Options are:
#' \itemize{
#' \item{\code{"v"}: }{Vertical signal track orientation.}
#' \item{\code{"h"}: }{Horizontal signal track orientation.}
#' }
#' @param width A numeric or unit object specifying signal plot width.
#' @param height A numeric or unit object specifying signal plot height.
#' @param gapdistance A numeric object specifying space between plots
#' @param just Justification of signal plot relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal justification
#' and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use
#' if \code{x} or \code{y} are only given as numeric vectors.
#' Default value is \code{default.units = "inches"}.
#' #' @param binSize A numeric specifying the length of each data
#' bin in basepairs. Default value is \code{binSize = NA}.
#' @param binCap A logical value indicating whether the function will
#' limit the number of data bins to 8,000.
#' Default value is \code{binCap = TRUE}.
#' @param negData A logical value indicating whether the data has both
#' positive and negative scores and the y-axis should be split.
#' Default value is \code{negData = FALSE}.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param linecolor A character value or vector of character values specifying the
#' line color(s) outlining the signal track(s).
#' Default value is \code{linecolor = "#37a7db"}.
#' @param fill A character value or vector of length 2 specifying
#' the fill color(s) of the signal track(s). Default value is \code{fill = NA}.
#' @param ymax A numeric specifying the fraction of the max y-value
#' to set as the height of the plot. Default value is \code{ymax = 1}.
#' @param range A numeric vector of length 2 specifying the y-range
#' of data to plot (c(min, max)).
#' @param scale A logical value indicating whether to include a data
#' scale label in the top left corner of the plot.
#' Default value is \code{scale = FALSE}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a
#' baseline along the x-axis. Default value is \code{baseline = TRUE}.
#' @param baseline.color Baseline color.
#' Default value is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width.
#' Default value is \code{baseline.lwd = 1}.
#' @param draw A logical value indicating whether graphics output should be
#' produced. Default value \code{draw = TRUE}.
#' @param params An optional \link[plotgardener]{pgParams} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns several \code{polygon} object containing relevant
#' placement and \link[grid]{grob} information.
#'
#' @examples
#' library("plotgardener")
#'library("org.Hs.eg.db")
#'library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#'library("plotgardenerData")
#'data("GM12878_HiC_10kb")
#'data("IMR90_HiC_10kb")
#'data("GM12878_ChIP_CTCF_signal")
#'data("IMR90_ChIP_CTCF_signal")
#'data("GM12878_ChIP_H3K27ac_signal")
#'data("IMR90_ChIP_H3K27ac_signal")
#'
#'
#'## Load libraries and datasets
#'testList <- list(GM12878_ChIP_CTCF_signal, GM12878_ChIP_H3K27ac_signal, IMR90_ChIP_CTCF_signal, IMR90_ChIP_H3K27ac_signal)
#'## Set genomic and dimension parameters in a `params` object
#'params <- pgParams(chrom = "chr21", chromstart = 28150000, chromend = 29150000, 
#'                     assembly = "hg19", x = 3.5, width = 1.5, default.units = "inches")
#'
#'pageCreate(width = 7, height = 6, default.units = "inches")
#'
#'plotMultiSignal(testList, chrom = "chr21", assembly = "hg19", params = c(params_c, ctcf_range),
#'                fill = "#253494",y = 0.2, height = 4,  x = 0.2, width = 2, gapdistance = 0.1, orientation = "h")


#' @export
plotMultiSignal<- function(data, 
                           binSize = NA, 
                           binCap = TRUE, 
                           negData = FALSE,
                           chrom, 
                           chromstart = NULL, 
                           chromend = NULL,
                           assembly = "hg38", 
                           linecolor= "#37a7db",
                           fill = NA,  
                           ymax = 1, 
                           range = NULL, 
                           scale = FALSE,
                           bg = NA, 
                           baseline = TRUE, 
                           baseline.color = "grey",
                           baseline.lwd = 1, 
                           orientation = "h",
                           x = NULL, 
                           y = NULL, 
                           width = NULL,
                           height = NULL, 
                           just = c("left", "top"),
                           default.units = "inches", gapdistance = .2,
                           draw = TRUE,
                           params = NULL, ...) {
  
  # =========================================================================
  # PARSE PARAMETERS
  # =========================================================================
  
  multisigInternal <- parseParams(
    params = params,
    defaultArgs = formals(eval(match.call()[[1]])),
    declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
    class = "multisigInternal"
  )
  
  ## Set gp
  multisigInternal$gp <- setGP(
    gpList = gpar(),
    params = multisigInternal, ...
  )
  
  ## Justification
  multisigInternal$just <- justConversion(just = multisigInternal$just)
  
  # =========================================================================
  # PARSE ASSEMBLY
  # =========================================================================
  
  multisigInternal$assembly <- parseAssembly(assembly = 
                                               multisigInternal$assembly)
  
  # =========================================================================
  # READ IN FILES/DATAFRAMES
  # =========================================================================
  
  data <- lapply(multisigInternal$data, read_rangeData,
                 assembly = multisigInternal$assembly,
                 chrom = multisigInternal$chrom,
                 start = multisigInternal$chromstart,
                 end = multisigInternal$chromend)
  
  # =========================================================================
  # SET RANGE, PLACEMENT COORDINATES, AND COLORS
  # =========================================================================
  nTracks <- length(multisigInternal$data)
  range <- findSignalRange(multisigInternal$data,negData = negData )
  
  xList <- getXCoordinates(x = multisigInternal$x, 
                           nTracks = nTracks, 
                           width = multisigInternal$width, 
                           orientation = multisigInternal$orientation, 
                           gapdistance = multisigInternal$gapdistance)
  
  yList <- getYCoordinates(y = multisigInternal$y,  
                           nTracks = nTracks, 
                           height = multisigInternal$height, 
                           orientation = multisigInternal$orientation, 
                           gapdistance = multisigInternal$gapdistance)
  
  if (multisigInternal$orientation == "h"){
    height <- (multisigInternal$height - 
                 (multisigInternal$gapdistance * (nTracks-1)))/nTracks
  } else if (multisigInternal$orientation == "v"){
    width <- (multisigInternal$width - 
                (multisigInternal$gapdistance * (nTracks-1)))/nTracks
  }
  else{
    stop("argument \" orientation\" is missing, ","with no default.", call. = FALSE)
  }
  
  linecolorList <- setColors(multisigInternal$linecolor, nTracks)
  fillList <- setColors(multisigInternal$fill, nTracks)
  
  # =========================================================================
  # CALL PLOTSIGNAL
  # =========================================================================
  
  ## Save to a group object?
  pmap(list(multisigInternal$data, xList, yList, linecolorList, fillList), 
       \(d, x, y, l, f){
         plotSignal(data = d, x = x, y = y, 
                    linecolor = l,
                    height = height, 
                    width = width,
                    range = range, 
                    orientation = multisigInternal$orientation, 
                    scale = multisigInternal$scale,
                    binSize = multisigInternal$binSize, 
                    binCap = multisigInternal$binCap, 
                    negData = multisigInternal$negData,
                    chrom = multisigInternal$chrom, 
                    chromstart = multisigInternal$chromstart, 
                    chromend = multisigInternal$chromend,
                    assembly = multisigInternal$assembly, 
                    fill = f, 
                    ymax = multisigInternal$ymax, 
                    bg = multisigInternal$bg, 
                    baseline = multisigInternal$baseline,
                    baseline.color = multisigInternal$baseline.color, 
                    baseline.lwd = multisigInternal$baseline.lwd, 
                    just = multisigInternal$just,
                    default.units = multisigInternal$default.units, 
                    draw = multisigInternal$draw)
       })
}  


# =========================================================================
# Helper functions
# =========================================================================

# =========================================================================
# Finding the range of data processed with read_rangeData()
# =========================================================================

findSignalRange <- function(data, negData = negData){
  #calculating the max from the processed data scores
  rangeMax<- lapply(data, dplyr::select, "score") %>%
    lapply(max) %>%
    unlist() %>%
    max
  #if no negative data, lower bound of range = 0
  if(negData == FALSE)
    rangeMin <- 0 
  #if there are negative values, lower bound of range is the min value 
  else{
    rangeMin<- lapply(data, dplyr::select, "score") %>%
      lapply(min) %>%
      unlist() %>%
      min
  }
  return(c(rangeMin,rangeMax))
}

# =========================================================================
# Fill in colors
# =========================================================================

setColors <- function(colorList, nTracks){
  #exit recursive with list of colors of the correct length
  if(length(colorList) == nTracks)
    return(colorList)
  
  #if there are too few colors for the number of tracks, repeat the colors provided
  else if(length(colorList) < nTracks){
    setColors(append(colorList, colorList), nTracks)
  }
  #if there are too many colors for the number of tracks, truncate the color list
  else if(length(colorList) > nTracks)
    return(colorList[1:nTracks])
}

# =========================================================================
# Calculate x coordinates of each track
# =========================================================================

getXCoordinates<- function(x, nTracks, width, orientation, gapdistance){
  #create a list for the x coordinates
  xList<- rep(x, nTracks)
  #if the orientation is h, the x coordinates stay the same, so return a list of the same x
  if(orientation == "h"){
   return(xList)
  } 
  
  else if(orientation == "v"){
    #calculating the width of each plot 
    internalWidth<- (width - (gapdistance * (nTracks-1)))/nTracks
    # first plot won't change x, others are shifting to the right 
    for (i in 2:nTracks)
      xList[i]<- (xList[i-1] + internalWidth + gapdistance)
    return(xList)
  } 
  
  else{
    stop("argument \" orientation\" is missing, ","with no default.", call. = FALSE)
  }
}


# =========================================================================
# Calculate y coordinates of each track
# =========================================================================

getYCoordinates<- function(y, nTracks, height, orientation, gapdistance){
  #create a list for the x coordinates
  yList<- rep(y, nTracks)
  #if the orientation is v, the y coordinates stay the same, so return a list of the same x
  if(orientation == "v"){
    return(yList)
  }
  
  else if(orientation == "h"){
    #calculating the height of each plot
    internalHeight<- (height - (gapdistance * (nTracks-1)))/nTracks
    # first plot won't change y, others are shifting to the right
    for (i in 2:nTracks)
      yList[i]<- (yList[i-1] + internalHeight + gapdistance)
    return(yList)
  } 
  
  else{
    stop("argument \" orientation\" is missing, ","with no default.", call. = FALSE)
  }
}

