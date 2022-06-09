#' Plot multiple signal tracks in line with each other
#' 
#' @usage plotMultiSignal(
#'    data,
#'    binSize = NA, 
#'    binCap = TRUE, 
#'    negData = FALSE,
#'    chrom, 
#'    chromstart = NULL, 
#'    chromend = NULL,
#'    assembly = "hg38", 
#'    linecolor= "#37a7db",
#'    fill = NA, 
#'    ymax = 1, 
#'    range = NULL, 
#'    scale = FALSE,
#'    label = NULL,
#'    bg = NA, 
#'    baseline = TRUE, 
#'    baseline.color = "grey",
#'    baseline.lwd = 1, 
#'    orientation = "h", 
#'    x = NULL, 
#'    y = NULL,
#'    width = NULL,
#'    height = NULL,
#'    just = c("left", "top"),
#'    gapdistance = .2,
#'    default.units = "inches", 
#'    draw = TRUE,
#'    params = NULL, ...
#' )
#'
#' @param data List of data to be plotted as character values specifying 
#' multiple bigwig file paths, dataframes in BED format, or
#' \link[GenomicRanges]{GRanges} objects with metadata column \code{score}.
#' @param binSize A numeric specifying the length of each data
#' bin in basepairs. Default value is \code{binSize = NA}.
#' @param binCap A logical value indicating whether the function will
#' limit the number of data bins to 8,000.
#' Default value is \code{binCap = TRUE}.
#' @param negData A logical value indicating whether any of the data has both
#' positive and negative scores and the y-axis of each signal track
#' should be split. Default value is \code{negData = FALSE}.
#' @param chrom Chromosome of region to be plotted, as a string.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param linecolor A character value or vector of character values specifying 
#' the line color(s) outlining the signal tracks.
#' Default value is \code{linecolor = "#37a7db"}.
#' @param fill A character value or vector specifying
#' the fill color(s) of the signal tracks. Default value is \code{fill = NA}.
#' @param ymax A numeric specifying the fraction of the max y-value
#' to set as the height of each plot. Default value is \code{ymax = 1}.
#' @param range A numeric vector of length 2 specifying the y-range
#' of data to plot (c(min, max)) in each signal track. If \code{range = NULL},
#' an optimal range for all signal tracks will be calculated.
#' @param scale A logical value indicating whether to include a data
#' scale label in the top left corner of each plot.
#' Default value is \code{scale = FALSE}.
#' @param label An optional character vector to conveniently add text labels
#' to signal tracks. If \code{scale = TRUE}, the labels will be drawn in the 
#' top right of the signal tracks. Otherwise, the label will be drawn in the 
#' top left of the plot. For more customizable labels, 
#' use \link[plotgardener]{plotText}.
#' @param bg Character value indicating background color.
#' Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a
#' baseline along the x-axis. Default value is \code{baseline = TRUE}.
#' @param baseline.color Baseline color.
#' Default value is \code{baseline.color = "grey"}.
#' @param baseline.lwd Baseline line width.
#' Default value is \code{baseline.lwd = 1}.
#' @param orientation A string specifying signal track orientations.
#' Default value is \code{orientation = "h"}. Options are:
#' \itemize{
#' \item{\code{"v"}: }{Vertical signal track orientations, where signal tracks
#' will be stacked from left to right.}
#' \item{\code{"h"}: }{Horizontal signal track orientations, where signal tracks
#' will be stacked from top to bottom.}
#' }
#' @param x A numeric vector or unit object specifying the overall multisignal 
#' x-location.
#' @param y A numeric, unit object, or character containing a "b"
#' combined with a numeric value specifying overall multisignal plot y-location.
#' The character value will
#' place the multisignal plot y relative to the bottom of the most recently
#' plotted plot according to the units of the plotgardener page.
#' @param width A numeric or unit object specifying overall multisignal plot 
#' width.
#' @param height A numeric or unit object specifying overall multisignal plot 
#' height.
#' @param just Justification of overall multisignal plot relative to 
#' its (x, y) location. If there are two values, the first value specifies 
#' horizontal justification and the second value specifies vertical 
#' justification. Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = c("left", "top")}.
#' @param gapdistance A numeric or unit object 
#' specifying space between plots. Default value is \code{gapdistance = 0.2}.
#' @param default.units A string indicating the default units to use
#' if \code{x} or \code{y} are only given as numerics.
#' Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be
#' produced. Default value \code{draw = TRUE}.
#' @param params An optional \link[plotgardener]{pgParams} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a list of \code{signal} objects containing relevant
#' genomic region, placement, and \link[grid]{grob} information for each signal
#' track.
#'
#' @examples
#' library("plotgardenerData")
#' data("GM12878_ChIP_CTCF_signal")
#' data("IMR90_ChIP_CTCF_signal")
#' data("GM12878_ChIP_H3K27ac_signal")
#' data("IMR90_ChIP_H3K27ac_signal")
#'
#' ## List of multiple signal datasets
#' signalList <- list(GM12878_ChIP_CTCF_signal, GM12878_ChIP_H3K27ac_signal, 
#'     IMR90_ChIP_CTCF_signal, IMR90_ChIP_H3K27ac_signal)
#' 
#' ## Create page
#' pageCreate(width = 6.9, height = 3.5, default.units = "inches")
#'
#' ## Plot multiple signals
#' multisignal <- plotMultiSignal(signalList, chrom = "chr21",
#'                                chromstart = 28150000, chromend = 29150000,
#'                                linecolor = c(brewer.pal(n = 9,"YlGnBu")[4],
#'                                              brewer.pal(n = 9,"YlGnBu")[5],
#'                                              brewer.pal(n = 9,"YlGnBu")[6],
#'                                              brewer.pal(n = 9,"YlGnBu")[7]),
#'                                label = c("GM12878 CTCF", "GM12878 H3K27ac",
#'                                          "IMR90 CTCF", "IMR90 H3K27ac"),
#'                                assembly = "hg19",
#'                                x = 0.2, y = 0.2,
#'                                width = 6.5, height = 3,
#'                                default.units = "inches",
#'                                gapdistance = 0.1)
#'                                
#' ## Plot genome label
#' plotGenomeLabel(
#'     chrom = "chr21",
#'     chromstart = 28150000, chromend = 29150000,
#'     assembly = "hg19",
#'     scale = "Kb",
#'     x = 0.2, y = 3.25, length = 6.5,
#'     default.units = "inches"
#' )                                
#' 
#' ## Hide page guides
#' pageGuideHide()              
#' @export
plotMultiSignal<- function(data, 
                           binSize = NA, 
                           binCap = TRUE, 
                           negData = FALSE,
                           chrom, 
                           chromstart = NULL, 
                           chromend = NULL,
                           assembly = "hg38", 
                           linecolor = "#37a7db",
                           fill = NA,  
                           ymax = 1, 
                           range = NULL, 
                           scale = FALSE,
                           label = NULL,
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
                           gapdistance = 0.2,
                           default.units = "inches", 
                           draw = TRUE,
                           params = NULL, ...) {
  
    # =========================================================================
    # FUNCTIONS
    # =========================================================================
    
    ## Define a function to fill in colors
    setColors <- function(colorList, nTracks){
        # exit recursive with list of colors of the correct length
        if (length(colorList) == nTracks){
            return(colorList)
        }
            
        # if there are too few colors for the number of tracks, 
        # repeat the colors provided
        else if (length(colorList) < nTracks){
            setColors(append(colorList, colorList), nTracks)
        }
        # if there are too many colors for the number of tracks, 
        # truncate the color list
        else if(length(colorList) > nTracks){
            return(colorList[seq(1, nTracks)])
        }
            
    }
    
    ## Define a function to calculate x coordinates of each track
    getXCoordinates <- function(x, nTracks, width, orientation, gapdistance){
        # create a list for the x coordinates
        xList <- rep(x, nTracks)
        # if the orientation is h, the x coordinates stay the same, so return 
        # a list of the same x
        if (orientation == "h"){
            return(xList)
        } 
        
        else if (orientation == "v"){
            # calculating the width of each plot 
            internalWidth <- (width - (gapdistance * (nTracks-1)))/nTracks
            # first plot won't change x, others are shifting to the right 
            for (i in 2:nTracks){
                xList[i]<- (xList[i-1] + internalWidth + gapdistance)
            }
                
            return(xList)
        } 
        
    }
    
    ## Define a function to calculate y coordinates of each track
    getYCoordinates <- function(y, nTracks, height, orientation, gapdistance){
        # create a list for the y coordinates
        yList <- rep(y, nTracks)
        # if the orientation is v, the y coordinates stay the same, so return 
        # a list of the same y
        if (orientation == "v"){
            return(yList)
        }
        
        else if (orientation == "h"){
            # calculating the height of each plot
            internalHeight <- (height - (gapdistance * (nTracks-1)))/nTracks
            # first plot won't change y, others are shifting down
            for (i in 2:nTracks){
                yList[i]<- (yList[i-1] + internalHeight + gapdistance) 
            }
                
            return(yList)
        } 

    }
    
    ## Define a function to get and add a signal gtree to one multisignal gtree 
    getGrobs <- function(object){
        
        assign("pg_multiSignals",
               addGrob(
                   gTree = get("pg_multiSignals", envir = pgEnv),
                   child = object$grobs
               ),
               envir = pgEnv
        )
        
    }
    
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
    # CATCH ERRORS
    # =========================================================================
  
    if (is.null(multisigInternal$data)) stop("argument \"orientation\" is",
                                     "missing, with no default.", call. = FALSE)
  
    # =========================================================================
    # PARSE ASSEMBLY
    # =========================================================================
  
    multisigInternal$assembly <- parseAssembly(assembly = 
                                               multisigInternal$assembly)
    
    # =========================================================================
    # PARSE UNITS
    # =========================================================================
    
    multisigInternal <- defaultUnits(
        object = multisigInternal,
        default.units = multisigInternal$default.units
    )
  
    multisigInternal$gapdistance <- 
        misc_defaultUnits(value = multisigInternal$gapdistance,
                          name = "gapdistance",
                          default.units = multisigInternal$default.units)
    
    # =========================================================================
    # READ IN FILES/DATAFRAMES AND SET RANGE
    # =========================================================================
  
    signal <- lapply(multisigInternal$data, read_rangeData,
                     assembly = multisigInternal$assembly,
                     chrom = multisigInternal$chrom,
                     start = multisigInternal$chromstart,
                     end = multisigInternal$chromend)
    
    
    range <- calcSignalRange(data = signal,
                             chrom = multisigInternal$chrom,
                             chromstart = multisigInternal$chromstart,
                             chromend = multisigInternal$chromend,
                             assembly = multisigInternal$assembly,
                             negData = multisigInternal$negData)
    
    # =========================================================================
    # SET PLACEMENT COORDINATES, COLORS, AND LABELS
    # =========================================================================
  
    nTracks <- length(multisigInternal$data)

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

    linecolorList <- setColors(multisigInternal$linecolor, nTracks)
    fillList <- setColors(multisigInternal$fill, nTracks)
    
    if (is.null(multisigInternal$label)){
        #multisigInternal$label <- rep(multisigInternal$label, nTracks)
        multisigInternal$label <- vector(mode = "list", length = nTracks)
    }
    
    # =========================================================================
    # CALL PLOTSIGNAL
    # =========================================================================
  
    ## Save to a group object?
    sigObjects <- purrr::pmap(list(multisigInternal$data, xList, yList, 
                    linecolorList, fillList, rep(FALSE, nTracks), 
                    multisigInternal$label), 
        \(d, x, y, l, f, draw, lab){
            suppressMessages(plotSignal(data = d, x = x, y = y, 
                        linecolor = l,
                        height = height, 
                        width = width,
                        range = range, 
                        orientation = multisigInternal$orientation, 
                        scale = multisigInternal$scale,
                        label = lab,
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
                        draw = draw))
        })
    
    # =========================================================================
    # GET AND DRAW GROBS
    # =========================================================================
    
    name <- paste0(
        "multisignal",
        length(grep(
            pattern = "multisignal",
            x = grid.ls(
                print = FALSE,
                recursive = FALSE
            )
        )) + 1
    )
    
    assign("pg_multiSignals", gTree(name = name), envir = pgEnv)
    invisible(lapply(sigObjects, getGrobs))
    
    sigGrobs <- get("pg_multiSignals", envir = pgEnv)

    # Call grid.draw to draw all these grobs at once
    if (multisigInternal$draw == TRUE){
        grid.draw(sigGrobs)
    }
    
    # =========================================================================
    # RETURN OBJECT
    # =========================================================================
    
    message("multisignal[", name, "]")
    invisible(sigObjects)
}  

#' Calculate a score range for multiple signals
#' 
#' @usage calcSignalRange(
#'     data,
#'     chrom = NULL,
#'     chromstart = 1,
#'     chromend = .Machine$integer.max,
#'     assembly = "hg38",
#'     negData = FALSE)
#' 
#' @param data List of data to be plotted as character values specifying 
#' multiple bigwig file paths, dataframes in BED format, or
#' \link[GenomicRanges]{GRanges} objects with metadata column \code{score}.
#' @param chrom Chromosome of data region ragne as a string, if range for a 
#' specific chromosome is desired.
#' @param chromstart Integer start position on chromosome to get data range.
#' @param chromend Integer end position on chromosome to get data range.
#' @param assembly Default genome assembly as a string or a
#' \link[plotgardener]{assembly} object.
#' Default value is \code{assembly = "hg38"}.
#' @param negData A logical value indicating whether any of the data has both
#' positive and negative scores and the signal range should be adjusted
#' accordingly. Default value is \code{negData = FALSE}. 
#' 
#' @return Returns a vector of length 2 with the calculated c(min, max) range.
#' 
#' @examples
#' library("plotgardenerData")
#' data("GM12878_ChIP_CTCF_signal")
#' data("IMR90_ChIP_CTCF_signal")
#' data("GM12878_ChIP_H3K27ac_signal")
#' data("IMR90_ChIP_H3K27ac_signal")
#' 
#' calcSignalRange(data = list(GM12878_ChIP_CTCF_signal, 
#'                             GM12878_ChIP_H3K27ac_signal, 
#'                             IMR90_ChIP_CTCF_signal, 
#'                             IMR90_ChIP_H3K27ac_signal),
#'                 chrom = "chr21",
#'                 chromstart = 28150000, chromend = 29150000,
#'                 assembly = "hg38", negData = FALSE)
#' 
#' @export
calcSignalRange <- function(data, chrom = NULL, chromstart = 1, 
                            chromend = .Machine$integer.max, 
                            assembly = "hg38", negData = FALSE){
    
    # Read and process data with read_rangeData()
    data <- lapply(data, read_rangeData,
                   assembly = assembly,
                   chrom = chrom,
                   start = chromstart,
                   end = chromend)
    
    # calculating the max from the processed data scores
    rangeMax <- lapply(data, dplyr::select, "score") %>%
        lapply(max) %>%
        unlist() %>%
        max
    # if no negative data, lower bound of range = 0
    if (negData == FALSE){
        rangeMin <- 0 
    }
        
    # if there are negative values, lower bound of range is the min value 
    else {
        rangeMin <- lapply(data, dplyr::select, "score") %>%
        lapply(min) %>%
        unlist() %>%
        min
    }
    return(c(rangeMin, rangeMax))
}