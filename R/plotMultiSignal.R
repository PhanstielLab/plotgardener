#' Plot multiple signal tracks inline with eachother
#' 
#' @usage plotMultiSignal(
#'     x,
#'     y,
#'     default.units = "inches",
#'     linecolor = "black",
#'     lwd = 1,
#'     lty = 1,
#'     fill = NA,
#'     alpha = 1,
#'     id = NULL,
#'     id.lengths = NULL,
#'     params = NULL,
#'     ...
#' )
#'
#' @param x A numeric vector or unit object specifying polygon
#' vertex x-locations.
#' @param y A numeric vector, unit object, or a character vector
#' of values containing a "b" combined with a numeric value specifying
#' polygon vertex y-locations.
#' The character vector will place polygon vertex y-locations relative
#' to the bottom of the most recently plotted plot according
#' to the units of the plotgardener page.
#' @param default.units A string indicating the default units to use
#' if \code{x} or \code{y} are only given as numeric vectors.
#' Default value is \code{default.units = "inches"}.
#' @param linecolor A character value specifying polygon line color.
#' Default value is \code{linecolor = "black"}.
#' @param lwd A numeric specifying polygon line width.
#'  Default value is \code{lwd = 1}.
#' @param lty A numeric specifying polygon line type.
#' Default value is \code{lty = 1}.
#' @param fill A character value specifying polygon fill color.
#' Default value is \code{fill = NA}.
#' @param alpha Numeric value specifying color transparency.
#' Default value is \code{alpha = 1}.
#' @param id A numeric vector used to separate locations in \code{x} and
#' \code{y} into multiple polygons. All locations with the same \code{id}
#' belong to the same polygon.
#' @param id.lengths A numeric vector used to separate locations in
#' \code{x} and \code{y} into multiple polygons. Specifies consecutive
#' blocks of locations which make up separate polygons.
#' @param params An optional \link[plotgardener]{pgParams} object containing
#' relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{polygon} object containing relevant
#' placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a page
#' pageCreate(width = 7.5, height = 6, default.units = "inches")
#'
#' ## Plot complex polygons one at a time
#' plotPolygon(
#'     x = c(2.6, 4.65, 4.75, 6.05, 1.4, 1.3),
#'     y = c(2.5, 3.1, 3.5, 4, 3.15, 2.8),
#'     fill = "#4a168e", linecolor = NA
#' )
#'
#' plotPolygon(
#'     x = c(4.65, 4.75, 6.05, 5.05, 4.4),
#'     y = c(3.1, 3.5, 4, 1.45, 1.2),
#'     fill = "#9d28b0", linecolor = NA
#' )
#'
#' ## Plot multiple triangles with different id's and colors
#' plotPolygon(
#'     x = c(
#'         0.45, 6.05, 3, 3, 6.05, 5.25, 4.4, 5.05, 4.95,
#'         1.3, 2.6, 1, 4.4, 4.95, 5, 4.95, 5, 6.25
#'     ),
#'     y = c(
#'         2.85, 4, 5.55, 5.55, 4, 5.55, 1.2, 1.45, 1.1,
#'         2.8, 2.5, 2.1, 1.2, 1.1, 0.45, 1.1, 0.45, 1.1
#'     ),
#'     id = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
#'     fill = c(
#'         "#ce93d9", "#bb6ac9", "#4a168e",
#'         "#7b1fa0", "#bb6ac9", "#ce93d9"
#'     ),
#'     linecolor = NA
#' )
#'
#########################################################################################
#########################################################################################
#########################################################################################
findSignalRange <- function(data){
   
  rangeMax<- lapply(data, select, "score") %>%
      lapply(max) %>%
      unlist() %>%
      max
   return(c(0,rangeMax))
}

setColors <- function(colorList, nTracks){
  if(length(colorList) == nTracks)
    return(colorList)
  else if(length(colorList) < nTracks)
    setColors(append(colorList, colorList), nTracks)
  else if(length(colorList) > nTracks)
    setColors(colorList[1:nTracks], nTracks)
}

## calculating new x,y coordinates of each track
getXCoordinates<- function(x, nTracks, width, orientation, gapdistance){
  xList<- rep(x, nTracks)
  if(orientation == "h"){
   return(xList)
  } else if(orientation == "v"){
    internalWidth<- (width - (gapdistance * (nTracks-1)))/nTracks
    for (i in 2:nTracks)
      xList[i]<- (xList[i-1] + internalWidth + gapdistance)
    return(xList)
  } else{
    stop("argument \" orientation\" is inco, ","with no default.", call. = FALSE)
  }
}
## calculating new x,y coordinates of each track
getYCoordinates<- function(y, nTracks, height, orientation, gapdistance){
  yList<- rep(y, nTracks)
  if(orientation == "h"){
    internalHeight<- (height - (gapdistance * (nTracks-1)))/nTracks
    for (i in 2:nTracks)
      yList[i]<- (yList[i-1] + internalHeight + gapdistance)
    return(yList)
  } else if(orientation == "v"){
      return(yList)
  } else{
    stop("argument \" orientation\" is missing, ","with no default.", call. = FALSE)
  }
}

#' @export
plotMultiSignal<- function(data, binSize = NA, binCap = TRUE, negData = FALSE,
                           chrom, chromstart = NULL, chromend = NULL,
                           assembly = "hg38", linecolor,
                           fill = NA, ymax = 1, range = NULL, scale = FALSE,
                           bg = NA, baseline = TRUE, baseline.color = "grey",
                           baseline.lwd = 1, orientation = "Teresa",
                           x = NULL, y = NULL, width = NULL,
                           height = NULL, just = c("left", "top"),
                           default.units = "inches", gapdistance = .2 ,draw = TRUE,
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
  
  
  
  
  #nTracks<- length(data)
  nTracks <- length(multisigInternal$data)
  #range<- findSignalRange(data)
  range <- findSignalRange(multisigInternal$data)
  
  #xList<- getXCoordinates(x = x, nTracks = length(data), width = width, orientation = orientation, gapdistance = gapdistance)
  xList <- getXCoordinates(x = multisigInternal$x, 
                           nTracks = length(multisigInternal$data), 
                           width = multisigInternal$width, 
                           orientation = multisigInternal$orientation, 
                           gapdistance = multisigInternal$gapdistance)
  #yList <- getYCoordinates(y = y,  nTracks = length(data), height = height, orientation = orientation, gapdistance = gapdistance)
 
  yList <- getYCoordinates(y = multisigInternal$y,  
                           nTracks = length(multisigInternal$data), 
                           height = multisigInternal$height, 
                           orientation = multisigInternal$orientation, 
                           gapdistance = multisigInternal$gapdistance)
  
  
  # if(orientation == "h"){
  #   height<- (height - (gapdistance * (nTracks-1)))/nTracks
  # } else if(orientation == "v"){
  #   width<- (width - (gapdistance * (nTracks-1)))/nTracks
  # }
  # else{
  #   #throw error
  # }
  
  if (multisigInternal$orientation == "h"){
    height <- (multisigInternal$height - (multisigInternal$gapdistance * (nTracks-1)))/nTracks
  } else if (multisigInternal$orientation == "v"){
    width <- (multisigInternal$width - (multisigInternal$gapdistance * (nTracks-1)))/nTracks
  }
  else{
    #throw error
  }
  
  colorList<- setColors(linecolor, nTracks)
  
  pmap(list(multisigInternal$data, xList, yList, colorList), \(d,x,y,c){
       # plotSignal(data = d, x = x, y = y, linecolor = c,
       #            params = params, height = height, width = width, 
       #            range = range, orientation = orientation, scale = scale, 
       #            binSize = binSize, binCap = binCap, negData = negData, 
       #            chrom = chrom, chromstart = chromstart, chromend = chromend, 
       #            assembly = assembly, fill = fill, ymax = ymax, bg = bg, baseline = baseline, 
       #            baseline.color = baseline.color, baseline.lwd = baseline.lwd, just = just, 
       #            default.units = default.units, draw = draw)
        plotSignal(data = d, x = x, y = y, linecolor = c,
                   height = height, width = width,
                   range = range, orientation = multisigInternal$orientation, 
                   scale = multisigInternal$scale,
                   binSize = multisigInternal$binSize, binCap = multisigInternal$binCap, 
                   negData = multisigInternal$negData,
                   chrom = multisigInternal$chrom, chromstart = multisigInternal$chromstart, 
                   chromend = multisigInternal$chromend,
                   assembly = multisigInternal$assembly, fill = multisigInternal$fill, 
                   ymax = multisigInternal$ymax, bg = multisigInternal$bg, baseline = multisigInternal$baseline,
                   baseline.color = multisigInternal$baseline.color, baseline.lwd = multisigInternal$baseline.lwd, 
                   just = multisigInternal$just,
                   default.units = multisigInternal$default.units, draw = multisigInternal$draw)
     })
}  

