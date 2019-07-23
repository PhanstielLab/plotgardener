#' @export

gridplotBedgraph <- function(signal, chrom, chromstart, chromend, range=NULL, color, transparency = 1.0, lwd =1, linecolor = NA, binSize = NA, addscale = FALSE, overlay =FALSE, binCap = TRUE, flip = FALSE, ymax = 1.04, width = 3.25, height = .625, x=4.25, y =5.5, ...  ){
  require(grid)
  if(is.na(binSize) == TRUE){
    binSize = (chromend - chromstart)/2000
  }
  
  if (is.na(linecolor) == TRUE){
    linecolor = color
  }
  
  #define function to normalize data to 0 - 1 scale
  normalize <- function(val, minval, maxval){
    return((val-minval)/(maxval-minval))
  }
  
  #ensure that the chromosome is a character
  signal[,1] = as.character(signal[,1])
  
  #filter for desired region
  signaltrack = signal[which(signal[,1] == chrom & ((signal[,2] > chromstart & signal[,2] < chromend | signal[,3] > chromstart & signal[,3] < chromend ))), (2:4)]
  
  #remove any duplicate rows
  signaltrack = signaltrack[!duplicated(signaltrack),]
  
  #exit if overlay is TRUE and there isn't enough data
  if(overlay == TRUE && nrow(signaltrack) < 2){
    return ("not enough data within range to plot")
  }
  
  #exit if overlay is FALSE and there isn't enough data
  if (nrow(signaltrack) < 2){
    
    return ("not enough data within range to plot")
  
  }
  
  binNum = (chromend-chromstart)/binSize
  
  #scale back binNum and print warning if binNum is greater than 8000
  if (binNum > 8000 && binCap == TRUE){
    binNum = 8000
    binSize = (chromend-chromstart)/binNum
    warning(paste0("Too many bins: adjusting to 8000 bins of size ", binSize, ". To override try binCap = FALSE."))
  }
  
  #scale bin size to 1 if binNum is larger than span
  if (binNum > (chromend - chromstart)){
    binNum = (chromend - chromstart)
    binSize = 1
    warning(paste0("Number of bins larger than plot length: adjusting to ", binNum, " bins of size 1."))
  }
  

  bin.dataframe <- data.frame(seq(chromstart, chromend-binSize, binSize), seq(chromstart+binSize, chromend, binSize), rep(0, times=binNum))
  
  
  # add column names
  colnames(bin.dataframe) = c("V1", "V2", "V3")
  
  #find the max signal values for each bin
  bin.signal <- function(line){ 
    signal=signaltrack
    line = as.integer(line)
    list=c(0)
    list = append(list, signal[,3][which((signal[,1] >= line[1] & signal[,1] < line[2]) |
                                           signal[,2] > line[1] & signal[,2] <= line[2] |
                                           signal[,1] < line[1] & signal[,2] > line[2])])
    return(max(list))
  }
  
  bin.dataframe[,3] = apply(bin.dataframe, 1, bin.signal)
  
  #use binned data as signal track
  signaltrack = bin.dataframe
  
  # make linking regions if neccesary
  linkingregions = cbind(signaltrack[1:(nrow(signaltrack)-1),2], signaltrack[2:nrow(signaltrack),1])
  
  linkingregions = matrix(linkingregions[which(linkingregions[,1] != linkingregions[,2]),],ncol=2)
  
  if (nrow(linkingregions) > 0)
  {
    linkingregions = cbind(linkingregions,0)
    
    # make col names the same
    names(linkingregions)[(1:3)] = c("V1","V2","V3")
    
    # add linking regions to signaltrack
    signaltrack = rbind(signaltrack,linkingregions)
  }
  
  # sort data
  signaltrack = signaltrack[order(signaltrack[,1]),]
  
  # convert two columns to one
  signaltrack = cbind(as.vector(t(signaltrack[,c(1,2)])),as.vector(t(signaltrack[,c(3,3)])))
  
  # add slighltly negative value to both ends to ensure proper polygon plotting
  #signaltrack = rbind(c(min(signaltrack[,1]),-.00001),signaltrack)
  #signaltrack = rbind(signaltrack, c(max(signaltrack[,1]),-.00001))
  
  if (flip == TRUE){
    signaltrack[,2] = signaltrack[,2]*-1
  }

  #determine the y-limits
  if (is.null(range) == TRUE){
    range = c(0,ymax*max(signaltrack[,2]))
    if (flip == TRUE){
      range = c(ymax*min(signaltrack[,2]),0)
    }
    
  }
  

  #original range before normalization
  y0 = range[1]
  y1 = range[2]

  
  
  ##normalize x from 0 to 1 based on chromstart/chromend
  signaltrack[,1] = (signaltrack[,1]-chromstart)/(chromend-chromstart)
  ##normalize y from 0 to 1 based on min/max
  maxsignal = max(signaltrack[,2])
  minsignal = min(signaltrack[,2])
  signaltrack[,2] = (signaltrack[,2]-minsignal)/(maxsignal-minsignal)

  
  #set the transparency
  rgbcol = col2rgb(color)
  finalcolor = rgb(rgbcol[1],rgbcol[2],rgbcol[3],alpha=transparency * 255,maxColorValue = 255)
  
  ###PLOT THE SIGNAL TRACK###
  if (overlay == FALSE){
    if (is.null(current.vpPath()) == FALSE){
      upViewport()
    }
    converted_coords = convert_coordinates(height, width, x, y)
    vp <- viewport(width = unit(width,"in"), height = unit(height, "in"), x = unit(converted_coords[1],"in"), y=unit(converted_coords[2],"in"))
    pushViewport(vp)
  } 
  if (overlay == TRUE){
    vp <- current.vpPath()
  }
  
    #grid.polygon(x=signaltrack[,1],y=signaltrack[,2],gp=gpar(fill=NA,lwd=lwd, col=linecolor))
    grid.segments(x0=signaltrack[c(1:1-length(signaltrack[,1])),1],y0=signaltrack[c(1:1-length(signaltrack[,2])),2],x1=signaltrack[c(2:length(signaltrack[,1])),1],y1 = signaltrack[c(2:length(signaltrack[,2])),2],gp=gpar(col=linecolor,lwd=lwd),name="plotBedgraph track")
    
  
  #add scale to upper right corner or bottom right corner
  if (addscale == TRUE){

    if (flip == FALSE){
      grid.text(paste(y0,y1,sep="-"),x=max(signaltrack[,1]),y=max(signaltrack[,2]), gp=gpar(col="grey"),name="bedgraph scale")
    }
    if (flip == TRUE){
      grid.text(paste(y1,y0*-1,sep="-"),x=max(signaltrack[,1]),y=min(signaltrack[,2]),gp=gpar(col="grey"),name=" bedgraph scale")
    }
  }
  
    return(list("viewport" = vp, "maxsignal" = maxsignal, "minsignal"= minsignal,"flip" = flip)) 
}
#gridplotBedgraph(Sushi_ChIPSeq_CTCF.bedgraph,chr="chr11",chromstart = 1955000,chromend = 1965000,color="blue", x=6.125, y =7.3687,transparency = 0.5)



