#' Maps a numeric or character vector to a color palette and returns
#' the vector of colors
#' @param vector Vector to map to color.
#' @param palette Color palette function.
#' @param range Range of values to map for a numerical value.
#' 
#' @examples 
#' ## Load paired ranges data in BEDPE format
#' library(plotgardenerData)
#' data("IMR90_DNAloops_pairs")
#' 
#' ## Add a length column
#' IMR90_DNAloops_pairs$length <- (IMR90_DNAloops_pairs$start2 - 
#'         IMR90_DNAloops_pairs$start1) / 1000
#' 
#' ## Map length column to a vector of colors
#' colors <- mapColors(vector = IMR90_DNAloops_pairs$length,
#'             palette = colorRampPalette(c("dodgerblue2", "firebrick2")))
#'             
#' ## Pass color vector into bbPlotPairsArches
#' heights <- IMR90_DNAloops_pairs$length / max(IMR90_DNAloops_pairs$length)    
#' pageCreate(width = 7.5, height = 2.1, default.units = "inches",
#'             showGuides = FALSE, xgrid = 0, ygrid = 0)       
#' params <- pgParams(
#'     chrom = "chr21",
#'     chromstart = 27900000, chromend = 30700000,
#'     assembly = "hg19",
#'     width = 7
#' )
#' 
#' archPlot <- plotPairsArches(
#'     data = IMR90_DNAloops_pairs, params = params,
#'     fill = colors,
#'     linecolor = "fill",
#'     archHeight = heights, alpha = 1,
#'     x = 0.25, y = 0.25, height = 1.5,
#'     just = c("left", "top"),
#'     default.units = "inches"
#' )            
#' 
#' annoGenomeLabel(plot = archPlot, x = 0.25, y = 1.78, scale = "Mb") 
#' annoHeatmapLegend(
#'     plot = archPlot, fontcolor = "black",
#'     x = 7.0, y = 0.25,
#'     width = 0.10, height = 1, fontsize = 10
#' )
#' plotText(
#'     label = "Kb", rot = 90, x = 6.9, y = 0.75,
#'     just = c("center", "center"),
#'     fontsize = 10
#' )
#' 
#' @details 
#' This function allows for the manual mapping of a numerical or factor
#' vector to a palette of colors. For a more automatic implementation 
#' of this functionality in plotgardener functions,  
#' \link[plotgardener]{colorby} objects can be used.
#' 
#' @return 
#' Returns a character vector of color values. If the input vector is 
#' numerical, this vector will have additional `palette` and `range`
#' attributes.
#' 
#' @seealso \link[plotgardener]{colorby}
#' 
#' @export
mapColors <- function(vector, palette, range = NULL){
    
    ## Define a function to catch errors
    error_mapColors <- function(vector, palette, range){
        
        ## palette errors
        if (!is(palette, "function")){
            stop("Please provide a palette function.", call. = FALSE)
        }
        
        ## range errors
        rangeErrors(range = range)
        
        ## Numerical vector for breaks
        if (is(vector, "numeric") | is(vector, "integer")){
            if (length(unique(vector)) == 1){
                warning("Not enough numerical values to map ",
                "to colors.", call. = FALSE)
                vector <- NULL
            }
        }
        
        return(vector)
        
    }
    
    ## Catch errors and update to NULL vector if necessary
    vector <- error_mapColors(vector = vector,
                            palette = palette,
                            range = range)
    
    if (is(vector, "numeric") | is(vector, "integer")){
        
        ## Update range, if necessary
        if (is.null(range)){
            breaks <- seq(min(vector), max(vector), length.out = 100)
            range <- c(min(vector), max(vector))
        } else {
            vector[which(vector < range[1])] <- range[1]
            vector[which(vector > range[2])] <- range[2]
            breaks <- seq(range[1], range[2], length.out = 100)
        }
        
        ## Map numbers to colors    
        colors <- palette(length(breaks) + 1)    
        colorVector <- as.character(cut(vector, c(-Inf, breaks, Inf), 
                                        labels = colors))
        attr(colorVector, "range") <- range
        
    } else {
        
        ## Convert if not a factor
        if (!is(vector, "factor")){
            vector <- as.factor(vector)
        }
        
        ## Map color palette to factor levels
        colors <- palette(length(levels(vector)))
        names(colors) <- levels(vector)
        colorVector <- colors[vector]
        
    } 
    attr(colorVector, "palette") <- palette
    return(colorVector)
    
}

## Color palette and range default assignments for mapColors
## Returns an updated object with the color palette and range
colorDefaults <- function(vector, palette = NULL, range = NULL, object){
    if (is(vector, "numeric") | is(vector, "integer")){
        
        if (is.null(palette)){
            palette <- colorRampPalette(brewer.pal(
                n = 9, "YlGnBu"
            ))
        }
        
        if (is.null(range)){
            
            range <- c(min(vector), max(vector))
            
            ## Switch back to NULL if invalid range
            if (range[1] >= range[2]){
                range <- NULL
            } 
            
        }
        
        object$zrange <- range
        
    } else {
        
        ## Convert if not a factor
        if (!is(vector, "factor")){
            vector <- as.factor(vector)
        }
        
        if (is.null(palette)){
            palette <- colorRampPalette(suppressWarnings(
                brewer.pal(n = length(levels(vector)),
                        "Paired")))
        }
        
    }
    
    object$color_palette <- palette
    
    return(object)
}

## Define a function that will parse a vector of colors vs. a colorby object
## Returns the final vector of colors and an updated plot object
# @param data Associated data, for finding `colorby` column
# @param fill Input fill - will either be a single value, a vector, or 
# a colorby object
# @param object The plot object, to be updated with any color_palette
# and zrange information
# @param subset A string describing the type of data, which will determine
# how to subset it. Options are ranges, pairs, pairs_clip, or manhattan.
parseColors <- function(data, fill, object, subset = NULL){
    
    ## `colorby` class
    if (is(fill, "colorby")){

        colorbyCol <- data[, fill$column]
        
        ## Scale numeric colorby data by the subsetted plotted region
        if ((is(colorbyCol, "numeric") | is(colorbyCol, "integer")) 
                & is.null(object$zrange)){
            
            if (fill$scalePerRegion == TRUE){
                if (subset == "ranges"){
                    subData <- data[which(data[,"chrom"] == object$chrom &
                                    data[,"start"] <= object$chromend &
                                    data[,"end"] >= object$chromstart),]
                    
                } else if (subset == "pairs"){
                    subData <- data[which(data[,"chrom1"] == object$chrom &
                                        data[,"chrom2"] == object$chrom &
                                        ((data[,"end1"] >= object$chromstart &
                                        data[,"end1"] <= object$chromend) |
                                        (data[,"start2"] <= object$chromstart &
                                        data[,"start2"] >= object$chromend))),]
                    subData <- data[which(data[, "chrom1"] == object$chrom &
                                            data[, "chrom2"] == object$chrom),]
                    overlappingRanges <- 
                        as.data.frame(subsetByOverlaps(ranges = 
                                IRanges(start = object$chromstart, 
                                        end = object$chromend),
                                x = IRanges(start = subData[,"start1"], 
                                            end = subData[,"end2"])))
                    subData <- subData[which(subData[,"start1"] %in% 
                                                overlappingRanges$start &
                                            subData[,"end2"] %in% 
                                                overlappingRanges$end),]
                } else if (subset == "pairs_clip"){
                    subData <- data[which(data[,"chrom1"] == object$chrom &
                                data[,"chrom2"] == object$chrom &
                                data[,"start1"] >= object$chromstart &
                                data[,"end1"] <= object$chromend &
                                data[,"start2"] >= object$chromstart &
                                data[,"end2"] <= object$chromend),]
                    
                } else if (subset == "pairs_noachor"){
                    subData <- data[which(data[,"chrom1"] == object$chrom &
                                        data[,"chrom2"] == object$chrom &
                                        ((data[,"start1"] >= object$chromstart &
                                        data[,"start1"] <= object$chromend) |
                                        (data[,"end2"] >= object$chromstart &
                                        data[,"end2"] <= object$chromend))), ]
                } else if (subset == "manhattan"){
                    subData <- data[which(data[,"chrom"] == object$chrom &
                                            data[,"pos"] >= object$chromstart &
                                            data[,"pos"] <= object$chromend),]
                } else {
                    subData <- data
                }
                
                range <- range(subData[,fill$column])
                if (range[1] >= range[2]){
                    range <- NULL
                }
                
                fill$range <- range
            }
            
        }
        
        ## Default palette and range
        object <- colorDefaults(vector = colorbyCol,
                                palette = fill$palette,
                                range = fill$range,
                                object = object)
        
        ## Pass into mapColors
        colors <- mapColors(vector = colorbyCol,
                            palette = object$color_palette,
                            range = object$zrange)
    } else {
        
        if (length(fill) == 1){
            colors <- as.character(rep(fill, nrow(data)))
        } else {
            colors <- as.character(rep(fill,
                        ceiling(nrow(data) / length(fill))
            )[seq(1, nrow(data))])
        }
        
        if (!is.null(attributes(fill))){
            object$color_palette <- attr(fill, "palette")
            object$zrange <- attr(fill, "range")
        }
        
    }
    
    return(list(colors, object))
    
}

## Define a function that will parse the objects for linecolors and 
## return a vector of linecolors
# @param linecolor The linecolor parameter value
# @param fillcolors The mapped vector of fillcolors, for use if 
# linecolor == "fill"
# @param data Associated data, for finding `colorby` column
# @param object The plot object
# @param subset A string describing the type of data, which will determine
# how to subset it. Options are ranges, pairs, or pairs_clip.
lineColors <- function(linecolor, fillcolors, data, object, subset = NULL){
    
    if (is.null(linecolor)){
        linecolor <- NA
    }
    
    if (!all(is.na(linecolor))){
        if (all(linecolor == "fill")){
            linecolors <- fillcolors
        } else {
            linecolors <- parseColors(data = data,
                                        fill = linecolor,
                                        object = object,
                                        subset = subset)[[1]]
        }
    } else {
        linecolors <- NA
    }
    
    return(linecolors)
    
}

## Define a function that makes a color transparent
# @param color color string
# @param alpha Alpha value of color
makeTransparent <- function(color, alpha) {
    if (is.null(alpha)) {
        alpha <- 1
    }

    rgb <- grDevices::col2rgb(color)
    transp <- rgb(rgb[1], rgb[2], rgb[3],
        alpha = alpha * 255,
        maxColorValue = 255
    )
    return(transp)
}
