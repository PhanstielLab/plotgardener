#' plots HiC interaction matrix
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param resolution the width in bp of each pixel
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param height height of plot in inches
#' @param width width of plot in inches
#' @param x x-coordinate of where to place plot relative to top left of plot
#' @param y y-coordinate of where to place plot relative to top left of plot
#' @param palette palette to use for representing interaction scores
#' @param norm if giving .hic file, hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#' @param half what sides of square plot; options are "both", top", or "bottom"
#' @param raster allows for rasterization of plot, which results in quicker plotting
#' @param addlegend add legend representing scores for color range
#' @param legendlocation if addlegend == TRUE, where relative to plot to place legend; options are "right", "left", top", or "bottom"
#' @param legendoffset if addlegend == TRUE, how much to offset the legend relative to the plot in inches
#'
#' @return Function will plot a HiC interaction matrix and return a list of the zrange min, zrange max, and sequential color palette
#'
#' @author Nicole Kramer
#' @export

bb_hic <- function(hic, chrom, chromstart, chromend, resolution = 10000, zrange = NULL, height = 3, width = 3, x = 1.75, y = 3,
                   palette = 'reds', norm = "KR", half = "both", raster = TRUE,
                   addlegend = FALSE, legendlocation = "right", legendoffset = 0.25, ...){


  # ERRORS
  # ======================================================================================================================================================================================
  if (chromend < chromstart){
    stop("Incorrect \"chromstart\" and \"chromend\".  \"chromstart\" cannot be larger than \"chromend\".")
  }

  # DEFINE A FUNCTION TO PLOT SQUARES
  # ======================================================================================================================================================================================
  drawpoly <- function(df, resolution, chromstart, chromend){

    ## Define the color
    col = rgb(df[4], df[5], df[6], maxColorValue = 255)
    x = df[1]
    y = df[2]

    ## Get coordinates for the points of the square, normalized to 0 to 1 based on chromstart and chromend
    xleft = x - .5 * resolution
    xleft.normalized = normalize(xleft, chromstart, chromend)
    xright = x + .5 * resolution
    xright.normalized = normalize(xright, chromstart, chromend)
    ytop = y + .5 * resolution
    ytop.normalized = normalize(ytop, chromstart, chromend)
    ybottom = y - .5 * resolution
    ybottom.normalized = normalize(ybottom, chromstart, chromend)

    ## Plot all the squares
    grid.polygon(x = c(xleft.normalized, xleft.normalized, xright.normalized, xright.normalized),
                 y = c(ybottom.normalized, ytop.normalized, ytop.normalized, ybottom.normalized), gp = gpar(col = col, fill = col))
  }


  # CHECK FOR DATA TYPE, EXTRACT FROM .HIC FILE IF NECESSARY, SUBSET DATAFRAME/ADJUST DATAFRAME
  # ======================================================================================================================================================================================
  # Dataframe
  if (class(hic) == "data.frame"){
    ## Check for correct dataframe format
    if(ncol(hic) > 3 | ncol(hic) < 3){
      stop("Incorrect dataframe format.  Input a dataframe with 3 columns: x, y, counts.")
    } else {

        ## Make sure dataframe is subsetted and has correct header names for processing
        hicregion <- hic[which(hic[ ,1] >= chromstart & hic[ ,1] <= chromend &
                               hic[ ,2] >= chromstart & hic[ ,2] <= chromend), ]
        colnames(hicregion) <- c("x", "y", "counts")

        ## Scale lower and upper bounds using zrange
        if(is.null(zrange)){
          ## If counts vector only has one value, keep everything as is and zrange min = zrange max
          if(length(unique(hicregion$counts)) == 1){
            zrange <- c(unique(hicregion$counts), unique(hicregion$counts))
          } else {
            zrange <- c(min(hicregion$counts), max(hicregion$counts))
            hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
            hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
            }
        } else {
          stopifnot(is.vector(zrange), length(zrange) == 2, zrange[2] > zrange[1])
          hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
          hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
        }

        ## Detect what format of data we have and adjust if necessary based on plottype
        data_matrix <- as.matrix(reshape::cast(hicregion, formula = x ~ y, value = "counts"))

        ## If data is sparse but we want a square plot, we need to make the data symmetric
        if(all(is.na(data_matrix[lower.tri(data_matrix)])) & half == "both"){
          lower <- hicregion[ ,c(2,1,3)]
          colnames(lower) <- c("x", "y", "counts")
          combined <- unique(rbind(hicregion, lower))
          combinedComplete <- tidyr::complete(combined, x, y)
          combinedComplete$counts[is.na(combinedComplete$counts)] <- 0
          hicregion <- as.data.frame(combinedComplete)

          }

    }
  } else if (class(hic) == "character"){
    ## Make sure the inputted file is a .hic file
    if(file_ext(hic) == "hic"){

      if(half == "both"){

        hicregion <- bb_rhic(hic = hic, format = "full", chrom = chrom, chromstart = chromstart, chromend = chromend,
                                resolution = resolution, zrange = zrange, norm = norm)

        ## Change correct header names for processing
        colnames(hicregion) <- c("x", "y", "counts")

        ## bb_rhic will do this, but this will give the values for returning later
        if(is.null(zrange)){
          zrange <- c(min(hicregion$counts), max(hicregion$counts))
          hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
          hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
        }
      } else if (half == "bottom" | half == "top"){
        hicregion <- bb_rhic(hic = hic, format = "sparse", chrom = chrom, chromstart = chromstart, chromend = chromend,
                                resolution = resolution, zrange = zrange, norm = norm)

        ## Change correct header names for processing
        colnames(hicregion) <- c("x", "y", "counts")

        ## extractHiC will do this, but this will give the values for returning later
        if(is.null(zrange)){
          zrange <- c(0, max(hicregion$counts))
          hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
          hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
        }
      } else {
        stop("Incorrect \"half\" argument.  Options are \"both\", \"top\", or \"bottom\".")
      }
    } else {
      stop("Incorrect input. Input a .hic file or a dataframe.")
    }
  } else {
    stop("Incorrect input. Input a .hic file or a dataframe.")
  }

  # VIEWPORTS
  # ======================================================================================================================================================================================
  ## Go up a viewport if not at the root of all viewports


  if(is.null(current.vpPath()) == FALSE){
    upViewport()
  }

  ## Get page_height and margins from bbEnv
  page_height <- get("height", envir = bbEnv)
  top_margin <- get("top_margin", envir = bbEnv)
  left_margin <- get("left_margin", envir = bbEnv)

  ## Start x and y based on given margins
  x <- x + left_margin
  y <- y + top_margin

  ## Make viewport of desired size at specified location
  converted_coords = convert_coordinates(height = height, width = width, x = x, y = y, pageheight = page_height)
  vp <- viewport(height = unit(height, "in"), width = unit(width, "in"), x = unit(converted_coords[1], "in"), y = unit(converted_coords[2], "in"))
  pushViewport(vp)

  # CONVERT NUMBERS TO COLORS
  # ======================================================================================================================================================================================
  ## Use colour_values function from "colourvalues" package to convert numbers to colors
  color_vector <- colourvalues::colour_values(hicregion[ ,3], palette = palette)

  # Sorted color vector for use in legend and get lowest color to fill in for triangular plot
  sorted_color_vector <- colour_values(sort(hicregion[ ,3]), palette = palette)
  sorted_colors <- unique(sorted_color_vector)


  # RASTERIZED PLOT
  # ======================================================================================================================================================================================
  if(raster == TRUE){

    ## Add color vector to hicregion dataframe
    hicregion <- cbind(hicregion, color_vector)

    ## Get the lowest color in the range to fill in matrix for triangle plots
    lowest_color <- sorted_colors[1]

    ## Remove unnecessary "counts" column
    hicregion = hicregion[ ,c(1, 2, 4)]

    ## Cast dataframe into a matrix
    reshapen <- as.matrix(reshape::cast(hicregion, formula = x ~ y, value = "color_vector"))


    if(half == "bottom" | half == "top"){
      ## Fill in all NA's with the lowest color in the range (essentially a 0)
      reshapen[is.na(reshapen)] <- lowest_color

      ## Replace the lower triangular part of the matrix with NA's to not have anything plot there
      reshapen[lower.tri(reshapen)] <- NA
      if (half == "bottom"){
        reshapen <- apply(reshapen, 2, rev)
      } else if (half == "top"){
        reshapen <- apply(reshapen, 1, rev)
      }
    } else if (half == "both"){
      ## Matrix already complete from extraction, reverse orientation based on columns
      reshapen <- apply(reshapen, 2, rev)
    } else {
      stop("Invalid plot type.")
    }

    grid.raster(reshapen)

    ## Go back to root viewport
    upViewport()
  }

  # NON-RASTERIZED PLOT #
  # ======================================================================================================================================================================================
  if (raster == "FALSE") {

    ## Append colors to hicdata and convert to rgb
    hicregion <- cbind(hicregion, t(col2rgb(color_vector)))

    ## Need to get lower part of dataframe to plot bottom triangle
    if(half == "bottom"){
      lower <- hicregion[ , c(2, 1, 3, 4, 5, 6)]
      colnames(lower) <- c("x", "y", "counts", "red", "green", "blue")
      hicregion <- lower
    }

    ## Plot squares with drawpoly function defined above
    invisible(apply(hicregion, 1, drawpoly, resolution = resolution, chromstart = chromstart, chromend = chromend))

    ## Go back to root viewport
    upViewport()

  }

  # LEGEND
  # ======================================================================================================================================================================================
  if (addlegend == TRUE){

    if (is.null(legendlocation)){
      stop("Cannot add legend. Please specify a legend location.")
    }

    ## Get max and min labels based on zrange
    min_z <- zrange[1]
    max_z <- zrange[2]

    if(legendlocation == "right"){
      legend_width <- 1/32 * width
      legend_height <- 1/4 * height
      legend_xcoord <- x + width + legendoffset - left_margin
      legend_ycoord <- y - top_margin
      bb_legend(color_vector = sorted_colors, min_label = min_z, max_label = max_z, height= legend_height,
                    width = legend_width, x = legend_xcoord, y = legend_ycoord, orientation = "vertical")

    }

    else if(legendlocation == "left"){
      legend_width <- 1/32 * width
      legend_height <- 1/4 * height
      legend_xcoord <- x - legendoffset - left_margin - legend_width
      legend_ycoord <- y - top_margin
      bb_legend(color_vector = sorted_colors, min_label = min_z, max_label = max_z, height= legend_height,
                    width = legend_width, x = legend_xcoord, y = legend_ycoord, orientation = "vertical")
    }

    else if(legendlocation == "top"){
      legend_width <- 1/4 * height
      legend_height <- 1/32 * width
      legend_xcoord <- x + (width - legend_width)/2 - left_margin
      legend_ycoord <- y - legendoffset - legend_height - top_margin

      bb_legend(color_vector = sorted_colors, min_label = min_z, max_label = max_z, height= legend_height,
                    width = legend_width, x = legend_xcoord, y = legend_ycoord, orientation = "horizontal")
    }

    else if(legendlocation == "bottom"){
      legend_width <- 1/4 * height
      legend_height <- 1/32 * width
      legend_xcoord <- x + (width - legend_width)/2 - left_margin
      legend_ycoord <- y + legendoffset + height - top_margin

      bb_legend(color_vector = sorted_colors, min_label = min_z, max_label = max_z, height= legend_height,
                    width = legend_width, x = legend_xcoord, y = legend_ycoord, orientation = "horizontal")
    }

  }


  #return sorted_color vector and zrange if choosing to add legend later
  return(list(c(zrange[1], zrange[2]), sorted_colors))
}



