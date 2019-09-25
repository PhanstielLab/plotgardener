#' plots HiC interaction matrix
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param half what sides of square plot; options are "both", top", or "bottom"
#' @param resolution the width in bp of each pixel
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param raster allows for rasterization of plot, which results in quicker plotting
#' @param width width of plot in inches
#' @param height height of plot in inches
#' @param x x-coordinate of where to place plot relative to top left of plot
#' @param y y-coordinate of where to place plot relative to top left of plot
#' @param units units of height, width, and x and y location of the plot
#' @param altchrom alternate chromosome for off-diagonal plotting or interchromosomal plotting
#' @param altchromstart alternate chromosome start for off-diagonal plotting or interchromosomal plotting
#' @param altchromend alternate chromosome end for off-diagonal plotting or interchromosomal plotting
#' @param althalf if plotting altchrom region, which side off diagonal to plot; options are "both", "top", or "bottom"
#' @param norm if giving .hic file, hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#'
#' @details macOS Preview anti-aliases images and will make rasterized plot appear blurry.
#'
#' @return Function will plot a HiC interaction matrix and return a list of the zrange min, zrange max, and sequential color palette
#'
#' @author Nicole Kramer
#'
#' @examples
#' data(bb_hic_data)
#' # Full Rasterized Plot
#' bb_makepage()
#' bb_hic_plot <- bb_hic(bb_hic_data)
#' # Half Rasterized Plot
#' bb_makepage()
#' bb_hic_plot <- bb_hic(bb_hic_data, half = "top")
#' # Non-rasterized Plot
#' bb_makepage()
#' bb_hic_plot <- bb_hic(bb_hic_data, raster = FALSE)
#'
#' @export
#'
#'
bb_plothic <- function(hic, chrom = "chr8", chromstart = 133600000, chromend = 134800000, half = "both", resolution = 10000, zrange = NULL,
                       palette = colorRampPalette(c("white", "dark red")), raster = TRUE, width = 3, height = 3, x = 1, y = 1,
                       units = "inches", altchrom = NULL, altchromstart = NULL, altchromend = NULL, althalf = "bottom", norm = "KR", ...){

  # ERRORS
  # ======================================================================================================================================================================================
  if (chromend < chromstart){

    stop("Invalid \"chromstart\" and \"chromend\".  \"chromstart\" cannot be larger than \"chromend\".")

  }

  if (!is.null(altchrom)){

    if(is.null(altchromstart) | is.null(altchromend)){

      stop("If specifying alternate chromosome, need to give altchromstart and altchromend.")

    }

    if (altchromend < altchromstart){

      stop("Invalid \"altchromstart\" and \"altchromend\".  \"altchromstart\" cannot be larger than \"altchromend\".")

    }

  }


  # DEFINE FUNCTIOND TO PLOT SQUARES/TRIANGLES
  # ======================================================================================================================================================================================
  drawpoly <- function(df, resolution, chrom = NULL, half, altchrom = NULL, althalf = NULL){

    col = df[3]
    x = as.numeric(df[1])
    y = as.numeric(df[2])

    xleft = x
    xright = x + resolution
    ybottom = y
    ytop = y + resolution

    if (!is.null(altchrom)){

      grid.polygon(x = c(xleft, xleft, xright, xright),
                   y = c(ybottom, ytop, ytop, ybottom), gp = gpar(col = NA, fill = col), default.units = "native")

    } else {

      if (half == "both"){

        ## Plot all squares
        grid.polygon(x = c(xleft, xleft, xright, xright),
                     y = c(ybottom, ytop, ytop, ybottom), gp = gpar(col = NA, fill = col), default.units = "native")
      } else if (half == "top"){

        ## Plot triangles along diagonal and squares above
        if (y > x){
          grid.polygon(x = c(xleft, xleft, xright, xright),
                       y = c(ybottom, ytop, ytop, ybottom), gp = gpar(col = NA, fill = col), default.units = "native")
        } else if (y == x) {
          grid.polygon(x = c(xleft, xleft, xright),
                       y = c(ybottom, ytop, ytop), gp = gpar(col = NA, fill = col), default.units = "native")
        }

      } else if (half == "bottom"){
        ## Plot triangles along diagonal and squares below
        if (y < x){
          grid.polygon(x = c(xleft, xleft, xright, xright),
                       y = c(ybottom, ytop, ytop, ybottom), gp = gpar(col = NA, fill = col), default.units = "native")
        } else if (y == x) {
          grid.polygon(x = c(xleft, xright, xright),
                       y = c(ybottom, ybottom, ytop), gp = gpar(col = NA, fill = col), default.units = "native")
        }
      }

    }

  }



  drawpoly_diagonal <- function(df, resolution, half){

    col = df[3]
    x = as.numeric(df[1])
    y = as.numeric(df[2])

    xleft = x - .5 * resolution
    xright = x + resolution + .5 * resolution
    ytop = y + resolution + .5 * resolution
    ybottom = y - .5 * resolution

    if (half == "top"){

      if (y == x){

        grid.polygon(x = c(xleft, xleft, xright),
                     y = c(ybottom, ytop, ytop), gp = gpar(col = NA, fill = col), default.units = "native")
      }

    } else if (half == "bottom"){

      if (y == x){

        grid.polygon(x = c(xleft, xright, xright),
                     y = c(ybottom, ybottom, ytop), gp = gpar(col= NA, fill = col), default.units = "native")
      }
    }
  }

  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  # Numeric chromosome
  chromNumber <- as.numeric(gsub(pattern = "chr", replacement = "", x = chrom))

  if (!is.null(altchrom)){

    # Numeric altchromosome
    altchromNumber <- as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))

    # Name of larger chromosome
    max_chrom <- paste0("chr", max(as.numeric(gsub(pattern = "chr", replacement = "", x = chrom)),
                                   as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))))

    # Name of smaller chromosome
    min_chrom <- paste0("chr", min(as.numeric(gsub(pattern = "chr", replacement = "", x = chrom)),
                                   as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))))

    if (chromNumber == altchromNumber){

      ## If same chromsome, "x" column of Straw hic counts dataframe is from minchromstart to minchromend and "y" columns is
      ## maxchromstart to maxchromend

      minchromstart <- min(chromstart, altchromstart)
      minchromend <- min(chromend, altchromend)
      maxchromstart <- max(chromstart, altchromstart)
      maxchromend <- max(chromend, altchromend)

    }

  }


  # CHECK FOR DATA TYPE, EXTRACT FROM .HIC FILE IF NECESSARY, SUBSET DATAFRAME/ADJUST DATAFRAME
  # ======================================================================================================================================================================================

  if (class(hic) %in% "data.frame"){
    ## DATAFRAME

    ## Check for correct dataframe format
    if(ncol(hic) > 3 | ncol(hic) < 3){

      stop("Invalid dataframe format.  Input a dataframe with 3 columns: chrA, chrB, counts.")

    } else {

      if(!is.null(altchrom)){
        ## ALTERNATE CHROMOSOME

        if(chromNumber == altchromNumber){
          ## SAME CHROM AND ALTCHROM (OFF DIAGONAL)

          ## Subset data
          hicregion <- hic[which(hic[,1] >= minchromstart & hic[,1] <= minchromend &
                                   hic[,2] >= maxchromstart & hic[,2] <= maxchromend),]

        } else {
          ## DIFFERENT CHROM AND ALTCHROM (INTERCHROMOSOMAL)

          ## Make sure smaller chrom is in "x" column and larger chrom is in "y" column
          if (min_chrom != colnames(hic)[1] | max_chrom != colnames(hic)[2]){

            ## Swap columns if not in proper order
            hic <- hic[ ,c(2, 1, 3)]
          }

          ## Subset data
          if (chromNumber < altchromNumber){
            ## If chrom < altchrom, column 1 = chrom and column 2 = altchrom

            hicregion <- hic[which(hic[,1] >= chromstart & hic[,1] <= chromend &
                                     hic[,2] >= altchromstart & hic[,2] <= altchromend),]

          } else if (chromNumber > altchromNumber){
            ## If chrom > altchrom, column 1 = altchrom and column 2 = chrom

            hicregion <- hic[which(hic[,1] >= altchromstart & hic[,1] <= altchromend &
                                     hic[,2] >= chromstart & hic[,2] <= chromend),]
          }

        }

        ## Rename columns for later processing
        colnames(hicregion) <- c("x", "y", "counts")

      } else {
        ## NO ALTERNATE CHROMOSOME

        ## Subset data
        hicregion <- hic[which(hic[ ,1] >= chromstart & hic[ ,1] <= chromend &
                                 hic[ ,2] >= chromstart & hic[ ,2] <= chromend), ]

        ## Rename columns for later processing
        colnames(hicregion) <- c("x", "y", "counts")

        ## If half == "both", need to make data symmetric
        if (half == "both"){

          lower <- hicregion[, c(2, 1, 3)]
          colnames(lower) <- c("x", "y", "counts")
          combined <- unique(rbind(hicregion, lower))
          hicregion <- combined

        }
      }

      ## Scale lower and upper bounds using zrange
      if(is.null(zrange)){
        ## If counts vector only has one value, keep everything as is and zrange min = zrange max
        if(length(unique(hicregion$counts)) == 1){
          zrange <- c(unique(hicregion$counts), unique(hicregion$counts))
        } else {
          zrange <- c(min(hicregion$counts), max(hicregion$counts))
        }
      } else {
        stopifnot(is.vector(zrange), length(zrange) == 2, zrange[2] > zrange[1])
        hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
        hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
      }

    }

  } else if (class(hic) == "character"){
    ## HIC FILE

    ## Make sure the inputted file is a .hic file
    if(file_ext(hic) == "hic"){

      if(!is.null(altchrom)){
        ## ALTERNATE CHROMOSOME

        hicregion <- bb_rhic(hic = hic, chrom = chrom, chromstart = chromstart, chromend = chromend,
                             resolution = resolution, zrange = zrange, norm = norm, altchrom = altchrom,
                             altchromstart = altchromstart, altchromend = altchromend)

        ## Rename columns for later processing
        colnames(hicregion) <- c("x", "y", "counts")

      } else {
        ## NO ALTERNATE CHROMOSOME

        hicregion <- bb_rhic(hic = hic, chrom = chrom, chromstart = chromstart, chromend = chromend,
                             resolution = resolution, zrange = zrange, norm = norm)

        ## Rename columns for later processing
        colnames(hicregion) <- c("x", "y", "counts")

        ## If half == "both", need to make data symmetric
        if(half == "both"){

          lower <- hicregion[ , c(2, 1, 3)]
          colnames(lower) <- c("x", "y", "counts")
          combined <- unique(rbind(hicregion, lower))
          hicregion <- combined

        }
      }

      ## Get zrange of data to return later
      if(is.null(zrange)){

        zrange <- c(min(hicregion$counts), max(hicregion$counts))

      }

    } else {

      stop("Invalid input. Input a .hic file or a dataframe.")

    }
  } else {

    stop("Invalid input. Input a .hic file or a dataframe.")

  }

  # DIMENSIONS AND COORDINATES
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)


  ## Convert x and y coordinates and height and width to same page_units
  old_x <- unit(x, units = units)
  old_y <- unit(y, units = units)
  old_height <- unit(height, units = units)
  old_width <- unit(width, units = units)
  new_x <- convertX(old_x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(old_y, unitTo = page_units, valueOnly = TRUE)
  new_height <- convertHeight(old_height, unitTo = page_units, valueOnly = TRUE)
  new_width <- convertWidth(old_width, unitTo = page_units, valueOnly = TRUE)

  ## Convert coordinates for viewport
  converted_coords = convert_coordinates(height = new_height, width = new_width, x = new_x, y = new_y, pageheight = page_height)

  # CONVERT NUMBERS TO COLORS
  # ======================================================================================================================================================================================
  ## Use bb_maptocolors to convert numbers to colors

  color_vector <- bb_maptocolors(hicregion[ ,3], col = palette, num = 100, range = zrange)

  # Sorted color vector for use in legend
  sorted_color_vector <- bb_maptocolors(sort(hicregion[ ,3]), col = palette, num = 100, range = zrange)
  sorted_colors <- unique(sorted_color_vector)

  ## Add color vector to hicregion dataframe
  hicregion <- cbind(hicregion, color_vector)

  ## Remove unnecessary "counts" column
  hicregion = hicregion[ ,c(1, 2, 4)]

  ## VIEWPORTS
  # ======================================================================================================================================================================================

  if (!is.null(altchrom)){

    if(chromNumber == altchromNumber){

      if (althalf == "bottom"){

        xscale = c(maxchromstart, maxchromend + resolution)
        yscale = c(minchromstart, minchromend + resolution)

      } else if (althalf == "top"){

        xscale = c(minchromstart, minchromend + resolution)
        yscale = c(maxchromstart, maxchromend + resolution)

      }


    } else if (chromNumber > altchromNumber){


      if (althalf == "bottom"){

        xscale = c(chromstart, chromend + resolution)
        yscale = c(altchromstart, altchromend + resolution)

      } else if (althalf == "top"){

        xscale = c(altchromstart, altchromend + resolution)
        yscale = c(chromstart, chromend + resolution)

      }

    } else if (chromNumber < altchromNumber){

      if (althalf == "bottom"){

        xscale = c(altchromstart, altchromend + resolution)
        yscale = c(chromstart, chromend + resolution)

      } else if (althalf == "top"){

        xscale = c(chromstart, chromend + resolution)
        yscale = c(altchromstart, altchromend + resolution)

      }

    }

  } else {

    xscale = c(chromstart, chromend + resolution)
    yscale = c(chromstart, chromend + resolution)

  }

  vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                 y = unit(converted_coords[2], units = page_units), clip = "on",
                 xscale = xscale, yscale = yscale)

  pushViewport(vp)

  ## PLOTTING
  # ======================================================================================================================================================================================

  if (raster == TRUE){

    ## Cast dataframe in raster matrices

    if (!is.null(altchrom)){

      if (althalf == "top"){

        reshapen <- as.matrix(reshape::cast(hicregion, formula = y ~ x, value = "color_vector"))

      } else if (althalf == "bottom"){

        reshapen <- as.matrix(reshape::cast(hicregion, formula = x ~ y, value = "color_vector"))
      }

    } else {

      if (half == "both"){

        reshapen <- as.matrix(reshape::cast(hicregion, formula = x ~ y, value = "color_vector"))

      } else if (half == "top"){

        reshapen <- as.matrix(reshape::cast(hicregion, formula = y ~ x, value = "color_vector"))

        ## Skip diagonal
        diag(reshapen) <- NA

      } else if (half == "bottom"){

        reshapen <- as.matrix(reshape::cast(hicregion, formula = x ~ y, value = "color_vector"))

        ## Skip diagonal
        diag(reshapen) <- NA

      }

    }

    ## Get matrix in proper orientation
    reshapen <- apply(reshapen, 2, rev)

    ## Plot matrix of colors
    grid.raster(reshapen, interpolate = FALSE)

    if (half == "bottom" | half == "top"){

      ## Plot diagonal with triangles
      invisible(apply(hicregion, 1, drawpoly_diagonal, resolution = resolution, half = half))

    }

  } else {
    ## NO RASTER

    if (!is.null(altchrom)){

      if (althalf == "bottom"){

        reverse <- hicregion[,c(2, 1, 3)]
        colnames(reverse) <- c("x", "y", "color_vector")
        hicregion <- reverse

      }

      ## Plot squares with drawpoly function defined above
      invisible(apply(hicregion, 1, drawpoly, resolution = resolution, chrom = chrom, half = half, altchrom = altchrom, althalf = althalf))

    } else {

      if (half == "bottom"){

        reverse <- hicregion[ , c(2, 1, 3)]
        colnames(reverse) <- c("x", "y", "color_vector")
        hicregion <- reverse

      }

      ## Plot squares with drawpoly function defined above
      invisible(apply(hicregion, 1, drawpoly, resolution = resolution, half = half))

    }

  }

  ## Go back up a viewport
  upViewport()

  return(list("chrom" = chrom, "chromstart" = chromstart, "chromend" = chromend, "x" = x, "y" = y,
              "height" = height, "width" = width, "units" = units, "color_palette" = sorted_colors, "zrange" = c(zrange[1], zrange[2])))

  # if (!is.null(altchrom)){
  #   hic_plot <- new("hic_plot", chrom = chrom, chromstart = chromstart, chromend = chromend, color_palette = sorted_colors, zrange = c(zrange[1], zrange[2]),
  #                   x = x, y = y, height = height, width = width, units = units, altchrom = altchrom, altchromstart = altchromstart, altchromend = altchromend)
  # } else {
  #   hic_plot <- new("hic_plot", chrom = chrom, chromstart = chromstart, chromend = chromend, color_palette = sorted_colors, zrange = c(zrange[1], zrange[2]),
  #                   x = x, y = y, height = height, width = width, units = units)
  # }


  #return(hic_plot)
  }
