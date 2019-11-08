#' plots HiC interaction matrix
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param resolution the width in bp of each pixel
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param width width of plot in inches
#' @param height height of plot in inches
#' @param x x-coordinate of where to place plot relative to top left of plot
#' @param y y-coordinate of where to place plot relative to top left of plot
#' @param units units of height, width, and x and y location of the plot
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param norm if giving .hic file, hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#' @param half what sides of square plot; options are "both", top", or "bottom"
#' @param raster allows for rasterization of plot, which results in quicker plotting
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

bb_hic <- function(hic, chrom = "chr8", chromstart = 133600000, chromend = 134800000,
                   altchrom = NULL, altchromstart = NULL, altchromend = NULL, althalf = "bottom",
                   resolution = 10000, zrange = NULL, width = 3, height = 3, x = 1, y = 1, units = "inches",
                   palette = colorRampPalette(c("white", "dark red")), norm = "KR", half = "both", raster = TRUE, ...){


  # ERRORS
  # ======================================================================================================================================================================================
  if (chromend < chromstart){
    stop("Incorrect \"chromstart\" and \"chromend\".  \"chromstart\" cannot be larger than \"chromend\".")
  }

  # DEFINE A FUNCTION TO PLOT SQUARES
  # ======================================================================================================================================================================================
  drawpoly <- function(df, resolution, chrom = NULL, half, altchrom = NULL, althalf = NULL){

    ## Define the color
    col = rgb(df[4], df[5], df[6], maxColorValue = 255)
    x = df[1]
    y = df[2]

    xleft = x - .5 * resolution
    xright = x + .5 * resolution
    ytop = y + .5 * resolution
    ybottom = y - .5 * resolution

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

    col = rgb(df[4], df[5], df[6], maxColorValue = 255)
    x = df[1]
    y = df[2]

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


  # CHECK FOR DATA TYPE, EXTRACT FROM .HIC FILE IF NECESSARY, SUBSET DATAFRAME/ADJUST DATAFRAME
  # ======================================================================================================================================================================================
  # Dataframe
  if (class(hic) %in% "data.frame"){
    ## Check for correct dataframe format
    if(ncol(hic) > 3 | ncol(hic) < 3){
      stop("Incorrect dataframe format.  Input a dataframe with 3 columns: chrA, chrB, counts.")
    } else {

        ## If we have an altchrom, need to subset columns separately based on that
        if(!is.null(altchrom)){

          if(is.null(altchromstart) | is.null(altchromend)){
            stop("If specifying alternate chromosome, need to give altchromstart and altchromend.")
          }

          ## Check if altchrom is same as chrom (plotting off the diagonal)
          if(chrom == altchrom){

            ## Get minchromstart to minchromend and maxchromstart to maxchromend to subset columns appropriately (the same way as Straw)
            minchromstart <- min(chromstart, altchromstart)
            minchromend <- min(chromend, altchromend)
            maxchromstart <- max(chromstart, altchromstart)
            maxchromend <- max(chromend, altchromend)

            hicregion <- hic[which(hic[,1] >= minchromstart & hic[,1] <= minchromend &
                               hic[,2] >= maxchromstart & hic[,2] <= maxchromend),]

          } else {
            ## Make sure smaller chrom is in "x" column and larger chrom is in "y" column
            max_chrom <- paste0("chr", max(as.numeric(gsub(pattern = "chr", replacement = "", x = chrom)),
                                           as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))))
            min_chrom <- paste0("chr", min(as.numeric(gsub(pattern = "chr", replacement = "", x = chrom)),
                                           as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))))

            if (min_chrom != colnames(hic)[1] | max_chrom != colnames(hic)[2]){
              hic <- hic[ ,c(2, 1, 3)]
            }
            ## Subset region

            ## If chrom is smaller than altchrom, it will be in column 1, and vice versa
            if (as.numeric(gsub(pattern = "chr", replacement = "", x = chrom)) < as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))){
              hicregion <- hic[which(hic[,1] >= chromstart & hic[,1] <= chromend &
                                       hic[,2] >= altchromstart & hic[,2] <= altchromend),]
            } else if (as.numeric(gsub(pattern = "chr", replacement = "", x = chrom)) > as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))){
              hicregion <- hic[which(hic[,1] >= altchromstart & hic[,1] <= altchromend &
                                       hic[,2] >= chromstart & hic[,2] <= chromend),]
            }
          }

          ## Rename columns for later processing
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

          ## Make sure data is symmetric for plotting off the diagonal
          if (chrom == altchrom){
            lower <- hicregion[ ,c(2,1,3)]
            colnames(lower) <- c("x", "y", "counts")
            combined <- unique(rbind(hicregion, lower))
            combinedComplete <- tidyr::complete(combined, x, y)
            combinedComplete$counts[is.na(combinedComplete$counts)] <- 0
            hicregion <- as.data.frame(combinedComplete)
          }

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

          ## Detect what format of data we have and adjust if necessary based on half
          data_matrix <- as.matrix(reshape::cast(hicregion, formula = x ~ y, value = "counts"))

          ## If data is sparse but we want a half = "both" plot, we need to make the data symmetric
          if(all(is.na(data_matrix[lower.tri(data_matrix)])) & half == "both"){
            lower <- hicregion[ ,c(2,1,3)]
            colnames(lower) <- c("x", "y", "counts")
            combined <- unique(rbind(hicregion, lower))
            combinedComplete <- tidyr::complete(combined, x, y)
            combinedComplete$counts[is.na(combinedComplete$counts)] <- 0
            hicregion <- as.data.frame(combinedComplete)

          }

        }
        }


  } else if (class(hic) == "character"){

    ## Make sure the inputted file is a .hic file
    if(file_ext(hic) == "hic"){

      ## Straw is going to give output in sequential order, no matter how you order chrom and altchrom
      if(!is.null(altchrom)){

        hicregion <- bb_rhic(hic = hic, chrom = chrom, chromstart = chromstart, chromend = chromend,
                             resolution = resolution, zrange = zrange, norm = norm, altchrom = altchrom, altchromstart = altchromstart, altchromend = altchromend)

        ## Change header names for processing
        colnames(hicregion) <- c("x", "y", "counts")

        ## bb_rhic will do this, but this will give the values for returning later
        if(is.null(zrange)){
          zrange <- c(min(hicregion$counts), max(hicregion$counts))
          hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
          hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
        }

      } else {
        if(half == "both"){

          upper <- bb_rhic(hic = hic, chrom = chrom, chromstart = chromstart, chromend = chromend,
                               resolution = resolution, zrange = zrange, norm = norm)

          ## Get symmetric region
          colnames(upper) <- c("x", "y", "counts")
          lower <- upper[ , c(2, 1, 3)]
          colnames(lower) <- c("x", "y", "counts")
          combined <- unique(rbind(upper, lower))
          hicregion <- combined

          ## bb_rhic will do this, but this will give the values for returning later
          if(is.null(zrange)){
            zrange <- c(min(hicregion$counts), max(hicregion$counts))
            hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
            hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
          }
        } else if (half == "bottom" | half == "top"){

          hicregion <- bb_rhic(hic = hic, chrom = chrom, chromstart = chromstart, chromend = chromend,
                               resolution = resolution, zrange = zrange, norm = norm)

          ## Change correct header names for processing
          colnames(hicregion) <- c("x", "y", "counts")

          ## bb_rhic will do this, but this will give the values for returning later
          if(is.null(zrange)){
            zrange <- c(0, max(hicregion$counts))
            hicregion$counts[hicregion$counts <= zrange[1]] <- zrange[1]
            hicregion$counts[hicregion$counts >= zrange[2]] <- zrange[2]
          }
        } else {
          stop("Invalid \"half\" argument. Options are \"both\", \"top\", or \"bottom\".")
        }
      }
    } else {
      stop("Incorrect input. Input a .hic file or a dataframe.")
    }
  } else {
    stop("Incorrect input. Input a .hic file or a dataframe.")
  }

  # VIEWPORT NAVIGATION
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

  # RASTERIZED PLOT
  # ======================================================================================================================================================================================
  if(raster == TRUE){

    ## Add color vector to hicregion dataframe
    hicregion1 <- cbind(hicregion, color_vector)

    hicregion2 <- cbind(hicregion, t(col2rgb(color_vector)))

    ## Remove unnecessary "counts" column
    hicregion1 = hicregion1[ ,c(1, 2, 4)]

    ## Make clipped viewport
    if(!is.null(altchrom)){

      ## Off diagonal plotting
      if (chrom == altchrom){

        minchromstart <- min(chromstart, altchromstart)
        minchromend <- min(chromend, altchromend)
        maxchromstart <- max(chromstart, altchromstart)
        maxchromend <- max(chromend, altchromend)

        if (althalf == "bottom"){

          ## X AXIS = LARGER NUMBERS, Y AXIS = SMALLER NUMBERS

          vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                         y = unit(converted_coords[2], units = page_units), clip = "on", xscale = c(maxchromstart, maxchromend), yscale = c(minchromstart, minchromend))
          pushViewport(vp)

          reshapen <- as.matrix(reshape::cast(hicregion1, formula = x ~ y, value = "color_vector"))


        } else if (althalf == "top"){

          ## X AXIS = SMALLER NUMBERS, Y AXIS = LARGER NUMBERS

          vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                         y = unit(converted_coords[2], units = page_units), clip = "on", xscale = c(minchromstart, minchromend), yscale = c(maxchromstart, maxchromend))
          pushViewport(vp)

          reshapen <- as.matrix(reshape::cast(hicregion1, formula = y ~ x, value = "color_vector"))


        } else{
          stop("Invalid \"althalf\" argument. Options are \"top\" or \"bottom\".")
        }

      } else {

        ## Interactions between different chromosomes

        chrom <- as.numeric(gsub(pattern = "chr", replacement = "", x = chrom))
        altchrom <- as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))


        if (althalf == "top"){

          ## X AXIS = SMALLER CHROM, Y AXIS = LARGER CHROM

          if (chrom > altchrom){
            vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                           y = unit(converted_coords[2], units = page_units), clip = "on", xscale = c(altchromstart, altchromend), yscale = c(chromstart, chromend))
          } else if (altchrom > chrom){
            vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                           y = unit(converted_coords[2], units = page_units), clip = "on", xscale = c(chromstart, chromend), yscale = c(altchromstart, altchromend))
          }

          pushViewport(vp)
          reshapen <- as.matrix(reshape::cast(hicregion1, formula = y ~ x, value = "color_vector"))

        } else if (althalf == "bottom"){

          ## X AXIS = LARGER CHROM, Y AXIS = SMALLER CHROM

          if (chrom > altchrom){
            vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                           y = unit(converted_coords[2], units = page_units), clip = "on", xscale = c(chromstart, chromend), yscale = c(altchromstart, altchromend))
          } else if (altchrom > chrom){
            vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                           y = unit(converted_coords[2], units = page_units), clip = "on", xscale = c(altchromstart, altchromend), yscale = c(chromstart, chromend))
          }

          pushViewport(vp)
          reshapen <- as.matrix(reshape::cast(hicregion1, formula = x ~ y, value = "color_vector"))

        } else {
          stop ("Invalid \"althalf\" argument. Options are \"top\" or \"bottom\".")
        }
      }


      ## Get data in proper orientation
      reshapen <- apply(reshapen, 2, rev)

      ## Plot matrix of colors
      grid.raster(reshapen, interpolate = FALSE)


    } else {

      vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                     y = unit(converted_coords[2], units = page_units), clip = "on", xscale = c(chromstart, chromend), yscale = c(chromstart, chromend))
      pushViewport(vp)



      reshapen <- as.matrix(reshape::cast(hicregion1, formula = x ~ y, value = "color_vector"))


      if (half == "bottom"){

        ## Skip diagonal
        diag(reshapen) <- NA
        reshapen <- apply(reshapen, 2, rev)

      } else if (half == "top"){

        ## Skip diagonal
        diag(reshapen) <- NA
        reshapen <- apply(reshapen, 1, rev)

      } else if (half == "both"){
        reshapen <- apply(reshapen, 2, rev)

      } else {
        stop("Invalid plot type.")
      }

      assign("TEST", reshapen, envir = globalenv())
      ## Plot rasterized version of everything (if "half" plot, not plotting diagonal)
      grid.raster(reshapen, interpolate = FALSE)

      ## Go back and plot triangles along diagonal for "half" plots
      if (half == "bottom" | half == "top"){

        invisible(apply(hicregion2, 1, drawpoly_diagonal, resolution = resolution, half = half))
      }

    }

    ## Go back to root viewport
    upViewport()
  }

  # NON-RASTERIZED PLOT #
  # ======================================================================================================================================================================================
  if (raster == "FALSE") {

    ## Append colors to hicdata and convert to rgb

    hicregion <- cbind(hicregion, t(col2rgb(color_vector)))

    ## Plotting altchrom
    if(!is.null(altchrom)){

      if (chrom == altchrom){

        minchromstart <- min(chromstart, altchromstart)
        minchromend <- min(chromend, altchromend)
        maxchromstart <- max(chromstart, altchromstart)
        maxchromend <- max(chromend, altchromend)


        if (althalf == "bottom"){

          ## X AXIS = LARGER NUMBERS, Y AXIS = SMALLER NUMBERS

          vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                         y = unit(converted_coords[2], units = page_units), xscale = c(maxchromstart, maxchromend), yscale = c(minchromstart, minchromend))
          pushViewport(vp)

          reverse <- hicregion[,c(2, 1, 3, 4, 5, 6)]
          colnames(reverse) <- c("x", "y", "counts", "red", "green", "blue")
          hicregion <- reverse

        } else if (althalf == "top"){

          ## X AXIS = SMALLER NUMBERS, Y AXIS = LARGER NUMBERS

          vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                         y = unit(converted_coords[2], units = page_units), xscale = c(minchromstart, minchromend), yscale = c(maxchromstart, maxchromend))
          pushViewport(vp)


        }else{
          stop("Invalid \"althalf\" argument. Options are \"top\" or \"bottom\".")
        }

      } else {

        chrom <- as.numeric(gsub(pattern = "chr", replacement = "", x = chrom))
        altchrom <- as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))

        if (althalf == "bottom"){

          if (chrom > altchrom){
            vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                           y = unit(converted_coords[2], units = page_units), xscale = c(chromstart, chromend), yscale = c(altchromstart, altchromend))
          } else if (altchrom > chrom){
            vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                           y = unit(converted_coords[2], units = page_units), xscale = c(altchromstart, altchromend), yscale = c(chromstart, chromend))
          }

          pushViewport(vp)

          reverse <- hicregion[, c(2, 1, 3, 4, 5, 6)]
          colnames(reverse) <- c("x", "y", "counts", "red", "green", "blue")
          hicregion <- reverse

        } else if (althalf == "top"){
          if (chrom > altchrom){
            vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                           y = unit(converted_coords[2], units = page_units), xscale = c(altchromstart, altchromend), yscale = c(chromstart, chromend))
          } else if (altchrom > chrom){
            vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                           y = unit(converted_coords[2], units = page_units), xscale = c(chromstart, chromend), yscale = c(altchromstart, altchromend))
          }

          pushViewport(vp)

        }

        chrom <- paste0("chr", chrom)
        altchrom <- paste0("chr", altchrom)
      }

      invisible(apply(hicregion, 1, drawpoly, resolution = resolution, chrom = chrom, half = half, altchrom = altchrom, althalf = althalf))

    } else {

      vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                     y = unit(converted_coords[2], units = page_units), xscale = c(chromstart, chromend), yscale = c(chromstart, chromend))
      pushViewport(vp)

      ## Need to get lower part of dataframe to plot bottom triangle
      if(half == "bottom"){
        lower <- hicregion[ , c(2, 1, 3, 4, 5, 6)]
        colnames(lower) <- c("x", "y", "counts", "red", "green", "blue")
        hicregion <- lower
      }

      ## Plot squares with drawpoly function defined above
      invisible(apply(hicregion, 1, drawpoly, resolution = resolution, half = half))
    }

    ## Go back up a viewport
    upViewport()

  }

  # if (!is.null(altchrom)){
  #   hic_plot <- new("hic_plot", chrom = chrom, chromstart = chromstart, chromend = chromend, color_palette = sorted_colors, zrange = c(zrange[1], zrange[2]),
  #                   x = x, y = y, height = height, width = width, units = units, altchrom = altchrom, altchromstart = altchromstart, altchromend = altchromend)
  # } else {
  #   hic_plot <- new("hic_plot", chrom = chrom, chromstart = chromstart, chromend = chromend, color_palette = sorted_colors, zrange = c(zrange[1], zrange[2]),
  #                   x = x, y = y, height = height, width = width, units = units)
  # }


  #return(hic_plot)
}



