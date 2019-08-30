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
  drawpoly <- function(df, resolution, chrom = NULL, chromstart, chromend, half, altchrom = NULL, altchromstart = NULL, altchromend = NULL, althalf = NULL){

    ## Define the color
    col = rgb(df[4], df[5], df[6], maxColorValue = 255)
    x = df[1]
    y = df[2]

    ## If altchrom, normalize x and y scales differently
    if (!is.null(altchrom)){

      ## If chrom == altchrom, we subsetted the first column with the lower range and the second column with the higher range
      if(chrom == altchrom){

        if(althalf == "top"){

          x.min <- min(chromstart, altchromstart)
          x.max <- min(chromend, altchromend)
          y.min <- max(chromstart, altchromstart)
          y.max <- max(chromend, altchromend)


        } else if (althalf == "bottom"){

          x.min <- max(chromstart, altchromstart)
          x.max <- max(chromend, altchromend)
          y.min <- min(chromstart, altchromstart)
          y.max <- min(chromend, altchromend)

        }

      } else {

        ## Separate altchrom/chrom to just get number and determine which is larger
        chrom <- as.numeric(gsub(pattern = "chr", replacement = "", x = chrom))
        altchrom <- as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))

        ## "TOP": smaller chrom on x and larger chrom on y
        if (althalf == "top"){

          if (chrom > altchrom){
            x.min <- altchromstart
            x.max <- altchromend
            y.min <- chromstart
            y.max <- chromend
          } else {
            x.min <- chromstart
            x.max <- chromend
            y.min <- altchromstart
            y.max <- altchromend
          }
          ## "BOTTOM": larger chrom on x and smaller chrom on y
        } else if (althalf == "bottom"){

          if (chrom > altchrom){
            x.min <- chromstart
            x.max <- chromend
            y.min <- altchromstart
            y.max <- altchromend

          } else {

            x.min <- altchromstart
            x.max <- altchromend
            y.min <- chromstart
            y.max <- chromend

          }

        }

      }

      xleft = x - .5 * resolution
      xleft.normalized = normalize(xleft, x.min, x.max)
      xright = x + .5 * resolution
      xright.normalized = normalize(xright, x.min, x.max)
      ytop = y + .5 * resolution
      ytop.normalized = normalize(ytop, y.min, y.max)
      ybottom = y - .5 * resolution
      ybottom.normalized = normalize(ybottom, y.min, y.max)
      grid.polygon(x = c(xleft.normalized, xleft.normalized, xright.normalized, xright.normalized),
                   y = c(ybottom.normalized, ytop.normalized, ytop.normalized, ybottom.normalized), gp = gpar(col = NA, fill = col))

    } else {

      ## Get coordinates for the points of the square/triangle, normalized to 0 to 1 based on chromstart and chromend
      xleft = x - .5 * resolution
      xleft.normalized = normalize(xleft, chromstart, chromend)
      xright = x + .5 * resolution
      xright.normalized = normalize(xright, chromstart, chromend)
      ytop = y + .5 * resolution
      ytop.normalized = normalize(ytop, chromstart, chromend)
      ybottom = y - .5 * resolution
      ybottom.normalized = normalize(ybottom, chromstart, chromend)

      if (half == "both"){

        ## Plot all squares
        grid.polygon(x = c(xleft.normalized, xleft.normalized, xright.normalized, xright.normalized),
                     y = c(ybottom.normalized, ytop.normalized, ytop.normalized, ybottom.normalized), gp = gpar(col = NA, fill = col))
      } else if (half == "top"){

        ## Plot triangles along diagonal and squares above
        if (y > x){
          grid.polygon(x = c(xleft.normalized, xleft.normalized, xright.normalized, xright.normalized),
                       y = c(ybottom.normalized, ytop.normalized, ytop.normalized, ybottom.normalized), gp = gpar(col = NA, fill = col))
        } else if (y == x) {
          grid.polygon(x = c(xleft.normalized, xleft.normalized, xright.normalized),
                       y = c(ybottom.normalized, ytop.normalized, ytop.normalized), gp = gpar(col = NA, fill = col))
        }

      } else if (half == "bottom"){
        ## Plot triangles along diagonal and squares below
        if (y < x){
          grid.polygon(x = c(xleft.normalized, xleft.normalized, xright.normalized, xright.normalized),
                       y = c(ybottom.normalized, ytop.normalized, ytop.normalized, ybottom.normalized), gp = gpar(col = NA, fill = col))
        } else if (y == x) {
          grid.polygon(x = c(xleft.normalized, xright.normalized, xright.normalized),
                       y = c(ybottom.normalized, ybottom.normalized, ytop.normalized), gp = gpar(col = NA, fill = col))
        }
      }

    }

  }

  drawpoly_diagonal <- function(df, resolution, chromstart, chromend, half){

    col = rgb(df[4], df[5], df[6], maxColorValue = 255)
    x = df[1]
    y = df[2]

    xleft = x - .5 * resolution
    xleft.normalized = normalize(xleft, chromstart, chromend)
    xright = x + resolution + .5 * resolution
    xright.normalized = normalize(xright, chromstart, chromend)
    ytop = y + resolution + .5 * resolution
    ytop.normalized = normalize(ytop, chromstart, chromend)
    ybottom = y - .5 * resolution
    ybottom.normalized = normalize(ybottom, chromstart, chromend)

    if (half == "top"){
      if (y == x){
        grid.polygon(x = c(xleft.normalized, xleft.normalized, xright.normalized),
                     y = c(ybottom.normalized, ytop.normalized, ytop.normalized), gp = gpar(col = NA, fill = col))
      }

    } else if (half == "bottom"){
      if (y == x){
        grid.polygon(x = c(xleft.normalized, xright.normalized, xright.normalized),
                     y = c(ybottom.normalized, ybottom.normalized, ytop.normalized), gp = gpar(col= NA,fill = col))
      }
    }
  }


  # CHECK FOR DATA TYPE, EXTRACT FROM .HIC FILE IF NECESSARY, SUBSET DATAFRAME/ADJUST DATAFRAME
  # ======================================================================================================================================================================================
  # Dataframe
  if (class(hic) == "data.frame"){
    ## Check for correct dataframe format
    if(ncol(hic) > 3 | ncol(hic) < 3){
      stop("Incorrect dataframe format.  Input a dataframe with 3 columns: x, y, counts.")
    } else {

        ## If we have an altchrom, need to subset columns separately based on that
        if(!is.null(altchrom)){

          if(is.null(altchromstart) | is.null(altchromend)){
            stop("If specifying alternate chromosome, need to give altchromstart and altchromend.")
          }

          ## Check if altchrom is same as chrom (plotting off the diagonal)
          if(chrom == altchrom){

            ## Compare chromstarts of chrom and altchrom to scale first column with lower chromstart and second column with higher chromstart
            if (chromstart > altchromstart){
              hicregion <- hic[which(hic[,1] >= altchromstart & hic[,1] <= altchromend &
                                       hic[,2] >= chromstart & hic[,2] <= chromend), ]
            } else if (altchromstart > chromstart){
              hicregion <- hic[which(hic[,1] >= chromstart & hic[,1] <= chromend &
                                       hic[,2] >= altchromstart & hic[,2] <= altchromend), ]
            } else {
              ## If they don't have the same start, but have different ending, subset first column based on chrom and second column based on altchrom
              hicregion <- hic[which(hic[,1] >= chromstart & hic[,1] <= chromend &
                                       hic[,2] >= altchromstart & hic[,2] <= altchromend), ]
            }

          } else {
            ## Subset region
            hicregion <- hic[which(hic$chrom >= chromstart & hic$chrom <= chromend &
                                     hic$altchrom >= altchromstart & hic$altchrom <= altchromend), ]

            ##
          }

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

          ## Make sure data is symmetric
          lower <- hicregion[ ,c(2,1,3)]
          colnames(lower) <- c("x", "y", "counts")
          combined <- unique(rbind(hicregion, lower))
          combinedComplete <- tidyr::complete(combined, x, y)
          combinedComplete$counts[is.na(combinedComplete$counts)] <- 0
          hicregion <- as.data.frame(combinedComplete)


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

        if (chrom == altchrom){

          ## Get symmetric region in case capturing section that includes diagonal for same chromosome
          hicregion <- bb_rhic(hic = hic, format = "full", chrom = chrom, chromstart = chromstart, chromend = chromend,
                               resolution = resolution, zrange = zrange, norm = norm, altchrom = altchrom, altchromstart = altchromstart, altchromend = altchromend)

        } else {

          hicregion <- bb_rhic(hic = hic, format = "sparse", chrom = chrom, chromstart = chromstart, chromend = chromend,
                               resolution = resolution, zrange = zrange, norm = norm, altchrom = altchrom, altchromstart = altchromstart, altchromend = altchromend)

        }

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

  # Sorted color vector for use in legend and for getting lowest color to fill in for triangular plot
  sorted_color_vector <- bb_maptocolors(sort(hicregion[ ,3]), col = palette, num = 100, range = zrange)
  sorted_colors <- unique(sorted_color_vector)

  lowest_color <- sorted_colors[1]

  # RASTERIZED PLOT
  # ======================================================================================================================================================================================
  if(raster == TRUE){

    ## Make clipped viewport
    vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units),
                   y = unit(converted_coords[2], units = page_units), clip = "on")
    pushViewport(vp)
    grid.rect()


    ## Add color vector to hicregion dataframe
    hicregion1 <- cbind(hicregion, color_vector)
    hicregion2 <- cbind(hicregion, t(col2rgb(color_vector)))

    ## Remove unnecessary "counts" column
    hicregion1 = hicregion1[ ,c(1, 2, 4)]


    ## Plotting altchromosome
    if(!is.null(altchrom)){

      ## Off digonal plotting
      if (chrom == altchrom){

        ## Need to subset for bottom or top
        if (althalf == "bottom"){

          ## Cast dataframe into a matrix with x vs y
          reshapen <- as.matrix(reshape::cast(hicregion1, formula = y ~ x, value = "color_vector"))


          ## Subset data for bottom side
          reshapen <- reshapen[which(rownames(reshapen) >= min(chromstart, altchromstart)  & rownames(reshapen) <= min(chromend, altchromend)),
                               which(colnames(reshapen) >= max(chromstart, altchromstart) & colnames(reshapen) <= max(chromend, altchromend))]


        } else if (althalf == "top"){

          ## Cast dataframe into a matrix with y vs x
          reshapen <- as.matrix(reshape::cast(hicregion1, formula = y ~ x, value = "color_vector"))

          ## Subset data for top side
          reshapen <- reshapen[which(rownames(reshapen) >= max(chromstart, altchromstart) & rownames(reshapen) <= max(chromend, altchromend)),
                               which(colnames(reshapen) >= min(chromstart, altchromstart) & colnames(reshapen) <= min(chromend, altchromend))]

        } else{
          stop("Invalid \"althalf\" argument. Options are \"top\" or \"bottom\".")
        }


      } else {

        ## Interactions between different chromosomes
        if (althalf == "top"){

          ## Cast dataframe into a matrix with x vs y (lower chrom on x and higher chrom on y)
          reshapen <- as.matrix(reshape::cast(hicregion1, formula = y ~ x, value = "color_vector"))

        } else if (althalf == "bottom"){

          ## Cast dataframe into a matrix with y vs x (higher chrom on x and lower chrom on x)
          reshapen <- as.matrix(reshape::cast(hicregion1, formula = x ~ y, value = "color_vector"))

        } else {
          stop ("Invalid \"althalf\" argument. Options are \"top\" or \"bottom\".")
        }
      }

      ## Fill in NA's with lowest color
      reshapen[is.na(reshapen)] <- lowest_color

      ## Get data in proper orientation
      reshapen <- apply(reshapen, 2, rev)

      ## Plot matrix of colors
      grid.raster(reshapen, interpolate = FALSE)

    } else {

      reshapen <- as.matrix(reshape::cast(hicregion1, formula = x ~ y, value = "color_vector"))
      if(half == "bottom" | half == "top"){

        ## Fill in all NA's with lowest color
        reshapen[is.na(reshapen)] <- lowest_color

        ## Replace the lower triangular part of the matrix with NA's to not have anything plot there
        reshapen[lower.tri(reshapen)] <- NA

        if (half == "bottom"){
          ## Skip diagonal
          diag(reshapen) <- NA
          reshapen <- apply(reshapen, 2, rev)
        } else if (half == "top"){
          ## Skip diagonal
          diag(reshapen) <- NA
          reshapen <- apply(reshapen, 1, rev)
        }
      } else if (half == "both"){
        ## Fill in all NA's with lowest color
        reshapen[is.na(reshapen)] <- lowest_color
        ## Matrix already complete from extraction, reverse orientation based on columns
        reshapen <- apply(reshapen, 2, rev)
      } else {
        stop("Invalid plot type.")
      }

      ## Plot rasterized version of everything (if "half" plot, not plotting diagonal)
      grid.raster(reshapen, interpolate = FALSE)

      ## Go back and plot triangles along diagonal for "half" plots
      if (half == "bottom" | half == "top"){
        invisible(apply(hicregion2, 1, drawpoly_diagonal, resolution = resolution, chromstart = chromstart, chromend = chromend, half = half))
      }

    }

    ## Go back to root viewport
    upViewport()
  }

  # NON-RASTERIZED PLOT #
  # ======================================================================================================================================================================================
  if (raster == "FALSE") {

    ## Create unclipped viewport
    vp <- viewport(height = unit(new_height, page_units), width = unit(new_width, page_units), x = unit(converted_coords[1], units = page_units), y = unit(converted_coords[2], units = page_units))
    pushViewport(vp)

    ## Append colors to hicdata and convert to rgb
    hicregion <- cbind(hicregion, t(col2rgb(color_vector)))

    ## Plotting altchrom
    if(!is.null(altchrom)){

      if (chrom == altchrom){

        ## Have to subset data again since we grabbed stuff beyond the diagonal
        if (althalf == "bottom"){
          hicregion <- hicregion[which(hicregion[,1] >= max(chromstart, altchromstart) & hicregion[,1] <= max(chromend, altchromend)
                                       & hicregion[,2] >= min(chromstart, altchromstart) & hicregion[,2] <= min(chromend, altchromend)),]

        } else if (althalf == "top"){
          hicregion <- hicregion[which(hicregion[,1] >= min(chromstart, altchromstart) & hicregion[,1] <= min(chromend, altchromend)
                                       & hicregion[,2] >= max(chromstart, altchromstart) & hicregion[,2] <= max(chromend, altchromend)),]

        }else{
          stop("Invalid \"althalf\" argument. Options are \"top\" or \"bottom\".")
        }

      } else {

        if (althalf == "bottom"){

          reverse <- hicregion[, c(2, 1, 3, 4, 5, 6)]
          colnames(reverse) <- c("x", "y", "counts", "red", "green", "blue")
          hicregion <- reverse

        }

      }


      ## Plot squares with drawpoly function defined above
      invisible(apply(hicregion, 1, drawpoly, resolution = resolution, chrom = chrom, chromstart = chromstart,
                      chromend = chromend, half = half, altchrom = altchrom, altchromstart = altchromstart, altchromend = altchromend, althalf = althalf))

    } else {

      ## Need to get lower part of dataframe to plot bottom triangle
      if(half == "bottom"){
        lower <- hicregion[ , c(2, 1, 3, 4, 5, 6)]
        colnames(lower) <- c("x", "y", "counts", "red", "green", "blue")
        hicregion <- lower
      }

      ## Plot squares with drawpoly function defined above
      invisible(apply(hicregion, 1, drawpoly, resolution = resolution, chromstart = chromstart, chromend = chromend, half = half))
    }

    ## Go back up a viewport
    upViewport()

  }

  return(list("chrom" = chrom, "chromstart" = chromstart, "chromend" = chromend, "x" = x, "y" = y,
              "height" = height, "width" = width, "units" = units, "color_palette" = sorted_colors, "zrange" = c(zrange[1], zrange[2])))
}



