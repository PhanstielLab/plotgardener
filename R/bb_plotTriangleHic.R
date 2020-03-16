#' plots HiC interaction matrix in triangular format
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted, based on build (i.e. for hg19 just a number, for hg38 string like "chr8")
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param resolution the width in bp of each pixel; options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param width A numeric or unit object specifying the bottom width of the triangle
#' @param height A numeric or unit object specifying the height of the triangle
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param just a string or numeric vector specifying the justification of the viewport relative to its (x, y) location
#' @param draw A logical value indicating whether graphics output should be produced
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numerics
#' @param norm if giving .hic file, hic data normalization; must be found in hic file
#'
#' @export
#'
#'
bb_plotTriangleHic <- function(hic, chrom, chromstart, chromend, resolution = 10000, zrange = NULL,
                               palette = colorRampPalette(c("white", "dark red")), width = NULL, height = NULL,
                               x = NULL, y = NULL, just = c("left", "top"), norm = "KR", default.units = "inches", draw = T, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## For more accurate calculation of sqrt(2)
  two <- mpfr(2, 120)

  ## Define a function that resets the just based on if the final plot will be a triangle or a trapezoid
  reset_just <- function(just, x, y, width, height){

    if (!is.null(x) & !is.null(y)){

      two <- mpfr(2, 120)
      desired_height <- convertHeight(height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
      calc_height <- convertWidth(width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)*0.5
      side_length <- (convertWidth(width, unitTo = get("page_units", envir = bbEnv), valueOnly = T))/sqrt(two)

      if (calc_height <= desired_height){

        ## here we'll have a triangle
        if (length(just == 2)){

          if (identical(just, c("left", "top")) | identical(just, c("right", "top"))){

            just <- "top"
            message("Entire triangle will be plotted.  Auto-adjusting plot justifiction to top.")
          }

        }

      }

    }
    return(just)
  }

  ## Define a function that catches errors for bb_plotTriangleHic
  errorcheck_bb_plotTriangleHic <- function(hic, hic_plot, norm){

    ###### hic/norm #####

    ## if it's a dataframe or datatable, it needs to be properly formatted
    if ("data.frame" %in% class(hic) && ncol(hic) != 3){

      stop("Invalid dataframe format.  Input a dataframe with 3 columns: chrA, chrB, counts.", call. = FALSE)

    }

    if (!"data.frame" %in% class(hic)){

      ## if it's a file path, it needs to be a .hic file
      if (file_ext(hic) != "hic"){

        stop("Invalid input. File must have a \".hic\" extension", call. = FALSE)

      }

      ## if it's a file path, it needs to exist
      if (!file.exists(hic)){

        stop(paste("File", hic, "does not exist."), call. = FALSE)

      }

      ## if it's a valid .hic file, it needs to have a valid norm parameter
      if (is.null(norm)){

        stop("If providing .hic file, please specify \'norm\'.", call. = FALSE)

      }

    }

    ##### chrom/chromstart/chromend #####

    if (is.null(hic_plot$chrom)){

      stop("Please specify \'chrom\'.", call. = FALSE)

    }

    ## Can't have only one NULL chromstart or chromend
    if (any(is.null(hic_plot$chromstart), is.null(hic_plot$chromend))){

      stop("Please specify \'chromstart\' and \'chromend\'.", call. = FALSE)

    }

    ## chromstart should be smaller than chromend
    if (hic_plot$chromstart > hic_plot$chromend){

      stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)

    }

    ##### resolution #####

    if (is.null(hic_plot$resolution)){

      stop("Invalid \'resolution\' value.  Options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000.", call. = FALSE)

    } else {

      if (!(hic_plot$resolution %in% c(2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000))){

        stop("Invalid \'resolution\' value.  Options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000.", call. = FALSE)

      }

    }

    ##### zrange #####

    ## Ensure properly formatted zrange
    if (!is.null(hic_plot$zrange)){

      ## zrange needs to be a vector
      if (!is.vector(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be length 2
      if (length(hic_plot$zrange) != 2){

        stop("\'zrange\' must be a vector of length 2.", call. = FALSE)

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of two numbers.", call. = FALSE)

      }

      ## second value should be larger than the first value
      if (hic_plot$zrange[1] >= hic_plot$zrange[2]){

        stop("\'zrange\' must be a vector of two numbers in which the 2nd value is larger than the 1st.", call. = FALSE)

      }

    }


    ##### height #####
    if (!is.null(hic_plot$height)){

      ## convert height to inches
      height <- convertHeight(hic_plot$height, unitTo = "inches", valueOnly = T)
      if (height < 0.05){
        stop("Height is too small for a valid triangle Hi-C plot.", call. = FALSE)

      }

    }

  }

  ## Define a function to check range of data in dataframe
  check_dataframe <- function(hic, hic_plot){

    if (min(hic[,1]) > hic_plot$chromstart | max(hic[,1]) < hic_plot$chromend | min(hic[,2]) > hic_plot$chromstart | max(hic[,2]) < hic_plot$chromend){

      warning("Data is incomplete for the specified range.", call. = FALSE)

    }

  }

  ## Define a function that reads in hic data for bb_plothic
  read_data <- function(hic, hic_plot, norm){

    ## if .hic file, read in with bb_rhic
    if (!("data.frame" %in% class(hic))){

      message(paste("Reading in hic file with", norm, "normalization."))

      readchromstart <- hic_plot$chromstart - hic_plot$resolution
      readchromend <- hic_plot$chromend + hic_plot$resolution

      hic <- bb_readHic(hic = hic, chrom = hic_plot$chrom, chromstart = readchromstart, chromend = readchromend,
                        resolution = hic_plot$resolution, zrange = hic_plot$zrange, norm = norm)

    } else {

      message("Reading in dataframe.")

      ## check range of data in dataframe
      check_dataframe(hic = hic, hic_plot = hic_plot)

    }

    ## Rename columns for later processing
    colnames(hic) <- c("x", "y", "counts")

    return(hic)

  }

  ## Define a function that subsets data
  subset_data <- function(hic, hic_plot){

    hic <- hic[which(hic[,1] >= hic_plot$chromstart - hic_plot$resolution &
                       hic[,1] < hic_plot$chromend &
                       hic[,2] >= hic_plot$chromstart - hic_plot$resolution &
                       hic[,2] < hic_plot$chromend),]

    return(hic)
  }

  ## Define a function that sets the zrange
  set_zrange <- function(hic, hic_plot){

    ## no zrange, only one value
    if (is.null(hic_plot$zrange) & length(unique(hic$counts)) == 1){

      zrange <- c(unique(hic$counts), unique(hic$counts))
      hic_plot$zrange <- zrange

    }

    ## no zrange, multiple values
    if (is.null(hic_plot$zrange) & length(unique(hic$counts)) > 1){

      zrange <- c(0, max(hic$counts))
      hic_plot$zrange <- zrange

    }

    return(hic_plot)

  }

  ## Define a function that converts the location to the bottom left of the triangle based on justification
  convert_just <- function(hic_plot){

    height <- convertHeight(hic_plot$height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    width <- convertWidth(hic_plot$width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    x <- convertX(hic_plot$x, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    y <- get("page_height", envir = bbEnv) - convertY(hic_plot$y, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    just <- hic_plot$justification

    ## Calculate height of triangle/trapezoid
    two <- mpfr(2, 120)
    desired_height <- height
    calc_height <- width*0.5
    side_length <- width/sqrt(two)

    if (calc_height > desired_height){
    ## here we'll have a trapezoid
      trap_top <- 2*(calc_height - desired_height)

      #height_diff <- calc_height - desired_height

      if (length(just) == 2){

        if (identical(just, c("left", "bottom"))){
          new_x <- x
          new_y <- y
        } else if (identical(just, c("right", "bottom"))){
          new_x <- x - width
          new_y <- y
        } else if (identical(just, c("left", "center"))){
          new_x <- x - (0.25*(width - trap_top))
          new_y <- y - (0.5*desired_height)
        } else if (identical(just, c("right", "center"))){
          new_x <- x - ((0.75*(width - trap_top)) + trap_top)
          new_y <- y - (0.5*desired_height)
        } else if (identical(just, c("center", "bottom"))){
          new_x <- x - ((0.5*trap_top) + (0.5*(width - trap_top)))
          new_y <- y
        } else if (identical(just, c("center", "top"))){
          new_x <- x - ((0.5*trap_top) + (0.5*(width - trap_top)))
          new_y <- y - desired_height
        } else if (identical(just, c("left", "top"))){
          new_x <- x - (0.5*(width - trap_top))
          new_y <- y - desired_height
        } else if (identical(just, c("right", "top"))){
          new_x <- x - ((0.5*(width - trap_top)) + trap_top)
          new_y <- y - desired_height
        } else {
          new_x <-  x - ((0.5*trap_top) + (0.5*(width - trap_top)))
          new_y <- y - (0.5*desired_height)
        }

      } else if (length(just) == 1){

        if (just == "left"){
          new_x <- x - (0.25*(width - trap_top))
          new_y <- y - (0.5*desired_height)
        } else if (just == "right"){
          new_x <- x - ((0.75*(width - trap_top)) + trap_top)
          new_y <- y - (0.5*desired_height)
        } else if (just == "bottom"){
          new_x <- x - ((0.5*trap_top) + (0.5*(width - trap_top)))
          new_y <- y
        } else if (just == "top"){
          new_x <- x - ((0.5*trap_top) + (0.5*(width - trap_top)))
          new_y <- y - desired_height
        } else {
          new_x <-  x - ((0.5*trap_top) + (0.5*(width - trap_top)))
          new_y <- y - (0.5*desired_height)
        }
      }

    } else {
      ## here we'll just have the triangle
      if (length(just) == 2){

        if (identical(just, c("left", "bottom"))){
          new_x <- x
          new_y <- y
        } else if (identical(just, c("right", "bottom"))){
          new_x <- x - width
          new_y <- y
        } else if (identical(just, c("left", "center"))){
          new_x <- x - (0.25*width)
          new_y <- y - (0.5*calc_height)
        } else if (identical(just, c("right", "center"))){
          new_x <- x - (0.75*width)
          new_y <- y - (0.5*calc_height)
        } else if (identical(just, c("center", "bottom"))){
          new_x <- x - (0.5*width)
          new_y <- y
        } else if (identical(just, c("center", "top"))){
          new_x <- x - (0.5*width)
          new_y <- y - calc_height
        } else {
          new_x <- x - (0.5*width)
          new_y <- y - (0.5*calc_height)
        }

      } else if (length(just) == 1){

        if (just == "left"){
          new_x <- x - (0.25*width)
          new_y <- y - (0.5*calc_height)
        } else if (just == "right"){
          new_x <- x - (0.75*width)
          new_y <- y - (0.5*calc_height)
        } else if (just == "bottom"){
          new_x <- x - (0.5*width)
          new_y <- y
        } else if (just == "top"){
          new_x <- x - (0.5*width)
          new_y <- y - calc_height
        } else {
          new_x <- x - (0.5*width)
          new_y <- y - (0.5*calc_height)
        }

      }

    }

    return(list(new_x, new_y))

  }

  ## Define a function that manually "clips" squares/triangles along edges
  manual_clip <- function(hic, hic_plot){

    ## make sure appropriate coordinates are in both the sides of the hic data (i.e. no missing pixels)
    all_coords <- seq(floor(hic_plot$chromstart/hic_plot$resolution)*hic_plot$resolution,
                      floor((hic_plot$chromend - hic_plot$resolution)/hic_plot$resolution)*hic_plot$resolution, hic_plot$resolution)

    hic_side <- hic[which(hic[,1] == min(hic[,1])),]
    hic_top <- hic[which(hic[,2] == max(hic[,2])),]

    hic <- hic[which(hic[,1] != min(hic[,1])),]
    hic <- hic[which(hic[,2] != max(hic[,2])),]

    side_missing <- dplyr::setdiff(all_coords, hic_side$y)
    top_missing <- dplyr::setdiff(all_coords, hic_top$x)

    side_missing_comp <- rep(hic_side[1,1], length(side_missing))
    top_missing_comp <- rep(hic_top[1,2], length(top_missing))

    add_side <- data.frame("x" = side_missing_comp, "y" = side_missing, "counts" = rep(NA, length(side_missing)),
                           "color" = rep(NA, length(side_missing)),
                           "width" = rep(hic_plot$resolution, length(side_missing)), "height" = rep(hic_plot$resolution, length(side_missing)))

    add_top <- data.frame("x" = top_missing, "y" = top_missing_comp, "counts" = rep(NA, length(top_missing)),
                          "color" = rep(NA, length(top_missing)),
                          "width" = rep(hic_plot$resolution, length(top_missing)), "height" = rep(hic_plot$resolution, length(top_missing)))

    sideTotal <- rbind(hic_side, add_side)
    sideTotal <- sideTotal[-which(sideTotal[,1] == min(sideTotal[,1]) & sideTotal[,2] == max(sideTotal[,2])),]
    topTotal <- rbind(hic_top, add_top)


    hic <- rbind(hic, sideTotal, topTotal)
    hic <- hic[order(as.numeric(rownames(hic))),]

    ## lowest most genomic coordinate of data set
    left_min <- min(hic[,1])
    ## highest most genomic coordinate of data set
    top_max <- max(hic[,2]) + hic_plot$resolution

    ############# Squares
    squares <- hic[which(hic[,2] > hic[,1]),]

    top_left <- squares[which(squares[,1] == min(squares[,1]) & squares[,2] == max(squares[,2])),]
    squares <- squares[-which(squares[,1] == min(squares[,1]) & squares[,2] == max(squares[,2])),]

    left_squares <- squares[which(squares[,1] == min(squares[,1])),]
    squares <- squares[-which(squares[,1] == min(squares[,1])),]

    top_squares <- squares[which(squares[,2] == max(squares[,2])),]
    squares <- squares[-which(squares[,2] == max(squares[,2])),]

    ############# Triangles
    triangles <- hic[which(hic[,2] == hic[,1]),]

    bottom_left <- triangles[which(triangles[,1] == min(triangles[,1])),]
    triangles <- triangles[-which(triangles[,1] == min(triangles[,1])),]

    top_right <- triangles[which(triangles[,2] == max(triangles[,2])),]
    triangles <- triangles[-which(triangles[,2] == max(triangles[,2])),]

    if (left_min < hic_plot$chromstart & top_max > hic_plot$chromend){

      new_width <- hic_plot$resolution - (hic_plot$chromstart - left_min)
      new_height <- hic_plot$resolution - (top_max - hic_plot$chromend)

      ## Adjust top left square from both left and top
      top_left$x <- chromstart
      top_left$width <- new_width
      top_left$height <- new_height

      ## Adjust left squares from left
      left_squares$x <- chromstart
      left_squares$width <- new_width

      ## Adjust top squares from top
      top_squares$height <- new_height

      ## Adjust bottom left triangle
      bottom_left$x <- chromstart
      bottom_left$y <- chromstart
      bottom_left$width <- new_width
      bottom_left$height <- new_width

      ## Adjust top right triangle
      top_right$height <- new_height
      top_right$width <- new_height

    } else {

      if (left_min < hic_plot$chromstart){

        new_width <- hic_plot$resolution - (hic_plot$chromstart - left_min)

        ## Adjust top left square from left only
        top_left$x <- chromstart
        top_left$width <- new_width

        ## Adjust left squares from left
        left_squares$x <- chromstart
        left_squares$width <- new_width

        ## Adjust bottom left triangle
        bottom_left$x <- chromstart
        bottom_left$y <- chromstart
        bottom_left$width <- new_width
        bottom_left$height <- new_width

      }

      else if (top_max > hic_plot$chromend){

        new_height <- hic_plot$resolution - (top_max - hic_plot$chromend)

        ## Adjust top left square from top only
        top_left$height <- new_height

        ## Adjust top squares from top
        top_squares$height <- new_height

        ## Adjust top triangle
        top_right$height <- new_height
        top_right$width <- new_height

      }

    }

    ## Recombine squares
    all_squares <- rbind(squares, top_left, left_squares, top_squares)

    ## Recombine triangles
    all_triangles <- rbind(triangles, top_right, bottom_left)

    ##vRecombine everything
    total <- rbind(all_squares, all_triangles)

    return(total)

  }

  ## Define a function that will manually "crop" the point of the triangle plot if the height is too large
  manual_cropTop <- function(hic, hic_plot){

    ## function that subsets regions of matrices
    matrix_data <- function(hic, data, upper_tri, lower_tri){

      if (nrow(data) <= 1){

        return(data)

      } else {

        cast <- as.matrix(reshape::cast(data, formula = y ~ x, value = "counts"))

        if (upper_tri == TRUE){

          cast[upper.tri(cast)] <- NA
        }

        if (lower_tri == TRUE){

          cast[lower.tri(cast)] <- NA
        }

        chop <- setNames(reshape2::melt(cast), c("y", "x", "counts"))
        chop <- chop[c(2,1,3)]
        chop <- subset(chop, !is.na(counts))
        final <- merge(chop, hic)

        return(final)

      }


    }

    ## function that makes shapes for bottom diagonal
    hic_bottomDiagonal <- function(diag, side_difference){

      x <- as.numeric(diag[1])
      y <- as.numeric(diag[2])
      width <- as.numeric(diag[5])
      height <- as.numeric(diag[6])

      x1 <- x
      x2 <- x + width - side_difference
      x3 <- x + width
      x4 <- x3
      x5 <- x1

      y1 <- y + side_difference
      y2 <- y + height
      y3 <- y2
      y4 <- y
      y5 <- y4

      col <- diag[4]

      hic_pentagon <- polygonGrob(x = c(x1, x2, x3, x4, x5),
                                  y = c(y1, y2, y3, y4, y5),
                                  gp = gpar(col = NA, fill = col),
                                  default.units = "native")

      assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_pentagon), envir = bbEnv)


    }

    ## function that makes shapes for top diagonal
    hic_topDiagonal <- function(diag, side_difference){

      x <- as.numeric(diag[1])
      y <- as.numeric(diag[2])
      width <- as.numeric(diag[5])
      height <- as.numeric(diag[6])

      x1 <- x + width - side_difference
      x2 <- x + width
      x3 <- x2

      y1 <- y
      y2 <- y + side_difference
      y3 <- y

      col <- diag[4]

      hic_triangle <- polygonGrob(x = c(x1, x2, x3),
                                  y = c(y1, y2, y3),
                                  gp = gpar(col = NA, fill = col),
                                  default.units = "native")

      assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_triangle), envir = bbEnv)

    }

    if (!is.null(hic_plot$x) & !is.null(hic_plot$y)){

      two <- mpfr(2, 120)
      desired_height <- convertHeight(hic_plot$height, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
      calc_height <- convertWidth(hic_plot$width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)*0.5
      width <- convertWidth(hic_plot$width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
      side_length <- width/sqrt(two)

      ## here we need to chop pixels out
      if (calc_height > desired_height){

        ## Extract the pixels that fall along the sides
        hic_side <- hic[which(hic[,1] == min(hic[,1])),]
        hic_top <- hic[which(hic[,2] == max(hic[,2])),]

        ## This is the length (in page units) to chop
        side_chop <- (calc_height - desired_height) * sqrt(two)

        ## Get the top left and bottom left to take pixel cropping into account
        topPixel <- hic[which(hic[,1] == min(hic[,1]) & hic[,2] == max(hic[,2])),]
        bottomLeft <- hic_side[which(hic_side[,2] == min(hic_side[,2])),]

        ## Subtract the top left and bottom left pixel from the other hic_side pixels
        sidePixels <- hic_side[-which(hic_side[,1] == min(hic_side[,1]) & hic_side[,2] == max(hic_side[,2])),]
        sidePixels <- sidePixels[-which(sidePixels[,2] == min(sidePixels[,2])),]

        ## Get the total number of pixels on a side (both sides should be the same)
        pixelNo <- nrow(sidePixels) + (topPixel$height/hic_plot$resolution) + (bottomLeft$height/hic_plot$resolution)

        ## This is length of a full pixel in page units
        pixelLength <- side_length/pixelNo

        ## This is the number of pixels to chop (being overconservative)
        pixelsChop <- side_chop/pixelLength - min(topPixel$height/hic_plot$resolution, topPixel$width/hic_plot$resolution)

        chop_data_side <- tail(hic_side, n = as.numeric(ceiling(pixelsChop)) + 1)
        chop_data_top <- head(hic_top, n = as.numeric(ceiling(pixelsChop)) + 1)

        ## Leftover pixel distances:
        ## coming from the bottom of a full pixel
        bottom_difference <- hic_plot$resolution - (as.numeric(pixelsChop) - as.numeric(floor(pixelsChop)))*hic_plot$resolution
        ## coming from the top of a full pixel
        top_difference <- (as.numeric(pixelsChop) - as.numeric(floor(pixelsChop)))*hic_plot$resolution

        ## Getting coordinates of square region of where to chop
        x_coords <- chop_data_top$x
        y_coords <- chop_data_side$y
        x_coords2 <- x_coords[-length(x_coords)]
        y_coords2 <- y_coords[-1]

        # Get 2 square regions of where to chop (one for larger diagonal and other for inner diagonal)
        square_chop <- hic[which(hic[,1] %in% x_coords & hic[,2] %in% y_coords),]
        square_chop2 <- hic[which(hic[,1] %in% x_coords2 & hic[,2] %in% y_coords2),]

        ## Get entire triangular region to eliminate
        removed <- matrix_data(hic = hic, data = square_chop, upper_tri = TRUE, lower_tri = FALSE)

        ## Get main diagonal to make pentagons
        bdiagonal <- matrix_data(hic = hic, data = square_chop, upper_tri = TRUE, lower_tri = TRUE)

        ## Get slightly higher diagonal to make tiny triangles
        tdiagonal <- matrix_data(hic = hic, data = square_chop2, upper_tri = TRUE, lower_tri = TRUE)

        ## BOTTOM DIAGONAL PENTAGONS
        ## Extract first and last pixels in case of pixel cropping

        ## First
        if (nrow(bdiagonal) > 1){

          bdiagonal_first <- bdiagonal[which(bdiagonal$x == min(bdiagonal$x) & bdiagonal$y == min(bdiagonal$y)),]
          if (bdiagonal_first$x >= (bdiagonal_first$x + bdiagonal_first$width - bottom_difference)){
            hic_pentagon1 <- polygonGrob(x = c(bdiagonal_first$x, bdiagonal_first$x,
                                               bdiagonal_first$x + bdiagonal_first$width, bdiagonal_first$x + bdiagonal_first$width),
                                         y = c(bdiagonal_first$y, bdiagonal_first$y + bdiagonal_first$height,
                                               bdiagonal_first$y + bdiagonal_first$height, bdiagonal_first$y),
                                         gp = gpar(col = NA, fill = bdiagonal_first$color),
                                         default.units = "native")
          } else {
            hic_pentagon1 <- polygonGrob(x = c(bdiagonal_first$x, bdiagonal_first$x + bdiagonal_first$width - bottom_difference,
                                               bdiagonal_first$x + bdiagonal_first$width, bdiagonal_first$x + bdiagonal_first$width,
                                               bdiagonal_first$x),
                                         y = c(bdiagonal_first$y + bdiagonal_first$height - (bdiagonal_first$width - bottom_difference),
                                               bdiagonal_first$y + bdiagonal_first$height,
                                               bdiagonal_first$y + bdiagonal_first$height, bdiagonal_first$y, bdiagonal_first$y),
                                         gp = gpar(col = NA, fill = bdiagonal_first$color),
                                         default.units = "native")
          }
          assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_pentagon1), envir = bbEnv)

          ## Last
          bdiagonal_last <- bdiagonal[which(bdiagonal$x == max(bdiagonal$x) & bdiagonal$y == max(bdiagonal$y)),]
          if ((bdiagonal_last$y + bdiagonal_last$height) <= (bdiagonal_last$y + bottom_difference)){
            hic_pentagon2 <- polygonGrob(x = c(bdiagonal_last$x, bdiagonal_last$x, bdiagonal_last$x + bdiagonal_last$width,
                                               bdiagonal_last$x + bdiagonal_last$width),
                                         y = c(bdiagonal_last$y, bdiagonal_last$y + bdiagonal_last$height, bdiagonal_last$y + bdiagonal_last$height,
                                               bdiagonal_last$y),
                                         gp = gpar(col = NA, fill = bdiagonal_last$color),
                                         default.units = "native")
          } else {
            hic_pentagon2 <- polygonGrob(x = c(bdiagonal_last$x, bdiagonal_last$x + bdiagonal_last$height - bottom_difference, bdiagonal_last$x + bdiagonal_last$width,
                                               bdiagonal_last$x + bdiagonal_last$width, bdiagonal_last$x),
                                         y = c(bdiagonal_last$y + bottom_difference, bdiagonal_last$y + bdiagonal_last$height, bdiagonal_last$y + bdiagonal_last$height,
                                               bdiagonal_last$y, bdiagonal_last$y),
                                         gp = gpar(col = NA, fill = bdiagonal_last$color),
                                         default.units = "native")
          }
          assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_pentagon2), envir = bbEnv)

          ## Remove first and last pixel
          bdiagonal <- bdiagonal[-which(bdiagonal$x == min(bdiagonal$x)),]
          bdiagonal <- bdiagonal[-which(bdiagonal$x == max(bdiagonal$x)),]

          if (nrow(bdiagonal) >= 1){

            invisible(apply(bdiagonal, 1, hic_bottomDiagonal, side_difference = bottom_difference))

          }

        } else if (nrow(bdiagonal) == 1){

          hic_pentagon <- polygonGrob(x = c(bdiagonal$x, bdiagonal$x, bdiagonal$x + bdiagonal$height - bottom_difference, bdiagonal$x + bdiagonal$width, bdiagonal$x + bdiagonal$width),
                                      y = c(bdiagonal$y, bdiagonal$y + bottom_difference, bdiagonal$y + bdiagonal$height, bdiagonal$y + bdiagonal$height, bdiagonal$y),
                                      gp = gpar(col = NA, fill = bdiagonal$color),
                                      default.units = "native")
          assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_pentagon), envir = bbEnv)

        }

        ## TOP DIAGONAL TRIANGLES
        ## Extract first and last pixels in case of pixel cropping
        if (nrow(tdiagonal) > 1){

          tdiagonal_first <- tdiagonal[which(tdiagonal$x == min(tdiagonal$x) & tdiagonal$y == min(tdiagonal$y)),]
          tdiagonal_last <- tdiagonal[which(tdiagonal$x == max(tdiagonal$x) & tdiagonal$y == max(tdiagonal$y)),]

          if (tdiagonal_first$x > (tdiagonal_first$x + tdiagonal_first$width - bottom_difference)){
            hic_triangle1 <- polygonGrob(x = c(tdiagonal_first$x, tdiagonal_first$x, tdiagonal_first$x + tdiagonal_first$width,
                                               tdiagonal_first$x + tdiagonal_first$width),
                                         y = c(tdiagonal_first$y, tdiagonal_first$y + bottom_difference - tdiagonal_first$width,
                                               tdiagonal_first$y + bottom_difference, tdiagonal_first$y),
                                         gp = gpar(col = NA, fill = tdiagonal_first$color),
                                         default.units = "native")
          } else {
            hic_triangle1 <- polygonGrob(x = c(tdiagonal_first$x + tdiagonal_first$width - bottom_difference,
                                               tdiagonal_first$x + tdiagonal_first$width, tdiagonal_first$x + tdiagonal_first$width),
                                         y = c(tdiagonal_first$y, tdiagonal_first$y + bottom_difference, tdiagonal_first$y),
                                         gp = gpar(col = NA, fill = tdiagonal_first$color),
                                         default.units = "native")
          }
          assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_triangle1), envir = bbEnv)
          if ((tdiagonal_last$y + tdiagonal_last$height) < (tdiagonal_last$y + bottom_difference)){
            hic_triangle2 <- polygonGrob(x = c(tdiagonal_last$x + tdiagonal_last$width - bottom_difference,
                                               tdiagonal_last$x + tdiagonal_last$width - (bottom_difference - tdiagonal_last$height),
                                               tdiagonal_last$x + tdiagonal_last$width,
                                               tdiagonal_last$x + tdiagonal_last$width),
                                         y = c(tdiagonal_last$y, tdiagonal_last$y + tdiagonal_last$height,
                                               tdiagonal_last$y + tdiagonal_last$height, tdiagonal_last$y),
                                         gp = gpar(col = NA, fill = tdiagonal_last$color),
                                         default.units = "native")
          } else {
            hic_triangle2 <- polygonGrob(x = c(tdiagonal_last$x + tdiagonal_last$width - bottom_difference,
                                               tdiagonal_last$x + tdiagonal_last$width,
                                               tdiagonal_last$x + tdiagonal_last$width),
                                         y = c(tdiagonal_last$y, tdiagonal_last$y + bottom_difference,
                                               tdiagonal_last$y),
                                         gp = gpar(col = NA, fill = tdiagonal_last$color),
                                         default.units = "native")
          }
          assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_triangle2), envir = bbEnv)

          ## Remove first and last pixel
          tdiagonal <- tdiagonal[-which(tdiagonal$x == min(tdiagonal$x)),]
          tdiagonal <- tdiagonal[-which(tdiagonal$x == max(tdiagonal$x)),]

          if (nrow(tdiagonal) >= 1){
            ## Triangles for the rest of the pixels
            invisible(apply(tdiagonal, 1, hic_topDiagonal, side_difference = bottom_difference))

          }

        } else if (nrow(tdiagonal) == 1) {

            if (((tdiagonal$y + tdiagonal$height) < (tdiagonal$y + bottom_difference)) & (tdiagonal$x > (tdiagonal$x + tdiagonal$width - bottom_difference))){

              hic_triangle <- polygonGrob(x = c(tdiagonal$x, tdiagonal$x, tdiagonal$x + tdiagonal$width - (bottom_difference - tdiagonal$height),
                                                tdiagonal$x + tdiagonal$width, tdiagonal$x + tdiagonal$width),
                                          y = c(tdiagonal$y, tdiagonal$y + bottom_difference - tdiagonal$width, tdiagonal$y + tdiagonal$height,
                                                tdiagonal$y + tdiagonal$height, tdiagonal$y),
                                          gp = gpar(col = NA, fill = tdiagonal$color),
                                          default.units = "native")
            } else {

              if ((tdiagonal$y + tdiagonal$height) < (tdiagonal$y + bottom_difference)){
                ## right side pixel chopping
                hic_triangle <- polygonGrob(x = c(tdiagonal$x + tdiagonal$width - bottom_difference,
                                                  tdiagonal$x + tdiagonal$width - (bottom_difference - tdiagonal$height),
                                                  tdiagonal$x + tdiagonal$width,
                                                  tdiagonal$x + tdiagonal$width),
                                            y = c(tdiagonal$y, tdiagonal$y + tdiagonal$height,
                                                  tdiagonal$y + tdiagonal$height, tdiagonal$y),
                                            gp = gpar(col = NA, fill = tdiagonal$color),
                                            default.units = "native")
              } else if (tdiagonal$x > (tdiagonal$x + tdiagonal$width - bottom_difference)){
                ## left side pixel chopping
                hic_triangle <- polygonGrob(x = c(tdiagonal$x, tdiagonal$x, tdiagonal$x + tdiagonal$width,
                                                  tdiagonal$x + tdiagonal$width),
                                            y = c(tdiagonal$y, tdiagonal$y + bottom_difference - tdiagonal$width,
                                                  tdiagonal$y + bottom_difference, tdiagonal$y),
                                            gp = gpar(col = NA, fill = tdiagonal$color),
                                            default.units = "native")

              } else {

                hic_triangle <- polygonGrob(x = c(tdiagonal$x + tdiagonal$width - bottom_difference, tdiagonal$x + tdiagonal$width, tdiagonal$x + tdiagonal$width),
                                            y = c(tdiagonal$y, tdiagonal$y + bottom_difference, tdiagonal$y),
                                            gp = gpar(col = NA, fill = tdiagonal$color),
                                            default.units = "native")

              }

            }

          assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_triangle), envir = bbEnv)

        }

        ## Remove pixels to be completely deleted or that have already been plotted with other shapes
        hic_removed <- dplyr::setdiff(hic, removed)

        return(list(hic_removed, as.numeric(pixelsChop)))

        } else if (calc_height < desired_height){

        units <- gsub(pattern = "[0-9]|\\.", replacement = "", as.character(hic_plot$width))
        converted_height <- convertHeight(unit(calc_height, get("page_units", envir = bbEnv)), unitTo = units, valueOnly = T)
        pixelsChop = 0

        warning(paste0("Specified height is larger than the height of a right triangle with the specified width. Adjusting height to ", converted_height, " ", units, "."), call. = FALSE)
        return(list(hic, pixelsChop))

      } else {

        pixelsChop = 0
        return(list(hic, pixelsChop))
      }

    } else {
      pixelsChop = 0
      return(list(hic, pixelsChop))

    }

  }

  ## Define a function that makes grobs for the hic diagonal
  hic_diagonal <- function(hic){

    col <- hic[4]
    x <- as.numeric(hic[1])
    y <- as.numeric(hic[2])
    width <- as.numeric(hic[5])
    height <- as.numeric(hic[6])

    xleft = x
    xright = x + width
    ybottom = y
    ytop = y + height

    hic_triangle <- polygonGrob(x = c(xleft, xleft, xright),
                                y = c(ybottom, ytop, ytop),
                                gp = gpar(col = NA, fill = col),
                                default.units = "native")

    assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_triangle), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  hic_plot <- structure(list(chrom = chrom, chromstart = as.numeric(chromstart), chromend = as.numeric(chromend), x = x, y = y, width = width, height = height, justification = NULL,
                             zrange = zrange, altchrom = chrom, altchromstart = chromstart, altchromend = chromend, resolution = resolution,
                             color_palette = NULL, grobs = NULL), class = "bb_trianglehic")
  attr(x = hic_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CHECK PLACEMENT
  # ======================================================================================================================================================================================

  check_placement(object = hic_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  hic_plot <- defaultUnits(object = hic_plot, default.units = default.units)

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotTriangleHic(hic = hic, hic_plot = hic_plot, norm = norm)

  # ======================================================================================================================================================================================
  # JUSTIFICATION OF PLOT
  # ======================================================================================================================================================================================

  new_just <- reset_just(just = just, x = hic_plot$x, y = hic_plot$y, width = hic_plot$width, height = hic_plot$height)
  hic_plot$justification <- new_just

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  hic <- read_data(hic = hic, hic_plot = hic_plot, norm = norm)

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  hic <- subset_data(hic = hic, hic_plot = hic_plot)

  # ======================================================================================================================================================================================
  # SET ZRANGE AND SCALE DATA
  # ======================================================================================================================================================================================

  hic_plot <- set_zrange(hic = hic, hic_plot = hic_plot)
  hic$counts[hic$counts <= hic_plot$zrange[1]] <- hic_plot$zrange[1]
  hic$counts[hic$counts >= hic_plot$zrange[2]] <- hic_plot$zrange[2]

  # ======================================================================================================================================================================================
  # CONVERT NUMBERS TO COLORS
  # ======================================================================================================================================================================================

  ## if we don't have an appropriate zrange (even after setting it based on a null zrange), can't scale to colors
  if (!is.null(hic_plot$zrange) & length(unique(hic_plot$zrange)) == 2){

    hic$color <- bb_maptocolors(hic$counts, col = palette, num = 100, range = hic_plot$zrange)
    sorted_colors <- unique(hic[order(hic$counts),]$color)
    hic_plot$color_palette <- sorted_colors

  } else {

    ## If we still have a null zrange or a length(unique(zrange)) == 1, means we couldn't do it in setzrange above (empty data or data with only 1 value)
    warning("Can't scale data to colors.", call. = FALSE)

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_trianglehic", length(grep(pattern = "bb_trianglehic", x = currentViewports)) + 1)

  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(1/sqrt(two), "snpc"), width = unit(1/sqrt(2), "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.25, "npc"),
                   xscale = c(chromstart, chromend), yscale = c(chromstart, chromend),
                   just = "center",
                   name = vp_name,
                   angle = -45)

    if (draw == TRUE){

      vp$name <- "bb_trianglehic1"
      grid.newpage()

    }

  } else {

    ## Get sides of viewport based on input width
    vp_side <- (convertWidth(hic_plot$width, unitTo = get("page_units", envir = bbEnv), valueOnly = T))/sqrt(two)

    ## Get bottom left point of triangle (hence bottom left of actual viewport) based on just
    bottom_coords <- convert_just(hic_plot = hic_plot)

    ## Make viewport
    vp <- viewport(height = unit(vp_side, get("page_units", envir = bbEnv)), width = unit(vp_side, get("page_units", envir = bbEnv)),
                   x = unit(bottom_coords[[1]], get("page_units", envir = bbEnv)),
                   y = unit(bottom_coords[[2]], get("page_units", envir = bbEnv)),
                   xscale = c(chromstart, chromend), yscale = c(chromstart, chromend),
                   just = c("left", "bottom"),
                   name = vp_name,
                   angle = -45)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("hic_grobs2", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  hic$width <- resolution
  hic$height <- resolution

  ## Manually "clip" the grobs that fall out of the desired chromstart to chromend region
  hic_total <- manual_clip(hic = hic, hic_plot = hic_plot)
  hic_total <- hic_total[which(hic_total$width != 0 | hic_total$height != 0),]
  hic_total <- hic_total[order(as.numeric(rownames(hic_total))),]

  ## Manually crop the grobs into a trapezoid if the input height is larger than the calculated height for the triangle
  cropped <- manual_cropTop(hic = hic_total, hic_plot = hic_plot)
  hic_total2 <- cropped[[1]]
  attr(x = hic_plot, which = "choppedPixels") <- cropped[[2]]

  ## Separate into squares and triangles
  squares <- hic_total2[which(hic_total2[,2] > hic_total2[,1]),]
  triangles <- hic_total2[which(hic_total2[,2] == hic_total2[,1]),]

  if (nrow(squares) > 0){

    ## Make square grobs and add to grob gTree
    hic_squares <- rectGrob(x = squares$x,
                            y = squares$y,
                            just = c("left", "bottom"),
                            width = squares$width,
                            height = squares$height,
                            gp = gpar(col = NA, fill = squares$color),
                            default.units = "native")
    assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_squares), envir = bbEnv)
  }

  if (nrow(triangles) > 0){
    ## Make triangle grobs and add to grob gTree
    invisible(apply(triangles, 1, hic_diagonal))

  }

  # ======================================================================================================================================================================================
  # IF DRAW == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(get("hic_grobs2", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  hic_plot$grobs <- get("hic_grobs2", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(hic_plot)

}
