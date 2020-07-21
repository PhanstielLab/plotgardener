#' plots HiC interaction matrix in triangular format
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted, based on build (i.e. for hg19 just a number, for hg38 string like "chr8")
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param resolution the width in bp of each pixel; for hic files, "auto" will attempt to choose a resolution based on the size of the region; for
#' dataframes, "auto" will attempt to detect the resolution the dataframe contains
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param assembly desired genome assembly
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
bb_plotTriangleHic <- function(hic, chrom, chromstart = NULL, chromend = NULL, resolution = "auto", zrange = NULL,
                               palette = colorRampPalette(c("white", "dark red")), assembly = "hg19", width = NULL, height = NULL,
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


    ## Can't have only one NULL chromstart or chromend
    if ((is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)) | (is.null(hic_plot$chromend) & !is.null(hic_plot$chromstart))){

      stop("Cannot have one \'NULL\' \'chromstart\' or \'chromend\'.", call. = FALSE)

    }

    if (!is.null(hic_plot$chromstart) & !is.null(hic_plot$chromend)){
      ## Chromstart should be smaller than chromend
      if (hic_plot$chromstart > hic_plot$chromend){

        stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)

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

  ## Define a function to adjust/detect resolution based on .hic file/dataframe
  adjust_resolution <- function(hic, hic_plot){

    if (!("data.frame" %in% class(hic))){
      ## Get range of data and try to pick a resolution to extract from hic file
      dataRange <- hic_plot$chromend - hic_plot$chromstart
      if (dataRange >= 150000000){
        bestRes <- 500000
      } else if (dataRange >= 75000000 & dataRange < 150000000){
        bestRes <- 250000
      } else if (dataRange >= 35000000 & dataRange < 75000000){
        bestRes <- 100000
      } else if (dataRange >= 20000000 & dataRange < 35000000){
        bestRes <- 50000
      } else if (dataRange >= 5000000 & dataRange < 20000000){
        bestRes <- 25000
      } else if (dataRange >= 3000000 & dataRange < 5000000){
        bestRes <- 10000
      } else {
        bestRes <- 5000
      }

      hic_plot$resolution <- as.integer(bestRes)

    } else {

      ## Try to detect resolution from data
      offDiag <- hic[which(hic[,1] != hic[,2]),]
      bpDiffs <- abs(offDiag[,2] - offDiag[,1])
      predRes <- min(bpDiffs)

      hic_plot$resolution <- as.integer(predRes)

    }

    return(hic_plot)
  }

  ## Define a function that reads in hic data
  read_data <- function(hic, hic_plot, norm, assembly){

    parse_chrom <- function(assembly, chrom){

      if (assembly == "hg19"){
        strawChrom <- as.numeric(gsub("chr", "", chrom))
      }

      return(strawChrom)
    }

    ## if .hic file, read in with bb_rhic
    if (!("data.frame" %in% class(hic))){

      strawChrom <- parse_chrom(assembly = assembly, chrom = hic_plot$chrom)

      readchromstart <- hic_plot$chromstart - hic_plot$resolution
      readchromend <- hic_plot$chromend + hic_plot$resolution

      hic <- bb_readHic(hic = hic, chrom = strawChrom, chromstart = readchromstart, chromend = readchromend,
                        resolution = hic_plot$resolution, zrange = hic_plot$zrange, norm = norm)

    } else {

      message(paste("Read in dataframe.", hic_plot$resolution, "BP resolution detected."))
      ## check range of data in dataframe
      check_dataframe(hic = hic, hic_plot = hic_plot)

    }
    ## Rename columns for later processing
    colnames(hic) <- c("x", "y", "counts")
    hic <- na.omit(hic)

    return(hic)

  }

  ## Define a function that subsets data
  subset_data <- function(hic, hic_plot){

    hic <- hic[which(hic[,1] >= floor(hic_plot$chromstart/hic_plot$resolution)*hic_plot$resolution &
                       hic[,1] < hic_plot$chromend &
                       hic[,2] >= floor(hic_plot$chromstart/hic_plot$resolution)*hic_plot$resolution &
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
      top_left$x <- hic_plot$chromstart
      top_left$width <- new_width
      top_left$height <- new_height

      ## Adjust left squares from left
      left_squares$x <- hic_plot$chromstart
      left_squares$width <- new_width

      ## Adjust top squares from top
      top_squares$height <- new_height

      ## Adjust bottom left triangle
      bottom_left$x <- hic_plot$chromstart
      bottom_left$y <- hic_plot$chromstart
      bottom_left$width <- new_width
      bottom_left$height <- new_width

      ## Adjust top right triangle
      top_right$height <- new_height
      top_right$width <- new_height

    } else {

      if (left_min < hic_plot$chromstart){

        new_width <- hic_plot$resolution - (hic_plot$chromstart - left_min)

        ## Adjust top left square from left only
        top_left$x <- hic_plot$chromstart
        top_left$width <- new_width

        ## Adjust left squares from left
        left_squares$x <- hic_plot$chromstart
        left_squares$width <- new_width

        ## Adjust bottom left triangle
        bottom_left$x <- hic_plot$chromstart
        bottom_left$y <- hic_plot$chromstart
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

  hic_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, x = x, y = y, width = width, height = height, justification = NULL,
                             zrange = zrange, altchrom = chrom, altchromstart = chromstart, altchromend = chromend, resolution = resolution,
                             color_palette = NULL, grobs = NULL, assembly = assembly), class = "bb_trianglehic")
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
  # WHOLE CHROM
  # ======================================================================================================================================================================================
  if (is.null(chromstart) & is.null(chromend)){
    if (assembly == "hg19"){
      genome <- bb_hg19
    }

    hic_plot$chromstart <- 1
    hic_plot$chromend <- genome[which(genome$chrom == chrom),]$length
    hic_plot$altchromstart <- 1
    hic_plot$altchromend <- genome[which(genome$chrom == chrom),]$length

  }

  # ======================================================================================================================================================================================
  # ADJUST RESOLUTION
  # ======================================================================================================================================================================================

  if (resolution == "auto"){
    hic_plot <- adjust_resolution(hic = hic, hic_plot = hic_plot)
  }

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  hic <- read_data(hic = hic, hic_plot = hic_plot, norm = norm, assembly = assembly)

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

    }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_trianglehic", length(grep(pattern = "bb_trianglehic", x = currentViewports)) + 1)

  if (is.null(x) & is.null(y)){

    inside_vp <- viewport(height = unit(1, "npc"), width = unit(0.5, "npc"),
                          x = unit(0, "npc"), y = unit(0, "npc"),
                          xscale = c(hic_plot$chromstart, hic_plot$chromend),
                          yscale = c(hic_plot$chromstart, hic_plot$chromend),
                          just = c("left", "bottom"),
                          name = paste0(vp_name, "_inside"),
                          angle = -45)

    outside_vp <- viewport(height = unit(0.75, "snpc"),
                           width = unit(1.5, "snpc"),
                           x = unit(0.125, "npc"),
                           y = unit(0.25, "npc"),
                           xscale = c(hic_plot$chromstart, hic_plot$chromend),
                           clip = "on",
                           just = c("left", "bottom"),
                           name = paste0(vp_name, "_outside"))


    if (draw == TRUE){

      inside_vp$name <- "bb_trianglehic1_inside"
      outside_vp$name <- "bb_trianglehic1_outside"
      grid.newpage()

    }

  } else {

    ## Get sides of viewport based on input width
    vp_side <- (convertWidth(hic_plot$width, unitTo = get("page_units", envir = bbEnv), valueOnly = T))/sqrt(two)

    ## Get bottom left point of triangle (hence bottom left of actual viewport) based on just
    bottom_coords <- convert_just(hic_plot = hic_plot)

    inside_vp <- viewport(height = unit(vp_side, get("page_units", envir = bbEnv)), width = unit(vp_side, get("page_units", envir = bbEnv)),
                          x = unit(0, "npc"),
                          y = unit(0, "npc"),
                          xscale = c(hic_plot$chromstart, hic_plot$chromend),
                          yscale = c(hic_plot$chromstart, hic_plot$chromend),
                          just = c("left", "bottom"),
                          name = paste0(vp_name, "_inside"),
                          angle = -45)

    ## Convert coordinates into same units as page for outside vp
    page_coords <- convert_page(object = hic_plot)

    outside_vp <- viewport(height = page_coords$height,
                           width = page_coords$width,
                           x = unit(bottom_coords[[1]], get("page_units", envir = bbEnv)),
                           y = unit(bottom_coords[[2]], get("page_units", envir = bbEnv)),
                           xscale = c(hic_plot$chromstart, hic_plot$chromend),
                           clip = "on",
                           just = c("left", "bottom"),
                           name = paste0(vp_name, "_outside"))
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS
  # ======================================================================================================================================================================================

  assign("hic_grobs2", gTree(vp = inside_vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROBS
  # ======================================================================================================================================================================================

  hic$width <- hic_plot$resolution
  hic$height <- hic_plot$resolution

  ## Manually "clip" the grobs that fall out of the desired chromstart to chromend region
  hic <- manual_clip(hic = hic, hic_plot = hic_plot)
  hic <- hic[which(hic$width != 0 | hic$height != 0),]
  hic <- hic[order(as.numeric(rownames(hic))),]

  ## Separate into squares for upper region and triangle shapes for the diagonal
  squares <- hic[which(hic[,2] > hic[,1]),]
  triangles <- hic[which(hic[,2] == hic[,1]),]

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

  if (nrow(squares) == 0 & nrow(triangles) == 0){
    warning("Warning: no data found in region.  Suggestions: check chromosome, check region.", call. = FALSE)
  }
  # ======================================================================================================================================================================================
  # IF DRAW == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    pushViewport(outside_vp)
    grid.draw(get("hic_grobs2", envir = bbEnv))
    upViewport()
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
