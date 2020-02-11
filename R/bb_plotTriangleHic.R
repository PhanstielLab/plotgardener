#' plots HiC interaction matrix in triangular format
#'
#' @param hic path to .hic file or 3 column dataframe of counts
#' @param chrom chromosome of region to be plotted, based on build (i.e. for hg19 just a number, for hg38 string like "chr8")
#' @param chromstart chromosome start of region to be plotted
#' @param chromend chromosome end of region to be plotted
#' @param resolution the width in bp of each pixel; options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param width A unit object specifying the bottom width of the triangle
#' @param x A unit object specifying x-location
#' @param y A unit object specifying y-location
#' @param just a string or numeric vector specifying the justification of the viewport relative to its (x, y) location
#' @param draw A logical value indicating whether graphics output should be produced
#' @param norm if giving .hic file, hic data normalization; options are "NONE", "VC", "VC_SQRT", and "KR"
#'
#' @export
#'
#'
bb_plotTriangleHic <- function(hic, chrom, chromstart, chromend, resolution = 10000, zrange = NULL,
                               palette = colorRampPalette(c("white", "dark red")), width = NULL, x = NULL, y = NULL, just = "top", norm = "KR", draw = T, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## For more accurate calculation of sqrt(2)
  two <- mpfr(2, 120)

  ## Define a function that catches valid just values because there isn't a left top or right top option
  valid_just <- function(just){

    if (length(just == 2)) {

      if (identical(just, c("left", "top")) | identical(just, c("right", "top"))){

        stop("Please enter a valid just value.")

      }
    }

  }

  ## Define a function that catches errors for bb_plotTriangleHic
  errorcheck_bb_plotTriangleHic <- function(hic, hic_plot, resolution, norm){

    ###### hic/norm #####

    ## if it's a dataframe or datatable, it needs to be properly formatted
    if ("data.frame" %in% class(hic) && ncol(hic) != 3){

      stop("Invalid dataframe format.  Input a dataframe with 3 columns: chrA, chrB, counts.")

    }

    if (!"data.frame" %in% class(hic)){

      ## if it's a file path, it needs to be a .hic file
      if (file_ext(hic) != "hic"){

        stop("Invalid input. File must have a \".hic\" extension")

      }

      ## if it's a file path, it needs to exist
      if (!file.exists(hic)){

        stop(paste("File", hic, "does not exist."))

      }

      ## if it's a valid .hic file, it needs to have a valid norm parameter
      if (is.null(norm)){

        stop("If providing .hic file, please specify \'norm\'.  Options are \'NONE\', \'VC\', \'VC_SQRT\', or \'KR\'.")

      } else {

        if (!norm %in% c("NONE", "VC", "VC_SQRT", "KR")) {

          stop("If providing .hic file, please specify \'norm\'.  Options are \'NONE\', \'VC\', \'VC_SQRT\', or \'KR\'.")

        }

      }

    }

    ##### chrom/chromstart/chromend #####

    if (is.null(hic_plot$chrom)){

      stop("Please specify \'chrom\'.")

    }

    ## Can't have only one NULL chromstart or chromend
    if (any(is.null(hic_plot$chromstart), is.null(hic_plot$chromend))){

      stop("Please specify \'chromstart\' and \'chromend\'.")

    }

    ## chromstart should be smaller than chromend
    if (hic_plot$chromstart > hic_plot$chromend){

      stop("\'chromstart\' should not be larger than \'chromend\'.")

    }

    ##### resolution #####

    if (is.null(resolution)){

      stop("Invalid \'resolution\' value.  Options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000.")

    } else {

      if (!(resolution %in% c(2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000))){

        stop("Invalid \'resolution\' value.  Options are 2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, or 5000.")

      }

    }

    ##### zrange #####

    ## Ensure properly formatted zrange
    if (!is.null(hic_plot$zrange)){

      ## zrange needs to be a vector
      if (!is.vector(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of length 2.")

      }

      ## zrange vector needs to be length 2
      if (length(hic_plot$zrange) != 2){

        stop("\'zrange\' must be a vector of length 2.")

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(hic_plot$zrange)){

        stop("\'zrange\' must be a vector of two numbers.")

      }

      ## second value should be larger than the first value
      if (hic_plot$zrange[1] >= hic_plot$zrange[2]){

        stop("\'zrange\' must be a vector of two numbers in which the 2nd value is larger than the 1st.")

      }

    }

    ##### just #####

    valid_just(just = hic_plot$justification)

  }

  ## Define a function to check range of data in dataframe
  check_dataframe <- function(hic, hic_plot){

    if (min(hic[,1]) > hic_plot$chromstart | max(hic[,1]) < hic_plot$chromend | min(hic[,2]) > hic_plot$chromstart | max(hic[,2]) < hic_plot$chromend){

      warning("Data is incomplete for the specified range.")

    }

  }

  ## Define a function that reads in hic data for bb_plothic
  read_data <- function(hic, hic_plot, norm, resolution){

    ## if .hic file, read in with bb_rhic
    if (!("data.frame" %in% class(hic))){

      message(paste("Reading in hic file with", norm, "normalization."))

      readchromstart <- hic_plot$chromstart - resolution
      readchromend <- hic_plot$chromend + resolution

      hic <- bb_readHic(hic = hic, chrom = hic_plot$chrom, chromstart = readchromstart, chromend = readchromend,
                        resolution = resolution, zrange = hic_plot$zrange, norm = norm)

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
  subset_data <- function(hic, hic_plot, resolution){

    hic <- hic[which(hic[,1] >= hic_plot$chromstart - resolution &
                       hic[,1] <= hic_plot$chromend + resolution &
                       hic[,2] >= hic_plot$chromstart - resolution &
                       hic[,2] <= hic_plot$chromend + resolution),]

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
  convert_just <- function(width, x, y, just){
    ## Convert to the same unit and get values only
    width <- convertWidth(x = width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    x <- convertX(x = x, unitTo = get("page_units", envir = bbEnv), valueOnly = T)
    y <- get("page_height", envir = bbEnv) - convertY(x = y, unitTo = get("page_units", envir = bbEnv), valueOnly = T)

    ## Calculate height of triangle
    height <- 0.5*width

    if (length(just) == 2){

      if (identical(just, c("left", "bottom"))){
        new_x <- x
        new_y <- y
      } else if (identical(just, c("right", "bottom"))){
        new_x <- x - width
        new_y <- y
      } else if (identical(just, c("left", "center"))){
        new_x <- x - (0.25*width)
        new_y <- y - (0.5*height)
      } else if (identical(just, c("right", "center"))){
        new_x <- x - (0.75*width)
        new_y <- y - (0.5*height)
      } else if (identical(just, c("center", "bottom"))){
        new_x <- x - (0.5*width)
        new_y <- y
      } else if (identical(just, c("center", "top"))){
        new_x <- x - (0.5*width)
        new_y <- y - height
      } else {
        new_x <- x - (0.5*width)
        new_y <- y - (0.5*height)
      }

    } else if (length(just) == 1){

      if (just == "left"){
        new_x <- x - (0.25*width)
        new_y <- y - (0.5*height)
      } else if (just == "right"){
        new_x <- x - (0.75*width)
        new_y <- y - (0.5*height)
      } else if (just == "bottom"){
        new_x <- x - (0.5*width)
        new_y <- y
      } else if (just == "top"){
        new_x <- x - (0.5*width)
        new_y <- y - height
      } else {
        new_x <- x - (0.5*width)
        new_y <- y - (0.5*height)
      }

    }

    return(list(new_x, new_y))

  }

  ## Define a function that manually "clips" squares/triangles along edges
  manual_clip <- function(hic, hic_plot, resolution){

    ## lowest most genomic coordinate of data set
    left_min <- min(hic[,1])
    ## highest most genomic coordinate of data set
    top_max <- max(hic[,2]) + resolution

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

    if (left_min < hic_plot$chromstart & top_max > (hic_plot$chromend + resolution)){

      new_width <- resolution - (hic_plot$chromstart - left_min)
      new_height <- top_max - (hic_plot$chromend + resolution)

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

        new_width <- resolution - (hic_plot$chromstart - left_min)

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

      else if (top_max > (hic_plot$chromend + resolution)){

        new_height <- top_max - (hic_plot$chromend + resolution)

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

    return(list(all_squares, all_triangles))

  }

  ## Define a function that makes grobs for the hic diagonal
  hic_diagonal <- function(hic, resolution){

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

  hic_plot <- structure(list(chrom = chrom, chromstart = chromstart, chromend = chromend, x = x, y = y, width = width, height = NULL, justification = just,
                             zrange = zrange, color_palette = NULL, grobs = NULL), class = "bb_trianglehic")
  hic_plot$height <-unit(convertWidth(width, unitTo = get("page_units", envir = bbEnv), valueOnly = T)*0.5, get("page_units", envir = bbEnv))
  attr(x = hic_plot, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = hic_plot)
  errorcheck_bb_plotTriangleHic(hic = hic, hic_plot = hic_plot, resolution = resolution, norm = norm)

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  hic <- read_data(hic = hic, hic_plot = hic_plot, norm = norm, resolution = resolution)

  # ======================================================================================================================================================================================
  # SUBSET DATA
  # ======================================================================================================================================================================================

  hic <- subset_data(hic = hic, hic_plot = hic_plot, resolution = resolution)

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
    warning("Can't scale data to colors.")

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_hic", length(grep(pattern = "bb_trianglehic", x = currentViewports)) + 1)

  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(1/sqrt(two), "snpc"), width = unit(1/sqrt(2), "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.25, "npc"),
                   xscale = c(chromstart, chromend + resolution), yscale = c(chromstart, chromend + resolution),
                   just = "center",
                   name = vp_name,
                   angle = -45)

    if (draw == TRUE){

      grid.newpage()

    }

  } else {

    ## Get sides of viewport based on input width
    vp_side <- (convertWidth(width, unitTo = get("page_units", envir = bbEnv), valueOnly = T))/sqrt(two)


    ## Get bottom left point of triangle (hence bottom left of actual viewport) based on just
    bottom_coords <- convert_just(width = width, x = x, y = y, just = just)

    ## Make viewport
    vp <- viewport(height = unit(vp_side, get("page_units", envir = bbEnv)), width = unit(vp_side, get("page_units", envir = bbEnv)),
                   x = unit(bottom_coords[[1]], get("page_units", envir = bbEnv)),
                   y = unit(bottom_coords[[2]], get("page_units", envir = bbEnv)),
                   xscale = c(chromstart, chromend + resolution), yscale = c(chromstart, chromend + resolution),
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

  ## Manually "clip" the grobs that fall out of the desired region
  shapes <- manual_clip(hic = hic, hic_plot = hic_plot, resolution = resolution)
  squares <- shapes[[1]]
  triangles <- shapes[[2]]


  ## Make square grobs and add to grob gTree
  hic_squares <- rectGrob(x = squares$x,
                          y = squares$y,
                          just = c("left", "bottom"),
                          width = squares$width,
                          height = squares$height,
                          gp = gpar(col = NA, fill = squares$color),
                          default.units = "native")

  assign("hic_grobs2", addGrob(gTree = get("hic_grobs2", envir = bbEnv), child = hic_squares), envir = bbEnv)

  ## Make triangle grobs and add to grob gTree
  invisible(apply(triangles, 1, hic_diagonal, resolution = resolution))

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
