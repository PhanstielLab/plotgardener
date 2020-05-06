#' annotates loops in a Hi-C plot
#'
#' @param hic hic plot to annotate
#' @param loops bedpe file or dataframe in bedpe file format with loop positions (chr1 and chr2 must be numbers)
#' @param half which half of hic plots to annotate; default is "inherit", which will inherit whatever is plotted; other options are "both", "top", or "bottom"
#' @param shift number of pixels on either end of loop in box/circle; number of pixels for length of arrow
#' @param type type of annotation; options are "box", "circle", or "arrow"
#' @param lty line type; options are "solid" (1), "dashed" (2), "dotted" (3), "dotdash" (4), "longdash" (5), or "twodash" (6)
#' @param lwd line width;
#' @param col line color
#'
#' @export
bb_annotateLoops <- function(hic, loops, half = "inherit", shift = 4, type = "box", lty = "solid", lwd = 1, col = "black", rotation = 0){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function to catch errors for bb_annotateLoops
  errorcheck_bb_annotateLoops <- function(hic, loops, half, type){

    ###### hic #####

    ## check type of input for hic
    if (!class(hic) %in% c("bb_hic", "bb_trianglehic" )){

      stop("Input plot must be a plot of class \'bb_hic\' or \'bb_trianglehic\'.", call. = FALSE)

    }

    ###### loops #####

    ## if loops is a dataframe or datatable, it needs to be properly formatted
    if ("data.frame" %in% class(loops) && ncol(loops) < 6){

      stop("Invalid dataframe format. Dataframe must be in bedpe format.", call. = FALSE)

    }

    if ("data.frame" %in% class(loops) && nrow(loops) < 1){

      stop("Loop input contains no values.", call. = FALSE)

    }


    ## if it's a file path, it needs to exist
    if (!"data.frame" %in% class(loops)){

      # ## File extension
      # if (file_ext(loops) != "bedpe"){
      #
      #   stop("Invalid input. File must have a \".bedpe\" extension")
      #
      # }

      ## File existence
      if (!file.exists(loops)){

        stop(paste("File", loops, "does not exist."), call. = FALSE)

      }

    }

    ###### half #####

    ## half needs to be a valid option
    if (!half %in% c("inherit", "both", "top", "bottom")){

      stop("Invalid \'half\'.  Options are \'inherit\', \'both\', \'top\', or \'bottom\'.", call. = FALSE)

    }

    ## half needs to be able to align with what kind of hic plot is plotted
    if (class(hic) == "bb_hic"){

      if (is.null(hic$althalf)){

        if ((hic$half == "top" | hic$half == "bottom") && (half == "both")){

          stop("Invalid \'half\' of plot to annotate.", call. = FALSE)

        }

        if (hic$half == "top" & half == "bottom"){

          stop("Invalid \'half\' of plot to annotate.", call. = FALSE)

        }

        if (hic$half == "bottom" & half == "top"){

          stop("Invalid \'half\' of plot to annotate.", call. = FALSE)

        }

      } else {

        if (hic$althalf == "bottom"){

          message(paste("Attempting to annotate loops where", hic$chrom, "is on the x-axis and", hic$altchrom, "is on the y-axis."), call. = FALSE)

        } else if (hic$althalf == "top"){

          message(paste("Attempting to annotate loops where", hic$altchrom, "is on the x-axis and", hic$chrom, "is on the y-axis."), call. = FALSE)

        }

      }

    } else if (class(hic) == "bb_trianglehic"){

      if (half == "both" | half == "bottom"){

        warning("Plot of class \'bb_trianglehic\' detected.  Loops will automatically be annotated in the upper triangular of the plot.", call. = FALSE)

      }

    }

    ###### annotation #####

    ## Check type of annotation
    if (!type %in% c("box", "circle", "arrow")){

      stop("Invalid \'type\' of annotation.  Options are \'box\', \'circle\', or \'arrow\'.", call. = FALSE)

    }

  }

  ## Define a function that subsets loop data for hic region
  subset_loops <- function(loops, object){

    chrom <- paste0("chr", object$chrom)
    altchrom <- paste0("chr", object$altchrom)

    if (object$chrom == object$altchrom){

      loops_subset <- loops[which(loops[,1] == chrom & loops[,4] == chrom & loops[,2] >= object$chromstart & loops[,3] <= object$chromend
                                  & loops[,5] >= object$chromstart & loops[,6] <= object$chromend),]

    } else {

      if (object$chrom > object$altchrom){

        loops_subset <- loops[which(loops[,1] == altchrom & loops[,4] == chrom & loops[,2] >= object$altchromstart & loops[,3] <= object$altchromend
                                    & loops[,5] >= object$chromstart & loops[,6] <= object$chromend),]

      } else if (object$chrom < object$altchrom){

        loops_subset <- loops[which(loops[,1] == chrom & loops[,4] == altchrom & loops[,2] >= object$chromstart & loops[,3] <= object$chromend
                                    & loops[,5] >= object$altchromstart & loops[,6] <= object$altchromend),]

      }

    }

    return(loops_subset)

  }

  ## Define a function that parses an inherited half
  inherit_half <- function(hic){

    if (class(hic) == "bb_hic"){

      if (is.null(hic$althalf)){

        half <- hic$half

      } else {

        half <- hic$althalf

      }

    } else if (class(hic) == "bb_trianglehic"){

      half <- "top"

    }

    return(half)

  }

  ## Define a function that checks if any of the loop annotations fall out of the plotted region of a triangle hic plot
  annot_limits <- function(loop, hic, shift, type){

    if (class(hic) == "bb_trianglehic"){

      if (attributes(hic)$choppedPixels != 0){

        ## this is a function that will determine if two line segments intersect
        check_intersection <- function(point1, point2, point3, point4){
          ## point1 and point2 define one line segment
          ## point3 and point4 define other line segment
          onSegment<- function(pointp, pointq, pointr){

            if (pointq[[1]] <= max(pointp[[1]], pointr[[1]]) && pointq[[1]] >= min(pointp[[1]], pointr[[1]]) &&
                pointq[[2]] <= max(pointp[[2]], pointr[[2]]) && pointq[[2]] >= min(pointp[[2]], pointr[[2]])){
              return(TRUE)
            } else {

              return(FALSE)
            }

          }

          orientation <- function(pointp, pointq, pointr){

            value <- (pointq[[2]] - pointp[[2]])*(pointr[[1]] - pointq[[1]]) -
              (pointq[[1]] - pointp[[1]])*(pointr[[2]] - pointq[[2]])

            if (value == 0){
              return(0)
            } else {
              if(value > 0){
                return(1)
              } else {
                return(2)
              }
            }
          }

          o1 <- orientation(pointp = point1, pointq = point2, pointr = point3)
          o2 <- orientation(pointp = point1, pointq = point2, pointr = point4)
          o3 <- orientation(pointp = point3, pointq = point4, pointr = point1)
          o4 <- orientation(pointp = point3, pointq = point4, pointr = point2)

          if (o1 != o2 && o3 != o4){
            return(TRUE)
          }

          if (o1 == 0 && onSegment(pointp = point1, pointq = point3, pointr = point2)){
            return(TRUE)
          }

          if (o2 == 0 && onSegment(pointp = point1, pointq = point4, pointr = point2)){
            return(TRUE)
          }

          if (o3 ==0 && onSegment(pointp = point3, pointq = point1, pointr = point4)){
            return(TRUE)
          }

          if(o4 == 0 && onSegment(pointp = point3, pointq = point2, pointr = point4)){
            return(TRUE)
          }

          return(FALSE)

        }

        ## this is a function that will determine if a line intersects a circle
        checkCircle <- function(a, b, c, x, y, radius){

          dist <- (abs(a*x + b*y + c))/sqrt(a^2+b^2)

          if (radius == dist){
            return(FALSE)
          }  else if (radius > dist){
            return(TRUE)
          } else {
            return(FALSE)
          }

        }

        normx <- ((hic$chromstart + attributes(hic)$choppedPixels*hic$resolution) - hic$chromstart)/(hic$chromend - hic$chromstart)
        normy <- ((hic$chromend - attributes(hic)$choppedPixels*hic$resolution) - hic$chromstart)/(hic$chromend - hic$chromstart)

        slope <- (1 - normy)/(normx - 0)

        if (type == "box"){

          side <- (as.numeric(loop[6]) - as.numeric(loop[5])) + (2 * shift * hic$resolution)
          center_x <- 0.5 * (as.numeric(loop[2]) + as.numeric(loop[3]))
          center_y <- 0.5 * (as.numeric(loop[5]) + as.numeric(loop[6]))

          check1 <- check_intersection(point1 = list(hic$chromstart, hic$chromend - attributes(hic)$choppedPixels*hic$resolution),
                                   point2 = list(hic$chromstart + attributes(hic)$choppedPixels*hic$resolution, hic$chromend),
                                   point3 = list(center_x - 0.5*side, center_y - 0.5*side),
                                   point4 = list(center_x - 0.5*side, center_y + 0.5*side))

          check2 <- check_intersection(point1 = list(hic$chromstart, hic$chromend - attributes(hic)$choppedPixels*hic$resolution),
                                   point2 = list(hic$chromstart + attributes(hic)$choppedPixels*hic$resolution, hic$chromend),
                                   point3 = list(center_x - 0.5*side, center_y - 0.5*side),
                                   point4 = list(center_x + 0.5*side, center_y - 0.5*side))

          check3 <- check_intersection(point1 = list(hic$chromstart, hic$chromend - attributes(hic)$choppedPixels*hic$resolution),
                                       point2 = list(hic$chromstart + attributes(hic)$choppedPixels*hic$resolution, hic$chromend),
                                       point3 = list(center_x - 0.5*side, center_y + 0.5*side),
                                       point4 = list(center_x + 0.5*side, center_y + 0.5*side))

          check4 <- check_intersection(point1 = list(hic$chromstart, hic$chromend - attributes(hic)$choppedPixels*hic$resolution),
                                       point2 = list(hic$chromstart + attributes(hic)$choppedPixels*hic$resolution, hic$chromend),
                                       point3 = list(center_x + 0.5*side, center_y + 0.5*side),
                                       point4 = list(center_x + 0.5*side, center_y - 0.5*side))



          if (check1 == TRUE | check2 == TRUE | check3 == TRUE | check4 == TRUE){

            warning("Loop annotation falls out of plotted region.", call. = FALSE)

          }

          norm_x0 <- ((center_x + 0.5*side) - hic$chromstart)/(hic$chromend - hic$chromstart)
          solve_val <- slope*norm_x0 + normy
          norm_y0 <- ((center_y - 0.5*side) - hic$chromstart)/(hic$chromend - hic$chromstart)
          if (norm_y0 > solve_val){
            warning("Loop annotation falls out of plotted region.", call. = FALSE)
          }

        } else if (type == "circle"){

          radius <- (0.5 * (as.numeric(loop[6]) - as.numeric(loop[5]))) + (shift * hic$resolution)
          center_x <- 0.5 * (as.numeric(loop[2]) + as.numeric(loop[3]))
          center_y <- 0.5 * (as.numeric(loop[5]) + as.numeric(loop[6]))

          rad1 <- center_x - radius
          norm_rad1 <- (rad1 - hic$chromstart)/(hic$chromend - hic$chromstart)
          rad2 <- center_x + radius
          norm_rad2 <- (rad2 - hic$chromstart)/(hic$chromend - hic$chromstart)
          norm_radius <- norm_rad2 - norm_rad1

          norm_centerx <- (center_x - hic$chromstart)/(hic$chromend - hic$chromstart)
          norm_centery <- (center_y - hic$chromstart)/(hic$chromend - hic$chromstart)

          check1 <- checkCircle(a = -slope, b = 1, c = -normy, x = norm_centerx, y = norm_centery, radius = norm_radius)

          if (check1 == TRUE){

            warning("Loop annotation falls out of plotted region.", call. = FALSE)

          }

        } else if (type == "arrow"){

          x0 <- as.numeric(loop[2]) - (0.5 * (as.numeric(loop[6]) - as.numeric(loop[5])))
          y0 <- as.numeric(loop[6]) + (0.5 * (as.numeric(loop[6]) - as.numeric(loop[5])))
          x1 <-  x0 - (shift * hic$resolution)
          y1 <- y0 + (shift * hic$resolution)

          check1 <- check_intersection(point1 = list(hic$chromstart, hic$chromend - attributes(hic)$choppedPixels*hic$resolution),
                                       point2 = list(hic$chromstart + attributes(hic)$choppedPixels*hic$resolution, hic$chromend),
                                       point3 = list(x0, y0), point4 = list(x1, y1))
          if (check1 == TRUE){

            warning("Loop annotation falls out of plotted region.", call. = FALSE)
          }
          norm_x0 <- (x0 - hic$chromstart)/(hic$chromend - hic$chromstart)
          solve_val <- slope*norm_x0 + normy
          norm_y0 <- (y0 - hic$chromstart)/(hic$chromend - hic$chromstart)
          if (norm_y0 > solve_val){
            warning("Loop annotation falls out of plotted region.", call. = FALSE)
          }

        }

      }

    }

  }

  ## Define a function to add box annotation
  boxAnnotation <- function(df, hic, object, shift, half){

    side <- (as.numeric(df[6]) - as.numeric(df[5])) + (2 * shift * hic$resolution)

    if (half == "bottom"){

      center_x <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      rect1 <- rectGrob(x = center_x, y = center_y, width = side, height = side, default.units = "native",
                        gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(rect1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = rect1), envir = bbEnv)

    } else if (half == "top"){

      center_x <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      rect1 <- rectGrob(x = center_x, y = center_y, width = side, height = side, default.units = "native",
                        gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(rect1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = rect1), envir = bbEnv)

    } else if (half == "both"){

      ## BOTTOM
      center_x1 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y1 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))

      ## TOP
      center_x2 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y2 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))

      rect1 <- rectGrob(x = center_x1, y = center_y1, width = side, height = side, default.units = "native",
                        gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))
      rect2 <- rectGrob(x = center_x2, y = center_y2, width = side, height = side, default.units = "native",
                        gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(rect1)
      grid.draw(rect2)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = rect1), envir = bbEnv)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = rect2), envir = bbEnv)

    }

  }

  ## Define a function to add circle annotation
  circleAnnotation <- function(df, hic, object, shift, half){

    radius <- (0.5 * (as.numeric(df[6]) - as.numeric(df[5]))) + (shift * hic$resolution)

    if (half == "bottom"){

      center_x <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      circ1 <- circleGrob(x = center_x, y = center_y, r = radius, default.units = "native",
                          gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(circ1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = circ1), envir = bbEnv)

    } else if (half == "top"){

      center_x <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      circ1 <- circleGrob(x = center_x, y = center_y, r = radius, default.units = "native",
                          gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(circ1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = circ1), envir = bbEnv)

    } else if (half == "both"){

      ## BOTTOM
      center_x1 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y1 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))

      ## TOP
      center_x2 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y2 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))

      circ1 <- circleGrob(x = center_x1, y = center_y1, r = radius, default.units = "native",
                          gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))
      circ2 <- circleGrob(x = center_x2, y = center_y2, r = radius, default.units = "native",
                          gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(circ1)
      grid.draw(circ2)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = circ1), envir = bbEnv)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = circ2), envir = bbEnv)

    }

  }

  ## Define a function to add arrow annotation
  arrowAnnotation <- function(df, hic, object, shift, half){

    if (half == "bottom"){

     x0 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
     y0 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

     arrow1 <- segmentsGrob(x0 = x0, y0 = y0, x1 = x0 + (shift * hic$resolution), y1 = y0 - (shift * hic$resolution),
                   arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                   gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, fill = object$gpar$col))

     grid.draw(arrow1)
     assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow1), envir = bbEnv)

    } else if (half == "top"){

      x0 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
      y0 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

      arrow1 <- segmentsGrob(x0 = x0, y0 = y0, x1 = x0 - (shift * hic$resolution), y1 = y0 + (shift * hic$resolution),
                    arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                    gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, fill = object$gpar$col))

      grid.draw(arrow1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow1), envir = bbEnv)

    } else if (half == "both"){

      ## BOTTOM
      x01 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
      y01 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

      ## TOP
      x02 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
      y02 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

      arrow1 <- segmentsGrob(x0 = x01, y0 = y01, x1 = x01 + (shift * hic$resolution), y1 = y01 - (shift * hic$resolution),
                    arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                    gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, fill = object$gpar$col))
      arrow2 <- segmentsGrob(x0 = x02, y0 = y02, x1 = x02 - (shift * hic$resolution), y1 = y02 + (shift * hic$resolution),
                    arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                    gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, fill = object$gpar$col))

      grid.draw(arrow1)
      grid.draw(arrow2)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow1), envir = bbEnv)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow2), envir = bbEnv)

    }

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT: GET REGION/DIMENSIONS FROM HIC PLOT INPUT
  # ======================================================================================================================================================================================

  loop_annot <- structure(list(chrom = hic$chrom, chromstart = hic$chromstart, chromend = hic$chromend, altchrom = hic$altchrom,
                               altchromstart = hic$altchromstart, altchromend = hic$altchromend, x = hic$x, y = hic$y, width = hic$width,
                               height = hic$height, justification = hic$just, grobs = NULL,
                             gpar = list(lty = lty, lwd = lwd, col = col)), class = "bb_loop")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot annotate loops without a BentoBox page.")
  errorcheck_bb_annotateLoops(hic = hic, loops = loops, half = half, type = type)

  # ======================================================================================================================================================================================
  # PARSE INHERITED HALF
  # ======================================================================================================================================================================================

  if (half == "inherit"){

    half <- inherit_half(hic = hic)

  }

  if (class(hic) == "bb_trianglehic"){

    half <- "top"

  }

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  if (!"data.frame" %in% class(loops)){

    loops <- as.data.frame(data.table::fread(loops))
    if (nrow(loops) < 1){
      stop("Loop input contains no values.", call. = FALSE)
    }

  }

  # ======================================================================================================================================================================================
  # SUBSET FOR LOOPS IN REGION
  # ======================================================================================================================================================================================
  ## Assuming loops are in first six columns only
  loops <- loops[,1:6]

  loops_subset <- subset_loops(loops = loops, object = loop_annot)

  # ======================================================================================================================================================================================
  # CHECK THE ANNOTATION LIMITS OF THE LOOPS
  # ======================================================================================================================================================================================

  invisible(apply(loops_subset, 1, annot_limits, hic = hic, shift = shift, type = type))

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Name viewport
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_loopAnnotation", length(grep(pattern = "bb_loopAnnotation", x = currentViewports)) + 1)

  ## Make viewport based on hic input viewport
  if (class(hic) == "bb_hic"){

    vp <- viewport(height = hic$grobs$vp$height, width = hic$grobs$vp$width,
                   x = hic$grobs$vp$x, y = hic$grobs$vp$y,
                   clip = "on",
                   xscale = hic$grobs$vp$xscale,
                   yscale = hic$grobs$vp$yscale,
                   just = hic$grobs$vp$justification,
                   name = vp_name)
  } else if (class(hic) == "bb_trianglehic"){

    vp <- viewport(height = hic$grobs$vp$height, width = hic$grobs$vp$width,
                   x = hic$grobs$vp$x, y = hic$grobs$vp$y,
                   xscale = hic$grobs$vp$xscale,
                   yscale = hic$grobs$vp$yscale,
                   just = hic$grobs$vp$justification,
                   name = vp_name,
                   angle = -45)

  }

  pushViewport(vp)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE OF GROBS
  # ======================================================================================================================================================================================

  assign("annotation_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  if (type == "box"){

    invisible(apply(loops_subset, 1, boxAnnotation, hic = hic, object = loop_annot, shift = shift, half = half))

  } else if (type == "circle"){

    invisible(apply(loops_subset, 1, circleAnnotation, hic = hic, object = loop_annot, shift = shift, half = half))

  } else if (type == "arrow"){

    invisible(apply(loops_subset, 1, arrowAnnotation, hic = hic, object = loop_annot, shift = shift, half = half))

  }

  ## Go back to root viewport
  upViewport()

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  loop_annot$grobs <- get("annotation_grobs", envir = bbEnv)
  grid.draw(loop_annot$grobs)
  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(loop_annot)

}
