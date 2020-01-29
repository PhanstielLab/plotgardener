#' annotates loops in a Hi-C plot
#'
#' @param hic hic plot to annotate
#' @param loops bedpe file or dataframe in bedpe file format with loop positions
#' @param shift number of pixels on either end of loop in box/circle; number of pixels for length of arrow
#' @param type type of annotation; options are "box", "circle", or "arrow"
#' @param lty line type; options are "solid" (1), "dashed" (2), "dotted" (3), "dotdash" (4), "longdash" (5), or "twodash" (6)
#' @param lwd line width;
#' @param col line color
#'
#' @export
bb_annotateLoops <- function(hic, loops, shift = 4, type = "box", lty = "dashed", lwd = 1, col = "black"){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function to catch errors for bb_annotateLoops
  errorcheck_bb_annotateLoops <- function(hic, loops, object){

    ## Check that plot is actually plotted
    if (length(grep(pattern = hic$grobs$vp$name, x = grid.ls(grobs = FALSE, viewport = TRUE, print = FALSE)$name)) == 0){

      stop("Hi-C plot is not plotted.")

    }

    ## Check that input plot is a hic plot
    if (class(hic) != "bb_hic"){

      stop("Input plot is not a Hi-C plot.")

    }

    ## Check loops dataframe format
    if (class(loops) %in% "data.frame" && ncol(loops) < 6){

      stop("Invalid dataframe format. Dataframe must be in bedpe format.")

    }

    ## Check loops bedpe file
    if (!class(loops) %in% "data.frame"){

      ## File extension
      if (file_ext(loops) != "bedpe"){

        stop("Invalid input. File must have a \".bedpe\" extension")

      }

      ## File existence
      if (!file.exists(loops)){

        stop(paste("File", loops, "does not exist."))

      }

    }

    ## Check type of annotation
    if (!object$type %in% c("box", "circle", "arrow")){

      stop("Invalid \'type\'.  Options are \'box\', \'circle\', or \'arrow\'.")

    }

  }

  ## Define a function that subsets loop data for hic region
  subset_loops <- function(loops, object){

    if (object$chrom == object$altchrom){

      loops_subset <- loops[which(loops[,1] == object$chrom & loops[,4] == object$chrom & loops[,2] >= object$chromstart & loops[,3] <= object$chromend
                                  & loops[,5] >= object$chromstart & loops[,6] <= object$chromend),]

    } else {

      if (object$chrom > object$altchrom){

        loops_subset <- loops[which(loops[,1] == object$altchrom & loops[,4] == object$chrom & loops[,2] >= object$altchromstart & loops[,3] <= object$altchromend
                                    & loops[,5] >= object$chromstart & loops[,6] <= object$chromend),]

      } else if (object$chrom < object$altchrom){

        loops_subset <- loops[which(loops[,1] == object$chrom & loops[,4] == object$altchrom & loops[,2] >= object$chromstart & loops[,3] <= object$chromend
                                    & loops[,5] >= object$altchromstart & loops[,6] <= object$altchromend),]

      }

    }

    return(loops_subset)

  }

  ## Define a function to add box annotation
  boxAnnotation <- function(df, hic, object, shift){

    side <- hic$additional_parameters$resolution + (2 * shift * hic$additional_parameters$resolution)

    if (hic$additional_parameters$half == "bottom"){

      center_x <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      rect1 <- rectGrob(x = center_x, y = center_y, width = side, height = side, default.units = "native",
                        gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(rect1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = rect1), envir = bbEnv)

    } else if (hic$additional_parameters$half == "top"){

      center_x <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      rect1 <- rectGrob(x = center_x, y = center_y, width = side, height = side, default.units = "native",
                        gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(rect1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = rect1), envir = bbEnv)

    } else if (hic$additional_parameters$half == "both"){

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
  circleAnnotation <- function(df, hic, object, shift){

    radius <- (0.5 * hic$additional_parameters$resolution) + (shift * hic$additional_parameters$resolution)

    if (hic$additional_parameters$half == "bottom"){

      center_x <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      circ1 <- circleGrob(x = center_x, y = center_y, r = radius, default.units = "native",
                          gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(circ1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = circ1), envir = bbEnv)

    } else if (hic$additional_parameters$half == "top"){

      center_x <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      circ1 <- circleGrob(x = center_x, y = center_y, r = radius, default.units = "native",
                          gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, col = object$gpar$col, fill = NA))

      grid.draw(circ1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = circ1), envir = bbEnv)

    } else if (hic$additional_parameters$half == "both"){

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
  arrowAnnotation <- function(df, hic, object, shift){

    if (hic$additional_parameters$half == "bottom"){

     x0 <- as.numeric(df[6]) + (0.5 * hic$additional_parameters$resolution)
     y0 <- as.numeric(df[2]) - (0.5 * hic$additional_parameters$resolution)

     arrow1 <- segmentsGrob(x0 = x0, y0 = y0, x1 = x0 + (shift * hic$additional_parameters$resolution), y1 = y0 - (shift * hic$additional_parameters$resolution),
                   arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                   gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, fill = object$gpar$col))

     grid.draw(arrow1)
     assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow1), envir = bbEnv)

    } else if (hic$additional_parameters$half == "top"){

      x0 <- as.numeric(df[2]) - (0.5 * hic$additional_parameters$resolution)
      y0 <- as.numeric(df[6]) + (0.5 * hic$additional_parameters$resolution)

      arrow1 <- segmentsGrob(x0 = x0, y0 = y0, x1 = x0 - (shift * hic$additional_parameters$resolution), y1 = y0 + (shift * hic$additional_parameters$resolution),
                    arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                    gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, fill = object$gpar$col))

      grid.draw(arrow1)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow1), envir = bbEnv)

    } else if (hic$additional_parameters$half == "both"){

      ## BOTTOM
      x01 <- as.numeric(df[6]) + (0.5 * hic$additional_parameters$resolution)
      y01 <- as.numeric(df[2]) - (0.5 * hic$additional_parameters$resolution)

      ## TOP
      x02 <- as.numeric(df[2]) - (0.5 * hic$additional_parameters$resolution)
      y02 <- as.numeric(df[6]) + (0.5 * hic$additional_parameters$resolution)

      arrow1 <- segmentsGrob(x0 = x01, y0 = y01, x1 = x01 + (shift * hic$additional_parameters$resolution), y1 = y01 - (shift * hic$additional_parameters$resolution),
                    arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                    gp = gpar(lty = object$gpar$lty, lwd = object$gpar$lwd, fill = object$gpar$col))
      arrow2 <- segmentsGrob(x0 = x02, y0 = y02, x1 = x02 - (shift * hic$additional_parameters$resolution), y1 = y02 + (shift * hic$additional_parameters$resolution),
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

  loop_annot <- structure(list(type = type, chrom = hic$chrom, chromstart = hic$chromstart, chromend = hic$chromend, altchrom = hic$altchrom,
                               altchromstart = hic$altchromstart, altchromend = hic$altchromend, x = hic$x, y = hic$y, width = hic$width,
                               height = hic$height, justification = hic$just, grobs = NULL,
                             gpar = list(lty = lty, lwd = lwd, col = col)), class = "bb_loopAnnotation")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot annotate loops without a BentoBox page.")
  errorcheck_bb_annotateLoops(hic = hic, loops = loops, object = loop_annot)

  # ======================================================================================================================================================================================
  # READ IN BEDPE FILE
  # ======================================================================================================================================================================================

  if (!class(loops) %in% "data.frame"){

    loops <- data.table::fread(loops)

  }

  # ======================================================================================================================================================================================
  # SUBSET FOR LOOPS IN REGION
  # ======================================================================================================================================================================================
  ## Assuming loops are in first six columns only
  loops <- loops[,1:6]

  loops_subset <- subset_loops(loops = loops, object = loop_annot)

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Name viewport
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_loopAnnotation", length(grep(pattern = "bb_loopAnnotation", x = current_viewports)) + 1)

  ## Make viewport based on hic input viewport
  vp <- viewport(height = hic$grobs$vp$height, width = hic$grobs$vp$width,
                 x = hic$grobs$vp$x, y = hic$grobs$vp$y, clip = "on", xscale = hic$grobs$vp$xscale, yscale = hic$grobs$vp$yscale, just = hic$grobs$vp$justification,
                 name = vp_name)
  pushViewport(vp)

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE OF GROBS
  # ======================================================================================================================================================================================

  assign("annotation_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  if (loop_annot$type == "box"){

    invisible(apply(loops_subset, 1, boxAnnotation, hic = hic, object = loop_annot, shift = shift))

  } else if (loop_annot$type == "circle"){

    invisible(apply(loops_subset, 1, circleAnnotation, hic = hic, object = loop_annot, shift = shift))

  } else if (loop_annot$type == "arrow"){

    invisible(apply(loops_subset, 1, arrowAnnotation, hic = hic, object = loop_annot, shift = shift))

  }

  ## Go back to root viewport
  upViewport()

  #seekViewport(name = vp$name)
  #popViewport()
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
