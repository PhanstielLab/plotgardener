#' annotates loops in a Hi-C plot
#'
#' @param hic hic plot to annotate
#' @param loops bedpe file or dataframe in bedpe file format with loop positions (chr1 and chr2 must be numbers)
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param half which half of hic plots to annotate; default is "inherit", which will inherit whatever is plotted; other options are "both", "top", or "bottom"
#' @param shift number of pixels on either end of loop in box/circle; number of pixels for length of arrow
#' @param type type of annotation; options are "box", "circle", or "arrow"
#'
#' @export
bb_annotateLoops <- function(hic, loops, params = NULL, half = "inherit", shift = 4, type = "box", ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================
  ## For more accurate calculation of sqrt(2)
  two <- mpfr(2, 120)

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

    numberChrom <- as.numeric(gsub("chr", "", object$chrom))
    numberAltChrom <- as.numeric(gsub("chr", "", object$altchrom))

    if (numberChrom == numberAltChrom){

      loops_subset <- loops[which(loops[,1] == object$chrom & loops[,4] == object$chrom & loops[,2] >= object$chromstart & loops[,3] <= object$chromend
                                  & loops[,5] >= object$chromstart & loops[,6] <= object$chromend),]

    } else {

      if (numberChrom > numberAltChrom){

        loops_subset <- loops[which(loops[,1] == object$altchrom & loops[,4] == object$chrom & loops[,2] >= object$altchromstart & loops[,3] <= object$altchromend
                                    & loops[,5] >= object$chromstart & loops[,6] <= object$chromend),]

      } else if (numberChrom < numberAltChrom){

        loops_subset <- loops[which(loops[,1] == object$chrom & loops[,4] == object$altchrom & loops[,2] >= object$chromstart & loops[,3] <= object$chromend
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

  ## Define a function to add box annotation
  boxAnnotation <- function(df, hic, object, shift, half){

    side <- (as.numeric(df[6]) - as.numeric(df[5])) + (2 * shift * hic$resolution)

    if (half == "bottom"){

      center_x <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      rect1 <- rectGrob(x = center_x, y = center_y, width = side, height = side, default.units = "native",
                        gp = object$gp)

      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = rect1), envir = bbEnv)

    } else if (half == "top"){

      center_x <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      rect1 <- rectGrob(x = center_x, y = center_y, width = side, height = side, default.units = "native",
                        gp = object$gp)

      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = rect1), envir = bbEnv)

    } else if (half == "both"){

      ## BOTTOM
      center_x1 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y1 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))

      ## TOP
      center_x2 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y2 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))

      rect1 <- rectGrob(x = center_x1, y = center_y1, width = side, height = side, default.units = "native",
                        gp = object$gp)
      rect2 <- rectGrob(x = center_x2, y = center_y2, width = side, height = side, default.units = "native",
                        gp = object$gp)

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
                          gp = object$gp)

      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = circ1), envir = bbEnv)

    } else if (half == "top"){

      center_x <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      circ1 <- circleGrob(x = center_x, y = center_y, r = radius, default.units = "native",
                          gp = object$gp)

      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = circ1), envir = bbEnv)

    } else if (half == "both"){

      ## BOTTOM
      center_x1 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))
      center_y1 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))

      ## TOP
      center_x2 <- 0.5 * (as.numeric(df[2]) + as.numeric(df[3]))
      center_y2 <- 0.5 * (as.numeric(df[5]) + as.numeric(df[6]))

      circ1 <- circleGrob(x = center_x1, y = center_y1, r = radius, default.units = "native",
                          gp = object$gp)
      circ2 <- circleGrob(x = center_x2, y = center_y2, r = radius, default.units = "native",
                          gp = object$gp)

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
                   gp = object$gp)

     assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow1), envir = bbEnv)

    } else if (half == "top"){

      x0 <- as.numeric(df[2]) - (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))
      y0 <- as.numeric(df[6]) + (0.5 * (as.numeric(df[6]) - as.numeric(df[5])))

      arrow1 <- segmentsGrob(x0 = x0, y0 = y0, x1 = x0 - (shift * hic$resolution), y1 = y0 + (shift * hic$resolution),
                    arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                    gp = object$gp)

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
                    gp = object$gp)
      arrow2 <- segmentsGrob(x0 = x02, y0 = y02, x1 = x02 - (shift * hic$resolution), y1 = y02 + (shift * hic$resolution),
                    arrow = arrow(length = unit(0.1, "inches"), ends = "first", type = "closed"), default.units = "native",
                    gp = object$gp)

      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow1), envir = bbEnv)
      assign("annotation_grobs", addGrob(gTree = get("annotation_grobs", envir = bbEnv), child = arrow2), envir = bbEnv)

    }

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(half)) half <- NULL
  if(missing(shift)) shift <- NULL
  if(missing(type)) type <- NULL

  ## Check if hic/loops arguments are missing (could be in object)
  if(!hasArg(hic)) hic <- NULL
  if(!hasArg(loops)) loops <- NULL

  ## Compile all parameters into an internal object
  bb_loopInternal <- structure(list(hic = hic, loops = loops, half = half, shift = shift, type = type), class = "bb_loopInternal")

  bb_loopInternal <- parseParams(bb_params = params, object_params = bb_loopInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_loopInternal$half)) bb_loopInternal$half <- "inherit"
  if(is.null(bb_loopInternal$shift)) bb_loopInternal$shift <- 4
  if(is.null(bb_loopInternal$type)) bb_loopInternal$type <- "box"

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT: GET REGION/DIMENSIONS FROM HIC PLOT INPUT
  # ======================================================================================================================================================================================

  loop_annot <- structure(list(chrom = bb_loopInternal$hic$chrom, chromstart = bb_loopInternal$hic$chromstart, chromend = bb_loopInternal$hic$chromend, altchrom = bb_loopInternal$hic$altchrom,
                               altchromstart = bb_loopInternal$hic$altchromstart, altchromend = bb_loopInternal$hic$altchromend, x = bb_loopInternal$hic$x, y = bb_loopInternal$hic$y,
                               width = bb_loopInternal$hic$width, height = bb_loopInternal$hic$height, just = bb_loopInternal$hic$just, grobs = NULL,
                               gp = gpar(...)), class = "bb_loop")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot annotate loops without a BentoBox page.")
  if(is.null(bb_loopInternal$hic)) stop("argument \"hic\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_loopInternal$loops)) stop("argument \"loops\" is missing, with no default.", call. = FALSE)


  errorcheck_bb_annotateLoops(hic = bb_loopInternal$hic, loops = bb_loopInternal$loops,
                              half = bb_loopInternal$half, type = bb_loopInternal$type)

  # ======================================================================================================================================================================================
  # PARSE INHERITED HALF
  # ======================================================================================================================================================================================

  half <- bb_loopInternal$half
  if (half == "inherit"){

    half <- inherit_half(hic = bb_loopInternal$hic)

  }

  if (class(bb_loopInternal$hic) == "bb_trianglehic"){

    half <- "top"

  }

  # ======================================================================================================================================================================================
  # READ IN FILE OR DATAFRAME
  # ======================================================================================================================================================================================

  loops <- bb_loopInternal$loops
  if (!"data.frame" %in% class(loops)){

    loops <- as.data.frame(data.table::fread(loops))
    if (nrow(loops) < 1){
      warning("Loop input contains no values.", call. = FALSE)
    }

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
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_loopAnnotation", length(grep(pattern = "bb_loopAnnotation", x = currentViewports)) + 1)

  ## Make viewport based on hic input viewport
  if (class(bb_loopInternal$hic) == "bb_hic"){

    vp <- viewport(height = bb_loopInternal$hic$grobs$vp$height, width = bb_loopInternal$hic$grobs$vp$width,
                   x = bb_loopInternal$hic$grobs$vp$x, y = bb_loopInternal$hic$grobs$vp$y,
                   clip = "on",
                   xscale = bb_loopInternal$hic$grobs$vp$xscale,
                   yscale = bb_loopInternal$hic$grobs$vp$yscale,
                   just = bb_loopInternal$hic$grobs$vp$justification,
                   name = vp_name)
  } else if (class(bb_loopInternal$hic) == "bb_trianglehic"){


    width <- convertUnit(bb_loopInternal$hic$outsideVP$width, unitTo = get("page_units", bbEnv), valueOnly = T)

    vp <- viewport(height = unit(width/sqrt(two), get("page_units", bbEnv)), width = unit(width/sqrt(two), get("page_units", bbEnv)),
                   x = bb_loopInternal$hic$outsideVP$x, y = bb_loopInternal$hic$outsideVP$y,
                   xscale = bb_loopInternal$hic$grobs$vp$xscale,
                   yscale = bb_loopInternal$hic$grobs$vp$yscale,
                   just = bb_loopInternal$hic$outsideVP$justification,
                   name = vp_name,
                   angle = -45)

  }


  # ======================================================================================================================================================================================
  # INITIALIZE GTREE OF GROBS
  # ======================================================================================================================================================================================

  assign("annotation_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  if (nrow(loops_subset) > 0){

    if (bb_loopInternal$type == "box"){

      loop_annot$gp$fill <- NA
      invisible(apply(loops_subset, 1, boxAnnotation, hic = bb_loopInternal$hic, object = loop_annot, shift = bb_loopInternal$shift, half = half))

    } else if (bb_loopInternal$type == "circle"){
      loop_annot$gp$fill <- NA
      invisible(apply(loops_subset, 1, circleAnnotation, hic = bb_loopInternal$hic, object = loop_annot, shift = bb_loopInternal$shift, half = half))

    } else if (bb_loopInternal$type == "arrow"){
      if (is.null(loop_annot$gp$col) & is.null(loop_annot$gp$fill)){
        loop_annot$gp$fill <- "black"
      } else {
        if(is.null(loop_annot$gp$fill)){
          loop_annot$gp$fill <- loop_annot$gp$col
        }
      }

      invisible(apply(loops_subset, 1, arrowAnnotation, hic = bb_loopInternal$hic, object = loop_annot, shift = bb_loopInternal$shift, half = half))

    }

  } else {

    warning("No loops found in region.", call. = FALSE)
  }


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
