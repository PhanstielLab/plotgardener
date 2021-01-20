#' Place a BentoBox plot that has been previously created but not drawn
#'
#' @usage
#' bb_pagePlotPlace(plot, x, y, width, height, just = c("left", "top"), default.units = "inches")
#'
#' @param plot BentoBox plot object to be placed, defined by the output of a BentoBox plotting function.
#' @param x A numeric or unit object specifying plot x-location.
#' @param y A numeric or unit object specifying plot y-location.
#' @param width A numeric or unit object specifying plot width.
#' @param height A numeric or unit object specifying plot height.
#' @param just Justification of plot relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{x}, \code{y}, \code{width}, or \code{height} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#'
#' @return Function will update dimensions of an input plot and return an updated BentoBox plot object.
#'
#' @examples
#' ## Load Hi-C data
#' data("bb_hicData")
#'
#' ## Create, but do not plot, square Hi-C plot
#' hicPlot <- bb_plotHicSquare(hicData = bb_hicData, resolution = 10000, zrange = c(0, 70), chrom = "chr21", chromstart = 28000000, chromend = 30300000, draw = FALSE)
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 4, height = 3.5, default.units = "inches", xgrid = 0, ygrid = 0)
#'
#' ## Place Hi-C plot on BentoBox page
#' bb_pagePlotPlace(plot = hicPlot, x = 0.5, y = 0.5, width = 2.5, height = 2.5, just = c("left", "top"), default.units = "inches", draw = TRUE)
#'
#' @export
bb_pagePlotPlace <- function(plot, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches",
                             draw = TRUE, params = NULL){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  set_x <- function(object, val){

    if (is.null(val)){

      object['x'] <- list(NULL)

    } else {

      object$x <- val

    }

    return(object)

  }

  set_y <- function(object, val){

    if (is.null(val)){

      object['y'] <- list(NULL)

    } else {

      object$y <- val

    }

    return(object)

  }

  set_width <- function(object, val){

    if (is.null(val)){

      object['width'] <- list(NULL)

    } else {

      object$width <- val

    }

    return(object)

  }

  set_height <- function(object, val){

    if (is.null(val)){

      object['height'] <- list(NULL)

    } else {

      object$height <- val

    }

    return(object)

  }

  set_values <- function(object, x, y, width, height){

    object <- set_x(object = object, val = x)
    object <- set_y(object = object, val = y)
    object <- set_width(object = object, val = width)
    object <- set_height(object = object, val = height)

    return(object)
  }

  replace_value <- function(val, new){

    return(new[val])

  }

  ## Define a function to parse coordinates
  parse_coordinates <- function(input_plot, output_plot){

    ## Make sublists of the dimensions and coordinates of the input and output plots
    inputCoords <- list(x = input_plot$x, y = input_plot$y, width = input_plot$width, height = input_plot$height, justification = input_plot$jusification)
    outputCoords <- list(x = output_plot$x, y = output_plot$y, width = output_plot$width, height = output_plot$height, justification = output_plot$justification)

    ## Determine which values in the output plot are NULL
    to_replace <- names(outputCoords[sapply(outputCoords, is.null)])
    not_replace <- outputCoords[!sapply(outputCoords, is.null)]

    ## Get corresponding values for those that are NULL from the input plot
    replaced <- unlist(lapply(to_replace, replace_value, new = inputCoords), recursive = F)

    ## Recombine values that weren't replaced and those that were
    new_coords <- c(not_replace, replaced)

    ## Assign new values to object
    # output_plot$x <- new_coords$x
    # output_plot$y <- new_coords$y
    # output_plot$width <- new_coords$width
    # output_plot$height <- new_coords$height

    output_plot <- set_values(object = output_plot, x = new_coords$x, y = new_coords$y, width = new_coords$width, height = new_coords$height)
    output_plot$justification <- new_coords$justification

    ## Return object
    return(output_plot)

  }

  ## Define a function that renames grobs and copies to a new gtree
  copy_grobs <- function(grob){

    new_name <- grobName(grob)

    grob$name <- new_name
    assign("new_gtree", addGrob(gTree = get("new_gtree", envir = bbEnv), child = grob), envir = bbEnv)

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if plot argument is missing (could be in object)
  if(!hasArg(plot)) plot <- NULL

  ## Compile all parameters into an internal object
  bb_place <- structure(list(plot = plot, x = x, y = y, width = width, height = height, draw = draw,
                                     just = just, default.units = default.units), class = "bb_place")

  bb_place <- parseParams(bb_params = params, object_params = bb_place)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_place$just)) bb_place$just <- c("left", "top")
  if(is.null(bb_place$default.units)) bb_place$default.units <- "inches"
  if(is.null(bb_place$draw)) bb_place$draw <- TRUE

  # ======================================================================================================================================================================================
  # ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_place$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # INITIALIZE PLOT OBJECT COPY
  # ======================================================================================================================================================================================

  object <- bb_place$plot

  # ======================================================================================================================================================================================
  # PARSE UNITS FOR INPUTS
  # ======================================================================================================================================================================================

  if (!is.null(bb_place$x)){

    if (!"unit" %in% class(bb_place$x)){

      if (!is.numeric(bb_place$x)){

        warning("x-coordinate is neither a unit object or a numeric value. Cannot parse x-coordinate.", call. = FALSE)
        bb_place$x <- NULL

      }

      if (is.null(bb_place$default.units)){

        warning("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
        bb_place$x <- NULL

      }

      bb_place$x <- unit(bb_place$x, bb_place$default.units)

    }
  }

  if (!is.null(bb_place$y)){

    if (!"unit" %in% class(bb_place$y)){

      if (!is.numeric(bb_place$y)){

        warning("y-coordinate is neither a unit object or a numeric value. Cannot parse y-coordinate.", call. = FALSE)
        bb_place$y <- NULL

      }

      if (is.null(bb_place$default.units)){

        warning("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
        bb_place$y <- NULL

      }

      bb_place$y <- unit(bb_place$y, bb_place$default.units)

    }
  }

  if (!is.null(bb_place$width)){

    if (!"unit" %in% class(bb_place$width)){

      if (!is.numeric(bb_place$width)){

        warning("width is neither a unit object or a numeric value. Cannot parse width.", call. = FALSE)
        bb_place$width <- NULL

      }

      if (is.null(bb_place$default.units)){

        warning("width detected as numeric.\'default.units\' must be specified.", call. = FALSE)
        bb_place$width <- NULL

      }

      bb_place$width <- unit(bb_place$width, bb_place$default.units)

    }
  }

  if (!is.null(bb_place$height)){

    if (!"unit" %in% class(bb_place$height)){

      if (!is.numeric(bb_place$height)){

        warning("height is neither a unit object or a numeric value. Cannot parse height.", call. = FALSE)
        bb_place$height <- NULL


      }

      if (is.null(bb_place$default.units)){

        warning("height detected as numeric.\'default.units\' must be specified.", call. = FALSE)
        bb_place$height <- NULL

      }

      bb_place$height <- unit(bb_place$height, bb_place$default.units)

    }
  }

  # ======================================================================================================================================================================================
  # UPDATE DIMENSIONS AND COORDINATES OF PLOT OBJECT BASED ON INPUTS
  # ======================================================================================================================================================================================

  object <- set_values(object = object, x = bb_place$x, y = bb_place$y, width = bb_place$width, height = bb_place$height)
  object$justification <- bb_place$just
  attr(x = object, which = "plotted") <- bb_place$draw

  # ======================================================================================================================================================================================
  # INHERIT DIMENSIONS/COOORDINATES WHERE NULL
  # ======================================================================================================================================================================================

  object <- parse_coordinates(input_plot = bb_place$plot, output_plot = object)

  # ======================================================================================================================================================================================
  # CALL ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = object)

  # ======================================================================================================================================================================================
  # DEFINE A NEW VIEWPORT
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0(gsub(pattern = "[0-9]", replacement = "", x = object$grobs$vp$name), length(grep(pattern = gsub(pattern = "[0-9]", replacement = "", x = object$grobs$vp$name), x = currentViewports)) + 1)


  ## If full placing information isn't provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(object$x) | is.null(object$y) | is.null(object$width) | is.null(object$height)){

    new_vp <- viewport(height = unit(1, "snpc"), width = unit(1, "snpc"),
                       x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                       clip = "on",
                       xscale = object$grobs$vp$xscale, yscale = object$grobs$vp$yscale,
                       just = "center",
                       name = vp_name)

    if (bb_place$draw == TRUE){

      grid.newpage()
      warning("Plot placement will only fill up the graphical device.", call. = FALSE)

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = object)

    ## Make viewport
    new_vp <- viewport(height = page_coords$height, width = page_coords$width,
                       x = page_coords$x, y = page_coords$y,
                       clip = "on",
                       xscale = object$grobs$vp$xscale, yscale = object$grobs$vp$yscale,
                       just = bb_place$just,
                       name = vp_name)
  }

  # ======================================================================================================================================================================================
  # RENAME GROBS
  # ======================================================================================================================================================================================

  assign("new_gtree", gTree(vp = new_vp), envir = bbEnv)
  invisible(lapply(object$grobs$children, copy_grobs))

  # ======================================================================================================================================================================================
  # ASSIGN NEW GROBS TO OBJECT
  # ======================================================================================================================================================================================

  object$grobs <- get("new_gtree", envir = bbEnv)

  # ======================================================================================================================================================================================
  # IF DRAW == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_place$draw == TRUE){

    grid.draw(object$grobs)

  }

  # ======================================================================================================================================================================================
  # RETURN UPDATED OBJECT
  # ======================================================================================================================================================================================

  return(object)


}
