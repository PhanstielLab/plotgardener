#' places a plot that has been previously created but not drawn
#'
#' @param plot plot to be placed
#' @param x A unit object specifying x-location.
#' @param y A unit object specifying y-location.
#' @param width A unit object specifying width.
#' @param height A unit object specifying height.
#' @param just A string or numeric vector specifying the justification of the plot relative to its (x, y) location
#'
#' @return Function will plot a HiC interaction matrix and return a bb_hicPlot object
#'
#'
#' @export
bb_placePlot <- function(plot, x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"), draw = T){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================
  ## Error function
  # if you give an x, need to give a y and vice versa

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
  # INITIALIZE PLOT OBJECT COPY
  # ======================================================================================================================================================================================

  object <- plot

  # ======================================================================================================================================================================================
  # UPDATE DIMENSIONS AND COORDINATES OF PLOT OBJECT BASED ON INPUTS
  # ======================================================================================================================================================================================
  object <- set_values(object = object, x = x, y = y, width = width, height = height)
  object$justification <- just
  attr(x = object, which = "plotted") <- draw

  # ======================================================================================================================================================================================
  # INHERIT DIMENSIONS/COOORDINATES WHERE NULL
  # ======================================================================================================================================================================================

  object <- parse_coordinates(input_plot = plot, output_plot = object)

  # ======================================================================================================================================================================================
  # CALL ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = object)

  # ======================================================================================================================================================================================
  # DEFINE A NEW VIEWPORT
  # ======================================================================================================================================================================================

  ## Get viewport name
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0(gsub(pattern = "[0-9]", replacement = "", x = object$grobs$vp$name), length(grep(pattern = gsub(pattern = "[0-9]", replacement = "", x = object$grobs$vp$name), x = current_viewports)) + 1)


  ## If full placing information isn't provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(object$x) | is.null(object$y) | is.null(object$width) | is.null(object$height)){

    new_vp <- viewport(height = unit(1, "snpc"), width = unit(1, "snpc"),
                       x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                       clip = "on",
                       xscale = object$grobs$vp$xscale, yscale = object$grobs$vp$yscale,
                       just = "center",
                       name = vp_name)

    if (draw == TRUE){

      grid.newpage()
      warning("Plot placement will only fill up the graphical device.")

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = object)

    ## Make viewport
    new_vp <- viewport(height = page_coords$height, width = page_coords$width,
                       x = page_coords$x, y = page_coords$y,
                       clip = "on",
                       xscale = object$grobs$vp$xscale, yscale = object$grobs$vp$yscale,
                       just = just,
                       name = vp_name)
  }





  # ## Convert coordinates into same units as page
  # page_coords <- convert_page(object = object)
  #
  #
  # ## Make viewport
  # new_vp <- viewport(height = page_coords$height, width = page_coords$width,
  #                x = page_coords$x, y = page_coords$y,
  #                clip = "on",
  #                xscale = object$grobs$vp$xscale, yscale = object$grobs$vp$yscale,
  #                just = just,
  #                name = vp_name)

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

  if (draw == TRUE){

    grid.draw(object$grobs)

  }

  # ======================================================================================================================================================================================
  # RETURN UPDATED OBJECT
  # ======================================================================================================================================================================================

  return(object)


}
