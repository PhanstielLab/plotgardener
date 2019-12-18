#' plots an Apa plot
#' @param apa path to APA txt file or an APA matrix
#' @param loopNumber number of DNA loops
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param zrangethe range of interaction scores to plot, where extreme values will be set to the max or min
#' @param x A unit object specifying x-location
#' @param y A unit object specifying y-location
#' @param width A unit object specifying width
#' @param height A unit object specifying height
#' @param just justification of x and y-coordinates
#' @param draw A logical value indicating whether graphics output should be produced

#' @export
bb_plotApa <- function(apa, loopNumber = 1, palette = colorRampPalette(c("white", "dark red")), zrange = NULL, x = NULL, y = NULL, width = NULL, height = NULL, just = c("center"),
                       draw = TRUE){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a functino that catches errors for bb_plotApa
  errorcheck_bb_plotApa <- function(apa, apa_plot){

    if (class(apa) != "matrix"){

      if (file_ext(apa) != "txt" ){

        stop("Invalid input. APA file must have a \".txt\" extension or be a matrix.")

      }

    }



    # zrange
    # ===================================================================================
    ## Ensure properly formatted zrange
    if (!is.null(apa_plot$zrange)){

      ## zrange needs to be a vector
      if (!is.vector(apa_plot$zrange)){

        stop("\'zrange\' must be a vector of length 2.")

      }

      ## zrange vector needs to be length 2
      if (length(apa_plot$zrange) != 2){

        stop("\'zrange\' must be a vector of length 2.")

      }

      ## zrange vector needs to be numbers
      if (!is.numeric(apa_plot$zrange)){

        stop("\'zrange\' must be a vector of two numbers.")

      }

      ## second value should be larger than the first value
      if (apa_plot$zrange[1] >= apa_plot$zrange[2]){

        stop("\'zrange\' must be a vector of two numbers in which the 2nd value is larger than the 1st.")

      }

    }

  }

  ## Define a function that reads in an apa file
  read_apa <- function(apa){

    if (class(apa) != "matrix"){

      ## Read in data from apa path
      data <- data.table::fread(apa)

      ## Remove brackets and convert to numeric
      data <- apply(data, 2, function(x) gsub("\\[|\\]", "", x) %>% as.numeric())

    }

    return(data)
  }

  ## Define a function that sets the apa zrange
  apa_zrange <- function(apa, apa_plot){

    ## no zrange, only one value
    if (is.null(apa_plot$zrange) & length(unique(apa)) == 1){

      zrange <- c(unique(apa), unique(apa))
      apa_plot$zrange <- zrange

    }

    ## no zrange, multiple values
    if (is.null(apa_plot$zrange) & length(unique(apa)) > 1){

      zrange <- c(min(apa), max(apa))
      apa_plot$zrange <- zrange

    }

    return(apa_plot)

  }

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  apa_plot <- structure(list(x = x, y = y, width = width, height = height, justification = just,
                             zrange = zrange, color_palette = NULL, grobs = NULL), class = "bb_apa")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_placement(object = apa_plot)
  errorcheck_bb_plotApa(apa = apa, apa_plot = apa_plot)

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  data <- read_apa(apa = apa)

  # ======================================================================================================================================================================================
  # GET AVERAGE LOOP STRENGTH
  # ======================================================================================================================================================================================

  ## Set scale to 0; divide aggregate values by number of loops to get average loop strength
  data <- data - min(data)
  data <- data/loopNumber

  # ======================================================================================================================================================================================
  # SET ZRANGE AND SCALE DATA
  # ======================================================================================================================================================================================

  apa_plot <- apa_zrange(apa = data, apa_plot = apa_plot)
  data[data <= apa_plot$zrange[1]] <- apa_plot$zrange[1]
  data[data >= apa_plot$zrange[2]] <- apa_plot$zrange[2]

  # ======================================================================================================================================================================================
  # CONVERT NUMBERS TO COLORS
  # ======================================================================================================================================================================================

  ## if we don't have an appropriate zrange (even after setting it based on a null zrange), can't scale to colors
  if (!is.null(apa_plot$zrange) & length(unique(apa_plot$zrange)) == 2){

    colors <- matrix(bb_maptocolors(data, col = palette, num = 100, range = apa_plot$zrange), nrow = 21, ncol = 21)
    sorted_colors <- bb_maptocolors(unique(data[order(data)]), col = palette, num = 100, range = apa_plot$zrange)
    apa_plot$color_palette <- sorted_colors

  } else {

    ## If we still have a null zrange or a length(unique(zrange)) == 1, means we couldn't do it in setzrange above (empty data or data with only 1 value)
    warning("Can't scale data to colors.")

  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_apa", length(grep(pattern = "bb_apa", x = current_viewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(x) & is.null(y)){

    vp <- viewport(height = unit(1, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   just = "center",
                   name = vp_name)

    if (draw == TRUE){

      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = apa_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   just = just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================

  ## Make raster grob
  apa_grob <- rasterGrob(colors, interpolate = F)
  apa_grobs <- gTree(vp = vp, children = gList(apa_grob))

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (draw == TRUE){

    grid.draw(apa_grobs)

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  apa_plot$grobs <- apa_grobs

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(apa_plot)

}
