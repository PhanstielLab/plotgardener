#' plots an Apa plot
#' @param apa path to APA txt file
#' @param loopNumber number of DNA loops
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param zrangethe range of interaction scores to plot, where extreme values will be set to the max or min
#' @param x x-coordinate of plot
#' @param y y-coordinate of plot
#' @param width width of plot
#' @param height height of plot
#' @param just justification of x and y-coordinates
#' @param units units of width, height, x, and y-coordinates

#' @export
bb_plotApa <- function(apa, loopNumber = 1, palette = colorRampPalette(c("white", "dark red")), zrange = NULL, x = 1, y = 1, width = 1, height = 1, just = c("center"),
                       units = "inches"){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a functino that catches errors for bb_plotApa
  errorcheck_bb_plotApa <- function(apa, apa_plot){

    if (file_ext(apa) != "txt"){

      stop("Invalid input. APA file must have a \".txt\" extension.")

    }


    if (apa_plot$width != apa_plot$height){

      warning("Width must equal height for a square APA plot.")

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

    ## Read in data from apa path
    data <- data.table::fread(apa)

    ## Remove brackets and convert to numeric
    data <- apply(data, 2, function(x) gsub("\\[|\\]", "", x) %>% as.numeric())

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

  apa_plot <- structure(list(x = x, y = y, width = width, height = height, units = units, justification = just,
                             zrange = zrange, color_palette = NULL, grobs = NULL, viewport = NULL), class = "bb_apa")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage()
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

  ## Convert coordinates into same units as page
  page_coords <- convert_page(object = apa_plot)

  ## Name viewport
  current_viewports <- lapply(current.vpTree()$children$bb_page$children, viewport_name)
  vp_name <- paste0("bb_apa", length(grep(pattern = "bb_apa", x = current_viewports)) + 1)

  ## Make viewport
  vp <- viewport(height = unit(page_coords[[1]]$height, page_coords[[3]]), width = unit(page_coords[[1]]$width, page_coords[[3]]),
                  x = unit(page_coords[[1]]$x, page_coords[[3]]), y = unit((page_coords[[2]]-page_coords[[1]]$y), page_coords[[3]]),
                  just = just, name = vp_name)
  apa_plot$viewport <- vp
  pushViewport(vp)

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  ## Make raster plot
  apa_grob <- rasterGrob(colors, interpolate = F)
  grid.draw(apa_grob)

  ## Go back to root viewport
  upViewport()

  ## Add grob to object
  apa_plot$grobs <- gList(apa_grob)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(apa_plot)

}
