#' plots an APA plot
#'
#' @param apa path to APA txt file or an APA matrix
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param loopNumber number of DNA loops
#' @param palette ColorRamp palette to use for representing interaction scores
#' @param zrange the range of interaction scores to plot, where extreme values will be set to the max or min
#' @param x A numeric or unit object specifying x-location
#' @param y A numeric or unit object specifying y-location
#' @param width A numeric or unit object specifying width
#' @param height A numeric or unit object specifying height
#' @param just justification of x and y-coordinates
#' @param default.units A string indicating the default units to use if x, y, width, or height are only given as numeric vectors
#' @param draw A logical value indicating whether graphics output should be produced

#' @export
bb_plotApa <- function(apa, params = NULL, loopNumber = 1, palette = colorRampPalette(c("white", "dark red")), zrange = NULL,
                       x = NULL, y = NULL, width = NULL, height = NULL, just = c("center"), default.units = "inches", draw = TRUE){

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
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(loopNumber)) loopNumber <- NULL
  if(missing(palette)) palette <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if apa argument is missing (could be in object)
  if(!hasArg(apa)) apa <- NULL

  ## Compile all parameters into an internal object
  bb_apaInternal <- structure(list(apa = apa, loopNumber = loopNumber, palette = palette, zrange = zrange,
                                   x = x, y = y, width = width, just = just, default.units = default.units, draw = draw), class = "bb_apaInternal")

  bb_apaInternal <- parseParams(bb_params = params, object_params = bb_apaInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_apaInternal$loopNumber)) bb_apaInternal$loopNumber <- 1
  if(is.null(bb_apaInternal$palette)) bb_apaInternal$palette <- colorRampPalette(c("white", "dark red"))
  if(is.null(bb_apaInternal$just)) bb_apaInternal$just <- c("center")
  if(is.null(bb_apaInternal$default.units)) bb_apaInternal$default.units <- "inches"
  if(is.null(bb_apaInternal$draw)) bb_apaInternal$draw <- TRUE

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  apa_plot <- structure(list(x = bb_apaInternal$x, y = bb_apaInternal$y, width = bb_apaInternal$width, height = bb_apaInternal$height, justification = bb_apaInternal$just,
                             zrange = bb_apaInternal$zrange, color_palette = NULL, grobs = NULL), class = "bb_apa")
  attr(x = apa_plot, which = "plotted") <- bb_apaInternal$draw

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_apaInternal$apa)) stop("argument \"apa\" is missing, with no default.", call. = FALSE)

  check_placement(object = apa_plot)
  errorcheck_bb_plotApa(apa = bb_apaInternal$apa, apa_plot = apa_plot)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  apa_plot <- defaultUnits(object = apa_plot, default.units = bb_apaInternal$default.units)

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  data <- read_apa(apa = bb_apaInternal$apa)

  # ======================================================================================================================================================================================
  # GET AVERAGE LOOP STRENGTH
  # ======================================================================================================================================================================================

  ## Set scale to 0; divide aggregate values by number of loops to get average loop strength
  data <- data - min(data)
  data <- data/bb_apaInternal$loopNumber

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

    colors <- matrix(bb_maptocolors(data, col = bb_apaInternal$palette, num = 100, range = apa_plot$zrange), nrow = 21, ncol = 21)
    sorted_colors <- bb_maptocolors(unique(data[order(data)]), col = bb_apaInternal$palette, num = 100, range = apa_plot$zrange)
    apa_plot$color_palette <- sorted_colors

  } else {
    colors <- NULL
  }

  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_apa", length(grep(pattern = "bb_apa", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(apa_plot$x) & is.null(apa_plot$y)){

    vp <- viewport(height = unit(1, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   just = "center",
                   name = vp_name)

    if (bb_apaInternal$draw == TRUE){

      vp$name <- "bb_apa1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = apa_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   just = bb_apaInternal$just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # GTREE
  # ======================================================================================================================================================================================

  assign("apa_grobs", gTree(vp = vp), envir = bbEnv)

  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================

  if (!is.null(colors)){

    ## Make raster grob
    apa_grob <- rasterGrob(colors, interpolate = F)
    assign("apa_grobs", addGrob(gTree = get("apa_grobs", envir = bbEnv), child = apa_grob), envir = bbEnv)

  } else {

    warning("APA data could not be plotted.", call. = FALSE)
  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  apa_plot$grobs <- get("apa_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_apaInternal$draw == TRUE){

    grid.draw(apa_plot$grobs)

  }

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(apa_plot)

}
