
#' wrapper to draw a rasterGrob based on BentoBox page coordinates and units
#'
#' @param image Any R object that can be coerced to a raster object
#' @param x A numeric vector or unit object specifying x-location
#' @param y A numeric vector or unit object specifying y-location
#' @param width A numeric vector or unit object specifying width
#' @param height A numeric vector or unit object specifying height
#' @param just The justification of the raster relative to its (x,y) location. If there are 2 values, the first specifies horizontal justification and the second specifies vertical justification. Options: "left","right","center","bottom", and "top". Numerically, 0 means left alignment and 1 means right alignment.
#' @param interpolate A logical value indicating whether to linearly interpolate the image
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param default.units A string indicating the default units to use if x or y are only given as numeric vectors
#'
#' @export
bb_plotRaster <- function(image, x, y, width, height, just = "center", interpolate = TRUE, params = NULL, default.units = "inches", ...){
  
  
  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================
  
  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(interpolate)) interpolate <- NULL
  if(missing(default.units)) default.units <- NULL
  
  ## Check if image/x/y arguments are missing (could be in object)
  if(!hasArg(image)) image <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(width)) width <- NULL
  if(!hasArg(height)) height <- NULL
  
  ## Compile all parameters into an internal object
  bb_rastInternal <- structure(list(image = image, x = x, y = y, width = width, height = height, just = just, interpolate = interpolate, default.units = default.units), class = "bb_rastInternal")
  
  bb_rastInternal <- parseParams(bb_params = params, object_params = bb_rastInternal)
  
  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_rastInternal$just)) bb_rastInternal$just <- "center"
  if(is.null(bb_rastInternal$interpolate)) bb_rastInternal$interpolate <- TRUE
  if(is.null(bb_rastInternal$default.units)) bb_rastInternal$default.units <- "inches"
  
  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================
  
  bb_rast <- structure(list(image = bb_rastInternal$image, x = bb_rastInternal$x, y = bb_rastInternal$y, width = bb_rastInternal$width, height = bb_rastInternal$height, just = bb_rastInternal$just, 
                            interpolate = bb_rastInternal$interpolate, grobs = NULL, gp = gpar()), class = "bb_rast")
  
  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================
  
  check_bbpage(error = "Cannot plot raster without a BentoBox page.")
  if(is.null(bb_rast$image)) stop ("arguement \"image\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rast$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rast$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rast$width)) stop("argument \"width\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rast$height)) stop("argument \"height\" is missing, with no default.", call. = FALSE)
  
  
  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================
  
  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)
  
  if (!"unit" %in% class(bb_rast$x)){
    
    if (!is.numeric(bb_rast$x)){
      
      stop("x-coordinate is neither a unit object or a numeric value. Cannot plot raster.", call. = FALSE)
      
    }
    
    if (is.null(bb_rastInternal$default.units)){
      
      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_rast$x <- unit(bb_rast$x, bb_rastInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_rast$y)){
    
    if (!is.numeric(bb_rast$y)){
      
      stop("y-coordinate is neither a unit object or a numeric value. Cannot plot raster.", call. = FALSE)
      
    }
    
    if (is.null(bb_rastInternal$default.units)){
      
      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_rast$y <- unit(bb_rast$y, bb_rastInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_rast$width)){
    
    if (!is.numeric(bb_rast$width)){
      
      stop("width is neither a unit object or a numeric value. Cannot plot raster.", call. = FALSE)
      
    }
    
    if (is.null(bb_rastInternal$default.units)){
      
      stop("width detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_rast$width <- unit(bb_rast$width, bb_rastInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_rast$height)){
    
    if (!is.numeric(bb_rast$height)){
      
      stop("height is neither a unit object or a numeric value. Cannot plot raster.", call. = FALSE)
      
    }
    
    if (is.null(bb_rastInternal$default.units)){
      
      stop("height detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_rast$height <- unit(bb_rast$height, bb_rastInternal$default.units)
    
  }
  ## Convert coordinates to page_units
  new_x <- convertX(bb_rast$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(bb_rast$y, unitTo = page_units, valueOnly = TRUE)
  
  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================
  
  rast <- grid.raster(image = bb_rast$image, x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), width = bb_rast$width, height = bb_rast$height, 
                    just = bb_rast$just, interpolate = bb_rast$interpolate, gp = bb_rast$gp)
  
  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================
  
  bb_rast$grobs <- rast
  
  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================
  
  return(bb_rast)
}



