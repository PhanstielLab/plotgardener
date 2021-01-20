
#' wrapper to draw a rectGrob based on BentoBox page coordinates and units
#'
#' @param x A numeric vector or unit object specifying x-location
#' @param y A numeric vector or unit object specifying y-location
#' @param width A numeric vector or unit object specifying width
#' @param height A numeric vector or unit object specifying height
#' @param just The justification of the rectangle relative to its (x,y) location
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param linecolor line color
#' @param fill fill color
#' @param lwd line width
#' @param lty line type
#' @param alpha color transparency 
#' @param default.units A string indicating the default units to use if x or y are only given as numeric vectors
#'
#' @export
bb_plotRect <- function(x, y, width, height, just = "center", params = NULL, linecolor = "black", fill = NA, lwd = 1, lty = 1, alpha = 0, default.units = "inches", ...){
  
  
  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================
  
  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(fill)) fill <- NULL
  if(missing(lwd)) lwd <- NULL
  if(missing(lty)) lty <- NULL
  if(missing(alpha)) alpha <- NULL
  if(missing(default.units)) default.units <- NULL
  
  ## Check if label/x/y arguments are missing (could be in object)
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(width)) width <- NULL
  if(!hasArg(height)) height <- NULL
  
  ## Compile all parameters into an internal object
  bb_rectInternal <- structure(list(x = x, y = y, width = width, height = height, just = just, linecolor = linecolor, fill = fill,
                                      lwd = lwd, lty = lty, alpha = alpha, default.units = default.units), class = "bb_rectInternal")
  
  bb_rectInternal <- parseParams(bb_params = params, object_params = bb_rectInternal)
  
  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_rectInternal$just)) bb_rectInternal$just <- "center"
  if(is.null(bb_rectInternal$linecolor)) bb_rectInternal$linecolor <- "black"
  if(is.null(bb_rectInternal$fill)) bb_rectInternal$fill <- NA
  if(is.null(bb_rectInternal$lwd)) bb_rectInternal$lwd <- 1
  if(is.null(bb_rectInternal$lty)) bb_rectInternal$lty <- 1
  if(is.null(bb_rectInternal$alpha)) bb_rectInternal$alpha <- 0
  if(is.null(bb_rectInternal$default.units)) bb_rectInternal$default.units <- "inches"
  
  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================
  
  bb_rect <- structure(list(x = bb_rectInternal$x, y = bb_rectInternal$y, width = bb_rectInternal$width, height = bb_rectInternal$height, just = bb_rectInternal$just, 
                            grobs = NULL, gp = gpar(col = bb_rectInternal$linecolor, fill = bb_rectInternal$fill, lwd = bb_rectInternal$lwd, lty = bb_rectInternal$lty, 
                                                    alpha = bb_rectInternal$alpha, ...)), class = "bb_rect")
  
  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================
  
  check_bbpage(error = "Cannot plot rectangle without a BentoBox page.")
  if(is.null(bb_rect$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rect$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rect$width)) stop("argument \"width\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_rect$height)) stop("argument \"height\" is missing, with no default.", call. = FALSE)
  
  
  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================
  
  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)
  
  if (!"unit" %in% class(bb_rect$x)){
    
    if (!is.numeric(bb_rect$x)){
      
      stop("x-coordinate is neither a unit object or a numeric value. Cannot plot rectangle.", call. = FALSE)
      
    }
    
    if (is.null(bb_rectInternal$default.units)){
      
      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_rect$x <- unit(bb_rect$x, bb_rectInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_rect$y)){
    
    if (!is.numeric(bb_rect$y)){
      
      stop("y-coordinate is neither a unit object or a numeric value. Cannot plot rectangle.", call. = FALSE)
      
    }
    
    if (is.null(bb_rectInternal$default.units)){
      
      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_rect$y <- unit(bb_rect$y, bb_rectInternal$default.units)
    
  }
  
  ## Convert coordinates to page_units
  new_x <- convertX(bb_rect$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(bb_rect$y, unitTo = page_units, valueOnly = TRUE)
  
  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================
  
  rect <- grid.rect(x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), width = bb_rect$width, height = bb_rect$height, 
                        just = bb_rect$just, gp = bb_rect$gp)
  
  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================
  
  bb_rect$grobs <- rect
  
  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================
  
  return(bb_rect)
}



