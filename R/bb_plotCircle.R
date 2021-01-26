
#' wrapper to draw a circleGrob based on BentoBox page coordinates and units
#'
#' @param x A numeric vector or unit object specifying x-location of the center of the circle
#' @param y A numeric vector or unit object specifying y-location of the center of the circle
#' @param r A numeric vector or unit object specifying radii
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param linecolor line color
#' @param fill fill color
#' @param lwd line width
#' @param lty line type
#' @param alpha color transparency 
#' @param default.units A string indicating the default units to use if x or y are only given as numeric vectors
#'
#' @export
bb_plotCircle <- function(x, y, r, params = NULL, linecolor = "black", fill = NA, lwd = 1, lty = 1, alpha = 1, default.units = "inches", ...){
  
  
  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================
  
  ## Check which defaults are not overwritten and set to NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(fill)) fill <- NULL
  if(missing(lwd)) lwd <- NULL
  if(missing(lty)) lty <- NULL
  if(missing(alpha)) alpha <- NULL
  if(missing(default.units)) default.units <- NULL
  
  ## Check if label/x/y arguments are missing (could be in object)
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL
  if(!hasArg(r)) r <- NULL
  
  ## Compile all parameters into an internal object
  bb_circleInternal <- structure(list(x = x, y = y, r = r, linecolor = linecolor, fill = fill,
                                      lwd = lwd, lty = lty, alpha = alpha, default.units = default.units), class = "bb_circleInternal")
  
  bb_circleInternal <- parseParams(bb_params = params, object_params = bb_circleInternal)
  
  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_circleInternal$linecolor)) bb_circleInternal$linecolor <- "black"
  if(is.null(bb_circleInternal$fill)) bb_circleInternal$fill <- NA
  if(is.null(bb_circleInternal$lwd)) bb_circleInternal$lwd <- 1
  if(is.null(bb_circleInternal$lty)) bb_circleInternal$lty <- 1
  if(is.null(bb_circleInternal$alpha)) bb_circleInternal$alpha <- 1
  if(is.null(bb_circleInternal$default.units)) bb_circleInternal$default.units <- "inches"
  
  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================
  
  bb_circle <- structure(list(x = bb_circleInternal$x, y = bb_circleInternal$y, r = bb_circleInternal$r, grobs = NULL,
                              gp = gpar(col = bb_circleInternal$linecolor, fill = bb_circleInternal$fill, lwd = bb_circleInternal$lwd, 
                                        lty = bb_circleInternal$lty, alpha = bb_circleInternal$alpha, ...)), class = "bb_circle")
  
  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================
  
  check_bbpage(error = "Cannot plot circle without a BentoBox page.")
  if(is.null(bb_circle$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_circle$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_circle$r)) stop("argument \"r\" is missing, with no default.", call. = FALSE)
  
  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================
  
  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)
  
  if (!"unit" %in% class(bb_circle$x)){
    
    if (!is.numeric(bb_circle$x)){
      
      stop("x-coordinate is neither a unit object or a numeric value. Cannot plot circle.", call. = FALSE)
      
    }
    
    if (is.null(bb_circleInternal$default.units)){
      
      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_circle$x <- unit(bb_circle$x, bb_circleInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_circle$y)){
    
    if (!is.numeric(bb_circle$y)){
      
      stop("y-coordinate is neither a unit object or a numeric value. Cannot plot circle.", call. = FALSE)
      
    }
    
    if (is.null(bb_circleInternal$default.units)){
      
      stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_circle$y <- unit(bb_circle$y, bb_circleInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_circle$r)){
    
    if (!is.numeric(bb_circle$r)){
      
      stop("Radius is neither a unit object or a numeric value. Cannot plot circle.", call. = FALSE)
      
    }
    
    if (is.null(bb_circleInternal$default.units)){
      
      stop("Radius detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_circle$r <- unit(bb_circle$r, bb_circleInternal$default.units)
    
  }
  
  ## Convert coordinates to page_units
  new_x <- convertX(bb_circle$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- convertY(bb_circle$y, unitTo = page_units, valueOnly = TRUE)
  
  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================
  
  circle <- grid.circle(x = unit(new_x, page_units), y = unit(page_height - new_y, page_units), r = bb_circle$r,
                        gp = bb_circle$gp)
  
  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================
  
  bb_circle$grobs <- circle
  
  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================
  
  return(bb_circle)
}



