
#' wrapper to draw a segmentsGrob based on BentoBox page coordinates and units
#'
#' @param x0 Numeric indicating the starting x-values of the line segments
#' @param y0 Numeric indicating the starting y-values of the line segments
#' @param x1 Numeric indicating the stopping x-values of the line segments
#' @param y1 Numeric indicating the stopping y-values of the line segments
#' @param params an optional "bb_params" object space containing relevant function parameters
#' @param arrow A list describing arrow heads to place at either end of the line segments, as produced by the arrow function
#' @param linecolor line color
#' @param lwd line width
#' @param lty line type
#' @param lineend Line end style (round, butt, square)
#' @param linejoin Line join style (round, mitre, bevel)
#' @param alpha color transparency 
#' @param default.units A string indicating the default units to use if x or y are only given as numeric vectors
#'
#' @export
bb_plotSegments <- function(x0, y0, x1, y1, params = NULL, arrow=NULL, linecolor = "black", lwd = 1, lty = 1, lineend = "butt", linejoin = "mitre", alpha = 1, default.units = "inches", ...){
  
  
  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================
  
  ## Check which defaults are not overwritten and set to NULL
  if(missing(arrow)) arrow <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(lwd)) lwd <- NULL
  if(missing(lty)) lty <- NULL
  if(missing(lineend)) lineend <- NULL
  if(missing(linejoin)) linejoin <- NULL
  if(missing(alpha)) alpha <- NULL
  if(missing(default.units)) default.units <- NULL
  
  ## Check if x0/y0/x1/y1 arguments are missing (could be in object)
  if(!hasArg(x0)) x0 <- NULL
  if(!hasArg(y0)) y0 <- NULL
  if(!hasArg(x1)) x1 <- NULL
  if(!hasArg(y1)) y1 <- NULL
  
  ## Compile all parameters into an internal object
  bb_segmentsInternal <- structure(list(x0 = x0, y0 = y0, x1 = x1, y0 = y0, y1 = y1, arrow = arrow, linecolor = linecolor,
                                      lwd = lwd, lty = lty, lineend = lineend, linejoin = linejoin, alpha = alpha, default.units = default.units), class = "bb_segmentsInternal")
  
  bb_segmentsInternal <- parseParams(bb_params = params, object_params = bb_segmentsInternal)
  
  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_segmentsInternal$arrow)) bb_segmentsInternal$arrow <- NULL
  if(is.null(bb_segmentsInternal$linecolor)) bb_segmentsInternal$linecolor <- "black"
  if(is.null(bb_segmentsInternal$lwd)) bb_segmentsInternal$lwd <- 1
  if(is.null(bb_segmentsInternal$lty)) bb_segmentsInternal$lty <- 1
  if(is.null(bb_segmentsInternal$lineend)) bb_segmentsInternal$lineend <- "butt"
  if(is.null(bb_segmentsInternal$linejoin)) bb_segmentsInternal$linejoin <- "mitre"
  if(is.null(bb_segmentsInternal$alpha)) bb_segmentsInternal$alpha <- 1
  if(is.null(bb_segmentsInternal$default.units)) bb_segmentsInternal$default.units <- "inches"
  
  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================
  
  bb_segments <- structure(list(x0 = bb_segmentsInternal$x0, y0 = bb_segmentsInternal$y0, x1 = bb_segmentsInternal$x1, y1 = bb_segmentsInternal$y1, 
                                arrow = bb_segmentsInternal$arrow, grobs = NULL, 
                                gp = gpar(col = bb_segmentsInternal$linecolor, lwd = bb_segmentsInternal$lwd,lty = bb_segmentsInternal$lty, lineend = bb_segmentsInternal$lineend, 
                                          linejoin = bb_segmentsInternal$linejoin, alpha = bb_segmentsInternal$alpha, ...)), class = "bb_segmentsInternal")
  
  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================
  
  check_bbpage(error = "Cannot plot segment without a BentoBox page.")
  if(is.null(bb_segments$x0)) stop("argument \"x0\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_segments$y0)) stop("argument \"y0\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_segments$x1)) stop("argument \"x1\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_segments$y1)) stop("argument \"y1\" is missing, with no default.", call. = FALSE)

  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================
  
  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)
  
  if (!"unit" %in% class(bb_segments$x0)){
    
    if (!is.numeric(bb_segments$x0)){
      
      stop("Starting x-coordinate is neither a unit object or a numeric value. Cannot plot segment.", call. = FALSE)
      
    }
    
    if (is.null(bb_segmentsInternal$default.units)){
      
      stop("Starting x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_segments$x0 <- unit(bb_segments$x0, bb_segmentsInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_segments$y0)){
    
    if (!is.numeric(bb_segments$y0)){
      
      stop("Starting y-coordinate is neither a unit object or a numeric value. Cannot plot segment.", call. = FALSE)
      
    }
    
    if (is.null(bb_segmentsInternal$default.units)){
      
      stop("Starting y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_segments$y0 <- unit(bb_segments$y0, bb_segmentsInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_segments$x1)){
    
    if (!is.numeric(bb_segments$x1)){
      
      stop("Stopping x-coordinate is neither a unit object or a numeric value. Cannot plot segment.", call. = FALSE)
      
    }
    
    if (is.null(bb_segmentsInternal$default.units)){
      
      stop("Stopping x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_segments$x1 <- unit(bb_segments$x1, bb_segmentsInternal$default.units)
    
  }
  
  if (!"unit" %in% class(bb_segments$y1)){
    
    if (!is.numeric(bb_segments$y1)){
      
      stop("Stopping y-coordinate is neither a unit object or a numeric value. Cannot plot segment.", call. = FALSE)
      
    }
    
    if (is.null(bb_segmentsInternal$default.units)){
      
      stop("Stopping y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)
      
    }
    
    bb_segments$y1 <- unit(bb_segments$y1, bb_segmentsInternal$default.units)
    
  }
  ## Convert coordinates to page_units
  new_x0 <- convertX(bb_segments$x0, unitTo = page_units, valueOnly = TRUE)
  new_y0 <- convertY(bb_segments$y0, unitTo = page_units, valueOnly = TRUE)
  new_x1 <- convertX(bb_segments$x1, unitTo = page_units, valueOnly = TRUE)
  new_y1 <- convertY(bb_segments$y1, unitTo = page_units, valueOnly = TRUE)
  
  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================
  
   segments <- grid.segments(x0 = unit(new_x0, page_units), y0 = unit(page_height - new_y0, page_units), x1 = unit(page_height - new_x1, page_units), 
                             y1 = unit(page_height - new_y1, page_units), gp = bb_segments$gp)
  
  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================
  
  bb_segments$grobs <- segments
  
  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================
  
  return(bb_segments)
}
