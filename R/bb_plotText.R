#' Plot text within a BentoBox layout
#'
#' @param label Character or expression of text to be plotted.
#' @param fontcolor A character value specifying text fontcolor.
#' Default value is \code{fontcolor = "black"}.
#' @param fontsize A numeric specifying text fontsize in points.
#' Default value is \code{fontsize = 12}.
#' @param rot A numeric specifying the angle to rotate the text.
#' Default value is \code{rot = 0}.
#' @param check.overlap A logical value to indicate whether to check
#' for and omit overlapping text.
#' Default value is \code{check.overlap = FALSE}.
#' @param x A numeric vector or unit object specifying text x-location.
#' @param y A numeric vector, unit object, or a character vector of
#' values containing a "b" combined with a numeric value
#' specifying text y-locations.
#' The character vector will place text y-locations relative to the
#' bottom of the most recently plotted BentoBox plot according to the
#' units of the BentoBox page.
#' @param just Justification of text relative to its (x, y) location.
#' If there are two values, the first value specifies horizontal
#' justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"},
#' \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}.
#' Default value is \code{just = "center"}.
#' @param default.units A string indicating the default units to use if
#' \code{x} or \code{y} are only given as numerics or numeric vectors.
#' Default value is \code{default.units = "inches"}.
#' @param params An optional \link[BentoBox]{bb_params} object
#' containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_text} object containing relevant
#' placement and \link[grid]{grob} information.
#'
#' @examples
#' ## Create a BentoBox page
#' bb_pageCreate(width = 4, height = 2, default.units = "inches")
#'
#' ## Plot text, adjusting fontsize and fontface
#' bb_plotText(label = "BentoBox", fontsize = 14, fontface = "bold",
#'             x = 1, y = 1, just = "center", default.units = "inches")
#'
#' ## Plot text, adjusting color, rotation, and fontfamily
#' bb_plotText(label = "coordinate-based", fontcolor = "#225EA8", rot = 90,
#'             fontfamily = "HersheyScript", x = 2, y = 1, just = "center",
#'             default.units = "inches")
#'
#' ## Plot a text label in multiple places at once
#' bb_plotText(label = "R", x = c(0.5, 1, 1.5), y = 1.5, just = "center",
#'             default.units = "inches")
#'
#' ## Plot a vector of text labels
#' bb_plotText(label = c("bb", "Bento", "Box"), x = 3, y = c(0.5, 1, 1.75),
#'             just = "center", default.units = "inches")
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @seealso \link[grid]{grid.text}
#'
#' @export
bb_plotText <- function(label, fontcolor = "black", fontsize = 12, rot = 0,
                        check.overlap = FALSE, x, y, just = "center",
                        default.units = "inches", params = NULL, ...){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(just)) just <- NULL
  if(missing(fontcolor)) fontcolor <- NULL
  if(missing(fontsize)) fontsize <- NULL
  if(missing(rot)) rot <- NULL
  if(missing(check.overlap)) check.overlap <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if label/x/y arguments are missing (could be in object)
  if(!hasArg(label)) label <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL

  ## Compile all parameters into an internal object
  bb_textInternal <- structure(list(label = label, x = x, y = y,
                                    just = just, fontcolor = fontcolor,
                                    fontsize = fontsize, rot = rot,
                                    check.overlap = check.overlap,
                                    default.units = default.units),
                               class = "bb_textInternal")

  bb_textInternal <- parseParams(bb_params = params,
                                 object_params = bb_textInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_textInternal$just)) bb_textInternal$just <- "center"
  if(is.null(bb_textInternal$fontcolor)) bb_textInternal$fontcolor <- "black"
  if(is.null(bb_textInternal$fontsize)) bb_textInternal$fontsize <- 12
  if(is.null(bb_textInternal$rot)) bb_textInternal$rot <- 0
  if(is.null(bb_textInternal$check.overlap)) bb_textInternal$check.overlap <- FALSE
  if(is.null(bb_textInternal$default.units)) bb_textInternal$default.units <- "inches"

  ## Set gp
  bb_textInternal$gp <- gpar(col = bb_textInternal$fontcolor,
                             fontsize = bb_textInternal$fontsize)
  bb_textInternal$gp <- setGP(gpList = bb_textInternal$gp,
                              params = bb_textInternal, ...)

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  bb_text <- structure(list(label = bb_textInternal$label,
                            x = bb_textInternal$x, y = bb_textInternal$y,
                            just = bb_textInternal$just, grobs = NULL),
                       class = "bb_text")

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  check_bbpage(error = "Cannot plot text without a BentoBox page.")
  if(is.null(bb_text$label)) stop("argument \"label\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_text$x)) stop("argument \"x\" is missing, with no default.",
                              call. = FALSE)
  if(is.null(bb_text$y)) stop("argument \"y\" is missing, with no default.",
                              call. = FALSE)

  # ======================================================================================================================================================================================
  # DEFINE PARAMETERS
  # ======================================================================================================================================================================================

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  if (!"unit" %in% class(bb_text$x)){

    if (!is.numeric(bb_text$x)){

      stop("x-coordinate is neither a unit object or a numeric value. Cannot plot text.", call. = FALSE)

    }

    if (is.null(bb_textInternal$default.units)){

      stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

    }

    bb_text$x <- unit(bb_text$x, bb_textInternal$default.units)

  }

  if (!"unit" %in% class(bb_text$y)){

    ## Check for "below" y-coords
    if (all(grepl("b", bb_text$y)) == TRUE){
      if (any(grepl("^[ac-zA-Z]+$", bb_text$y)) == TRUE){
        stop("\'below\' y-coordinate(s) detected with additional letters. Cannot parse y-coordinate(s).", call. = FALSE)
      }

      if(any(is.na(as.numeric(gsub("b","", bb_text$y))))){
        stop("\'below\' y-coordinate(s) does not have a numeric associated with it. Cannot parse y-coordinate(s).", call. = FALSE)
      }

      bb_text$y <- unit(unlist(lapply(bb_text$y, plot_belowY)),
                        get("page_units", envir = bbEnv))

    } else {

      if (!is.numeric(bb_text$y)){

        stop("y-coordinate is neither a unit object or a numeric value. Cannot plot text.", call. = FALSE)

      }

      if (is.null(bb_textInternal$default.units)){

        stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      bb_text$y <- unit(bb_text$y, bb_textInternal$default.units)

    }

  }

  ## Convert coordinates to page_units
  new_x <- convertX(bb_text$x, unitTo = page_units, valueOnly = TRUE)
  new_y <- page_height - convertY(bb_text$y, unitTo = page_units,
                                  valueOnly = TRUE)

  # ======================================================================================================================================================================================
  # MAKE GROB
  # ======================================================================================================================================================================================

  name <- paste0("bb_text",
                 length(grep(pattern = "bb_text",
                             x = grid.ls(print = FALSE,
                                         recursive = FALSE))) + 1)
  text <- grid.text(label = bb_text$label, x = unit(new_x, page_units),
                    y = unit(new_y, page_units), just = bb_text$just,
                    gp = bb_textInternal$gp, rot = bb_textInternal$rot,
                    check.overlap = bb_textInternal$check.overlap,
                    name = name)

  # ======================================================================================================================================================================================
  # ADD GROB TO OBJECT
  # ======================================================================================================================================================================================

  bb_text$grobs <- text

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_text[", text$name, "]"))
  invisible(bb_text)
}
