#' Handle BentoBox color scaling parameters
#'
#' \code{colorby} should be used to create a set of parameters
#' that specify color scaling for the functions \code{bb_plotBedpe},
#' \code{bb_plotBedpeArches}, and \code{bb_plotBed}.
#'
#' @param column String specifying name of data column to scale colors by.
#' @param range A numeric vector specifying the range of values to
#' apply a color scale to.
#'
#' @return Returns a "\code{bb_colorby}" object.
#' @export
colorby <- function(column, range = NULL){

  colorbyObject <- structure(list(column = column, range = range),
                             class = "bb_colorby")
  return(colorbyObject)

}
