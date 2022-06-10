#' Generate row positions for a number of plot elements with a specified height
#' and space between them
#' 
#' @param y A numeric or unit object specifying the starting row y-position.
#' @param height A numeric or unit object specifying the height of rows.
#' @param space A numeric or unit object specifying the space between rows.
#' @param n An integer specifying the number of elements to generate row 
#' positions for.
#' @param default.units A string indicating the default units to use
#' if \code{y}, \code{h}, or \code{s} are only given as numerics.
#' Default value is \code{default.units = "inches"}.
#' 
#' @examples 
#' # Starting at 0.5 units, return a vector of positions for 3 objects that
#' # are 2 units in height with 0.1 units of space between them
#' 
#' pageLayoutRow(y = 0.5, height = 2, space = 0.1, n = 3, 
#'             default.units = "inches")
#' 
#' @returns Returns a unit vector of page positions.
#' @export
pageLayoutRow <- function(y, height, space, n, default.units = "inches"){
    ## Parse default units
    y <- misc_defaultUnits(value = y,
                           name = "y",
                           default.units = default.units)
    height <- misc_defaultUnits(value = height,
                           name = "height",
                           default.units = default.units)
    space <- misc_defaultUnits(value = space,
                           name = "space",
                           default.units = default.units)
    
    ## Calculate positions
    y + ((height+space) * seq(0,(n-1)))
    
}

#' Generate column positions for a number of plot elements with a specified 
#' width and space between them
#' 
#' @param x A numeric or unit object specifying the starting column x-position.
#' @param width A numeric or unit object specifying the width of columns.
#' @param space A numeric or unit object specifying the space between columns.
#' @param n  An integer specifying the number of elements to generate column 
#' positions for.
#' @param default.units A string indicating the default units to use
#' if \code{x}, \code{w}, or \code{s} are only given as numerics.
#' Default value is \code{default.units = "inches"}
#' 
#' @examples 
#' # Starting at 0.5 units, return a vector of positions for 3 objects that
#' # are 2 units in width with 0.1 units of space between them
#' 
#' pageLayoutCol(x = 0.5, width = 2, space = 0.1, n = 3, 
#'             default.units = "inches")
#' 
#' @returns Returns a unit vector of page positions.
#' @export
pageLayoutCol <- \(x, width, space, n, default.units = "inches") 

    ## Calculate positions
    pageLayoutRow(y = x, height = width, space = space, n = n, 
                  default.units = default.units)


