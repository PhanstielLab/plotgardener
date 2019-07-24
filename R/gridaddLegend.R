#' @export
#'

gridaddLegend <- function(color_vector, min_label, max_label, height, width, x, y, pageheight, orientation = "vertical", fontsize = 8, fontcolor = "grey", ...){

  # Convert user-inputted coordinates
  legend_coordinates <- convert_coordinates(height, width, x, y, pageheight)

  # Viewports
  if(is.null(current.vpPath()) == FALSE){
    upViewport()
  }

  vp <- viewport(height=unit(height,"in"),width=unit(width,"in"),x=unit(legend_coordinates[1],"in"),y=unit(legend_coordinates[2],"in"))
  pushViewport(vp)

  ncolors <- length(color_vector)


  if (orientation == "vertical"){

    xcoords <- rep(0.5, ncolors)
    ycoords <- seq(from = (1/ncolors)/2, to = 1, by = 1/ncolors)
    grid.rect(xcoords, ycoords, height = 1/ncolors, width = 1, gp=gpar(col=color_vector,fill=color_vector))

    grid.text(min_label, x = 0.5, y = -0.1, gp = gpar(fontsize = fontsize, col = fontcolor))
    grid.text(max_label, x = 0.5, y = 1.1, gp = gpar(fontsize = fontsize, col = fontcolor))
  }



  if (orientation == "horizontal"){

    xcoords <- seq(from = (1/ncolors)/2, to = 1, by = 1/ncolors)
    ycoords <- rep(0.5, ncolors)
    grid.rect(xcoords, ycoords, height = 1, width = 1/ncolors, gp=gpar(col=color_vector,fill=color_vector))
    grid.text(min_label, x = -0.1, y = 0.5, gp = gpar(fontsize = fontsize, col = fontcolor))
    grid.text(max_label, x = 1.1, y = 0.5, gp = gpar(fontsize = fontsize, col = fontcolor))


  }


}
