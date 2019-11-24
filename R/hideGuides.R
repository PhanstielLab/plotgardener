hideGuides <- function(){

  # if(length(envir$guides) > 0){
  #   for(i in 1:length(envir$guides)) grid.remove(gPath(envir$guides[[i]]$name)$name)
  #   bbEnv$guides <- list()
  # }

  guides <- get("guides", envir = bbEnv)

  if (length(guides$children) > 0){

    ## convert the grobs into gpaths
    grob_list <- lapply(guides$children, convert_gpath)

    ## get the last grob
    last_grob <- grob_list[length(grob_list)]

    ## remove the last grob from the list of grobs
    grob_list <- grob_list[- length(grob_list)]

    ## remove all grobs except last one, not redrawing each time (to save time)
    invisible(lapply(grob_list, grid.remove, redraw = FALSE))

    ## remove last grob with redrawing, now removing all the grobs
    invisible(lapply(last_grob, grid.remove))

  }

}
