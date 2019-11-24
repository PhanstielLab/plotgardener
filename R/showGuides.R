showGuides <- function(){
  #for(i in 1:length(envir$guides)) grid.draw(envir$guides[[i]])

  guides <- get("guides", envir = bbEnv)

  invisible(lapply(guides$children, grid.draw))



}
