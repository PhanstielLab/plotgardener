showGuides <- function(envir=bbEnv){
  for(i in 1:length(envir$guides)) grid.draw(envir$guides[[i]])
}
