appendGuide <- function(x, envir=bbEnv){
  envir$guides[[length(envir$guides)+1]] <- x
  grid.draw(x)
}
