hideGuides <- function(envir=bbEnv){
  if(length(envir$guides) > 0){
    for(i in 1:length(envir$guides)) grid.remove(gPath(envir$guides[[i]]$name)$name)
    bbEnv$guides <- list()
  }
}
