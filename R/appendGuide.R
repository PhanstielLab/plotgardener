appendGuide <- function(x){
  #envir$guides[[length(envir$guides)+1]] <- x
  assign("guides", addGrob(gTree = get("guides", envir = bbEnv), child = x), envir = bbEnv)
  grid.draw(x)
}
