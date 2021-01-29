# Maps numeric vector to color palette
# vec: numeric vector to map to color
# col: color palette to map to
# num: number of bins of colors
# range: range of values to map
bb_maptocolors <- function(vec, col, num = 100, range = NULL){

  if (is.null(range) == TRUE){
    breaks <- seq(min(vec), max(vec), length.out = num)
  } else {
    vec[which(vec < range[1])] = range[1]
    vec[which(vec > range[2])] = range[2]
    breaks <- seq(range[1], range[2], length.out = num)
  }

  cols <- col(length(breaks) + 1)
  colvec <- as.character(cut(vec, c(-Inf, breaks, Inf), labels = cols))
  return(colvec)


}
