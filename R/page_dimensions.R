bbEnv <- new.env()
if("pageDim" %in% ls(envir = .GlobalEnv)){
  assign("height", unname(pageDim("height")), envir = bbEnv)
  assign("width", unname(pageDim("width")), envir = bbEnv)
  assign("top_margin", unname(pageDim("top_margin")), envir = bbEnv)
  assign("bottom_margin", unname(pageDim("bottom_margin")), envir = bbEnv)
  assign("left_margin", unname(pageDim("left_margin")), envir = bbEnv)
  assign("right_margin", unname(pageDim("right_margin")), envir = bbEnv)
} else {
  assign("height", 11, envir = bbEnv)
  assign("width", 8.5, envir = bbEnv)
  assign("top_margin", 1, envir = bbEnv)
  assign("bottom_margin", 1, envir = bbEnv)
  assign("left_margin", 1, envir = bbEnv)
  assign("right_margin", 1, envir = bbEnv)

}





