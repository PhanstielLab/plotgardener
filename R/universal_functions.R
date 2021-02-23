## Define a function to turn gTree child into a gPath
convert_gpath <- function(grob){

  ## Get the name of the grob
  name <- grob$name

  ## Turn it into a gPath
  gpath <- gPath(name)

  return(gpath)

}

## Define a function to grab the name of a viewport
viewport_name <- function(viewport){

  return(viewport$name)

}

## Define a function that gets the children of a group viewport
vp_children <- function(group_name){

  children <- unlist(current.vpTree()$children$bb_page$children[group_name], recursive = F)[[2]]
  return(children)

}

## Define a function to get a list of current viewports
current_viewports <- function(){

  if (!"bb_page" %in% names(lapply(current.vpTree()$children, viewport_name))){

    current <- as.list(names(lapply(current.vpTree()$children, viewport_name)))

  } else {

    ## Check for groups
    page_children <- names(lapply(current.vpTree()$children$bb_page$children, viewport_name))

    if (length(grep(pattern = "bb_group", x = page_children)) > 0){

      group_vps <- as.list(page_children[grep(pattern = "bb_group", x = page_children)])

      group_children <- unlist(lapply(group_vps, vp_children), recursive = F)

      children_vps <- lapply(group_children, viewport_name)

      current <- c(page_children, children_vps)

    } else {

      current <- as.list(page_children)

    }

  }

  return(current)
}

## Define a function to convert viewport x and y into center based on justification
adjust_vpCoords <- function(viewport){
  vp_y <- viewport$y

  if (length(viewport$justification == 2)){

    if ("left" %in% viewport$justification & "center" %in% viewport$justification){

      ## convert the x-coordinate only
      vp_x <- viewport$x + (0.5 * viewport$width)
    } else if ("right" %in% viewport$justification & "center" %in% viewport$justification){
      ## convert the x-coordinate only
      vp_x <- viewport$x - (0.5 * viewport$width)
    } else if ("center" %in% viewport$justification & "bottom" %in% viewport$justification){
      ## convert the y-coordinate only
      vp_x <- viewport$x
      vp_y <- vp_y + (0.5 * viewport$height)
    } else if ("center" %in% viewport$justification & "top" %in% viewport$justification){
      ## convert the y-coordinate only
      vp_x <- viewport$x
      vp_y <- vp_y - (0.5 * viewport$height)
    } else if ("left" %in% viewport$justification & "top" %in% viewport$justification){
      ## convert x-coordinate and y-coordinate
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- vp_y - (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "top" %in% viewport$justification){
      ## convert x-coordinate and y-coordinate
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- vp_y - (0.5 * viewport$height)
    } else if ("left" %in% viewport$justification & "bottom" %in% viewport$justification){
      ## convert x-coordinate and y-coordinate
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- vp_y + (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "bottom" %in% viewport$justification){
      ## convert x-coordinate and y-coordinate
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- vp_y + (0.5 * viewport$height)
    } else {
      ## no conversion
      vp_x <- viewport$x
    }

  } else if (length(viewport$justification == 1)){

    if (viewport$justification == "left"){
      ## convert the x-coordinate only
      vp_x <- viewport$x + (0.5 * viewport$width)
    } else if (viewport$justification == "right"){
      ## convert the x-coordinate only
      vp_x <- viewport$x - (0.5 * viewport$width)
    } else if (viewport$justification == "bottom"){
      ## convert the y-coordinate only
      vp_x <- viewport$x
      vp_y <- vp_y + (0.5 * viewport$height)
    } else if (viewport$justification == "top"){
      ## convert the y-coordinate only
      vp_x <- viewport$x
      vp_y <- vp_y - (0.5 * viewport$height)
    } else {
      ## no conversion
      vp_x <- viewport$x
    }

  }

  return(list(vp_x, vp_y))

}

## Define a function to change viewport x and y-coordinates to top left based on justification
vp_topLeft <- function(viewport){

  if (length(viewport$justification == 2)){

    if ("left" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y + (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y + (0.5 * viewport$height)
    } else if ("center" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y + (viewport$height)
    } else if ("center" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if ("left" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y
    } else if ("right" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y
    } else if ("left" %in% pviewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y + (viewport$height)
    } else if ("right" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y + (viewport$height)
    } else {
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y + (0.5 * viewport$height)
    }

  } else if (length(viewport$justification == 1)){

    if (viewport$justification == "left"){
      vp_x <- viewport$x
      vp_y <- viewport$y + (0.5 * viewport$height)
    } else if (viewport$justification == "right"){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y + (0.5 * viewport$height)
    } else if (viewport$justification == "bottom"){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y + (viewport$height)
    } else if (viewport$justification == "top"){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y
    } else {
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y + (0.5 * viewport$height)
    }

  }

  return(list(vp_x, vp_y))

}

## Define a function to change viewport x and y-coordinates to bottom left based on justification
vp_bottomLeft <- function(viewport){

  if (length(viewport$justification == 2)){

    if ("left" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if ("center" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if ("center" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else if ("left" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y - (viewport$height)
    } else if ("right" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else if ("left" %in% pviewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y
    } else if ("right" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y
    } else {
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    }

  } else if (length(viewport$justification == 1)){

    if (viewport$justification == "left"){
      vp_x <- viewport$x
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (viewport$justification == "right"){
      vp_x <- viewport$x - (viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (viewport$justification == "bottom"){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if (viewport$justification == "top"){
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else {
      vp_x <- viewport$x - (0.5 * viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    }

  }

  return(list(vp_x, vp_y))

}

## Define a function to change viewport x and y-coordinates to bottom right based on justification
vp_bottomRight <- function(viewport){

  if (length(viewport$justification == 2)){

    if ("left" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x + viewport$width
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if ("right" %in% viewport$justification & "center" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if ("center" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if ("center" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else if ("left" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x + viewport$width
      vp_y <- viewport$y - (viewport$height)
    } else if ("right" %in% viewport$justification & "top" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y - (viewport$height)
    } else if ("left" %in% pviewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x + viewport$width
      vp_y <- viewport$y
    } else if ("right" %in% viewport$justification & "bottom" %in% viewport$justification){
      vp_x <- viewport$x
      vp_y <- viewport$y
    } else {
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    }

  } else if (length(viewport$justification == 1)){

    if (viewport$justification == "left"){
      vp_x <- viewport$x + viewport$width
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (viewport$justification == "right"){
      vp_x <- viewport$x
      vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (viewport$justification == "bottom"){
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y
    } else if (viewport$justification == "top"){
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y - (viewport$height)
    } else {
      vp_x <- viewport$x + (0.5 * viewport$width)
      vp_y <- viewport$y - (0.5 * viewport$height)
    }

  }

  return(list(vp_x, vp_y))

}

## Define a function to convert to page units
convert_page <- function(object){

  ## Get page_height and its units from bbEnv through bb_makepage
  page_height <- get("page_height", envir = bbEnv)
  page_units <- get("page_units", envir = bbEnv)

  ## Convert x and y coordinates and height and width to same page_units
  old_x <- object$x
  old_y <- object$y
  old_height <- object$height
  old_width <- object$width
  new_x <- convertX(old_x, unitTo = page_units)
  new_y <- convertY(unit(page_height, units = page_units) - convertY(old_y, unitTo = page_units), unitTo = page_units)
  new_height <- convertHeight(old_height, unitTo = page_units)
  new_width <- convertWidth(old_width, unitTo = page_units)

  object$x <- new_x
  object$y <- new_y
  object$height <- new_height
  object$width <- new_width


  return(object)

}

## Define a function to make sure a bb_page viewport exists
check_bbpage <- function(error){

  if (!"bb_page" %in% current.vpPath()){

    stop(error, call. = FALSE)

  }

}

## Define a function to check dimensions/placing coordinates
check_placement <- function(object){

  if (attributes(object)$plotted == T){

    ## If giving placement coordinates
    if (!is.null(object$x) | !is.null(object$y)){

      ## 1. Need both an x and y coordinate
      if (!is.null(object$x) & is.null(object$y)){

        stop("Placement detected with y value missing.", call. = FALSE)

      }

      if (!is.null(object$y) & is.null(object$x)){

        stop("Placement detected with x value missing.", call. = FALSE)

      }

      ## 2. Need plot dimensions
      if (is.null(object$width)){

        stop("Placement detected with plot width missing.", call. = FALSE)

      }

      if (is.null(object$height)){

        stop("Placement detected with plot height missing.", call. = FALSE)

      }

      ## 3. Need a bb_page
      check_bbpage(error = "Must make a BentoBox page with bb_makePage() before placing a plot.")

    }

  }

}

## Define a character vector of valid coordinate systems to work in
validUnits <- c("npc", "native", "inches", "cm", "mm", "points", "bigpts", "picas", "dida",
                "cicero", "scaledpts", "char", "lines", "snpc")

## Define a function that converts coordinates/dimensions into default units
defaultUnits <- function(object, default.units){

  if (!(is.null(object$x) & is.null(object$y))){

    if (!"unit" %in% class(object$x)){

      if (!is.numeric(object$x)){

        stop("x-coordinate is neither a unit object nor a numeric value. Cannot place object.", call. = FALSE)

      }

      if (is.null(default.units)){

        stop("x-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      object$x <- unit(object$x, default.units)

    }


    if (!"unit" %in% class(object$y)){

      ## Check for "below" y-coord
      if (grepl("b", object$y) == TRUE){
        if (grepl("^[ac-zA-Z]+$", object$y) == TRUE){
          stop("\'below\' y-coordinate detected with additional letters. Cannot parse y-coordinate.", call. = FALSE)
        }

        if(is.na(as.numeric(gsub("b","", object$y)))){
          stop("\'below\' y-coordinate does not have a numeric associated with it. Cannot parse y-coordinate.", call. = FALSE)
        }

       object$y <- plot_belowY(y_coord = object$y)

      } else {

        if (!is.numeric(object$y)){

          stop("y-coordinate is neither a unit object nor a numeric value. Cannot place object.", call. = FALSE)

        }

        if (is.null(default.units)){

          stop("y-coordinate detected as numeric.\'default.units\' must be specified.", call. = FALSE)

        }

        object$y <- unit(object$y, default.units)

      }


    }

    if (!"unit" %in% class(object$width)){

      if (!is.numeric(object$width)){

        stop("Width is neither a unit object nor a numeric value. Cannot place object.", call. = FALSE)

      }

      if (is.null(default.units)){

        stop("Width detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      object$width <- unit(object$width, default.units)

    }

    if (!"unit" %in% class(object$height)){

      if (!is.numeric(object$height)){

        stop("Height is neither a unit object nor a numeric value. Cannot place object.", call. = FALSE)

      }

      if (is.null(default.units)){

        stop("Height detected as numeric.\'default.units\' must be specified.", call. = FALSE)

      }

      object$height <- unit(object$height, default.units)

    }

  }

  return(object)

}

## Define a function that determines the chromosome offsets for a plot with multiple chromosomes (manhattan plots)
spaceChroms <- function(assemblyData, space){

  ## space is the space in between each chromosome as a fraction of the width of the plot

  ## Determine the offset for each chromomse
  cumsums <- cumsum(as.numeric(assemblyData[,2]))
  spacer <- cumsums[length(cumsums)] * space
  additionalSpace <- (1:length(cumsums)-0) * spacer

  ## Start position
  startPos <- c(0, cumsums[1:length(cumsums)-1])
  startPos <- startPos + additionalSpace
  assemblyData[,3] <- startPos

  ## Stop Position
  stopPos <- cumsums + (1:(length(cumsums))) * spacer
  assemblyData[,4] <- stopPos

  colnames(assemblyData) <- c("chrom", "length", "start", "stop")

  return(assemblyData)

}

## Define a function that makes a color transparent
makeTransparent <- function(color, alpha){

  if (is.null(alpha)){
    alpha <- 1
  }

  rgb <- col2rgb(color)
  transp <- rgb(rgb[1], rgb[2], rgb[3], alpha = alpha*255, maxColorValue = 255)
  return(transp)

}

# Define a function that parses parameters from a bb_params object
parseParams <- function(bb_params, object_params){

  getNull <- function(name, params){

    val <- params[[name]]
    if (is.null(val)){
      return(name)
    }

  }

  if (!is.null(bb_params)){
    ## First make sure it's actually a "bb_params" object
    if (class(bb_params) != "bb_params"){
      warning("Input object ignored. Object must be a \'bb_params\' class object.", call. = FALSE)
    } else {

      objectNames <- names(object_params)

      ## Get the NULL object params that haven't been set with other parameter input
      nullNames <- unlist(lapply(objectNames, getNull, params = object_params))
      ## Get any of the object values that match those names
      objectMatches <- bb_params[which(names(bb_params) %in% nullNames)]
      ## Replace object_param values with bb_params values
      object_params[match(names(objectMatches), names(object_params))] <- objectMatches


      addParams <- bb_params[which(!names(bb_params) %in% objectNames)]

      object_params <- c(object_params, addParams)

    }

  }

  return(object_params)
}

## Define a function that parses gpar parameters
setGP <- function(gpList, params, ...){

  availGPs <- names(get.gpar())
  gpMatches <- params[which(names(params) %in% availGPs)]
  gpList[names(gpMatches)] <- gpMatches
  gpList[names(list(...))] <- list(...)

  ## Reset with fontface first
  if ("fontface" %in% names(gpList)){
    otherParams <- gpList
    gpList <- gpar(fontface = gpList$fontface)
    gpList[names(otherParams)] <- otherParams
  }

  return(gpList)
}


## Define a function to get the data packages of a genome assembly
getPackages <- function(genome, TxDb=NULL){
  switch(genome,
         "bosTau8" = structure(list("Genome" = "bosTau8", "TxDb" = ifelse(is.null(TxDb), "TxDb.Btaurus.UCSC.bosTau8.refGene", TxDb),
                          "OrgDb" = "org.Bt.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Btaurus.UCSC.bosTau8"), class = "bb_assembly"),
         "bosTau9" = structure(list("Genome" = "bosTau9", "TxDb" = ifelse(is.null(TxDb), "TxDb.Btaurus.UCSC.bosTau9.refGene", TxDb),
                          "OrgDb" = "org.Bt.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Btaurus.UCSC.bosTau9"), class = "bb_assembly"),
         "canFam3" = structure(list("Genome" = "canFam3","TxDb" = ifelse(is.null(TxDb), "TxDb.Cfamiliaris.UCSC.canFam3.refGene", TxDb),
                          "OrgDb" = "org.Cf.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Cfamiliaris.UCSC.canFam3"), class = "bb_assembly"),
         "ce6" = structure(list("Genome" = "ce6","TxDb" = ifelse(is.null(TxDb), "TxDb.Celegans.UCSC.ce6.ensGene", TxDb),
                      "OrgDb" = "org.Ce.eg.db", gene.id.column = "GENEID", display.column = "GENEID",
                      "BSgenome" = "BSgenome.Celegans.UCSC.ce6"), class = "bb_assembly"),
         "ce11" = structure(list("Genome" = "ce11","TxDb" = ifelse(is.null(TxDb), "TxDb.Celegans.UCSC.ce11.refGene", TxDb),  # TxDb.Celegans.UCSC.ce11.ensGene also available
                       "OrgDb" = "org.Ce.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                       "BSgenome" = "BSgenome.Celegans.UCSC.ce11"), class = "bb_assembly"),
         "danRer10" = structure(list("Genome" = "danRer10","TxDb" = ifelse(is.null(TxDb), "TxDb.Drerio.UCSC.danRer10.refGene", TxDb),
                           "OrgDb" = "org.Dr.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                           "BSgenome" = "BSgenome.Drerio.UCSC.danRer10"), class = "bb_assembly"),
         "danRer11" = structure(list("Genome" = "danRer11","TxDb" = ifelse(is.null(TxDb), "TxDb.Drerio.UCSC.danRer11.refGene", TxDb),
                           "OrgDb" = "org.Dr.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                           "BSgenome" = "BSgenome.Drerio.UCSC.danRer11"), class = "bb_assembly"),
         "dm3" = structure(list("Genome" = "dm3","TxDb" = ifelse(is.null(TxDb), "TxDb.Dmelanogaster.UCSC.dm3.ensGene", TxDb),
                      "OrgDb" = "org.Dm.eg.db", gene.id.column = "ENSEMBL", display.column = "SYMBOL",
                      "BSgenome" = "BSgenome.Dmelanogaster.UCSC.dm3"), class = "bb_assembly"),
         "dm6" = structure(list("Genome" = "dm6","TxDb" = ifelse(is.null(TxDb), "TxDb.Dmelanogaster.UCSC.dm6.ensGene", TxDb),
                      "OrgDb" = "org.Dm.eg.db", gene.id.column = "ENSEMBL", display.column = "SYMBOL",
                      "BSgenome" = "BSgenome.Dmelanogaster.UCSC.dm6"), class = "bb_assembly"),
         "galGal4" = structure(list("Genome" = "galGal4","TxDb" = ifelse(is.null(TxDb), "TxDb.Ggallus.UCSC.galGal4.refGene", TxDb),
                          "OrgDb" = "org.Gg.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Ggallus.UCSC.galGal4"), class = "bb_assembly"),
         "galGal5" = structure(list("Genome" = "galGal5", "TxDb" = ifelse(is.null(TxDb), "TxDb.Ggallus.UCSC.galGal5.refGene", TxDb),
                          "OrgDb" = "org.Gg.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Ggallus.UCSC.galGal5"), class = "bb_assembly"),
         "galGal6" = structure(list("Genome" = "galGal6","TxDb" = ifelse(is.null(TxDb), "TxDb.Ggallus.UCSC.galGal6.refGene", TxDb),
                          "OrgDb" = "org.Gg.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Ggallus.UCSC.galGal6"), class = "bb_assembly"),
         "hg18" = structure(list("Genome" = "hg18","TxDb" = ifelse(is.null(TxDb), "TxDb.Hsapiens.UCSC.hg18.knownGene", TxDb),
                       "OrgDb" = "org.Hs.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                       "BSgenome" = "BSgenome.Hsapiens.UCSC.hg18"), class = "bb_assembly"),
         "hg19" = structure(list("Genome" = "hg19","TxDb" = ifelse(is.null(TxDb), "TxDb.Hsapiens.UCSC.hg19.knownGene", TxDb),
                       "OrgDb" = "org.Hs.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                       "BSgenome" = "BSgenome.Hsapiens.UCSC.hg19"), class = "bb_assembly"),
         "hg38" = structure(list("Genome" = "hg38","TxDb" = ifelse(is.null(TxDb), "TxDb.Hsapiens.UCSC.hg38.knownGene", TxDb),
                       "OrgDb" = "org.Hs.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                       "BSgenome" = "BSgenome.Hsapiens.UCSC.hg38"), class = "bb_assembly"),
         "mm9" = structure(list("Genome" = "mm9","TxDb" = ifelse(is.null(TxDb), "TxDb.Mmusculus.UCSC.mm9.knownGene", TxDb),
                      "OrgDb" = "org.Mm.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                      "BSgenome" = "BSgenome.Mmusculus.UCSC.mm9"), class = "bb_assembly"),
         "mm10" = structure(list("Genome" = "mm10","TxDb" = ifelse(is.null(TxDb), "TxDb.Mmusculus.UCSC.mm10.knownGene", TxDb),
                       "OrgDb" = "org.Mm.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                       "BSgenome" = "BSgenome.Mmusculus.UCSC.mm10"), class = "bb_assembly"),
         "rheMac3" = structure(list("Genome" = "rheMac3","TxDb" = ifelse(is.null(TxDb), "TxDb.Mmulatta.UCSC.rheMac3.refGene", TxDb),
                          "OrgDb" = "org.Mmu.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Mmulatta.UCSC.rheMac3"), class = "bb_assembly"),
         "rheMac8" = structure(list("Genome" = "rheMac8","TxDb" = ifelse(is.null(TxDb), "TxDb.Mmulatta.UCSC.rheMac8.refGene", TxDb),
                          "OrgDb" = "org.Mmu.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Mmulatta.UCSC.rheMac8"), class = "bb_assembly"),
         "rheMac10" = structure(list("Genome" = "rheMac10","TxDb" = ifelse(is.null(TxDb), "TxDb.Mmulatta.UCSC.rheMac10.refGene", TxDb),
                           "OrgDb" = "org.Mmu.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                           "BSgenome" = "BSgenome.Mmulatta.UCSC.rheMac10"), class = "bb_assembly"),
         "panTro5" = structure(list("Genome" = "panTro5","TxDb" = ifelse(is.null(TxDb), "TxDb.Ptroglodytes.UCSC.panTro5.refGene", TxDb),
                          "OrgDb" = "org.Pt.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Ptroglodytes.UCSC.panTro5"), class = "bb_assembly"),
         "panTro6" = structure(list("Genome" = "panTro6","TxDb" = ifelse(is.null(TxDb), "TxDb.Ptroglodytes.UCSC.panTro6.refGene", TxDb),
                          "OrgDb" = "org.Pt.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Ptroglodytes.UCSC.panTro6"), class = "bb_assembly"),
         "rn4" = structure(list("Genome" = "rn4","TxDb" = ifelse(is.null(TxDb), "TxDb.Rnorvegicus.UCSC.rn4.ensGene", TxDb),
                      "OrgDb" = "org.Rn.eg.db", gene.id.column = "ENSEMBL", display.column = "SYMBOL",
                      "BSgenome" = "BSgenome.Rnorvegicus.UCSC.rn4"), class = "bb_assembly"),
         "rn5" = structure(list("Genome" = "rn5","TxDb" = ifelse(is.null(TxDb), "TxDb.Rnorvegicus.UCSC.rn5.refGene", TxDb),
                      "OrgDb" = "org.Rn.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                      "BSgenome" = "BSgenome.Rnorvegicus.UCSC.rn5"), class = "bb_assembly"),
         "rn6" = structure(list("Genome" = "rn6","TxDb" = ifelse(is.null(TxDb), "TxDb.Rnorvegicus.UCSC.rn6.refGene", TxDb),
                      "OrgDb" = "org.Rn.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                      "BSgenome" = "BSgenome.Rnorvegicus.UCSC.rn6"), class = "bb_assembly"),
         "sacCer2" = structure(list("Genome" = "sacCer2","TxDb" = ifelse(is.null(TxDb), "TxDb.Scerevisiae.UCSC.sacCer2.sgdGene", TxDb),
                          "OrgDb" = "org.Sc.sgd.db", gene.id.column = "ORF", display.column = "GENENAME",
                          "BSgenome" = "BSgenome.Scerevisiae.UCSC.sacCer2"), class = "bb_assembly"),
         "sacCer3" = structure(list("Genome" = "sacCer3","TxDb" = ifelse(is.null(TxDb), "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", TxDb),
                          "OrgDb" = "org.Sc.sgd.db", gene.id.column = "ORF", display.column = "GENENAME",
                          "BSgenome" = "BSgenome.Scerevisiae.UCSC.sacCer3"), class = "bb_assembly"),
         "susScr3" = structure(list("Genome" = "susScr3","TxDb" = ifelse(is.null(TxDb), "TxDb.Sscrofa.UCSC.susScr3.refGene", TxDb),
                          "OrgDb" = "org.Ss.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                          "BSgenome" = "BSgenome.Sscrofa.UCSC.susScr3"), class = "bb_assembly"),
         "susScr11" = structure(list("Genome" = "susScr11","TxDb" = ifelse(is.null(TxDb), "TxDb.Sscrofa.UCSC.susScr11.refGene", TxDb),
                           "OrgDb" = "org.Ss.eg.db", gene.id.column = "ENTREZID", display.column = "SYMBOL",
                           "BSgenome" = "BSgenome.Sscrofa.UCSC.susScr11"), class = "bb_assembly"),
  )
}

## Define a function to get the assembly info based on a string (ie default) or bb_assembly object
parse_bbAssembly <- function(assembly){

  availDefaults <- c("bosTau8", "bosTau9", "canFam3", "ce6", "ce11", "danRer10",
                     "danRer11", "dm3", "dm6", "galGal4", "galGal5", "galGal6",
                     "hg18", "hg19", "hg38", "mm9", "mm10", "rheMac3", "rheMac8",
                     "rehMac10", "panTro5", "panTro6", "rn4", "rn5", "rn6", "sacCer2",
                     "sacCer3", "susScr3", "susScr11")

  ## If it's just a string, get the default
  if (class(assembly) == "character"){

    if (!assembly %in% availDefaults){
      stop("\'assembly\' not available as a default. Please make a bb_assembly object with `bb_assembly()` or pick an assembly from the defaults listed with `bb_genomes()`.", call. = FALSE)
    }

    assemblyData <- getPackages(genome = assembly)

    ## If it's a bb_assembly object, use those
  } else if (class(assembly) == "bb_assembly"){

    assemblyData <- assembly

  } else {
    stop("Invalid \'assembly\' type. Please make a bb_assembly object with `bb_assembly()` or input an assembly string from the defaults listed with `bb_genomes()`.", call. = FALSE)
  }


  return(assemblyData)


}

## Define a function to check that a TxDb, Org, or BSgenome package is loaded
check_loadedPackage <- function(package, message){
  if (!package %in% (.packages())){
    warning(message, call. = FALSE)
    return(FALSE)
  } else {
    return(TRUE)
  }

}

## Define a function that will by default prioritize genes by citations (if available) or gene lengths
defaultGenePriorities <- function(data, assembly, transcript = FALSE){

  ## Define our list of available defaults that have citations
  availCitations<- list(bosTau8 = "Citations.Btaurus.NCBI.bosTau8", bosTau9 = "Citations.Btaurus.NCBI.bosTau9",
                        canFam3 = "Citations.Cfamiliaris.NCBI.canFam3", ce11 = "Citations.Celegans.NCBI.ce11",
                        danRer10 = "Citations.Drerio.NCBI.danRer10", danRer11 = "Citations.Drerio.NCBI.danRer11",
                        dm3 = "Citations.Dmelanogaster.NCBI.dm3", dm6 = "Citations.Dmelanogaster.NCBI.dm6",
                        galGal4 = "Citations.Ggallus.NCBI.galGal4", galGal5 = "Citations.Ggallus.NCBI.galGal5", galGal6 = "Citations.Ggallus.NCBI.galGal6",
                        hg18 = "Citations.Hsapiens.NCBI.hg18", hg19 = "Citations.Hsapiens.NCBI.hg19", hg38 = "Citations.Hsapiens.NCBI.hg38",
                        mm9 = "Citations.Mmusculus.NCBI.mm9", mm10 = "Citations.Mmusculus.NCBI.mm10",
                        rheMac3 = "Citations.Mmulatta.NCBI.rheMac3", rheMac8 = "Citations.Mmulatta.NCBI.rheMac8", rehMac10 = "Citations.Mmulatta.NCBI.rheMac10",
                        panTro5 = "Citations.Ptroglodytes.NCBI.panTro5", panTro6 = "Citations.Ptroglodytes.NCBI.panTro6",
                        rn4 = "Citations.Rnorvegicus.NCBI.rn4", rn5 = "Citations.Rnorvegicus.NCBI.rn5", rn6 = "Citations.Rnorvegicus.NCBI.rn6",
                        sacCer2 = "Citations.Scerevisiae.NCBI.sacCer2", sacCer3 = "Citations.Scerevisiae.NCBI.sacCer3",
                        susScr3 = "Citations.Sscrofa.NCBI.susScr3", susScr11 = "Citations.Sscrofa.NCBI.susScr11")

  ## Define assemblies whose TxDb IDs will need to be converted to ENTREZID from a different ID
  convertIDs <- list(dm3 = "ENSEMBL", dm6 = "ENSEMBL",
                     rn4 = "ENSEMBL", sacCer2 = "ORF",
                     sacCer3 = "ORF")

  ## Split TxDb assembly name
  assemblyName <- unlist(strsplit(assembly$TxDb, split = "[.]"))

  ## If assembly is included in package, access citations
  if (any(names(availCitations) %in% assemblyName)){

    name <- names(availCitations)[which(names(availCitations) %in% assemblyName)]

    ## Convert necessary builds to ENTREZID
    if (name %in% names(convertIDs)){

      ## Load associated OrgDb
      org_db <- eval(parse(text = assembly$OrgDb))

      ## Convert gene ids in data to ENTREZID based on previous keytype
      entrezIDs  <- suppressMessages(select(org_db, keys = data$GENEID, columns = "ENTREZID", keytype = convertIDs[[name]]))

      data$ENTREZID <- entrezIDs

    } else {
      data$ENTREZID <- data$GENEID
    }

    ## Get internal citation data and match based on ENTREZID
    citationData <- eval(parse(text = availCitations[[name]]))
    updatedData <- suppressMessages(dplyr::left_join(x = data, y = citationData, by = "ENTREZID"))

    ## Set any missing citations to 0
    updatedData[is.na(updatedData$Citations),]$Citations <- rep(0, nrow(updatedData[is.na(updatedData$Citations),]))

    if (transcript == TRUE){
      updatedData <- updatedData[duplicated(updatedData$TXNAME) == F,]
    }


    ## Order based on citation number
    updatedData <- updatedData[order(updatedData$Citations, decreasing = TRUE),]

  } else {

    ## With no internal citation data, set default priority based on gene/transcript length
    #data$GeneLength <- data$TXEND - data$TXSTART
    updatedData <- data[order(data$length, decreasing = TRUE),]

  }

  ## Return data with priorities
  return(updatedData)

}

## Parse y-coordinates for plots relative to bottom of last plot
plot_belowY <- function(y_coord){

  currentViewports <- current_viewports()
  if (length(currentViewports) == 0){
    stop("No previous plot detected. Cannot define a \'below\' y-coordinate.", call. = FALSE)
  }

  lastPlot <- currentViewports[length(currentViewports)]
  ## Go into viewport to get coordinates/dimensions
  seekViewport(name = lastPlot)
  x <- current.viewport()$x
  y <- current.viewport()$y
  width <- current.viewport()$width
  height <- current.viewport()$height
  just <- current.viewport()$just

  ## Exit viewport
  upViewport()

  ## Convert to bottom
  bottomCoords <- vp_bottomLeft(viewport(x = x, y = y, width = width, height = height, just = just))

  ## Parse number and convert to page units
  space <- as.numeric(gsub("b", "", y_coord))
  space <- unit(space, get("page_units", envir = bbEnv))

  new_y <- (unit(get("page_height", envir = bbEnv), get("page_units", envir = bbEnv)) - bottomCoords[[2]]) + space
  new_y <- convertY(new_y, unitTo = get("page_units", envir = bbEnv))

  return(new_y)
}
