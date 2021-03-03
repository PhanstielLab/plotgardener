#' Plot a Manhattan plot
#'
#' @param data Data to be plotted, as a character value specifying a file path of GWAS data, a dataframe, or a \link[GenomicRanges]{GRanges} object. Each of these data types must have the following columns:
#' \itemize{
#' \item{\code{"chr"}: }{Chromosome names. This column must be a character.}
#' \item{\code{"pos"}: }{Chromosomal position. This column must be an integer or numeric.}
#' \item{\code{"p"}: }{p-value. This column must be numeric. p-values will be converted to -log(10) space.}
#' \item{\code{"snp"}(optional): }{SNP name or rsid. This column should be a character.}
#' }
#' @param sigVal A numeric specifying the significance level of p-values. Along with data p-values, this value will be converted to -log1=(10) space.
#' Default value is \code{sigVal = 5e-08}.
#' @param chrom Chromosome of region to be plotted, as a string. If left \code{NULL}, all chromosomes found in data will be plotted.
#' @param chromstart Integer start position on chromosome to be plotted.
#' @param chromend Integer end position on chromosome to be plotted.
#' @param assembly Default genome assembly as a string or a \link[BentoBox]{bb_assembly} object. Default value is \code{assembly = "hg19"}.
#' @param fill Character value(s) as a single value, vector, or palette specifying fill colors of data points.
#' If \code{scaleLD} is supplied, colors will be mapped to \code{scaleLD} values. If \code{scaleLD} is not supplied, color vectors and palettes will
#' only be mapped to different chromosomes of a multi-chromosomal plot. Default value is \code{fill = "black"}.
#' @param pch A numeric value or numeric vector specifying point symbols.
#' If \code{scaleLD} is supplied, point symbols will be mapped to \code{scaleLD} values. Default value is \code{pch = 19}.
#' @param cex A numeric indicating the amount by which points should be scaled relative to the default. Default value is \code{cex = 0.25}.
#' @param leadSNP A list specifying the lead SNP in the desired region and any associated aesthetic features of the lead SNP data point and text label.
#' The lead SNP should be specified as a character with the name slot \code{"snp"} in the list. Accepted lead SNP aesthetic features in the list include
#' \code{fill}, \code{pch}, \code{cex}, \code{fontcolor}, and \code{fontsize}.
#' @param scaleLD A character value specifying the \code{data} column name of linkage disequilibrium (LD) scores to apply \code{fill} and/or
#' \code{pch} vectors or functions to. LD scores will be grouped into the following ranges: 0-0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.
#' @param sigLine Logical value indicating whether to draw a line at the significance level indicated with \code{sigVal}. Default value is \code{sigLine = FALSE}.
#' @param sigCol Single character value specifying the color of significant data points. If \code{scaleLD} is supplied,
#' \code{sigCol} will be ignored.
#' @param ymax A numeric specifying the fraction of the max y-value to set as the height of the plot. Default value is \code{ymax = 1}.
#' @param range A numeric vector of length 2 specifying the y-range of p-values to plot (c(min, max)).
#' @param space A numeric value indicating the space between each chromsome as a fraction of the width of the plot, if plotting multiple chromosomes. Default value is \code{space = 0.01}.
#' @param bg Character value indicating background color. Default value is \code{bg = NA}.
#' @param baseline Logical value indicating whether to include a baseline along the x-axis. Default value is \code{baseline = FALSE}.
#' @param x A numeric or unit object specifying Manhattan plot x-location.
#' @param y A numeric, unit object, or character containing a "b" combined with a numeric value specifying Manhattan plot y-location. The character value will
#' place the Manhattan plot y relative to the bottom of the most recently plotted BentoBox plot according to the units of the BentoBox page.
#' @param width A numeric or unit object specifying Manhattan plot width.
#' @param height A numeric or unit object specifying Manhattan plot height.
#' @param just Justification of Manhattan plot relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{x}, \code{y}, \code{width}, or \code{height} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param draw A logical value indicating whether graphics output should be produced. Default value is \code{draw = TRUE}.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_manhattan} object containing relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load genomic assembly information
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' ## Load GWAS data
#' data("bb_gwasData")
#'
#' ## Create a page
#' bb_pageCreate(width = 7.5, height = 4.5, default.units = "inches")
#'
#' ## Plot all GWAS data
#' manhattanPlot <- bb_plotManhattan(data = bb_gwasData, assembly = "hg19",
#'                                   fill = c("grey", "#37a7db"), sigLine = TRUE,
#'                                   col = "grey", lty = 2, range = c(0, 14),
#'                                   x = 0.5, y = 0, width = 6.5, height = 2,
#'                                   just = c("left", "top"), default.units = "inches")
#' ## Annotate genome label
#' bb_annoGenomeLabel(plot = manhattanPlot, x = 0.5, y = 2, fontsize = 8, just = c("left", "top"),
#'                    default.units = "inches" )
#' bb_plotText(label = "Chromosome", fontsize = 8,
#'             x = 3.75, y = 2.20, just = "center", default.units = "inches")
#'
#' ## Annotate y-axis
#' bb_annoYaxis(plot = manhattanPlot, at = c(0, 2, 4, 6, 8, 10, 12, 14),
#'              axisLine = TRUE, fontsize = 8)
#'
#' ## Plot y-axis label
#' bb_plotText(label = "-log10(p-value)", x = 0.15, y = 1, rot = 90,
#'             fontsize = 8, fontface = "bold", just = "center", default.units = "inches")
#'
#'
#' ## Plot GWAS data zooming in on chromosome 11, highlighting a lead SNP, and coloring by LD score
#' leadSNP_p <- min(bb_gwasData[which(bb_gwasData$chr == "chr11"),]$p)
#' leadSNP <- bb_gwasData[which(bb_gwasData$p == leadSNP_p),]$snp
#' chr11_manhattanPlot <- bb_plotManhattan(data = bb_gwasData, chrom = "chr11",
#'                                         chromstart = 60000000, chromend = 130000000,
#'                                         fill = c("#1f4297", "#37a7db", "green", "orange", "red"),
#'                                         sigLine = TRUE, col = "grey", lty = 2, range = c(0, 16),
#'                                         leadSNP = list("snp" = leadSNP, pch = 18, cex = 0.75, fill = "#7ecdbb", fontsize = 8),
#'                                         scaleLD = "LD",
#'                                         x = 0.5, y = 2.5, width = 6.5, height = 1.5,
#'                                         just = c("left", "top"), default.units = "inches")
#'
#' ## Plot legend for LD scores
#' bb_plotLegend(legend = c("LD Ref Var",
#'                          expression(paste("0.4", ">", "r"^{paste("2")}, "", ">=", "0.2")),
#'                          expression(paste("0.2", ">", "r"^{paste("2")}, "", ">=", "0")),
#'                          "no LD data"),
#'               fill = c("#7ecdbb", "#37a7db", "#1f4297", "grey"), cex = 0.75,
#'               pch = c(18, 19, 19, 19), border = FALSE, x = 7, y = 2.5,
#'               width = 1.5, height = 0.6, just = c("right", "top"), default.units = "inches")
#'
#'
#' ## Annotate genome label
#' bb_annoGenomeLabel(plot = chr11_manhattanPlot, x = 0.5, y = 4.01, fontsize = 8, scale = "Mb",
#'                    just = c("left", "top"), default.units = "inches" )
#'
#' ## Annotate y-axis
#' bb_annoYaxis(plot = chr11_manhattanPlot, at = c(0, 2, 4, 6, 8, 10, 12, 14, 16),
#'              axisLine = TRUE, fontsize = 8)
#'
#' ## Plot y-axis label
#' bb_plotText(label = "-log10(p-value)", x = 0.15, y = 3.25, rot = 90,
#'             fontsize = 8, fontface = "bold", just = "center",
#'             default.units = "inches")
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @details
#' A Manhattan plot can be placed on a BentoBox coordinate page by providing plot placement parameters:
#' \preformatted{
#' bb_plotManhattan(data,
#'                  chrom = NULL,
#'                  chromstart = NULL, chromend = NULL,
#'                  x, y, width, height, just = c("left", "top"),
#'                  default.units = "inches")
#' }
#' This function can also be used to quickly plot an unannotated Manhattan plot by ignoring plot placement parameters:
#' \preformatted{
#' bb_plotManhattan(data,
#'                  chrom = NULL,
#'                  chromstart = NULL, chromend = NULL)
#' }
#'
#' @export
bb_plotManhattan <- function(data, sigVal = 5e-08, chrom = NULL, chromstart = NULL, chromend = NULL, assembly = "hg19", fill = "black", pch = 19, cex = 0.25,
                             leadSNP = NULL, scaleLD = NULL, sigLine = FALSE, sigCol = NULL, ymax = 1, range = NULL, space = 0.01, bg = NA, baseline = FALSE,
                             x = NULL, y = NULL, width = NULL, height = NULL, just = c("left", "top"), default.units = "inches", draw = TRUE, params = NULL, ...){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  ## Define a function that checks for errors in bb_plotManhattan
  errorcheck_bb_plotmanhattan <- function(bedfile, chrom, chromstart, chromend, object, leadSNP, scaleLD){

    ## check bedfile columns
    if (!"chr" %in% colnames(bedfile)){
      stop("\'chr\' column not found in data.", call. = FALSE)
    } else {
      if(class(bedfile$chr) != "character"){
        stop("\'chr\' column must be a character.", call. = FALSE)
      }
    }
    if (!"pos" %in% colnames(bedfile)){
      stop("\'pos\' column not found in data.", call. = FALSE)
    } else {
      if(class(bedfile$pos) != "numeric" & class(bedfile$pos) != "integer"){
        stop("\'pos\' column must be an integer or numeric.", call. = FALSE)
      }

    }
    if(!"p" %in% colnames(bedfile)){
      stop("\'p\' column not found in data.", call. = FALSE)
    } else {
      if(class(bedfile$p) != "numeric"){
        stop("\'p\' column must be numeric.", call. = FALSE)
      }
    }


    if (!is.null(chrom)){

      ## Need both chromstart and chromend if trying to do a region within a chrom

      if (!is.null(chromstart) & is.null(chromend)){

        stop("If specifying \'chromstart\', need to provide \'chromend\'.", call. = FALSE)

      }

      if (!is.null(chromend) & is.null(chromstart)){

        stop("If specifying \'chromend\', need to provide \'chromstart\'.", call. = FALSE)

      }

      if (!is.null(chromstart) & !is.null(chromend)){

        ## chromstart cannot be larger than chromend

        if (chromstart == chromend){
          stop("Genomic region is 0 bp long.", call. = FALSE)
        }

        if (chromstart > chromend){

          stop("\'chromstart\' should not be larger than \'chromend\'.", call. = FALSE)
        }

      }

    } else {

      if (!is.null(chromstart) | !is.null(chromend)){
        warning("Plotting multiple chromosomes. \'chromstart\' and \'chromend\' inputs will be ignored.", call. = FALSE)
      }

    }

    ## range
    if (!is.null(object$range)){

      ## range needs to be a vector
      if (!is.vector(object$range)){

        stop("\'range\' must be a vector of length 2.", call. = FALSE)

      }

      ## range vector needs to be length 2
      if (length(man_plot$range) != 2){

        stop("\'range\' must be a vector of length 2.", call. = FALSE)

      }

      ## range vector needs to be numbers
      if (!is.numeric(man_plot$range)){

        stop("\'range\' must be a vector of two numbers.", call. = FALSE)

      }

      ## second value should be larger than the first value
      if (man_plot$range[1] >= man_plot$range[2]){

        stop("\'range\' must be a vector of two numbers in which the 2nd value is larger than the 1st.", call. = FALSE)

      }

    }

    ## lead SNP
    if (!is.null(leadSNP)){
      if (class(leadSNP) != "list"){
        stop("\'leadSNP\' must be a list with a \'snp\' name slot and any other aesthetic options for that SNP,
             like \'fill\', \'pch\', \'cex\', \'fontcolor\', and \'fontsize\'.", call. = FALSE)
      }

      if (!"snp" %in% names(leadSNP)){
        stop("\'leadSNP\' must be a list with a \'snp\' name slot and any other aesthetic options for that SNP,
             like \'fill\', \'pch\', \'cex\', \'fontcolor\', and \'fontsize\'.", call. = FALSE)
      }

    }

    ## scaleLD
    if(!is.null(scaleLD)){
      if(class(scaleLD) != "character"){
        stop("\'scaleLD\' input must be a character specifying the name of a column in data.", call. = FALSE)
      }
      if (!scaleLD %in% colnames(bedfile)){
        stop(paste0(scaleLD, "not found in data."), call. = FALSE)
      }


    }


  }

  ## Define a function that parses the data into an internal format
  parse_data <- function(bedfile, chrom, chromstart, chromend, scaleLD){

    ## Parse 'chr', 'pos', and 'p' columns
    chrCol <- which(colnames(bedfile) == "chr")
    posCol <- which(colnames(bedfile) == "pos")
    pCol <- which(colnames(bedfile) == "p")

    ## Rearrange/reformat columns in case not in order
    data <- data.frame("chr" = bedfile[,chrCol], "pos" = bedfile[,posCol], "p" = bedfile[,pCol])
    data$chr <- as.character(data$chr)
    data$pos <- as.numeric(data$pos)
    data$p <- as.numeric(data$p)

    ## Parse optional 'snp' and scaleLD columns
    if ("snp" %in% colnames(bedfile)){
      snpCol <- which(colnames(bedfile) == "snp")
      data$snp <- bedfile[,snpCol]
    }

    if (!is.null(scaleLD)){
      ldCol <- which(colnames(bedfile) == scaleLD)
      data$ld <- bedfile[,ldCol]
    }


    ## Subset data
    if (!is.null(chrom)){
      data <- data[which(data$chr == chrom),]
      if (!is.null(chromstart) & !is.null(chromend)){
        data <- data[which(data$pos >= chromstart & data$pos <= chromend),]
      }

    }
    return(data)
  }

  ## Define a function that adds offsets to the beddata based on the offsets of the genome assembly
  bed_offset <- function(bedData, offsetAssembly){

    offsetChrom <- function(offsetChrom, bedData){

      chromMatch <- as.character(offsetChrom[1])
      offset <- as.numeric(offsetChrom[3])
      bedMatches <- bedData[which(bedData[,1] == chromMatch),]

      bedMatches[,2] <- bedMatches[,2] + offset

      return(bedMatches)
    }

    updatedBed <- apply(offsetAssembly, 1, offsetChrom, bedData = bedData)
    updatedBed <- dplyr::bind_rows(updatedBed)

    return(updatedBed)

  }

  ## Define a function that adjusts the yrange of the plot
  manhattan_range <- function(bedData, object){

    if (is.null(object$range)){

      object$range <- c(0, object$ymax * max(-log10(bedData[,3])))

    }

    return(object)
  }

  ## Define a function that parses colors
  parse_color <- function(fillcolor, offsetAssembly, bedData){

    if (!is.null(offsetAssembly)){

      ## parse type of color input
      if (class(fillcolor) == "function"){

        newCol <- fillcolor(nrow(offsetAssembly))

      } else {

        newCol <- rep(fillcolor, ceiling(nrow(offsetAssembly)/length(fillcolor)))

      }

      ## Assign associated color number
      colNum <- length(newCol)
      colNum_vector <- rep_len(1:colNum, length.out = nrow(offsetAssembly))

      ## Get color based on number
      colVec <- newCol[colNum_vector]
      offsetAssembly <- cbind(offsetAssembly, colVec)

      ## Get associated chroms in bedfile and assign the color
      chromColor <- function(offsetChrom, bedData){

        chromMatch <- as.character(offsetChrom[1])
        chromCol <- as.character(offsetChrom[5])
        bedMatches <- bedData[which(bedData[,1] == chromMatch),]
        bedCols <- colnames(bedMatches)
        bedMatches <- cbind(bedMatches, rep(chromCol, nrow(bedMatches)))
        colnames(bedMatches) <- c(bedCols, "color")
        return(bedMatches)
      }


      colorBed <- apply(offsetAssembly, 1, chromColor, bedData = bedData)
      colorBed <- dplyr::bind_rows(colorBed)


    } else {


      if (class(fillcolor) == "function") fillcolor <- fillcolor(1)
      bedCols <- colnames(bedData)
      colorBed <- cbind(bedData, rep(fillcolor[1], nrow(bedData)))
      colnames(colorBed) <- c(bedCols, "color")
    }

    return(colorBed)

  }

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(assembly)) assembly <- NULL
  if(missing(fill)) fill <- NULL
  if(missing(pch)) pch <- NULL
  if(missing(space)) space <- NULL
  if(missing(cex)) cex <- NULL
  if(missing(ymax)) ymax <- NULL
  if(missing(sigVal)) sigVal <- NULL
  if(missing(sigLine)) sigLine <- NULL
  if(missing(bg)) bg <- NULL
  if(missing(baseline)) baseline <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL
  if(missing(draw)) draw <- NULL

  ## Check if bed/pVals arguments are missing (could be in object)
  if(!hasArg(data)) data <- NULL

  ## Compile all parameters into an internal object
  bb_manInternal <- structure(list(data = data, leadSNP = leadSNP, chrom = chrom, chromstart = chromstart, chromend = chromend, assembly = assembly,
                                   fill = fill, pch = pch, space = space, cex = cex, ymax = ymax, range = range, sigVal = sigVal, scaleLD = scaleLD,
                                   sigLine = sigLine, sigCol = sigCol, bg = bg, baseline = baseline, x = x, y = y, width = width, height = height, just = just,
                                   default.units = default.units, draw = draw), class = "bb_manInternal")

  bb_manInternal <- parseParams(bb_params = params, object_params = bb_manInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_manInternal$assembly)) bb_manInternal$assembly <- "hg19"
  if(is.null(bb_manInternal$fill)) bb_manInternal$fill <- "black"
  if(is.null(bb_manInternal$pch)) bb_manInternal$pch <- 19
  if(is.null(bb_manInternal$space)) bb_manInternal$space <- 0.01
  if(is.null(bb_manInternal$cex)) bb_manInternal$cex <- 0.25
  if(is.null(bb_manInternal$ymax)) bb_manInternal$ymax <- 1
  if(is.null(bb_manInternal$sigVal)) bb_manInternal$sigVal <- 5e-08
  if(is.null(bb_manInternal$sigLine)) bb_manInternal$sigLine <- FALSE
  if(is.null(bb_manInternal$bg)) bb_manInternal$bg <- NA
  if(is.null(bb_manInternal$baseline)) bb_manInternal$baseline <- FALSE
  if(is.null(bb_manInternal$just)) bb_manInternal$just <- c("left", "top")
  if(is.null(bb_manInternal$default.units)) bb_manInternal$default.units <- "inches"
  if(is.null(bb_manInternal$draw)) bb_manInternal$draw <- TRUE

  ## Set gp
  bb_manInternal$gp <- gpar(cex = bb_manInternal$cex)
  bb_manInternal$gp <- setGP(gpList = bb_manInternal$gp, params = bb_manInternal, ...)

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  man_plot <- structure(list(chrom = bb_manInternal$chrom, chromstart = bb_manInternal$chromstart, chromend = bb_manInternal$chromend, assembly = NULL,
                             range = bb_manInternal$range, ymax = bb_manInternal$ymax, space = bb_manInternal$space, x = bb_manInternal$x, y = bb_manInternal$y,
                             width = bb_manInternal$width, height = bb_manInternal$height, just = bb_manInternal$just, grobs = NULL), class = "bb_manhattan")
  attr(x = man_plot, which = "plotted") <- bb_manInternal$draw

  # ======================================================================================================================================================================================
  # CATCH MISSING ARGUMENT AND PLACEMENT ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_manInternal$data)) stop("argument \"data\" is missing, with no default.", call. = FALSE)

  check_placement(object = man_plot)

  # ======================================================================================================================================================================================
  # READ IN DATA
  # ======================================================================================================================================================================================

  bedfile <- bb_manInternal$data
  ## Read in data if it's not a dataframe or data.table
  if (!"data.frame" %in% class(bedfile)){
    if (!"GRanges" %in% class(bedfile)){
      bedfile <- fread(bedfile)
    }

  }

  bedfile <- as.data.frame(bedfile)

  # ======================================================================================================================================================================================
  # CATCH MORE ERRORS
  # ======================================================================================================================================================================================

  errorcheck_bb_plotmanhattan(bedfile = bedfile, chrom = man_plot$chrom, chromstart = man_plot$chromstart, chromend = man_plot$chromend,
                              object = man_plot, leadSNP = bb_manInternal$leadSNP, scaleLD = bb_manInternal$scaleLD)

  # ======================================================================================================================================================================================
  # PARSE ASSEMBLY
  # ======================================================================================================================================================================================

  man_plot$assembly <- parse_bbAssembly(assembly = bb_manInternal$assembly)

  # ======================================================================================================================================================================================
  # PARSE UNITS
  # ======================================================================================================================================================================================

  man_plot <- defaultUnits(object = man_plot, default.units = bb_manInternal$default.units)

  # ======================================================================================================================================================================================
  # READ AND SUBSET DATA
  # ======================================================================================================================================================================================

  bed_data <- parse_data(bedfile = bedfile, chrom = man_plot$chrom, chromstart = man_plot$chromstart, chromend = man_plot$chromend, scaleLD = bb_manInternal$scaleLD)

  if (nrow(bed_data) > 0){

    # ======================================================================================================================================================================================
    # MULTIPLE CHROMOSOMES
    # ======================================================================================================================================================================================

    if (is.null(man_plot$chrom)){

      man_plot$chromstart <- NULL
      man_plot$chromend <- NULL
      chroms <- as.character(unique(bed_data$chr))

      ## Get chrom sizes based on assembly data
      txdbChecks <- check_loadedPackage(package = man_plot$assembly$TxDb, message = paste(paste0("`", man_plot$assembly$TxDb,"`"), "not loaded. Please install and load to generate full genome assembly Manhattan plot."))
      if (txdbChecks == TRUE){

        tx_db <- eval(parse(text = man_plot$assembly$TxDb))
        assembly_data <- as.data.frame(setDT(as.data.frame(seqlengths(tx_db)), keep.rownames = TRUE))
        assembly_data <- assembly_data[which(assembly_data[,1] %in% chroms),]
        man_plot$chrom <- assembly_data[,1]

        if (any(!chroms %in% assembly_data[,1])){
          non_txdb <- chroms[which(!chroms %in% assembly_data[,1])]
          warning(paste("Chromosome(s)", paste0("'", non_txdb, "'", collapse = ", "), "not found in", paste0("`", man_plot$assembly$TxDb, "`"), "and will be ignored."), call. = FALSE)
        }

        ## get the offsets based on spacer for the assembly
        offsetAssembly <- spaceChroms(assemblyData = assembly_data, space = bb_manInternal$space)

        ## remove bed_data data that aren't in the genome assembly
        bed_data <- bed_data[bed_data$chr %in% offsetAssembly[,1],]

        ## Add chromosome offsets to bed_data
        bed_data <- bed_offset(bedData = bed_data, offsetAssembly = offsetAssembly)

        ## Set viewport xscale
        cumsums <- cumsum(as.numeric(assembly_data[,2]))
        spacer <- cumsums[length(cumsum(as.numeric(assembly_data[,2])))] * bb_manInternal$space

        xscale <- c(0, max(offsetAssembly[,4]) + spacer)

      } else {
        xscale <- c(0, 1)
      }



    } else {

      # ======================================================================================================================================================================================
      # SINGLE CHROMOSOME
      # ======================================================================================================================================================================================
      offsetAssembly <- NULL

      ## Whole single chromosome needs chromosome length information
      if (is.null(man_plot$chromstart) & is.null(man_plot$chromend)){

        txdbChecks <- check_loadedPackage(package = man_plot$assembly$TxDb, message = paste(paste0("`", man_plot$assembly$TxDb,"`"), "not loaded. Please install and load to generate full chromosome Manhattan plot."))
        if (txdbChecks == TRUE){

          tx_db <- eval(parse(text = man_plot$assembly$TxDb))
          assembly_data <- seqlengths(tx_db)

          if (!man_plot$chrom %in% names(assembly_data)){
            warning(paste("Chromosome", paste0("'", man_plot$chrom, "'"), "not found in", paste0("`", man_plot$assembly$TxDb, "`"), "and data for entire chromosome cannot be plotted."), call. = FALSE)
            xscale <- c(0, 1)
          } else {
            man_plot$chromstart <- 1
            man_plot$chromend <- assembly_data[[man_plot$chrom]]
            xscale <- c(man_plot$chromstart, man_plot$chromend)

          }


        }



      } else {
        xscale <- c(man_plot$chromstart, man_plot$chromend)
        txdbChecks <- TRUE
      }


    }

    if (txdbChecks != FALSE){

      # ======================================================================================================================================================================================
      # Y-LIMITS
      # ======================================================================================================================================================================================

      man_plot <- manhattan_range(bedData = bed_data, object = man_plot)

      # ======================================================================================================================================================================================
      # LD COLORING VS SIGNFICANCE COLORING
      # ======================================================================================================================================================================================

      if (!is.null(bb_manInternal$scaleLD)){
        ## Make sure LD column is numeric
        bed_data$ld <- as.numeric(bed_data$ld)

        leadSNP_data <- NULL
        if ("snp" %in% colnames(bed_data)){
          ## Remove lead SNP from data
          leadSNP_data <- bed_data[which(bed_data$snp == bb_manInternal$leadSNP$snp),]
          bed_data <- suppressMessages(dplyr::anti_join(bed_data, leadSNP_data))
        }

        ## Group LD column into LD ranges
        bed_data <- dplyr::group_by(bed_data, LDgrp = cut(bed_data$ld, c(0, 0.2, 0.4, 0.6, 0.8, 1)))

        ## Apply fill and pch to LD groups
        if (class(bb_manInternal$fill) == "function") bb_manInternal$fill <- bb_manInternal$fill(5)
        bed_data$color <- bb_manInternal$fill[bed_data$LDgrp]

        if (length(bb_manInternal$pch) > 1){
          bed_data$pch <- bb_manInternal$pch[bed_data$LDgrp]
        } else{
          bed_data$pch <- rep(bb_manInternal$pch, nrow(bed_data))
        }

        ## If any are NA, set to "grey"
        bed_data[which(is.na(bed_data$color)),]$color <- "grey"

        colorBed <- rbind(bed_data, leadSNP_data)


      } else {

        # ======================================================================================================================================================================================
        # SPLIT DATA INTO SIG AND NON-SIG VALUES FOR ADDITIONAL COLORING
        # ======================================================================================================================================================================================

        sigBed <- bed_data[bed_data$p <= bb_manInternal$sigVal,]
        nonsigBed <- bed_data[bed_data$p > bb_manInternal$sigVal,]

        # ======================================================================================================================================================================================
        # SIGNFICANT COLORING
        # ======================================================================================================================================================================================

        color_nonsig <- parse_color(fillcolor = bb_manInternal$fill, offsetAssembly = offsetAssembly, bedData = nonsigBed)

        if (!is.null(bb_manInternal$sigCol)){

          color_sig <- parse_color(fillcolor = bb_manInternal$sigCol, offsetAssembly = offsetAssembly, bedData = sigBed)

        } else {

          color_sig <- parse_color(fillcolor = bb_manInternal$fill, offsetAssembly = offsetAssembly, bedData = sigBed)
        }


        colorBed <- rbind(color_nonsig, color_sig)
        colorBed$pch <- rep(bb_manInternal$pch[1], nrow(colorBed))



      }



    } else {
      xscale <- c(0, 1)
      man_plot$range <- c(0, 1)
    }

  } else {
    txdbChecks <- TRUE
    xscale <- c(0, 1)
    man_plot$range <- c(0, 1)
    warning("No data found in region.", call. = FALSE)

  }


  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Get viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_manhattan", length(grep(pattern = "bb_manhattan", x = currentViewports)) + 1)

  ## If placing information is provided but plot == TRUE, set up it's own viewport separate from bb_makepage
  ## Not translating into page_coordinates
  if (is.null(man_plot$x) & is.null(man_plot$y)){

    vp <- viewport(height = unit(0.25, "snpc"), width = unit(1, "snpc"),
                   x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                   clip = "on",
                   xscale = xscale, yscale = c(man_plot$range[1], man_plot$range[2]),
                   just = "center",
                   name = vp_name)

    if (bb_manInternal$draw == TRUE){

      vp$name <- "bb_manhattan1"
      grid.newpage()

    }

  } else {

    ## Convert coordinates into same units as page
    page_coords <- convert_page(object = man_plot)

    ## Make viewport
    vp <- viewport(height = page_coords$height, width = page_coords$width,
                   x = page_coords$x, y = page_coords$y,
                   clip = "on",
                   xscale = xscale, yscale = c(man_plot$range[1], man_plot$range[2]),
                   just = bb_manInternal$just,
                   name = vp_name)
  }

  # ======================================================================================================================================================================================
  # INITIALIZE GTREE FOR GROBS WITH BACKGROUND
  # ======================================================================================================================================================================================
  backgroundGrob <- rectGrob(gp = gpar(fill = bb_manInternal$bg, col = NA), name = "background")
  assign("manhattan_grobs", gTree(vp = vp, children = gList(backgroundGrob)), envir = bbEnv)

  if (nrow(bed_data) > 0 & txdbChecks == TRUE){

    # ======================================================================================================================================================================================
    # SUBSET DATA FOR LEAD SNP
    # ======================================================================================================================================================================================
    leadSNP_row <- data.frame(matrix(nrow = 0, ncol = 1))
    if (!is.null(bb_manInternal$leadSNP)){

      if ("snp" %in% colnames(colorBed)){
        # Find index SNP in data
        leadSNP_row <- colorBed[which(colorBed$snp == bb_manInternal$leadSNP$snp),]
        if (nrow(leadSNP_row) > 0){

          ## Remove from data to plot separately
          colorBed <- suppressMessages(dplyr::anti_join(colorBed, leadSNP_row))

        } else {

          warning("Specified lead SNP not found in data.", call. = FALSE)
        }

      }

    }
    # ======================================================================================================================================================================================
    # POINTS
    # ======================================================================================================================================================================================

    if ("col" %in% names(bb_manInternal$gp)){
      bb_manInternal$gp$linecolor <- bb_manInternal$gp$col
    }

    bb_manInternal$gp$col <- colorBed$color


    points <- pointsGrob(x = colorBed$pos, y = -log10(colorBed$p), pch = colorBed$pch,
                         gp = bb_manInternal$gp,
                         default.units = "native")
    assign("manhattan_grobs", addGrob(gTree = get("manhattan_grobs", envir = bbEnv), child = points), envir = bbEnv)

    # ======================================================================================================================================================================================
    # LEAD SNP
    # ======================================================================================================================================================================================

    if (nrow(leadSNP_row) > 0){
      bb_manInternal$gp$col <- bb_manInternal$leadSNP$fill
      bb_manInternal$gp$cex <- bb_manInternal$leadSNP$cex
      if (is.null(bb_manInternal$leadSNP$pch)){
        bb_manInternal$leadSNP$pch <- bb_manInternal$pch[1]
      }
      if (is.null(bb_manInternal$leadSNP$cex)){
        bb_manInternal$gp$cex <- bb_manInternal$cex
      }

      point <- pointsGrob(x = leadSNP_row$pos, y = -log10(leadSNP_row$p), pch = bb_manInternal$leadSNP$pch,
                          gp = bb_manInternal$gp,
                          default.units = "native")
      snp <- textGrob(label = leadSNP_row$snp, x = leadSNP_row$pos, y = unit(-log10(leadSNP_row$p), "native") + unit(1.5, "mm"),
                      just = "bottom", gp = gpar(fontsize = bb_manInternal$leadSNP$fontsize, col = bb_manInternal$leadSNP$fontcolor), default.units = "native")
      assign("manhattan_grobs", addGrob(gTree = get("manhattan_grobs", envir = bbEnv), child = point), envir = bbEnv)
      assign("manhattan_grobs", addGrob(gTree = get("manhattan_grobs", envir = bbEnv), child = snp), envir = bbEnv)

    }

    # ======================================================================================================================================================================================
    # SIGLINE
    # ======================================================================================================================================================================================

    if (bb_manInternal$sigLine == TRUE){

      bb_manInternal$gp$col <- bb_manInternal$gp$linecolor
      sigGrob <- segmentsGrob(x0 = unit(0, "npc"), y0 = unit(-log10(bb_manInternal$sigVal), "native"), x1 = unit(1, "npc"), y1 = unit(-log10(bb_manInternal$sigVal), "native"),
                              gp = bb_manInternal$gp)
      assign("manhattan_grobs", addGrob(gTree = get("manhattan_grobs", envir = bbEnv), child = sigGrob), envir = bbEnv)
    }

    if (bb_manInternal$baseline == TRUE){
      bb_manInternal$gp$col <- bb_manInternal$gp$linecolor
      bb_manInternal$gp$lwd <- 1.5
      baselineGrob <- segmentsGrob(x0 = unit(0, "npc"), y0 = 0, x1 = unit(1, "npc"), y1 = 0,
                              gp = bb_manInternal$gp, default.units = "native")
      assign("manhattan_grobs", addGrob(gTree = get("manhattan_grobs", envir = bbEnv), child = baselineGrob), envir = bbEnv)

    }
  }

  # ======================================================================================================================================================================================
  # IF PLOT == TRUE, DRAW GROBS
  # ======================================================================================================================================================================================

  if (bb_manInternal$draw == TRUE){

    grid.draw(get("manhattan_grobs", envir = bbEnv))

  }

  # ======================================================================================================================================================================================
  # ADD GROBS TO OBJECT
  # ======================================================================================================================================================================================

  man_plot$grobs <- get("manhattan_grobs", envir = bbEnv)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_manhattan[", vp$name, "]"))
  invisible(man_plot)

}
