#' Annotate genomic coordinates along the x or y-axis of a BentoBox plot
#'
#' @param plot Input BentoBox plot to annotate genomic coordinates. Genomic coordinates and
#' assembly will be inherited from \code{plot}.
#' @param fontsize A numeric specifying text fontsize in points. Default value is \code{fontsize = 10}.
#' @param fontcolor A character value indicating the color for text. Default value is \code{fontcolor = "black"}.
#' @param linecolor A character value indicating the color of the genome label axis. Default value is \code{linecolor = "black"}.
#' @param scale A character value indicating the scale of the coordinates along the genome label. Default value is \code{scale = "bp"}. Options are:
#' \itemize{
#' \item{\code{"bp"}: }{base pairs.}
#' \item{\code{"Kb"}: }{kilobase pairs. 1 kilobase pair is equal to 1000 base pairs.}
#' \item{\code{"Mb"}: }{megabase pairs. 1 megabase pair is equal to 1000000 base pairs.}
#' }
#' @param commas A logical value indicating whether to include commas in start and stop labels. Default value is \code{commas = TRUE}.
#' @param sequence A logical value indicating whether to include sequence information above the label of an x-axis (only at appropriate resolutions).
#' @param boxWidth A numeric value indicating the width of the boxes representing sequence information at appropriate resolutions. Default value is \code{boxWidth = 0.5}.
#' @param axis A character value indicating along which axis to add genome label. Sequence information will not be displayed along a y-axis. Default value is \code{axis = "x"}.
#' Options are:
#' \itemize{
#' \item{\code{"x"}: }{Genome label will be plotted along the x-axis.}
#' \item{\code{"y"}: }{Genome label will be plotted along the y-axis. This is typically used for a square Hi-C plot made with \code{bb_plotHicSquare}.}
#' }
#' @param at A numeric vector of x-value locations for tick marks.
#' @param tcl A numeric specifying the length of tickmarks as a fraction of text height. Default value is \code{tcl = 0.5}.
#' @param x A numeric or unit object specifying genome label x-location.
#' @param y A numeric or unit object specifying genome label y-location.
#' @param just Justification of genome label relative to its (x, y) location. If there are two values, the first value specifies horizontal justification and the second value specifies vertical justification.
#' Possible string values are: \code{"left"}, \code{"right"}, \code{"centre"}, \code{"center"}, \code{"bottom"}, and \code{"top"}. Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if \code{x} or \code{y} are only given as numerics. Default value is \code{default.units = "inches"}.
#' @param params An optional \link[BentoBox]{bb_assembly} object containing relevant function parameters.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Returns a \code{bb_genomeLabel} object containing relevant genomic region, placement, and \link[grid]{grob} information.
#'
#' @examples
#' ## Load hg19 genomic annotation packages
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 5, height = 4, default.units = "inches")
#'
#' ## Plot and place gene track on a BentoBox page
#' genesPlot <- bb_plotGenes(chrom = "chr8", chromstart = 1000000, chromend = 2000000,
#'                           assembly = "hg19",
#'                           x = 0.5, y = 0.5, width = 4, height = 1.5, just = c("left", "top"),
#'                           default.units = "inches")
#'
#' ## Annotate x-axis genome label
#' bb_annoGenomeLabel(plot = genesPlot, scale = "Kb",
#'                    x = 0.5, y = 2, just = c("left", "top"), default.units = "inches")
#'
#' @export
bb_annoGenomeLabel <- function(plot, fontsize = 10, fontcolor = "black", linecolor = "black", scale = "bp",
                               commas = TRUE, sequence = TRUE, boxWidth = 0.5, axis = "x", at = NULL,
                               tcl = 0.5, x, y, just = c("left", "top"), default.units = "inches", params = NULL, ...){


  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  if(missing(fontsize)) fontsize <- NULL
  if(missing(fontcolor)) fontcolor <- NULL
  if(missing(linecolor)) linecolor <- NULL
  if(missing(scale)) scale <- NULL
  if(missing(commas)) commas <- NULL
  if(missing(sequence)) sequence <- NULL
  if(missing(boxWidth)) boxWidth <- NULL
  if(missing(axis)) axis <- NULL
  if(missing(tcl)) tcl <- NULL
  if(missing(just)) just <- NULL
  if(missing(default.units)) default.units <- NULL

  ## Check if plot/x/y arguments are missing (could be in object)
  if(!hasArg(plot)) plot <- NULL
  if(!hasArg(x)) x <- NULL
  if(!hasArg(y)) y <- NULL

  ## Compile all parameters into an internal object
  bb_genomeLabelInternal <- structure(list(plot = plot, x = x, y = y, just = just, scale = scale,
                                           fontsize = fontsize, fontcolor = fontcolor, linecolor = linecolor, commas = commas,
                                           sequence = sequence, axis = axis, boxWidth = boxWidth, at = at, tcl = tcl, default.units = default.units), class = "bb_genomeLabelInternal")
  bb_genomeLabelInternal <- parseParams(bb_params = params, object_params = bb_genomeLabelInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_genomeLabelInternal$fontsize)) bb_genomeLabelInternal$fontsize <- 10
  if(is.null(bb_genomeLabelInternal$fontcolor)) bb_genomeLabelInternal$fontcolor <- "black"
  if(is.null(bb_genomeLabelInternal$linecolor)) bb_genomeLabelInternal$linecolor <- "black"
  if(is.null(bb_genomeLabelInternal$scale)) bb_genomeLabelInternal$scale <- "bp"
  if(is.null(bb_genomeLabelInternal$commas)) bb_genomeLabelInternal$commas <- TRUE
  if(is.null(bb_genomeLabelInternal$sequence)) bb_genomeLabelInternal$sequence <- TRUE
  if(is.null(bb_genomeLabelInternal$boxWidth)) bb_genomeLabelInternal$boxWidth <- 0.5
  if(is.null(bb_genomeLabelInternal$axis)) bb_genomeLabelInternal$axis <- "x"
  if(is.null(bb_genomeLabelInternal$tcl)) bb_genomeLabelInternal$tcl <- 0.5
  if(is.null(bb_genomeLabelInternal$just)) bb_genomeLabelInternal$just <- c("left", "top")
  if(is.null(bb_genomeLabelInternal$default.units)) bb_genomeLabelInternal$default.units <- "inches"

  # ======================================================================================================================================================================================
  # CATCH ARGUMENT/PLOT INPUT ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_genomeLabelInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_genomeLabelInternal$x)) stop("argument \"x\" is missing, with no default.", call. = FALSE)
  if(is.null(bb_genomeLabelInternal$y)) stop("argument \"y\" is missing, with no default.", call. = FALSE)

  ## Check that input plot is a valid type of plot to be annotated
  if (class(bb_genomeLabelInternal$plot) != "bb_manhattan"){

    ## Manhattan plots can do whole genome assembly but other plots can't
    inputNames <- attributes(bb_genomeLabelInternal$plot)$names
    if (!("chrom" %in% inputNames) | !("chromstart" %in% inputNames) | !("chromend" %in% inputNames)){

      stop("Invalid input plot. Please input a plot that has genomic coordinates associated with it.", call. = FALSE)

    }

  }

  # ======================================================================================================================================================================================
  # ASSIGN PARAMETERS BASED ON PLOT INPUT
  # ======================================================================================================================================================================================

  chrom <- bb_genomeLabelInternal$plot$chrom
  chromstart <- bb_genomeLabelInternal$plot$chromstart
  chromend <- bb_genomeLabelInternal$plot$chromend
  assembly <- bb_genomeLabelInternal$plot$assembly
  x <- bb_genomeLabelInternal$x
  y <- bb_genomeLabelInternal$y
  length <- bb_genomeLabelInternal$plot$width
  if (bb_genomeLabelInternal$axis == "y"){
    length <- bb_genomeLabelInternal$plot$height
  }
  just <- bb_genomeLabelInternal$just

  ## Whole genome Manhattan plot spacing
  space <- bb_genomeLabelInternal$plot$space

  # ======================================================================================================================================================================================
  # CALL BB_PLOTGENOMELABEL
  # ======================================================================================================================================================================================

  bb_genomeLabel <- bb_plotGenomeLabel(chrom = chrom, chromstart = chromstart, chromend = chromend, assembly = assembly,
                                       fontsize = bb_genomeLabelInternal$fontsize, fontcolor = bb_genomeLabelInternal$fontcolor,
                                       linecolor = bb_genomeLabelInternal$linecolor, scale = bb_genomeLabelInternal$scale,
                                       commas = bb_genomeLabelInternal$commas, sequence = bb_genomeLabelInternal$sequence,
                                       boxWidth = bb_genomeLabelInternal$boxWidth, axis = bb_genomeLabelInternal$axis,
                                       at = bb_genomeLabelInternal$at, tcl = bb_genomeLabelInternal$tcl,
                                       x = x, y = y, length = length, just = just, default.units = bb_genomeLabelInternal$default.units,
                                       space = space, ...)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  return(bb_genomeLabel)

}
