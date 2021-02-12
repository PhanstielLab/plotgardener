#' Add an x-axis to a plot
#'
#' @param plot Plot object to annotate with x-axis.
#' @param at A numeric vector of x-value locations for tick marks.
#' @param label A logical value indicating whether to draw the labels on the tick marks, or an expression or character vector which specify the labels to use. If not logical, must be the same length as the \code{at} argument.
#' @param main A logical value indicating whether to draw the x-axis at the bottom of the plot. Default value is \code{main = TRUE}. Options are:
#' \itemize{
#' \item{\code{TRUE}: }{x-axis is drawn at the bottom of the plot.}
#' \item{\code{FALSE}: }{x-axis is drawn at the top of the plot.}
#' }
#' @param gp Grid graphical parameters. See \link[grid]{gpar}.
#' @param params An optional \link[BentoBox]{bb_params} object containing relevant function parameters.
#'
#' @return Returns a \code{bb_xaxis} object containing relevant \link[grid]{grob} information.
#'
#' @examples
#' ## Load transcript information
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("org.Hs.eg.db")
#'
#' ## Create BentoBox page
#' bb_pageCreate(width = 3, height = 2, default.units = "inches")
#'
#' ## Plot gene transcripts
#' transcriptPlot <- bb_plotTranscripts(chrom = "chr1", chromstart = 1000000, chromend = 2000000,
#'                                      x = 0, y = 0, width = 3, height = 1.5,
#'                                      just = c("left", "top"), default.units = "inches")
#'
#' ## Add standard x-axis to transcript plot
#' bb_annoXaxis(plot = transcriptPlot, at = c(1000000, 1250000, 1500000, 1750000, 2000000),
#'              gp = gpar(fontsize = 8))
#'
#' ## Hide page guides
#' bb_pageGuideHide()
#'
#' @export
bb_annoXaxis <- function(plot, at = NULL, label = TRUE, main = TRUE, gp = gpar(), scipen = 999, axisLine = FALSE, params = NULL){

  # ======================================================================================================================================================================================
  # PARSE PARAMETERS
  # ======================================================================================================================================================================================

  ## Check which defaults are not overwritten and set to NULL
  if(missing(label)) label <- NULL
  if(missing(main)) main <- NULL

  ## Check if plot argument is missing (could be in object)
  if(!hasArg(plot)) plot <- NULL

  ## Compile all parameters into an internal object
  bb_xInternal <- structure(list(plot = plot, at = at, label = label, main = main, gp = gp), class = "bb_xInternal")

  bb_xInternal <- parseParams(bb_params = params, object_params = bb_xInternal)

  ## For any defaults that are still NULL, set back to default
  if(is.null(bb_xInternal$label)) bb_xInternal$label <- TRUE
  if(is.null(bb_xInternal$main)) bb_xInternal$main <- TRUE

  # ======================================================================================================================================================================================
  # CATCH ERRORS
  # ======================================================================================================================================================================================

  if(is.null(bb_xInternal$plot)) stop("argument \"plot\" is missing, with no default.", call. = FALSE)
  check_bbpage(error = "Cannot add an x-axis without a BentoBox page.")

  # ======================================================================================================================================================================================
  # SAVE USER'S SCIPEN AND SET OURS TEMPORARILY
  # ======================================================================================================================================================================================

  usrscipen = getOption("scipen")
  options(scipen = scipen)

  # ======================================================================================================================================================================================
  # INITIALIZE OBJECT
  # ======================================================================================================================================================================================

  xAxis <- structure(list(grobs = NULL), class = "bb_xaxis")

  # ======================================================================================================================================================================================
  # CREATE GROB WITHOUT DRAWING
  # ======================================================================================================================================================================================

  xGrob <- xaxisGrob(at = bb_xInternal$at, label = bb_xInternal$label, main = bb_xInternal$main, gp = bb_xInternal$gp, vp = bb_xInternal$plot$grobs$vp)

  # ======================================================================================================================================================================================
  # GET CENTER OF INPUT PLOT VIEWPORT BASED ON INPUT PLOT TYPE AND JUSTIFICATION
  # ======================================================================================================================================================================================

  if (class(bb_xInternal$plot) == "bb_genes"){

    plotVP <- bb_xInternal$plot$grobs$children$background$vp

  } else if (class(bb_xInternal$plot) == "bb_hicTriangle"){

    plotVP <- bb_xInternal$plot$outsideVP

  } else {

    plotVP <- bb_xInternal$plot$grobs$vp

  }
  adjusted_vp <- adjust_vpCoords(viewport = plotVP)
  # ======================================================================================================================================================================================
  # VIEWPORTS
  # ======================================================================================================================================================================================

  ## Make viewport name
  currentViewports <- current_viewports()
  vp_name <- paste0("bb_xaxis", length(grep(pattern = "bb_xaxis", x = currentViewports)) + 1)

  ## Define viewport
  if (bb_xInternal$main == TRUE){

    vp <- viewport(width = plotVP$width, height = heightDetails(xGrob),
                   x = adjusted_vp[[1]] - 0.5 * (plotVP$width), y = adjusted_vp[[2]] - 0.5 * (plotVP$height),
                   just = c("left", "top"), xscale = plotVP$xscale, name = vp_name)


  } else if (bb_xInternal$main == FALSE){

    vp <- viewport(width = plotVP$width, height = heightDetails(xGrob),
                   x = adjusted_vp[[1]] - 0.5 * (plotVP$width), y = adjusted_vp[[2]] + 0.5 * (plotVP$height),
                   just = c("left", "bottom"), xscale = plotVP$xscale, name = vp_name)

  }

  # ======================================================================================================================================================================================
  # PLOT
  # ======================================================================================================================================================================================

  xGrob <- grid.xaxis(at = bb_xInternal$at, label = bb_xInternal$label, main = bb_xInternal$main, gp = bb_xInternal$gp, vp = vp)
  if(!axisLine) grid.remove(paste0(xGrob$name, "::major"))
  xaxis_grobs <- gTree(vp = vp, children = gList(xGrob))
  xAxis$grobs <- xaxis_grobs

  # ======================================================================================================================================================================================
  # RESTORE USER'S SCIPEN
  # ======================================================================================================================================================================================

  options(scipen = usrscipen)

  # ======================================================================================================================================================================================
  # RETURN OBJECT
  # ======================================================================================================================================================================================

  message(paste0("bb_xaxis[", vp_name, "]"))
  invisible(xAxis)

}
