#' Draw an SV neo-chromosome assembly backbone
#'
#' Draw a directional fragment backbone for a `neoAssembly` object. Each
#' fragment is drawn as an arrow pointing in its neo-chromosome orientation.
#'
#' @param neo A `neoAssembly` object from \code{readNeoAssembly()}.
#' @param x A numeric or unit object specifying the backbone x-location.
#' @param y A numeric, unit object, or character containing a "b" combined
#' with a numeric value specifying the backbone y-location.
#' @param width A numeric or unit object specifying the backbone width.
#' @param height A numeric or unit object specifying the backbone height.
#' @param just Justification of the backbone relative to its (x, y) location.
#' Default value is \code{just = c("left", "top")}.
#' @param default.units A string indicating the default units to use if
#' \code{x}, \code{y}, \code{width}, or \code{height} are only given as
#' numerics. Default value is \code{default.units = "inches"}.
#' @param fill Fill color(s) for fragments. A single value, a vector of
#' length equal to the number of fragments, or a named vector keyed by
#' chromosome name. Default value is \code{fill = "#D4EEF9"}.
#' @param linecolor Fragment outline color. Default value is
#' \code{linecolor = "white"}.
#' @param lwd Fragment outline line width. Default value is \code{lwd = 0.5}.
#' @param fragment.height Fraction of \code{height} used by the fragment
#' arrows. Default value is \code{fragment.height = 0.6}.
#' @param arrow.size Arrow-tip size in page units. If \code{NULL}, the tip is
#' sized relative to each fragment.
#' @param labels Logical value indicating whether to draw fragment labels, or
#' a character vector of labels. Default value is \code{labels = TRUE}.
#' @param label.fontsize Fragment label font size.
#' Default value is \code{label.fontsize = 7}.
#' @param label.fontcolor Fragment label color.
#' Default value is \code{label.fontcolor = "black"}.
#' @param show.seams Logical value indicating whether to draw fragment
#' boundaries. Default value is \code{show.seams = TRUE}.
#' @param seam.color Fragment boundary color.
#' Default value is \code{seam.color = "grey40"}.
#' @param seam.lwd Fragment boundary line width. Default value is
#' \code{seam.lwd = 0.5}.
#' @param seam.lty Fragment boundary line type. Default value is
#' \code{seam.lty = 1}.
#' @param params An optional \link[plotgardener]{pgParams} object.
#' @param ... Additional grid graphical parameters. See \link[grid]{gpar}.
#'
#' @return Invisibly returns a list containing the fragment polygons, labels,
#' and boundary segments.
#'
#' @examples
#' \dontrun{
#' neo <- readNeoAssembly("HM1_merged.assemblies.txt", "C1", extend = 200000)
#' pageCreate(width = 6, height = 2, default.units = "inches")
#' plotNeoAssembly(neo, x = 0.5, y = 0.5, width = 3, height = 0.2)
#' }
#' @export
plotNeoAssembly <- function(neo, x, y, width, height,
                            just = c("left", "top"),
                            default.units = "inches",
                            fill = "#D4EEF9", linecolor = "white", lwd = 0.5,
                            fragment.height = 0.6, arrow.size = NULL,
                            labels = TRUE, label.fontsize = 7,
                            label.fontcolor = "black",
                            show.seams = TRUE, seam.color = "grey40",
                            seam.lwd = 0.5, seam.lty = 1,
                            params = NULL, ...) {

    if (!inherits(neo, "neoAssembly")) {
        stop("`neo` must be a `neoAssembly` object.", call. = FALSE)
    }
    if (missing(x) || missing(y) || missing(width) || missing(height)) {
        stop("`x`, `y`, `width`, and `height` must be supplied.",
             call. = FALSE)
    }

    check_page(error = "Cannot plot a neo assembly without a `plotgardener` page.")

    assemblyInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "neoAssemblyInternal"
    )
    assemblyInternal$just <- justConversion(just = assemblyInternal$just)

    page_units <- get("page_units", envir = pgEnv)

    xNum <- convertX(
        misc_defaultUnits(assemblyInternal$x, "x", assemblyInternal$default.units),
        unitTo = page_units, valueOnly = TRUE
    )
    yNum <- convertY(
        misc_defaultUnits(assemblyInternal$y, "y", assemblyInternal$default.units),
        unitTo = page_units, valueOnly = TRUE
    )
    widthNum <- convertWidth(
        misc_defaultUnits(assemblyInternal$width, "width",
                          assemblyInternal$default.units),
        unitTo = page_units, valueOnly = TRUE
    )
    heightNum <- convertHeight(
        misc_defaultUnits(assemblyInternal$height, "height",
                          assemblyInternal$default.units),
        unitTo = page_units, valueOnly = TRUE
    )

    just <- assemblyInternal$just
    xLeft <- xNum - just[1] * widthNum
    yTop <- yNum - (1 - just[2]) * heightNum

    fragments <- neo$fragments
    nFragments <- nrow(fragments)
    totalLength <- neo$total_length

    fragWidths <- fragments$length / totalLength * widthNum
    fragX1 <- xLeft + fragments$offset / totalLength * widthNum
    fragX2 <- xLeft + (fragments$offset + fragments$length) /
        totalLength * widthNum

    fragmentHeight <- assemblyInternal$fragment.height * heightNum
    if (fragmentHeight <= 0) {
        stop("`fragment.height` must produce a positive fragment height.",
             call. = FALSE)
    }
    yBottom <- yTop + (heightNum - fragmentHeight) / 2
    yCenter <- yTop + heightNum / 2
    yTopArrow <- yBottom + fragmentHeight

    ## Resolve fill colors
    fillColors <- assemblyInternal$fill
    if (length(fillColors) == 1) {
        fillColors <- rep(fillColors, nFragments)
    } else if (!is.null(names(fillColors)) &&
               all(fragments$chrom %in% names(fillColors))) {
        fillColors <- unname(fillColors[fragments$chrom])
    } else {
        fillColors <- rep_len(fillColors, nFragments)
    }

    polygons <- vector("list", nFragments)
    for (i in seq_len(nFragments)) {
        x1 <- fragX1[i]
        x2 <- fragX2[i]
        tip <- if (is.null(assemblyInternal$arrow.size)) {
            min(0.1, (x2 - x1) * 0.15)
        } else {
            assemblyInternal$arrow.size
        }
        tip <- min(tip, (x2 - x1) * 0.45)

        if (fragments$flipped[i]) {
            px <- c(x2, x1 + tip, x1, x1 + tip, x2)
            py <- c(yBottom, yBottom, yCenter, yTopArrow, yTopArrow)
        } else {
            px <- c(x1, x2 - tip, x2, x2 - tip, x1)
            py <- c(yBottom, yBottom, yCenter, yTopArrow, yTopArrow)
        }

        polygons[[i]] <- plotPolygon(
            x = px, y = py,
            fill = fillColors[i],
            linecolor = assemblyInternal$linecolor,
            lwd = assemblyInternal$lwd,
            default.units = page_units, ...
        )
    }

    labelGrobs <- list()
    if (!identical(assemblyInternal$labels, FALSE)) {
        if (isTRUE(assemblyInternal$labels)) {
            labelText <- fragments$chrom_raw
        } else {
            labelText <- rep_len(assemblyInternal$labels, nFragments)
        }
        for (i in seq_len(nFragments)) {
            labelGrobs[[i]] <- plotText(
                label = labelText[i],
                x = (fragX1[i] + fragX2[i]) / 2,
                y = yCenter,
                just = c("center", "center"),
                fontsize = assemblyInternal$label.fontsize,
                fontcolor = assemblyInternal$label.fontcolor,
                default.units = page_units
            )
        }
    }

    seamGrobs <- list()
    if (isTRUE(assemblyInternal$show.seams)) {
        boundaries <- xLeft + cumsum(fragments$length) / totalLength * widthNum
        boundaries <- boundaries[-length(boundaries)]
        for (i in seq_along(boundaries)) {
            seamGrobs[[i]] <- plotSegments(
                x0 = boundaries[i], y0 = yTop,
                x1 = boundaries[i], y1 = yTop + heightNum,
                linecolor = assemblyInternal$seam.color,
                lwd = assemblyInternal$seam.lwd,
                lty = assemblyInternal$seam.lty,
                default.units = page_units
            )
        }
    }

    message("neoAssembly[", neo$id, "]")
    invisible(list(
        fragments = polygons,
        labels = labelGrobs,
        seams = seamGrobs
    ))
}

#' Add flexible labels to an SV neo-chromosome assembly
#'
#' This is a label-only companion to \code{plotNeoAssembly()}. It is useful
#' when fragment labels need custom y-offsets, gaps, fonts, rotation, or
#' inversion highlighting without adding more arguments to the backbone
#' function.
#'
#' @param neo A `neoAssembly` object.
#' @param x A numeric or unit object specifying the label row x-location.
#' @param y A numeric or unit object specifying the label row y-location.
#' @param width A numeric or unit object specifying the label row width.
#' @param height A numeric or unit object specifying the label row height.
#' Default value is \code{height = 0}.
#' @param just Justification of the label row relative to (x, y).
#' Default value is \code{just = c("left", "top")}.
#' @param default.units Default units for numeric \code{x}, \code{y},
#' \code{width}, and \code{height}. Default value is
#' \code{default.units = "inches"}.
#' @param labels Logical value or character vector of labels.
#' If \code{TRUE}, labels are generated from \code{label.style}.
#' Default value is \code{labels = TRUE}.
#' @param label.style Label style when \code{labels = TRUE}: \code{"chrom"}
#' (chromosome name), \code{"coordinate"} (chromosome and Mb coordinates), or
#' \code{"fragment"} (fragment number). Default value is
#' \code{label.style = "chrom"}.
#' @param label.fontsize Label font size. Default value is
#' \code{label.fontsize = 7}.
#' @param label.fontcolor Label color. Can be a vector recycled across
#' fragments. Default value is \code{label.fontcolor = "black"}.
#' @param label.colors Optional explicit per-fragment colors. Overrides
#' \code{label.fontcolor} and inversion coloring.
#' @param label.fontfamily Label font family. Default value is
#' \code{label.fontfamily = ""}.
#' @param label.fontface Label font face. Default value is
#' \code{label.fontface = 1}.
#' @param label.rot Label rotation in degrees. Default value is
#' \code{label.rot = 0}.
#' @param label.position Label position relative to the fragment row:
#' \code{"center"}, \code{"above"}, or \code{"below"}.
#' Default value is \code{label.position = "center"}.
#' @param label.offset Vertical offset in page units for \code{"above"} or
#' \code{"below"} positions. Default value is \code{label.offset = 0.1}.
#' @param show.inversion Logical value indicating whether flipped fragments
#' should be marked. Default value is \code{show.inversion = FALSE}.
#' @param inversion.suffix Suffix appended to flipped-fragment labels.
#' Default value is \code{inversion.suffix = " (inv)"}.
#' @param inversion.color Color for flipped-fragment labels. Default value is
#' \code{inversion.color = "red"}.
#' @param params An optional \link[plotgardener]{pgParams} object.
#' @param ... Additional arguments passed to \code{plotText()}.
#'
#' @return Invisibly returns a list of \code{plotText()} objects.
#'
#' @examples
#' \dontrun{
#' plotNeoAssembly(neo, x = 1, y = 3, width = 3, height = 0.15,
#'                 labels = FALSE)
#' annoNeoFragmentLabels(
#'     neo, x = 1, y = 3.25, width = 3,
#'     label.style = "coordinate",
#'     label.position = "below",
#'     label.offset = 0.05,
#'     label.fontface = "italic",
#'     show.inversion = TRUE
#' )
#' }
#' @export
annoNeoFragmentLabels <- function(neo, x, y, width, height = 0,
                                  just = c("left", "top"),
                                  default.units = "inches",
                                  labels = TRUE, label.style = "chrom",
                                  label.fontsize = 7,
                                  label.fontcolor = "black",
                                  label.colors = NULL,
                                  label.fontfamily = "",
                                  label.fontface = 1,
                                  label.rot = 0,
                                  label.position = "center",
                                  label.offset = 0.1,
                                  show.inversion = FALSE,
                                  inversion.suffix = " (inv)",
                                  inversion.color = "red",
                                  params = NULL, ...) {

    if (!inherits(neo, "neoAssembly")) {
        stop("`neo` must be a `neoAssembly` object.", call. = FALSE)
    }
    if (missing(x) || missing(y) || missing(width)) {
        stop("`x`, `y`, and `width` must be supplied.", call. = FALSE)
    }
    label.position <- match.arg(label.position,
                                c("center", "above", "below"))
    label.style <- match.arg(label.style,
                             c("chrom", "coordinate", "fragment"))
    check_page(error = "Cannot annotate a neo assembly without a `plotgardener` page.")

    labelInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "neoLabelInternal"
    )
    labelInternal$just <- justConversion(just = labelInternal$just)

    page_units <- get("page_units", envir = pgEnv)
    xNum <- convertX(
        misc_defaultUnits(labelInternal$x, "x", labelInternal$default.units),
        unitTo = page_units, valueOnly = TRUE
    )
    yNum <- convertY(
        misc_defaultUnits(labelInternal$y, "y", labelInternal$default.units),
        unitTo = page_units, valueOnly = TRUE
    )
    widthNum <- convertWidth(
        misc_defaultUnits(labelInternal$width, "width",
                          labelInternal$default.units),
        unitTo = page_units, valueOnly = TRUE
    )
    heightNum <- convertHeight(
        misc_defaultUnits(labelInternal$height, "height",
                          labelInternal$default.units),
        unitTo = page_units, valueOnly = TRUE
    )

    just <- labelInternal$just
    xLeft <- xNum - just[1] * widthNum
    yTop <- yNum - (1 - just[2]) * heightNum
    yCenter <- yTop + heightNum / 2

    fragments <- neo$fragments
    nFragments <- nrow(fragments)
    fragWidths <- fragments$length / neo$total_length * widthNum
    fragCenters <- xLeft + (fragments$offset + fragments$length / 2) /
        neo$total_length * widthNum

    if (isTRUE(labels)) {
        labelText <- switch(label.style,
            chrom = fragments$chrom_raw,
            coordinate = {
                left <- ifelse(fragments$flipped,
                               fragments$end, fragments$start) / 1e6
                right <- ifelse(fragments$flipped,
                                fragments$start, fragments$end) / 1e6
                paste0(fragments$chrom_raw, ": ",
                       sprintf("%.2f", left), "-",
                       sprintf("%.2f", right), " Mb")
            },
            fragment = paste0("fragment ", fragments$id)
        )
    } else {
        labelText <- rep_len(labels, nFragments)
    }

    if (isTRUE(show.inversion) && !isTRUE(labels)) {
        labelText[fragments$flipped] <- paste0(
            labelText[fragments$flipped], inversion.suffix
        )
    } else if (isTRUE(show.inversion) && isTRUE(labels)) {
        labelText[fragments$flipped] <- paste0(
            labelText[fragments$flipped], inversion.suffix
        )
    }

    if (!is.null(label.colors)) {
        labelColors <- rep_len(label.colors, nFragments)
    } else {
        labelColors <- rep_len(label.fontcolor, nFragments)
        if (isTRUE(show.inversion)) {
            labelColors[fragments$flipped] <- inversion.color
        }
    }

    yLabel <- switch(label.position,
        center = rep(yCenter, nFragments),
        above = rep(yTop - label.offset, nFragments),
        below = rep(yTop + heightNum + label.offset, nFragments)
    )

    grobs <- lapply(seq_len(nFragments), function(i) {
        plotText(
            label = labelText[i],
            x = fragCenters[i],
            y = yLabel[i],
            just = c("center", "center"),
            fontsize = label.fontsize,
            fontcolor = labelColors[i],
            fontfamily = label.fontfamily,
            fontface = label.fontface,
            rot = label.rot,
            default.units = page_units, ...
        )
    })

    invisible(grobs)
}

#'
#' This mirrors the visual convention used by the neo-triangle examples:
#' a solid domain outline around the full neo-chromosome and dotted domain
#' outlines at each internal fragment seam.
#'
#' @param plot A Hi-C plot object from \code{plotHicSquare()},
#' \code{plotHicTriangle()}, or \code{plotHicRectangle()}.
#' @param neo A `neoAssembly` object.
#' @param outer Logical value indicating whether to draw the outer
#' neo-chromosome domain outline. Default value is \code{outer = TRUE}.
#' @param seams Logical value indicating whether to draw internal fragment
#' seams. Default value is \code{seams = TRUE}.
#' @param linecolor Domain line color. Default value is
#' \code{linecolor = "black"}.
#' @param lwd Domain line width. Default value is \code{lwd = 0.5}.
#' @param lty Outer domain line type. Default value is \code{lty = 1}.
#' @param seam.color Seam line color. Defaults to \code{linecolor}.
#' @param seam.lwd Seam line width. Defaults to \code{lwd}.
#' @param seam.lty Seam line type. Default value is \code{seam.lty = 2}
#' (dotted).
#' @param ... Additional arguments passed to \code{annoDomains()}.
#'
#' @return Invisibly returns a list with \code{outer} and \code{seams}
#' \code{annoDomains} objects.
#'
#' @examples
#' \dontrun{
#' neoHic <- readNeoHic("HM1_merged.mcool", neo, resolution = 5000)
#' h <- plotNeoHic(neoHic, type = "triangle", x = 1, y = 1,
#'                 width = 2, height = 1)
#' annoNeoDomains(h, neo)
#' }
#' @export
annoNeoDomains <- function(plot, neo, outer = TRUE, seams = TRUE,
                           linecolor = "black", lwd = 0.5, lty = 1,
                           seam.color = NULL, seam.lwd = NULL,
                           seam.lty = 2, ...) {

    if (!inherits(neo, "neoAssembly")) {
        stop("`neo` must be a `neoAssembly` object.", call. = FALSE)
    }
    if (!is(plot, "hicSquare") && !is(plot, "hicTriangle") &&
        !is(plot, "hicRectangle")) {
        stop("`plot` must be a Hi-C plot object.", call. = FALSE)
    }

    totalLength <- neo$total_length
    if (is.null(seam.color)) seam.color <- linecolor
    if (is.null(seam.lwd)) seam.lwd <- lwd
    outerGrobs <- NULL
    seamGrobs <- list()

    if (isTRUE(outer)) {
        masterDomain <- data.frame(
            chrom = "neo",
            start = 0,
            end = totalLength
        )
        outerGrobs <- annoDomains(
            plot = plot, data = masterDomain,
            linecolor = linecolor, lwd = lwd, lty = lty, fill = NA, ...
        )
    }

    if (isTRUE(seams)) {
        boundaries <- cumsum(neo$fragments$length)
        boundaries <- boundaries[-length(boundaries)]
        ## Triangle plots mirror the pixel data rather than the viewport, so
        ## seam annotations must be mirrored explicitly to stay aligned.
        if (is(plot, "hicTriangle") &&
            identical(plot$xaxisDir, "decrease")) {
            boundaries <- totalLength - boundaries
        }
        seamGrobs <- lapply(boundaries, function(boundary) {
            seamDomain <- data.frame(
                chrom = c("neo", "neo"),
                start = c(0, boundary),
                end = c(boundary, totalLength)
            )
            annoDomains(
                plot = plot, data = seamDomain,
                linecolor = seam.color, lwd = seam.lwd, lty = seam.lty,
                fill = NA, ...
            )
        })
    }

    invisible(list(outer = outerGrobs, seams = seamGrobs))
}

#' Plot IGV-style gene models in neo-chromosome coordinates
#'
#' Draws gene models through the neo-chromosome topology, using IGV-style
#' display modes. At whole-chromosome scale genes collapse to solid blocks;
#' when zoomed in, gene bodies and exon blocks are shown. Input can be a
#' GTF/GFF file, a `GRanges` object, or a data frame. Using a GTF/GFF file
#' allows arbitrary GENCODE releases without requiring a matching TxDb
#' package.
#'
#' @param data A GTF/GFF path, `GRanges` object, data frame, or a
#' `TxDb`/`EnsDb` annotation object. Data frames need `chrom`, `start`,
#' `end`, and optionally `strand`, `gene`, `transcript`, and `feature`
#' columns.
#' @param neo A `neoAssembly` object.
#' @param txdb Optional `TxDb` or `EnsDb` annotation. If `data` is a
#' `TxDb`/`EnsDb`, it is used automatically.
#' @param orgdb Optional `OrgDb` used to map gene IDs to symbols.
#' @param assembly Assembly name or object used to resolve default
#' `TxDb`/`OrgDb` packages. Defaults to `neo$assembly`.
#' @param gene.id.column,display.column Column names used when reading a
#' `TxDb` annotation.
#' @param geneOrder Optional character vector giving label priority order.
#' @param labelClusters Logical value indicating whether labels should be
#' rotated to reduce overlap. Default value is `FALSE`.
#' @param labelClusterRot Rotation in degrees for clustered labels.
#' Default value is 45.
#' @param geneBackground Color for non-highlighted genes when
#' `geneHighlights` is supplied. Default value is `"grey"`.
#' @param geneMode Gene layout mode. `"expand"` draws all genes at full
#' opacity in packed rows. `"focus"` draws priority genes (`geneHighlights`,
#' `geneOrder`, or named `labels`) at full opacity and all other genes at
#' `focus.alpha`. `"compact"` draws all genes on a single line.
#' Default value is `"expand"`.
#' @param focus.alpha Alpha applied to non-priority genes in focus mode.
#' Default value is 0.5.
#' @param reader GTF/GFF reader: `"auto"`, `"readGFF"`, `"import"`, or
#' `"tabix"`. Default value is `"auto"`.
#' @param cache Optional cache directory for parsed annotations, or `TRUE`
#' for the plotgardener user cache. Default value is `NULL` (no cache).
#' @param x,y,width,height Placement parameters.
#' @param just Justification. Default value is \code{c("left", "top")}.
#' @param default.units Default units. Default value is \code{"inches"}.
#' @param displayMode Gene display mode. \code{"auto"} chooses
#' \code{"collapsed"} or \code{"expanded"} from basepairs per inch.
#' \code{"squished"} uses expanded detail with shorter rows.
#' Default value is \code{displayMode = "auto"}.
#' @param chromstart,chromend Neo-coordinate window to display. Defaults to
#' the full neo-chromosome.
#' @param color Gene color. Can be a single color or a named vector by gene.
#' Default value is \code{"black"}.
#' @param highlight.color Color for highlighted genes.
#' Default value is \code{"red"}.
#' @param geneHighlights Optional character vector of gene names, a named
#' color vector, or a data frame with `gene` and `color` columns. Genes named
#' without colors use \code{highlight.color}.
#' @param gene.type Gene biotypes to retain when a \code{gene_type} column is
#' available, e.g. \code{"protein_coding"}. Use \code{NULL} to keep all
#' biotypes. Default value is \code{"protein_coding"}.
#' @param labels Logical value or character vector of gene names to label.
#' Default value is \code{FALSE} so the track stays clean.
#' @param label.fontsize Gene label font size. Default value is 5.
#' @param label.fontcolor Gene label color for non-highlighted genes.
#' Default value is \code{"black"}.
#' @param label.position Label position: \code{"above"} or \code{"center"}.
#' Default value is \code{"above"}.
#' @param label.overlap Logical value indicating whether labels should be
#' stacked into extra levels to avoid collisions. Default value is
#' \code{TRUE}.
#' @param show.strand Strand display: \code{"arrows"}, \code{"labels"}, or
#' \code{"none"}. Logical \code{TRUE}/\code{FALSE} are accepted for backward
#' compatibility. Default value is \code{"arrows"}.
#' @param strand.color Strand-label color. Default value is \code{"black"}.
#' @param arrow.spacing Spacing between strand chevrons in points.
#' Default value is 30.
#' @param body.lwd Gene-body line width in expanded mode. Default value is
#' 0.8.
#' @param exon.height Exon block height as a fraction of row height.
#' Default value is 0.7.
#' @param feature.min.px Minimum feature width in pixels at 72 dpi.
#' Default value is 3.
#' @param exon.min.px Minimum exon width in pixels at 72 dpi.
#' Default value is 1.
#' @param exon.min.width Deprecated point-based exon width override. When
#' supplied, it takes precedence over \code{exon.min.px}.
#' @param row.height Optional explicit row height in \code{default.units}.
#' @param auto.collapse.bp.per.inch Basepairs per inch above which
#' \code{displayMode = "auto"} collapses genes. Default value is 50000.
#' @param max.rows Maximum number of gene rows. Default value is 4.
#' @param params An optional \link[plotgardener]{pgParams} object.
#' @param ... Additional arguments passed to \code{plotText()}.
#'
#' @return Invisibly returns a list with the gene rows, grobs, display mode,
#' and displayed window.
#' @export
plotNeoGenes <- function(data = NULL, neo, x, y, width, height,
                         txdb = NULL, orgdb = NULL,
                         assembly = NULL,
                         gene.id.column = NULL, display.column = NULL,
                         just = c("left", "top"),
                         default.units = "inches",
                         displayMode = c("auto", "collapsed", "expanded",
                                         "squished"),
                         chromstart = 0, chromend = NULL,
                         color = "black", highlight.color = "red",
                         geneHighlights = NULL, gene.type = "protein_coding",
                         geneOrder = NULL, labelClusters = FALSE,
                         labelClusterRot = 45,
                         geneBackground = "grey",
                         geneMode = c("expand", "focus", "compact"),
                         focus.alpha = 0.5,
                         labels = FALSE, label.fontsize = 5,
                         label.fontcolor = "black",
                         label.position = c("above", "center"),
                         label.overlap = TRUE,
                         show.strand = c("arrows", "labels", "none"),
                         strand.color = "black",
                         arrow.spacing = 30,
                         body.lwd = 0.8, exon.height = 0.7,
                         feature.min.px = 3, exon.min.px = 1,
                         exon.min.width = NULL,
                         row.height = NULL,
                         auto.collapse.bp.per.inch = 50e3,
                         max.rows = 4,
                         reader = c("auto", "readGFF", "import", "tabix"),
                         cache = NULL,
                         params = NULL, ...) {
    if (!inherits(neo, "neoAssembly")) {
        stop("`neo` must be a `neoAssembly` object.", call. = FALSE)
    }

    ## Parse parameters ------------------------------------------------------
    geneInternal <- parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n = 2),
        class = "neoGeneInternal"
    )
    geneInternal$just <- justConversion(just = geneInternal$just)
    if (is.null(geneInternal$assembly)) {
        geneInternal$assembly <- neo$assembly
    }
    assembly <- geneInternal$assembly

    displayMode <- match.arg(
        geneInternal$displayMode,
        c("auto", "collapsed", "expanded", "squished")
    )
    label.position <- match.arg(geneInternal$label.position,
                                c("above", "center"))
    reader <- match.arg(geneInternal$reader,
                        c("auto", "readGFF", "import", "tabix"))
    geneMode <- match.arg(geneInternal$geneMode,
                          c("expand", "focus", "compact"))
    showStrand <- geneInternal$show.strand
    if (is.logical(showStrand)) {
        showStrand <- if (isTRUE(showStrand)) "labels" else "none"
    }
    showStrand <- match.arg(showStrand, c("arrows", "labels", "none"))
    if (is.null(geneInternal$chromend)) {
        geneInternal$chromend <- neo$total_length
    }

    ## Resolve the annotation source ----------------------------------------
    ## A TxDb/EnsDb can be passed as `data`, via `txdb=`, or left to the
    ## default assembly packages. A GTF/GFF path uses the fast reader.
    if (methods::is(data, "TxDb") || inherits(data, "EnsDb")) {
        if (is.null(txdb)) {
            txdb <- data
        }
        data <- NULL
    }

    featureTypes <- c("gene", "exon", "CDS", "UTR",
                      "five_prime_UTR", "three_prime_UTR")
    featureTags <- c("gene_id", "gene_name", "gene_type",
                     "transcript_id", "transcript_type")

    if (is.null(data)) {
        if (is.null(txdb)) {
            stop("Either `data` or `txdb` must be supplied.", call. = FALSE)
        }
        df <- .readNeoGenesTxDb(
            txdb = txdb, orgdb = orgdb, neo = neo,
            assemblyName = geneInternal$assembly,
            gene.id.column = gene.id.column,
            display.column = display.column
        )
    } else if (is.character(data) && length(data) == 1 &&
               file.exists(data)) {
        df <- .readNeoGenesGFF(
            file = data, neo = neo,
            type = featureTypes, tags = featureTags,
            reader = reader, cache = cache
        )
    } else if (inherits(data, "GRanges")) {
        df <- as.data.frame(data)
        df$chrom <- as.character(df$seqnames)
        df$seqnames <- NULL
        if ("type" %in% colnames(df)) {
            df$feature <- as.character(df$type)
        } else if ("feature" %in% colnames(df)) {
            df$feature <- as.character(df$feature)
        } else {
            df$feature <- "exon"
        }
        if (!"gene" %in% colnames(df)) {
            df$gene <- if ("gene_name" %in% colnames(df)) {
                as.character(df$gene_name)
            } else if ("gene_id" %in% colnames(df)) {
                as.character(df$gene_id)
            } else {
                paste0("gene", seq_len(nrow(df)))
            }
        }
        if (!"transcript" %in% colnames(df)) {
            df$transcript <- df$gene
        }
    } else {
        df <- as.data.frame(data, stringsAsFactors = FALSE)
        if (!all(c("chrom", "start", "end") %in% colnames(df))) {
            stop("`data` must contain `chrom`, `start`, and `end` columns.",
                 call. = FALSE)
        }
        if (!"feature" %in% colnames(df)) {
            if ("type" %in% colnames(df)) {
                df$feature <- as.character(df$type)
            } else {
                df$feature <- "exon"
            }
        }
        if (!"gene" %in% colnames(df)) {
            df$gene <- if ("name" %in% colnames(df)) {
                df$name
            } else {
                paste0("gene", seq_len(nrow(df)))
            }
        }
        if (!"transcript" %in% colnames(df)) {
            df$transcript <- df$gene
        }
    }
    if (!"strand" %in% colnames(df)) {
        df$strand <- "*"
    }
    if (!"gene_type" %in% colnames(df) && "gene_biotype" %in% colnames(df)) {
        df$gene_type <- df$gene_biotype
    }
    df$gene <- as.character(df$gene)
    df$feature <- as.character(df$feature)
    df$strand <- as.character(df$strand)

    ## Restrict to requested gene biotypes (GENCODE `gene_type`). TxDb
    ## annotations often have no biotype column; only filter when the column
    ## actually carries values.
    if (!is.null(geneInternal$gene.type) &&
        "gene_type" %in% colnames(df) &&
        any(!is.na(df$gene_type) & nzchar(as.character(df$gene_type)))) {
        df <- df[as.character(df$gene_type) %in% geneInternal$gene.type, ,
                 drop = FALSE]
    }
    if (nrow(df) == 0) {
        return(invisible(list(genes = NULL, grobs = list(),
                              displayMode = displayMode,
                              chromstart = geneInternal$chromstart,
                              chromend = geneInternal$chromend)))
    }

    feature <- toupper(df$feature)
    exonRows <- df[feature == "EXON", , drop = FALSE]
    cdsRows <- df[feature == "CDS", , drop = FALSE]
    geneRows <- df[feature == "GENE", , drop = FALSE]
    if (nrow(exonRows) == 0) {
        exonRows <- geneRows
    }
    if (nrow(exonRows) == 0) {
        exonRows <- df
    }

    mapRows <- function(rows) {
        if (nrow(rows) == 0) {
            return(NULL)
        }
        mapped <- mapNeoFeatures(rows, neo)
        if (nrow(mapped) == 0) {
            return(NULL)
        }
        mapped$gene <- as.character(mapped$gene)
        mapped$strand <- as.character(mapped$strand)
        mapped$key <- paste(mapped$gene, mapped$fragment, mapped$strand,
                            sep = "|")
        mapped
    }
    mappedExons <- mapRows(exonRows)
    mappedCDS <- mapRows(cdsRows)
    if (is.null(mappedExons)) {
        return(invisible(list(genes = NULL, grobs = list(),
                              displayMode = displayMode,
                              chromstart = geneInternal$chromstart,
                              chromend = geneInternal$chromend)))
    }

    exonSplit <- split(mappedExons, mappedExons$key)
    cdsSplit <- if (!is.null(mappedCDS)) split(mappedCDS, mappedCDS$key) else list()
    keys <- names(exonSplit)
    genes <- lapply(keys, function(key) {
        g <- exonSplit[[key]]
        exons <- IRanges::reduce(IRanges::IRanges(g$start, g$end))
        cds <- NULL
        if (!is.null(cdsSplit[[key]])) {
            cg <- cdsSplit[[key]]
            cds <- IRanges::reduce(IRanges::IRanges(cg$start, cg$end))
        }
        data.frame(
            gene = g$gene[1],
            strand = g$strand[1],
            fragment = g$fragment[1],
            start = min(g$start),
            end = max(g$end),
            exonStart = I(list(IRanges::start(exons))),
            exonEnd = I(list(IRanges::end(exons))),
            cdsStart = I(list(if (is.null(cds)) numeric(0) else
                IRanges::start(cds))),
            cdsEnd = I(list(if (is.null(cds)) numeric(0) else
                IRanges::end(cds))),
            nExons = length(exons),
            stringsAsFactors = FALSE
        )
    })
    genes <- do.call(rbind, genes)
    genes <- genes[order(genes$start), , drop = FALSE]

    ## Clip to the displayed neo-coordinate window.
    chromstart <- max(0, geneInternal$chromstart)
    chromend <- min(neo$total_length, geneInternal$chromend)
    if (chromend <= chromstart) {
        stop("`chromend` must be greater than `chromstart`.", call. = FALSE)
    }
    genes <- genes[genes$end >= chromstart & genes$start <= chromend, ,
                   drop = FALSE]
    if (nrow(genes) == 0) {
        return(invisible(list(genes = genes, grobs = list(),
                              displayMode = displayMode,
                              chromstart = chromstart,
                              chromend = chromend)))
    }
    genes$start <- pmax(genes$start, chromstart)
    genes$end <- pmin(genes$end, chromend)

    ## Greedy row packing with explicit overflow handling.
    packGenes <- function(genes, maxRows, overflow = c("drop", "last")) {
        overflow <- match.arg(overflow)
        genes$row <- 0L
        rowEnds <- numeric(0)
        for (i in seq_len(nrow(genes))) {
            placed <- FALSE
            for (r in seq_along(rowEnds)) {
                if (genes$start[i] >= rowEnds[r]) {
                    genes$row[i] <- r
                    rowEnds[r] <- genes$end[i]
                    placed <- TRUE
                    break
                }
            }
            if (!placed) {
                rowEnds <- c(rowEnds, genes$end[i])
                genes$row[i] <- length(rowEnds)
            }
        }
        if (any(genes$row > maxRows)) {
            if (overflow == "last") {
                genes$row[genes$row > maxRows] <- maxRows
            } else {
                warning("Not enough rows for all genes; increase `max.rows`.",
                        call. = FALSE)
                genes <- genes[genes$row <= maxRows, , drop = FALSE]
            }
        }
        genes
    }
    if (geneMode == "compact") {
        genes$row <- 1L
    } else {
        genes <- packGenes(
            genes, geneInternal$max.rows,
            overflow = if (geneMode == "focus") "last" else "drop"
        )
    }
    if (nrow(genes) == 0) {
        return(invisible(list(genes = genes, grobs = list(),
                              displayMode = displayMode,
                              chromstart = chromstart,
                              chromend = chromend)))
    }

    ## Placement in page units ----------------------------------------------
    page_units <- get("page_units", envir = pgEnv)
    xNum <- convertX(misc_defaultUnits(geneInternal$x, "x",
                                       geneInternal$default.units),
                     unitTo = page_units, valueOnly = TRUE)
    yNum <- convertY(misc_defaultUnits(geneInternal$y, "y",
                                       geneInternal$default.units),
                     unitTo = page_units, valueOnly = TRUE)
    widthNum <- convertWidth(misc_defaultUnits(geneInternal$width, "width",
                                               geneInternal$default.units),
                             unitTo = page_units, valueOnly = TRUE)
    heightNum <- convertHeight(misc_defaultUnits(geneInternal$height,
                                                 "height",
                                                 geneInternal$default.units),
                               unitTo = page_units, valueOnly = TRUE)
    just <- geneInternal$just
    xLeft <- xNum - just[1] * widthNum
    yTop <- yNum - (1 - just[2]) * heightNum

    ## Choose display mode from actual basepairs-per-inch.
    widthIn <- convertWidth(unit(widthNum, page_units), unitTo = "inches",
                            valueOnly = TRUE)
    if (widthIn <= 0) {
        stop("`width` must be positive.", call. = FALSE)
    }
    bpPerInch <- (chromend - chromstart) / widthIn
    if (displayMode == "auto") {
        displayMode <- if (bpPerInch > geneInternal$auto.collapse.bp.per.inch) {
            "collapsed"
        } else {
            "expanded"
        }
    }

    nRows <- max(genes$row)
    if (!is.null(geneInternal$row.height)) {
        rowHeight <- convertHeight(
            misc_defaultUnits(geneInternal$row.height, "height",
                              geneInternal$default.units),
            unitTo = page_units, valueOnly = TRUE
        )
    } else {
        rowHeight <- heightNum / nRows
    }
    if (displayMode == "squished") {
        rowHeight <- rowHeight * 0.6
    }
    if (rowHeight * nRows > heightNum) {
        warning("`row.height` exceeds the track height; using equal rows.",
                call. = FALSE)
        rowHeight <- heightNum / nRows
    }
    yOffset <- (heightNum - rowHeight * nRows) / 2

    toX <- function(pos) {
        xLeft + (pos - chromstart) / (chromend - chromstart) * widthNum
    }
    featureMinWidth <- convertWidth(unit(geneInternal$feature.min.px,
                                         "points"),
                                    unitTo = page_units, valueOnly = TRUE)
    if (!is.null(geneInternal$exon.min.width)) {
        exonMinWidth <- convertWidth(unit(geneInternal$exon.min.width,
                                          "points"),
                                     unitTo = page_units, valueOnly = TRUE)
    } else {
        exonMinWidth <- convertWidth(unit(geneInternal$exon.min.px,
                                          "points"),
                                     unitTo = page_units, valueOnly = TRUE)
    }
    exonMinBp <- exonMinWidth / widthNum * (chromend - chromstart)
    ## Gene colors ----------------------------------------------------------
    colorMap <- NULL
    if (!is.null(geneInternal$geneHighlights)) {
        geneHighlights <- geneInternal$geneHighlights
        if (is.data.frame(geneHighlights)) {
            colorMap <- stats::setNames(
                as.character(geneHighlights$color),
                as.character(geneHighlights$gene)
            )
        } else if (!is.null(names(geneHighlights))) {
            colorMap <- as.character(geneHighlights)
            names(colorMap) <- names(geneHighlights)
        } else {
            colorMap <- stats::setNames(
                rep_len(geneInternal$highlight.color, length(geneHighlights)),
                as.character(geneHighlights)
            )
        }
    }
    geneColor <- function(g) {
        if (!is.null(colorMap) && g %in% names(colorMap)) {
            return(unname(colorMap[g]))
        }
        if (!is.null(colorMap)) {
            return(geneInternal$geneBackground)
        }
        if (length(geneInternal$color) > 1 &&
            !is.null(names(geneInternal$color)) &&
            g %in% names(geneInternal$color)) {
            return(unname(geneInternal$color[g]))
        }
        if (length(geneInternal$color) == 1) {
            geneInternal$color
        } else {
            geneInternal$color[1]
        }
    }
    isHighlighted <- vapply(genes$gene, function(g) {
        (!is.null(colorMap) && g %in% names(colorMap)) ||
            (length(geneInternal$color) > 1 &&
             !is.null(names(geneInternal$color)) &&
             g %in% names(geneInternal$color))
    }, logical(1))
    geneCols <- unname(vapply(genes$gene, geneColor, character(1)))

    ## Priority genes keep full exon detail in compact mode.
    priorityNames <- character(0)
    if (!is.null(colorMap)) {
        priorityNames <- c(priorityNames, names(colorMap))
    }
    if (is.character(geneInternal$labels)) {
        priorityNames <- c(priorityNames, geneInternal$labels)
    }
    if (!is.null(geneInternal$geneOrder)) {
        priorityNames <- c(priorityNames, geneInternal$geneOrder)
    }
    priorityNames <- unique(priorityNames)
    priorityIndex <- genes$gene %in% priorityNames
    if (geneMode == "focus" && !any(priorityIndex)) {
        priorityIndex <- rep(TRUE, nrow(genes))
    }

    ## Label layout ---------------------------------------------------------
    ## Priority: explicit `geneOrder`, then highlighted genes, then the rest.
    labelPriority <- seq_len(nrow(genes))
    ord <- if (is.null(geneInternal$geneOrder)) {
        rep(NA_integer_, nrow(genes))
    } else {
        match(genes$gene, geneInternal$geneOrder)
    }
    labelPriority <- order(is.na(ord), ord, !isHighlighted, labelPriority)
    labelPriority <- match(seq_len(nrow(genes)), labelPriority)
    drawLabel <- if (isTRUE(geneInternal$labels)) {
        rep(TRUE, nrow(genes))
    } else if (is.character(geneInternal$labels)) {
        genes$gene %in% geneInternal$labels
    } else {
        rep(FALSE, nrow(genes))
    }
    if (geneMode == "focus") {
        drawLabel <- drawLabel & priorityIndex
    }
    labelLevel <- integer(nrow(genes))
    labelHeight <- convertHeight(
        grobHeight(textGrob("Hg", gp = gpar(
            fontsize = geneInternal$label.fontsize
        ))),
        unitTo = page_units, valueOnly = TRUE
    )
    if (isTRUE(geneInternal$label.overlap) &&
        !isTRUE(geneInternal$labelClusters) &&
        any(drawLabel) && label.position == "above") {
        labelIdx <- which(drawLabel)
        labelCx <- vapply(labelIdx, function(i) {
            toX((genes$start[i] + genes$end[i]) / 2)
        }, numeric(1))
        labelW <- vapply(labelIdx, function(i) {
            convertWidth(
                grobWidth(textGrob(
                    genes$gene[i],
                    gp = gpar(fontsize = geneInternal$label.fontsize)
                )),
                unitTo = page_units, valueOnly = TRUE
            )
        }, numeric(1))
        for (r in unique(genes$row[labelIdx])) {
            sub <- which(genes$row[labelIdx] == r)
            sub <- sub[order(labelPriority[labelIdx[sub]], labelCx[sub])]
            placedCx <- list()
            placedW <- list()
            for (j in sub) {
                lev <- 0L
                repeat {
                    key <- as.character(lev)
                    ok <- TRUE
                    if (!is.null(placedCx[[key]])) {
                        ok <- all(abs(labelCx[j] - placedCx[[key]]) >=
                            (labelW[j] + placedW[[key]]) / 2 + 0.005)
                    }
                    if (ok) break
                    lev <- lev + 1L
                }
                placedCx[[as.character(lev)]] <-
                    c(placedCx[[as.character(lev)]], labelCx[j])
                placedW[[as.character(lev)]] <-
                    c(placedW[[as.character(lev)]], labelW[j])
                labelLevel[labelIdx[j]] <- lev
            }
        }
    }

    ## Drawing helpers ------------------------------------------------------
    clipIntervals <- function(starts, ends) {
        if (length(starts) == 0) {
            return(list(start = numeric(0), end = numeric(0)))
        }
        keep <- ends >= chromstart & starts <= chromend
        list(
            start = pmax(starts[keep], chromstart),
            end = pmin(ends[keep], chromend)
        )
    }
    mergeIntervals <- function(starts, ends, minGap) {
        if (length(starts) == 0) {
            return(list(start = numeric(0), end = numeric(0)))
        }
        ord <- order(starts)
        starts <- starts[ord]
        ends <- ends[ord]
        ms <- starts[1]
        me <- ends[1]
        if (length(starts) > 1) {
            for (e in 2:length(starts)) {
                if (starts[e] - me[length(me)] < minGap) {
                    me[length(me)] <- max(me[length(me)], ends[e])
                } else {
                    ms <- c(ms, starts[e])
                    me <- c(me, ends[e])
                }
            }
        }
        list(start = ms, end = me)
    }

    grobs <- list()
    for (i in seq_len(nrow(genes))) {
        g <- genes[i, ]
        center <- yTop + yOffset + (g$row - 0.5) * rowHeight
        col <- geneCols[i]
        alpha <- if (geneMode == "focus" && !priorityIndex[i]) {
            geneInternal$focus.alpha
        } else {
            1
        }
        x1 <- toX(g$start)
        x2 <- toX(g$end)
        labelTop <- center

        if (displayMode == "collapsed") {
            blockHeight <- rowHeight * 0.45
            geneWidthPx <- x2 - x1
            blockWidth <- max(geneWidthPx, featureMinWidth)
            blockX <- if (blockWidth > geneWidthPx) {
                (x1 + x2) / 2 - blockWidth / 2
            } else {
                x1
            }
            grobs[[length(grobs) + 1]] <- plotRect(
                x = blockX, y = center - blockHeight / 2,
                width = blockWidth, height = blockHeight,
                just = c("left", "top"),
                fill = col, linecolor = col, lwd = 0.2, alpha = alpha,
                default.units = page_units
            )
            labelTop <- center - blockHeight / 2
            if (showStrand == "arrows" && blockWidth > blockHeight) {
                tipSize <- min(blockHeight * 0.7, blockWidth * 0.35)
                tipShift <- tipSize * 0.5
                if (g$strand == "-") {
                    tipX <- blockX - tipShift
                    arrowX <- c(tipX + tipSize, tipX, tipX + tipSize)
                } else {
                    tipX <- blockX + blockWidth + tipShift
                    arrowX <- c(tipX - tipSize, tipX, tipX - tipSize)
                }
                grobs[[length(grobs) + 1]] <- plotPolygon(
                    x = arrowX,
                    y = c(center - blockHeight / 2, center,
                          center + blockHeight / 2),
                    fill = col, linecolor = col, lwd = 0.2, alpha = alpha,
                    default.units = page_units
                )
            }
        } else {
            exons <- clipIntervals(g$exonStart[[1]], g$exonEnd[[1]])
            if (g$nExons <= 1) {
                blockHeight <- rowHeight * geneInternal$exon.height
                blockWidth <- max(x2 - x1, featureMinWidth)
                blockX <- if (blockWidth > x2 - x1) {
                    (x1 + x2) / 2 - blockWidth / 2
                } else {
                    x1
                }
                grobs[[length(grobs) + 1]] <- plotRect(
                    x = blockX, y = center - blockHeight / 2,
                    width = blockWidth, height = blockHeight,
                    just = c("left", "top"),
                    fill = col, linecolor = col, lwd = 0.2, alpha = alpha,
                    default.units = page_units
                )
                labelTop <- center - blockHeight / 2
                if (showStrand == "arrows" && blockWidth > blockHeight) {
                    tipSize <- min(blockHeight * 0.7, blockWidth * 0.35)
                    tipShift <- tipSize * 0.5
                    if (g$strand == "-") {
                        tipX <- blockX - tipShift
                        arrowX <- c(tipX + tipSize, tipX, tipX + tipSize)
                    } else {
                        tipX <- blockX + blockWidth + tipShift
                        arrowX <- c(tipX - tipSize, tipX, tipX - tipSize)
                    }
                    grobs[[length(grobs) + 1]] <- plotPolygon(
                        x = arrowX,
                        y = c(center - blockHeight / 2, center,
                              center + blockHeight / 2),
                        fill = col, linecolor = col, lwd = 0.2, alpha = alpha,
                        default.units = page_units
                    )
                }
            } else {
                grobs[[length(grobs) + 1]] <- plotSegments(
                    x0 = x1, y0 = center, x1 = x2, y1 = center,
                    linecolor = col, lwd = geneInternal$body.lwd,
                    alpha = alpha, default.units = page_units
                )
                fullHalf <- rowHeight * geneInternal$exon.height / 2
                utrHalf <- fullHalf / 2
                exons <- mergeIntervals(exons$start, exons$end, exonMinBp)
                cds <- clipIntervals(g$cdsStart[[1]], g$cdsEnd[[1]])
                cds <- mergeIntervals(cds$start, cds$end, exonMinBp)
                hasCDS <- length(cds$start) > 0

                exonHalf <- if (hasCDS) utrHalf else fullHalf
                for (e in seq_along(exons$start)) {
                    exonWidth <- max(toX(exons$end[e]) - toX(exons$start[e]),
                                     exonMinWidth)
                    grobs[[length(grobs) + 1]] <- plotRect(
                        x = toX(exons$start[e]), y = center - exonHalf,
                        width = exonWidth, height = exonHalf * 2,
                        just = c("left", "top"),
                        fill = col, linecolor = col, lwd = 0.2, alpha = alpha,
                        default.units = page_units
                    )
                }
                if (hasCDS) {
                    for (e in seq_along(cds$start)) {
                        cdsWidth <- max(toX(cds$end[e]) - toX(cds$start[e]),
                                        exonMinWidth)
                        grobs[[length(grobs) + 1]] <- plotRect(
                            x = toX(cds$start[e]), y = center - fullHalf,
                            width = cdsWidth, height = fullHalf * 2,
                            just = c("left", "top"),
                            fill = col, linecolor = col, lwd = 0.2,
                            alpha = alpha, default.units = page_units
                        )
                    }
                }
                labelTop <- center - fullHalf

                if (showStrand == "arrows") {
                    step <- convertWidth(unit(geneInternal$arrow.spacing,
                                              "points"),
                                         unitTo = page_units,
                                         valueOnly = TRUE)
                    if (step > 0 && (x2 - x1) > step) {
                        arrowHalf <- rowHeight * 0.12
                        arrowSize <- min(step * 0.32, rowHeight * 0.24)
                        arrowXs <- seq(x1 + step / 2, x2 - step / 2,
                                       by = step)
                        for (ax in arrowXs) {
                            if (g$strand == "-") {
                                ax0 <- ax + arrowSize / 2
                            } else {
                                ax0 <- ax - arrowSize / 2
                            }
                            grobs[[length(grobs) + 1]] <- plotSegments(
                                x0 = ax0, y0 = center - arrowHalf,
                                x1 = ax, y1 = center,
                                linecolor = col, lwd = 0.3,
                                alpha = alpha, default.units = page_units
                            )
                            grobs[[length(grobs) + 1]] <- plotSegments(
                                x0 = ax, y0 = center,
                                x1 = ax0, y1 = center + arrowHalf,
                                linecolor = col, lwd = 0.3,
                                alpha = alpha, default.units = page_units
                            )
                        }
                    }
                }
            }
        }

        if (showStrand == "labels" &&
            !(geneMode == "focus" && !priorityIndex[i])) {
            strandLabel <- if (g$strand == "-") {
                "-"
            } else if (g$strand == "+") {
                "+"
            } else {
                "*"
            }
            grobs[[length(grobs) + 1]] <- plotText(
                label = strandLabel, x = x1 - 0.02, y = center,
                just = c("right", "center"),
                fontsize = geneInternal$label.fontsize,
                fontcolor = geneInternal$strand.color,
                default.units = page_units, ...
            )
        }

        if (drawLabel[i]) {
            labelCx <- toX((g$start + g$end) / 2)
            if (label.position == "center") {
                labelY <- center
                labelJust <- c("center", "center")
            } else {
                labelY <- labelTop - 0.015 -
                    labelLevel[i] * (labelHeight + 0.005)
                labelJust <- c("center", "bottom")
                if (labelLevel[i] > 0) {
                    grobs[[length(grobs) + 1]] <- plotSegments(
                        x0 = labelCx, y0 = labelTop,
                        x1 = labelCx, y1 = labelY,
                        linecolor = "grey60", lwd = 0.3,
                        default.units = page_units
                    )
                }
            }
            labelCol <- if (isHighlighted[i]) col else
                geneInternal$label.fontcolor
            grobs[[length(grobs) + 1]] <- plotText(
                label = g$gene, x = labelCx, y = labelY,
                just = labelJust,
                fontsize = geneInternal$label.fontsize,
                fontcolor = labelCol,
                rot = if (isTRUE(geneInternal$labelClusters)) {
                    geneInternal$labelClusterRot
                } else {
                    0
                },
                check.overlap = FALSE,
                default.units = page_units, ...
            )
        }
    }

    invisible(list(
        genes = genes, grobs = grobs, displayMode = displayMode,
        chromstart = chromstart, chromend = chromend,
        bpPerInch = bpPerInch, geneMode = geneMode
    ))
}
