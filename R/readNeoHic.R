#' @keywords internal
.neoFileType <- function(file) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "hic") {
        "hic"
    } else if (ext %in% c("cool", "mcool")) {
        "cool"
    } else {
        stop("`data` must be a .hic, .cool, or .mcool file.", call. = FALSE)
    }
}

#' @keywords internal
.neoFileResolutions <- function(file, fileType) {
    if (fileType == "hic") {
        strawr::readHicBpResolutions(file)
    } else {
        readCoolBpResolutions(file)
    }
}

#' @keywords internal
.neoAutoResolution <- function(file, fileType, neo, targetBins = 1000) {
    resolutions <- .neoFileResolutions(file, fileType)
    target <- neo$total_length / targetBins
    resolutions[which.min(abs(resolutions - target))]
}

#' @keywords internal
.neoMapPosition <- function(position, fragments, resolution, chrom = NULL) {
    candidates <- fragments
    if (!is.null(chrom)) {
        candidates <- candidates[candidates$chrom == chrom, , drop = FALSE]
    }
    if (nrow(candidates) == 0) {
        return(list(neo = NA_real_, fragment = NA_integer_))
    }

    exact <- which(position >= candidates$start & position < candidates$end)
    if (length(exact) == 0) {
        return(list(neo = NA_real_, fragment = NA_integer_))
    }
    exact <- exact[1]

    frag <- candidates[exact, ]
    if (frag$flipped) {
        neoPosition <- frag$offset + (frag$end - (position + resolution))
    } else {
        neoPosition <- frag$offset + (position - frag$start)
    }

    neoPosition <- max(frag$offset,
                       min(neoPosition, frag$offset + frag$length - resolution))
    list(neo = neoPosition, fragment = frag$id)
}

#' @keywords internal
.neoReadPair <- function(data, fileType, fi, fj, resolution, norm, matrix,
                         assembly) {
    if (fileType == "hic") {
        start1 <- max(1, fi$start)
        end1 <- max(start1 + 1, fi$end)
        start2 <- max(1, fj$start)
        end2 <- max(start2 + 1, fj$end)

        if (fi$id == fj$id) {
            readHic(
                file = data, chrom = fi$chrom,
                chromstart = start1, chromend = end1,
                resolution = resolution, norm = norm, matrix = matrix,
                assembly = assembly, zrange = NULL, quiet = TRUE
            )
        } else {
            readHic(
                file = data, chrom = fi$chrom,
                chromstart = start1, chromend = end1,
                altchrom = fj$chrom,
                altchromstart = start2, altchromend = end2,
                resolution = resolution, norm = norm, matrix = matrix,
                assembly = assembly, zrange = NULL, quiet = TRUE
            )
        }
    } else {
        start1 <- max(0, fi$start)
        end1 <- max(start1 + 1, fi$end)
        start2 <- max(0, fj$start)
        end2 <- max(start2 + 1, fj$end)

        if (fi$id == fj$id) {
            readCool(
                file = data, chrom = fi$chrom,
                chromstart = start1, chromend = end1,
                resolution = resolution, norm = norm,
                zrange = NULL, quiet = TRUE
            )
        } else {
            readCool(
                file = data, chrom = fi$chrom,
                chromstart = start1, chromend = end1,
                altchrom = fj$chrom,
                altchromstart = start2, altchromend = end2,
                resolution = resolution, norm = norm,
                zrange = NULL, quiet = TRUE
            )
        }
    }
}

#' @keywords internal
.neoEstimateResolution <- function(data) {
    starts <- sort(unique(c(data$start1, data$start2)))
    diffs <- diff(starts)
    diffs <- diffs[diffs > 0]
    if (length(diffs) == 0) {
        stop("Cannot estimate `resolution` from the data frame; supply a ",
             "numeric resolution.", call. = FALSE)
    }
    as.numeric(min(diffs))
}

#' @keywords internal
.neoNormalizePairData <- function(data) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    cn <- colnames(data)

    if (all(c("chrom1", "start1", "chrom2", "start2") %in% cn)) {
        countCol <- if ("counts" %in% cn) {
            "counts"
        } else if ("count" %in% cn) {
            "count"
        } else {
            stop("`data` must contain a `counts` or `count` column.",
                 call. = FALSE)
        }
        out <- data.frame(
            chrom1 = as.character(data$chrom1),
            start1 = as.numeric(data$start1),
            chrom2 = as.character(data$chrom2),
            start2 = as.numeric(data$start2),
            counts = as.numeric(data[[countCol]]),
            stringsAsFactors = FALSE
        )
    } else if (all(c("chrom", "start1", "start2") %in% cn)) {
        countCol <- if ("counts" %in% cn) {
            "counts"
        } else if ("count" %in% cn) {
            "count"
        } else {
            stop("`data` must contain a `counts` or `count` column.",
                 call. = FALSE)
        }
        out <- data.frame(
            chrom1 = as.character(data$chrom),
            start1 = as.numeric(data$start1),
            chrom2 = as.character(data$chrom),
            start2 = as.numeric(data$start2),
            counts = as.numeric(data[[countCol]]),
            stringsAsFactors = FALSE
        )
    } else if (ncol(data) >= 5 &&
               is.character(data[[1]]) && is.numeric(data[[2]]) &&
               is.character(data[[3]]) && is.numeric(data[[4]]) &&
               is.numeric(data[[5]])) {
        ## Accept `readCool()`/`readHic()` interchromosomal-like output:
        ## chrom1, start1, chrom2, start2, counts
        out <- data.frame(
            chrom1 = as.character(data[[1]]),
            start1 = as.numeric(data[[2]]),
            chrom2 = as.character(data[[3]]),
            start2 = as.numeric(data[[4]]),
            counts = as.numeric(data[[5]]),
            stringsAsFactors = FALSE
        )
    } else {
        stop("`data` must be a data frame with columns ",
             "`chrom1`, `start1`, `chrom2`, `start2`, `counts`, or ",
             "`chrom`, `start1`, `start2`, `counts`.", call. = FALSE)
    }

    out
}

#' @keywords internal
.neoMakeHic <- function(matrixDf, neo, resolution, norm, matrix, file,
                        zrange, quiet, pairCounts = NULL) {
    if (!is.null(zrange)) {
        matrixDf <- scale_data(upper = matrixDf, zrange = zrange)
    }

    neoHic <- structure(
        list(
            matrix = matrixDf,
            neo = neo,
            resolution = resolution,
            norm = norm,
            matrix_type = matrix,
            file = file,
            pair_counts = pairCounts
        ),
        class = c("neoHic", "list")
    )

    if (!quiet) {
        message("neoHic[", nrow(matrixDf), " pixels at ",
                resolution, " bp]")
    }

    neoHic
}


#' Reconstruct an SV-impacted contact matrix in neo-chromosome coordinates
#'
#' @param data Path to a `.hic`, `.cool`, or `.mcool` file, or a data frame
#' of precomputed genomic contacts with columns \code{chrom1}, \code{start1},
#' \code{chrom2}, \code{start2}, and \code{counts}. A data frame with
#' \code{chrom}, \code{start1}, \code{start2}, and \code{counts} is also
#' accepted for intrachromosomal data.
#' @param neo A `neoAssembly` object from \code{readNeoAssembly()}, or a path
#' to an assembly file. If a path is supplied, \code{id} is required.
#' @param id Assembly ID when \code{neo} is an assembly file path.
#' @param resolution Resolution in basepairs, or \code{"auto"} to choose
#' a resolution close to 1000 bins across the neo-chromosome. For data frame
#' input, \code{"auto"} estimates the smallest positive spacing between bins.
#' @param norm Normalization to pass to \code{readHic()} or \code{readCool()}.
#' Default value is \code{norm = "NONE"}.
#' @param matrix Matrix type for `.hic` files. Default value is
#' \code{matrix = "observed"}.
#' @param zrange Optional numeric vector of length 2 used to clip the final
#' reconstructed counts. Default value is \code{zrange = NULL}.
#' @param assembly Genome assembly for `.hic` files and for reading an
#' assembly file path. Default value is \code{assembly = "hg38"}.
#' @param extend Basepairs to extend the first and last fragments when
#' \code{neo} is a file path. Default value is \code{extend = 0}.
#' @param prefix Prefix for chromosome names when \code{neo} is a file path.
#' Default value is \code{prefix = "chr"}.
#' @param quiet A logical indicating whether to print progress messages.
#'
#' @return A `neoHic` object: a list with \code{matrix} (a data frame of
#' \code{x}, \code{y}, \code{counts} in neo coordinates), \code{neo},
#' \code{resolution}, \code{norm}, and \code{matrix_type}.
#'
#' @details
#' For file input, the reconstruction reads each fragment pair from the source
#' file and maps genomic bin starts onto the concatenated neo-chromosome.
#' Contacts are assigned only to the fragment pair being read, so pixels from
#' adjacent fragments are not double-counted. For data frame input, each row's
#' genomic coordinates are mapped to neo coordinates using the same rules.
#'
#' @examples
#' \dontrun{
#' neo <- readNeoAssembly("HM1_merged.assemblies.txt", "C1", extend = 200000)
#' neoHic <- readNeoHic("HM1_merged.mcool", neo, resolution = 5000)
#' plotNeoHic(neoHic, type = "triangle")
#'
#' ## Or pass the assembly path and ID directly
#' neoHic2 <- readNeoHic(
#'     "HM1_merged.mcool",
#'     "HM1_merged.assemblies.txt", id = "C1",
#'     extend = 200000, resolution = 5000
#' )
#' }
#' @export
readNeoHic <- function(data, neo, id = NULL, resolution = "auto",
                       norm = "NONE", matrix = "observed", zrange = NULL,
                       assembly = "hg38", extend = 0, prefix = "chr",
                       quiet = FALSE) {

    if (missing(data) || is.null(data)) {
        stop("argument \"data\" is missing, with no default.", call. = FALSE)
    }
    if (is.character(neo)) {
        if (is.null(id)) {
            stop("`id` must be supplied when `neo` is an assembly file path.",
                 call. = FALSE)
        }
        neo <- readNeoAssembly(
            file = neo, id = id, assembly = assembly,
            extend = extend, prefix = prefix, quiet = quiet
        )
    }
    if (!inherits(neo, "neoAssembly")) {
        stop("`neo` must be a `neoAssembly` object or an assembly file path.",
             call. = FALSE)
    }

    if (is.data.frame(data)) {
        pairData <- .neoNormalizePairData(data)
        if (identical(resolution, "auto")) {
            resolution <- .neoEstimateResolution(pairData)
        }
        if (!is.numeric(resolution) || length(resolution) != 1 ||
            resolution <= 0) {
            stop("`resolution` must be a positive number or \"auto\".",
                 call. = FALSE)
        }
        resolution <- as.numeric(resolution)

        fragments <- neo$fragments
        mx <- lapply(seq_len(nrow(pairData)), function(k) {
            .neoMapPosition(pairData$start1[k], fragments,
                            resolution = resolution,
                            chrom = pairData$chrom1[k])
        })
        my <- lapply(seq_len(nrow(pairData)), function(k) {
            .neoMapPosition(pairData$start2[k], fragments,
                            resolution = resolution,
                            chrom = pairData$chrom2[k])
        })
        keep <- vapply(seq_len(nrow(pairData)), function(k) {
            !is.na(mx[[k]]$fragment) && !is.na(my[[k]]$fragment)
        }, logical(1))

        if (!any(keep)) {
            matrixDf <- data.frame(
                x = numeric(0), y = numeric(0), counts = numeric(0)
            )
        } else {
            nx <- vapply(mx[keep], `[[`, numeric(1), "neo")
            ny <- vapply(my[keep], `[[`, numeric(1), "neo")
            matrixDf <- data.frame(
                x = pmin(nx, ny),
                y = pmax(nx, ny),
                counts = pairData$counts[keep],
                stringsAsFactors = FALSE
            )
            matrixDf <- aggregate(counts ~ x + y, data = matrixDf, FUN = sum)
        }

        return(.neoMakeHic(
            matrixDf = matrixDf, neo = neo, resolution = resolution,
            norm = norm, matrix = matrix, file = NA_character_,
            zrange = zrange, quiet = quiet
        ))
    }

    if (!is.character(data) || length(data) != 1) {
        stop("`data` must be a file path or data frame.", call. = FALSE)
    }
    if (!file.exists(data)) {
        stop("Data file does not exist: ", data, call. = FALSE)
    }

    fileType <- .neoFileType(data)
    if (identical(resolution, "auto")) {
        resolution <- .neoAutoResolution(
            file = data, fileType = fileType, neo = neo
        )
    }
    if (!is.numeric(resolution) || length(resolution) != 1 ||
        resolution <= 0) {
        stop("`resolution` must be a positive number or \"auto\".",
             call. = FALSE)
    }
    resolution <- as.numeric(resolution)

    fragments <- neo$fragments
    nFragments <- nrow(fragments)
    parts <- list()
    pairStats <- list()

    for (i in seq_len(nFragments)) {
        for (j in i:nFragments) {
            fi <- fragments[i, ]
            fj <- fragments[j, ]
            if (!quiet) {
                message("neo pair ", i, "-", j, ": ",
                        fi$chrom, ":", fi$start, "-", fi$end, " / ",
                        fj$chrom, ":", fj$start, "-", fj$end)
            }

            pairData <- .neoReadPair(
                data = data, fileType = fileType, fi = fi, fj = fj,
                resolution = resolution, norm = norm, matrix = matrix,
                assembly = assembly
            )

            if (nrow(pairData) == 0) {
                pairStats[[length(pairStats) + 1]] <- data.frame(
                    fragment1 = i, fragment2 = j,
                    pixels = 0L, counts = 0, stringsAsFactors = FALSE
                )
                next
            }

            x <- as.numeric(pairData[[1]])
            y <- as.numeric(pairData[[2]])
            counts <- as.numeric(pairData[[3]])

            mx <- lapply(x, .neoMapPosition,
                         fragments = fragments, resolution = resolution,
                         chrom = fi$chrom)
            my <- lapply(y, .neoMapPosition,
                         fragments = fragments, resolution = resolution,
                         chrom = fj$chrom)

            keep <- vapply(seq_along(x), function(k) {
                fx <- mx[[k]]$fragment
                fy <- my[[k]]$fragment
                if (is.na(fx) || is.na(fy)) {
                    FALSE
                } else if (i == j) {
                    fx == i && fy == i
                } else {
                    (fx == i && fy == j) || (fx == j && fy == i)
                }
            }, logical(1))

            if (!any(keep)) {
                pairStats[[length(pairStats) + 1]] <- data.frame(
                    fragment1 = i, fragment2 = j,
                    pixels = 0L, counts = 0, stringsAsFactors = FALSE
                )
                next
            }

            nx <- vapply(mx[keep], `[[`, numeric(1), "neo")
            ny <- vapply(my[keep], `[[`, numeric(1), "neo")

            parts[[length(parts) + 1]] <- data.frame(
                x = pmin(nx, ny),
                y = pmax(nx, ny),
                counts = counts[keep],
                stringsAsFactors = FALSE
            )
            pairStats[[length(pairStats) + 1]] <- data.frame(
                fragment1 = i, fragment2 = j,
                pixels = sum(keep), counts = sum(counts[keep]),
                stringsAsFactors = FALSE
            )
        }
    }

    if (length(parts) == 0) {
        matrixDf <- data.frame(
            x = numeric(0), y = numeric(0), counts = numeric(0)
        )
    } else {
        matrixDf <- do.call(rbind, parts)
        matrixDf <- aggregate(counts ~ x + y, data = matrixDf, FUN = sum)
    }

    pairCounts <- if (length(pairStats) > 0) {
        do.call(rbind, pairStats)
    } else {
        NULL
    }
    if (!quiet && !is.null(pairCounts)) {
        emptyInter <- pairCounts$fragment1 != pairCounts$fragment2 &
            pairCounts$pixels == 0
        if (any(emptyInter)) {
            warning(
                "No inter-fragment contacts found for fragment pairs: ",
                paste0(pairCounts$fragment1[emptyInter], "-",
                       pairCounts$fragment2[emptyInter], collapse = ", "),
                ". The source Hi-C file may have no contacts in those ",
                "fragment intervals.", call. = FALSE
            )
        }
    }

    .neoMakeHic(
        matrixDf = matrixDf, neo = neo, resolution = resolution,
        norm = norm, matrix = matrix,
        file = normalizePath(data, mustWork = TRUE),
        zrange = zrange, quiet = quiet,
        pairCounts = pairCounts
    )
}

#' @export
print.neoHic <- function(x, ...) {
    cat("neoHic\n")
    cat("  matrix:", nrow(x$matrix), "pixels\n")
    cat("  resolution:", x$resolution, "bp\n")
    cat("  norm:", x$norm, "\n")
    cat("  neo length:", x$neo$total_length, "bp\n")
    if (!is.null(x$pair_counts)) {
        cat("  fragment-pair pixels:\n")
        print(x$pair_counts, row.names = FALSE)
    }
    invisible(x)
}

#' Plot a reconstructed neo-chromosome contact matrix
#'
#' @param neoHic A `neoHic` object from \code{readNeoHic()}.
#' @param type Hi-C plot type: \code{"triangle"}, \code{"square"}, or
#' \code{"rectangle"}. Default value is \code{type = "triangle"}.
#' @param flip Orientation. \code{"none"} or \code{FALSE} leaves the plot
#' unchanged, \code{"horizontal"} mirrors the x-axis so the opposite
#' chromosome is on the left, and \code{TRUE} retains the legacy vertical
#' triangle flip. Default value is \code{flip = "none"}.
#' @param ... Additional arguments passed to the corresponding
#' \code{plotHic*} function.
#'
#' @return The plot object returned by \code{plotHicTriangle()},
#' \code{plotHicSquare()}, or \code{plotHicRectangle()}.
#'
#' @examples
#' \dontrun{
#' neoHic <- readNeoHic("HM1_merged.mcool", neo, resolution = 5000)
#' plotNeoHic(neoHic, type = "square", x = 1, y = 1, width = 2, height = 2)
#' }
#' @export
plotNeoHic <- function(neoHic, type = c("triangle", "square", "rectangle"),
                        flip = c("none", "horizontal"), ...) {
    if (!inherits(neoHic, "neoHic")) {
        stop("`neoHic` must be a `neoHic` object.", call. = FALSE)
    }
    type <- match.arg(type)
    flip <- flip[1]
    verticalFlip <- FALSE
    if (is.logical(flip)) {
        verticalFlip <- isTRUE(flip)
        flip <- "none"
    } else {
        flip <- match.arg(flip, c("none", "horizontal"))
    }

    args <- list(
        data = neoHic$matrix,
        chrom = "neo",
        chromstart = 0,
        chromend = neoHic$neo$total_length,
        resolution = neoHic$resolution,
        ...
    )
    if (verticalFlip) {
        args$flip <- TRUE
    }
    if (identical(flip, "horizontal")) {
        if (type == "rectangle") {
            warning("Horizontal flipping is currently supported for triangle ",
                    "and square neo-Hi-C plots.", call. = FALSE)
        } else {
            args$xaxisDir <- "decrease"
        }
    }

    fun <- switch(type,
        triangle = plotHicTriangle,
        square = plotHicSquare,
        rectangle = plotHicRectangle
    )

    do.call(fun, args)
}
