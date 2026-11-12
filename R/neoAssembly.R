#' Read an SV neo-chromosome assembly topology
#'
#' Parse an assembly file from a neoloop/SV workflow and build a
#' `neoAssembly` object describing the concatenated neo-chromosome.
#'
#' @param file Path to an `.assemblies.txt` file.
#' @param id Assembly ID to extract, for example `"C1"`.
#' @param assembly Genome assembly name. Default value is
#' \code{assembly = "hg38"}.
#' @param extend Number of basepairs to extend the first and last fragment.
#' Default value is \code{extend = 0}.
#' @param prefix Prefix to add to chromosome names that do not already have
#' one. Default value is \code{prefix = "chr"}.
#' @param quiet A logical indicating whether to print messages.
#'
#' @return A `neoAssembly` object: a list with \code{id}, \code{fragments},
#' \code{total_length}, and \code{assembly}. \code{fragments} is a data frame
#' with columns \code{id}, \code{chrom}, \code{chrom_raw}, \code{start},
#' \code{end}, \code{flipped}, \code{offset}, and \code{length}.
#'
#' @details
#' Assembly lines are usually produced by NeoLoopFinder's
#' \code{assemble-complexSVs}. Each line contains an assembly ID, zero or more
#' SV records (\code{type,chrA,b1,strandA,chrB,b2,strandB}), and two span
#' coordinates (\code{chr,position}). Consecutive node pairs define fragments;
#' a fragment is flipped when its entry position is greater than its exit
#' position. See \code{vignette("neo_reconstruction")} for a full description
#' and examples.
#'
#' @references
#' Wang X, Xu J, Zhang B, Hou Y, Song F, Lyu H, Yue F. Genome-wide detection of
#' enhancer-hijacking events from chromatin interaction data in rearranged
#' genomes. Nature Methods (2021).
#'
#' @examples
#' \dontrun{
#' neo <- readNeoAssembly(
#'     file = "HM1_merged.assemblies.txt",
#'     id = "C1",
#'     assembly = "hg38",
#'     extend = 200000
#' )
#' neo
#' }
#' @export
readNeoAssembly <- function(file, id, assembly = "hg38", extend = 0,
                            prefix = "chr", quiet = FALSE) {

    if (missing(file) || is.null(file)) {
        stop("argument \"file\" is missing, with no default.", call. = FALSE)
    }
    if (missing(id) || is.null(id)) {
        stop("argument \"id\" is missing, with no default.", call. = FALSE)
    }
    if (!file.exists(file)) {
        stop("Assembly file does not exist: ", file, call. = FALSE)
    }
    if (!is.numeric(extend) || length(extend) != 1 || extend < 0) {
        stop("`extend` must be a single non-negative number.", call. = FALSE)
    }

    lines <- readLines(file, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    id <- as.character(id)
    target <- lines[startsWith(lines, paste0(id, "\t")) |
                    startsWith(lines, paste0(id, " "))]

    if (length(target) == 0) {
        stop("Assembly ID '", id, "' not found in ", file, call. = FALSE)
    }

    parts <- strsplit(target[1], "\t", fixed = TRUE)[[1]]
    parts <- trimws(parts)
    parts <- parts[nzchar(parts)]
    if (length(parts) < 3) {
        stop("Malformed assembly line for ID '", id, "'.", call. = FALSE)
    }

    parseSpan <- function(x) {
        pieces <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
        if (length(pieces) < 2) {
            stop("Malformed assembly span: ", x, call. = FALSE)
        }
        list(chrom = pieces[1], pos = as.numeric(pieces[2]))
    }

    span1 <- parseSpan(parts[length(parts) - 1])
    span2 <- parseSpan(parts[length(parts)])

    svStrings <- character(0)
    if (length(parts) > 3) {
        svStrings <- parts[2:(length(parts) - 2)]
    }

    nodesChrom <- span1$chrom
    nodesPos <- span1$pos

    if (length(svStrings) > 0) {
        for (sv in svStrings) {
            svParts <- trimws(strsplit(sv, ",", fixed = TRUE)[[1]])
            if (length(svParts) < 6) {
                stop("Malformed SV entry: ", sv, call. = FALSE)
            }
            nodesChrom <- c(nodesChrom, svParts[2], svParts[5])
            nodesPos <- c(nodesPos, as.numeric(svParts[3]),
                          as.numeric(svParts[6]))
        }
    }

    nodesChrom <- c(nodesChrom, span2$chrom)
    nodesPos <- c(nodesPos, span2$pos)

    if (length(nodesChrom) %% 2 != 0) {
        stop("Assembly topology has an odd number of nodes.", call. = FALSE)
    }

    nFragments <- length(nodesChrom) / 2
    fragments <- vector("list", nFragments)
    currentOffset <- 0

    formatChrom <- function(chrom) {
        chrom <- as.character(chrom)
        if (startsWith(chrom, "chr")) {
            chrom
        } else {
            paste0(prefix, chrom)
        }
    }

    for (i in seq_len(nFragments)) {
        idx1 <- (i - 1) * 2 + 1
        idx2 <- idx1 + 1
        rawChrom <- nodesChrom[idx1]
        posEnter <- nodesPos[idx1]
        posExit <- nodesPos[idx2]
        flipped <- posEnter > posExit
        start <- min(posEnter, posExit)
        end <- max(posEnter, posExit)

        if (extend > 0 && i == 1) {
            if (flipped) {
                end <- end + extend
            } else {
                start <- max(0, start - extend)
            }
        }
        if (extend > 0 && i == nFragments) {
            if (flipped) {
                start <- max(0, start - extend)
            } else {
                end <- end + extend
            }
        }

        fragLength <- end - start
        if (fragLength <= 0) {
            stop("Fragment ", i, " has a non-positive length.", call. = FALSE)
        }

        fragments[[i]] <- data.frame(
            id = i,
            chrom = formatChrom(rawChrom),
            chrom_raw = as.character(rawChrom),
            start = start,
            end = end,
            flipped = flipped,
            offset = currentOffset,
            length = fragLength,
            stringsAsFactors = FALSE
        )
        currentOffset <- currentOffset + fragLength
    }

    fragments <- do.call(rbind, fragments)
    rownames(fragments) <- NULL

    neo <- structure(
        list(
            id = id,
            file = normalizePath(file, mustWork = TRUE),
            assembly = assembly,
            extend = extend,
            fragments = fragments,
            total_length = currentOffset
        ),
        class = c("neoAssembly", "list")
    )

    if (!quiet) {
        message("neoAssembly[", id, "] ", nFragments, " fragments, ",
                currentOffset, " bp")
    }

    neo
}

#' @export
print.neoAssembly <- function(x, ...) {
    cat("neoAssembly", x$id, "\n")
    cat("  fragments:", nrow(x$fragments), "\n")
    cat("  total length:", x$total_length, "bp\n")
    cat("  assembly:", x$assembly, "\n")
    invisible(x)
}

#' Map genomic features to neo-chromosome coordinates
#'
#' @param data A data frame with columns \code{chrom}, \code{start}, and
#' \code{end}, or a \link[GenomicRanges]{GRanges} object.
#' @param neo A `neoAssembly` object from \code{readNeoAssembly()}.
#'
#' @return A data frame with \code{chrom = "neo"}, neo-coordinate
#' \code{start}/\code{end}, and a \code{fragment} column. Features that span
#' multiple fragments are split across fragments. If a \code{strand} column is
#' present, it is flipped for reverse-complemented fragments.
#'
#' @examples
#' \dontrun{
#' neo <- readNeoAssembly("HM1_merged.assemblies.txt", "C1", extend = 200000)
#' features <- data.frame(
#'     chrom = "chr12", start = 65822483, end = 65968291,
#'     name = "HMGA2", strand = "+"
#' )
#' mapNeoFeatures(features, neo)
#' }
#' @export
mapNeoFeatures <- function(data, neo) {

    if (!inherits(neo, "neoAssembly")) {
        stop("`neo` must be a `neoAssembly` object.", call. = FALSE)
    }

    if (inherits(data, "GRanges")) {
        df <- as.data.frame(data)
        if ("seqnames" %in% colnames(df)) {
            df$chrom <- as.character(df$seqnames)
            df$seqnames <- NULL
        }
    } else {
        df <- as.data.frame(data, stringsAsFactors = FALSE)
    }

    required <- c("chrom", "start", "end")
    if (!all(required %in% colnames(df))) {
        stop("`data` must contain `chrom`, `start`, and `end` columns.",
             call. = FALSE)
    }

    fragments <- neo$fragments
    mappedParts <- list()

    for (i in seq_len(nrow(fragments))) {
        frag <- fragments[i, ]
        subsetDf <- df[df$chrom == frag$chrom | df$chrom == frag$chrom_raw, ,
                       drop = FALSE]
        if (nrow(subsetDf) == 0) {
            next
        }

        subsetDf <- subsetDf[subsetDf$start <= frag$end &
                             subsetDf$end >= frag$start, , drop = FALSE]
        if (nrow(subsetDf) == 0) {
            next
        }

        originalStart <- subsetDf$start
        originalEnd <- subsetDf$end
        clippedStart <- pmax(originalStart, frag$start)
        clippedEnd <- pmin(originalEnd, frag$end)

        if (frag$flipped) {
            subsetDf$start <- frag$offset + (frag$end - clippedEnd)
            subsetDf$end <- frag$offset + (frag$end - clippedStart)
            if ("strand" %in% colnames(subsetDf)) {
                subsetDf$strand <- ifelse(
                    subsetDf$strand == "+", "-",
                    ifelse(subsetDf$strand == "-", "+", subsetDf$strand)
                )
            }
        } else {
            subsetDf$start <- frag$offset + (clippedStart - frag$start)
            subsetDf$end <- frag$offset + (clippedEnd - frag$start)
        }

        subsetDf$chrom <- "neo"
        subsetDf$fragment <- frag$id
        mappedParts[[length(mappedParts) + 1]] <- subsetDf
    }

    if (length(mappedParts) == 0) {
        return(data.frame(
            chrom = character(0), start = numeric(0), end = numeric(0),
            fragment = integer(0)
        ))
    }

    mapped <- do.call(rbind, mappedParts)
    rownames(mapped) <- NULL
    mapped
}
