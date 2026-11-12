#' Map a signal track onto neo-chromosome coordinates
#'
#' Maps bigWig, `GRanges`, or data-frame signal data through the
#' neo-chromosome topology and optionally bins it in neo coordinates. This
#' handles RNA-seq, ChIP-seq, and other genome-wide signal tracks in the same
#' way, including fragments that are inverted.
#'
#' @param data A bigWig path, `GRanges` object, or data frame with
#' `chrom`, `start`, `end`, and a signal column.
#' @param neo A `neoAssembly` object.
#' @param binSize Bin size in basepairs for the neo-coordinate output. Use
#' `NULL` to keep the input intervals. Default value is `1000`.
#' @param fun Aggregation function for binning: `"mean"`, `"median"`,
#' `"sum"`, `"max"`, or `"min"`. Default value is `"mean"`.
#' @param value.col Name of the signal column. If `NULL`, the function looks
#' for `"score"`, `"balanced"`, `"value"`, or the last numeric column.
#' @param chromstart,chromend Optional neo-coordinate window. Only intervals
#' overlapping this window are returned.
#' @param ... Additional arguments passed to `readBigwig()`.
#'
#' @return A data frame with `chrom`, `start`, `end`, and `score`.
#' @export
getNeoSignal <- function(data, neo, binSize = 1000,
                         fun = c("mean", "median", "sum", "max", "min"),
                         value.col = NULL,
                         chromstart = NULL, chromend = NULL, ...) {
    if (!inherits(neo, "neoAssembly")) {
        stop("`neo` must be a `neoAssembly` object.", call. = FALSE)
    }
    fun <- match.arg(fun)

    ## Read a bigWig fragment by fragment, or accept pre-read data --------
    if (is.character(data) && length(data) == 1) {
        if (!file.exists(data)) {
            stop("Signal file does not exist: ", data, call. = FALSE)
        }
        isBigWig <- grepl("\\.(bw|bigwig|bigWig)(\\.gz)?$", data,
                          ignore.case = TRUE)
        if (isBigWig) {
            parts <- lapply(seq_len(nrow(neo$fragments)), function(i) {
                frag <- neo$fragments[i, ]
                bw <- readBigwig(
                    file = data, chrom = frag$chrom,
                    chromstart = max(1, frag$start), chromend = frag$end,
                    ...
                )
                if (nrow(bw) == 0) {
                    return(NULL)
                }
                df <- as.data.frame(bw)
                if ("seqnames" %in% colnames(df)) {
                    df$chrom <- as.character(df$seqnames)
                    df$seqnames <- NULL
                }
                mapNeoFeatures(df, neo)
            })
            parts <- Filter(NROW, parts)
            if (length(parts) == 0) {
                return(data.frame(
                    chrom = character(0), start = numeric(0),
                    end = numeric(0), score = numeric(0)
                ))
            }
            raw <- do.call(rbind, parts)
        } else {
            ## Text table input (for example a virtual-4C table).
            if (!requireNamespace("data.table", quietly = TRUE)) {
                stop("`data.table` is required to read tabular signal files.",
                     call. = FALSE)
            }
            raw <- as.data.frame(data.table::fread(data))
            raw <- .neoSignalFrame(raw, neo, value.col)
        }
    } else if (inherits(data, "GRanges")) {
        raw <- as.data.frame(data)
        raw <- .neoSignalFrame(raw, neo, value.col)
    } else {
        raw <- .neoSignalFrame(as.data.frame(data), neo, value.col)
    }

    if (nrow(raw) == 0) {
        return(data.frame(
            chrom = character(0), start = numeric(0),
            end = numeric(0), score = numeric(0)
        ))
    }

    if (is.null(binSize)) {
        out <- data.frame(
            chrom = "neo", start = raw$start, end = raw$end,
            score = raw$score, stringsAsFactors = FALSE
        )
    } else {
        if (!is.numeric(binSize) || length(binSize) != 1 || binSize <= 0) {
            stop("`binSize` must be a positive number or NULL.",
                 call. = FALSE)
        }
        raw$bin <- floor(raw$start / binSize) * binSize
        binned <- stats::aggregate(
            score ~ bin, data = raw, FUN = fun, na.rm = TRUE
        )
        out <- data.frame(
            chrom = "neo",
            start = binned$bin,
            end = binned$bin + binSize - 1,
            score = binned$score,
            stringsAsFactors = FALSE
        )
    }

    out <- out[order(out$start), , drop = FALSE]
    rownames(out) <- NULL
    if (!is.null(chromstart) || !is.null(chromend)) {
        lo <- if (is.null(chromstart)) -Inf else chromstart
        hi <- if (is.null(chromend)) Inf else chromend
        out <- out[out$end >= lo & out$start <= hi, , drop = FALSE]
    }
    out
}

.neoSignalFrame <- function(df, neo, value.col = NULL) {
    if (!all(c("chrom", "start", "end") %in% colnames(df)) &&
        !all(c("seqnames", "start", "end") %in% colnames(df))) {
        stop("Signal data must contain `chrom`, `start`, and `end` columns.",
             call. = FALSE)
    }
    if ("seqnames" %in% colnames(df)) {
        df$chrom <- as.character(df$seqnames)
        df$seqnames <- NULL
    }
    if (!"chrom" %in% colnames(df)) {
        stop("Signal data must contain a `chrom` column.", call. = FALSE)
    }

    if (!is.null(value.col)) {
        if (!value.col %in% colnames(df)) {
            stop("`value.col` not found in signal data: ", value.col,
                 call. = FALSE)
        }
        score <- df[[value.col]]
    } else {
        candidates <- c("score", "balanced", "value")
        hit <- candidates[candidates %in% colnames(df)][1]
        if (is.na(hit)) {
            numericCols <- names(df)[vapply(df, is.numeric, logical(1))]
            numericCols <- setdiff(numericCols, c("start", "end"))
            hit <- utils::tail(numericCols, 1)
        }
        if (is.na(hit) || length(hit) == 0) {
            stop("Could not determine a signal column; supply `value.col`.",
                 call. = FALSE)
        }
        score <- df[[hit]]
    }

    mapped <- mapNeoFeatures(
        data.frame(
            chrom = as.character(df$chrom),
            start = as.numeric(df$start),
            end = as.numeric(df$end),
            score = as.numeric(score),
            stringsAsFactors = FALSE
        ),
        neo
    )
    mapped[, c("chrom", "start", "end", "score"), drop = FALSE]
}

#' Map a virtual-4C table onto neo-chromosome coordinates
#'
#' Convenience wrapper around `getNeoSignal()` for cooltools virtual-4C TSV
#' files. The table is mapped through the neo topology and adjacent intervals
#' are clipped so they do not overlap.
#'
#' @param file Path to a virtual-4C TSV file.
#' @param neo A `neoAssembly` object.
#' @param value.col Signal column. Default value is `"balanced"`.
#' @param binSize Optional bin size. Default value is `NULL`, which keeps the
#' table's own intervals.
#' @param fun Aggregation function when `binSize` is supplied.
#' @param ... Additional arguments passed to `getNeoSignal()`.
#'
#' @return A data frame with `chrom`, `start`, `end`, and `score`.
#' @export
getNeoV4C <- function(file, neo, value.col = "balanced", binSize = NULL,
                      fun = "mean", ...) {
    out <- getNeoSignal(
        data = file, neo = neo, binSize = binSize,
        fun = fun, value.col = value.col, ...
    )
    if (nrow(out) == 0) {
        return(out)
    }
    out <- out[order(out$start), , drop = FALSE]
    out$end <- pmin(
        out$end - 1,
        c(out$start[-1] - 1, Inf)
    )
    out <- out[out$start <= out$end, , drop = FALSE]
    rownames(out) <- NULL
    out
}
