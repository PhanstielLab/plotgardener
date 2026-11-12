## Internal gene-annotation readers for plotNeoGenes().
##
## The GTF/GFF path is intentionally separate from the renderer so that the
## expensive parse can be cached and reused across plots/samples.

.neoGeneSessionCache <- new.env(parent = emptyenv())

.neoGeneRegions <- function(neo) {
    GenomicRanges::GRanges(
        seqnames = neo$fragments$chrom,
        ranges = IRanges::IRanges(
            start = neo$fragments$start,
            end = neo$fragments$end
        )
    )
}

.neoGeneRegionFrame <- function(neo) {
    data.frame(
        chrom = as.character(neo$fragments$chrom),
        chrom_raw = as.character(neo$fragments$chrom_raw),
        start = as.numeric(neo$fragments$start),
        end = as.numeric(neo$fragments$end),
        stringsAsFactors = FALSE
    )
}

.neoGeneCacheKey <- function(file, regions, type, tags) {
    info <- file.info(file)
    key <- paste(
        normalizePath(file, mustWork = FALSE),
        ifelse(is.na(info$size), "", info$size),
        ifelse(is.na(info$mtime), "", as.numeric(info$mtime)),
        paste(regions$chrom, regions$chrom_raw, regions$start, regions$end,
              sep = ":", collapse = ";"),
        paste(type, collapse = ","),
        paste(tags, collapse = ","),
        sep = "|"
    )
    if (requireNamespace("digest", quietly = TRUE)) {
        return(digest::digest(key, algo = "xxhash64"))
    }
    key <- gsub("[^A-Za-z0-9._-]+", "_", key)
    substr(key, 1, 180)
}

.neoGeneCacheGet <- function(cache, key) {
    if (is.null(cache)) {
        return(NULL)
    }
    if (isTRUE(cache)) {
        dir <- tools::R_user_dir("plotgardener", which = "cache")
    } else {
        dir <- cache
    }
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    path <- file.path(dir, paste0("neo_genes_", key, ".rds"))
    if (file.exists(path)) {
        return(readRDS(path))
    }
    NULL
}

.neoGeneCacheSet <- function(cache, key, value) {
    if (is.null(cache)) {
        return(invisible(NULL))
    }
    if (isTRUE(cache)) {
        dir <- tools::R_user_dir("plotgardener", which = "cache")
    } else {
        dir <- cache
    }
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    path <- file.path(dir, paste0("neo_genes_", key, ".rds"))
    saveRDS(value, path)
    invisible(path)
}

.neoGeneSubsetRegions <- function(df, regions) {
    if (nrow(df) == 0) {
        return(df)
    }
    keep <- rep(FALSE, nrow(df))
    for (i in seq_len(nrow(regions))) {
        keep <- keep | (
            (df$chrom == regions$chrom[i] |
             df$chrom == regions$chrom_raw[i]) &
            df$start <= regions$end[i] &
            df$end >= regions$start[i]
        )
    }
    df[keep, , drop = FALSE]
}

.neoGeneNormalize <- function(df, regions = NULL) {
    if (is.null(df) || nrow(df) == 0) {
        return(data.frame(
            chrom = character(0), start = numeric(0), end = numeric(0),
            strand = character(0), gene = character(0),
            transcript = character(0), feature = character(0),
            gene_type = character(0), stringsAsFactors = FALSE
        ))
    }
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    chromCol <- intersect(c("chrom", "seqid", "seqnames"), colnames(df))[1]
    featureCol <- intersect(c("feature", "type"), colnames(df))[1]
    geneCol <- intersect(c("gene", "gene_name", "name"), colnames(df))[1]
    idCol <- intersect(c("gene_id", "geneId", "GENEID"), colnames(df))[1]
    txCol <- intersect(c("transcript", "transcript_id", "TXNAME", "TXID"),
                       colnames(df))[1]
    typeCol <- intersect(c("gene_type", "gene_biotype", "TXTYPE"),
                         colnames(df))[1]

    chrom <- if (!is.na(chromCol)) as.character(df[[chromCol]]) else NA_character_
    gene <- if (!is.na(geneCol)) as.character(df[[geneCol]]) else NA_character_
    if (!is.na(idCol)) {
        ids <- as.character(df[[idCol]])
        gene[is.na(gene) | !nzchar(gene)] <- ids[is.na(gene) | !nzchar(gene)]
    }
    gene[is.na(gene) | !nzchar(gene)] <- paste0("gene", seq_len(nrow(df)))[
        is.na(gene) | !nzchar(gene)
    ]
    transcript <- if (!is.na(txCol)) as.character(df[[txCol]]) else gene
    transcript[is.na(transcript) | !nzchar(transcript)] <-
        gene[is.na(transcript) | !nzchar(transcript)]

    out <- data.frame(
        chrom = chrom,
        start = as.numeric(df$start),
        end = as.numeric(df$end),
        strand = if ("strand" %in% colnames(df)) {
            as.character(df$strand)
        } else {
            "*"
        },
        gene = gene,
        transcript = transcript,
        feature = if (!is.na(featureCol)) {
            as.character(df[[featureCol]])
        } else {
            "exon"
        },
        gene_type = if (!is.na(typeCol)) {
            as.character(df[[typeCol]])
        } else {
            NA_character_
        },
        stringsAsFactors = FALSE
    )
    out$strand[is.na(out$strand) | !nzchar(out$strand)] <- "*"
    out <- out[!is.na(out$chrom) & !is.na(out$start) & !is.na(out$end), ,
               drop = FALSE]

    ## GENCODE carries gene_type on every row, but some GTF/GFF files only
    ## put it on the gene rows. Fill it down by gene_id/gene name.
    hasType <- !is.na(out$gene_type) & nzchar(as.character(out$gene_type))
    if (any(hasType) && any(!hasType)) {
        typeMap <- tapply(
            as.character(out$gene_type[hasType]),
            out$gene[hasType],
            function(x) x[1]
        )
        fillIdx <- !hasType & out$gene %in% names(typeMap)
        out$gene_type[fillIdx] <- unname(typeMap[out$gene[fillIdx]])
    }

    ## Drop duplicate feature rows (common in TxDb joins).
    key <- paste(out$chrom, out$start, out$end, out$strand, out$gene,
                 out$transcript, out$feature, sep = "|")
    out <- out[!duplicated(key), , drop = FALSE]

    if (!is.null(regions)) {
        out <- .neoGeneSubsetRegions(out, regions)
    }
    out
}

.neoGeneReadGFF <- function(file, regions, type, tags, reader) {
    regionsGR <- GenomicRanges::GRanges(
        seqnames = regions$chrom,
        ranges = IRanges::IRanges(regions$start, regions$end)
    )
    if (reader == "readGFF") {
        if (!requireNamespace("rtracklayer", quietly = TRUE)) {
            stop("`rtracklayer` is required to read GTF/GFF files.",
                 call. = FALSE)
        }
        ## `readGFF` applies the seqid/type filters during the native scan and
        ## only materializes the requested columns/tags. It is much lighter
        ## than import(..., which = regions), which parses the whole file into
        ## a GRanges first.
        df <- suppressWarnings(rtracklayer::readGFF(
            file,
            filter = list(
                seqid = unique(c(regions$chrom, regions$chrom_raw)),
                type = type
            ),
            columns = c("seqid", "start", "end", "strand", "type"),
            tags = tags
        ))
        return(.neoGeneNormalize(df, regions))
    }

    if (reader == "tabix") {
        if (!requireNamespace("Rsamtools", quietly = TRUE) ||
            !requireNamespace("rtracklayer", quietly = TRUE)) {
            stop("`Rsamtools` and `rtracklayer` are required for the tabix ",
                 "reader.", call. = FALSE)
        }
        if (!file.exists(paste0(file, ".tbi"))) {
            stop("Tabix index not found: ", paste0(file, ".tbi"),
                 call. = FALSE)
        }
        format <- if (grepl("\\.gff3?(\\.gz)?$", file, ignore.case = TRUE)) {
            "GFF3"
        } else {
            "GTF"
        }
        tbx <- Rsamtools::TabixFile(file)
        present <- Rsamtools::seqnamesTabix(tbx)
        queryGR <- regionsGR[
            as.character(GenomicRanges::seqnames(regionsGR)) %in% present
        ]
        if (length(queryGR) == 0) {
            return(.neoGeneNormalize(NULL))
        }
        gr <- suppressWarnings(rtracklayer::import(
            tbx, format = format, which = queryGR,
            feature.type = type
        ))
        return(.neoGeneNormalize(as.data.frame(gr), regions))
    }

    ## Fallback: the original full-file import.
    if (!requireNamespace("rtracklayer", quietly = TRUE)) {
        stop("`rtracklayer` is required to read GTF/GFF files.",
             call. = FALSE)
    }
    gr <- suppressMessages(rtracklayer::import(
        file, feature.type = type, which = regionsGR
    ))
    .neoGeneNormalize(as.data.frame(gr), regions)
}

.readNeoGenesGFF <- function(file, neo, type, tags,
                             reader = c("auto", "readGFF", "import",
                                        "tabix"),
                             cache = NULL) {
    reader <- match.arg(reader)
    regions <- .neoGeneRegionFrame(neo)
    if (reader == "auto") {
        if (file.exists(paste0(file, ".tbi")) &&
            requireNamespace("Rsamtools", quietly = TRUE)) {
            reader <- "tabix"
        } else if (requireNamespace("rtracklayer", quietly = TRUE)) {
            reader <- "readGFF"
        } else {
            reader <- "import"
        }
    }

    key <- .neoGeneCacheKey(file, regions, type, tags)
    cached <- .neoGeneCacheGet(cache, key)
    if (!is.null(cached)) {
        return(cached)
    }
    sessionKey <- paste(reader, key, sep = "|")
    if (exists(sessionKey, envir = .neoGeneSessionCache, inherits = FALSE)) {
        cached <- get(sessionKey, envir = .neoGeneSessionCache,
                      inherits = FALSE)
        .neoGeneCacheSet(cache, key, cached)
        return(cached)
    }

    df <- .neoGeneReadGFF(file, regions, type, tags, reader)
    assign(sessionKey, df, envir = .neoGeneSessionCache)
    .neoGeneCacheSet(cache, key, df)
    df
}

.readNeoGenesTxDb <- function(txdb, orgdb, neo, assemblyName,
                              gene.id.column = NULL,
                              display.column = NULL) {
    if (!requireNamespace("GenomicFeatures", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
        stop("`GenomicFeatures` and `AnnotationDbi` are required for TxDb ",
             "gene annotation.", call. = FALSE)
    }
    if (is.null(display.column)) {
        display.column <- if (is.null(orgdb)) "GENEID" else "SYMBOL"
    }
    if (is.null(gene.id.column)) {
        gene.id.column <- if (identical(display.column, "GENEID")) {
            "GENEID"
        } else {
            "ENTREZID"
        }
    }
    asm <- assembly(
        Genome = assemblyName,
        TxDb = txdb,
        OrgDb = if (is.null(orgdb)) NA else orgdb,
        gene.id.column = gene.id.column,
        display.column = display.column
    )
    regions <- .neoGeneRegionFrame(neo)
    chroms <- unique(regions$chrom)
    parts <- lapply(chroms, function(ch) {
        rr <- regions[regions$chrom == ch, , drop = FALSE]
        getExons(
            assembly = asm,
            chromosome = ch,
            start = min(rr$start),
            stop = max(rr$end)
        )
    })
    parts <- Filter(function(x) !is.null(x) && nrow(x) > 0, parts)
    if (length(parts) == 0) {
        return(.neoGeneNormalize(NULL))
    }
    d <- do.call(rbind, parts)

    displayCol <- if (display.column %in% colnames(d)) {
        display.column
    } else if ("SYMBOL" %in% colnames(d)) {
        "SYMBOL"
    } else {
        "GENEID"
    }
    gene <- as.character(d[[displayCol]])
    gene[is.na(gene) | !nzchar(gene)] <- as.character(d$GENEID)[
        is.na(gene) | !nzchar(gene)
    ]
    transcript <- if ("TXNAME" %in% colnames(d)) {
        as.character(d$TXNAME)
    } else {
        as.character(d$TXID)
    }
    transcript[is.na(transcript) | !nzchar(transcript)] <-
        gene[is.na(transcript) | !nzchar(transcript)]

    exons <- data.frame(
        chrom = as.character(d$TXCHROM),
        start = as.numeric(d$EXONSTART),
        end = as.numeric(d$EXONEND),
        strand = as.character(d$TXSTRAND),
        gene = gene,
        transcript = transcript,
        feature = "exon",
        gene_type = if ("TXTYPE" %in% colnames(d)) {
            as.character(d$TXTYPE)
        } else {
            NA_character_
        },
        stringsAsFactors = FALSE
    )

    if (all(c("CDSSTART", "CDSEND") %in% colnames(d))) {
        cdsKeep <- !is.na(d$CDSSTART) & !is.na(d$CDSEND)
        cds <- data.frame(
            chrom = as.character(d$TXCHROM)[cdsKeep],
            start = as.numeric(d$CDSSTART)[cdsKeep],
            end = as.numeric(d$CDSEND)[cdsKeep],
            strand = as.character(d$TXSTRAND)[cdsKeep],
            gene = gene[cdsKeep],
            transcript = transcript[cdsKeep],
            feature = "CDS",
            gene_type = exons$gene_type[cdsKeep],
            stringsAsFactors = FALSE
        )
        exons <- rbind(exons, cds)
    }
    .neoGeneNormalize(exons, regions)
}
