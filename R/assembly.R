#' Make a \code{assembly} object for alternate TxDb, OrgDb,and BSgenome
#' genomic annotation packages
#' 
#' @usage assembly(
#'     Genome,
#'     TxDb,
#'     OrgDb,
#'     gene.id.column = "ENTREZID",
#'     display.column = "SYMBOL",
#'     BSgenome = NULL
#' )
#'
#' @param Genome String indicating the name of the genome assembly.
#' @param TxDb String of existing TxDb package name, a TxDb object, or an
#' \link[ensembldb]{EnsDb} object. EnsDb objects are converted to a TxDb for
#' plotting; provide \code{OrgDb} and set \code{gene.id.column = "ENSEMBL"}
#' if gene symbols are desired.
#' @param OrgDb String of the desired OrgDb package name.
#' @param gene.id.column String of the TxDb column name that refers to
#' the given TxDb gene IDs.
#' Default value is \code{gene.id.column = "ENTREZID"}.
#' @param display.column String of the OrgDb column name that refers to
#' the type of gene symbol to be displayed in plots.
#' Default value is \code{display.column = "SYMBOL"}.
#' @param BSgenome String of the desired BSgenome package name.
#'
#' @return Returns a \code{assembly} object containing all input parameters.
#'
#' @examples
#' ## Create a custom `assembly` object for hg38/GRCh38 packages
#' newAssembly <- assembly(
#'     Genome = "hg38_GRCh38",
#'     TxDb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
#'     OrgDb = "org.Hs.eg.db",
#'     BSgenome = "BSgenome.Hsapiens.NCBI.GRCh38"
#' )
#' @seealso \link[GenomicFeatures]{TxDb}, \link[ensembldb]{EnsDb},
#' \link[AnnotationDbi]{OrgDb-class}, \link[BSgenome]{BSgenome}
#'
#' @export
assembly <- function(Genome, TxDb, OrgDb, gene.id.column = "ENTREZID",
                        display.column = "SYMBOL", BSgenome = NULL) {
    object <- structure(list(
        Genome = Genome, TxDb = TxDb, OrgDb = OrgDb,
        gene.id.column = gene.id.column,
        display.column = display.column
    ),
    class = "assembly"
    )

    if (is(object$TxDb, "TxDb")) {
        object$TxDb$packageName <- paste0("TxDb.", Genome)
    } else if (is(object$TxDb, "EnsDb")) {
        if (!requireNamespace("ensembldb", quietly = TRUE) ||
            !requireNamespace("GenomicFeatures", quietly = TRUE)) {
            stop("`ensembldb` and `GenomicFeatures` are required to use an ",
                 "EnsDb object in `assembly()`.", call. = FALSE)
        }
        object$EnsDb <- object$TxDb
        tx <- ensembldb::transcripts(object$TxDb)
        object$TxDb <- GenomicFeatures::makeTxDbFromGRanges(tx)
        object$TxDb$packageName <- paste0("EnsDb.", Genome)
        if (is.null(object$OrgDb) ||
            (is.character(object$OrgDb) && length(object$OrgDb) == 1 &&
             is.na(object$OrgDb))) {
            object$gene.id.column <- "GENEID"
            object$display.column <- "GENEID"
        } else if (identical(object$gene.id.column, "ENTREZID")) {
            object$gene.id.column <- "ENSEMBL"
        }
    }

    if (!is.null(BSgenome)) {
        object$BSgenome <- BSgenome
    }

    invisible(object)
}
