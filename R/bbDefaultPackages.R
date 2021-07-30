#' Display the default genomic annotation packages associated with a
#' genome build
#'
#' @param Genome String indicating the name of the genome assembly.
#'
#' @return Returns a list of the default data packages for a genome build.
#'
#' @examples
#' ## View default genomic annotation packages associated with "hg19"
#' bbDefaultPackages(Genome = "hg19")
#'
#' ## View default genomic annotation packages associated with "mm9"
#' bbDefaultPackages(Genome = "mm9")
#' @export
bbDefaultPackages <- function(Genome) {

    availDefaults <- default_genomePackages$Genome

    if (!Genome %in% availDefaults) {
        stop("Inputted genome not an available default. To see the included ",
            "defaults, use `bbGenomes()`.", call. = FALSE)
    }

    defaults <- default_genomePackages[
        default_genomePackages$Genome == Genome,c("Genome",
                                                "TxDb", 
                                                "OrgDb",
                                                "gene.id.column",
                                                "display.column",
                                                "BSgenome")]
    
    return(str(defaults))
}
