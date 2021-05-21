#' Display the default genomic annotation packages associated with a
#' genome build
#'
#' @param Genome String indicating the name of the genome assembly.
#'
#' @return Returns a list of the default data packages for a genome build.
#'
#' @examples
#' ## View default genomic annotation packages associated with "hg19"
#' bb_defaultPackages(Genome = "hg19")
#'
#' ## View default genomic annotation packages associated with "mm9"
#' bb_defaultPackages(Genome = "mm9")
#' @export
bb_defaultPackages <- function(Genome) {
    availDefaults <- c(
        "bosTau8", "bosTau9", "canFam3", "ce6", "ce11",
        "danRer10", "danRer11", "dm3", "dm6", "galGal4",
        "galGal5", "galGal6", "hg18", "hg19", "hg38", "mm9",
        "mm10", "rheMac3", "rheMac8", "rehMac10", "panTro5",
        "panTro6", "rn4", "rn5", "rn6", "sacCer2",
        "sacCer3", "susScr3", "susScr11"
    )

    if (!Genome %in% availDefaults) {
        stop("Inputted genome not an available default. To see the included
            defaults, use `bb_genomes()`.", call. = FALSE)
    }

    defaults <- getPackages(genome = Genome)
    return(str(defaults))
}
