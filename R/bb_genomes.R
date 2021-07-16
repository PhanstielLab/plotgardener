#' Display the included available default genome assemblies
#'
#' @return Returns the included available default genome assemblies
#'
#' @examples
#' bb_genomes()
#' @export
bb_genomes <- function() {
    
    availDefaults <- default_genomePackages$Genome

    return(cat(availDefaults, sep = "\n"))
}
