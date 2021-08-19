#' Display the included available default genome assemblies
#'
#' @return Returns the included available default genome assemblies
#'
#' @examples
#' genomes()
#' @export
genomes <- function() {
    
    availDefaults <- default_genomePackages$Genome

    return(cat(availDefaults, sep = "\n"))
}
