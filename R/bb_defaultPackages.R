#' Displays the default data packages associated with a genome build
#'
#' @param genome string indicating the genome assembly
#' @export
bb_defaultPackages <- function(genome){


  availDefaults <- c("bosTau8", "bosTau9", "canFam3", "ce6", "ce11", "danRer10",
                     "danRer11", "dm3", "dm6", "galGal4", "galGal5", "galGal6",
                     "hg18", "hg19", "hg38", "mm9", "mm10", "rheMac3", "rheMac8",
                     "rehMac10", "panTro5", "panTro6", "rn4", "rn5", "rn6", "sacCer2",
                     "sacCer3", "susScr3", "susScr11")

  if (!genome %in% availDefaults){
    stop("Inputted genome not an available default. To see the included defaults, use `bb_genomes()`.", call. = FALSE)
  }

  defaults <- getPackages(genome = genome)
  return(str(defaults))

}
