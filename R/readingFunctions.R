# Define a function to read in various kinds of genomic range data
#' @importFrom plyranges %>%
read_rangeData <- function(data, assembly, chrom = NULL, 
                        start = NULL, end = NULL, type = "range"){
    
    if (!"data.frame" %in% class(data)) {
        if (!is(data, "GRanges")) {
            if (file_ext(data) == "bam") {
                indexFile <- paste0(data, ".bai")
                if (!file.exists(indexFile)) {
                    stop("Cannot read in bam file without a ",
                        "corresponding bam index file (.bai) in the ",
                        "same directory.", call. = FALSE)
                }
                data <- plyranges::read_bam(data) %>%
                    plyranges::filter_by_overlaps(GenomicRanges::GRanges(
                        seqnames = chrom,
                        ranges = IRanges::IRanges(
                            start = start,
                            end = end
                        )
                    )) %>%
                    dplyr::mutate()
            } else if (file_ext(data) %in% c("bw", "bigWig",
                                            "bigwig", "bedgraph")) {

                data <- bbReadBigwig(
                    file = data,
                    chrom = chrom,
                    chromstart = start,
                    chromend = end
                )
            } else {
                data <- data.table::fread(data)
            }
        } else {
            
            ## check GRanges genome with assembly input
            checkAssemblyMatch(data = data, assembly = assembly)
            
        }
    }
    
    data <- as.data.frame(data)
    
    if (type == "range"){
        ## Rename columns if necessary
        colnames(data)[seq(1, 3)] <- c("chrom", "start", "end")
        
        ## Check column input types
        columnClasses <- list("start" = TRUE, "end" = TRUE)
        if (!is(data[, "start"], "integer") &
            !is(data[, "start"], "numeric")){
            columnClasses["start"] <- FALSE
        }
        if (!is(data[, "end"], "integer") &
            !is(data[, "end"], "numeric")){
            columnClasses["end"] <- FALSE
        }
        
        if (any(columnClasses == FALSE)){
            
            ## Get the ones that are wrong
            wrongClasses <- names(which(columnClasses == FALSE))
            stop(cat(wrongClasses, sep = ", "), 
            " are not the correct input type. ",
            "Chromosome column must be a character and chromstart and",
            " chromend columns must be integers or numerics.", call. = FALSE)
        }
        
        
    }
    
    return(data)
    
}

# Define a function to read in various kinds of genomic paired range data
read_pairedData <- function(data, assembly, warning = FALSE){
    
    if (!"data.frame" %in% class(data)) {
        if (!is(data, "GInteractions")) {
            data <- as.data.frame(data.table::fread(data))
        } else {
            ## check GInteractions genome with assembly input
            checkAssemblyMatch(data = data, assembly = assembly)
            
            ## Reorder GInteractions columns
            data <- as.data.frame(data)
            dataSubset <- data[, c(
                "seqnames1", "start1", "end1",
                "seqnames2", "start2", "end2"
            )]
            
            data <- data[, which(!colnames(data) %in%
                                colnames(dataSubset))]
            data <- cbind(dataSubset, data)
        }
    } else {
        data <- as.data.frame(data)
    }
    
    ## Rename columns
    colnames(data)[seq(1, 6)] <- c("chrom1", "start1", "end1",
                                "chrom2", "start2", "end2")
    
    ## Check column input types
    columnClasses <- list("start1" = TRUE, "end1" = TRUE,
                        "start2" = TRUE, "end2" = TRUE)
    if (!is(data[, "start1"], "integer") &
        !is(data[, "start1"], "numeric")){
        columnClasses["start1"] <- FALSE
    }
    if (!is(data[, "start2"], "integer") &
        !is(data[, "start2"], "numeric")){
        columnClasses["start2"] <- FALSE
    }
    if (!is(data[, "end1"], "integer") &
        !is(data[, "end1"], "numeric")){
        columnClasses["end1"] <- FALSE
    }
    if (!is(data[, "end2"], "integer") &
        !is(data[, "end2"], "numeric")){
        columnClasses["end2"] <- FALSE
    }
    
    if (any(columnClasses == FALSE)){
        
        ## Get the ones that are wrong
        wrongClasses <- names(which(columnClasses == FALSE))
        stop(wrongClasses, " are not the correct input type. ",
            "Chromosome columns must be a character and chromstart and",
            " chromend columns must be integers or numerics.", call. = FALSE)
    }
    
    if (warning == TRUE){
        if (nrow(data) < 1){
            warning("\'data\' input contains no values.", call. = FALSE)
        }
        
    }
    
    return(data)
    
    }
    
# Define a function that checks matching for GRanges/GInteractions
# objects and declared BentoBox assembly
checkAssemblyMatch <- function(data, assembly){
    
    genome <- unique(GenomeInfoDb::genome(data))
    if (!is.na(genome)){
        if (genome != assembly$Genome){
            warning("Input data assembly detected as ",
                    genome, " and BentoBox assembly ",
                    "detected as ", assembly$Genome, ".",
                    .call = FALSE)
        }
    }    
}