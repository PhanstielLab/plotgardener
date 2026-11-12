# Define a function to read in various kinds of genomic range data
#' @importFrom plyranges %>%
#' @importFrom GenomicRanges GRanges
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

                data <- readBigwig(
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
        if (!all(c("chrom", "start", "end") %in% colnames(data))) {
            colnames(data)[seq(1, 3)] <- c("chrom", "start", "end")
        }
        ## Put the expected columns first without changing other columns
        data <- data[, c("chrom", "start", "end",
                         setdiff(colnames(data),
                                 c("chrom", "start", "end"))), drop = FALSE]

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
        
        
    } else if (type == "GWAS") {
        
        ## Rename chrom and position columns if necessary
        if (!all(c("chrom", "pos") %in% colnames(data))) {
            colnames(data)[seq(1, 2)] <- c("chrom", "pos")
        }
        data <- data[, c("chrom", "pos",
                         setdiff(colnames(data),
                                 c("chrom", "pos"))), drop = FALSE]
        ## column data types will be checked in plotManhattan errorchecks
        
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
    
    ## Check anchor orders and try to reorder if necessary
    anchor1 <- GRanges(seqnames = data[,"chrom1"], 
                       ranges = IRanges(start = data[,"start1"], 
                                        end = data[,"end1"]))
    
    anchor2 <- GRanges(seqnames = data[,"chrom2"], 
                       ranges = IRanges(start = data[,"start2"], 
                                        end = data[,"end2"]))
    
    swap <- anchor2 < anchor1
    if (any(swap)){
        # Swap any anchors where anchor2 < anchor1
        temp <- anchor1[swap]
        anchor1[swap] <- anchor2[swap]
        anchor2[swap] <- temp
        # Convert back to 6 column bedpe format
        swapped_data <- cbind(as.data.frame(anchor1)[, c("seqnames", "start", "end")],
                              as.data.frame(anchor2)[,c("seqnames", "start", "end")])
        colnames(swapped_data) <- c("chrom1", "start1", "end1", 
                                    "chrom2", "start2", "end2")
        # Join back with any trailing data columns
        data <- cbind(swapped_data, data[,-c(1:6)])
        warning("One or more pairs contained unordered anchors and were",
                " swapped for plotting.", call. = FALSE)
    }    
    
    if (warning == TRUE){
        if (nrow(data) < 1){
            warning("\'data\' input contains no values.", call. = FALSE)
            }
    }
    return(data)
    }
    
# Define a function that checks matching for GRanges/GInteractions
# objects and declared assembly
checkAssemblyMatch <- function(data, assembly){
    
    genome <- unique(GenomeInfoDb::genome(data))
    if (!is.na(genome)){
        if (genome != assembly$Genome){
            warning("Input data assembly detected as ",
                    genome, " and plotgardener assembly ",
                    "detected as ", assembly$Genome, ".",
                    .call = FALSE)
        }
    }    
}