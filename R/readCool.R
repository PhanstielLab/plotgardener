#' Check for .(m)cool file and contents
#' @author Sarah Parker
#'
#' @param file Path to .(m)cool file
#'
#' @return A character string, either \code{".cool"} or \code{".mcool"},
#'   indicating the file type. Aborts with an error if the file is not a
#'   valid .(m)cool file.
#'
#' @importFrom glue glue glue_collapse
#' @importFrom rlang abort
#' @importFrom rhdf5 H5Fis_hdf5 h5ls
#' @importFrom dplyr pull
.checkCool <- function(file){
    
    ## Check if file is hdf5
    isHDF5 <- H5Fis_hdf5(file)
    
    if(!isHDF5){
        abort(c("File must be a `.cool` or `.mcool` file",
                "x" = glue("{file} is not in HDF5 format")))
    }
    
    lvl1contents <- h5ls(file, recursive = 1)
    expectedContents <- c("bins", "chroms", "indexes", "pixels")
    
    ## Check if it is mcool or cool
    if (nrow(lvl1contents) == 1){
        ## mcool file has 1 row-resolutions    
        
        ## Get deeper contents of mcool
        lvl2contents <- h5ls(file, recursive = 3) |>
            pull(name) |>
            unique()
        
        ## Check for expected datasets and resolutions
        
        if (all(expectedContents %in% lvl2contents) & 
           lvl1contents$name == "resolutions"){
            return(".mcool")
        } 
        
    } else {
        ## check that cool file contains expected datasets
        if (all(expectedContents %in% lvl1contents$name)){
            return(".cool")
        } 
    }
    
    ## Abort with message if not a cool file and doesn't contain expected contents
    expectedNames <- glue_collapse(glue("`{expectedContents}`"), sep = ", ")
    abort(c("File must be a `.cool` or `.mcool` file",
            "x" = glue("{file} does not contain expected datasets"),
            "i" = glue("Expected datasets {expectedNames}")))
    
}

#' Read basepair resolutions from an .(m)cool file
#' @author Sarah Parker
#'
#' @param file A character value specifying the path to the .(m)cool file
#' @importFrom rhdf5 h5ls h5read
#'
#' @returns Vector of basepair resolutions
#' 
#' @export
readCoolBpResolutions <- function(file){
    
    ## Check if file is `.mcool` or `.cool`
    fileType <- .checkCool(file)
    
    ## find all valid resolutions
    if(fileType == ".mcool"){
        ## get contents from cool file
        contents <- h5ls(file, recursive = 2)
        
        res <- contents[contents$group == "/resolutions","name"] |>
            as.integer() |>
            sort()
    } else {
        ## .cool files only have one resolution, calculate from bins
        start <- h5read(file, name = "/bins/start", index = list(1))
        end <- h5read(file, name = "/bins/end", index = list(1))
        res <- as.integer(end - start)
    }
    
    return(res)
}


#' Read normalizations included in .(m)cool files
#' @author Sarah Parker, Nicole Kramer
#'
#' @details
#' The "BALANCE" normalization refers to applying the pre-calculated matrix 
#' balancing weights in the `weight` dataset of `file`, typically present
#' in files created using cooler. VC is vanilla coverage, 
#' VC_SQRT is square root of vanilla coverage, 
#' and KR is Knight-Ruiz normalization.
#' 
#' Please note that if using a file from HiC-Pro, ICE normalizations will come
#' from files stored in a separate folder, and thus will not contain any 
#' normalizations explicitly called "ICE" or "BALANCE" since values are 
#' included already normalized. This normalization can be specified as "NONE".
#' 
#' @param file A character value specifying the path to the .(m)cool file
#' @param resolution optional, specify which resolution(s) to read 
#' normalization types from. Default is all resolutions in `file`.
#' 
#' @importFrom rlang arg_match
#' @importFrom rhdf5 h5ls
#'
#' @returns A vector or list of vectors of available normalizations
#' 
#' @export
readCoolNorms <- function(file, resolution = NULL){
    
    readBpNorm <- function(resolutionPath, file){
        bpNorm <- h5read(file, name = paste0(resolutionPath, "/bins")) |>
            names()
        bpNorm <- bpNorm[!(bpNorm %in% c("chrom", "start", "end"))] |>
            gsub("weight","BALANCE", x = _) |> # replace "weight" with "BALANCE"
            c("NONE") # add "NONE" for raw counts
        
        return(bpNorm)
    }
    
    fileType <- .checkCool(file)
    
    if (is.null(resolution)){
        if (fileType == ".mcool"){
            ## Get all resolutions from mcool file
            resolutions <- readCoolBpResolutions(file)
            
            ## Set dataset paths for each resolution
            datasetPaths <- paste0("/resolutions/", as.integer(resolutions))
            
            ## Get normalizations from the names of datasets in bins
            fileNorms <- lapply(datasetPaths, readBpNorm, file)
            names(fileNorms) <- resolutions
            
        } else {
            
            ## Get normalizations just from /bins
            fileNorms <- readBpNorm("", file)
        }
        
    } else {
        
        if (!resolution %in% readCoolBpResolutions(file)) {
            abort(c(
                glue("resolution={resolution} is not found in {file}."),
                "i"= glue("Use `readCoolBpResolutions({file})` to \\
                          see available resolutions.")
            ))
        }
        
        datasetPath <- ""
        ## if `.mcool` file, add the proper resolution to path
        if (fileType == ".mcool"){
            datasetPath <- paste0("/resolutions/", as.integer(resolution))
        }
        
        fileNorms <- readBpNorm(datasetPath, file)
        
    }
    
    return(fileNorms)
}


#' Read chromosomes included in .(m)cool files
#' @author Sarah Parker, Nicole Kramer
#' @param file A character value specifying the path to the .(m)cool file
#' @param resolution optional, specify which resolution(s) to read 
#' chromosomes from. Default is all resolutions in `file`.
#' 
#' @importFrom rhdf5 h5read
#'
#' @returns Data frame or list of data frames of chromosome names and lengths
#'
#' @export
readCoolChroms <- function(file, resolution = NULL){
    
    readBpChrom <- function(resolutionPath, file){
        bpChrom <- h5read(file, name = paste0(resolutionPath, "/chroms")) |>
            as.data.frame()
        
        ## add index column of the order of chromosomes
        bpChrom$index <- 1:nrow(bpChrom)
        return(bpChrom)
    }
    
    ## check if input is a cool file
    fileType <- .checkCool(file)
    
    if (is.null(resolution)){
        if (fileType == ".mcool"){
            
            ## Get all resolutions from mcool file
            resolutions <- readCoolBpResolutions(file)
            
            ## Set dataset paths for each resolution
            datasetPaths <- paste0("/resolutions/", as.integer(resolutions))
            
            ## Get chroms for each resolution
            chromInfo <- lapply(datasetPaths, readBpChrom, file)
            names(chromInfo) <- resolutions
            
        } else {
            chromInfo <- readBpChrom("", file)
        }
    } else {
        
        if (!resolution %in% readCoolBpResolutions(file)) {
            abort(c(
                glue("resolution={resolution} is not found in {file}."),
                "i"= glue("Use `readCoolBpResolutions({file})` to \\
                          see available resolutions.")
            ))
        }

        datasetPath <- ""
        
        ## if `.mcool` file, add the proper resolution to path
        if (fileType == ".mcool"){
            datasetPath <- paste0("/resolutions/", as.integer(resolution))
        }
    
        chromInfo <- readBpChrom(datasetPath, file)
        
    }
    
    return(chromInfo)
}


#' Determine best resolution for size of region for .(m)cool files
#' 
#' @param file Path to .(m)cool file
#' @param chromstart Chromstart of region
#' @param chromend Chromend of region
#'
#' @return A numeric value specifying the automatically determined resolution
#'   in base pairs.
.coolAutoResolution <- function(file, chromstart, chromend){
    
    fileResolutions <- readCoolBpResolutions(file)
    
    if (is.null(chromstart) & is.null(chromend)){
        autoRes <- max(fileResolutions)
    } else {
        dataRange <- chromend - chromstart
        if (dataRange >= 150000000) {
            autoRes <- max(fileResolutions)
        } else if (dataRange >= 75000000 & dataRange < 150000000) {
            autoRes <- 250000
            autoRes <- fileResolutions[which(
                abs(fileResolutions - autoRes) == min(
                    abs(fileResolutions - autoRes)
                )
            )]
        } else if (dataRange >= 35000000 & dataRange < 75000000) {
            autoRes <- 100000
            autoRes <- fileResolutions[which(
                abs(fileResolutions - autoRes) == min(
                    abs(fileResolutions - autoRes)
                )
            )]
        } else if (dataRange >= 20000000 & dataRange < 35000000) {
            autoRes <- 50000
            autoRes <- fileResolutions[which(
                abs(fileResolutions - autoRes) == min(
                    abs(fileResolutions - autoRes)
                )
            )]
        } else if (dataRange >= 5000000 & dataRange < 20000000) {
            autoRes <- 25000
            autoRes <- fileResolutions[which(
                abs(fileResolutions - autoRes) == min(
                    abs(fileResolutions - autoRes)
                )
            )]
        } else if (dataRange >= 3000000 & dataRange < 5000000) {
            autoRes <- 10000
            autoRes <- fileResolutions[which(
                abs(fileResolutions - autoRes) == min(
                    abs(fileResolutions - autoRes)
                )
            )]
        } else {
            autoRes <- 5000
            autoRes <- fileResolutions[which(
                abs(fileResolutions - autoRes) == min(
                    abs(fileResolutions - autoRes)
                )
            )]
        }
    }
    
    return(as.integer(autoRes))
}

#' Add (alt)chromstart and (alt)chromend for NULL (alt)chrom region of 
#' .(m)cool files
#' 
#' @param file Path to .(m)cool file
#' @param chrom Chromosome of region; can also be altchromosome
#' @param resolution Resolution to read chromsome info from
#'
#' @return A list of length two containing the chromosome start position (1)
#'   and the chromosome length.
.coolRegion <- function(file, chrom, resolution){
    
    chromInfo <- readCoolChroms(file, resolution = resolution)
    chromLength <- chromInfo[chromInfo$name == chrom, "length"]
    
    return(list(1, chromLength))
}

#' Read in data for a bin chunk
#' 
#' @param binChunk The binChunk indeces to read
#' @param file Path to .(m)cool file
#' @param bin_offsets Read in bin1 offsets
#' @param binChunkSize Size of bin chunk, for comparison against the end of 
#' the bin chunk
#' @param datasetPath Dataset path, for specifying resolution in .mcool file
#' @param end1bin Bin where end1 starts
#' @param start2bin Bin for chr2 starts
#' @param end2bin Bin for end2 starts
#'
#' @return An integer vector of pixel indices corresponding to interactions
#'   within the specified bin chunk, or \code{NA} if none are found.
#'
#' @importFrom rhdf5 h5read
.pullBinChunks <- function(binChunk, file, bin_offsets, binChunkSize,
                           datasetPath, end1bin,
                           start2bin, end2bin){
    
    binChunkEnd <- binChunk + binChunkSize - 1
    
    if (binChunkEnd > bin_offsets[end1bin+1]+1){
        binChunkEnd <- bin_offsets[end1bin+1]+1
    }
    if (binChunkEnd > binChunk){
        bin2s <- h5read(file, 
                        name = paste0(datasetPath,"/pixels/bin2_id"),
                        index = list(binChunk:binChunkEnd)) 
        
        ## Find indexes for interactions with all bin2s between start2 and end2
        count_idx <- which(bin2s %in% start2bin:end2bin) + binChunk - 1
    } else {
        count_idx <- NA
    }
    
    return(count_idx)
}


#' Error checking function for .(m)cool files
#' 
#' @param file Path to .(m)cool file
#' @param chrom User-inputted chromosome
#' @param chromstart User-inputted chromstart, can still be NULL at this point.
#' @param chromend User-inputted chromend, can still be NULL at this point.
#' @param zrange User-inputted zrange.
#' @param altchrom User-inputted alt chromosome.
#' @param altchromstart User-inputted alt chromstart.
#' @param altchromend User-inputted alt chromend.
#' @param norm User-inputted normalization.
#' @param resolution Resolution, either user-inputted or determined by 'auto'.
#'
#' @return Called for side effects (input validation). Aborts with an
#'   informative error if any inputs are invalid.
#'
#' @importFrom glue glue
#' @importFrom rlang abort
.checkCoolErrors <- function(file, chrom, chromstart, chromend, zrange,
                       altchrom, altchromstart, altchromend, norm, resolution){
    ## File input
    .checkCool(file)
    
    ## Norm
    fileNorms <- readCoolNorms(file, resolution = resolution)
    if (!norm %in% fileNorms){
        abort(c(
            glue("norm={norm} is not found in {file}."),
            "i" = glue("Use `readCoolNorms({file})` to see available norms.")
        ))
    }
    
    ## Chroms
    chromInfo <- readCoolChroms(file, resolution = resolution)
    
    if (!chrom %in% chromInfo$name){
        abort(c(
            glue("chrom={chrom} not found in {file}."),
            "i" = glue("Check file chromosome names with \\
                       `readCoolChroms({file}).")
        ))
    }
    
    ## Genomic region
    regionErrors(chromstart = chromstart, chromend = chromend)
    
    ## Check that chromstart is in bounds
    if (!is.null(chromstart)){
        chromLength <- chromInfo[chromInfo$name == chrom, "length"]
        if (chromstart < 0 | chromstart > chromLength){
            abort(c(glue("chromstart must be between 0 and chromosome end."),
                    "i" = "{chrom} length is {chromLength} bp."))
        }
    }
    
    if (!is.null(altchrom)){
        if (!altchrom %in% chromInfo$name){
            abort(c(
                glue("altchrom={altchrom} not found in {file}."),
                "i" = glue("Check file chromosome names with \\
                           `readCoolChroms({file}).")
            ))
        }
        
        ## Can't specify altchrom without a chrom
        if (is.null(chrom)){
            abort(c("Specified altchrom={altchrom} but did not give chrom."))
        }
        
        regionErrors(chromstart = altchromstart,
                     chromend = altchromend)
        
        altchromLength <- chromInfo[chromInfo$name == altchrom, "length"]
        if (altchromstart < 0 | altchromstart > altchromLength){
            abort(c(glue("altchromstart must be between 0 and alt chromosome end."),
                    "i" = "{altchrom} length is {altchromLength} bp."))
        }
        
        ## If giving same chrom and altchrom, need to specify
        ## chromstart/chromend and altchromstart/altchromend
        
        if (chrom == altchrom) {
            if (is.null(chromstart) |
                is.null(chromend) |
                is.null(altchromstart) |
                is.null(altchromend)) {
                
                abort(c(glue("No chromstart, chromend, altchromstart, and \\
                             altchromend given for same chrom and altchrom."),
                        "i" = glue("If trying to get all interactions between \\
                                   one chromsome, just specify chrom.")))
                
            }
        }
        
    }
    
    ## zrange errors
    rangeErrors(range = zrange)
    
}

#' Read a .(m)cool file and return Hi-C data as a dataframe
#' @author Sarah Parker, Nicole Kramer
#' @usage readCool(
#'     file,
#'     chrom,
#'     chromstart = NULL,
#'     chromend = NULL,
#'     altchrom = NULL,
#'     altchromstart = NULL,
#'     altchromend = NULL,
#'     resolution = "auto",
#'     zrange = NULL,
#'     norm = "NONE",
#'     binChunkSize = 5e6,
#'     params = NULL,
#'     quiet = FALSE
#' )
#' 
#' @param file A character value specifying the path to the .(m)cool file.
#' @param chrom Chromosome of data, as a string.
#' @param chromstart Integer start position on chromosome.
#' @param chromend Integer end position on chromosome.
#' @param altchrom Alternate chromosome for interchromosomal data,
#' as a string.
#' @param altchromstart Alternate chromosome integer start position
#' for interchromosomal data.
#' @param altchromend Alternate chromosome integer end position
#' for interchromosomal data.
#' @param resolution A numeric specifying the width of each pixel.
#' "auto" will attempt to choose a resolution in basepairs based on
#' the size of the region.
#' @param zrange A numeric vector of length 2 specifying the range of
#' interaction scores, where extreme values will be set to the max or min.
#' @param norm Character value specifying hic data normalization method.
#' This value must be found in the .(m)cool file.
#' Default value is \code{norm = "NONE"}.
#' @param binChunkSize A numeric specifying the number of bin indices to read 
#' from a file for a given region at a given resolution. If the total amount of 
#' data is larger than the \code{binChunkSize}, data will be read in multiple
#' chunks. Default value is \code{binChunkSize = 5e6}.
#' @param params An optional \link[plotgardener]{pgParams} object
#' containing relevant function parameters.
#' @param quiet A logical indicating whether or not to print messages.
#' 
#' 
#' @return Returns a 3-column dataframe in sparse upper triangular
#' format with the following columns: \code{chrom}, \code{altchrom},
#' \code{counts}.
#' 
#' @examples
#' 
#' ## .cool file
#' coolFile <- file.path(tempdir(), "Rao2014-IMR90-MboI-allreps-filtered.1000kb.cool")
#' download.file(url = "https://usgs2.osn.mghpcc.org/cooler01/examples/hg19/Rao2014-IMR90-MboI-allreps-filtered.1000kb.cool",
#'     destfile = coolFile, mode = "wb")
#' 
#' ## Read in region `chr2:10000000-22000000` at 1000Kb cool file resolution
#' coolData <- readCool(file = coolFile, chrom = "chr2", chromstart = 10000000,
#'                      chromend = 22000000,
#'                      resolution = 1000000)
#' 
#' ## .mcool file
#' mcoolFile <- file.path(tempdir(), "LEUK_HEK_PJA27_inter_30.mcool")
#' download.file(url = "https://zenodo.org/records/10906240/files/LEUK_HEK_PJA27_inter_30.mcool?download=1",
#'     destfile = mcoolFile, mode = "wb")
#'
#' ## Read in region `chr2:1000000-5000000` at 100Kb resolution
#' mcoolData_100Kb <- readCool(file = mcoolFile, chrom = "2",
#'                             chromstart = 1000000, chromend = 5000000,
#'                             resolution = 100000)
#' 
#' ## Read in data for chr2 at 2500Kb resolution 
#' mcoolData_2500Kb <- readCool(file = mcoolFile, chrom = "2",
#'                              resolution = 2500000)
#' @seealso \link[plotgardener]{readHic}
#'
#' @importFrom rlang inform warn
#' @export
readCool <- function(file, chrom, chromstart = NULL, chromend = NULL, 
                     altchrom = NULL, altchromstart = NULL, altchromend = NULL, 
                     resolution = "auto", zrange = NULL, norm = "NONE",
                     binChunkSize = 5e6, params = NULL, quiet = FALSE){
    
    # =========================================================================
    # PARSE PARAMETERS
    # =========================================================================
    
    rcool <- parseParams(params = params, 
                        defaultArgs = formals(eval(match.call()[[1]])),
                        declaredArgs = lapply(match.call()[-1], 
                                              eval.parent, n = 2),
                        class = "rcool")
    
    if (is.null(rcool$file)){
        abort(c("argument \"file\" is missing, with no default."))}

    if (is.null(rcool$chrom)){
        abort(c("argument \"chrom\" is missing, with no default."))}
    
    # =========================================================================
    # ADJUST RESOLUTION
    # =========================================================================

    if (rcool$resolution == "auto") {
        rcool$resolution <- .coolAutoResolution(
            file = rcool$file,
            chromstart = rcool$chromstart,
            chromend = rcool$chromend
        )
    } 

    # =========================================================================
    # CATCH ERRORS
    # =========================================================================

    .checkCoolErrors(
        file = rcool$file, chrom = rcool$chrom,
        chromstart = rcool$chromstart,
        chromend = rcool$chromend, zrange = rcool$zrange,
        altchrom = rcool$altchrom,
        altchromstart = rcool$altchromstart,
        altchromend = rcool$altchromend, norm = rcool$norm,
        resolution = rcool$resolution
    )

    # =========================================================================
    # SET REGION PARAMETERS
    # =========================================================================
    
    if (is.null(rcool$chromstart) & is.null(rcool$chromend)){
        chrRegion <- .coolRegion(file = rcool$file,
                                 chrom = rcool$chrom,
                                 resolution = rcool$resolution)
        
        rcool$chromstart <- chrRegion[[1]]
        rcool$chromend <- chrRegion[[2]]
    }

    
    ## For off diagonal plotting, grabbing whole symmetric region
    if (!is.null(rcool$altchrom)) {
        
        if (is.null(rcool$altchromstart) & is.null(rcool$altchromend)){
            altchrRegion <- .coolRegion(file = rcool$file,
                                        chrom = rcool$altchrom,
                                        resolution = rcool$resolution)
            rcool$altchromstart <- altchrRegion[[1]]
            rcool$altchromend <- altchrRegion[[2]]
        }
        
        
        if (rcool$chrom == rcool$altchrom) {
            rcool$chromstart <- min(rcool$chromstart, rcool$altchromstart)
            rcool$chromend <- max(rcool$chromend, rcool$altchromend)
        }
    } else {
        rcool$altchrom <- rcool$chrom
        rcool$altchromstart <- rcool$chromstart
        rcool$altchromend <- rcool$chromend
        
    }
    # =========================================================================
    # EXTRACT SPARSE UPPER TRIANGULAR COUNT MATRIX
    # =========================================================================

    fileType <- .checkCool(rcool$file)
    chrInfo <- readCoolChroms(rcool$file, rcool$resolution)
    datasetPath <- ""
    if (fileType == ".mcool"){
        datasetPath <- paste0("/resolutions/", as.integer(rcool$resolution))
    }

    ## Find bin IDs for both chrom locations -----------------------------------
    ## adjust starts and ends to start of bins
    
    start1 <- (rcool$chromstart %/% rcool$resolution) * rcool$resolution
    start2 <- (rcool$altchromstart %/% rcool$resolution) * rcool$resolution

    end1 <- (rcool$chromend %/% rcool$resolution) * rcool$resolution
    end2 <- (rcool$altchromend %/% rcool$resolution) * rcool$resolution
    
    ## Get vector of chromosome offsets 
    chrom_offsets <- h5read(rcool$file, 
                            name = paste0(datasetPath, "/indexes/chrom_offset"))
    
    ## find bin id for chr1:start1:end1
    chr1indx <- chrInfo[chrInfo$name == rcool$chrom, "index"]
    
    ## get start1 bin
    chr1_starts <- h5read(rcool$file, name = paste0(datasetPath, "/bins/start"),
                          index = list((chrom_offsets[chr1indx]+1):
                                           chrom_offsets[chr1indx+1]))
    start1bin <- which(chr1_starts == start1) + chrom_offsets[chr1indx] - 1
    
    ## get end1 bin (i.e. bin where end1 starts)
    end1bin <- which(chr1_starts == end1) + chrom_offsets[chr1indx] - 1
    
    ## find bin id for chr2:start2:end2
    chr2indx <- chrInfo[chrInfo$name == rcool$altchrom, "index"]
    
    ## get start2 bin
    chr2_starts <- h5read(rcool$file, name = paste0(datasetPath,"/bins/start"),
                          index = list((chrom_offsets[chr2indx]+1):
                                           chrom_offsets[chr2indx+1]))
    start2bin <- which(chr2_starts == start2) + chrom_offsets[chr2indx] - 1

    ## get end2 bin
    end2bin <- which(chr2_starts == end2) + chrom_offsets[chr2indx] - 1
    
    ## Extract counts and match locations --------------------------------------
    
    ## Get bin offsets for counts
    bin_offsets <- h5read(rcool$file, 
                          name = paste0(datasetPath,"/indexes/bin1_offset"))
                          
    ## Pull all bin2 ids for interactions with all bin1s between start1 and end1
    ## Read 5 million indices at a time for more efficient calling
    allBin1s <- (bin_offsets[start1bin+1]+1):(bin_offsets[end1bin+1]+1)
    
    if (length(allBin1s) > rcool$binChunkSize){
        binChunks <- seq(bin_offsets[start1bin+1]+1,
                         bin_offsets[end1bin+1]+1,
                         binChunkSize)
        
        
        count_idx <- na.omit(unlist(lapply(binChunks, .pullBinChunks, 
                                           file = rcool$file,
               bin_offsets = bin_offsets, binChunkSize = rcool$binChunkSize,
               datasetPath = datasetPath, end1bin = end1bin,
               start2bin = start2bin, end2bin = end2bin)))
        
    } else {
        bin2s <- h5read(rcool$file, 
                        name = paste0(datasetPath,"/pixels/bin2_id"),
                        index = list(allBin1s)) 
        
        ## Find indexes for interactions with all bin2s between start2 and end2
        count_idx <- which(bin2s %in% start2bin:end2bin) + 
            bin_offsets[start1bin+1]
    }
    
    ## Pull out bin1 ids, bin2 ids, and counts for all interactions in slice
    bin1ids <- h5read(rcool$file, name = paste0(datasetPath,"/pixels/bin1_id"), 
                      index = list(count_idx))
    
    bin2ids <- h5read(rcool$file, name = paste0(datasetPath,"/pixels/bin2_id"), 
                      index = list(count_idx))
    
    counts <- h5read(rcool$file, name = paste0(datasetPath,"/pixels/count"), 
                     index = list(count_idx)) |>
        #defaults to integer, needs to be numeric to convert to NA correctly
        as.numeric()
    
    ## Multiply by normalization factors, if applicable
    if (rcool$norm == "BALANCE"){
        bin1norm <- h5read(rcool$file, 
                           name = paste0(datasetPath,"/bins/weight"),
                           list(as.numeric(bin1ids+1)))
        bin2norm <- h5read(rcool$file, 
                           name = paste0(datasetPath,"/bins/weight"),
                           list(as.numeric(bin2ids+1)))
        counts <- counts / (bin1norm * bin2norm)
    } else if (rcool$norm != "NONE"){
        bin1norm <- h5read(rcool$file, 
                           name = paste0(datasetPath,"/bins/", rcool$norm),
                           list(as.numeric(bin1ids+1)))
        bin2norm <- h5read(rcool$file, 
                           name = paste0(datasetPath,"/bins/", rcool$norm),
                           list(as.numeric(bin2ids+1)))
        counts <- counts / (bin1norm * bin2norm)
    }
    
    ## Get genomic locations for each bin id 
    xs <- h5read(rcool$file, 
                 name = paste0(datasetPath,"/bins/start"),
                 index = list(as.numeric(bin1ids+1)))
    ys <- h5read(rcool$file, 
                 name = paste0(datasetPath,"/bins/start"),
                 index = list(as.numeric(bin2ids+1)))
    
    ## create sparse upper triangular matrix
    upper <- data.frame(x = xs, y = ys, counts = counts)
    
    # =========================================================================
    # SCALE DATA WITH ZRANGE
    # =========================================================================
    upper <- scale_data(upper = upper, zrange = rcool$zrange)
    
    # =========================================================================
    # FORMAT DATA IN PROPER ORDER AND WITH LABELS
    # =========================================================================
    
    ## Rename columns based on chrom/altchrom
    if (rcool$chrom == rcool$altchrom){
        colnames(upper) <- c(paste0(rcool$chrom, "_A"),
                             paste0(rcool$altchrom, "_B"),
                             "counts")
    } else {
        colnames(upper) <- c(rcool$chrom, rcool$altchrom, "counts")
    }
    
    # =========================================================================
    # REMOVE NAN VALUES
    # =========================================================================
    
    upper <- na.omit(upper)
    
    # =========================================================================
    # RETURN DATAFRAME
    # =========================================================================
    if (nrow(upper) == 0) {
        warn(c(glue("No data found in region.")))
    } else {
        if (!rcool$quiet) {
            inform(c(glue("Read in {fileType} file with {rcool$norm} \\
                          normalization at {rcool$resolution} BP resolution.")))
        }
    }
    return(upper)
    
    
}