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
    rootNames <- unique(lvl1contents$name)
    
    ## Check if it is mcool or cool. An mcool file has a root-level
    ## `resolutions` group; it can contain additional root groups, so do not
    ## rely on the number of root rows.
    if ("resolutions" %in% rootNames){
        ## Get deeper contents of mcool
        lvl2contents <- h5ls(file, recursive = 3) |>
            pull(name) |>
            unique()
        
        ## Check for expected datasets and resolutions
        if (all(expectedContents %in% lvl2contents)){
            return(".mcool")
        }
    } else {
        ## check that cool file contains expected datasets
        if (all(expectedContents %in% rootNames)){
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
#' The `weight` column is the standard cooler balancing weight, typically
#' present in files created using `cooler balance`.
#' VC is vanilla coverage, VC_SQRT is square root of vanilla coverage,
#' and KR is Knight-Ruiz normalization. Other columns such as `sweight`,
#' `CNV`, or `obj_weight` are returned by their stored names. Use
#' `readCoolNorms()` to see the exact column names available in a file.
#' 
#' Cooler weights are multiplicative: balanced counts are calculated as
#' \code{observed * weight1 * weight2}. Contacts involving bins with
#' non-finite weights are dropped; the bins remain in the coordinate space.
#' 
#' Please note that if using a file from HiC-Pro, ICE normalizations will come
#' from files stored in a separate folder, and thus will not contain any 
#' normalizations explicitly called "ICE" or "weight" since values are 
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
        bpNorm <- bpNorm[!(bpNorm %in% c("chrom", "start", "end"))]
        bpNorm <- unique(c(bpNorm, "NONE")) # add "NONE" for raw counts
        
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


#' Choose an available normalization from a .(m)cool file
#'
#' @param file Path to .(m)cool file
#' @param resolution Resolution to inspect
#'
#' @return A character value naming an available normalization.
.coolAutoNorm <- function(file, resolution = NULL){
    fileNorms <- readCoolNorms(file, resolution = resolution)

    if ("KR" %in% fileNorms){
        return("KR")
    }
    if ("weight" %in% fileNorms){
        return("weight")
    }
    if ("sweight" %in% fileNorms){
        return("sweight")
    }
    available <- setdiff(fileNorms, "NONE")
    if (length(available) > 0){
        return(available[1])
    }
    return("NONE")
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
    chrom <- .coolMatchChrom(chrom, chromInfo)
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
#' @param endPixel Last pixel index to read (1-based, inclusive)
#' @param start2bin Bin for chr2 starts
#' @param end2bin Bin for end2 starts
#'
#' @return An integer vector of pixel indices corresponding to interactions
#'   within the specified bin chunk, or \code{NA} if none are found.
#'
#' @importFrom rhdf5 h5read
.pullBinChunks <- function(binChunk, file, bin_offsets, binChunkSize,
                           datasetPath, endPixel,
                           start2bin, end2bin){
    
    binChunkEnd <- min(binChunk + binChunkSize - 1, endPixel)
    
    if (binChunkEnd >= binChunk){
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


#' Match a user chromosome name to a chromosome name in a .(m)cool file
#'
#' @param chrom User-supplied chromosome name.
#' @param chromInfo Chromosome data frame from \code{readCoolChroms()}.
#'
#' @return The matching chromosome name, or \code{NULL} if no match is found.
.coolMatchChrom <- function(chrom, chromInfo){
    if (is.null(chrom) || length(chrom) == 0){
        return(NULL)
    }
    chrom <- as.character(chrom)
    if (chrom %in% chromInfo$name){
        return(chrom)
    }

    ## Accept `chr1` for files that use `1` and vice versa.
    altChrom <- if (startsWith(chrom, "chr")){
        sub("^chr", "", chrom)
    } else {
        paste0("chr", chrom)
    }
    if (altChrom %in% chromInfo$name){
        return(altChrom)
    }

    return(NULL)
}


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
    
    if (is.null(.coolMatchChrom(chrom, chromInfo))){
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
        fileChrom <- .coolMatchChrom(chrom, chromInfo)
        chromLength <- chromInfo[chromInfo$name == fileChrom, "length"]
        if (chromstart < 0 | chromstart > chromLength){
            abort(c(glue("chromstart must be between 0 and chromosome end."),
                    "i" = "{chrom} length is {chromLength} bp."))
        }
    }
    
    if (!is.null(altchrom)){
        if (is.null(.coolMatchChrom(altchrom, chromInfo))){
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
        
        fileAltchrom <- .coolMatchChrom(altchrom, chromInfo)
        altchromLength <- chromInfo[chromInfo$name == fileAltchrom, "length"]
        if (!is.null(altchromstart) &&
            (altchromstart < 0 | altchromstart > altchromLength)){
            abort(c(glue("altchromstart must be between 0 and alt chromosome end."),
                    "i" = "{altchrom} length is {altchromLength} bp."))
        }
        
        ## If giving same chrom and altchrom, need to specify
        ## chromstart/chromend and altchromstart/altchromend
        
        if (.coolMatchChrom(chrom, chromInfo) ==
            .coolMatchChrom(altchrom, chromInfo)) {
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
#' @param norm Character value specifying the exact name of a column in the
#' `bins` group of the .(m)cool file (for example `"weight"`, `"sweight"`,
#' `"KR"`, or `"NONE"`). Use \code{readCoolNorms()} to see the actual column
#' names in a file. \code{"NONE"} returns raw counts; \code{"auto"} chooses
#' \code{"KR"} if available, then \code{"weight"}, then \code{"sweight"}, then
#' the first available weight, then \code{"NONE"}.
#' Default value is \code{norm = "NONE"}.
#' @details
#' Cooler normalization weights are multiplicative:
#' \code{balanced[i, j] = observed[i, j] * weight[i] * weight[j]}.
#' Contacts involving bins with non-finite weights are dropped; the bins
#' themselves remain in the coordinate space.
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
#' ## .mcool file (large download, skipped in automated checks)
#' \dontrun{
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
#' }
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

    if (identical(rcool$norm, "auto")) {
        rcool$norm <- .coolAutoNorm(
            file = rcool$file,
            resolution = rcool$resolution
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
    # FILE TYPE AND CHROMOSOME NAMES
    # =========================================================================

    fileType <- .checkCool(rcool$file)
    chrInfo <- readCoolChroms(rcool$file, rcool$resolution)

    ## Match user chromosome names to the names used by the file. `.mcool`
    ## files from some pipelines use `1` while plotgardener users often use
    ## `chr1` (and vice versa).
    rcool$chrom <- .coolMatchChrom(rcool$chrom, chrInfo)
    if (!is.null(rcool$altchrom)){
        rcool$altchrom <- .coolMatchChrom(rcool$altchrom, chrInfo)
    }

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
            ## For same-chromosome off-diagonal regions, read the full
            ## symmetric span on both coordinates. Otherwise bin1 spans the
            ## full region while bin2 is restricted to the smaller region.
            rcool$chromstart <- min(rcool$chromstart, rcool$altchromstart)
            rcool$chromend <- max(rcool$chromend, rcool$altchromend)
            rcool$altchromstart <- rcool$chromstart
            rcool$altchromend <- rcool$chromend
        }
    } else {
        rcool$altchrom <- rcool$chrom
        rcool$altchromstart <- rcool$chromstart
        rcool$altchromend <- rcool$chromend
        
    }
    # =========================================================================
    # EXTRACT SPARSE UPPER TRIANGULAR COUNT MATRIX
    # =========================================================================

    datasetPath <- ""
    if (fileType == ".mcool"){
        datasetPath <- paste0("/resolutions/", as.integer(rcool$resolution))
    }

    ## Find bin IDs for both chrom locations -----------------------------------
    ## adjust starts and ends to start of bins

    outputChrom <- rcool$chrom
    outputAltchrom <- rcool$altchrom

    ## Get vector of chromosome offsets
    chrom_offsets <- h5read(rcool$file,
                            name = paste0(datasetPath, "/indexes/chrom_offset"))

    ## Cooler stores only the upper triangle, ordered by chromosome index.
    ## When the requested `chrom` comes after `altchrom`, read the rectangle in
    ## storage order and swap x/y back afterward. Without this, interchromosomal
    ## queries in reverse chromosome order silently return zero pixels.
    swappedChroms <- FALSE
    if (rcool$chrom != rcool$altchrom) {
        chr1Order <- chrom_offsets[chrInfo[chrInfo$name == rcool$chrom,
                                           "index"]]
        chr2Order <- chrom_offsets[chrInfo[chrInfo$name == rcool$altchrom,
                                           "index"]]
        if (chr1Order > chr2Order) {
            swappedChroms <- TRUE
            oldChrom <- rcool$chrom
            oldStart <- rcool$chromstart
            oldEnd <- rcool$chromend
            rcool$chrom <- rcool$altchrom
            rcool$chromstart <- rcool$altchromstart
            rcool$chromend <- rcool$altchromend
            rcool$altchrom <- oldChrom
            rcool$altchromstart <- oldStart
            rcool$altchromend <- oldEnd
        }
    }

    start1 <- (rcool$chromstart %/% rcool$resolution) * rcool$resolution
    start2 <- (rcool$altchromstart %/% rcool$resolution) * rcool$resolution

    end1 <- (rcool$chromend %/% rcool$resolution) * rcool$resolution
    end2 <- (rcool$altchromend %/% rcool$resolution) * rcool$resolution

    ## Find the bin containing a coordinate. `findInterval()` is robust to
    ## regions that do not start exactly on a bin boundary.
    findBin <- function(coord, starts, chromOffset){
        bin <- findInterval(coord, starts)
        bin <- max(1L, min(bin, length(starts)))
        return(as.integer(chromOffset + bin - 1))
    }
    
    ## find bin id for chr1:start1:end1
    chr1indx <- chrInfo[chrInfo$name == rcool$chrom, "index"]
    
    ## get start1 and end1 bins
    chr1_starts <- h5read(rcool$file, name = paste0(datasetPath, "/bins/start"),
                          index = list((chrom_offsets[chr1indx]+1):
                                           chrom_offsets[chr1indx+1]))
    start1bin <- findBin(start1, chr1_starts, chrom_offsets[chr1indx])
    end1bin <- findBin(end1, chr1_starts, chrom_offsets[chr1indx])
    if (end1bin < start1bin){
        end1bin <- start1bin
    }
    
    ## find bin id for chr2:start2:end2
    chr2indx <- chrInfo[chrInfo$name == rcool$altchrom, "index"]
    
    ## get start2 and end2 bins
    chr2_starts <- h5read(rcool$file, name = paste0(datasetPath,"/bins/start"),
                          index = list((chrom_offsets[chr2indx]+1):
                                           chrom_offsets[chr2indx+1]))
    start2bin <- findBin(start2, chr2_starts, chrom_offsets[chr2indx])
    end2bin <- findBin(end2, chr2_starts, chrom_offsets[chr2indx])
    if (end2bin < start2bin){
        end2bin <- start2bin
    }
    
    ## Extract counts and match locations --------------------------------------
    
    ## Get bin offsets for counts. `bin_offsets` has length nbins + 1; the
    ## pixels for bin `b` (0-based) are `bin_offsets[b + 1] + 1` through
    ## `bin_offsets[b + 2]`. The old code stopped one bin early and dropped
    ## most of the last bin's pixels.
    bin_offsets <- h5read(rcool$file, 
                          name = paste0(datasetPath,"/indexes/bin1_offset"))
                          
    firstPixel <- bin_offsets[start1bin + 1] + 1
    lastPixel <- bin_offsets[end1bin + 2]
    
    ## Pull all bin2 ids for interactions with all bin1s between start1 and end1
    ## Read `binChunkSize` indices at a time for more efficient calling
    if (lastPixel >= firstPixel) {
        if ((lastPixel - firstPixel + 1) > rcool$binChunkSize){
            binChunks <- seq(firstPixel, lastPixel, by = rcool$binChunkSize)
            
            count_idx <- na.omit(unlist(lapply(binChunks, .pullBinChunks, 
                                               file = rcool$file,
                   bin_offsets = bin_offsets, binChunkSize = rcool$binChunkSize,
                   datasetPath = datasetPath, endPixel = lastPixel,
                   start2bin = start2bin, end2bin = end2bin)))
        } else {
            bin2s <- h5read(rcool$file, 
                            name = paste0(datasetPath,"/pixels/bin2_id"),
                            index = list(firstPixel:lastPixel)) 
            
            ## Find indexes for interactions with all bin2s between start2 and end2
            count_idx <- which(bin2s %in% start2bin:end2bin) + 
                firstPixel - 1
        }
    } else {
        count_idx <- integer(0)
    }
    
    ## Pull out bin1 ids, bin2 ids, and counts for all interactions in slice
    if (length(count_idx) > 0){
        bin1ids <- h5read(rcool$file, 
                          name = paste0(datasetPath,"/pixels/bin1_id"), 
                          index = list(count_idx))
        
        bin2ids <- h5read(rcool$file, 
                          name = paste0(datasetPath,"/pixels/bin2_id"), 
                          index = list(count_idx))
        
        counts <- h5read(rcool$file, 
                         name = paste0(datasetPath,"/pixels/count"), 
                         index = list(count_idx)) |>
            #defaults to integer, needs to be numeric to convert to NA correctly
            as.numeric()
    } else {
        bin1ids <- integer(0)
        bin2ids <- integer(0)
        counts <- numeric(0)
    }
    
    ## Multiply by normalization factors, if applicable. Cooler weights are
    ## multiplicative: balanced = observed * weight_i * weight_j.
    if (rcool$norm != "NONE" && length(counts) > 0){
        weights <- h5read(rcool$file, 
                          name = paste0(datasetPath,"/bins/", rcool$norm))
        bin1norm <- weights[as.numeric(bin1ids) + 1]
        bin2norm <- weights[as.numeric(bin2ids) + 1]
        counts <- counts * bin1norm * bin2norm
        counts[!is.finite(counts)] <- NA_real_
    }
    
    ## Get genomic locations for each bin id 
    if (length(count_idx) > 0){
        xs <- h5read(rcool$file, 
                     name = paste0(datasetPath,"/bins/start"),
                     index = list(as.numeric(bin1ids+1)))
        ys <- h5read(rcool$file, 
                     name = paste0(datasetPath,"/bins/start"),
                     index = list(as.numeric(bin2ids+1)))
    } else {
        xs <- numeric(0)
        ys <- numeric(0)
    }
    
    if (swappedChroms) {
        tmp <- xs
        xs <- ys
        ys <- tmp
    }

    ## create sparse upper triangular matrix
    upper <- data.frame(x = xs, y = ys, counts = counts)
    
    # =========================================================================
    # SCALE DATA WITH ZRANGE
    # =========================================================================
    upper <- scale_data(upper = upper, zrange = rcool$zrange)
    
    # =========================================================================
    # FORMAT DATA IN PROPER ORDER AND WITH LABELS
    # =========================================================================
    
    ## Rename columns based on original chrom/altchrom
    if (outputChrom == outputAltchrom){
        colnames(upper) <- c(paste0(outputChrom, "_A"),
                             paste0(outputAltchrom, "_B"),
                             "counts")
    } else {
        colnames(upper) <- c(outputChrom, outputAltchrom, "counts")
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