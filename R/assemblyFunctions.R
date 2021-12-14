## Define a function that checks formatting agreement between given chrom
## and the chromosome in provided data 
# @param data Input data
# @param chrom Input chrom
# @param Type of data, either "ranges" or "pairs"
chromDataAgreement <- function(data, chrom, type){
    chrom_chrCheck <- grepl("chr", chrom)
    if (type == "ranges"){
        ## Just check first column 
        data_chrCheck <- all(grepl("chr", data[,"chrom"]))
    } else if (type == "pairs"){
        ## Check chr1 and chr2 columns
        data_chrCheck <- all(grepl("chr", data[,"chrom1"])) &
            all(grepl("chr", data[,"chrom2"]))
    }
    
    
    if (chrom_chrCheck != data_chrCheck){
        warning("Format of chromosome in data does not match ",
                "format of `chrom`.", call. = FALSE)
    }
    
}

## Define a function that checks for whole chromosome data and
## sets the plot xscale accordingly
# @param object plot object
# @param objectInternal internal plot object
# @param plotType string of plot type to show up in error message
genomicScale <- function(object, objectInternal, plotType) {
    if (is.null(object$chromstart) & is.null(object$chromend)) {
        if (is(object$assembly$TxDb, "TxDb")) {
            txdbChecks <- TRUE
        } else {
            
            if (!requireNamespace(object$assembly$TxDb, quietly = TRUE)){
                txdbChecks <- FALSE
                warning("`", object$assembly$TxDb, "` not available. Please",
                        " load to generate ", plotType, ".", call. = FALSE)
            } else {
                txdbChecks <- TRUE
            }
            
        }
        objectInternal$xscale <- c(0, 1)
        if (txdbChecks == TRUE) {
            if (is(object$assembly$TxDb, "TxDb")) {
                tx_db <- object$assembly$TxDb
            } else {
                tx_db <- eval(parse(text = paste0(as.name(object$assembly$TxDb),
                                            "::",
                                            as.name(object$assembly$TxDb))))
            }
            assembly_data <- GenomeInfoDb::seqlengths(tx_db)
            if (!object$chrom %in% names(assembly_data)) {
                warning("Chromosome ",
                        "'", object$chrom, "' ",
                        "not found in ",
                        "`", tx_db$packageName, "`",
                        " and data for entire chromosome cannot be plotted.",
                        call. = FALSE
                )
            } else {
                object$chromstart <- 1
                object$chromend <- assembly_data[[object$chrom]]
                if (class(object) %in% c("hicTriangle", "hicRectangle")) {
                    object$altchromstart <- 1
                    object$altchromend <- assembly_data[[object$chrom]]
                }
                objectInternal$xscale <- c(object$chromstart, object$chromend)
            }    
            
        }

    } else {
        txdbChecks <- TRUE
        objectInternal$xscale <- c(object$chromstart, object$chromend)
    }

    objectInternal$txdbChecks <- txdbChecks

    return(list(object, objectInternal))
}

## Define a function that checks for and gets gene/transcript data
# @param object plot object
# @param objectInternal internal plot object
geneData <- function(object, objectInternal) {

    ## TxDb
    if (is(object$assembly$TxDb, "TxDb")) {
        txdbChecks <- TRUE
    } else {
        
        if (!requireNamespace(object$assembly$TxDb, quietly = TRUE)){
            txdbChecks <- FALSE
            warning("`", object$assembly$TxDb, "` not available. Please",
                    " load to plot genes or transcripts.", call. = FALSE)
        } else {
            txdbChecks <- TRUE
        }
        
    }

    ## orgDb
    if (!requireNamespace(object$assembly$OrgDb, quietly = TRUE)){
        orgdbChecks <- FALSE
        warning("`", object$assembly$OrgDb, "` not available. Please load",
                " to plot genes or transcripts.", call. = FALSE)
    } else {
        orgdbChecks <- TRUE
    }
    

    ## Data
    data <- data.frame(matrix(ncol = 22, nrow = 0))
    xscale <- c(0, 1)

    if (txdbChecks == TRUE & orgdbChecks == TRUE) {

        ## Load txdb
        if (is(object$assembly$TxDb, "TxDb")) {
            tx_db <- object$assembly$TxDb
        } else {
            tx_db <- eval(parse(text = paste0(as.name(object$assembly$TxDb),
                                            "::",
                                            as.name(object$assembly$TxDb))))
        }

        genome <- GenomeInfoDb::seqlengths(tx_db)

        if (object$assembly$gene.id.column ==
            object$assembly$display.column) {
            objectInternal$displayCol <- "GENEID"
        } else {
            objectInternal$displayCol <- object$assembly$display.column
        }

        if (!object$chrom %in% names(genome)) {
            warning("Chromosome ", "'", object$chrom, "'",
                "not found in ", "`", tx_db$packageName, "`",
                " and data for entire chromosome cannot be plotted.",
                call. = FALSE
            )
        } else {
            if (is.null(object$chromstart) & is.null(object$chromend)) {
                object$chromstart <- 1
                object$chromend <- genome[[object$chrom]]
            }

            data <- getExons(
                assembly = object$assembly,
                chromosome = object$chrom,
                start = object$chromstart,
                stop = object$chromend
            )
            xscale <- c(object$chromstart, object$chromend)
        }
    }

    objectInternal$xscale <- xscale
    objectInternal$data <- data
    return(list(object, objectInternal))
}

## Define a function that will by default prioritize genes by citations
# (if available) or gene lengths
# @param data data frame of gene or transcript data
# @param assembly assembly associated with gene data
# @param transcript a logical indicating whether or not we're
# plotting transcripts or not
defaultGenePriorities <- function(data, assembly, transcript = FALSE, 
                                geneHighlights = NULL, displayCol = "SYMBOL") {

    availCitations <- default_genomePackages[
        !is.na(default_genomePackages$Citations),]$Genome

    ## Define assemblies whose TxDb IDs will need to be converted to
    ## ENTREZID from a different ID
    convertIDs <- list(
        dm3 = "ENSEMBL", dm6 = "ENSEMBL",
        rn4 = "ENSEMBL", sacCer2 = "ORF",
        sacCer3 = "ORF"
    )

    assemblyName <- assembly$Genome
    ## If assembly is included in package, access citations
    if (any(availCitations %in% assemblyName)) {
        
        name <- default_genomePackages[
            default_genomePackages$Genome %in% assemblyName,]$Genome
        
        ## Convert necessary builds to ENTREZID
        if (name %in% names(convertIDs)) {

            ## Load associated OrgDb
            org_db <- eval(parse(text = paste0(as.name(assembly$OrgDb),
                                            "::",
                                            as.name(assembly$OrgDb))))

            ## Convert gene ids in data to ENTREZID based on previous keytype
            entrezIDs <- suppressMessages(
                AnnotationDbi::select(org_db,
                    keys = data$GENEID,
                    columns = "ENTREZID",
                    keytype = convertIDs[[name]]
                )
            )

            data$ENTREZID <- entrezIDs
        } else {
            data$ENTREZID <- data$GENEID
        }

        ## Get internal citation data and match based on ENTREZID
        citationData <- default_genomePackages[
            default_genomePackages$Genome == name,]$Citations
        citationData <- eval(parse(text = citationData))
        updatedData <- suppressMessages(dplyr::left_join(
            x = data,
            y = citationData,
            by = "ENTREZID"
        ))

        ## Set any missing citations to 0
        updatedData[is.na(updatedData$Citations), ]$Citations <-
            rep(0, nrow(updatedData[is.na(updatedData$Citations), ]))

        if (transcript == TRUE) {
            updatedData <-
                updatedData[duplicated(updatedData$TXNAME) == FALSE, ]
        }


        ## Order based on citation number
        updatedData <- updatedData[order(updatedData$Citations,
            decreasing = TRUE
        ), ]
    } else {

        ## With no internal citation data, set default priority
        ## based on gene/transcript length
        updatedData <- data[order(data$length, decreasing = TRUE), ]
    }

    ## Put any gene highlights at the top of the priority
    if (transcript == FALSE){
        if (!is.null(geneHighlights)){
            updatedData <- updatedData[c(which(updatedData[[displayCol]] 
                                        %in% geneHighlights), 
                                        which(!updatedData[[displayCol]] 
                                              %in% geneHighlights)),]
        }
    }
    ## Return data with priorities
    return(updatedData)
}

## Define a function that determines the chromosome offsets for a plot
## with multiple chromosomes (manhattan plots)
# @param assemblyData data.frame of an assembly's chromosomes and lengths
# @param space the space between each chrom as a fraction of plot width
spaceChroms <- function(assemblyData, space) {

    ## Determine the offset for each chromsome
    cumsums <- cumsum(as.numeric(assemblyData[, "length"]))
    spacer <- cumsums[length(cumsums)] * space
    additionalSpace <- (seq(1, length(cumsums) - 0)) * spacer

    ## Start position
    startPos <- c(0, cumsums[seq(1, length(cumsums) - 1)])
    startPos <- startPos + additionalSpace
    assemblyData[, "start"] <- startPos

    ## Stop Position
    stopPos <- cumsums + (seq(1, (length(cumsums)))) * spacer
    assemblyData[, "end"] <- stopPos

    return(assemblyData)
}
