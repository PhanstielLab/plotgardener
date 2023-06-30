## Define a function for `pgParams` and parameter parsing logic
# @param params pgParams object to override default arguments of
# parent function
# @param defaultArgs List of defaults for each argument of parent function
# @param declaredArgs List of arguments to override all others
# @param class Name of internal function class
parseParams <- function(params = params,
                        defaultArgs = formals(eval(match.call()[[1]])),
                        declaredArgs = lapply(match.call()[-1],
                            eval.parent,
                            n = 2
                        ),
                        class) {
    ## Remove 'params' and '...' from defaultArgs and declaredArgs
    defaultArgs[["params"]] <- NULL
    declaredArgs[["params"]] <- NULL

    if ("..." %in% names(defaultArgs)) {
        defaultArgs[["..."]] <- NULL
    }
    if ("..." %in% names(declaredArgs)) {
        declaredArgs[["..."]] <- NULL
    }

    ## If params are supplied override matching defaultArguments
    if (!is.null(params)) {
        if (is(params, "pgParams")) {
            ## Replace matching defaultArgs with params
            matchedParams <- params[na.omit(sort(match(
                names(defaultArgs),
                names(params)
            )))]
            defaultArgs[na.omit(match(
                names(params),
                names(defaultArgs)
            ))] <- matchedParams
        } else {
            warning("Input object ignored. Object must be a",
                " \'pgParams\' class object.",
                call. = FALSE
            )
        }
    }

    ## Replace default args with declared args
    if (length(declaredArgs) != 0) {
        suppressWarnings(defaultArgs[names(defaultArgs)
        %in% names(declaredArgs)] <- declaredArgs)
    }
    ## Set arguments without default to NULL
    unset <- unlist(lapply(defaultArgs, is.name))
    defaultArgs[unset] <- lapply(lapply(defaultArgs[unset], deparse), as.null)

    ## Add arguments to object and evaluate
    object <- structure(
        .Data = defaultArgs,
        class = class
    )
    object <- lapply(object, eval, rlang::ns_env("plotgardener"))

    return(object)
}

## Define a function that parses gpar parameters
# @param gpList Input gpar() class list
# @param params Internal function object
# @param ... Additional arguments passed in from parent function call
setGP <- function(gpList, params, ...) {
    availGPs <- c("fill", "col", "lty", "lwd", "cex", "fontsize", 
                  "lineheight", "font", "fontfamily", "alpha", "lineend",
                  "linejoin", "linemitre", "lex")
    gpMatches <- params[which(names(params) %in% availGPs)]
    gpList[names(gpMatches)] <- gpMatches
    gpList[names(list(...))] <- list(...)
    ## Reset with fontface first
    if ("fontface" %in% names(gpList)) {
        otherParams <- gpList
        gpList <- gpar(fontface = gpList$fontface)
        gpList[names(otherParams)] <- otherParams
        gpList$font <- NULL
    }

    return(gpList)
}

## Define a function that checks chromstart/chromend 
## region errors
regionErrors <- function(chromstart, chromend){
    
    ## Can't have only one NULL chromstart or chromend
    if ((is.null(chromstart) &
            !is.null(chromend)) |
            (is.null(chromend) &
            !is.null(chromstart))) {
        stop("Cannot have one \'NULL\' \'chromstart\' or ",
            "\'chromend\'.", call. = FALSE)
    }
    
    if (!is.null(chromstart) & !is.null(chromend)){
        
        ## Can't have a 0 bp wide region
        if (chromstart == chromend) {
            stop("Genomic region is 0 bp wide.", call. = FALSE)
        }
        
        ## start can't be larger than end
        if (chromstart > chromend) {
            stop("\'chromstart\' should not be larger than \'chromend\'.",
                call. = FALSE
            )
        }
        
    }
}

## Define a function that checks just and gets corresponding numeric value
justConversion <- function(just){
    
    ## Convert input just to comma-separated string for finding in dataframe
    just <- paste(just, collapse = ", ")
    
    ## Valid combinations of string and numeric justifications
    validJusts <- c(plotJusts$just, plotJusts$numeric)
    
    ## Throw error if not a valid Just
    if (!just %in% validJusts){
        stop("Invalid `just` option.", call. = FALSE)
    }
    
    ## Convert to robust numeric values
    if (!just %in% plotJusts$numeric){
        just <- plotJusts[plotJusts$just == just,]$numeric
    }
    
    just <- as.numeric(unlist(strsplit(just, split = ", ")), 
                                as.is = TRUE)
    
    return(just)
}

## Define a function that checks the format of range parameters
rangeErrors <- function(range){
    
    if (!is.null(range)){
        
        ## Needs to be a vector
        if (!is.vector(range)) {
            stop("\'(z)range\' must be a vector of length 2.", call. = FALSE)
        }
        
        ## Vector needs to be length 2
        if (length(range) != 2) {
            stop("\'(z)range\' must be a vector of length 2.", call. = FALSE)
        }
        
        ## Vector needs to be numbers
        if (!is.numeric(range)) {
            stop("\'(z)range\' must be a vector of two numbers.",
                call. = FALSE
            )
        }
        
        ## second value should be larger than the first value
        if (range[1] >= range[2]) {
            stop("\'(z)range\' must be a vector of two numbers ",
                "in which the 2nd value is larger than the 1st.",
                call. = FALSE
            )
        }
        
    }
    
}

## Define a function that checks for hic input errors
hicErrors <- function(hic, norm){
    
    ## if it's a dataframe or datatable, it needs to be properly formatted
    if (is(hic, "data.frame") && ncol(hic) != 3) {
        stop("Invalid dataframe format.  Input a dataframe with 3 ",
            "columns: chrA, chrB, counts.", call. = FALSE)
    }
    
    if (!is(hic, "data.frame")) {
        
        ## if it's a file path, it needs to be a .hic file
        if (file_ext(hic) != "hic") {
            stop("Invalid input. File must have a \".hic\" extension",
                call. = FALSE
            )
        }
        
        ## if it's a valid .hic file, it needs to have a valid norm
        if (is.null(norm)) {
            stop("If providing .hic file, please specify \'norm\'.",
                call. = FALSE
            )
        }
    }
}

## Define a function to get the assembly info based on a string (ie default)
## or assembly object
# @param assembly assembly input from a plot object
parseAssembly <- function(assembly) {
    
    availDefaults <- default_genomePackages$Genome

    ## If it's just a string, get the default
    if (is(assembly, "character")) {
        if (!assembly %in% availDefaults) {
            stop("\'assembly\' not available as a default. Please make a ",
                "assembly object with `assembly()` or pick an assembly ",
                "from the defaults listed with `genomes()`.",
                call. = FALSE
            )
        }

        assemblyData <- default_genomePackages[
            default_genomePackages$Genome == assembly,]
        assemblyData <- assembly(
            Genome = assemblyData$Genome,
            TxDb = assemblyData$TxDb,
            OrgDb = assemblyData$OrgDb,
            gene.id.column = assemblyData$gene.id.column,
            display.column = assemblyData$display.column,
            BSgenome = assemblyData$BSgenome
        )

        ## If it's a assembly object, use those
    } else if (is(assembly, "assembly")) {
        assemblyData <- assembly
    } else {
        stop("Invalid \'assembly\' type. Please make a assembly object ",
            "with `assembly()` or input an assembly string from the ",
            "defaults listed with `genomes()`.",
            call. = FALSE
        )
    }


    return(assemblyData)
}

## Define a function that converts coordinates/dimensions
## for standard plot objects into default units
# @param object Function object containing x, y, width, height valus
# @param default.units String value of default.units
defaultUnits <- function(object, default.units) {
    if (!(is.null(object$x)) & !(is.null(object[["y"]]))) {
        if (!"unit" %in% class(object$x)) {
            if (!is.numeric(object$x)) {
                stop("x-coordinate is neither a unit object nor a ",
                    "numeric value. Cannot place object.",
                    call. = FALSE
                )
            }

            if (is.null(default.units)) {
                stop("x-coordinate detected as numeric.\'default.units\' ",
                    "must be specified.",
                    call. = FALSE
                )
            }

            object$x <- unit(object$x, default.units)
        }


        if (!"unit" %in% class(object$y)) {

            ## Check for "below" y-coord
            if (grepl("b", object$y) == TRUE) {
                if (grepl("^[ac-zA-Z]+$", object$y) == TRUE) {
                    stop("\'below\' y-coordinate detected with additional ",
                        "letters. Cannot parse y-coordinate.",
                        call. = FALSE
                    )
                }

                if (is.na(suppressWarnings(as.numeric(gsub("b", "", 
                                                    object$y))))) {
                    stop("\'below\' y-coordinate does not have a numeric ",
                        "associated with it. Cannot parse y-coordinate.",
                        call. = FALSE
                    )
                }

                object$y <- plot_belowY(y_coord = object$y)
            } else {
                if (!is.numeric(object$y)) {
                    stop("y-coordinate is neither a unit object nor a ",
                        "numeric value. Cannot place plot.",
                        call. = FALSE
                    )
                }

                if (is.null(default.units)) {
                    stop("y-coordinate detected as numeric. \'default.units\' ",
                        "must be specified.",
                        call. = FALSE
                    )
                }

                object$y <- unit(object$y, default.units)
            }
        }

        if (!"unit" %in% class(object$width)) {
            if (!is.numeric(object$width)) {
                stop("width is neither a unit object nor a numeric value. ",
                    "Cannot place plot.",
                    call. = FALSE
                )
            }

            if (is.null(default.units)) {
                stop("width detected as numeric. \'default.units\' must ",
                    "be specified.",
                    call. = FALSE
                )
            }

            object$width <- unit(object$width, default.units)
        }

        if (!"unit" %in% class(object$height)) {
            if (!is.numeric(object$height)) {
                stop("height is neither a unit object nor a numeric ",
                    "value. Cannot place plot.",
                    call. = FALSE
                )
            }

            if (is.null(default.units)) {
                stop("height detected as numeric.v\'default.units\' ",
                    "must be specified.",
                    call. = FALSE
                )
            }

            object$height <- unit(object$height, default.units)
        }
    }

    return(object)
}

## Define a function that converts coordinates/dimensions
## for other miscellaneous elements into default units
# @param value value to be checked/converted
# @param name string name identifying kind of value
# @param default.units String value of default.units
# @param funName String value of function name called within
# @param yBelow Logical indicating whether to allow "below" y coords
misc_defaultUnits <- function(value, name, default.units,
                            funName = NULL, yBelow = TRUE) {
    if (!"unit" %in% class(value)) {
        ## Check for y
        if (grepl("y", name)) {
            ## Check for "below" y-coord
            if (all(grepl("b", value))) {
                if (yBelow == FALSE) {
                    stop("\'below\' ", name, "-coordinate detected. Cannot ",
                        "parse \'below\' ", name, "-coordinate for ",
                        funName, ".",
                        call. = FALSE
                    )
                } else {
                    ## "below" y-coord additional letters
                    if (any(grepl("^[ac-zA-Z]+$", value))) {
                        stop("\below\' ", name, "-coordinate detected with ",
                            "additional letters. Cannot parse ",
                            name, "-coordinate.",
                            call. = FALSE
                        )
                    }
                    ## "below" y-coord missing a number
                    if (any(is.na(suppressWarnings(as.numeric(gsub("b", "", 
                                                            value)))))) {
                        stop("\below\' ", name, "-coordinate does not have ",
                            "a numeric associated with it. Cannot parse ",
                            name, "-coordinate.",
                            call. = FALSE
                        )
                    }

                    value <- unit(
                        unlist(lapply(value, plot_belowY)),
                        get("page_units", envir = pgEnv)
                    )
                }
            } else {
                ## y-coord, not detected as "below"
                if (!is.numeric(value)) {
                    stop(name, "-coordinate is neither a unit object nor a ",
                        "numeric value. Cannot place plot.",
                        call. = FALSE
                    )
                }
                if (is.null(default.units)) {
                    stop(name, "-coordinate detected as numeric. ",
                        "\'default.units\' must be specified.",
                        call. = FALSE
                    )
                }
                value <- unit(value, default.units)
            }
        } else {
            if (grepl("x", name)) {
                name <- paste0(name, "-coordinate")
            }

            if (!is.numeric(value)) {
                stop(name, " is neither a unit object nor a ",
                    "numeric value. Cannot place plot.",
                    call. = FALSE
                )
            }
            if (is.null(default.units)) {
                stop(name, " detected as numeric. ",
                    "\'default.units\' must be specified.",
                    call. = FALSE
                )
            }
            value <- unit(value, default.units)
        }
    }

    return(value)
}

## Define a function to catch errors for functions that don't handle
## colorby and for colorby related errors
checkColorby <- function(fill, colorby = TRUE, data = NULL){
    
    ## colorby not allowed for some functions
    if (colorby == FALSE){
        if (is(fill, "colorby")){
            stop("`fill` cannot be a `colorby` object for this function.",
                call. = FALSE)
        }
    } else {
        ## If colorby allowed, check for appropriate colorby class,
        ## column in data, and number of colorby column in data
        if (!is(fill, "character") & !is(fill, "factor")){
            
            if (!is(fill, "colorby")){
                stop("`colorby` not of class \"colorby\". Input colorby ",
                    "information with `colorby()`.", call. = FALSE)
            }
            
            ## Check for `colorby` column
            if (!any(colnames(data) == fill$column)) {
                stop("`colorby` column not found in data. Check ",
                    "`colorby` column name.",
                    call. = FALSE
                )
            }
            
            if (length(which(colnames(data) == fill$column)) > 1) {
                stop("Multiple matching `colorby` columns found in data. ",
                    "Please provide `colorby` column name with only ",
                    "one occurrence.", call. = FALSE)
            }
            
        }
        
    }
    
}
