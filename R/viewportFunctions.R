## Define a function to grab the name of a viewport
viewport_name <- function(viewport) {
    return(viewport$name)
}

## Define a function to add a plot viewport to internal list of
## valid "last plotted" viewports
addViewport <- function(vpName) {
    pg_vpTree <- get("pg_vpTree", envir = pgEnv)
    pg_vpTree[length(pg_vpTree) + 1] <- vpName
    assign("pg_vpTree", pg_vpTree, envir = pgEnv)
}

## Define a function to get a list of current viewports
current_viewports <- function() {
    if (!"page" %in% names(lapply(
        current.vpTree()$children,
        viewport_name
    ))) {
        current <- as.list(names(lapply(
            current.vpTree()$children,
            viewport_name
        )))
    } else {
        page_children <- names(lapply(
            current.vpTree()$children$page$children,
            viewport_name
        ))
        current <- as.list(page_children)
    }

    return(current)
}

## Define a function to convert viewport x and y into center just
adjust_vpCoords <- function(viewport) {
    vp_y <- viewport$y
    ## left
    if (identical(viewport$justification, c(0, 0.5))) {
        vp_x <- viewport$x + (0.5 * viewport$width)
    } else if (identical(viewport$justification, c(1, 0.5))) {
        ## right
        vp_x <- viewport$x - (0.5 * viewport$width)
    } else if (identical(viewport$justification, c(0.5, 0))) {
        ## bottom
        vp_x <- viewport$x
        vp_y <- vp_y + (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(0.5, 1))) {
        ## top
        vp_x <- viewport$x
        vp_y <- vp_y - (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(0, 1))) {
        ## left top
        vp_x <- viewport$x + (0.5 * viewport$width)
        vp_y <- vp_y - (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(1, 1))) {
        ## right top
        vp_x <- viewport$x - (0.5 * viewport$width)
        vp_y <- vp_y - (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(0, 0))) {
        ## left bottom
        vp_x <- viewport$x + (0.5 * viewport$width)
        vp_y <- vp_y + (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(1, 0))) {
        ## right bottom
        vp_x <- viewport$x - (0.5 * viewport$width)
        vp_y <- vp_y + (0.5 * viewport$height)
    } else {
        ## center
        vp_x <- viewport$x
    }

    return(list(vp_x, vp_y))
}

## Define a function to change viewport x and y-coordinates to top left
vp_topLeft <- function(viewport) {
    if (identical(viewport$justification, c(0, 0.5))) {
        ## left
        vp_x <- viewport$x
        vp_y <- viewport$y + (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(1, 0.5))) {
        ## right
        vp_x <- viewport$x - (viewport$width)
        vp_y <- viewport$y + (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(0.5, 0))) {
        ## bottom
        vp_x <- viewport$x - (0.5 * viewport$width)
        vp_y <- viewport$y + (viewport$height)
    } else if (identical(viewport$justification, c(0.5, 1))) {
        ## top
        vp_x <- viewport$x - (0.5 * viewport$width)
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(0, 1))) {
        ## left top
        vp_x <- viewport$x
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(1, 1))) {
        ## right top
        vp_x <- viewport$x - (viewport$width)
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(0, 0))) {
        ## left bottom
        vp_x <- viewport$x
        vp_y <- viewport$y + (viewport$height)
    } else if (identical(viewport$justification, c(1, 0))) {
        ## right bottom
        vp_x <- viewport$x - (viewport$width)
        vp_y <- viewport$y + (viewport$height)
    } else {
        ## center
        vp_x <- viewport$x - (0.5 * viewport$width)
        vp_y <- viewport$y + (0.5 * viewport$height)
    }

    return(list(vp_x, vp_y))
}

## Define a function to change viewport x and y-coordinates to bottom left
vp_bottomLeft <- function(viewport) {
    if (identical(viewport$justification, c(0, 0.5))) {
        ## left
        vp_x <- viewport$x
        vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(1, 0.5))) {
        ## right
        vp_x <- viewport$x - (viewport$width)
        vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(0.5, 0))) {
        ## bottom
        vp_x <- viewport$x - (0.5 * viewport$width)
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(0.5, 1))) {
        ## top
        vp_x <- viewport$x - (0.5 * viewport$width)
        vp_y <- viewport$y - (viewport$height)
    } else if (identical(viewport$justification, c(0, 1))) {
        ## left top
        vp_x <- viewport$x
        vp_y <- viewport$y - (viewport$height)
    } else if (identical(viewport$justification, c(1, 1))) {
        ## right top
        vp_x <- viewport$x - (viewport$width)
        vp_y <- viewport$y - (viewport$height)
    } else if (identical(viewport$justification, c(0, 0))) {
        ## left bottom
        vp_x <- viewport$x
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(1, 0))) {
        ## right bottom
        vp_x <- viewport$x - (viewport$width)
        vp_y <- viewport$y
    } else {
        ## center
        vp_x <- viewport$x - (0.5 * viewport$width)
        vp_y <- viewport$y - (0.5 * viewport$height)
    }

    return(list(vp_x, vp_y))
}

## Define a function to change viewport x and y-coordinates to bottom right
vp_bottomRight <- function(viewport) {
    if (identical(viewport$justification, c(0, 0.5))) {
        ## left
        vp_x <- viewport$x + viewport$width
        vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(1, 0.5))) {
        ## right
        vp_x <- viewport$x
        vp_y <- viewport$y - (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(0.5, 0))) {
        ## bottom
        vp_x <- viewport$x + (0.5 * viewport$width)
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(0.5, 1))) {
        ## top
        vp_x <- viewport$x + (0.5 * viewport$width)
        vp_y <- viewport$y - (viewport$height)
    } else if (identical(viewport$justification, c(0, 1))) {
        ## left top
        vp_x <- viewport$x + viewport$width
        vp_y <- viewport$y - (viewport$height)
    } else if (identical(viewport$justification, c(1, 1))) {
        ## right top
        vp_x <- viewport$x
        vp_y <- viewport$y - (viewport$height)
    } else if (identical(viewport$justification, c(0, 0))) {
        ## left bottom
        vp_x <- viewport$x + viewport$width
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(1, 0))) {
        ## right bottom
        vp_x <- viewport$x
        vp_y <- viewport$y
    } else {
        ## center
        vp_x <- viewport$x + (0.5 * viewport$width)
        vp_y <- viewport$y - (0.5 * viewport$height)
    }

    return(list(vp_x, vp_y))
}

vp_topRight <- function(viewport) {
    if (identical(viewport$justification, c(0, 0.5))) {
        ## left
        vp_x <- viewport$x + viewport$width
        vp_y <- viewport$y + (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(1, 0.5))) {
        ## right
        vp_x <- viewport$x
        vp_y <- viewport$y + (0.5 * viewport$height)
    } else if (identical(viewport$justification, c(0.5, 0))) {
        ## bottom
        vp_x <- viewport$x + (0.5 * viewport$width)
        vp_y <- viewport$y + viewport$height
    } else if (identical(viewport$justification, c(0.5, 1))) {
        ## top
        vp_x <- viewport$x + (0.5 * viewport$width)
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(0, 1))) {
        ## left top
        vp_x <- viewport$x + viewport$width
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(1, 1))) {
        ## right top
        vp_x <- viewport$x
        vp_y <- viewport$y
    } else if (identical(viewport$justification, c(0, 0))) {
        ## left bottom
        vp_x <- viewport$x + viewport$width
        vp_y <- viewport$y + viewport$height
    } else if (identical(viewport$justification, c(1, 0))) {
        ## right bottom
        vp_x <- viewport$x
        vp_y <- viewport$y + viewport$height
    } else {
        ## center
        vp_x <- viewport$x + (0.5 * viewport$width)
        vp_y <- viewport$y + (0.5 * viewport$height)
    }

    return(list(vp_x, vp_y))
}

## Define a function to convert to page units
convert_page <- function(object) {

    ## Get page_height and its units from pgEnv
    page_height <- get("page_height", envir = pgEnv)
    page_units <- get("page_units", envir = pgEnv)

    ## Convert x and y coordinates and height and width to same page_units
    old_x <- object$x
    old_y <- object$y
    old_height <- object$height
    old_width <- object$width
    new_x <- convertX(old_x, unitTo = page_units)
    new_y <- convertY(unit(page_height,
        units = page_units
    ) - convertY(old_y,
        unitTo = page_units
    ),
    unitTo = page_units
    )
    new_height <- convertHeight(old_height, unitTo = page_units)
    new_width <- convertWidth(old_width, unitTo = page_units)

    object$x <- new_x
    object$y <- new_y
    object$height <- new_height
    object$width <- new_width


    return(object)
}

## Parse y-coordinates for plots relative to bottom of last plot
plot_belowY <- function(y_coord) {
    prevPlots <- unlist(get("pg_vpTree", envir = pgEnv))

    if (length(prevPlots) == 0) {
        stop("No previous plot detected. Cannot define ",
            "a \'below\' y-coordinate.",
            call. = FALSE
        )
    }

    lastPlot <- prevPlots[length(prevPlots)]
    ## Go into viewport to get coordinates/dimensions
    seekViewport(name = lastPlot)
    x <- current.viewport()$x
    y <- current.viewport()$y
    width <- current.viewport()$width
    height <- current.viewport()$height
    just <- current.viewport()$just

    ## Exit viewport
    upViewport()

    ## Convert to bottom
    bottomCoords <- vp_bottomLeft(viewport(
        x = x, y = y,
        width = width, height = height,
        just = just
    ))

    ## Parse number and convert to page units
    space <- utils::type.convert(gsub("b", "", y_coord), as.is = TRUE)
    space <- unit(space, get("page_units", envir = pgEnv))

    new_y <- (unit(
        get("page_height", envir = pgEnv),
        get("page_units", envir = pgEnv)
    ) - bottomCoords[[2]]) + space
    new_y <- convertY(new_y, unitTo = get("page_units", envir = pgEnv))

    return(new_y)
}

## Define a function to distinguish between plots with 
## multiple viewports for certain annotation functions
getAnnoViewport <- function(plot){
    
    if (is(plot, "genes")){
        plotVP <- plot$grobs$children$background$vp
    } else if (is(plot, "hicTriangle") |
            is(plot, "hicRectangle")){
        plotVP <- plot$outsideVP
    } else {
        plotVP <- plot$grobs$vp
    }
    
    return(plotVP)
    
}

## Define a function that converts chromstart/chromend to page units for
## multi-chromosome manhattan plot
convertManhattan <- function(object, manhattanPlot) {

    ## Get assembly data
    if (is(object$assembly$TxDb, "TxDb")) {
        txdbChecks <- TRUE
    } else {
        
        if (!requireNamespace(object$assembly$TxDb, quietly = TRUE)){
            txdbChecks <- FALSE
        } else {
            txdbChecks <- TRUE
        }
        
    }

    if (txdbChecks == TRUE) {
        if (is(object$assembly$TxDb, "TxDb")) {
            tx_db <- object$assembly$TxDb
        } else {
            tx_db <- eval(parse(text = paste0(as.name(object$assembly$TxDb), 
                                            "::", 
                                            as.name(object$assembly$TxDb))))
        }

        assembly_data <- as.data.frame(setDT(as.data.frame(
            GenomeInfoDb::seqlengths(tx_db)
        ),
        keep.rownames = TRUE
        ))
        colnames(assembly_data) <- c("chrom", "length")
        assembly_data <- assembly_data[which(assembly_data[, "chrom"] %in%
            manhattanPlot$chrom), ]

        ## Get the offset based on spacer for the assembly
        offsetAssembly <- spaceChroms(
            assemblyData = assembly_data,
            space = manhattanPlot$space
        )
        offsetAssembly <- offsetAssembly[which(
            offsetAssembly$chrom == object$chrom
        ), ]

        ## Convert chromstart and chromend to chrom offsetAssembly range
        oldRange <- offsetAssembly[, "length"] - 1
        newRange <- offsetAssembly[, "end"] - offsetAssembly[, "start"]
        newStart <- (((object$chromstart - 1) * newRange) / oldRange) +
            offsetAssembly[, "start"]
        newEnd <- (((object$chromend - 1) * newRange) / oldRange) +
            offsetAssembly[, "start"]

        ## Convert new chromstart and chromend to page units
        start <- convertX(unit(newStart, "native"),
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )
        end <- convertX(unit(newEnd, "native"),
            unitTo = get("page_units", envir = pgEnv),
            valueOnly = TRUE
        )

        return(list(start, end))
    }
}
