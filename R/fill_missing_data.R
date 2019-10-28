## Function to fill the missing gaps of data
fill_missing_data <- function(dataframe, chrom, chromstart, chromend, altchrom, altchromstart, altchromend, resolution, fill_missing){

  # ======================================================================================================================================================================================
  # FUNCTIONS
  # ======================================================================================================================================================================================

  check_coords_chrom <- function(chromstart, chromend){

    if (is.null(chromstart) | is.null(chromend)){

      return(TRUE)

    }

    return(FALSE)
  }

  check_coords_altchrom <- function(chromstart, chromend, altchromstart, altchromend){

    if (is.null(chromstart) | is.null(chromend) | is.null(altchromstart) | is.null(altchromend)){

      return(TRUE)

    }

    return(FALSE)
  }

  get_chromstart_chromend <- function(dataframe, chrom, altchrom){

    if (is.null(altchrom)){

      new_chromstart <- min(dataframe[,1], dataframe[,2])
      new_chromend <- max(dataframe[,1], dataframe[,2])

      return(list("new_chromstart" = new_chromstart, "new_chromend" = new_chromend))

    } else {

        if (is.null(chromstart) & is.null(chromend) & !is.null(altchromstart) & !is.null(altchromend)){

          new_chromstart <- min(dataframe[,1])
          new_chromend <- max(dataframe[,1])
          new_altchromstart <- altchromstart
          new_altchromend <- altchromend

        } else if (!is.null(chromstart) & !is.null(chromend) & is.null(altchromstart) & is.null(altchromstart)){

          new_chromstart <- chromstart
          new_chromend <- chromend
          new_altchromstart <- min(dataframe[,2])
          new_altchromend <- max(dataframe[,2])

        } else if (is.null(chromstart) & is.null(chromend) & is.null(altchromstart) & is.null(altchromend)){

          new_chromstart <- min(dataframe[,1])
          new_chromend <- max(dataframe[,1])
          new_altchromstart <- min(dataframe[,2])
          new_altchromend <- max(dataframe[,2])

        }

      return(list("new_chromstart" = new_chromstart, "new_chromend" = new_chromend, "new_altchromstart" = new_altchromstart,
                  "new_altchromend" = new_altchromend))
    }

  }

  fill_gaps <- function(dataframe, chrom, chromstart, chromend, altchrom, altchromstart, altchromend, resolution, fill_missing){

    ## Get necessary sequence of coordinates for chrom and altchrom based on resolution spacing
    if (is.null(altchrom)){

      chrom_coords <- seq(chromstart, chromend, resolution)
      altchrom_coords <- chrom_coords

    } else {

      if (chrom == altchrom){

        min_chromstart <- min(chromstart, altchromstart)
        max_chromend <- max(chromend, altchromend)

        chrom_coords <- seq(min_chromstart, max_chromend, resolution)
        altchrom_coords <- chrom_coords

      } else {

        chrom_coords <- seq(chromstart, chromend, resolution)
        altchrom_coords <- seq(altchromstart, altchromend, resolution)

      }

    }

    ## Get what is missing
    chrom_missing <- setdiff(chrom_coords, unique(dataframe[,1]))
    altchrom_missing <- setdiff(altchrom_coords, unique(dataframe[,2]))

    ## Get a complement coordinate for missing coordinates
    chrom_missing_comp <- rep(dataframe[,2][1], length(chrom_missing))
    altchrom_missing_comp <- rep(dataframe[,1][1], length(altchrom_missing))

    ## Add 0's for these new complementary pairs
    if (fill_missing == T){

      add_chrom <- data.frame(chrom = chrom_missing, altchrom = chrom_missing_comp, "counts" = rep(0, length(chrom_missing)))
      add_altchrom <- data.frame(chrom = altchrom_missing_comp, altchrom = altchrom_missing, "counts" = rep(0, length(altchrom_missing)))

    } else if (fill_missing == F){

      add_chrom <- data.frame(chrom = chrom_missing, altchrom = chrom_missing_comp, "counts" = rep(NA, length(chrom_missing)))
      add_altchrom <- data.frame(chrom = altchrom_missing_comp, altchrom = altchrom_missing, "counts" = rep(NA, length(altchrom_missing)))

    }


    ## Rename columns to combine dataframes
    colnames(add_chrom) <- c("x", "y", "counts")
    colnames(add_altchrom) <- c("x", "y", "counts")

    ## Combine dataframes
    added_dataframe <- rbind(dataframe, add_chrom, add_altchrom)

    return(added_dataframe)

  }

  upper_matrix_complete <- function(dataframe){

    ## Cast into matrix, filling NA values with 0's
    data.table::setDT(dataframe)
    matrix <- as.data.frame(data.table::dcast(dataframe, x ~ y, value.var = "counts", fill = 0))

    ## First column is rownames
    rownames(matrix) <- matrix[,1]
    matrix <- matrix[,-1]

    ## Fill lower triangulars with NA's again
    matrix <- data.matrix(matrix)
    matrix[lower.tri(matrix)] <- NA

    ## Put matrix back into dataframe format
    complete_dataframe <- setNames(melt(matrix), c("x", "y", "counts"))

    ## Remove NA values
    data.table::setDT(complete_dataframe)
    subset_dataframe <- subset(complete_dataframe, !is.na(counts))

    return(as.data.frame(subset_dataframe))

  }

  upper_matrix_incomplete <- function(dataframe){

    ## Cast into matrix, filling NA values
    data.table::setDT(dataframe)
    matrix <- as.data.frame(data.table::dcast(dataframe, x ~ y, value.var = "counts"))

    ## First column is rownames
    rownames(matrix) <- matrix[,1]
    matrix <- matrix[,-1]

    ## Fill lower triangulars with 'remove' tag
    matrix <- data.matrix(matrix)
    matrix[lower.tri(matrix)] <- "remove"

    ## Put matrix back into dataframe format
    complete_dataframe <- setNames(melt(matrix), c("x", "y", "counts"))

    ## Remove 'remove' values
    data.table::setDT(complete_dataframe)
    subset_dataframe <- subset(complete_dataframe, counts != "remove" | is.na(counts))

    return(as.data.frame(subset_dataframe))

  }

  # ======================================================================================================================================================================================
  # SET MISSING CHROMSTART/CHROMEND AND/OR ALTCHROMSTART/ALTCHROMEND
  # ======================================================================================================================================================================================

  if (is.null(altchrom)){

    missing <- check_coords_chrom(chromstart = chromstart, chromend = chromend)

    if (missing == T){

      new_coords <- get_chromstart_chromend(dataframe = dataframe, chrom = chrom, altchrom = altchrom)
      chromstart <- new_coords$new_chromstart
      chromend <- new_coords$new_chromend

    }

  } else {

    missing <- check_coords_altchrom(chromstart = chromstart, chromend = chromend, altchromstart = altchromstart, altchromend = altchromend)

    if (missing == T){

      new_coords <- get_chromstart_chromend(dataframe = dataframe, chrom = chrom, altchrom = altchrom)
      chromstart <- new_coords$new_chromstart
      chromend <- new_coords$new_chromend
      altchromstart <- new_coords$new_altchromstart
      altchromend <- new_coords$new_altchromend

    }

  }

  # ======================================================================================================================================================================================
  # GET A DATAFRAME WITH LARGE GAPS FILLED
  # ======================================================================================================================================================================================

    filled_gaps <- fill_gaps(dataframe = dataframe, chrom = chrom, chromstart = chromstart, chromend = chromend, altchrom = altchrom,
                             altchromstart = altchromstart, altchromend = altchromend, resolution = resolution, fill_missing = fill_missing)


  # ======================================================================================================================================================================================
  # MAKE FURTHER COMPLETIONS TO DATA BASED ON CHROM/ALTCHROM INTERACTION
  # ======================================================================================================================================================================================

    if (is.null(altchrom)){
    ## No altchrom, just need to cast into matrix and complete upper triangular values

      if (fill_missing == T){

        filled_data <- upper_matrix_complete(filled_gaps)

      } else if (fill_missing == F){

        filled_data <- upper_matrix_incomplete(filled_gaps)
      }


    } else {

      if (chrom == altchrom){

        ## Cast into matrix and complete upper triangular values
        if (fill_missing == T){

          filled_data <- upper_matrix_complete(filled_gaps)

        } else if (fill_missing == F){

          filled_data <- upper_matrix_incomplete(filled_gaps)

        }

        ## Make symmetric
        lower_dataframe <- filled_data[,c(2, 1, 3)]
        colnames(lower_dataframe) <- c("x", "y", "counts")
        symmetric_dataframe <- unique(rbind(filled_data, lower_dataframe))

        ## Subset data for chromstart/chromend in column 1 and altchromstart/altchromend in column 2
        filled_data <- symmetric_dataframe[which(symmetric_dataframe[,1] >= chromstart & symmetric_dataframe[,1] <= chromend
                                               & symmetric_dataframe[,2] >= altchromstart & symmetric_dataframe[,2] <= altchromend),]

      } else {
        ## different chrom and altchrom, need to complete data
        filled_data <- tidyr::complete(filled_gaps, x, y)

        ## if fill_missing == T, fill NA's with 0's
        if (fill_missing == T){

          filled_data$counts[is.na(filled_data$counts)] <- 0

        }

      }

    }

    return(filled_data)

}
