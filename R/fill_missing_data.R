fill_missing_data <- function(dataframe, chrom, chromstart, chromend, altchrom = NULL, altchromstart = NULL, altchromend = NULL, resolution){

  if (!is.null(altchrom)){

    if (chrom == altchrom){

      ## SAME CHROM, OFF DIAGONAL ##

      if (!is.null(chromstart) & !is.null(chromend) & is.null(altchromstart) & is.null(altchromend)){

        ## dataframe[,1] is whole chrom/altchrom
        ## dataframe[,2] is chromstart to chromend

        minchromstart <- min(dataframe[,1])
        minchromend <- max(dataframe[,1])
        maxchromstart <- chromstart
        maxchromend <- chromend

      } else if (is.null(chromstart) & is.null(chromend) & !is.null(altchromstart) & !is.null(altchromend)){

        ## dataframe[,1] is whole chrom/altchrom
        ## dataframe[,2] is altchromstart to altchromend

        minchromstart <- min(dataframe[,1])
        minchromend <- max(dataframe[,1])
        maxchromstart <- altchromstart
        maxchromend <- altchromend

      } else if (is.null(chromstart) & is.null(chromend) & is.null(altchromstart) & is.null(altchromend)){

        ## dataframe[,1] is whole chrom/altchrom
        ## dataframe[,2] is whole chrom/altchrom

        minchromstart <- min(dataframe[,1])
        minchromend <- max(dataframe[,1])
        maxchromstart <- min(dataframe[,2])
        maxchromend <- max(dataframe[,2])

      } else {

        ## dataframe[,1] is min(chromstart, altchromstart) to min(chromend, altchromend)
        ## dataframe[,2] is max(chromstart, altchromstart) to max(chromend, altchromend)

        minchromstart <- min(chromstart, altchromstart)
        minchromend <- min(chromend, altchromend)
        maxchromstart <- max(chromstart, altchromstart)
        maxchromend <- max(chromend, altchromend)

      }

      ## Get necessary sequence of coordinates for col1 and col2 based on resolution spacing
      col1_coords <- seq(minchromstart, minchromend, resolution)
      col2_coords <- seq(maxchromstart, maxchromend, resolution)

      ## Get what is missing
      col1_missing <- setdiff(col1_coords, unique(dataframe[,1]))
      col2_missing <- setdiff(col2_coords, unique(dataframe[,2]))

      ## Get a complement coordinate for missing coordinate
      col1_missing_comp <- rep(dataframe[,2][1], length(col1_missing))
      col2_missing_comp <- rep(dataframe[,1][1], length(col2_missing))

      ## Add 0's for these new complementary pairs
      add_col1 <- data.frame(chrom = col1_missing, altchrom = col1_missing_comp, "counts" = rep(0, length(col1_missing)))
      add_col2 <- data.frame(chrom = col2_missing_comp, altchrom = col2_missing, "counts" = rep(0, length(col2_missing)))

      ## Rename columns to combine dataframes
      colnames(add_col1) <- c("x", "y", "counts")
      colnames(add_col2) <- c("x", "y", "counts")

      ## Combine dataframes
      added_dataframe <- rbind(dataframe, add_col1, add_col2)

      ## Make data symmetric just in case need to include diagonal
      symmetric_dataframe <- added_dataframe[,c(2, 1, 3)]
      colnames(symmetric_dataframe) <- c("x", "y", "counts")

      complete_dataframe <- unique(rbind(added_dataframe, symmetric_dataframe))

      ## Subset again
      complete_dataframe <- complete_dataframe[which(complete_dataframe[,1] >= minchromstart & complete_dataframe[,1] <= minchromend
                                   & complete_dataframe[,2] >= maxchromstart & complete_dataframe[,2] <= maxchromend),]

      ## Cast into matrix and complete upper values
      matrix <- as.matrix(reshape::cast(complete_dataframe, formula = x ~ y, value = "counts"))
      matrix[is.na(matrix)] <- 0

      ## Get back into dataframe format
      added_dataframe_complete <- data.frame(x = as.numeric(rownames(matrix)[row(matrix)]), y = as.numeric(colnames(matrix)[col(matrix)]), counts = as.numeric(c(matrix)))


    } else {

      ## Separate altchrom/chrom to just get number and determine which is larger
      chrom <- as.numeric(gsub(pattern = "chr", replacement = "", x = chrom))
      altchrom <- as.numeric(gsub(pattern = "chr", replacement = "", x = altchrom))

      if (chrom > altchrom){

        ## ALTCHROM IS DATAFRAME[,1] AND CHROM IS DATAFRAME[,2]

        if (is.null(chromstart) & is.null(chromend) & !is.null(altchromstart) & !is.null(altchromend)){

          chromstart <- min(dataframe[,2])
          chromend <- max(dataframe[,2])

        } else if (!is.null(chromstart) & !is.null(chromend) & is.null(altchromstart) & is.null(altchromstart)){

          altchromstart <- min(dataframe[,1])
          altchromend <- max(dataframe[,1])

        } else if (is.null(chromstart) & is.null(chromend) & is.null(altchromstart) & is.null(altchromend)){

          chromstart <- min(dataframe[,2])
          chromend <- max(dataframe[,2])
          altchromstart <- min(dataframe[,1])
          altchromend <- max(dataframe[,1])

        }

        ## Get necessary sequence of coordinates for chrom and altchrom based on resolution spacing
        chrom_coords <- seq(chromstart, chromend, resolution)
        altchrom_coords <- seq(altchromstart, altchromend, resolution)

        ## Get what is missing
        chrom_missing <- setdiff(chrom_coords, unique(dataframe[,2]))
        altchrom_missing <- setdiff(altchrom_coords, unique(dataframe[,1]))

        ## Get a complement coordinate for missing coordinate
        chrom_missing_comp <- rep(dataframe[,1][1], length(chrom_missing))
        altchrom_missing_comp <- rep(dataframe[,2][1], length(altchrom_missing))

        ## Add 0's for these new complementary pairs
        add_chrom <- data.frame(altchrom = chrom_missing_comp, chrom = chrom_missing, "counts" = rep(0, length(chrom_missing)))
        add_altchrom <- data.frame(altchrom = altchrom_missing, chrom = altchrom_missing_comp, "counts" = rep(0, length(altchrom_missing)))


      } else if (altchrom > chrom){

        ## CHROM IS DATAFRAME[,1] AND ALTCHROM IS DATAFRAME[,2]

        if (is.null(chromstart) & is.null(chromend) & !is.null(altchromstart) & !is.null(altchromend)){

          chromstart <- min(dataframe[,1])
          chromend <- max(dataframe[,1])

        } else if (!is.null(chromstart) & !is.null(chromend) & is.null(altchromstart) & is.null(altchromstart)){

          altchromstart <- min(dataframe[,2])
          altchromend <- max(dataframe[,2])

        } else if (is.null(chromstart) & is.null(chromend) & is.null(altchromstart) & is.null(altchromend)){

          chromstart <- min(dataframe[,1])
          chromend <- max(dataframe[,1])
          altchromstart <- min(dataframe[,2])
          altchromend <- max(dataframe[,2])

        }

        ## Get necessary sequence of coordinates for chrom and altchrom based on resolution spacing
        chrom_coords <- seq(chromstart, chromend, resolution)
        altchrom_coords <- seq(altchromstart, altchromend, resolution)

        ## Get what is missing
        chrom_missing <- setdiff(chrom_coords, unique(dataframe[,1]))
        altchrom_missing <- setdiff(altchrom_coords, unique(dataframe[,2]))

        ## Get a complement coordinate for missing coordinate
        chrom_missing_comp <- rep(dataframe[,2][1], length(chrom_missing))
        altchrom_missing_comp <- rep(dataframe[,1][1], length(altchrom_missing))

        ## Add 0's for these new complementary pairs
        add_chrom <- data.frame(chrom = chrom_missing, altchrom = chrom_missing_comp, "counts" = rep(0, length(chrom_missing)))
        add_altchrom <- data.frame(chrom = altchrom_missing_comp, altchrom = altchrom_missing, "counts" = rep(0, length(altchrom_missing)))

      }

      ## Rename columns to combine dataframes
      colnames(add_chrom) <- c("x", "y", "counts")
      colnames(add_altchrom) <- c("x", "y", "counts")

      ## Combine and complete dataframes
      added_dataframe <- rbind(dataframe, add_chrom, add_altchrom)
      added_dataframe_complete <- tidyr::complete(added_dataframe, x, y)
      added_dataframe_complete$counts[is.na(added_dataframe_complete$counts)] <- 0

    }

  } else {

    if (is.null(chromstart) & is.null(chromend)){

      chromstart <- min(dataframe[,1])
      chromend <- max(dataframe[,1])


    }

    ## Get necessary sequence of coordinates based on resolution spacing
    chrom_coords <- seq(chromstart, chromend, resolution)

    ## Get what is missing
    chrom_missing <- setdiff(chrom_coords, unique(dataframe[,1]))

    ## Get a complement coordinate for missing coordinate
    chrom_missing_comp <- rep(dataframe[,2][1], length(chrom_missing))

    ## Add 0's for these new complementary pairs
    add_chrom <- data.frame(chrom = chrom_missing, chrom = chrom_missing_comp, "counts" = rep(0, length(chrom_missing)))

    ## Rename columns to combine dataframes
    colnames(add_chrom) <- c("x", "y", "counts")

    ## Combine dataframes
    added_dataframe <- rbind(dataframe, add_chrom)

    ## Cast into matrix and complete upper values
    matrix <- as.matrix(reshape::cast(added_dataframe, formula = x ~ y, value = "counts"))
    matrix[is.na(matrix)] <- 0
    matrix[lower.tri(matrix)] <- NA

    ## Get back into dataframe format
    added_dataframe_complete <- data.frame(x = rownames(matrix)[row(matrix)], y = colnames(matrix)[col(matrix)], counts = c(matrix))

    ## Remove NA values
    added_dataframe_complete <- added_dataframe_complete[complete.cases(added_dataframe_complete),]

  }

  return(as.data.frame(added_dataframe_complete))
}
