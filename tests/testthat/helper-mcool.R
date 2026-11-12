## Tests for the .cool/.mcool reader. These use a small synthetic HDF5 file so
## they do not depend on network access or an external cooler executable.

make_test_mcool <- function(path, chromNames = c("chr1", "chr2")) {
    skip_if_not_installed("rhdf5")

    ## Two chromosomes: chr1 has 10 bins, chr2 has 8 bins (100 bp bins)
    nBins1 <- 10L
    nBins2 <- 8L
    starts1 <- seq(0L, by = 100L, length.out = nBins1)
    starts2 <- seq(0L, by = 100L, length.out = nBins2)
    starts <- c(starts1, starts2)
    ends <- starts + 100L
    chromIds <- c(rep(0L, nBins1), rep(1L, nBins2))
    nbins <- length(starts)

    ## Sparse upper-triangular pixels
    pixels <- expand.grid(
        bin1_id = seq_len(nbins) - 1L,
        bin2_id = seq_len(nbins) - 1L
    )
    pixels <- pixels[pixels$bin1_id <= pixels$bin2_id, ]
    pixels <- pixels[order(pixels$bin1_id, pixels$bin2_id), ]
    rownames(pixels) <- NULL
    pixels$count <- as.integer(seq_len(nrow(pixels)) %% 7L + 1L)

    bin1Tab <- tabulate(pixels$bin1_id + 1L, nbins = nbins)
    bin1Offset <- as.integer(c(0L, cumsum(bin1Tab)))
    chromOffset <- as.integer(c(0L, nBins1, nBins1 + nBins2))

    ## Deliberately include a masked bin (NaN weight) and a second weight type
    weight <- c(1, 0.5, 2, 4, 0.25, 1, 1, 1, 1, NaN,
                1, 1, 0.5, 0.5, 2, 2, 1, 1)
    sweight <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                 1, 1, 1, 1, 1, 1, 1, 1)

    h5 <- rhdf5::h5createFile(path)
    rhdf5::h5createGroup(path, "resolutions")
    rhdf5::h5createGroup(path, "resolutions/100")
    rhdf5::h5createGroup(path, "resolutions/100/bins")
    rhdf5::h5createGroup(path, "resolutions/100/chroms")
    rhdf5::h5createGroup(path, "resolutions/100/indexes")
    rhdf5::h5createGroup(path, "resolutions/100/pixels")
    ## Extra root groups must not confuse .checkCool()
    rhdf5::h5createGroup(path, "extra_root_group")

    rhdf5::h5write(chromIds, path, "resolutions/100/bins/chrom")
    rhdf5::h5write(starts, path, "resolutions/100/bins/start")
    rhdf5::h5write(ends, path, "resolutions/100/bins/end")
    rhdf5::h5write(weight, path, "resolutions/100/bins/weight")
    rhdf5::h5write(sweight, path, "resolutions/100/bins/sweight")

    rhdf5::h5write(chromNames, path, "resolutions/100/chroms/name")
    rhdf5::h5write(c(max(ends[chromIds == 0]), max(ends[chromIds == 1])),
                   path, "resolutions/100/chroms/length")

    rhdf5::h5write(bin1Offset, path, "resolutions/100/indexes/bin1_offset")
    rhdf5::h5write(chromOffset, path, "resolutions/100/indexes/chrom_offset")

    rhdf5::h5write(pixels$bin1_id, path, "resolutions/100/pixels/bin1_id")
    rhdf5::h5write(pixels$bin2_id, path, "resolutions/100/pixels/bin2_id")
    rhdf5::h5write(pixels$count, path, "resolutions/100/pixels/count")
    rhdf5::h5closeAll()

    list(
        path = path,
        pixels = pixels,
        starts = starts,
        weight = weight,
        sweight = sweight,
        nBins1 = nBins1
    )
}
