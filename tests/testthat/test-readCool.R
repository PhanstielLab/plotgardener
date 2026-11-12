test_that(".checkCool recognizes mcool files with extra root groups", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)
    expect_equal(.checkCool(tmp), ".mcool")
})

test_that("readCool returns every pixel in a multi-bin region", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    got <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 300,
        resolution = 100, norm = "NONE"
    )

    expected <- fx$pixels[
        fx$pixels$bin1_id %in% 0:3 & fx$pixels$bin2_id %in% 0:3,
    ]
    expect_equal(nrow(got), nrow(expected))
    expect_setequal(got$counts, expected$count)

    ## The last bin (bin 3) must actually be present
    expect_true(any(got[[1]] == 300 | got[[2]] == 300))
})

test_that("readCool uses multiplicative cooler weights", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    got <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 300,
        resolution = 100, norm = "weight"
    )

    expected <- fx$pixels[
        fx$pixels$bin1_id %in% 0:3 & fx$pixels$bin2_id %in% 0:3,
    ]
    expectedCounts <- expected$count *
        fx$weight[expected$bin1_id + 1] *
        fx$weight[expected$bin2_id + 1]

    expect_equal(got$counts, expectedCounts)
})

test_that("readCool supports sweight and auto normalization", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    expect_equal(readCoolBpResolutions(tmp), 100)
    expect_true("sweight" %in% readCoolNorms(tmp, resolution = 100))
    expect_true("weight" %in% readCoolNorms(tmp, resolution = 100))
    
    expect_equal(.coolAutoNorm(tmp, resolution = 100), "weight")

    got <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 200,
        resolution = 100, norm = "sweight"
    )
    expect_true(nrow(got) > 0)

    gotAuto <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 200,
        resolution = 100, norm = "auto"
    )
    expect_true(nrow(gotAuto) > 0)

    gotWeight <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 200,
        resolution = 100, norm = "weight"
    )
    expect_true(nrow(gotWeight) > 0)
})

test_that("readCool matches chromosome names with and without chr", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    got <- readCool(
        file = tmp, chrom = "1",
        chromstart = 0, chromend = 200,
        resolution = 100, norm = "NONE"
    )
    expect_true(nrow(got) > 0)
    expect_true(all(grepl("chr1", colnames(got)[1:2])))
})

test_that("readCool handles regions that are not bin-aligned", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    got <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 51, chromend = 349,
        resolution = 100, norm = "NONE"
    )
    expected <- fx$pixels[
        fx$pixels$bin1_id %in% 0:3 & fx$pixels$bin2_id %in% 0:3,
    ]
    expect_equal(nrow(got), nrow(expected))
})

test_that("readCool returns interchromosomal contacts", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    got <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 300,
        altchrom = "chr2", altchromstart = 0, altchromend = 300,
        resolution = 100, norm = "NONE"
    )

    expected <- fx$pixels[
        fx$pixels$bin1_id %in% 0:3 & fx$pixels$bin2_id %in% 10:13,
    ]
    expect_equal(nrow(got), nrow(expected))
    expect_setequal(got$counts, expected$count)
})

test_that("readCool chunked reading matches unchunked reading", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    full <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 500,
        resolution = 100, norm = "NONE"
    )
    chunked <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 500,
        resolution = 100, norm = "NONE",
        binChunkSize = 5
    )

    expect_equal(chunked, full)
})

test_that("readCool drops contacts with non-finite weights", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    got <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 900, chromend = 1000,
        resolution = 100, norm = "weight"
    )

    ## bin 9 has a NaN weight, so all of its contacts should be dropped
    expect_equal(nrow(got), 0)
})

test_that("readCool reads interchromosomal pixels in reverse chromosome order", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    ## chr2 comes after chr1 in the file, so this query is in reverse order.
    got <- readCool(
        file = tmp, chrom = "chr2",
        chromstart = 0, chromend = 300,
        altchrom = "chr1", altchromstart = 0, altchromend = 300,
        resolution = 100, norm = "NONE"
    )
    expected <- fx$pixels[
        fx$pixels$bin1_id %in% 0:3 & fx$pixels$bin2_id %in% 10:13,
    ]
    expect_equal(nrow(got), nrow(expected))
    expect_setequal(got$counts, expected$count)
})

test_that("readCool expands same-chromosome off-diagonal regions symmetrically", {
    tmp <- withr::local_tempfile(fileext = ".mcool")
    fx <- make_test_mcool(tmp)

    got <- readCool(
        file = tmp, chrom = "chr1",
        chromstart = 0, chromend = 300,
        altchrom = "chr1", altchromstart = 500, altchromend = 800,
        resolution = 100, norm = "NONE"
    )
    expected <- fx$pixels[
        fx$pixels$bin1_id %in% 0:8 & fx$pixels$bin2_id %in% 0:8,
    ]
    expect_equal(nrow(got), nrow(expected))
})
