test_that("read_rangeData respects named range columns in any order", {
    df <- data.frame(
        start = c(1, 10),
        chrom = c("chr1", "chr2"),
        end = c(5, 20),
        score = c(1, 2)
    )
    out <- plotgardener:::read_rangeData(
        data = df, assembly = "hg38", type = "range"
    )
    expect_equal(colnames(out)[1:3], c("chrom", "start", "end"))
    expect_equal(out$chrom, c("chr1", "chr2"))
    expect_equal(out$start, c(1, 10))
    expect_equal(out$end, c(5, 20))
})

test_that("read_rangeData respects named GWAS columns in any order", {
    df <- data.frame(
        p = c(0.1, 0.2),
        chrom = c("chr1", "chr2"),
        pos = c(100, 200),
        snp = c("rs1", "rs2")
    )
    out <- plotgardener:::read_rangeData(
        data = df, assembly = "hg38", type = "GWAS"
    )
    expect_equal(colnames(out)[1:2], c("chrom", "pos"))
    expect_equal(out$chrom, c("chr1", "chr2"))
    expect_equal(out$pos, c(100, 200))
})
