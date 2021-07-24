test_that("convert_gpath", {
    
    testGrob <- rectGrob()
    ## Grob becomes a gpath
    expect_true(is(BentoBox:::convert_gpath(testGrob), "gPath"))
})

test_that("check_placement", {
    
    testObject <- list("x" = NULL,
                       "y" = NULL)
    
    attr(testObject, "plotted") <- FALSE
    ## doesn't check if not plotted
    expect_silent(BentoBox:::check_placement(testObject))
    
    ## error one null x or y
    testObject <- list("x" = unit(3, "inches"),
                       "y" = NULL)
    
    attr(testObject, "plotted") <- TRUE
    expect_error(BentoBox:::check_placement(testObject))
    
    ## error for 0 width
    testObject <- list("x" = unit(2, "inches"),
                       "y" = unit(0.5, "npc"),
                        "width" = unit(3, "inches"),
                       "height" = unit(0, "inches"))
    
    attr(testObject, "plotted") <- TRUE
    expect_error(BentoBox:::check_placement(testObject))
    
    ## No errors with all checks
    bb_pageCreate()
    testObject <- list("x" = unit(2, "inches"),
                       "y" = unit(0.5, "npc"),
                       "width" = unit(3, "inches"),
                       "height" = unit(3, "inches"))
    attr(testObject, "plotted") <- TRUE
    expect_silent(BentoBox:::check_placement(testObject))
    
})

test_that("assignRows", {
    library(BentoBoxData)
    data("IMR90_ChIP_CTCF_reads")
    
    ## set seed
    set.seed(nrow(IMR90_ChIP_CTCF_reads))
    
    IMR90_ChIP_CTCF_reads <- IMR90_ChIP_CTCF_reads[
        which(IMR90_ChIP_CTCF_reads[, "chrom"] == "chr21"
              & IMR90_ChIP_CTCF_reads[, "start"] <=
                  29074000 & IMR90_ChIP_CTCF_reads[, "end"] >=
                  29073000), ]
    
    ## Message without limitLabel 
    expect_warning(BentoBox:::assignRows(data = IMR90_ChIP_CTCF_reads[, c(2,3)],
                          maxRows = 3,
                          wiggle = 10,
                          rowCol = 2, 
                          limitLabel = FALSE,
                          gTree = NULL),
                   regexp = "Not enough plotting space for all provided elements.")
    
    ## Message with limitLabel
    assign("pileup_grobs", gTree(), envir = bbEnv)
    expect_warning(BentoBox:::assignRows(data = IMR90_ChIP_CTCF_reads[, c(2,3)],
                                         maxRows = 3,
                                         wiggle = 10,
                                         rowCol = 2, 
                                         limitLabel = TRUE,
                                         gTree = "pileup_grobs"))
})