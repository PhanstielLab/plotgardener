test_that("pgParams object created", {
    object <- pgParams(chrom = "chr1")
    expect_equal(class(object), "pgParams")
})

test_that("pgParams object concatenates", {
    object1 <- pgParams(chrom = "chr1")
    object2 <- pgParams(chromstart = 1000000)
    object3 <- pgParams(chrom = "chr1", chromstart = 1000000)

    expect_equal(class(c(object1, object2)), "pgParams")
    expect_equal(c(object1, object2), object3)
    expect_warning(c(object1, "gene"))
})

test_that("pgParams parameter parsing", {
    
    testFun <- function(data, chrom, chromstart = NULL, chromend = NULL,
                        assembly = "hg19", 
                        params = NULL, ...){
        newObject <- plotgardener:::parseParams(params = params,
                                            defaultArgs = formals(eval(match.call()[[1]])),
                                            declaredArgs = lapply(match.call()[-1], eval.parent,
                                                                  n = 2),
                                            class = "test_object")
        return(newObject)
    }
    
    ## All default parameters
    expect_mapequal(
        testFun(),
        list(data = NULL,
             chrom = NULL,
             chromstart = NULL,
             chromend = NULL,
             assembly = "hg19"
        )
    )
    
    ## Parameter replaces default
    expect_mapequal(
        testFun(chromstart = 1000000),
        list(data = NULL,
             chrom = NULL,
             chromstart = 1000000,
             chromend = NULL,
             assembly = "hg19"
        )
    )
    
    ## pgParams replaces default
    params <- pgParams(assembly = "hg38")
    expect_mapequal(
        testFun(params = params),
        list(data = NULL,
             chrom = NULL,
             chromstart = NULL,
             chromend = NULL,
             assembly = "hg38"
        )
    )
    
    ## Input overrides pgParams
    expect_mapequal(
        testFun(params = params, assembly = "mm9"),
        list(data = NULL,
             chrom = NULL,
             chromstart = NULL,
             chromend = NULL,
             assembly = "mm9"
        )
    )

    ## Multiple replacements
    expect_mapequal(
        testFun(params = params, assembly = "mm9", chromstart = 1000000),
        list(data = NULL,
             chrom = NULL,
             chromstart = 1000000,
             chromend = NULL,
             assembly = "mm9"
        )
    )
    
})
