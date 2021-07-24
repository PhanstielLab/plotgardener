test_that("parseParams", {
    
    test <- function(params = NULL, alpha, hic, altchrom = 1, altchromend = NULL, 
                     arch = "blah", ...) {
        ## Access default arguments of parent function like this:
        defArgs <- formals(eval(match.call()[[1]]))
        ## Access user-declared arguments of parent function like this:
        decArgs <- lapply(match.call()[-1], eval)
        ## Call parseParams like this:
        x <- BentoBox:::parseParams(params = params, defaultArgs = defArgs, 
                                    declaredArgs = decArgs,
                                    class = "object")
        ## Crete internal object like this:
        object <- structure(
            .Data = x,
            class = "object"
        )
        return(object)
    }
    
    ## All default parameters
    expect_equal(test(),
                 structure(
                     .Data = list("alpha" = NULL,
                     "hic" = NULL,
                     "altchrom" = 1,
                     "altchromend" = NULL,
                     "arch" = "blah"),
                     class = "object"
                 ))
    
    ## Default argument replaced
    expect_equal(test(arch = "not_blah"),
                 structure(
                     .Data = list("alpha" = NULL,
                                  "hic" = NULL,
                                  "altchrom" = 1,
                                  "altchromend" = NULL,
                                  "arch" = "not_blah"),
                     class = "object"
                 ))
    
    ## Multiple default replacements
    expect_equal(test(alpha = 10, altchrom = 2),
                 structure(
                     .Data = list("alpha" = 10,
                                  "hic" = NULL,
                                  "altchrom" = 2,
                                  "altchromend" = NULL,
                                  "arch" = "blah"),
                     class = "object"
                 ))
    ## Overriding bb_params parameter
    expect_equal(test(params = bb_params(alpha = 100), arch = "not_blah", alpha = 10),
                 structure(
                     .Data = list("alpha" = 10,
                                  "hic" = NULL,
                                  "altchrom" = 1,
                                  "altchromend" = NULL,
                                  "arch" = "not_blah"),
                     class = "object"
                 ))
    
    ## Warning without a bb_params object
    expect_warning(test(params = c(hic = "no")))
})

test_that("setGP", {
    
    ## Empty gpar
    expect_equal(BentoBox:::setGP(
        gpList = gpar(),
        params = NULL
    ), gpar())
    
    ## One gpar parameter replaced
    expect_equal(BentoBox:::setGP(
        gpList = gpar(),
        params = NULL,
        lwd = 2
    ), gpar(lwd = 2))

    ## Multiple gpar parameters replaced
    expect_equal(BentoBox:::setGP(
        gpList = gpar(),
        params = NULL,
        lwd = 2, lty = 2
    ), gpar(lwd = 2, lty = 2))
    
})

test_that("bb_regionErrors", {
    
    ## No errors
    expect_silent(BentoBox:::bb_regionErrors(chromstart = 1000000,
                               chromend = 2000000))
    ## 0 bp wide region error
    expect_error(BentoBox:::bb_regionErrors(chromstart = 1000000,
                                             chromend = 1000000))
    ## One null chromstart/chromend
    expect_error(BentoBox:::bb_regionErrors(chromstart = 1000000,
                                            chromend = NULL))
    
    ## chromstart bigger than chromend
    expect_error(BentoBox:::bb_regionErrors(chromstart = 2000000,
                                            chromend = 1000000))
})

test_that("bb_justConversion", {
    
    ## left, top should be 0, 1
    expect_equal(BentoBox:::bb_justConversion(just = c("left", "top")),
                 c(0, 1))
    
    ## top, left should be 0, 1
    expect_equal(BentoBox:::bb_justConversion(just = c("top", "left")),
                 c(0, 1))
    
    ## Invalid just option error
    expect_error(BentoBox:::bb_justConversion(just = "invalid"))
    
})

test_that("parse_bbAssembly", {
    
    ## Invalid default assembly
    expect_error(BentoBox:::parse_bbAssembly(assembly = "none"))
    
    ## Default hg19 assembly
    expect_equal(BentoBox:::parse_bbAssembly(assembly = "hg19"),
                 structure(.Data = list("Genome" = "hg19",
                                        "TxDb" = "TxDb.Hsapiens.UCSC.hg19.knownGene",
                                        "OrgDb" = "org.Hs.eg.db",
                                        "gene.id.column" = "ENTREZID",
                                        "display.column" = "SYMBOL",
                                        "BSgenome" = "BSgenome.Hsapiens.UCSC.hg19"),
                           class = "bb_assembly"))
    
    ## Input bb_assembly object just returns itself
    assemblyobject <- bb_assembly(Genome = "testing",
                                  TxDb = "TxDb",
                                  OrgDb = "OrgDb")
    
    expect_equal(BentoBox:::parse_bbAssembly(assembly = assemblyobject),
                 structure(.Data = list("Genome" = "testing",
                                        "TxDb" = "TxDb",
                                        "OrgDb" = "OrgDb",
                                        "gene.id.column" = "ENTREZID",
                                        "display.column" = "SYMBOL"),
                           class = "bb_assembly"))
})

test_that("defaultUnits", {
    
    ## No coordinates or dimensions need to be converted
    testObject <- list("x" = unit(1, "inches"),
                       "y" = unit(2, "inches"),
                       "width" = unit(3, "npc"),
                       "height" = unit(4, "cm"))
    expect_equal(BentoBox:::defaultUnits(object = testObject, 
                                         default.units = "inches"),
                 testObject)
    expect_equal(BentoBox:::misc_defaultUnits(value = unit(1, "inches"),
                                              name = "x0",
                                              default.units = "npc"),
                 unit(1, "inches"))
    
    ## Some coordinates converted to default units
    testObject <- list("x" = 1,
                       "y" = unit(2, "inches"),
                       "width" = 3,
                       "height" = unit(4, "cm"))
    expect_equal(BentoBox:::defaultUnits(object = testObject, 
                                         default.units = "inches"),
                 list("x" = unit(1, "inches"),
                      "y" = unit(2, "inches"),
                      "width" = unit(3, "inches"),
                      "height" = unit(4, "cm")))
    expect_equal(BentoBox:::misc_defaultUnits(value = 1,
                                              name = "x0",
                                              default.units = "npc"),
                 unit(1, "npc"))
    
    ## below y-coordinate errors
    testObject <- list("x" = 1,
                       "y" = "plb",
                       "width" = 3,
                       "height" = unit(4, "cm"))
    expect_error(BentoBox:::defaultUnits(object = testObject,
                                         default.units = "inches"))
    expect_error(BentoBox:::misc_defaultUnits(value = "plb",
                                              name = "y0",
                                              default.units = "npc"))
})