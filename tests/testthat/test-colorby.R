test_that("colorby object created", {
    object1 <- colorby(column = "foo")
    object2 <- colorby(column = "foo", range = c(0, 1))

    expect_equal(class(object1), "bb_colorby")
    expect_equal(class(object2), "bb_colorby")
})

test_that("mapped colors", {
    ## Numeric vector color mappings
    numVector <- c(1, 2, 3, 2, 1)
    colorVector <- bbMapColors(vector = numVector,
                palette = colorRampPalette(c("blue", "red", "black")),
                range = NULL)
    expect_setequal(colorVector, c("#0000FF", "#FF0000", "#050000", "#FF0000", 
                                "#0000FF"))
    colorVector <- bbMapColors(vector = numVector,
                                palette = colorRampPalette(c("blue", "red", "black")),
                                range = c(1, 2))
    
    expect_setequal(colorVector, c("#0000FF", "#050000", "#050000", 
                                   "#050000", "#0000FF"))
    
    
    ## Factor vector where levels are defaulted to: cat, dog, pig
    factVector <- factor(c("dog", "cat", "pig", "cat"))
    ## cat mapped to blue, dog mapped to red, pig mapped to black
    colorVector <- bbMapColors(vector = factVector,
                                palette = colorRampPalette(c("blue", "red", "black")))
    expect_mapequal(colorVector[levels(factVector)], c("dog" = "#FF0000", "cat" = "#0000FF",
                                      "pig" = "#000000"))
    
    ## Factor vector where levels are: dog, cat, pig
    factVector <- factor(c("dog", "cat", "pig", "cat"))
    levels(factVector) <- c("dog", "cat", "pig")
    ## dog mapped to blue, cat mapped to red, pig mapped to black
    colorVector <- bbMapColors(vector = factVector,
                                palette = colorRampPalette(c("blue", "red", "black")))
    expect_mapequal(colorVector[levels(factVector)],
                    c("dog" = "#0000FF", "cat" = "#FF0000", "pig" = "#000000"))
    
    
    ## Subset of factor vector with four levels, but only three values present
    factVector <- factor(c("dog", "cat", "pig", "cat", "giraffe"))
    factVector <- factVector[1:4]
    # levels are cat, dog, giraffe, pig
    ## cat mapped to blue, dog mapped to red, giraffe mapped to black, pig mapped to yellow
    colorVector <- bbMapColors(vector = factVector,
                                palette = colorRampPalette(c("blue", "red", "black", "yellow")))
    expect_mapequal(colorVector[c("dog", "cat", "pig")],
                    c("cat" = "#0000FF", "dog" = "#FF0000", "pig" = "#FFFF00"))
    
    ## Non-numeric, non-factor vector
    nonnumVector <- c("dog", "cat", "pig", "pig")
    # cat mapped to blue, dog mapped to red, pig mapped to black
    colorVector <- bbMapColors(vector = nonnumVector,
                                palette = colorRampPalette(c("blue", "red", "black")))
    expect_mapequal(colorVector[c("dog", "cat", "pig")],
                    c("cat" = "#0000FF", "dog" = "#FF0000", "pig" = "#000000"))
    
})