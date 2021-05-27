test_that("colorby object created", {
    object1 <- colorby(column = "foo")
    object2 <- colorby(column = "foo", range = c(0, 1))

    expect_equal(class(object1), "bb_colorby")
    expect_equal(class(object2), "bb_colorby")
})
