# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### apply methods
###
context("GTuplesList apply methods")

test_that("GTuplesList endoapply works", {
    expect_identical(gtl0, endoapply(gtl0, function(x) {x}))
    expect_identical(gtl1, endoapply(gtl1, function(x) {x}))
    expect_identical(gtl2, endoapply(gtl2, function(x) {x}))
    expect_identical(gtl4, endoapply(gtl4, function(x) {x}))
})
