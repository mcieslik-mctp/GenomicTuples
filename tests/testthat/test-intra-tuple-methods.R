# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###
context("GTuples intra-tuple methods")

test_that("GTuples shift", {
    expect_equal(shift(gt0), gt0)
    expect_equal(shift(gt0, 2), gt0)
    expect_equal(shift(gt1), gt1)
    # proper dispatch to GRanges
    expect_equal(granges(shift(gt1, 2)), shift(granges(gt1), 2))
    expect_equal(granges(shift(gt2, 2)), shift(granges(gt2), 2))
    expect_equal(granges(shift(gt3, 2)), shift(granges(gt3), 2))
    ## is the internal field updated?
    expect_equal(gt3@internalPos + 2, shift(gt3, 2)@internalPos)
    expect_equal(gt4@internalPos + 2, shift(gt4, 2)@internalPos)

})

test_that("GTuplesList shift", {
    expect_equal(shift(gtl0), gtl0)
    expect_equal(shift(gtl0, 2), gtl0)
    expect_equal(shift(gtl1), gtl1)
    # proper dispatch to GRanges
    expect_equal(granges(unlist(shift(gtl1,2))), shift(granges(unlist(gtl1)),2))
    expect_equal(granges(unlist(shift(gtl2,2))), shift(granges(unlist(gtl2)),2))
    expect_equal(granges(unlist(shift(gtl3,2))), shift(granges(unlist(gtl3)),2))
    ## is the internal field updated?
    expect_equal(gtl3@unlistData@internalPos + 2, shift(gtl3, 2)@unlistData@internalPos)
    expect_equal(gtl4@unlistData@internalPos + 2, shift(gtl4, 2)@unlistData@internalPos)
})

test_that("GTuples narrow", {
    ## GRanges can be narrowed down to zero-width
    expect_true(all(width(narrow(granges(gt2), 3)) == 0))
    ## if the narrow method is inherited from GRanges the following will
    ## violate GTuples constraints pos1 < pos0
    expect_error(narrow(gt2, 3), "GTuples do not currently support the 'flank' method.")
})

test_that("GTuplesList narrow", {
    expect_error(narrow(gtl2, 3), "GTuplesList do not currently support the 'flank' method.")
})
