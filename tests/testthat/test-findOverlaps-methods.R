# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###
context("findOverlaps methods")

test_that("GTuples,GTuples overlaps", {
    ## empty
    hits <- findOverlaps(gt0, gt0)
    expect_true(inherits(hits, "Hits"))
    expect_equal(length(hits), 0)
    ## identical 1
    hits <- findOverlaps(gt1, gt1)
    expect_true(inherits(hits, "Hits"))
    expect_equal(length(hits), 10)
    expect_equal(hits@queryHits, 1:10)
    expect_equal(hits@queryHits, hits@subjectHits)
    ## identical 2
    tmp2 <- GTuples(seqnames = c("chr1", "chr1"),
                    tuples = matrix(c(1L,2L,3L,4L), nrow=2),
                    strand = "+")
    hits = findOverlaps(tmp2, tmp2, type="equal")
    expect_equal(length(hits), 2)
    expect_equal(hits@queryHits, 1:2)
    hits = findOverlaps(tmp2, tmp2, type="any")
    expect_equal(length(hits), 4)
    expect_equal(hits@queryHits, c(1,1,2,2))
    expect_equal(hits@subjectHits, c(1,2,1,2))
    ## identical 3
    tmp3 <- GTuples(seqnames = c("chr1", "chr1"),
                    tuples = matrix(c(1L,2L,3L,4L,5L,6L), nrow=2),
                    strand = "+")
    hits = findOverlaps(tmp3, tmp3, type="equal")
    expect_equal(length(hits), 2)
    expect_equal(hits@queryHits, 1:2)
    hits = findOverlaps(tmp3, tmp3, type="any")
    expect_equal(length(hits), 4)
    expect_equal(hits@queryHits, c(1,1,2,2))
    expect_equal(hits@subjectHits, c(1,2,1,2))
    
})
