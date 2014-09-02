# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### comparison
###
context("GTuples comparison methods")

test_that("GTuples,GTuples comparison errors work", {
    # empty fails
    expect_error(.GTuples.compare(gt0, gt0), "Cannot compare empty GTuples.")
    expect_error(gt0 == gt0, "Cannot compare empty GTuples.")
    expect_error(.GTuples.compare(gt0, gt1), "Cannot compare 'GTuples' with different length.")
    expect_error(.GTuples.compare(gt1, gt2), "Cannot compare 'GTuples' objects of different 'size'.")
    expect_error(.GTuples.compare(gt2, gt3), "Cannot compare 'GTuples' objects of different 'size'.")
    expect_error(.GTuples.compare(gt3, gt4), "Cannot compare 'GTuples' objects of different 'size'.")
    expect_error(.GTuples.compare(gt3[1], gt4[1]), "Cannot compare 'GTuples' objects of different 'size'.")

    # switching chromosome names
    seqinfo <- Seqinfo(paste0("chr", 3:1), c(1000, 2000, 1500), NA, "mock1")
    gt3_fake <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"),
                     c(1, 3, 2, 4)),
                   tuples = matrix(c(1:10, 2:11, 3:12), ncol = 3),
                   strand = Rle(strand(c("-", "+", "*", "+", "-")),
                     c(1, 2, 2, 3, 2)),
                   score = 1:10, seqinfo = seqinfo)
    # error is too long and detailed to capture
    expect_error(gt3 == gt3_fake)
    expect_error(granges(gt3) == granges(gt3_fake))
    # Granges fails the same way
    expect_equal(
      tryCatch(gt3 == gt3_fake, error=function(e) as.character(e)),
      tryCatch(granges(gt3) == granges(gt3_fake), error=function(e) as.character(e))
    )
    
})

test_that("GTuples,GTuples internal equal works", {
    # length 1
    expect_equal(.GTuples.compare(gt1[1], gt1[1]), 0)
    expect_equal(.GTuples.compare(gt2[1], gt2[1]), 0)
    expect_equal(.GTuples.compare(gt3[1], gt3[1]), 0)
    expect_equal(.GTuples.compare(gt4[1], gt4[1]), 0)
    
})
