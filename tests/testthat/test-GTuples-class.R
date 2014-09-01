# NB: Several objects used in testing are defined in 
# tests/testthat/helper-make-test-data.R

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###
context("GTuples validity methods")

test_that(".valid.GTuples.pos for 1-tuples", {
  expect_that(GTuples('chr1', tuples = matrix(12:1, ncol = 1)), 
               not(throws_error()))
})
test_that(".valid.GTuples.pos for 2-tuples", {
    expect_error(GTuples('chr1', tuples = matrix(12:1, ncol = 2)), 
               "negative widths are not allowed")
    expect_error(GTuples('chr1', tuples = cbind(11:20, 11:20)), 
                 "positions in each tuple must be sorted")
    expect_error(GTuples('chr1', tuples = cbind(11:20, c(12:20, 20L))), 
                 "positions in each tuple must be sorted")
})
test_that(".valid.GTuples.pos for m-tuples (m > 2)", {
  expect_error(GTuples('chr1', tuples = matrix(12:1, ncol = 3)), 
               "negative widths are not allowed")
  expect_error(GTuples('chr1', tuples = cbind(11:20, 1:10, 31:40)), 
               "positions in each tuple must be sorted")
  expect_error(GTuples('chr1', tuples = cbind(11:20, c(12:20, 20L), 31:40)), 
               "positions in each tuple must be sorted")
  expect_error(GTuples('chr1', tuples = matrix(12:1, ncol = 4)), 
               "negative widths are not allowed")
  expect_error(GTuples('chr1', tuples = cbind(11:20, 1:10, 31:40, 41:50)), 
               "positions in each tuple must be sorted")
  expect_error(GTuples('chr1', tuples = cbind(11:20, 31:40, 1:10, 41:50)), 
               "positions in each tuple must be sorted")
  expect_error(GTuples('chr1', tuples = cbind(11:20, c(12:20, 20L), 31:40, 
                                              41:50)), 
               "positions in each tuple must be sorted")
})

test_that(".valid.GTuples.mcols", {
  c("seqnames", "ranges", "strand",
    "seqlevels", "seqlengths", "isCircular",
    #"genome",
    "start", "end", "width", "element",
    "tuples", "internalPos", "size")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       seqnames = letters[1:10]), 
               "formal argument \"seqnames\" matched by multiple")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       ranges = letters[1:10]), 
               "formal argument \"ranges\" matched by multiple")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       strand = letters[1:10]), 
               "formal argument \"strand\" matched by multiple")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       strand = letters[1:10]), 
               "formal argument \"strand\" matched by multiple")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       seqlevels = letters[1:10]), 
               "names of metadata columns cannot be one of")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       seqlengths = letters[1:10]), 
               "length of supplied 'seqlengths' must equal the number of")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       isCircular = letters[1:10]), 
               "names of metadata columns cannot be one of")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       start = letters[1:10]), 
               "names of metadata columns cannot be one of")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       end = letters[1:10]), 
               "names of metadata columns cannot be one of")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       width = letters[1:10]), 
               "names of metadata columns cannot be one of")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       element = letters[1:10]), 
               "names of metadata columns cannot be one of")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       tuples = letters[1:10]), 
               "formal argument \"tuples\" matched by multiple")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       internalPos = letters[1:10]), 
               "names of metadata columns cannot be one of")
  expect_error(GTuples(seqnames = 'chr1', tuples = matrix(1:10), strand = '+',
                       size = letters[1:10]), 
               "names of metadata columns cannot be one of")
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###
context("GTuples constructor")

test_that("GTuples constructor returns a valid GTuples object when m = 0", {
  expect_true(validObject(new("GTuples")))
  expect_true(validObject(gt0))
})

test_that("GTuples constructor returns a GTuples object when m = 1", {
  expect_true(validObject(gt1))
})

test_that("GTuples constructor returns a GTuples object when m = 2", {
  expect_true(validObject(gt4))
})

test_that("GTuples constructor returns a GTuples object when m >= 3", {
  expect_true(validObject(gt3))
  expect_true(validObject(gt4))
})

test_that("GTuples constructor returns warnings on unexpected input", {
  expect_warning(GTuples('chr1', tuples = matrix(c(1.1, 2, 3), ncol = 1)), 
                 "Converting 'tuples' to integer mode")
})

test_that("GTuples constructor returns errors on bad input", {
  expect_error(GTuples('chr1', tuples = 1:10), 
               "'tuples' must be an integer matrix")
  expect_error(GTuples('chr1', tuples = as.matrix(letters)), 
               "'tuples' must be an integer matrix")
  expect_error(GTuples('chr1', tuples = matrix(c(1, NA), ncol = 1)), 
              "'NA' detected in 'tuples'")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###
context("GTuples coercion")

test_that("as.data.frame works", {
  expect_error(as.data.frame(gt0), "arguments imply differing number of rows")
  expect_identical(as.data.frame(gt1), 
                   data.frame(seqnames = factor(as.character(seqnames(gt1))),
                              pos1 = start(gt1),
                              strand = factor(as.character(strand(gt1)), 
                                              levels = c('+', '-', '*')),
                              score = mcols(gt1)$score
                   )
  )
  expect_identical(as.data.frame(gt2), 
                   data.frame(seqnames = factor(as.character(seqnames(gt4))),
                              pos1 = start(gt1),
                              pos2 = end(gt2),
                              strand = factor(as.character(strand(gt2)), 
                                              levels = c('+', '-', '*')),
                              score = mcols(gt2)$score
                   )
  )
  expect_identical(as.data.frame(gt3), 
                   data.frame(seqnames = factor(as.character(seqnames(gt3))),
                              pos1 = start(gt3),
                              pos2 = gt3@internalPos,
                              pos3 = end(gt3),
                              strand = factor(as.character(strand(gt3)), 
                                              levels = c('+', '-', '*')),
                              score = mcols(gt1)$score
                   )
  )
  expect_identical(as.data.frame(gt4), 
                   data.frame(seqnames = factor(as.character(seqnames(gt4))),
                              pos1 = start(gt4),
                              pos2 = gt4@internalPos[, 1],
                              pos3 = gt4@internalPos[, 2],
                              pos4 = end(gt4),
                              strand = factor(as.character(strand(gt4)), 
                                              levels = c('+', '-', '*')),
                              score = mcols(gt4)$score
                   )
  )
})

test_that("granges works", {
  expect_identical(granges(gt0), GRanges())
  expect_identical(granges(gt1, use.mcols = FALSE), 
                   GRanges(seqnames(gt1), IRanges(start(gt1), end(gt1)),
                           strand(gt1), 
                           seqinfo = seqinfo(gt1)))
  expect_identical(granges(gt1, use.mcols = TRUE), 
                   GRanges(seqnames(gt1), IRanges(start(gt1), end(gt1)),
                           strand(gt1), score = mcols(gt1)$score, 
                           seqinfo = seqinfo(gt1)))
  expect_identical(granges(gt4, use.mcols = FALSE), 
                   GRanges(seqnames(gt4), IRanges(start(gt4), end(gt4)),
                           strand(gt4), 
                           seqinfo = seqinfo(gt4)))
  expect_identical(granges(gt4, use.mcols = TRUE), 
                   GRanges(seqnames(gt4), IRanges(start(gt4), end(gt4)),
                           strand(gt4), score = mcols(gt4)$score, 
                           seqinfo = seqinfo(gt4)))
  expect_identical(granges(gt3, use.mcols = FALSE), 
                   GRanges(seqnames(gt3), IRanges(start(gt3), end(gt3)),
                           strand(gt3), 
                           seqinfo = seqinfo(gt3)))
  expect_identical(granges(gt3, use.mcols = TRUE), 
                   GRanges(seqnames(gt3), IRanges(start(gt3), end(gt3)),
                           strand(gt3), score = mcols(gt3)$score, 
                           seqinfo = seqinfo(gt3)))
  expect_identical(granges(gt4, use.mcols = FALSE), 
                   GRanges(seqnames(gt4), IRanges(start(gt4), end(gt4)),
                           strand(gt4), 
                           seqinfo = seqinfo(gt4)))
  expect_identical(granges(gt4, use.mcols = TRUE), 
                   GRanges(seqnames(gt4), IRanges(start(gt4), end(gt4)),
                           strand(gt4), score = mcols(gt4)$score, 
                           seqinfo = seqinfo(gt4)))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning
###

context("GTuples updating and cloning")

test_that("update works on all relevant slots", {
  x <- gt1
  x <- update(x, seqnames = rev(seqnames(x)))
  expect_identical(x, GTuples(rev(seqnames(gt1)), tuples(gt1), strand(gt1),
                              score = mcols(gt1)$score, seqinfo = seqinfo(gt1)))
  x <- gt1
  x <- update(x, ranges = rev(ranges(x)))
  expect_identical(x, GTuples(seqnames(gt1), as.matrix(rev(tuples(gt1))), 
                              strand(gt1), score = mcols(gt1)$score, 
                              seqinfo = seqinfo(gt1)))
  x <- gt1
  x <- update(x, strand = rev(strand(x)))
  expect_identical(x, GTuples(seqnames(gt1), tuples(gt1), rev(strand(gt1)),
                              score = mcols(gt1)$score, seqinfo = seqinfo(gt1)))
  x <- gt1
  x <- update(x, elementMetadata = DataFrame(score = Rle(0L, 10)))
  expect_identical(x, GTuples(seqnames(gt1), tuples(gt1), strand(gt1),
                              score = Rle(0L, 10), seqinfo = seqinfo(gt1)))
  x <- gt1
  seqinfo <- Seqinfo(seqnames = c("chr1", "chr2", "chr3"), 
                     seqlengths = c(10000L, 20000L, 15000L), 
                     isCircular = c(NA, NA, NA), 
                     genome = c("mock1", "mock1", "mock1"))
  x <- update(x, seqinfo = seqinfo)
  expect_identical(x, GTuples(seqnames(gt1), tuples(gt1), strand(gt1),
                              score = mcols(gt1)$score, seqinfo = seqinfo))
  x <- gt3
  x <- update(x, ranges = IRanges(start(x) + 10L, end(x) + 10L), 
                                  internalPos = x@internalPos + 10L)
  expect_identical(x, GTuples(seqnames(gt3), unname(tuples(gt3)) + 10L, 
                              strand(gt3), score = mcols(gt3)$score, 
                              seqinfo = seqinfo(gt3)))
  x <- gt4
  x <- update(x, ranges = IRanges(start(x) + 10L, end(x) + 10L), 
              internalPos = x@internalPos + 10L)
  expect_identical(x, GTuples(seqnames(gt4), unname(tuples(gt4)) + 10L, 
                              strand(gt4), score = mcols(gt4)$score, 
                              seqinfo = seqinfo(gt4)))
  
  
  
})

test_that("clone works", {
  
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###
test_that("GRanges inherited getters work", {
    expect_identical(seqnames(gt1), gt1@seqnames)
    expect_identical(ranges(gt2), gt2@ranges)
    expect_identical(strand(gt3), gt3@strand)
    expect_identical(seqinfo(gt4), gt4@seqinfo)
    expect_identical(seqinfo(gt4), gt4@seqinfo)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###
test_that("inherited split works", {
    ## by integer
    gt2_s = split(gt2, 1:10)
    expect_equal(length(gt2_s), 10)
    expect_true(inherits(gt2_s, class(gtl2)))
    ## by Rle
    gt2_s = split(gt2, seqnames(gt2))
    expect_equal(length(gt2_s), 3)
    expect_true(inherits(gt2_s, class(gtl2)))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Concatenating
###
test_that("concatenation works", {
    expect_equal(c(gt3[1:5], gt3[6:10]), gt3)
    expect_error(c(gt3, granges(gt3)))
    expect_error(c(gt3, gt4))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

context("GTuples accessors")

test_that("size works", {
  expect_identical(size(gt0), NA_integer_)
  expect_identical(size(gt1), 1L)
  expect_identical(size(gt2), 2L)
  expect_identical(size(gt3), 3L)
  expect_identical(size(gt4), 4L)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### IPD
###
test_that("IPD works", {
    expect_error(IPD(gt1))
    expect_equal(IPD(gt2), matrix(1, nrow=10, ncol=1))
    expect_equal(IPD(gt3), matrix(1, nrow=10, ncol=2))
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine and split
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###



