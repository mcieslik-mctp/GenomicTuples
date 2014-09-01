### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GTuples objects used in tests
###
### Copied from examples in man/GTuples-class.Rd

seqinfo <- Seqinfo(paste0("chr", 1:3), c(1000, 2000, 1500), NA, "mock1")
# Empty GTuples object
gt0 <- GTuples()
# 1-tuples
gt1 <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"),
                              c(1, 3, 2, 4)),
               tuples = matrix(c(1:10), ncol = 1),
               strand = Rle(strand(c("-", "+", "*", "+", "-")),
                            c(1, 2, 2, 3, 2)),
               score = 1:10, seqinfo = seqinfo)
# 2-tuples
gt2 <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"),
                              c(1, 3, 2, 4)),
               tuples = matrix(c(1:10, 2:11), ncol = 2),
               strand = Rle(strand(c("-", "+", "*", "+", "-")),
                            c(1, 2, 2, 3, 2)),
               score = 1:10, seqinfo = seqinfo)
# 3-tuples
gt3 <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"),
                              c(1, 3, 2, 4)),
               tuples = matrix(c(1:10, 2:11, 3:12), ncol = 3),
               strand = Rle(strand(c("-", "+", "*", "+", "-")),
                            c(1, 2, 2, 3, 2)),
               score = 1:10, seqinfo = seqinfo)
# 4-tuples
gt4 <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"),
                              c(1, 3, 2, 4)),
               tuples = matrix(c(1:10, 2:11, 3:12, 4:13), ncol = 4),
               strand = Rle(strand(c("-", "+", "*", "+", "-")),
                            c(1, 2, 2, 3, 2)),
               score = 1:10, seqinfo = seqinfo)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GTuplesList objects used in tests
###
### Copied from examples in man/GTuplesList-class.Rd

gtl0 <- GTuplesList(A = gt0, B = gt0)
gtl1 <- GTuplesList(A = gt1[1:5], B = gt1[6:10])
gtl2 <- GTuplesList(A = gt2[1:5], B = gt2[6:10])
gtl3 <- GTuplesList(A = gt3[1:5], B = gt3[6:10])
gtl4 <- GTuplesList(A = gt4[1:5], B = gt4[6:10])
