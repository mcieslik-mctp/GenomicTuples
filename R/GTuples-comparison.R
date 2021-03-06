### =========================================================================
### Comparing and ordering genomic tuples
### -------------------------------------------------------------------------
### I. UNIQUE AND DUPLICATED ELEMENTS WITHIN A GTuples OBJECT
### ---------------------------------------------------------------
### Two elements of a GTuples object (i.e. two genomic tuples) are
### considered equal iff they are on the same underlying sequence and strand,
### and have the same positions. duplicated() and unique() on a
### GTuples object are conforming to this.
###
### II. ORDERING THE ELEMENTS OF A GTuples OBJECT
### ---------------------------------------------------
### The "natural order" for the elements of a GTuples object is to order
### them (a) first by sequence level, (b) then by strand, (c) then by pos_{1}, 
### ..., pos_{m}. This way, the space of genomic tuples is totally ordered.
### order(), sort(), and rank() on a GTuples object are using this
### "natural order".
###
### III. ELEMENT-WISE (AKA "PARALLEL") COMPARISON OF 2 GTuples OBJECTS
### ------------------------------------------------------------------------
### We want the "==", "!=", "<=", ">=", "<" and ">" operators between 2
### GTuples objects to be compatible with the "natural order" defined
### previously. Defining those operators when the 2 objects have *identical*
### seqlevels() is straighforward but we can in fact extend this comparison
### to the following situation:
###   (A) 'e1' and 'e2' have compatible sets of underlying sequences, that is,
###       'seqinfo(e1)' and 'seqinfo(e2)' can be merged.
###   (B) 'seqlevels(e1)' and 'seqlevels(e2)' are in the same order. Note that
###       (A) guarantees that the seqlevels of one is a subset of the seqlevels 
###       of the other. (B) is saying that this subset should be a subsequence.
### Pre-comparison step: if (A) and (B) are satisfied, then the 2 seqinfo() are
### merged and the seqlevels() of the result is assigned back to each object
### to compare. This is a way to have 2 objects with identical seqlevels()
### before the comparison can actually be performed and meaningful.
### The reason (B) is required for the pre-comparison step is because we want
### this step to preserve the original order of the seqlevels() in *both*
### objects. Without this precaution, the expected anti-symetric property of
### some operators would not be satisfied e.g. 'any(e1 < e2 & e2 < e1)' could
### be TRUE.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compare() and related methods.
###
### compare() is based on the method defined for GRanges objects but cannot 
### explicity inherit the method due to the internalPos slot in GTuples 
### objects. However, unlike the compare() method defined for GRanges, the 
### compare() method for GTuples requires that both x and y have the same 
### length. I also define the element wise (aka "parallel") operators '<=' and 
### '=='. The other element wise operators (`!=`, `>=`, `<`, `>`) work 
### out-of-the-box on GTuples objects via inheritance from GRanges -> Vector.

# .GTuples.compare is adapted from GenomicRanges:::.GenomicRanges.compare.
# On 2 GTuples objects, it returns one of the 3 codes: If the first tuple in 
# the pair is "<" than the second tuple then the return value for that element 
# is < 0, if the first tuple in the pair is "==" the second tuple then the 
# return value is 0, and if the first tuple is ">" that the second tuple then 
# the return value is > 0.
.GTuples.compare <- function(x, y) {
  
  # Different to comparing GRanges, I only allow comparison if x, y have same
  # length.
  # shortened error message because a long error trigger line formatting
  # that breaks the testthat error parser.
  if (length(x) != length(y)){
    stop("Cannot compare 'GTuples' with different length.")
  }

  # This is where .GTuples.compare really differs from .GenomicRanges.compare
  # NOTE: moved this up because the next 'if' will fail on NA != NA
  if (is.na(size(x)) || is.na(size(y))) {
    stop("Cannot compare empty GTuples.")
  }
  
  # Check 'size' is identical
  if (size(x) != size(y)) {
    stop("Cannot compare 'GTuples' objects of different 'size'.")
  }
  
  # Pre-comparison step (see above for details).
  # merge() will fail if 'x' and 'y' don't have compatible underlying
  # sequences.
  seqinfo <- merge(seqinfo(x), seqinfo(y))
  seqlevels <- seqlevels(seqinfo)
  if (any(diff(match(seqlevels(y), seqlevels)) < 0L)){
    stop("the 2 objects to compare have seqlevels in incompatible orders")
  }
  # This should only insert new seqlevels in the existing ones i.e. it
  # should NEVER drop or reorder existing levels
  seqlevels(x) <- seqlevels(y) <- seqlevels
  
  if (size(x) == 1L) {
    val <- .compareGTuplesCpp(int_seqnames = as.integer(seqnames(x)) - 
                                as.integer(seqnames(y)), 
                              int_strand = as.integer(strand(x)) - 
                                as.integer(strand(y)),
                              int_pos = as.matrix(start(x) - start(y)))
  } else if (size(x) > 1L) {
    val <- .compareGTuplesCpp(int_seqnames = as.integer(seqnames(x)) - 
                                as.integer(seqnames(y)), 
                              int_strand = as.integer(strand(x)) - 
                                as.integer(strand(y)),
                              int_pos = cbind(start(x) - start(y), 
                                              x@internalPos - y@internalPos, 
                                              end(x) - end(y)))
  }
  
  return(val)
}

#' @export
setMethod("compare", 
          c("GTuples", "GTuples"), 
          function(x, y) {
            .GTuples.compare(x, y)
          }
)

#' @export
setMethod("<=", 
          c("GTuples", "GTuples"), 
          function(e1, e2) {
            .GTuples.compare(e1, e2) <= 0L
          }
)

#' @export
setMethod("==", 
          c("GTuples", "GTuples"), 
          function(e1, e2) {
            .GTuples.compare(e1, e2) == 0L
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### duplicated()  
###
### Can't rely on using duplicated() via inheritance from Vector because 
### GTuples inherits from GRanges, rather than from `Vector` directly.
### Furthermore, simply deferring to the duplicated() method for GRanges would 
### only check pos1, posm of of each tuple, which is okay provided size < 3 but 
### not if size >= 3.
###
### unique() will work out-of-the-box on an GTuples object thanks to the 
### method for GRanges objects, which inherits from Vector objects.
.duplicated.GTuples <- function(x, incomparables = FALSE, fromLast = FALSE, 
                                method = c("hash", "base")) {
  
  if (!identical(incomparables, FALSE)) {
    stop("\"duplicated\" method for GTuples objects ",
         "only accepts 'incomparables=FALSE'")
  }
  
  method <- match.arg(method)
  
  # Store tuples as integer matrix
  if (size(x) == 1L) { 
    mat <- c(as.integer(seqnames(x)), as.integer(strand(x)), start(x))
    dim(mat) <- c(length(x), 2 + size(x))
  } else {
    mat <- c(as.integer(seqnames(x)), as.integer(strand(x)), start(x),
             x@internalPos, end(x))
    dim(mat) <- c(length(x), 2 + size(x))
  }
  
  if (method == "base") {
    # base::duplicated.array is slow for large matrices.
    val <- duplicated.array(mat, incomparables = incomparables, MARGIN = 1, 
                            fromLast = fromLast)
  } else {
    # Create a hash of each row of mat by computing the inner product of mat 
    # with a vector of prime numbers.
    # The choice of prime numbers is arbitrary and a better hash function could 
    # be used.
    PRIMES <- c(139L, 7919L, 2129L, 557L, 6857L, 3677L, 761L, 3023L, 6863L, 
                1361L, 6733L, 1811L, 1979L, 5573L, 1129L, 3659L, 3389L, 2383L, 
                211L, 7603L, 7487L, 1171L, 3877L, 7649L, 3767L, 683L, 2819L, 
                2161L, 5237L, 6737L, 1279L, 7547L, 1913L, 4787L, 2593L, 1063L, 
                17L, 4657L, 5659L, 6269L, 2141L, 5693L, 5L, 5399L, 3529L, 
                5189L, 509L, 3329L, 6197L, 4817L, 1741L, 6869L)
    
    # Use the first "2 + size(x)" elements of PRIMES in the inner product
    # If "size(x)" is large then resample from PRIMES; don't expect this to be 
    # required since "length(PRIMES)" is 52, which can handle up to 50-tuples.
    if (size(x) <= length(PRIMES) - 2) {
      primes <- PRIMES[seq_len(2 + size(x))]
    } else {
      primes <- c(PRIMES, sample(PRIMES, 
                                 size(x) - length(PRIMES), replace = TRUE))
    }
    
    hash <- mat %*% primes
    
    # Check the hash for duplicates, which identified **candidate** duplicate 
    # tuples. Need to check from start and end to make sure I get all copies 
    # of any duplicates.
    val <- duplicated(hash) | duplicated(hash, fromLast = TRUE)
    
    # Check whether the candidates are truly duplicates
    val[val] <- duplicated.array(mat[val, ], incomparables = incomparables, 
                                 MARGIN = 1, fromLast = fromLast)
  }
  
  return(as.vector(val))
}

# TODO: Is there is a need for an S3/S4 combo? If so, why?
# S3/S4 combo for duplicated.GTuples
#' @export
duplicated.GTuples <- function(x, incomparables = FALSE, ...) {
  .duplicated.GTuples(x, incomparables = incomparables, ...)
}
#' @export
setMethod("duplicated", 
          "GTuples", 
          .duplicated.GTuples
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### match()
###
### %in%, findMatches(), countMatches() will work out-of-the-box on GTuples 
### objects thanks to the method for Vector objects.

# Effectively just calls findOverlaps with type = equal.
#' @export
setMethod("match", 
          c("GTuples", "GTuples"), 
          function(x, table, nomatch = NA_integer_, incomparables = NULL, 
                   ignore.strand = FALSE) {
            
            if (!isSingleNumberOrNA(nomatch)) {
              stop("'nomatch' must be a single number or NA")
            }
            if (!is.integer(nomatch)) {
              nomatch <- as.integer(nomatch)
            }
            if (!is.null(incomparables)) {
              stop("\"match\" method for GenomicRanges objects only accepts,",
                   "'incomparables=NULL'")
            }
            if (!isTRUEorFALSE(ignore.strand)) {
              stop("'ignore.strand' must be TRUE or FALSE")
            }
            ## Calling merge() is the way to check that 'x' and 'table' are 
            ## based on the same reference genome.
            merge(seqinfo(x), seqinfo(table))
            
            findOverlaps(x, table, type = "equal", select = "first")
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### order() and related methods.
###
### The order() and rank() methods for GTuples objects are consistent with the 
### order implied by compare().
### sort is defined via inheritance to GRanges

#' @export
setMethod(order, 
          "GTuples", 
          function(..., na.last = TRUE, decreasing = FALSE) {
            
            if (!isTRUEorFALSE(decreasing)){
              stop("'decreasing' must be TRUE or FALSE")
            }
            
            args <- list(...)
            
            size <- sapply(args, size)
            
            if (!.zero_range(size)) {
              stop("All GTuples objects must have the same 'size' value.")
            } else{
              size <- size[1]
            }
            
            if (size < 3L){
              # If size < 3 just defer to the order method defined for GRanges
              # Can't simply use callNextMethod() because (oddly) there is no 
              # order method defined for GRanges (rather it is defined for 
              # GenomicRanges).
              order(do.call("c", lapply(args, function(x){as(x, "GRanges")})), 
                    na.last = na.last, decreasing = decreasing)
            } else {
              # If size >= 3 then need to define an order method that takes 
              # note of the "internal positions in each tuple.
              order_args <- vector("list", (size + 2L) * length(args))
              idx <- (size + 2L) * seq_len(length(args))
              order_args[seq.int(from = 1, to = max(idx), by = size + 2)] <- 
                lapply(args, function(x) {
                  as.factor(seqnames(x))
                })
              order_args[seq.int(from = 2, to = max(idx), by = size + 2)] <- 
                lapply(args, function(x) {
                  as.factor(strand(x))
                })
              order_args[seq.int(from = 3, to = max(idx), by = size + 2)] <- 
                lapply(args, start)
              order_args[seq.int(from = 4, to = max(idx), by = size + 2) + 
                           rep(seq(0, size - 3, by = 1))] <- 
                lapply(args, function(x) {
                  x@internalPos
                })
              
              order_args[idx] <- lapply(args, function(x){end(x)})
              do.call(order, c(order_args, list(na.last = na.last, 
                                                decreasing = decreasing)))
            }
          }
)
