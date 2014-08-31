### =========================================================================
### Intra-tuple methods
### -------------------------------------------------------------------------
###
### The methods documented in this page on this page override those in 
### GenomicRanges intra-range-methods.R. Basically, I allow some of these 
### methods (with modifications) and don't allow others, at least for now.
### shift()
### narrow()
### flank()
### promoters()
### reflect()
### resize()
### restrict()
### trim()
### Zooming (symmetrically scales the width).
###
### Some of these methods could be defined for GTuples via inheritance but they 
### would effectively just treat them as GRanges, which may not play nice with
### the internalPos slot. But, and more importantly, it's not clear to me that 
### these methods are sensible or necessary for intra-tuples; I'm happy to 
### implement these if there is a good use case.
### TODO: Are there any other intra-tuple methods that make sense?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### shift()
###

#' @export
setMethod("shift", 
          "GTuples",
          function(x, shift = 0L, use.names = TRUE) {
            new_ranges <- shift(ranges(x), shift = shift, use.names = use.names)
            new_internalPos <- x@internalPos + shift
            GenomicRanges:::clone(x, ranges = new_ranges, 
                                  internalPos = new_internalPos) 
          }
)

#' @export
setMethod("shift", 
          "GTuplesList",
          function(x, shift = 0L, use.names = TRUE) {
            if (is(shift, "IntegerList")) {
              if (length(shift) != length(x) || 
                    any(elementLengths(shift) != elementLengths(x))) {
                stop("IntegerList 'shift' not of same dimension as 'x'")
              }
              shift <- unlist(shift, use.names = FALSE)
            }
            ranges(x@unlistData) <-
              shift(x@unlistData@ranges, shift, use.names = use.names)
              x@unlistData@internalPos <- x@unlistData@internalPos + shift
            x
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### narrow()
###
### Method for GTuples defined via inheritance to GRanges

# TODO: Should I explicitly define this via callNextMethod()?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### flank()
###

#' @export
setMethod("flank", 
          "GTuples", 
          function(x, width, start = TRUE, both = FALSE, use.names = TRUE, 
                   ignore.strand = FALSE) {
            stop(paste0(class(x), " do not currently support the 'flank' ", 
                        "method."))
          }
)

#' @export
setMethod("flank", 
          "GTuplesList", 
          function(x, width, start = TRUE, both = FALSE, use.names = TRUE, 
                   ignore.strand=FALSE) { 
            stop(paste0(class(x), " do not currently support the 'flank' ", 
                        "method."))
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### promoters()
###

#' @export
setMethod("promoters", 
          "GTuples", 
          function(x, upstream = 2000, downstream = 200, ...) {
            stop(paste0(class(x), " do not currently support the 'promoters' ", 
                        "method."))
          }
)

#' @export
setMethod("promoters", 
          "GTuplesList", 
          function(x, upstream = 2000, downstream = 200, ...) { 
            stop(paste0(class(x), " do not currently support the 'promoters' ", 
                        "method."))
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reflect()
###

### TODO: Investigate "reflect" method for GTuples objects once it's added for 
### GenomicRanges

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### resize()
###

#' @export
setMethod("resize", 
          "GTuples", 
          function(x, width, fix = "start", use.names = TRUE, 
                   ignore.strand = FALSE) {
            stop(paste0(class(x), " do not currently support the 'resize' ", 
                        "method."))
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### restrict()
###

#' @export
setMethod("restrict", 
          "GenomicRanges",
          function(x, start = NA, end = NA, keep.all.ranges = FALSE, 
                   use.names = TRUE) {
            stop(paste0(class(x), " do not currently support the 'restrict' ", 
                        "method."))
          }
)

#' @export
setMethod("restrict", 
          "GTuplesList", 
          function(x, start = NA, end = NA, keep.all.ranges = FALSE, 
                   use.names = TRUE) { 
            stop(paste0(class(x), " do not currently support the 'restrict' ", 
                        "method."))
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### trim()
###

### Method for GTuples defined via inheritance to GRanges

# TODO: Should I explicitly define this via callNextMethod()?
# TODO: If trim is required, i.e. the ranges are out of bounds, then the 
# warning messages references GRanges rather than GTuples. Can this be changed 
# to reference the actual class?
# TODO: Currently trim is not defined for GRangesList but if this changes then 
# trim will also be defined via inheritance for GTuplesList.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Zooming (symmetrically scales the width).
###
#' @export
setMethod("Ops", 
          c("GTuples", "numeric"),
          function(e1, e2) {
            stop(paste0(class(e1), " do not currently support the 'zoom' ", 
                        "method."))
          }
)