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
### TODO: I need to decide on any inter-tuple methods.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### shift()
###

setMethod("shift", 
          "GTuples",
          function(x, shift = 0L, use.names = TRUE) {
            new_ranges <- shift(ranges(x), shift = shift, use.names = use.names)
            new_internalPos <- x@internalPos + shift
            GenomicRanges:::clone(x, ranges = new_ranges, 
                                  internalPos = new_internalPos) 
          }
)

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

# TODO: flank seems to be broken in BioC devel. 
# TODO: Test and understand what it does, particularly with GTuples and 
# GTuplesList objects.
# TODO: Could define flank for GTuples and GTuplesList via inheritance from
# GRanges and GRangesList but don't know if this really makes sense.

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

# TODO: promoters seems to be broken in BioC devel. 
# TODO: Test and understand what it does, particularly with GTuples and 
# GTuplesList objects.
# TODO: Could define promoters for GTuples and GTuplesList via inheritance from
# GRanges and GRangesList but don't know if this really makes sense.

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

### TODO: Add "reflect" method for GTuples objects once it's added for 
### GenomicRanges

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### resize()
###

# TODO: resize seems to be broken in BioC devel. 
# TODO: Test and understand what it does, particularly with GTuples and 
# GTuplesList objects.
# TODO: Could define resize for GTuples and GTuplesList via inheritance from
# GRanges and GRangesList but don't know if this really makes sense.

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

# TODO: Could define restrict for GTuples and GTuplesList via inheritance from
# GRanges and GRangesList but don't know if this really makes sense.
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

setMethod("Ops", 
          c("GTuples", "numeric"),
          function(e1, e2) {
            stop(paste0(class(x), " do not currently support the 'zoom' ", 
                        "method."))
          }
)