### =========================================================================
### GTuples objects: tuples of genomic positions
### -------------------------------------------------------------------------
###

setClassUnion(name = "matrixOrNULL", members = c("matrix", "NULL"))

#' @export
setClass("GTuples",
         contains = "GRanges",
         representation(
           internalPos = "matrixOrNULL", 
           size = "integer"),
         prototype(
           internalPos = NULL,
           size = NA_integer_)
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GTuples.pos <- function(object) {
  
  msg <- NULL
  
  # Check tuples are sorted; only required if m > 1.
  if (isTRUE(object@size > 2) && length(object) != 0L) {
    if (!.allTuplesSorted(pos1 = object@ranges@start, 
                          internal_pos = object@internalPos, 
                          posm = object@ranges@start + 
                            object@ranges@width - 1)) {
      msg <- validMsg(msg, paste0("positions in each tuple must be sorted in ", 
                                  "strictly increasing order, i.e. 'pos1' < ", 
                                  "'pos2'  < ...) < ", 
                                  paste0('pos', object@size)))
    }
  } else if (isTRUE(object@size == 2)) {
    if (isTRUE(any(object@ranges@width == 0))) {
      msg <- validMsg(msg, 
                      paste0("positions in each tuple must be sorted in ", 
                             "strictly increasing order, i.e. 'pos1' < ", 
                             "'pos2'"))
    }
  }
  
  # Check all tuples have positive positions.
  # Only need to check pos1 and posm because already checked 
  # pos1 < internalPos < posm
  # NB: min(x) < 0 is faster than any(x < 0)
  if (!is.na(object@size) && length(object) != 0L) {
    if (min(object@ranges@start) < 0 || min(object@ranges@start + 
                                              object@ranges@width - 1) < 0) {
      msg <- validMsg(msg, paste0("positions in each tuple must be positive ",
                             "integers."))
    }
  }
  
  return(msg)
}

INVALID.GT.COLNAMES <- c("seqnames", "ranges", "strand",
                         "seqlevels", "seqlengths", "isCircular",
                         #"genome",
                         "start", "end", "width", "element",
                         "tuples", "internalPos", "size")

.valid.GTuples.mcols <- function(object) {
  if (any(INVALID.GT.COLNAMES %in% colnames(mcols(object)))) {
    msg <- c("names of metadata columns cannot be one of ",
             paste0("\"", INVALID.GT.COLNAMES, "\"", collapse=", "))
    return(paste(msg, collapse=" "))
  }
  NULL
}

.valid.GTuples <- function(object) {
  
  # Include all .valid.GTuples.* functions in this vector
  msg <- c(.valid.GTuples.pos(object), .valid.GTuples.mcols(object))
  
  if (is.null(msg)){
    return(TRUE)
  } else{
    return(msg)
  }
}

setValidity2("GTuples", .valid.GTuples)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
GTuples <- function(seqnames = Rle(), tuples = matrix(), 
                    strand = Rle("*", length(seqnames)), ..., 
                    seqlengths = NULL, seqinfo = NULL) {
  
  # Only need to check the tuples, all others get checked by the GRanges 
  # constructor
  
  # Check tuples
  if (!is.matrix(tuples)) {
    stop("'tuples' must be an integer matrix") 
  }
  if (!is.integer(tuples)) {
    if (!all(is.na(tuples))) {
      warning("Converting 'tuples' to integer mode")
    }
    mode(tuples) <- "integer"
  }
  
  # Get size of tuples
  if (all(is.na(tuples))) {
    size <- NA_integer_
  } else {
    size <- ncol(tuples)
  }
  
  # Create IRanges
  if (!is.na(size)) {
    ranges <- IRanges(start = tuples[, 1], end = tuples[, size])
  } else {
    ranges <- IRanges()
  }
  
  # Create internalPos
  if (is.na(size) || size < 3) {
    internalPos <- NULL
  } else {
    internalPos <- tuples[, seq(from = 2L, to = size - 1, by = 1L), 
                          drop = FALSE]
  }
  
  # Create GRanges
  gr <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand, 
                seqlengths = seqlengths, seqinfo = seqinfo, ...)
  
  new("GTuples", gr, internalPos = internalPos, size = size)
}

# TODO: Test
#' @export
setMethod("updateObject", 
          "GTuples", 
          function(object, ..., verbose=FALSE) { 
            if (verbose) {
              message("updateObject(object = 'GTuples')")
            }
            if (is(try(object@seqinfo, silent = TRUE), "try-error")) {
              object <- new(class(object),
                            seqnames = object@seqnames,
                            ranges = object@ranges,
                            strand = object@strand,
                            elementMetadata = object@elementMetadata,
                            metadata = object@metadata,
                            seqinfo = Seqinfo(seqnames = 
                                                names(object@seqlengths),
                                              seqlengths = object@seqlengths),
                            size = object@size,
                            internalPos = object@internalPos)
              return(object)
            }
            if (is(try(validObject(object@seqinfo, complete = TRUE), 
                       silent=TRUE), "try-error")) {
              object@seqinfo <- updateObject(object@seqinfo)
              return(object)
            }
            object
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setMethod("as.data.frame", 
          "GTuples", 
          function(x, row.names = NULL, optional = FALSE, ...) {
            tuples <- tuples(x)
            if (missing(row.names)) {
              row.names <- names(x)
            }
            if (!is.null(names(x))) {
              names(x) <- NULL
            }
            mcols_df <- as.data.frame(mcols(x), ...)
            extraColumnNames <- GenomicRanges:::extraColumnSlotNames(x)
            extraColumnNames <- extraColumnNames[extraColumnNames != 
                                                   'internalPos']
            if (length(extraColumnNames) > 0L) {
              mcols_df <- cbind(as.data.frame(extraColumnSlotsAsDF, ...), 
                                mcols_df)
            }
            # TODO: S4 dispatch seems to be going awry so I have to force
            # the conversion of seqnames(x) and strand(x) to factor.
            # This is done using the exact code called by as.factor,Rle-method
#             data.frame(seqnames = as.factor(seqnames(x)), 
#                        as.data.frame(tuples), strand = as.factor(strand(x)), 
#                        mcols_df, row.names = row.names, 
#                        stringsAsFactors = FALSE)
            seqnames <- rep.int(as.factor(runValue(seqnames(x))), 
                                runLength(seqnames(x)))
            strand <- rep.int(as.factor(runValue(strand(x))), 
                              runLength(strand(x)))
            data.frame(seqnames = seqnames,
                       as.data.frame(tuples),
                       strand = strand,
                       mcols_df,
                       row.names = row.names,
                       stringsAsFactors = FALSE)
          }
)

#' @export
setMethod("granges", 
          "GTuples",
          function(x, use.mcols = FALSE) {
            callNextMethod()
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.
###
### From GenomicRanges-class.R "For an object with a pure S4 slot 
### representation, these both map to initialize. Reference classes will want 
### to override 'update'. Other external representations need further 
### customization." Note, however, that these are not exported from 
### GenomicRanges. I think I can safely use these for GTuples via inheritance
### to GenomicRanges, but should be careful whenever using them and test well.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

### Not exported. 'x' *must* be an unnamed list of length >= 1 (not checked).
.unlist_list_of_GTuples <- function(x, ignore.mcols = FALSE) {
  if (!isTRUEorFALSE(ignore.mcols)) {
    stop("'ignore.mcols' must be TRUE or FALSE")
  }
  ans_class <- class(x[[1L]])
  ans_seqinfo <- do.call(merge, lapply(x, seqinfo))
  ans_seqnames <- do.call(c, lapply(x, seqnames))
  ans_ranges <- do.call(c, lapply(x, ranges))
  ans_strand <- do.call(c, lapply(x, strand))
  ans_internalPos <- do.call(rbind, lapply(x, function(xx) {
    xx@internalPos
  }))
  ans_size <- sapply(x, size)[1]
  if (ignore.mcols) {
    ans_mcols <- new("DataFrame", nrows = length(ans_ranges))
  } else {
    ans_mcols <- do.call(rbind, lapply(x, mcols, FALSE))
  }
  new(ans_class, seqnames = ans_seqnames, ranges = ans_ranges, 
      strand = ans_strand, elementMetadata = ans_mcols, seqinfo = ans_seqinfo, 
      size = ans_size, internalPos = ans_internalPos)
}

#' @export
setMethod("c", 
          "GTuples", 
          function(x, ..., ignore.mcols = FALSE, recursive = FALSE) {
            if (!identical(recursive, FALSE)) {
              stop("'recursive' argument not supported")
            }
            if (missing(x)) {
              args <- unname(list(...))
            } else {
              args <- unname(list(x, ...))
            }
            if (!.zero_range(sapply(args, size))) {
              stop("Cannot concatenate GTuples containing tuples of ", 
                   "different 'size'.")
            }
            .unlist_list_of_GTuples(args, ignore.mcols = ignore.mcols)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###
# TODO: Examples in GTuples-class.Rd

#' @include AllGenerics.R
#' @export
setMethod("size", 
          "GTuples", 
          function(x) {
            x@size
          }
)

# TODO: Examples in GTuples-class.Rd
#' @include AllGenerics.R
#' @export
setMethod("tuples", 
          "GTuples", 
          function(x, use.mcols = FALSE) {
            if (use.mcols) {
              stop("Sorry, only use.mcols = FALSE is currently supported.")
            }
            if (!isTRUEorFALSE(use.mcols)) {
              stop("'use.mcols' must be TRUE or FALSE")
            }
            if (is.na(size(x))) {
              ans <- matrix()
            } else if (size(x) == 1L) {
              ans <- as.matrix(start(x))
              colnames(ans) <- paste0('pos', seq_len(size(x)))
            } else{
              ans <- cbind(start(x), x@internalPos, end(x))
              colnames(ans) <- paste0('pos', seq_len(size(x)))
            }
            return(ans)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Splitting
###

# Via inheritance to split,Vector-method.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters
###
#' @export
setReplaceMethod("tuples", 
                 "GTuples", 
                 function(x, value) {
                   if (!is(value, "matrix")) {
                     value <- as(value, "matrix")
                   }
                   mode(value) <- "integer"
                   n <- length(x)
                   k <- length(value)
                   if (k != n) {
                     stop(k, " elements in value to replace ", n, " elements")
                   }
                   m <- ncol(value)
                   if (m != size(x)) {
                     stop(paste0("Size of tuples in 'x' ", size(x), " not ", 
                                 "equal to size of tuples in 'value' ", m))
                   }
                   if (is.na(m)) {
                     x
                   } else if (m == 1L) {
                     start(x) <- value
                     x
                   } else if (m == 2L) {
                     ranges(x) <- value
                     x
                   } else if (m > 2L) {
                     start(x) <- value[, 1]
                     x@internalPos <- value[, seq.int(2, m - 1, 1)]
                     end(x) <- value[, m]
                     x
                   }
                 }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Tuples methods
###

# TODO

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities
###

#' @include AllGenerics.R
#' @export
setMethod("IPD", "GTuples", function(x) {
  size <- size(x)
  if (is.na(size)) {
    stop("Cannot compute IPD from an empty GTuples.")
  } else if (isTRUE(size == 1L)) {
    stop("It does not make sense to compute IPD when size = 1.")
  } else if (isTRUE(size == 2L)) {
    ipd <- width(x)
  } else {
    ipd <- .IPD(start(x), as.matrix(x@internal_pos), end(x))
  }
  
  return(ipd)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###
# TODO: Test the subsetting methods defined in GenomicRanges-class.R

# TODO: Should I explicitly define this via callNextMethod() 
# extractROWS works via inheritance because it handles extraColumnSlots 
# (internalPos)

# TODO: Should I explicitly define this via callNextMethod()
# "[" works via inheritance because it calls extractROWS

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

# TODO: Decide if I should support the with.classinfo argument?
# Ensure the extraPos column "sticks" during subsetting, etc.
setMethod(GenomicRanges:::extraColumnSlotNames, "GTuples",
          function(x) {
            c("internalPos")
          })

# The show method is adapted from that of GRanges
.makeNakedMatFromGTuples <- function(x) {
  lx <- length(x)
  nc <- ncol(mcols(x))
  if (!is.na(x@size)) {
    if (x@size == 1L) {
      ans <- cbind(as.character(x@seqnames), x@ranges@start, 
                   as.character(x@strand))
    } else if (x@size == 2L) {
      ans <- cbind(as.character(x@seqnames), x@ranges@start, 
                   x@ranges@start + x@ranges@width - 1, as.character(x@strand))
    } else {
      ans <- cbind(as.character(x@seqnames), x@ranges@start, 
                   x@internalPos, x@ranges@start + x@ranges@width - 1, 
                   as.character(x@strand))
    }
    colnames(ans) <- c("seqnames", paste0('pos', seq_len(x@size)), "strand")
  } else{
    ans <- cbind(as.character(x@seqnames), as.character(x@strand))
    colnames(ans) <- c("seqnames", "strand")
  }
  extraColumnNames <- GenomicRanges:::extraColumnSlotNames(x)
  extraColumnNames <- extraColumnNames[extraColumnNames != "internalPos"]
  if (length(extraColumnNames) > 0L) {
    ans <- do.call(cbind, c(list(ans), 
                            lapply(GenomicRanges:::extraColumnSlots(x), 
                                   showAsCell)))
  }
  if (nc > 0L) {
    tmp <- do.call(data.frame, c(lapply(mcols(x), S4Vectors:::showAsCell), 
                                 list(check.names = FALSE)))
    ans <- cbind(ans, `|` = rep.int("|", lx), as.matrix(tmp))
  }
  return(ans)
}

showGTuples <- function(x, margin = "", with.classinfo = FALSE, 
                        print.seqlengths = FALSE) {
  if (!identical(with.classinfo, FALSE)) {
    stop("'with.classinfo' not implemented")
  }
  lx <- length(x)
  nc <- ncol(mcols(x))
  
  if (!is.na(x@size)) {
    cat(class(x), " with ", lx, " x ", 
        ifelse(lx == 1L, paste0(x@size, "-tuple"), 
               paste0(x@size, "-tuples")), 
        " and ", nc, " metadata ", ifelse(nc == 1L, "column", "columns"), 
        ":\n", sep = "")
  } else{
    cat(class(x), " with 0 tuples and 0 metadata columns:\n", sep = "")
  }
  
  out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x, .makeNakedMatFromGTuples)
  ## These lines commented out because classinfo is more complicated for GTuples 
  ## objects than GRanges objects. For example, some of the `pos` information 
  ## is stored in an IRanges object while some is stored in a matrix.
  #if (with.classinfo) {
  #    .COL2CLASS <- c(seqnames = "Rle", ranges = "IRanges", 
  #        strand = "Rle")
  #    extraColumnNames <- extraColumnSlotNames(x)
  #    .COL2CLASS <- c(.COL2CLASS, getSlots(class(x))[extraColumnNames])
  #    classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
  #    stopifnot(identical(colnames(classinfo), colnames(out)))
  #    out <- rbind(classinfo, out)
  #}
  
  if (nrow(out) != 0L){ 
    rownames(out) <- paste0(margin, rownames(out))
  }
  print(out, quote = FALSE, right = TRUE)
  if (print.seqlengths) {
    cat(margin, "---\n", sep = "")
    GenomicRanges:::showSeqlengths(x, margin = margin)
  }
}

#' @export
setMethod("show", "GTuples", function(object){
  showGTuples(object, margin="  ", print.seqlengths = TRUE)
})