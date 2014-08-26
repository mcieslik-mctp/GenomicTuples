### =========================================================================
### GTuplesList objects
### -------------------------------------------------------------------------
###

#' @export
setClass("GTuplesList",
         contains = c("GRangesList"),
         representation(
           unlistData="GTuples",
           elementMetadata="DataFrame"
         ),
         prototype(
           elementType = "GTuples"
         )
)

# TODO: Needed?
setClassUnion("GTuplesORGTuplesList", c("GTuples", "GTuplesList"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

INVALID.GT.COLNAMES <- c("seqnames", "ranges", "strand",
                         "seqlevels", "seqlengths", "isCircular",
                         #"genome",
                         "start", "end", "width", "element",
                         "tuples", "internalPos", "size")
.valid.GTuplesList.mcols <- function(object) {
  msg <- NULL
  object_mcols <- object@elementMetadata
  if (nrow(object_mcols) != length(object)) {
    msg <- "'mcols(object)' has an incorrect number of rows"
  }
  if (any(INVALID.GT.COLNAMES %in% colnames(mcols(object)))) {
    msg <- c("names of metadata columns cannot be one of ",
             paste0("\"", INVALID.GT.COLNAMES, "\"", collapse=", "))
  }
  if (!is.null(rownames(object_mcols))) {
    msg <- c(msg, "'mcols(object)' cannot have row names")
  }
  msg
}

.valid.GTuplesList <- function(x) {
  c(.valid.GTuplesList.mcols(x))
}

setValidity2("GTuplesList", .valid.GTuplesList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @export
GTuplesList <- function(...) {
  listData <- list(...)
  if (length(listData) == 0L) {
    unlistData <- GTuples()
  } else {
    if (length(listData) == 1L && is.list(listData[[1L]])) {
      listData <- listData[[1L]]
    }
    if (!all(sapply(listData, is, "GTuples"))) {
      stop("all elements in '...' must be GTuples objects")
    }
    if (!.zero_range(sapply(listData, size))) {
      stop("all GTuples in '...' must have the same 'size'")
    }
    unlistData <- suppressWarnings(do.call("c", unname(listData)))
  }

  relist(unlistData, PartitioningByEnd(listData))
}

# TODO: Test
#' @export
setMethod("updateObject", 
          "GTuplesList", 
          function(object, ..., verbose = FALSE) { 
            if (verbose) {
              message("updateObject(object = 'GTuplesList')")
            }
            if (is(try(validObject(object@unlistData, complete = TRUE), 
                       silent=TRUE), "try-error")) {
              object@unlistData <- updateObject(object@unlistData)
              return(object)
            }
            object
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors.
###

#' @include AllGenerics.R
#' @export
setMethod("size", 
          "GTuplesList", 
          function(x) {
            size(x[[1]])
          }
)

# TODO: tuples, tuples<-

# TODO: None at this point. Coercion to DataFrameList might be useful.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

# TODO: Does c(gtl, gtl)  correctly deal with size and internalPos slots?
# TODO: grglist when method is implemented in GenomicRanges

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities.
###

#' @include AllGenerics.R
#' @export
setMethod("IPD", 
          "GTuplesList", 
          function(x) {
            unlisted_x <- unlist(x, use.names = FALSE)
            relist(IPD(unlisted_x), x)
          }
)

# TODO: stack (GenomicRangesList-class.R)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

# TODO: Won't be able to rely on subsetting method inherited from GRangesList 
# because this will ignore 'size' and 'internalPos' slots.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going from GRanges to GRangesList with extractList() and family.
###

#' @export
setMethod("relistToClass", 
          "GTuples", 
          function(x) {
            "GTuplesList"
          }
)

# TODO: splitAsListReturnedClass ?

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

# TODO: Decide if I should support the with.classinfo argument?
## The show method is adapted from that of GRangesList
showList <- function(object, showFunction, with.classinfo, ...) {
  k <- length(object)
  cumsumN <- cumsum(elementLengths(object))
  N <- tail(cumsumN, 1)
  cat(class(object), " of length ", k, ":\n", sep = "")
  if (k == 0L) {
    cat("<0 elements>\n\n")
  } else if ((k == 1L) || ((k <= 3L) && (N <= 20L))) {
    nms <- names(object)
    defnms <- paste0("[[", seq_len(k), "]]")
    if (is.null(nms)) {
      nms <- defnms
    } else {
      empty <- nchar(nms) == 0L
      nms[empty] <- defnms[empty]
      nms[!empty] <- paste0("$", nms[!empty])
    }
    for (i in seq_len(k)) {
      cat(nms[i], "\n")
      showFunction(object[[i]], margin="  ",
                   with.classinfo = with.classinfo)
      if (with.classinfo)
        with.classinfo <- FALSE
      cat("\n")
    }
  } else {
    sketch <- function(x) c(head(x, 3), "...", tail(x, 3))
    if (k >= 3 && cumsumN[3L] <= 20) {
      showK <- 3
    } else if (k >= 2 && cumsumN[2L] <= 20) {
      showK <- 2
    } else{
      showK <- 1
    }
    diffK <- k - showK
    nms <- names(object)[seq_len(showK)]
    defnms <- paste0("[[", seq_len(showK), "]]")
    if (is.null(nms)) {
      nms <- defnms
    } else {
      empty <- nchar(nms) == 0L
      nms[empty] <- defnms[empty]
      nms[!empty] <- paste0("$", nms[!empty])
    }
    for (i in seq_len(showK)) {
      cat(nms[i], "\n")
      showFunction(object[[i]], margin="  ",
                   with.classinfo=with.classinfo)
      if (with.classinfo)
        with.classinfo <- FALSE
      cat("\n")
    }
    if (diffK > 0) {
      cat("...\n<", k - showK,
          ifelse(diffK == 1, " more element>\n", " more elements>\n"),
          sep="")
    }
  }
  cat("---\n")
  GenomicRanges:::showSeqlengths(object)
}

setMethod("show", 
          "GTuplesList", 
          function(object) {
            showList(object, showGTuples, FALSE)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deconstruction/reconstruction of a GTuplesList into/from a GTuples
### object.
###
### For internal use only (not exported).
###

# TODO: Required or can inherit from GRangesList?
