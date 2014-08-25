### -------------------------------------------------------------------------
### size
###

#' Get the size of genomic tuples stored in a \code{GTuples} object.
#' 
#' @param x A \code{\link{GTuples}} or \code{\link{GTuplesList}} object.
#' 
#' @export
#' 
#' @return An integer.
setGeneric("size", function(x) {
  standardGeneric("size")
})

### -------------------------------------------------------------------------
### IPD
###

## TODO: Update @return in the case of a GTuplesList
#' Compute the IPD of genomic tuples stored in a \code{\link{GTuples}} or 
#' \code{\link{GTuplesList}} object.
#' 
#' 
#' @param x A \code{\link{GTuples}} object.
#' 
#' @export
#' 
#' @return If \code{size(x)} > 1, an integer matrix with the same number of 
#' rows as \code{length(x)} and number of columns equal to \code{size(x) - 1}. 
#' If \code{size(x)} == 1 or \code{size(x)} is \code{NA} then IPD is not 
#' defined and an error is returned.
setGeneric("IPD", function(x) {
  standardGeneric("IPD")
})

#' @export
setGeneric("tuples", function(x, use.mcols = FALSE) {
  standardGeneric("tuples")
})

#' @export
setGeneric("tuples<-", function(x, ..., value) {
  standardGeneric("tuples<-")
})