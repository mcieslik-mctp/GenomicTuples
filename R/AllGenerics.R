### -------------------------------------------------------------------------
### size
###

#' @export
setGeneric("size", function(x) {
  standardGeneric("size")
})

### -------------------------------------------------------------------------
### IPD
###

#' @export
setGeneric("IPD", function(x) {
  standardGeneric("IPD")
})

### -------------------------------------------------------------------------
### tuples
###

#' @export
setGeneric("tuples", function(x, use.mcols = FALSE) {
  standardGeneric("tuples")
})

#' @export
setGeneric("tuples<-", function(x, ..., value) {
  standardGeneric("tuples<-")
})