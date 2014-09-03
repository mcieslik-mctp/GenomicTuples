### =========================================================================
### Helper functions not exported.
### For emphasis, __no function in this file will be exported__, because 
### `@keywords internal` applies to the whole file 
### (https://github.com/ramhiser/sparsediscrim/issues/26).
### This means that these functions will not be documented by roxygen2, even
### though the functions have roxygen2 tags.
### =========================================================================

#' Check whether all elements of a numeric vector are identical (within machine precision)
#' @param x a numeric vector.
# 
#' @return TRUE if all elements of the vector are identical (within machine 
#' precision). FALSE in all other cases, including if the vector contains any 
#' NAs.
#' 
#' @export
#' @keywords internal
#' 
#' @note This function is based on Hadley and John's answer to 
#' http://stackoverflow.com/q/4752275
.zero_range <- function(x, tol = .Machine$double.eps ^ 0.5) {
  if (length(x) == 1) {
    val <- TRUE
  } 
  if (any(is.na(x)) & any(is.na(x))){
    val <- FALSE
  } else{
    val <- (abs(max(x) - min(x)) < tol)
  }
  
  return(val)
}

# TODO: This may not be in the latest version of S4Vectors available as a 
# binary for Bioconductor. It is included when S4Vectors is compiled from the 
# source. Figure out why.
# # Define a replaceROWS function with signature NULL. Required for when slots 
# # that are extraColumnSlots are NULL and calling the replaceROWS method for 
# # GTuples (via inheritance to the replaceROWS method for GenomicRanges).
# #' @export
# setMethod("replaceROWS", 
#           "NULL",
#           function(x, i, value) {
#             NULL
#           }
# )