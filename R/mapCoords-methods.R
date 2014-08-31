### =========================================================================
### mapCoords methods
### -------------------------------------------------------------------------
###

# TODO: These might be useful for tuples, e.g. for mapping between different 
# builds of the reference genome. Is this even what this method is designed for?
#' @export
setMethod("mapCoords", 
          c("GTuples", "GTuplesList"), 
          function(x, to, ..., ignore.strand = FALSE, elt.loc = FALSE, 
                   elt.hits = FALSE) {
            stop(paste0(class(x), " do not currently support the 'mapCoords' ", 
                        "method."))
          }
)