### =========================================================================
### "coverage" methods
### -------------------------------------------------------------------------
###

setMethod("coverage", 
          "GTuples",
          function(x, shift = 0L, width = NULL, weight = 1L,
                   method = c("auto", "sort", "hash"))
          {
            stop(paste0(class(x), " do not currently support the 'coverage' ", 
                        "method."))
          }
)

setMethod("coverage", 
          "GTuplesList",
          function(x, shift = 0L, width = NULL, weight = 1L,
                   method=c("auto", "sort", "hash"))
          {
            stop(paste0(class(x), " do not currently support the 'coverage' ", 
                        "method."))
          }
)
