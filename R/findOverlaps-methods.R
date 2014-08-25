### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------

#' @export
#' 
#' @keywords internal
.findIdentical.GTuples <- function(query, subject, select, ignore.strand) {
  
  # Assumes size is identical for query and subject
  # Create (size - 1) x 2-tuples (pairs) from the tuples as GRanges and 
  # run findOverlaps with parameters chosen to select only equal matches. Then, 
  # intersect the Hits objects from each pair to create the final Hits object.
  # It's clunky, but it works.
    
  # Create all pairs and run findOverlaps
  size <- size(query)
  tuples_idx <- seq_len(size - 1L)
  
  q_seqnames <- seqnames(query)
  s_seqnames <- seqnames(subject)
  if (ignore.strand) {
    q_strand <- Rle('*', length(query))
    s_strand <- Rle('*', length(subject))
  } else{
    q_strand <- strand(query)
    s_strand <- strand(subject)
  }
  q_tuples <- tuples(query)
  s_tuples <- tuples(subject)
  
  pair_hits <- lapply(tuples_idx, function(i, q_seqnames, s_seqnames, q_strand, 
                                   s_strand, q_tuples, s_tuples) {
    query <- GRanges(q_seqnames, IRanges(q_tuples[, i], q_tuples[, i + 1]), 
                     q_strand)
    subject <- GRanges(s_seqnames, IRanges(s_tuples[, i], s_tuples[, i + 1]), 
                       s_strand)
    findOverlaps(query, subject, maxgap = 0L, minoverlap = 1L, type = "equal", 
                 select = select)
    }, q_seqnames , s_seqnames, q_strand, s_strand, q_tuples, s_tuples)
  
  # Intersect the Hits objects
  tuples_hits <- Reduce(intersect, pair_hits)
  
  return(tuples_hits)
}


# TODO: What happens if findOverlaps (type = 'equal') with m-tuples (m > 2) 
# when positions are identical but strands differ, e.g. '+' vs. '+' (should be 
# a match) or '+' vs. '-' (should not be a match).
# TODO: Not quite sure why method breaks if 'ignore.strand = FALSE' rather than 
# '...' is included in setMethod().
# There is a specially defined method for findOverlaps when both the query and 
# the subject are GTuples objects. This is to allow for "equal" matching 
# between GTuples. All other methods are defined via inheritance.
# If either the subject or the query is not a GTuples object then it defers to the findOverlaps method defined for GRanges objects. 
#' @export
setMethod("findOverlaps", signature = c("GTuples", "GTuples"), 
          function(query, subject, maxgap = 0L, minoverlap = 1L, 
                   type = c("equal", "start", "end", "within", "any"), 
                   select = c("all", "first", "last", "arbitrary"), 
                   ...) {
            
            # Argument matching
            select <- match.arg(select)
            type <- match.arg(type)
            
            # Check GTuples are compatible (i.e. have the same size)
            q_size <- size(query)
            s_size <- size(subject)
            if (q_size != s_size) {
              stop("Cannot findOverlaps between '", class(query), "' and '", 
                   class(subject), "' if they have different 'size'.")
            }
            size <- q_size
            
            # TODO: Make this special case faster
            if (size >= 3 && type == 'equal') {
              
              # The 'maxgap' and 'minoverlap' parameters aren't used when 
              # type is 'equal'.
              # Report a warning if these differ from the defaults.
              if (maxgap != 0L || minoverlap != 1L) {
                warning(paste0("'maxgap' and 'minoverlap' ignored when ", 
                               "'type = equal'."))
                maxmap <- 0L
                minoverlap <- 1L
              }
              
              seqinfo <- merge(seqinfo(query), seqinfo(subject))
              # The internal function .findIdentical.GTuples() hasn't been 
              # tested with circular chromosomes.
              
              if (isTRUE(any(isCircular(seqinfo)))) {
                stop("Cannot handle circular chromosomes.")
              }
              
              .findIdentical.GTuples(query, subject, select, ignore.strand)
            } else{
              callNextMethod()
            }
          })

### =========================================================================
### findOverlaps-based methods
### -------------------------------------------------------------------------

setMethod("countOverlaps", signature = c("GTuples", "GTuples"), 
          function(query, subject, 
                   maxgap = 0L, minoverlap = 1L, 
                   type = c("equal", "start", "end", "within", "any"), 
                   ignore.strand = FALSE) {
            counts <- queryHits(findOverlaps(query, subject, maxgap = maxgap, 
                                             minoverlap = minoverlap, 
                                             type = match.arg(type), 
                                             ignore.strand = ignore.strand)) 
            structure(tabulate(counts, NROW(query)), names=names(query))
})

setMethod("overlapsAny", signature = c("GTuples", "GTuples"), 
          function(query, subject, 
                   maxgap = 0L, minoverlap = 1L, 
                   type = c("equal", "start", "end", "within", "any"), 
                   ignore.strand = FALSE) {
            !is.na(findOverlaps(query, subject, maxgap = maxgap,
                                minoverlap = minoverlap,
                                type = match.arg(type),
                                select = "first",
                                ignore.strand = ignore.strand))
          })

setMethod("subsetByOverlaps", signature = c("GTuples", "GTuples"), 
          function(query, subject, 
                   maxgap = 0L, minoverlap = 1L, 
                   type = c("equal", "start", "end", "within", "any"), 
                   ignore.strand = FALSE) {
            query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                      minoverlap = minoverlap,
                                      type = match.arg(type),
                                      select = "first",
                                      ignore.strand = ignore.strand))]
          })