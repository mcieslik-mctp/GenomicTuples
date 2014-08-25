\name{GTuples-comparison}

\alias{GTuples-comparison}

\alias{compare,GTuples,GTuples-method}

\alias{duplicated,GTuples-method}
\alias{duplicated.GTuples}

\alias{match,GTuples,GTuples-method}

\alias{order,GTuples-method}
\alias{sort,GTuples-method}
\alias{sort.GTuples}
\alias{rank,GTuples-method}


\title{Comparing and ordering genomic tuples}

\description{
  Methods for comparing and ordering the elements in one or more
  \link{GTuples} objects.
}

\usage{
## Element-wise (aka "parallel") comparison of 2 GTuples objects
## -------------------------------------------------------------------

\S4method{==}{GTuples,GTuples}(e1, e2)

\S4method{<=}{GTuples,GTuples}(e1, e2)

## duplicated()
## ------------

\S4method{duplicated}{GTuples}(x, incomparables=FALSE, fromLast=FALSE,
           method=c("hash", "base"))

## match()
## -------

\S4method{match}{GTuples,GTuples}(x, table, nomatch=NA_integer_, incomparables=NULL,
      method=c("auto"), ignore.strand=FALSE)

## order() and related methods
## ----------------------------

\S4method{order}{GTuples}(..., na.last=TRUE, decreasing=FALSE)

\S4method{sort}{GTuples}(x, decreasing=FALSE, ignore.strand=FALSE, by)

\S4method{rank}{GTuples}(x, na.last=TRUE,
     ties.method=c("average", "first", "random", "max", "min"))

## Generalized element-wise (aka "parallel") comparison of 2 GTuples
## objects
## ------------------------------------------------------------------------

\S4method{compare}{GTuples,GTuples}(x, y)
}

\arguments{
  \item{e1, e2, x, table, y}{
    \link{GTuples} objects.
  }
  \item{incomparables}{
    Not supported.
  }
  \item{method}{
    There are two methods implemented: \code{hash} (default) and \code{base}. 
    The \code{base} method is not recommended as it is much slower when the 
    \code{GTuples} object contains a large number of tuples.
  }
  \item{fromLast, nomatch}{
    See \code{?`\link[GRanges]{Ganges-comparison}`} in the IRanges
    package for a description of these arguments.
  }
  \item{ignore.strand}{
    Whether or not the strand should be ignored when comparing 2 genomic
    tuples.
  }
  \item{...}{
    Additional \link{GTuples} objects used for breaking ties.
  }
  \item{na.last}{
    Ignored.
  }
  \item{decreasing}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{ties.method}{
    A character string specifying how ties are treated. Only \code{"first"}
    is supported for now.
  }
  \item{by}{
    An optional formula that is resolved against \code{as.env(x)}; 
    the resulting variables are passed to \code{order} to generate the
    ordering permutation.
  }
}

\details{
  Two elements of a \link{GTuples} object (i.e. two genomic ranges) are
  considered equal iff they are on the same underlying sequence and strand,
  and have the same positions. \code{duplicated()} and \code{unique()}
  on a \link{GTuples} object are conforming to this.

  The "natural order" for the elements of a \link{GTuples} object is to
  order them (a) first by sequence level, (b) then by strand, (c) then by 
  \eqn{pos_{1}, \ldots, pos_{m}}.
  This way, the space of genomic tuples is totally ordered.
  % TODO: Check this applies for GTuples and not just GRanges
  %Note that the \code{reduce} method for \link{GTuples} uses this
  %"natural order" implicitly. 
  \code{order()}, \code{sort()}, and \code{rank()} on a \link{GTuples}
  object are using this "natural order".

  Also \code{==}, \code{!=}, \code{<=}, \code{>=}, \code{<} and \code{>}
  on \link{GTuples} objects are using this "natural order".
}

\author{Peter Hickey}

\seealso{
  \itemize{
    \item The \link{GTuples} class.
    \item \link[GRanges]{GRanges-comparison} in the \strong{GRanges}
          package for comparing and ordering genomic ranges.
    % TODO: Should there be intra-tuple and inter-tuple methods?
    % TODO: Should intra-range and inter-range methods fail for GTuples?
    %\item \link[GTuples]{intra-range-methods} and
    %      \link[GTuples]{inter-range-methods} for
    %      intra and inter range transformations.
    \item \link[GTuples]{setops-methods} for set operations on
          \link{GTuples} objects.
    \item \link[GTuples]{findOverlaps-methods} for finding
          overlapping genomic ranges.
  }
}

\examples{
% TODO: Check all examples work as intended
## Example 1-tuples
gt1a <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               tuples = matrix(c(1:9, 7L), ncol = 1),
               strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
               seqlengths=c(chr1 = 110, chr2 = 120, chr3 = 130))
gt1 <- c(gt1a, gt1a[7:3])

## Example 2-tuples
gt2a <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               tuples = matrix(c(1:9, 7L, 11:19, 17L), ncol = 2),
               strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
               seqlengths=c(chr1 = 110, chr2 = 120, chr3 = 130))
gt2 <- c(gt2a, gt2a[7:3])

## Example 3-tuples
gt3a <- GTuples(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               tuples = matrix(c(1:9, 7L, 11:19, 17L, 21:29, 27L), ncol = 3),
               strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
               seqlengths=c(chr1 = 110, chr2 = 120, chr3 = 130))
gt3 <- c(gt3a, gt3a[7:3])


## ---------------------------------------------------------------------
## A. ELEMENT-WISE (AKA "PARALLEL") COMPARISON OF 2 GTuples OBJECTS
## ---------------------------------------------------------------------
gt1[2] == gt1[2]  # TRUE
gt2[2] == gt2[2]  # TRUE
gt3[2] == gt3[2]  # TRUE

gt1[2] == gt1[5]  # FALSE
gt2[2] == gt2[5]  # FALSE
gt3[2] == gt3[5]  # FALSE

gt1 == gt1[4]
gt2 == gt2[4]
gt3 == gt3[4]

gt1 >= gt1[3]
gt2 >= gt2[3]
gt3 >= gt3[3]

## ---------------------------------------------------------------------
## B. duplicated(), unique()
## ---------------------------------------------------------------------
duplicated(gt1)
duplicated(gt2)
duplicated(gt3)

unique(gt1)
unique(gt2)
unique(gt3)

## ---------------------------------------------------------------------
## C. match(), %in%
## ---------------------------------------------------------------------
table1 <- gt1[1:7]
table2 <- gt2[1:7]
table3 <- gt3[1:7]

match(gt1, table1)
match(gt2, table2)
match(gt3, table3)

match(gt1, table1, ignore.strand = TRUE)
match(gt2, table2, ignore.strand = TRUE)
match(gt3, table3, ignore.strand = TRUE)

gt1 \%in\% table1
gt2 \%in\% table2
gt3 \%in\% table3

## ---------------------------------------------------------------------
## D. findMatches(), countMatches()
## ---------------------------------------------------------------------
findMatches(gt1, table1)
findMatches(gt2, table2)
findMatches(gt3, table3)

countMatches(gt1, table1)
countMatches(gt2, table2)
countMatches(gt3, table3)

findMatches(gt1, table1, ignore.strand = TRUE)
findMatches(gt2, table2, ignore.strand = TRUE)
findMatches(gt3, table3, ignore.strand = TRUE)

countMatches(gt1, table1, ignore.strand = TRUE)
countMatches(gt2, table2, ignore.strand = TRUE)
countMatches(gt3, table3, ignore.strand = TRUE)

gt1_levels <- unique(gt1)
gt2_levels <- unique(gt2)
gt3_levels <- unique(gt3)

countMatches(gt1_levels, gt1)
countMatches(gt2_levels, gt2)
countMatches(gt3_levels, gt3)

## ---------------------------------------------------------------------
## E. order() AND RELATED METHODS
## ---------------------------------------------------------------------
order(gt1)
order(gt2)
order(gt3)

sort(gt1)
sort(gt2)
sort(gt3)

sort(gt1, ignore.strand = TRUE)
sort(gt2, ignore.strand = TRUE)
sort(gt3, ignore.strand = TRUE)

% by is not yet implemented for GTuples
%## TODO: Broken. Please fix!
%#sort(gr, by = ~ seqnames + start + end) # equivalent to (but slower than) above
%score(gr) <- rev(seq_len(length(gr)))
%## TODO: Broken. Please fix!
%#sort(gr, by = ~ score)

rank(gt1)
rank(gt2)
rank(gt3)

## ---------------------------------------------------------------------
## F. GENERALIZED ELEMENT-WISE COMPARISON OF 2 GTuples OBJECTS
## ---------------------------------------------------------------------
% TODO
%gr2 <- GTuples(c(rep("chr1", 12), "chr2"), IRanges(c(1:11, 6:7), width=3))
%strand(gr2)[12] <- "+"
%gr3 <- GTuples("chr1", IRanges(5, 9))

%compare(gr2, gr3)
%rangeComparisonCodeToLetter(compare(gr2, gr3))
%}

\keyword{methods}