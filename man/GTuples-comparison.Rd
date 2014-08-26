\name{GTuples-comparison}

\alias{GTuples-comparison}

\alias{compare,GTuples,GTuples-method}

\alias{==,GTuples,GTuples-method}
\alias{<=,GTuples,GTuples-method}

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
## -----------------------------------------------------------------------------

\S4method{==}{GTuples,GTuples}(e1, e2)

\S4method{<=}{GTuples,GTuples}(e1, e2)

## duplicated()
## ------------

\S4method{duplicated}{GTuples}(x, incomparables = FALSE, fromLast = FALSE,
           method = c("hash", "base"))

## match()
## -------

\S4method{match}{GTuples,GTuples}(x, table, nomatch=NA_integer_, 
      incomparables = NULL, method = c("auto"), ignore.strand = FALSE)

## order() and related methods
## -----------------------------------------------------------------------------

\S4method{order}{GTuples}(..., na.last = TRUE, decreasing = FALSE)

\S4method{sort}{GTuples}(x, decreasing = FALSE, ignore.strand = FALSE, by)

\S4method{rank}{GTuples}(x, na.last = TRUE,
     ties.method = c("average", "first", "random", "max", "min"))

## Generalized element-wise (aka "parallel") comparison of 2 GTuples
## objects
## -----------------------------------------------------------------------------

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
  Two elements of a \link{GTuples} object (i.e. two genomic tuples) are
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
    \item \link{intra-tuple-methods} for intra-tuple transformations.
    \item \link[GTuples]{setops-methods} for set operations on
          \link{GTuples} objects.
    \item \link[GTuples]{findOverlaps-methods} for finding
          overlapping genomic ranges.
  }
}

\examples{
## TODO
}
\keyword{methods}