\name{intra-tuple-methods}

\alias{intra-tuple-methods}

\alias{shift}
\alias{shift,GTuples-method}
\alias{shift,GTuplesList-method}

\alias{narrow}
\alias{narrow,GTuples-method}

% TODO: Delete these aliases if not required
%\alias{flank}
%\alias{flank,GTuples-method}
%\alias{flank,GTuplesList-method}

%\alias{promoters}
%\alias{promoters,GenomicRanges-method}
%\alias{promoters,GRangesList-method}

%\alias{resize}
%\alias{resize,GenomicRanges-method}

%\alias{restrict}
%\alias{restrict,GenomicRanges-method}
%\alias{restrict,GRangesList-method}

\alias{trim,GenomicRanges-method}

%\alias{Ops,GenomicRanges,numeric-method}

\title{Intra range transformations of a GTuples or GTuplesLists object}

\description{
  This man page documents intra tuple transformations of a
  \code{\link{GTuples}} or a \code{\link{GTuplesList}} object.
  
  \strong{NOTE}: These are \strong{not} the same as intra-range methods defined 
  in the \pkg{GenomicRanges} package 
  (\code{\link[GenomicRanges]{intra-range-methods}}) or in the \pkg{IRanges} 
  package (\code{\link[IRanges]{intra-range-methods}}).
}

% TODO: Update if adding more intra-tuple methods
\usage{
\S4method{shift}{GTuples}(x, shift = 0L, use.names = TRUE)
\S4method{shift}{GTuplesList}(x, shift = 0L, use.names = TRUE)

\S4method{narrow}{GTuples}(x, start = NA, end = NA, width = NA, 
use.names = TRUE)

\S4method{trim}{GTuples}(x, use.names = TRUE)
}

% TODO: Update if adding more intra-tuple methods
\arguments{
  \item{x}{A \code{\link{GTuples}} or \code{\link{GTuplesList}} object.}
  \item{shift, use.names, start, end, width}{
    See \code{?`\link[IRanges]{intra-range-methods}`}.
  }
  \item{ignore.strand}{
    \code{TRUE} or \code{FALSE}. Whether the strand of the input tuples
    should be ignored or not. See details below.
  }
  \item{\dots}{
    Additional arguments to methods.
  }
}

% TODO: Update if adding more intra-tuple methods
\details{
  \itemize{
    \item{}{\code{shift} behaves like the \code{shift} method for 
      \link[GenomicRanges]{GRanges} objects, except that any \code{internalPos} 
      are also shifted. See \code{?`\link[GenomicRanges]{intra-range-methods}`} 
      for further details of the \code{shift} method.
    }
    \item(){\code{narrow} on a \code{\link{GTuples}} object behaves
      like on a \code{\link[GenomicRanges]{GRanges}} object. See 
      \code{?`\link[GenomicRanges]{intra-range-methods}`} for the details.

      A major difference though is that it returns a \code{\link{GTuples}}
      object instead of a \code{\link[GenomicRanges]{GRanges}} object.
      The returned object is \emph{parallel} (i.e. same length and names)
      to the original object \code{x} and any \code{internalPos} are unaltered.
      
      A \code{narrow} operation of a \code{\link{GTuples}} object can only be 
      applied if all tuples remained ordered, i.e. \eqn{pos_{1} < \ldots < 
      pos_{m}}, where \eqn{m} is the \code{\link{size}} of the 
      \code{\link{GTuples}} object.
    }
    \item{}{\code{trim} trims out-of-bound ranges located on non-circular
      sequences whose length is not NA.
  }
}

\author{Peter Hickey}

\seealso{
  \itemize{
    \item \code{\link{GTuples}} and \code{\link{GTuplesList}} objects.
    \item The \link[GenomicRanges]{intra-range-methods} man page in the
          \pkg{GenomicRanges} package.
  }
}

\examples{
## TODO
}

\keyword{methods}
\keyword{utilities}
