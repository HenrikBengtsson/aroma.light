%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Thu Jun 29 11:45:06 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2006.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{sampleTuples}
\alias{sampleTuples.default}
\alias{sampleTuples}

\title{Sample tuples of elements from a set}

\description{
  Sample tuples of elements from a set.
  The elements within a sampled tuple are unique, i.e. no two elements
  are the same.
}

\usage{\method{sampleTuples}{default}(x, size, length, ...)}

\arguments{
 \item{x}{A set of elements to sample from.}
 \item{size}{The number of tuples to sample.}
 \item{length}{The length of each tuple.}
 \item{...}{Additional arguments passed to \code{\link[base]{sample}}().}
}

\value{
  Returns a NxK \code{\link[base]{matrix}} where N = \code{size} and K = \code{length}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
pairs <- sampleTuples(1:10, size=5, length=2)
print(pairs)

triples <- sampleTuples(1:10, size=5, length=3)
print(triples)

# Allow tuples with repeated elements
quadruples <- sampleTuples(1:3, size=5, length=4, replace=TRUE)
print(quadruples)
}

\seealso{
  \code{\link[base]{sample}}().
}


\keyword{utilities}
