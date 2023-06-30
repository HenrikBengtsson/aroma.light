%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  averageQuantile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{averageQuantile}
\alias{averageQuantile}

\alias{averageQuantile.list}
\alias{averageQuantile.matrix}

\title{Gets the average empirical distribution}

\usage{
\method{averageQuantile}{list}(X, ...)
\method{averageQuantile}{matrix}(X, ...)
}

\description{
  Gets the average empirical distribution for a set of samples.
}

\arguments{
  \item{X}{A \code{\link[base]{list}} with K \code{\link[base]{numeric}} \code{\link[base]{vector}}s, or a \code{\link[base]{numeric}} NxK \code{\link[base]{matrix}}.
    If a \code{\link[base]{list}}, the \code{\link[base]{vector}}s may be of different lengths.}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length equal to the longest \code{\link[base]{vector}}
  in argument \code{X}.
}

\section{Missing values}{
  Missing values are excluded.
}

\seealso{
  \code{\link{normalizeQuantileRank}}().
  \code{\link{normalizeQuantileSpline}}().
  \code{\link[stats]{quantile}}.
}

\author{
  Parts adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002
  & 2006.  Original code by Ben Bolstad at Statistics Department,
  University of California.
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
