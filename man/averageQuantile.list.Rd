%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  averageQuantile.list.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{averageQuantile.list}
\alias{averageQuantile.list}

\title{Gets the average empirical distribution}

\usage{
\method{averageQuantile}{list}(X, ...)
}

\description{
  Gets the average empirical distribution for a set of samples of different sizes.
}

\arguments{
  \item{X}{a \code{\link[base]{list}} with \code{\link[base]{numeric}} \code{\link[base]{vector}}s.  The \code{\link[base]{vector}}s may be of
    different lengths.}
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
  \code{\link[aroma.light:normalizeQuantileRank.list]{*normalizeQuantileRank}()}.
  \code{\link[aroma.light:normalizeQuantileSpline.list]{*normalizeQuantileSpline}()}.
  \code{\link[stats]{quantile}}.
}

\author{
  Parts adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002
  \& 2006.  Original code by Ben Bolstad at Statistics Department,
  University of California.
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
