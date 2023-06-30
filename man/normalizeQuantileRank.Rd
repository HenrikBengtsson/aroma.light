%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileRank.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeQuantileRank}
\alias{normalizeQuantileRank}

\alias{normalizeQuantileRank.numeric}
\alias{normalizeQuantileRank.list}
\alias{normalizeQuantile}
\alias{normalizeQuantile.default}

\title{Normalizes the empirical distribution of one of more samples to a target distribution}

\usage{
\method{normalizeQuantileRank}{numeric}(x, xTarget, ties=FALSE, ...)
\method{normalizeQuantileRank}{list}(X, xTarget=NULL, ...)
\method{normalizeQuantile}{default}(x, ...)
}

\description{
  Normalizes the empirical distribution of one of more samples to a target distribution.

  The average sample distribution is calculated either robustly or not
  by utilizing either \code{weightedMedian()} or \code{weighted.mean()}.
  A weighted method is used if any of the weights are different from one.
}

\arguments{
  \item{x, X}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length N or a \code{\link[base]{list}} of length N
    with \code{\link[base]{numeric}} \code{\link[base]{vector}}s.
    If a \code{\link[base]{list}}, then the \code{\link[base]{vector}}s may be of different lengths.}
  \item{xTarget}{The target empirical distribution as a \emph{sorted}
    \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{M}.
    If \code{\link[base]{NULL}} and \code{X} is a \code{\link[base]{list}}, then the target distribution is
    calculated as the average empirical distribution of the samples.}
  \item{ties}{Should ties in \code{x} be treated with care or not?
    For more details, see "limma:normalizeQuantiles".}
  \item{...}{Not used.}
}

\value{
  Returns an object of the same shape as the input argument.
}

\section{Missing values}{
  Missing values are excluded when estimating the "common" (the baseline).
  Values that are \code{\link[base]{NA}} remain \code{\link[base]{NA}} after normalization.
  No new \code{\link[base]{NA}}s are introduced.
}

\section{Weights}{
  Currently only channel weights are support due to the way quantile
  normalization is done.
  If signal weights are given, channel weights are calculated from these
  by taking the mean of the signal weights in each channel.
}

\examples{
# Simulate ten samples of different lengths
N <- 10000
X <- list()
for (kk in 1:8) {
  rfcn <- list(rnorm, rgamma)[[sample(2, size=1)]]
  size <- runif(1, min=0.3, max=1)
  a <- rgamma(1, shape=20, rate=10)
  b <- rgamma(1, shape=10, rate=10)
  values <- rfcn(size*N, a, b)

  # "Censor" values
  values[values < 0 | values > 8] <- NA

  X[[kk]] <- values
}

# Add 20\% missing values
X <- lapply(X, FUN=function(x) {
  x[sample(length(x), size=0.20*length(x))] <- NA
  x
})

# Normalize quantiles
Xn <- normalizeQuantile(X)

# Plot the data
layout(matrix(1:2, ncol=1))
xlim <- range(X, na.rm=TRUE)
plotDensity(X, lwd=2, xlim=xlim, main="The original distributions")
plotDensity(Xn, lwd=2, xlim=xlim, main="The normalized distributions")
}

\author{
  Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 & 2006.
  Original code by Ben Bolstad at Statistics Department, University of
  California.
}

\seealso{
  To calculate a target distribution from a set of samples, see
  \code{\link{averageQuantile}}().
  For an alternative empirical density normalization methods, see
  \code{\link{normalizeQuantileSpline}}().
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
