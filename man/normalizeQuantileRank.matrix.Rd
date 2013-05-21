%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileRank.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantileRank.matrix}
\alias{normalizeQuantileRank.matrix}
\alias{matrix.normalizeQuantile}
\alias{normalizeQuantile.matrix}
\alias{normalizeQuantile,matrix-method}


\title{Weighted sample quantile normalization}

\usage{
\method{normalizeQuantileRank}{matrix}(X, ties=FALSE, robust=FALSE, weights=NULL,
  typeOfWeights=c("channel", "signal"), ...)
}

\description{
  Normalizes channels so they all have the same average sample
  distributions.

  The average sample distribution is calculated either robustly or not
  by utilizing either \code{weightedMedian()} or \code{weighted.mean()}.
  A weighted method is used if any of the weights are different from one.
}

\arguments{
  \item{X}{a numerical NxK \code{\link[base]{matrix}} with the K columns representing the
    channels and the N rows representing the data points.}
  \item{robust}{If \code{\link[base:logical]{TRUE}}, the (weighted) median function is used for
           calculating the average sample distribution, otherwise the
           (weighted) mean function is used.}
  \item{ties}{Should ties be specially treated or not?}
  \item{weights}{If \code{\link[base]{NULL}}, non-weighted normalization is done.
    If channel weights, this should be a \code{\link[base]{vector}} of length K specifying
    the weights for each channel.
    If signal weights, it should be an NxK \code{\link[base]{matrix}} specifying the
    weights for each signal.
  }
  \item{typeOfWeights}{A \code{\link[base]{character}} string specifying the type of
    weights given in argument \code{weights}.}
  \item{...}{Not used.}
}

\value{
  Returns an NxK \code{\link[base]{matrix}}.
}

\section{Missing values}{
  Missing values are excluded when estimating the "common" (the baseline)
  distribution. Values that are \code{\link[base]{NA}} before remain \code{\link[base]{NA}}. No new \code{\link[base]{NA}}s are
  introduced.
}

\section{Weights}{
  Currently only channel weights are support due to the way quantile
  normalization is done.
  If signal weights are given, channel weights are calculated from these
  by taking the mean of the signal weights in each channel.
}

\examples{
# Simulate three samples with on average 20\% missing values
N <- 10000
X <- cbind(rnorm(N, mean=3, sd=1),
           rnorm(N, mean=4, sd=2),
           rgamma(N, shape=2, rate=1))
X[sample(3*N, size=0.20*3*N)] <- NA

# Normalize quantiles
Xn <- normalizeQuantile(X)

# Plot the data
layout(matrix(1:2, ncol=1))
xlim <- range(X, Xn, na.rm=TRUE);
plotDensity(X, lwd=2, xlim=xlim, main="The three original distributions")
plotDensity(Xn, lwd=2, xlim=xlim, main="The three normalized distributions")
}

\author{
  Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 \& 2006.
  Original code by Ben Bolstad at Statistics Department, University of
  California.
  Support for calculating the average sample distribution using (weighted)
  mean or median was added by Henrik Bengtsson.
}

\seealso{
  \code{\link[stats]{median}}, \code{\link[matrixStats]{weightedMedian}},
  \code{\link[base]{mean}}() and \code{\link[stats]{weighted.mean}}.
  \code{\link[aroma.light:normalizeQuantileSpline.matrix]{*normalizeQuantileSpline}()}.
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
