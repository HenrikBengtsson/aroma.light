%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileSpline.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantileSpline.numeric}
\alias{normalizeQuantileSpline.numeric}

\title{Normalizes the empirical distribution of one or more samples to a target distribution}

\usage{
\method{normalizeQuantileSpline}{numeric}(x, w=NULL, xTarget, sortTarget=TRUE, ..., robust=TRUE)
}

\description{
  Normalizes the empirical distribution of one or more samples to a target distribution.
  After normalization, all samples have the same average empirical
  density distribution.
}

\arguments{
  \item{x}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N}.}
  \item{w}{an optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N} of weights.}
  \item{xTarget}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N}.}
  \item{sortTarget}{If \code{\link[base:logical]{TRUE}}, argument \code{xTarget} is sorted.}
  \item{...}{Arguments passed to (\code{\link[stats]{smooth.spline}}
     or \code{\link[aroma.light]{robustSmoothSpline}}), e.g. \code{w}
     for weights.}
  \item{robust}{If \code{\link[base:logical]{TRUE}}, the normalization function is
     estimated robustly.}
}

\value{
  Returns an object of the same type and dimensions as the input.
  Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N}.
}

\section{Missing values}{
  Both argument \code{x} and \code{xTarget} may contain non-finite values.
  These values do not affect the estimation of the normalization function.
  Missing values and other non-finite values in \code{x},
  remain in the output as is.  No new missing values are introduced.
}

\author{Henrik Bengtsson}

\seealso{
  %%The target distribution can be calculated as the average
  %%using \code{\link[aroma.light:averageQuantile.numeric]{*averageQuantile}()}.

  Internally either \code{\link[aroma.light]{robustSmoothSpline}} or
  \code{\link[stats]{smooth.spline}} is used.

  An alternative normalization method that is also normalizing the
  empirical densities of samples is \code{\link[aroma.light:normalizeQuantileRank.numeric]{*normalizeQuantileRank}()}.
  Contrary to this method, that method requires that all samples are
  based on the exact same set of data points and it is also more likely
  to over-correct in the tails of the distributions.
}

\references{
  [1] H. Bengtsson, R. Irizarry, B. Carvalho, and T. Speed, \emph{Estimation and assessment of raw copy numbers at the single locus level}, Bioinformatics, 2008.
 \cr
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
