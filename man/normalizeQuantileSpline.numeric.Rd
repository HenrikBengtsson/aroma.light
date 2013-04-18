%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileSpline.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantileSpline.numeric}
\alias{normalizeQuantileSpline.numeric}

\title{Normalizes the empirical distribution of a single sample to a target distribution}

\usage{
\method{normalizeQuantileSpline}{numeric}(x, w=NULL, xTarget, sortTarget=TRUE, ..., robust=TRUE)
}

\description{
  Normalizes the empirical distribution of a single sample to a target distribution.
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
  Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N}.
}

\section{Missing values}{
  Both argument \code{X} and \code{xTarget} may contain non-finite values.
  These values do not affect the estimation of the normalization function.
  Non-finite values in \code{X}, remain in the output.
}

\seealso{
  Internally either \code{\link[aroma.light]{robustSmoothSpline}} or
  \code{\link[stats]{smooth.spline}} is used.
  \code{\link{normalizeQuantileSpline.matrix}}().
  \code{\link[aroma.light:normalizeQuantileRank.numeric]{*normalizeQuantileRank}()}.
}

\author{Henrik Bengtsson}

\references{
  [1] H. Bengtsson, R. Irizarry, B. Carvalho, and T. Speed, \emph{Estimation and assessment of raw copy numbers at the single locus level}, Bioinformatics, 2008.
 \cr
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
