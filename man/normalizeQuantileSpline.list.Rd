%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileSpline.list.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantileSpline.list}
\alias{normalizeQuantileSpline.list}

\title{Normalizes the empirical distribution of a set of samples to a target distribution}

\usage{
\method{normalizeQuantileSpline}{list}(X, xTarget=NULL, ...)
}

\description{
  Normalizes the empirical distribution of a set of samples to a target distribution.  The samples may differ in size.
}

\arguments{
  \item{X}{a \code{\link[base]{list}} with \code{\link[base]{numeric}} \code{\link[base]{vector}}s.  The \code{\link[base]{vector}}s may be of
    different lengths.}
  \item{xTarget}{The target empirical distribution.  If \code{\link[base]{NULL}}, the target
    distribution is calculated as the average empirical distribution of
    the samples.}
  \item{...}{Passed to \code{\link{normalizeQuantileSpline.numeric}}().}
}

\value{
  Returns a \code{\link[base]{list}} of normalized \code{\link[base]{numeric}} \code{\link[base]{vector}} of the same lengths as the
  corresponding ones in the input matrix.
}

\section{Missing values}{
  Missing values are excluded.
  Values that are \code{\link[base]{NA}} remain \code{\link[base]{NA}} after normalization.
  No new \code{\link[base]{NA}}s are introduced.
}

\author{Henrik Bengtsson}

\seealso{
  The target empirical distribution is calculated as the average
  using \code{\link[aroma.light:averageQuantile.list]{*averageQuantile}()}.
  Each \code{\link[base]{vector}} is normalized toward this target disribution using
  \code{\link{normalizeQuantileSpline.numeric}}().
  \code{\link[aroma.light:normalizeQuantileRank.list]{*normalizeQuantileRank}()}.
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
