%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileRank.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{normalizeQuantileRank.numeric}
\alias{normalizeQuantileRank.numeric}
 \alias{numeric.normalizeQuantile}
\alias{normalizeQuantile.numeric}
\alias{normalizeQuantile,numeric-method}


 \title{Normalizes the empirical distribution of a single sample to a target distribution}

 \usage{
\method{normalizeQuantileRank}{numeric}(x, xTarget, ties=FALSE, ...)
}

 \description{
   Normalizes the empirical distribution of a single sample to a target distribution.
 }

 \arguments{
   \item{x}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N}.}
   \item{xTarget}{a \emph{sorted} \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{M}.}
   \item{ties}{Should ties in \code{x} be treated with care or not?
     For more details, see "limma:normalizeQuantiles".}
   \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N}.
 }

 \section{Missing values}{
   It is only the empirical distribution of the non-missing values that is
   normalized to the target distribution.  All \code{\link[base]{NA}} values remain \code{\link[base]{NA}} after
   normalization.  No new \code{\link[base]{NA}}s are introduced.
 }

 \author{
   Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 \& 2006.
   Original code by Ben Bolstad at Statistics Department, University of
   California.
 }

 \seealso{
   To calculate a target distribution from a set of samples, see
   \code{\link{averageQuantile.list}}().
   This method is used by \code{\link{normalizeQuantileRank.list}}().
   \code{\link[aroma.light:normalizeQuantileSpline.numeric]{*normalizeQuantileSpline}()}.
 }




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
