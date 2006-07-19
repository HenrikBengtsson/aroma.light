%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Thu Jun 29 11:45:04 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2006.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{normalizeQuantile.numeric}
\alias{normalizeQuantile.numeric}
\alias{normalizeQuantile.numeric}


 \title{Normalizes the empirical distribution of a single sample to a target distribution}

 \usage{\method{normalizeQuantile}{numeric}(x, xTarget, ties=FALSE, ...)}

 \description{
   Normalizes the empirical distribution of a single sample to a target distribution.
 }

 \arguments{
   \item{x}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N}.}
   \item{xTarget}{a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{M}.}
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

 \seealso{
   To calculate a target distribution from a set of samples, see
   \code{\link[aroma.light:getAverageQuantile.numeric]{*getAverageQuantile}()}.
   This method is used by \code{\link{normalizeQuantile.list}}().
 }

 \author{
   Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 \& 2006.
   Original code by Ben Bolstad at Statistics Department, University of
   California.
 }




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
