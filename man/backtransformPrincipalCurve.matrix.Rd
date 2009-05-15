%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  backtransformPrincipalCurve.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{backtransformPrincipalCurve.matrix}
\alias{backtransformPrincipalCurve.matrix}
\alias{backtransformPrincipalCurve.matrix}

\alias{backtransformPrincipalCurve.numeric}

\title{Reverse affine transformation}

\description{
  Reverse affine transformation.
}

\usage{\method{backtransformPrincipalCurve}{matrix}(X, fit, dimensions=NULL, targetDimension=NULL, ...)}

\arguments{
 \item{X}{An NxK \code{\link[base]{matrix}} containing data to be backtransformed.}
 \item{fit}{An object of class \code{principal.curve} as returned by
   \code{\link[aroma.light:fitPrincipalCurve.matrix]{*fitPrincipalCurve}()}.}
 \item{dimensions}{An (optional) subset of of D dimensions all in [1,K]
   to be returned (and backtransform).}
 \item{targetDimension}{An (optional) index specifying the dimension
   in [1,K] to be used as the target dimension.  All other
   dimensions will be normalized toward this dimension, which will
   not be changed.}
 \item{...}{Passed internally to \code{\link[stats]{smooth.spline}}.}
}

\value{
  The backtransformed NxK (or NxD) \code{\link[base]{matrix}}.
}

\examples{\dontrun{See help(fitPrincipalCurve.matrix).}}

\seealso{
  \code{\link[aroma.light:fitPrincipalCurve.matrix]{*fitPrincipalCurve}()}
}
\keyword{methods}
