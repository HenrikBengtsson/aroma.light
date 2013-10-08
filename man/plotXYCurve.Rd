%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plotXYCurve.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotXYCurve}
\alias{plotXYCurve}

\alias{plotXYCurve.numeric}
\alias{plotXYCurve.matrix}

\title{Plot the relationship between two variables as a smooth curve}

\usage{
\method{plotXYCurve}{numeric}(x, y, col=1L, lwd=2, dlwd=1, dcol=NA, xlim=NULL, ylim=xlim, xlab=NULL,
  ylab=NULL, curveFit=smooth.spline, ..., add=FALSE)
\method{plotXYCurve}{matrix}(X, Y, col=seq(length = nrow(X)), lwd=2, dlwd=1, dcol=NA, xlim=NULL,
  ylim=xlim, xlab=NULL, ylab=NULL, curveFit=smooth.spline, ..., add=FALSE)
}

\description{
  Plot the relationship between two variables as a smooth curve.
}

\arguments{
  \item{x, y, X, Y}{Two \code{\link[base]{numeric}} \code{\link[base]{vector}}s of length N for one curve (K=1),
    or two \code{\link[base]{numeric}} NxK \code{\link[base]{matrix}}:es for K curves.}
  \item{col}{The color of each curve.
    Either a scalar specifying the same value of all curves,
    or a \code{\link[base]{vector}} of K curve-specific values.}
  \item{lwd}{The line width of each curve.
    Either a scalar specifying the same value of all curves,
    or a \code{\link[base]{vector}} of K curve-specific values.}
  \item{dlwd}{The width of each density curve.}
  \item{dcol}{The fill color of the interior of each density curve.}
  \item{xlim, ylim}{The x and y plotting limits.}
  \item{xlab, ylab}{The x and y labels.}
  \item{curveFit}{The \code{\link[base]{function}} used to fit each curve.  The two first
    arguments of the function must take \code{x} and \code{y}, and the
    function must return a \code{\link[base]{list}} with fitted elements \code{x} and
    \code{y}.}
  \item{...}{Additional arguments passed to \code{\link[graphics]{lines}}
    used to draw each curve.}
  \item{add}{If \code{\link[base:logical]{TRUE}}, the graph is added to the current plot, otherwise
    a new plot is created.}
}

\value{
  Returns nothing.
}

\section{Missing values}{
  Data points (x,y) with non-finite values are excluded.
}

\author{Henrik Bengtsson}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
