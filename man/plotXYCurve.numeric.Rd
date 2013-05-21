%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plotXYCurve.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotXYCurve.numeric}
\alias{plotXYCurve.numeric}

\title{Plot the relationship between two variables as a smooth curve}

\usage{
\method{plotXYCurve}{numeric}(x, y, lwd=2, col=1, dlwd=1, dcol=NA, xlim=NULL, ylim=xlim, xlab=NULL,
  ylab=NULL, curveFit=smooth.spline, ..., add=FALSE)
}

\description{
  Plot the relationship between two variables as a smooth curve.
}

\arguments{
  \item{x, y}{Two \code{\link[base]{numeric}} \code{\link[base]{vector}} of length N.}
  \item{lwd}{The width of the curve.}
  \item{col}{The color of the curve.}
  \item{dlwd}{The width of the density curves.}
  \item{dcol}{The fill color of the interior of the density curves.}
  \item{xlim, ylim}{The x and y plotting limits.}
  \item{xlab, ylab}{The x and y labels.}
  \item{curveFit}{The \code{\link[base]{function}} used to fit the curve.  The two first
    arguments of the function must take \code{x} and \code{y}, and the
    function must return a \code{\link[base]{list}} with fitted elements \code{x} and
    \code{y}.}
  \item{...}{Additional arguments passed to \code{\link[graphics]{lines}}
    used to draw the curve.}
  \item{add}{If \code{\link[base:logical]{TRUE}}, the graph is added to the current plot, otherwise
    a new plot is created.}
}

\value{
  Returns (invisibly) the curve fit.
}

\section{Missing values}{
  Data points (x,y) with non-finite values are excluded.
}

\seealso{
  \code{\link{plotXYCurve.matrix}}().
}

\author{Henrik Bengtsson}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
