%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plotXYCurve.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotXYCurve.matrix}
\alias{plotXYCurve.matrix}

\title{Plot the relationship between two variables as a smooth curve}

\usage{
\method{plotXYCurve}{matrix}(X, Y, col=1:nrow(X), lwd=NULL, xlim=NULL, ylim=xlim, xlab=NULL, ylab=NULL,
  ..., add=FALSE)
}

\description{
  Plot the relationship between two variables as a smooth curve.
}

\arguments{
  \item{X, Y}{Two \code{\link[base]{numeric}} NxK \code{\link[base]{matrix}}.}
  \item{col}{A \code{\link[base]{vector}} of colors to be used for each of columns.}
  \item{lwd}{A \code{\link[base]{vector}} of line widths to be used for each of columns.}
  \item{xlim, ylim}{The x and y plotting limits.}
  \item{xlab, ylab}{The x and y labels.}
  \item{...}{Additional arguments passed to \code{\link{plotXYCurve.numeric}}().}
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
  Internally \code{\link{plotXYCurve.numeric}}() is used.
}

\author{Henrik Bengtsson}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
