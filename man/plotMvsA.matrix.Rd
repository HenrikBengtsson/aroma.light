%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plotMvsA.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotMvsA.matrix}
\alias{plotMvsA.matrix}

\title{Plot log-ratios vs log-intensities}

\description{
  Plot log-ratios vs log-intensities.
}

\usage{
\method{plotMvsA}{matrix}(X, Alab="A", Mlab="M", Alim=c(0, 16), Mlim=c(-1, 1) * diff(Alim) * aspectRatio, aspectRatio=1, pch=".", ..., add=FALSE)
}

\arguments{
 \item{X}{Nx2 \code{\link[base]{matrix}} with two channels and N observations.}
 \item{Alab,Mlab}{Labels on the x and y axes.}
 \item{Alim,Mlim}{Plot range on the A and M axes.}
 \item{aspectRatio}{Aspect ratio between \code{Mlim} and \code{Alim}.}
 \item{pch}{Plot symbol used.}
 \item{...}{Additional arguments accepted by \code{\link[graphics]{points}}.}
 \item{add}{If \code{\link[base:logical]{TRUE}}, data points are plotted in the current plot,
   otherwise a new plot is created.}
}

\details{
 Red channel is assumed to be in column one and green in column two.
 Log-ratio are calculated taking channel one over channel two.
}

\value{
  Returns nothing.
}

\author{Henrik Bengtsson}
\keyword{methods}
