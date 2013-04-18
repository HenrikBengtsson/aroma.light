%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plotMvsAPairs.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotMvsAPairs.matrix}
\alias{plotMvsAPairs.matrix}

\title{Plot log-ratios/log-intensities for all unique pairs of data vectors}

\description{
  Plot log-ratios/log-intensities for all unique pairs of data vectors.
}

\usage{
\method{plotMvsAPairs}{matrix}(X, Alab="A", Mlab="M", Alim=c(0, 16), Mlim=c(-1, 1) * diff(Alim), pch=".", ..., add=FALSE)
}

\arguments{
 \item{X}{NxK \code{\link[base]{matrix}} where N is the number of observations and
   K is the number of channels.}
 \item{Alab,Mlab}{Labels on the x and y axes.}
 \item{Alim,Mlim}{Plot range on the A and M axes.}
 \item{pch}{Plot symbol used.}
 \item{...}{Additional arguments accepted by \code{\link[graphics]{points}}.}
 \item{add}{If \code{\link[base:logical]{TRUE}}, data points are plotted in the current plot,
   otherwise a new plot is created.}
}

\details{
 Log-ratios and log-intensities are calculated for each neighboring pair
 of channels (columns) and plotted. Thus, in total there will be K-1
 data set plotted.

 The colors used for the plotted pairs are 1, 2, and so on. To change
 the colors, use a different color palette.
}

\value{
  Returns nothing.
}

\author{Henrik Bengtsson}
\keyword{methods}
