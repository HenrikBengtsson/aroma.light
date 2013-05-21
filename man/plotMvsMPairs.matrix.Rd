%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  plotMvsMPairs.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{plotMvsMPairs.matrix}
\alias{plotMvsMPairs.matrix}

\title{Plot log-ratios vs log-ratios for all pairs of columns}

\description{
  Plot log-ratios vs log-ratios for all pairs of columns.
}

\usage{
\method{plotMvsMPairs}{matrix}(X, xlab="M", ylab="M", xlim=c(-1, 1) * 6, ylim=xlim, pch=".", ...,
  add=FALSE)
}

\arguments{
 \item{X}{Nx2K \code{\link[base]{matrix}} where N is the number of observations and
   2K is an even number of channels.}
 \item{xlab,ylab}{Labels on the x and y axes.}
 \item{xlim,ylim}{Plot range on the x and y axes.}
 \item{pch}{Plot symbol used.}
 \item{...}{Additional arguments accepted by \code{\link[graphics]{points}}.}
 \item{add}{If \code{\link[base:logical]{TRUE}}, data points are plotted in the current plot,
   otherwise a new plot is created.}
}

\details{
 Log-ratio are calculated by over paired columns, e.g. column 1 and 2,
 column 3 and 4, and so on.
}

\value{
  Returns nothing.
}

\author{Henrik Bengtsson}
\keyword{methods}
