%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeAverage.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeAverage.matrix}
\alias{normalizeAverage.matrix}
\alias{normalizeAverage.list}

\title{Rescales channel vectors to get the same average}

\description{
 Rescales channel vectors to get the same average.
}

\usage{
\method{normalizeAverage}{matrix}(x, baseline=1, avg=median, targetAvg=2200, ...)
}

\arguments{
  \item{x}{A \code{\link[base]{numeric}} NxK \code{\link[base]{matrix}} (or \code{\link[base]{list}} of length K).}
  \item{baseline}{An \code{\link[base]{integer}} in [1,K] specifying which channel should be
     the baseline.}
  \item{avg}{A \code{\link[base]{function}} for calculating the average of one channel.}
  \item{targetAvg}{The average that each channel should have afterwards.
     If \code{\link[base]{NULL}}, the baseline column sets the target average.}
  \item{...}{Additional arguments passed to the \code{avg} \code{\link[base]{function}}.}
}

\value{
 Returns a normalized \code{\link[base]{numeric}} NxK \code{\link[base]{matrix}} (or \code{\link[base]{list}} of length K).
}

\author{Henrik Bengtsson}
\keyword{methods}
