%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  pairedAlleleSpecificCopyNumbers.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{pairedAlleleSpecificCopyNumbers}
\alias{pairedAlleleSpecificCopyNumbers}

\alias{pairedAlleleSpecificCopyNumbers.numeric}

\title{Calculating tumor-normal paired allele-specific copy number stratified on genotypes}

\description{
 Calculating tumor-normal paired allele-specific copy number stratified on genotypes.
 The method is a single-sample (single-pair) method.
 It requires paired tumor-normal parent-specific copy number signals.
}

\usage{
 \method{pairedAlleleSpecificCopyNumbers}{numeric}(thetaT, betaT, thetaN, betaN,
  muN=callNaiveGenotypes(betaN), ...)
}

\arguments{
 \item{thetaT, betaT}{Theta and allele-B fraction signals for the tumor.}
 \item{thetaN, betaN}{Total and allele-B fraction signals for the
    matched normal.}
 \item{muN}{An optional \code{\link[base]{vector}} of length J containing
    normal genotypes calls in (0,1/2,1,\code{\link[base]{NA}}) for (AA,AB,BB).}
 \item{...}{Not used.}
}

\value{
 Returns a \code{\link[base]{data.frame}} with elements \code{CT}, \code{betaT} and \code{muN}.
}

\seealso{
 This definition of calculating tumor-normal paired ASCN is related
 to how the \code{\link{normalizeTumorBoost}}() method calculates normalized
 tumor BAFs.
}

\author{Pierre Neuvial, Henrik Bengtsson}
\keyword{methods}
