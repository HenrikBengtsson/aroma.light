%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  sampleCorrelations.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{sampleCorrelations.matrix}
\alias{sampleCorrelations.matrix}

\title{Calculates the correlation for random pairs of observations}

\description{
  Calculates the correlation for random pairs of observations.
}

\usage{
\method{sampleCorrelations}{matrix}(X, MARGIN=1, pairs=NULL, npairs=max(5000, nrow(X)), ...)
}

\arguments{
 \item{X}{An NxK \code{\link[base]{matrix}} where N >= 2 and K >= 2.}
 \item{MARGIN}{The dimension (1 or 2) in which the observations are.
   If \code{MARGIN==1} (\code{==2}), each row (column) is an observation.}
 \item{pairs}{If a Lx2 \code{\link[base]{matrix}}, the L index pairs for which the
   correlations are calculated.
   If \code{\link[base]{NULL}}, pairs of observations are sampled.}
 \item{npairs}{The number of correlations to calculate.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{double}} \code{\link[base]{vector}} of length \code{npairs}.
}

\author{Henrik Bengtsson}

\examples{

# Simulate 20000 genes with 10 observations each
X <- matrix(rnorm(n=20000), ncol=10)

# Calculate the correlation for 5000 random gene pairs
cor <- sampleCorrelations(X, npairs=5000)
print(summary(cor))

}

\seealso{
  \code{\link[base]{sample}}().
}

\references{
 [1] A. Ploner, L. Miller, P. Hall, J. Bergh & Y. Pawitan.
     \emph{Correlation test to assess low-level processing of high-density
     oligonucleotide microarray data}. BMC Bioinformatics, 2005, vol 6.
}


\keyword{methods}
\keyword{utilities}
