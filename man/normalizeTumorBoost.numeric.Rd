%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeTumorBoost.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeTumorBoost.numeric}
\alias{normalizeTumorBoost.numeric}
\alias{normalizeTumorBoost}

\title{Normalizes allele B fractions for a tumor given a match normal}

\description{
 TumorBoost [1] is a normalization method that normalizes the allele B
 fractions of a tumor sample given the allele B fractions and genotypes
 of a matched normal.
 The method is a single-sample (single-pair) method.
 It does not require total copy-number estimates.
 The normalization is done such that the total copy number is
 unchanged afterwards.
}

\usage{
\method{normalizeTumorBoost}{numeric}(betaT, betaN, muN=callNaiveGenotypes(betaN),
  flavor=c("v4", "v3", "v2", "v1"), preserveScale=TRUE, ...)
}

\arguments{
 \item{betaT, betaN}{Two \code{\link[base]{numeric}} \code{\link[base]{vector}}s each of length J with
    tumor and normal allele B fractions, respectively.}
 \item{muN}{An optional \code{\link[base]{vector}} of length J containing
    normal genotypes calls in (0,1/2,1,\code{\link[base]{NA}}) for (AA,AB,BB).}
 \item{flavor}{A \code{\link[base]{character}} string specifying the type of
    correction applied.}
 \item{preserveScale}{If \code{\link[base:logical]{TRUE}}, SNPs that are heterozygous in the
   matched normal are corrected for signal compression using an estimate
   of signal compression based on the amount of correction performed
   by TumorBoost on SNPs that are homozygous in the matched normal.}
 \item{...}{Argument passed to \code{\link{callNaiveGenotypes}}(), if called.}
}

\value{
  Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J containing the normalized
  allele B fractions for the tumor.
  Attribute \code{modelFit} is a \code{\link[base]{list}} containing model fit parameters.
}

\details{
  Allele B fractions are defined as the ratio between the allele B signal
  and the sum of both (all) allele signals at the same locus.
  Allele B fractions are typically within [0,1], but may have a slightly
  wider support due to for instance negative noise.
  This is typically also the case for the returned normalized
  allele B fractions.
}

\section{Flavors}{
 This method provides a few different "flavors" for normalizing the
 data.  The following values of argument \code{flavor} are accepted:
 \itemize{
  \item{v4: (default) The TumorBoost method, i.e. Eqns. (8)-(9) in [1].}
  \item{v3: Eqn (9) in [1] is applied to both heterozygous and homozygous
            SNPs, which effectly is v4 where the normalized allele B
            fractions for homozygous SNPs becomes 0 and 1.}
  \item{v2: ...}
  \item{v1: TumorBoost where correction factor is forced to one, i.e.
            \eqn{\eta_j=1}.  As explained in [1], this is a suboptimal
            normalization method.  See also the discussion in the
            paragraph following Eqn (12) in [1].}
 }
}

\section{Preserving scale}{
 Allele B fractions are more or less compressed toward a half, e.g.
 the signals for homozygous SNPs are slightly away from zero and one.
 The TumorBoost method decreases the correlation in allele B fractions
 between the tumor and the normal \emph{conditioned on the genotype}.
 What it does not control for is the mean level of the allele B fraction
 \emph{conditioned on the genotype}.

 By design, most flavors of the method will correct the homozygous SNPs
 such that their mean levels get close to the expected zero and
 one levels.  However, the heterozygous SNPs will typically keep the
 same mean levels as before.
 One possibility is to adjust the signals such as the mean levels of
 the heterozygous SNPs relative to that of the homozygous SNPs is
 the same after as before the normalization.

 If argument \code{preserveScale=TRUE}, then SNPs that are heterozygous
 (in the matched normal) are corrected for signal compression using
 an estimate of signal compression based on the amount of correction
 performed by TumorBoost on SNPs that are homozygous
 (in the matched normal).

 The option of preserving the scale is \emph{not} discussed in the
 TumorBoost paper [1].
}

\examples{
library(R.utils)

# Load data
pathname <- system.file("data-ex/TumorBoost,fracB,exampleData.Rbin", package="aroma.light")
data <- loadObject(pathname)
attachLocally(data)
pos <- position/1e6
muN <- genotypeN

layout(matrix(1:4, ncol=1))
par(mar=c(2.5,4,0.5,1)+0.1)
ylim <- c(-0.05, 1.05)
col <- rep("#999999", length(muN))
col[muN == 1/2] <- "#000000"

# Allele B fractions for the normal sample
plot(pos, betaN, col=col, ylim=ylim)

# Allele B fractions for the tumor sample
plot(pos, betaT, col=col, ylim=ylim)

# TumorBoost w/ naive genotype calls
betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN)
plot(pos, betaTN, col=col, ylim=ylim)

# TumorBoost w/ external multi-sample genotype calls
betaTNx <- normalizeTumorBoost(betaT=betaT, betaN=betaN, muN=muN)
plot(pos, betaTNx, col=col, ylim=ylim)
}

\author{Henrik Bengtsson, Pierre Neuvial}

\references{
 [1] H. Bengtsson, P. Neuvial & T.P. Speed,
     \emph{TumorBoost: Normalization of allele-specific tumor copy numbers
     from a single pair of tumor-normal genotyping microarrays},
     BMC Bioinformatics, 2010, 11:245. [PMID 20462408]\cr
}
\keyword{methods}
