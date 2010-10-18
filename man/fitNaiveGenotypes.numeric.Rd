%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitNaiveGenotypes.numeric.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{fitNaiveGenotypes.numeric}
\alias{fitNaiveGenotypes.numeric}
\alias{fitNaiveGenotypes.numeric}

 \alias{fitNaiveGenotypes}

 \title{Fit naive genotype model from a normal sample}

 \description{
   Fit naive genotype model from a normal sample.
 }

 \usage{\method{fitNaiveGenotypes}{numeric}(y, cn=rep(2, length(y)), subsetToFit=NULL, flavor=c("density"), adjust=1.5, ..., censorAt=c(-0.1, 1.1), verbose=FALSE)}

 \arguments{
  \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J containing allele B fractions
    for a normal sample.}
  \item{cn}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J specifying the true
    total copy number in \eqn{\{0,1,2,NA\}} at each locus.  This can be
    used to specify which loci are diploid and which are not, e.g.
    autosomal and sex chromosome copy numbers.}
  \item{subsetToFit}{An optional \code{\link[base]{integer}} or \code{\link[base]{logical}} \code{\link[base]{vector}} specifying
    which loci should be used for estimating the model.
    If \code{\link[base]{NULL}}, all loci are used.}
  \item{flavor}{A \code{\link[base]{character}} string specifying the type of algorithm used.}
  \item{adjust}{A postive \code{\link[base]{double}} specifying the amount smoothing for
    the empirical density estimator.}
  \item{...}{Additional arguments passed to \code{\link{findPeaksAndValleys}}().}
  \item{censorAt}{A \code{\link[base]{double}} \code{\link[base]{vector}} of length two specifying the range
    for which values are considered finite.  Values below (above) this
    range are treated as -\code{\link[base:is.finite]{Inf}} (+\code{\link[base:is.finite]{Inf}}).}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
 }

 \value{
   Returns a \code{\link[base]{list}} of \code{\link[base]{list}}s.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   To call genotypes see \code{\link[aroma.light:callNaiveGenotypes.numeric]{*callNaiveGenotypes}()}.
   Internally \code{\link{findPeaksAndValleys}}() is used to identify the thresholds.
 }
\keyword{methods}
