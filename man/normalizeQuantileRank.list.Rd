%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileRank.list.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{normalizeQuantileRank.list}
\alias{normalizeQuantileRank.list}
 \alias{list.normalizeQuantile}
\alias{normalizeQuantile.list}
\alias{normalizeQuantile,list-method}


 \title{Normalizes the empirical distribution of a set of samples to a target distribution}

 \usage{
\method{normalizeQuantileRank}{list}(X, xTarget=NULL, ...)
}

 \description{
   Normalizes the empirical distribution of a set of samples to a target distribution.  The samples may differ in size.
 }

 \arguments{
   \item{X}{a \code{\link[base]{list}} with \code{\link[base]{numeric}} \code{\link[base]{vector}}s.  The \code{\link[base]{vector}}s may be of
     different lengths.}
   \item{xTarget}{The target empirical distribution.  If \code{\link[base]{NULL}}, the target
     distribution is calculated as the average empirical distribution of
     the samples.}
   \item{...}{Passed to \code{\link{normalizeQuantileRank.numeric}}().}
 }

 \value{
   Returns a \code{\link[base]{list}} of normalized \code{\link[base]{numeric}} \code{\link[base]{vector}} of the same lengths as the
   corresponding ones in the input matrix.
 }

 \section{Missing values}{
   Missing values are excluded.
   Values that are \code{\link[base]{NA}} remain \code{\link[base]{NA}} after normalization.
   No new \code{\link[base]{NA}}s are introduced.
 }

 \examples{
# Simulate ten samples of different lengths
N <- 10000
X <- list()
for (kk in 1:8) {
  rfcn <- list(rnorm, rgamma)[[sample(2, size=1)]]
  size <- runif(1, min=0.3, max=1)
  a <- rgamma(1, shape=20, rate=10)
  b <- rgamma(1, shape=10, rate=10)
  values <- rfcn(size*N, a, b)

  # "Censor" values
  values[values < 0 | values > 8] <- NA

  X[[kk]] <- values
}

# Add 20\% missing values
X <- lapply(X, FUN=function(x) {
  x[sample(length(x), size=0.20*length(x))] <- NA;
  x
})

# Normalize quantiles
Xn <- normalizeQuantile(X)

# Plot the data
layout(matrix(1:2, ncol=1))
xlim <- range(X, na.rm=TRUE);
plotDensity(X, lwd=2, xlim=xlim, main="The original distributions")
plotDensity(Xn, lwd=2, xlim=xlim, main="The normalized distributions")
}

 \author{
   Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 \& 2006.
   Original code by Ben Bolstad at Statistics Department, University of
   California.
 }

 \seealso{
   The target empirical distribution is calculated as the average
   using \code{\link[aroma.light:averageQuantile.list]{*averageQuantile}()}.
   Each \code{\link[base]{vector}} is normalized toward this target disribution using
   \code{\link{normalizeQuantileRank.numeric}}().
   \code{\link[aroma.light:normalizeQuantileSpline.list]{*normalizeQuantileSpline}()}.
 }




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
