%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileSpline.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeQuantileSpline}
\alias{normalizeQuantileSpline}

\alias{normalizeQuantileSpline.numeric}
\alias{normalizeQuantileSpline.matrix}
\alias{normalizeQuantileSpline.list}

\title{Normalizes the empirical distribution of one or more samples to a target distribution}

\usage{
\method{normalizeQuantileSpline}{numeric}(x, w=NULL, xTarget, sortTarget=TRUE, robust=TRUE, ...)
\method{normalizeQuantileSpline}{matrix}(X, w=NULL, xTarget=NULL, sortTarget=TRUE, robust=TRUE, ...)
\method{normalizeQuantileSpline}{list}(X, w=NULL, xTarget=NULL, sortTarget=TRUE, robust=TRUE, ...)
}

\description{
  Normalizes the empirical distribution of one or more samples to a target distribution.
  After normalization, all samples have the same average empirical
  density distribution.
}

\arguments{
  \item{x, X}{A single (\eqn{K=1}) \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N},
    a \code{\link[base]{numeric}} \eqn{NxK} \code{\link[base]{matrix}}, or a \code{\link[base]{list}} of length \eqn{K} with
    \code{\link[base]{numeric}} \code{\link[base]{vector}}s, where \eqn{K} represents the number of samples
    and \eqn{N} the number of data points.}
  \item{w}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{N} of weights
    specific to each data point.}
  \item{xTarget}{The target empirical distribution as a \emph{sorted}
    \code{\link[base]{numeric}} \code{\link[base]{vector}} of length \eqn{M}.
    If \code{\link[base]{NULL}} and \code{X} is a \code{\link[base]{list}}, then the target distribution is
    calculated as the average empirical distribution of the samples.}
  \item{sortTarget}{If \code{\link[base:logical]{TRUE}}, argument \code{xTarget} will be sorted,
    otherwise it is assumed to be already sorted.}
  \item{robust}{If \code{\link[base:logical]{TRUE}}, the normalization function is
     estimated robustly.}
  \item{...}{Arguments passed to (\code{\link[stats]{smooth.spline}}
     or \code{\link[aroma.light]{robustSmoothSpline}}).}
}

\value{
  Returns an object of the same type and dimensions as the input.
}

\section{Missing values}{
  Both argument \code{X} and \code{xTarget} may contain non-finite values.
  These values do not affect the estimation of the normalization function.
  Missing values and other non-finite values in \code{X},
  remain in the output as is.  No new missing values are introduced.
}

\examples{
# Simulate three samples with on average 20\% missing values
N <- 10000
X <- cbind(rnorm(N, mean=3, sd=1),
           rnorm(N, mean=4, sd=2),
           rgamma(N, shape=2, rate=1))
X[sample(3*N, size=0.20*3*N)] <- NA

# Plot the data
layout(matrix(c(1,0,2:5), ncol=2, byrow=TRUE))
xlim <- range(X, na.rm=TRUE);
plotDensity(X, lwd=2, xlim=xlim, main="The three original distributions")

Xn <- normalizeQuantile(X)
plotDensity(Xn, lwd=2, xlim=xlim, main="The three normalized distributions")
plotXYCurve(X, Xn, xlim=xlim, main="The three normalized distributions")

Xn2 <- normalizeQuantileSpline(X, xTarget=Xn[,1], spar=0.99)
plotDensity(Xn2, lwd=2, xlim=xlim, main="The three normalized distributions")
plotXYCurve(X, Xn2, xlim=xlim, main="The three normalized distributions")
}

\author{Henrik Bengtsson}

\seealso{
  The target distribution can be calculated as the average
  using \code{\link{averageQuantile}}().

  Internally either
  \code{\link[aroma.light]{robustSmoothSpline}} (\code{robust=TRUE}) or
  \code{\link[stats]{smooth.spline}} (\code{robust=FALSE}) is used.

  An alternative normalization method that is also normalizing the
  empirical densities of samples is \code{\link{normalizeQuantileRank}}().
  Contrary to this method, that method requires that all samples are
  based on the exact same set of data points and it is also more likely
  to over-correct in the tails of the distributions.
}

\references{
  [1] H. Bengtsson, R. Irizarry, B. Carvalho, and T. Speed, \emph{Estimation and assessment of raw copy numbers at the single locus level}, Bioinformatics, 2008.
 \cr
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
