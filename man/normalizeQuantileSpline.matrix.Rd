%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileSpline.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantileSpline.matrix}
\alias{normalizeQuantileSpline.matrix}

\title{Weighted sample quantile normalization}

\usage{
\method{normalizeQuantileSpline}{matrix}(X, xTarget, ...)
}

\description{
  Normalizes channels so they all have the same average sample
  distributions.
}

\arguments{
  \item{X}{A \code{\link[base]{numeric}} NxK \code{\link[base]{matrix}} with the K columns representing the
    channels and the N rows representing the data points.}
  \item{xTarget}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length N.}
  \item{...}{Additional arguments passed to
    \code{\link{normalizeQuantileSpline.numeric}}().}
}

\value{
  Returns an NxK \code{\link[base]{matrix}}.
}

\section{Missing values}{
  Both argument \code{X} and \code{xTarget} may contain non-finite values.
  These values do not affect the estimation of the normalization function.
  Non-finite values in \code{X}, remain in the output.
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

\seealso{
  Internally \code{\link{normalizeQuantileSpline.numeric}}() is used.
  \code{\link[aroma.light:normalizeQuantileRank.matrix]{*normalizeQuantileRank}()}.
}

\author{Henrik Bengtsson}

\references{
  [1] H. Bengtsson, R. Irizarry, B. Carvalho, and T. Speed, \emph{Estimation and assessment of raw copy numbers at the single locus level}, Bioinformatics, 2008.
 \cr
}




\keyword{methods}
\keyword{nonparametric}
\keyword{multivariate}
\keyword{robust}
