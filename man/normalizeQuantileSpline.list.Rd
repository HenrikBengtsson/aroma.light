%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeQuantileSpline.list.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantileSpline.list}
\alias{normalizeQuantileSpline.list}

\title{Normalizes the empirical distribution of one or more samples to a target distribution}

\usage{
\method{normalizeQuantileSpline}{list}(X, xTarget=NULL, ...)
}

\description{
  Normalizes the empirical distribution of one or more samples to a target distribution.
  After normalization, all samples have the same average empirical
  density distribution.
}

\arguments{
  \item{X}{a \code{\link[base]{list}} of length K with \code{\link[base]{numeric}} \code{\link[base]{vector}}s that may be off different lenghts.}
  \item{xTarget}{The target empirical distribution.
    If \code{\link[base]{NULL}}, the target distribution is calculated as the average
    empirical distribution of the samples.
  }
  \item{...}{Passed to \code{\link{normalizeQuantileSpline.numeric}}().}
}

\value{
  Returns an object of the same type and dimensions as the input.
  \code{\link[base]{list}} of normalized \code{\link[base]{numeric}} \code{\link[base]{vector}} of the same lengths as the
  corresponding ones in the input matrix.
}

\section{Missing values}{
  Both argument \code{X} and \code{xTarget} may contain non-finite values.
  These values do not affect the estimation of the normalization function.
  Missing values and other non-finite values in \code{X},
  remain in the output as is.  No new missing values are introduced.
}

\author{Henrik Bengtsson}

\seealso{
  The target distribution can be calculated as the average
  using \code{\link[aroma.light:averageQuantile.list]{*averageQuantile}()}.

  Internally either \code{\link[aroma.light]{robustSmoothSpline}} or
  \code{\link[stats]{smooth.spline}} is used.

  An alternative normalization method that is also normalizing the
  empirical densities of samples is \code{\link[aroma.light:normalizeQuantileRank.list]{*normalizeQuantileRank}()}.
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
