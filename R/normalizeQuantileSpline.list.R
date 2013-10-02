###########################################################################/**
# @set "class=list"
# @RdocMethod normalizeQuantileSpline
#
# @title "Normalizes the empirical distribution of one or more samples to a target distribution"
#
# @synopsis
#
# \description{
#   @get "title".
#   After normalization, all samples have the same average empirical
#   density distribution.
# }
#
# \arguments{
#   \item{X}{a @list of length K with @numeric @vectors that may be off different lenghts.}
#   \item{xTarget}{The target empirical distribution.
#     If @NULL, the target distribution is calculated as the average
#     empirical distribution of the samples.
#   }
#   \item{...}{Passed to @see "normalizeQuantileSpline.numeric".}
# }
#
# \value{
#   Returns an object of the same type and dimensions as the input.
#   @list of normalized @numeric @vector of the same lengths as the
#   corresponding ones in the input matrix.
# }
#
# \section{Missing values}{
#   Both argument \code{X} and \code{xTarget} may contain non-finite values.
#   These values do not affect the estimation of the normalization function.
#   Missing values and other non-finite values in \code{X},
#   remain in the output as is.  No new missing values are introduced.
# }
#
# @author "HB"
#
# \seealso{
#   The target distribution can be calculated as the average
#   using @seemethod "averageQuantile".
#
#   Internally either @see "aroma.light::robustSmoothSpline" or
#   @see "stats::smooth.spline" is used.
#
#   An alternative normalization method that is also normalizing the
#   empirical densities of samples is @seemethod "normalizeQuantileRank".
#   Contrary to this method, that method requires that all samples are
#   based on the exact same set of data points and it is also more likely
#   to over-correct in the tails of the distributions.
# }
#
# \references{
#   [1] @include "../incl/BengtssonH_etal_2008.bib.Rdoc" \cr
# }
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("normalizeQuantileSpline", "list", function(X, xTarget=NULL, ...) {
  # Get the target quantile for all channels (columns)?
  if (is.null(xTarget))
    xTarget <- averageQuantile(X);

  # Normalizes the data
  nTarget <- length(xTarget);
  X <- lapply(X, FUN=function(x) {
    normalizeQuantileSpline(x, xTarget=xTarget, ...);
  })

  X;
})




##############################################################################
# HISTORY:
# 2008-04-14
# o Created from normalizeQuantileRank.list.R.
##############################################################################
