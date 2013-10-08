###########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizeQuantileSpline
#
# @title "Normalizes the empirical distribution of one or more samples to a target distribution"
#
# \usage{
# @usage normalizeQuantileSpline,matrix
# }
#
# \description{
#   @get "title".
#   After normalization, all samples have the same average empirical
#   density distribution.
# }
#
# \arguments{
#   \item{X}{A @numeric NxK @matrix with the K columns representing the
#     channels and the N rows representing the data points.}
#   \item{xTarget}{A @numeric @vector of length N.}
#   \item{...}{Additional arguments passed to
#     @see "normalizeQuantileSpline" for @numeric:s.}
# }
#
# \value{
#   Returns an object of the same type and dimensions as the input.
#   Returns an NxK @matrix.
# }
#
# \section{Missing values}{
#   Both argument \code{X} and \code{xTarget} may contain non-finite values.
#   These values do not affect the estimation of the normalization function.
#   Missing values and other non-finite values in \code{X},
#   remain in the output as is.  No new missing values are introduced.
# }
#
# @examples "../incl/normalizeQuantileSpline.matrix.Rex"
#
# @author "HB"
#
# \seealso{
#   %%The target distribution can be calculated as the average
#   %%using @see "averageQuantile".
#
#   Internally either @see "aroma.light::robustSmoothSpline" or
#   @see "stats::smooth.spline" is used.
#
#   An alternative normalization method that is also normalizing the
#   empirical densities of samples is @see "normalizeQuantileRank".
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
setMethodS3("normalizeQuantileSpline", "matrix", function(X, xTarget, ...) {
  # Argument 'xTarget':
  if (!is.numeric(xTarget)) {
    throw("Argument 'xTarget' is not numeric: ", mode(xTarget));
  }

  if (length(xTarget) != nrow(X)) {
    throw("Argument 'xTarget' is of different length than the number of rows in 'X': ", length(xTarget) , " != ", nrow(X));
  }

  # Sort the target distribution once
  xTarget <- sort(xTarget, na.last=TRUE);

  # Normalize each of the columns towards the target distribution
  for (cc in seq(length=ncol(X))) {
    X[,cc] <- normalizeQuantileSpline(X[,cc], xTarget=xTarget,
                                                      sortTarget=FALSE, ...);
  }

  X;
}) # normalizeQuantileSpline()


##############################################################################
# HISTORY:
# 2013-05-25
# o SPEEDUP: Removed gc() called for every column.
# 2007-02-04
# o Created from normalizeQuantile.matrix.R.
##############################################################################
