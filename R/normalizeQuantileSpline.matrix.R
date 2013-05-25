###########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizeQuantileSpline
#
# @title "Weighted sample quantile normalization"
#
# @synopsis
#
# \description{
#   Normalizes channels so they all have the same average sample
#   distributions.
# }
#
# \arguments{
#   \item{X}{A @numeric NxK @matrix with the K columns representing the
#     channels and the N rows representing the data points.}
#   \item{xTarget}{A @numeric @vector of length N.}
#   \item{...}{Additional arguments passed to
#     @see "normalizeQuantileSpline.numeric".}
# }
#
# \value{
#   Returns an NxK @matrix.
# }
#
# \section{Missing values}{
#   Both argument \code{X} and \code{xTarget} may contain non-finite values.
#   These values do not affect the estimation of the normalization function.
#   Non-finite values in \code{X}, remain in the output.
# }
#
# @examples "../incl/normalizeQuantileSpline.matrix.Rex"
#
# \seealso{
#   Internally @see "normalizeQuantileSpline.numeric" is used.
#   @seemethod "normalizeQuantileRank".
# }
#
# @author
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
