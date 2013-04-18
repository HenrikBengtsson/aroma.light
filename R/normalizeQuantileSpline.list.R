###########################################################################/**
# @set "class=list"
# @RdocMethod normalizeQuantileSpline
#
# @title "Normalizes the empirical distribution of a set of samples to a target distribution"
#
# @synopsis
#
# \description{
#   @get "title".  The samples may differ in size.
# }
#
# \arguments{
#   \item{X}{a @list with @numeric @vectors.  The @vectors may be of
#     different lengths.}
#   \item{xTarget}{The target empirical distribution.  If @NULL, the target
#     distribution is calculated as the average empirical distribution of
#     the samples.}
#   \item{...}{Passed to @see "normalizeQuantileSpline.numeric".}
# }
#
# \value{
#   Returns a @list of normalized @numeric @vector of the same lengths as the
#   corresponding ones in the input matrix.
# }
#
# \section{Missing values}{
#   Missing values are excluded.
#   Values that are @NA remain @NA after normalization.
#   No new @NAs are introduced.
# }
#
# @author "HB"
#
# \seealso{
#   The target empirical distribution is calculated as the average
#   using @seemethod "averageQuantile".
#   Each @vector is normalized toward this target disribution using
#   @see "normalizeQuantileSpline.numeric".
#   @seemethod "normalizeQuantileRank".
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
