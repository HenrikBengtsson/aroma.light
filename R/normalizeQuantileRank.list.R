###########################################################################/**
# @set "class=list"
# @RdocMethod normalizeQuantileRank
# @aliasmethod normalizeQuantile
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
#   \item{...}{Passed to @see "normalizeQuantileRank.numeric".}
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
# @examples "../incl/normalizeQuantileRank.list.Rex"
#
# \author{
#   Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 \& 2006.
#   Original code by Ben Bolstad at Statistics Department, University of
#   California.
# }
#
# \seealso{
#   The target empirical distribution is calculated as the average 
#   using @seemethod "averageQuantile".
#   Each @vector is normalized toward this target disribution using
#   @seemethod "normalizeQuantileRank.numeric".
#   @seemethod "normalizeQuantileSpline".
# }
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("normalizeQuantileRank", "list", function(X, xTarget=NULL, ...) {
  # Get the target quantile for all channels (columns)?
  if (is.null(xTarget))
    xTarget <- averageQuantile(X);

  # Normalizes the data
  nTarget <- length(xTarget);
  X <- lapply(X, FUN=function(x) {
    normalizeQuantileRank(x, xTarget=xTarget, ...);
  })

  X;
})




##############################################################################
# HISTORY:
# 2008-04-14
# o Renamed normalizeQuantile() to normalizeQuantileRank().  Keeping the old
#   name for backward compatibility.
# 2006-05-12
# o Created from normalizeQuantile.matrix.R.  It has been optimized for 
#   memory. Hence, the normalization is done using a two-pass procedure.
##############################################################################
