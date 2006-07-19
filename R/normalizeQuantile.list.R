###########################################################################/**
# @set "class=list"
# @RdocMethod normalizeQuantile
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
#   \item{...}{Passed to @see "normalizeQuantile.numeric".}
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
# @examples "../incl/normalizeQuantile.list.Rex"
#
# \seealso{
#   The target empirical distribution is calculated as the average 
#   using @seemethod "averageQuantile".
#   Each @vector is normalized toward this target disribution using
#   @seemethod "normalizeQuantile.numeric".
# }
#
# \author{
#   Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 \& 2006.
#   Original code by Ben Bolstad at Statistics Department, University of
#   California.
# }
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("normalizeQuantile", "list", function(X, xTarget=NULL, ...) {
  # Get the target quantile for all channels (columns)?
  if (is.null(xTarget))
    xTarget <- averageQuantile(X);

  # Normalizes the data
  nTarget <- length(xTarget);
  X <- lapply(X, FUN=function(x) {
    normalizeQuantile(x, xTarget=xTarget, ...);
  })

  X;
}) # normalizeQuantile()




##############################################################################
# HISTORY:
# 2006-05-12
# o Created from normalizeQuantile.matrix.R.  It has been optimized for 
#   memory. Hence, the normalization is done using a two-pass procedure.
##############################################################################
