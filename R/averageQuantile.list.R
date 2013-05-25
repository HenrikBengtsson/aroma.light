###########################################################################/**
# @set "class=list"
# @RdocMethod averageQuantile
#
# @title "Gets the average empirical distribution"
#
# @synopsis
#
# \description{
#   @get "title" for a set of samples of different sizes.
# }
#
# \arguments{
#   \item{X}{a @list with @numeric @vectors.  The @vectors may be of
#     different lengths.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric @vector of length equal to the longest @vector
#   in argument \code{X}.
# }
#
# \section{Missing values}{
#   Missing values are excluded.
# }
#
# \seealso{
#   @seemethod "normalizeQuantileRank".
#   @seemethod "normalizeQuantileSpline".
#   @see "stats::quantile".
# }
#
# \author{
#   Parts adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002
#   \& 2006.  Original code by Ben Bolstad at Statistics Department,
#   University of California.
# }
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("averageQuantile", "list", function(X, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfChannels <- length(X);
  if(nbrOfChannels == 1)
    return(X);

  nbrOfObservations <- unlist(lapply(X, FUN=length), use.names=FALSE);
  maxNbrOfObservations <- max(nbrOfObservations);
  if(maxNbrOfObservations == 1)
    return(X);


  # Note: We use variable 'tt' as a temporary variable for anything.
  # Has NAs?
  # tt <- !is.na(Xcc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the sample quantile for all channels (columns)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # A vector specifying the number of observations in each column
  nbrOfFiniteObservations <- rep(maxNbrOfObservations, times=nbrOfChannels);

  # Construct the sample quantiles
  quantiles <- (0:(maxNbrOfObservations-1))/(maxNbrOfObservations-1);

  # Create a vector to hold the target distribution
  xTarget <- vector("double", maxNbrOfObservations);
  for (cc in 1:nbrOfChannels) {
    Xcc <- X[[cc]];

    # Order and sort the values
    Scc <- sort(Xcc);

    # The number of non-NA observations
    nobs <- length(Scc);

    # Too few data points?
    if(nobs < maxNbrOfObservations) {
      # Get the sample quantiles for those values
      bins <- (0:(nobs-1))/(nobs-1);

      # Interpolate to get the values at positions specified by
      # 'quantile' using data points given by 'bins' and 'Scc'.
      Scc <- approx(x=bins, y=Scc, xout=quantiles, ties="ordered")$y;
    }

    # Incremental mean
    xTarget <- xTarget + Scc;
  }
  # Not needed anymore
  Scc <- Xcc <- NULL;

  xTarget <- xTarget/nbrOfChannels;

  xTarget;
}) # averageQuantile()




##############################################################################
# HISTORY:
# 2007-01-22
# o BUG FIX: averageQuantile.list() did not deal with vectors of different
#   length correctly. Thanks Alicia Oshlack, WEHI.
# 2006-05-12
# o Created from normalizeQuantile.matrix.R.  It has been optimized for
#   memory. Hence, the normalization is done using a two-pass procedure.
##############################################################################
