###########################################################################/**
# @RdocGeneric averageQuantile
# @alias averageQuantile.list
# @alias averageQuantile.matrix
#
# @title "Gets the average empirical distribution"
#
# \usage{
# @usage averageQuantile,list
# @usage averageQuantile,matrix
# }
#
# \description{
#   @get "title" for a set of samples.
# }
#
# \arguments{
#   \item{X}{A @list with K @numeric @vectors, or a @numeric NxK @matrix.
#     If a @list, the @vectors may be of different lengths.}
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
#   @see "normalizeQuantileRank".
#   @see "normalizeQuantileSpline".
#   @see "stats::quantile".
# }
#
# \author{
#   Parts adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002
#   & 2006.  Original code by Ben Bolstad at Statistics Department,
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
  # Argument 'X':
  nbrOfChannels <- length(X);
  # Nothing to do?
  if(nbrOfChannels == 1L)
    return(X);

  nbrOfObservations <- unlist(lapply(X, FUN=length), use.names=FALSE);
  maxNbrOfObservations <- max(nbrOfObservations);
  # Nothing to do?
  if(maxNbrOfObservations == 1L)
    return(X);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the sample quantile for all channels (columns)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Construct the sample quantiles
  quantiles <- (0:(maxNbrOfObservations-1L))/(maxNbrOfObservations-1L);

  # Create a vector to hold the target distribution
  xTarget <- vector("double", length=maxNbrOfObservations);
  for (cc in 1:nbrOfChannels) {
    Xcc <- X[[cc]];

    # Order and sort the values
    Scc <- sort(Xcc, na.last=NA);

    # The number of non-NA observations
    nobs <- length(Scc);

    # Too few data points?
    if(nobs < maxNbrOfObservations) {
      # Get the sample quantiles for those values
      bins <- (0:(nobs-1L))/(nobs-1L);

      # Interpolate to get the values at positions specified by
      # 'quantile' using data points given by 'bins' and 'Scc'.
      Scc <- approx(x=bins, y=Scc, xout=quantiles, ties="ordered")$y;

      bins <- NULL; # Not needed anymore
    }

    # Incremental mean
    xTarget <- xTarget + Scc;

    Scc <- NULL; # Not needed anymore
  }
  Xcc <- NULL; # Not needed anymore

  xTarget <- xTarget/nbrOfChannels;

  xTarget;
}) # averageQuantile.list()


setMethodS3("averageQuantile", "matrix", function(X, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfChannels <- ncol(X);
  # Nothing to do?
  if(nbrOfChannels == 1L)
    return(X);

  nbrOfObservations <- nrow(X);
  # Nothing to do?
  if(nbrOfObservations == 1L)
    return(X);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the sample quantile for all channels (columns)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Construct the sample quantiles
  quantiles <- (0:(nbrOfObservations-1L))/(nbrOfObservations-1L);

  # Create a vector to hold the target distribution
  xTarget <- vector("double", length=nbrOfObservations);
  for (cc in 1:nbrOfChannels) {
    Xcc <- X[,cc,drop=TRUE];

    # Order and sort the values
    Scc <- sort(Xcc, na.last=NA);

    # The number of non-NA observations
    nobs <- length(Scc);

    # Too few data points?
    if(nobs < nbrOfObservations) {
      # Get the sample quantiles for those values
      bins <- (0:(nobs-1L))/(nobs-1L);

      # Interpolate to get the values at positions specified by
      # 'quantile' using data points given by 'bins' and 'Scc'.
      Scc <- approx(x=bins, y=Scc, xout=quantiles, ties="ordered")$y;

      bins <- NULL; # Not needed anymore
    }

    # Incremental mean
    xTarget <- xTarget + Scc;

    Scc <- NULL; # Not needed anymore
  }
  Xcc <- NULL; # Not needed anymore

  xTarget <- xTarget/nbrOfChannels;

  xTarget;
}) # averageQuantile.matrix()




##############################################################################
# HISTORY:
# 2013-10-08
# o Added averageQuantile() for matrices.
# 2007-01-22
# o BUG FIX: averageQuantile.list() did not deal with vectors of different
#   length correctly. Thanks Alicia Oshlack, WEHI.
# 2006-05-12
# o Created from normalizeQuantile.matrix.R.  It has been optimized for
#   memory. Hence, the normalization is done using a two-pass procedure.
##############################################################################
