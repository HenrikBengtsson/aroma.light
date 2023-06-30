###########################################################################/**
# @RdocGeneric normalizeQuantileRank
# @alias normalizeQuantileRank.numeric
# @alias normalizeQuantileRank.list
# @alias normalizeQuantile
# @alias normalizeQuantile.default
#
# @title "Normalizes the empirical distribution of one of more samples to a target distribution"
#
# \usage{
# @usage normalizeQuantileRank,numeric
# @usage normalizeQuantileRank,list
# @usage normalizeQuantile,default
# }
#
# \description{
#   @get "title".
#
#   The average sample distribution is calculated either robustly or not
#   by utilizing either \code{weightedMedian()} or \code{weighted.mean()}.
#   A weighted method is used if any of the weights are different from one.
# }
#
# \arguments{
#   \item{x, X}{a @numeric @vector of length N or a @list of length N
#     with @numeric @vectors.
#     If a @list, then the @vectors may be of different lengths.}
#   \item{xTarget}{The target empirical distribution as a \emph{sorted}
#     @numeric @vector of length \eqn{M}.
#     If @NULL and \code{X} is a @list, then the target distribution is
#     calculated as the average empirical distribution of the samples.}
#   \item{ties}{Should ties in \code{x} be treated with care or not?
#     For more details, see "limma:normalizeQuantiles".}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an object of the same shape as the input argument.
# }
#
# \section{Missing values}{
#   Missing values are excluded when estimating the "common" (the baseline).
#   Values that are @NA remain @NA after normalization.
#   No new @NAs are introduced.
# }
#
# \section{Weights}{
#   Currently only channel weights are support due to the way quantile
#   normalization is done.
#   If signal weights are given, channel weights are calculated from these
#   by taking the mean of the signal weights in each channel.
# }
#
# @examples "../incl/normalizeQuantileRank.list.Rex"
#
# \author{
#   Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 & 2006.
#   Original code by Ben Bolstad at Statistics Department, University of
#   California.
# }
#
# \seealso{
#   To calculate a target distribution from a set of samples, see
#   @see "averageQuantile".
#   For an alternative empirical density normalization methods, see
#   @see "normalizeQuantileSpline".
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


setMethodS3("normalizeQuantileRank", "numeric", function(x, xTarget, ties=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n <- length(x);

  # Argument 'xTarget':
  if (!is.numeric(xTarget)) {
    throw("Argument 'xTarget' is not numeric: ", mode(xTarget));
  }
  nTarget <- length(xTarget);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Different length of sample and target?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nDiff <- (nTarget - n);
  if (nDiff > 0L) {
    # Add hoc fix for differences in lengths.
    naValue <- NA;
    storage.mode(naValue) <- storage.mode(x);
    x <- c(x, rep(naValue, times=nDiff));
    n <- n + nDiff;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For all columns, get for each sample quantile the value of
  # average sample distribution at that quantile.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  quantiles <- (0:(nTarget-1))/(nTarget-1);

  ok <- !is.na(x);
  nok <- sum(ok);

  if(nok < n) {
    # Get the sample quantiles for those values
    if (ties) {
      r <- rank(x[ok]);
      xNew <- (r-1)/(nok-1);
    } else {
      xNew <- (0:(nok-1))/(nok-1);
    }

    # Interpolate to get the xTarget's at positions specified by
    # 'quantile' using data points given by 'xNew' and 'xTarget'.
    if (!ties) {
      # Order and sort the values
      ok <- ((1:n)[ok])[order(x[ok])];
    }

    x[ok] <- approx(x=quantiles, y=xTarget, xout=xNew, ties="ordered")$y;

    if (nDiff > 0L) {
      x <- x[1:(n-nDiff)];
    }
  } else {
    if (ties || n != nTarget) {
      r <- rank(x);
      xNew <- (r-1)/(n-1);
      x <- approx(x=quantiles, y=xTarget, xout=xNew, ties="ordered")$y;
    } else {
      ok <- order(x);
      x[ok] <- xTarget;
    }
  }

  x;
})


setMethodS3("normalizeQuantile", "default", function(x, ...) {
  normalizeQuantileRank(x, ...);
})



##############################################################################
# HISTORY:
# 2013-10-07
# o DOCUMENTATION: Merged the documentation for normalizeQuantileRank()
#   for numeric vectors and lists.
# 2011-04-12
# o Now using NAs of the correct storage type.
# 2008-04-14
# o Renamed normalizeQuantile() to normalizeQuantileRank().  Keeping the old
#   name for backward compatibility.
# 2006-05-21
# o Now 'x' and 'xTarget' may be of different lengths.
# 2006-05-15
# o Now the method can normalize vectors of length different from 'xTarget'.
# 2006-05-12
# o Created from normalizeQuantile.matrix.R.  It has been optimized for
#   memory. Hence, the normalization is done using a two-pass procedure.
##############################################################################
