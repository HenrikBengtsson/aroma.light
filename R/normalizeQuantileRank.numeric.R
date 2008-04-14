###########################################################################/**
# @set "class=numeric"
# @RdocMethod normalizeQuantileRank
# @aliasmethod normalizeQuantile
#
# @title "Normalizes the empirical distribution of a single sample to a target distribution"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{x}{a @numeric @vector of length \eqn{N}.}
#   \item{xTarget}{a \emph{sorted} @numeric @vector of length \eqn{M}.}
#   \item{ties}{Should ties in \code{x} be treated with care or not?  
#     For more details, see "limma:normalizeQuantiles".}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric @vector of length \eqn{N}.
# }
#
# \section{Missing values}{
#   It is only the empirical distribution of the non-missing values that is
#   normalized to the target distribution.  All @NA values remain @NA after
#   normalization.  No new @NAs are introduced.
# }
# 
# \author{
#   Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 \& 2006.
#   Original code by Ben Bolstad at Statistics Department, University of
#   California.
# }
#
# \seealso{
#   To calculate a target distribution from a set of samples, see
#   @seemethod "getAverageQuantile".
#   This method is used by @see "normalizeQuantileRank.list".
#   @seemethod "normalizeQuantileSpline".
# }
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
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
  if (nDiff > 0) {
    # Add hoc fix for differences in lengths.
    x <- c(x, rep(NA, times=nDiff));
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

    if (nDiff > 0) {
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
