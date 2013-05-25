###########################################################################/**
# @set "class=numeric"
# @RdocMethod normalizeQuantileSpline
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
#   \item{w}{an optional @numeric @vector of length \eqn{N} of weights.}
#   \item{xTarget}{a @numeric @vector of length \eqn{N}.}
#   \item{sortTarget}{If @TRUE, argument \code{xTarget} is sorted.}
#   \item{...}{Arguments passed to (@see "stats::smooth.spline"
#      or @see "aroma.light::robustSmoothSpline"), e.g. \code{w}
#      for weights.}
#   \item{robust}{If @TRUE, the normalization function is
#      estimated robustly.}
# }
#
# \value{
#   Returns a @numeric @vector of length \eqn{N}.
# }
#
# \section{Missing values}{
#   Both argument \code{X} and \code{xTarget} may contain non-finite values.
#   These values do not affect the estimation of the normalization function.
#   Non-finite values in \code{X}, remain in the output.
# }
#
# \seealso{
#   Internally either @see "aroma.light::robustSmoothSpline" or
#   @see "stats::smooth.spline" is used.
#   @see "normalizeQuantileSpline.matrix".
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
setMethodS3("normalizeQuantileSpline", "numeric", function(x, w=NULL, xTarget, sortTarget=TRUE, ..., robust=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n <- length(x);

  # Argument 'w':
  if (!is.null(w)) {
    if (!is.numeric(w)) {
      throw("Argument 'w' is not numeric: ", mode(w));
    }
    if (length(w) != n) {
      throw("Argument 'w' is of different length than 'x': ",
                                                       length(w), " != ", n);
    }
  }

  # Argument 'xTarget':
  if (!is.numeric(xTarget)) {
    throw("Argument 'xTarget' is not numeric: ", mode(xTarget));
  }
  if (length(xTarget) != n) {
    throw("Argument 'xTarget' is of different length than 'x': ",
                                               length(xTarget), " != ", n);
  }


  # Sort target distribution?
  if (sortTarget) {
    xTarget <- sort(xTarget, na.last=TRUE);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # A) Fit normalization function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sort signals (w/ weights) to be normalized
  o <- order(x, na.last=TRUE);
  xx <- x[o];
  if (!is.null(w))
    w <- w[o];

  # Not needed anymore
  o <- NULL;

  # Keep only finite values
  ok <- (is.finite(xx) & is.finite(xTarget));

  # Exclude data points with zero weight
  if (!is.null(w))
    ok <- (ok & w > 0);

  xx <- xx[ok];
  if (!is.null(w))
    w <- w[ok];
  xTarget <- xTarget[ok];

  # Not needed anymore
  ok <- NULL;

  if (robust) {
    # robustSmoothSpline() does not return 'data'.
    fit <- robustSmoothSpline(x=xx, w=w, y=xTarget, ...);
  } else {
    # smooth.spline() returns 'data' by default.
    fit <- smooth.spline(x=xx, w=w, y=xTarget, keep.data=FALSE, ...);
  }

  # Not needed anymore
  xx <- xTarget <- NULL;

  # Not needed below
  fit[c("x", "y", "w", "yin", "call")] <- NULL;  # Saves < 1MB, though.


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # B) Normalize the data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ok <- is.finite(x);
  x[ok] <- predict(fit, x=x[ok])$y;


  x;
}) # normalizeQuantileSpline.numeric()



##############################################################################
# HISTORY:
# 2013-05-25
# o SPEEDUP: Removed all three gc() calls.
# 2007-03-28
# o BUG FIX: Weights 'w' are now correctly ordered.
# o BUG FIX: Due to an incorrect if(), TRUE & FALSE was swapped for 'robust'.
# o Memory optimization; now the fitting is not keeping the data.
# o Renamed argument 'sort' to 'sortTarget'.
# 2007-03-22
# o Updated the Rdocs slightly.
# 2007-02-05
# o Now normalizeQuantileSpline() handles NAs too.
# 2007-02-04
# o Created from normalizeQuantile.numeric.R.
##############################################################################
