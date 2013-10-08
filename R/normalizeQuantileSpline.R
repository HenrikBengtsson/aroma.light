###########################################################################/**
# @RdocGeneric normalizeQuantileSpline
# @alias normalizeQuantileSpline.numeric
# @alias normalizeQuantileSpline.matrix
# @alias normalizeQuantileSpline.list
#
# @title "Normalizes the empirical distribution of one or more samples to a target distribution"
#
# \usage{
# @usage normalizeQuantileSpline,numeric
# @usage normalizeQuantileSpline,matrix
# @usage normalizeQuantileSpline,list
# }
#
# \description{
#   @get "title".
#   After normalization, all samples have the same average empirical
#   density distribution.
# }
#
# \arguments{
#   \item{x, X}{A single (\eqn{K=1}) @numeric @vector of length \eqn{N},
#     a @numeric \eqn{NxK} @matrix, or a @list of length \eqn{K} with
#     @numeric @vectors, where \eqn{K} represents the number of samples
#     and \eqn{N} the number of data points.}
#   \item{w}{An optional @numeric @vector of length \eqn{N} of weights
#     specific to each data point.}
#   \item{xTarget}{The target empirical distribution as a \emph{sorted}
#     @numeric @vector of length \eqn{M}.
#     If @NULL and \code{X} is a @list, then the target distribution is
#     calculated as the average empirical distribution of the samples.}
#   \item{sortTarget}{If @TRUE, argument \code{xTarget} will be sorted,
#     otherwise it is assumed to be already sorted.}
#   \item{robust}{If @TRUE, the normalization function is
#      estimated robustly.}
#   \item{...}{Arguments passed to (@see "stats::smooth.spline"
#      or @see "aroma.light::robustSmoothSpline").}
# }
#
# \value{
#   Returns an object of the same type and dimensions as the input.
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
#   The target distribution can be calculated as the average
#   using @see "averageQuantile".
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
setMethodS3("normalizeQuantileSpline", "list", function(X, w=NULL, xTarget=NULL, sortTarget=TRUE, robust=TRUE, ...) {
  # Argument 'xTarget':
  if (is.null(xTarget)) {
    # Get the target quantile for all channels?
    xTarget <- averageQuantile(X);
    sortTarget <- FALSE;
  } else if (!is.numeric(xTarget)) {
    throw("Argument 'xTarget' is not numeric: ", mode(xTarget));
  }

  # Sort target distribution?
  if (sortTarget) {
    xTarget <- sort(xTarget, na.last=TRUE);
  }

  # Normalizes the data
  nTarget <- length(xTarget);
  X <- lapply(X, FUN=function(x) {
    normalizeQuantileSpline(x, w=w, xTarget=xTarget, sortTarget=FALSE,
                            robust=TRUE, ...);
  })

  X;
})


setMethodS3("normalizeQuantileSpline", "matrix", function(X, w=NULL, xTarget=NULL, sortTarget=TRUE, robust=TRUE, ...) {
  # Argument 'xTarget':
  if (is.null(xTarget)) {
    # Get the target quantile for all channels?
    xTarget <- averageQuantile(X);
    sortTarget <- FALSE;
  } else if (!is.numeric(xTarget)) {
    throw("Argument 'xTarget' is not numeric: ", mode(xTarget));
  }
  if (length(xTarget) != nrow(X)) {
    throw("Argument 'xTarget' is of different length than the number of rows in 'X': ", length(xTarget) , " != ", nrow(X));
  }

  # Sort target distribution?
  if (sortTarget) {
    xTarget <- sort(xTarget, na.last=TRUE);
  }

  # Normalize each of the columns towards the target distribution
  for (cc in seq(length=ncol(X))) {
    X[,cc] <- normalizeQuantileSpline(X[,cc], w=w, xTarget=xTarget,
                              sortTarget=FALSE, robust=robust, ...);
  }

  X;
}) # normalizeQuantileSpline.matrix()


setMethodS3("normalizeQuantileSpline", "numeric", function(x, w=NULL, xTarget, sortTarget=TRUE, robust=TRUE, ...) {
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
# 2013-10-08
# o SPEEDUP: Now normalizeQuantileSpline(..., sortTarget=TRUE) sorts the
#   target only once for lists of vectors just as done for matrices.
# 2013-10-07
# o normalizeQuantileSpline() for list:s gained explicit argument 'w',
#   'sortTarget' and 'robust'.
# 2013-05-25
# o SPEEDUP: Removed all gc() calls.
# 2008-04-14
# o Added normalizeQuantileSpline() for list:s, which works analogously as
#   normalizeQuantileRank() for list:s.
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
