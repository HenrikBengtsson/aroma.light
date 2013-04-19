#########################################################################/**
# @set "class=matrix"
# @RdocMethod backtransformPrincipalCurve
# @alias backtransformPrincipalCurve.numeric
#
# @title "Reverse transformation of principal-curve fit"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An NxK @matrix containing data to be backtransformed.}
#  \item{fit}{An MxL principal-curve fit object of class
#    \code{principal.curve} as returned by @seemethod "fitPrincipalCurve".
#    Typically \eqn{L = K}, but not always.
#  }
#  \item{dimensions}{An (optional) subset of of D dimensions all in [1,L]
#    to be returned (and backtransform).}
#  \item{targetDimension}{An (optional) index specifying the dimension
#    in [1,L] to be used as the target dimension of the \code{fit}.
#    More details below.}
#  \item{...}{Passed internally to @see "stats::smooth.spline".}
# }
#
# \value{
#   The backtransformed NxK (or NxD) @matrix.
# }
#
# \details{
#   Each column in X ("dimension") is backtransformed independentently
#   of the others.
# }
#
# \section{Target dimension}{
#   By default, the backtransform is such that afterward the signals are
#   approximately proportional to the (first) principal curve as fitted
#   by @seemethod "fitPrincipalCurve".  This scale and origin of this
#   principal curve is not uniquely defined.
#   If \code{targetDimension} is specified, then the backtransformed signals
#   are approximately proportional to the signals of the target dimension,
#   and the signals in the target dimension are unchanged.
# }
#
# \section{Subsetting dimensions}{
#   Argument \code{dimensions} can be used to backtransform a subset of
#   dimensions (K) based on a subset of the fitted dimensions (L).
#   If \eqn{K = L}, then both \code{X} and \code{fit} is subsetted.
#   If \eqn{K <> L}, then it is assumed that \code{X} is already
#   subsetted/expanded and only \code{fit} is subsetted.
# }
#
# @examples "../incl/backtransformPrincipalCurve.matrix.Rex"
#
# \seealso{
#   @seemethod "fitPrincipalCurve"
# }
#*/#########################################################################
setMethodS3("backtransformPrincipalCurve", "matrix", function(X, fit, dimensions=NULL, targetDimension=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'X'
  if (!is.numeric(X)) {
    stop("Argument 'X' is not numeric: ", mode(X));
  }

  dimnamesX <- dimnames(X);
  dimX <- dim(X);
  K <- dimX[2];
  if (!is.matrix(X)) {
    X <- as.matrix(X);
  }

  # Argument 'fit'
  if (!inherits(fit, "principal.curve")) {
    stop("Argument 'fit' is not a principal.curve object: ", class(fit)[1]);
  }

  # Argument 'dimensions'
  dimS <- dim(fit$s);
  L <- dimS[2];
  if (!is.null(dimensions)) {
    dimensions <- as.integer(dimensions);
    if (any(dimensions < 1 | dimensions > L)) {
      stop("Argument 'dimensions' contains values out of range [1,", L, "]");
    }
  }

  # Argument 'targetDimension':
  if (!is.null(targetDimension)) {
    targetDimension <- as.integer(targetDimension);
    if (length(targetDimension) != 1) {
      stop("Argument 'targetDimension' should be a scalar or NULL.");
    }
    if (targetDimension < 1 | targetDimension > L) {
      stop("Argument 'targetDimension' is out of range [1,", L, "]: ",
                                                           targetDimension);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transform towards a target dimension?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hasTargetDimension <- (!is.null(targetDimension));
  if (hasTargetDimension) {
    lambda <- fit$s[,targetDimension];
  } else {
    lambda <- fit$lambda;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset dimensions?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  s <- fit$s;
  if (!is.null(dimensions)) {
    s <- s[,dimensions,drop=FALSE];
    if (K == L) {
      X <- X[,dimensions,drop=FALSE];
      dimX <- dim(X);
      dimnamesX <-   dimnames(X);
    }
    dimS <- dim(s);
    L <- dimS[2];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find backtransformations and backtransform data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  naValue <- NA;
  mode(naValue) <- mode(X);
  Xhat <- matrix(naValue, nrow=dimX[1], ncol=dimX[2]);

  okLambda <- is.finite(lambda);

  for (kk in seq(length=L)) {
    sKK <- s[,kk];
    ok <- (is.finite(sKK) & okLambda);
    fitKK <- smooth.spline(sKK[ok], lambda[ok], ...);

    Xkk <- X[,kk];
    keep <- which(is.finite(Xkk));
    Xkk <- Xkk[keep];
    XhatKK <- predict(fitKK, x=Xkk)$y;

    # Sanity check
    stopifnot(length(XhatKK) == length(keep));

    Xhat[keep,kk] <- XhatKK;
  }

  rm(sKK, lambda, fitKK, XhatKK, keep, s);

  dim(Xhat) <- dimX;
  dimnames(Xhat) <- dimnamesX;

  Xhat;
}) # backtransformPrincipalCurve()


setMethodS3("backtransformPrincipalCurve", "numeric", function(X, ...) {
  X <- as.matrix(X);
  backtransformPrincipalCurve(X, ...);
})

###########################################################################
# HISTORY:
# 2013-04-18
# o BUG FIX: backtransformPrincipalCurve() gave an error if the
#   pricipal curve was fitted using data with missing values.
#   Now backtransformPrincipalCurve() preserves dimension names.
# 2009-05-29
# o BUG FIX: Previous bug fix in backtransformPrincipalCurve() regarding
#   argument 'dimension' broke the initial purpose of this argument. Since
#   both use cases are still of interest, how the subsetting is done is now
#   based on whether the number of dimensions of the input data and the
#   model fit match or not. See help(backtransformPrincipalCurve.matrix).
#   Added several cases to the example code for testing this.
# o Added more Rdoc comments.
# 2009-05-12
# o BUG FIX: backtransformPrincipalCurve(..., dimensions) did not subset
#   the 'X' matrix. Also, the method now returns a matrix of the same
#   number of columns requested.  The Rd example now illustrates this.
#   Thanks to Pierre Neuvial, UC Berkeley for the troublshooting and fix.
# 2009-02-08
# o An error was thrown in backtransformPrincipalCurve() if argument
#   'dimensions' was specified.
# o BUG FIX:
# 2009-01-12
# o Updated validation of arguments such that it does not require R.utils.
# 2008-10-08
# o Added argument 'targetDimension' to backtransformPrincipalCurve().
# 2008-10-07
# o Created.
###########################################################################
