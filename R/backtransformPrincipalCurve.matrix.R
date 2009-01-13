#########################################################################/**
# @set "class=matrix"
# @RdocMethod backtransformPrincipalCurve
# @alias backtransformPrincipalCurve.numeric
#
# @title "Reverse affine transformation"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An NxK @matrix containing data to be backtransformed.}
#  \item{fit}{An object of class \code{principal.curve} as returned by
#    @seemethod "fitPrincipalCurve".}
#  \item{dimensions}{An (optional) subset of of dimensions all in [1,K]
#    to be returned (and backtransform).}
#  \item{targetDimension}{An (optional) index specifying the dimension
#    in [1,K] to be used as the target dimension.  All other
#    dimensions will be normalized toward this dimension, which will
#    not be changed.}
#  \item{...}{Passed internally to @see "stats::smooth.spline".}
# }
#
# \value{
#   The backtransformed NxK @matrix.
# }
#
# \examples{\dontrun{See help(fitPrincipalCurve.matrix).}}
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

  dim <- dim(X);
  if (!is.matrix(X)) {
    X <- as.matrix(X);
  }

  # Argument 'fit'
  if (!inherits(fit, "principal.curve")) {
    stop("Argument 'fit' is not a principal.curve object: ", class(fit)[1]);
  }

  # Argument 'dimensions'
  p <- ncol(fit$s);
  if (!is.null(dimensions)) {
    dimensions <- as.integer(dimensions);
    if (any(dimensions < 1 | dimension > p)) {
      stop("Argument 'dimensions' contains values out of range [1,", p, "]");
    }
  }

  if (!is.null(targetDimension)) {
    targetDimension <- as.integer(targetDimension);
    if (length(targetDimension) != 1) {
      stop("Argument 'targetDimension' should be a scalar or NULL.");
    }
    if (targetDimension < 1 | targetDimension > p) {
      stop("Argument 'targetDimension' is out of range [1,", p, "]: ", 
                                                           targetDimension);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Pre-allocate result matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  naValue <- NA;
  mode(naValue) <- mode(X);
  Xhat <- matrix(naValue, nrow=dim[1], ncol=dim[2]);


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
  # Find backtransformations and backtransform data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  s <- fit$s;
  if (!is.null(dimensions)) {
    s <- s[,dimensions,drop=FALSE];
  }

    
  for (kk in seq(length=ncol(s))) {
    sKK <- s[,kk];
    fitKK <- smooth.spline(sKK, lambda, ...);

    Xkk <- X[,kk];
    keep <- which(is.finite(Xkk));
    Xkk <- Xkk[keep];
    XhatKK <- predict(fitKK, x=Xkk)$y;

    # Sanity check
    stopifnot(length(XhatKK) == length(keep));

    Xhat[keep,kk] <- XhatKK;
  }

  rm(sKK, lambda, fitKK, XhatKK, keep, s);

  dim(Xhat) <- dim;
  Xhat;
}) # backtransformPrincipalCurve()


setMethodS3("backtransformPrincipalCurve", "numeric", function(X, ...) {
  X <- as.matrix(X);
  backtransformPrincipalCurve(X, ...);
})

###########################################################################
# HISTORY:
# 2009-01-12
# o Updated validation of arguments such that it does not require R.utils.
# 2008-10-08
# o Added argument 'targetDimension' to backtransformPrincipalCurve().
# 2008-10-07
# o Created.
###########################################################################
