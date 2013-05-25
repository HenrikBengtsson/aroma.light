#########################################################################/**
# @set "class=matrix"
# @RdocMethod fitXYCurve
# @alias fitXYCurve
#
# @title "Fitting a smooth curve through paired (x,y) data"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An Nx2 @matrix where the columns represent the two channels
#    to be normalized.}
#  \item{weights}{If @NULL, non-weighted normalization is done.
#    If data-point weights are used, this should be a @vector of length
#    N of data point weights used when estimating the normalization
#    function.
#  }
#  \item{typeOfWeights}{A @character string specifying the type of
#    weights given in argument \code{weights}.
#  }
#  \item{method}{@character string specifying which method to use when
#    fitting the intensity-dependent function.
#    Supported methods:
#     \code{"loess"} (better than lowess),
#     \code{"lowess"} (classic; supports only zero-one weights),
#     \code{"spline"} (more robust than lowess at lower and upper
#                      intensities; supports only zero-one weights),
#     \code{"robustSpline"} (better than spline).
#  }
#  \item{bandwidth}{A @double value specifying the bandwidth of the
#   estimator used.
#  }
#  \item{satSignal}{Signals equal to or above this threshold will not
#    be used in the fitting.
#  }
#  \item{...}{Not used.}
# }
#
# \value{
#   A named @list structure of class \code{XYCurve}.
# }
#
# \section{Missing values}{
#  The estimation of the function will only be made based on complete
#  non-saturated observations, i.e. observations that contains no @NA
#  values nor saturated values as defined by \code{satSignal}.
# }
#
# \section{Weighted normalization}{
#  Each data point, that is, each row in \code{X}, which is a
#  vector of length 2, can be assigned a weight in [0,1] specifying how much
#  it should \emph{affect the fitting of the normalization function}.
#  Weights are given by argument \code{weights}, which should be a @numeric
#  @vector of length N.
#
#  Note that the lowess and the spline method only support zero-one
#  \{0,1\} weights.
#  For such methods, all weights that are less than a half are set to zero.
# }
#
# \section{Details on loess}{
#  For @see "stats::loess", the arguments \code{family="symmetric"},
#  \code{degree=1}, \code{span=3/4},
#  \code{control=loess.control(trace.hat="approximate"},
#  \code{iterations=5}, \code{surface="direct")} are used.
# }
#
# @author "HB"
#
# \examples{
#  @include "../incl/fitXYCurve.matrix.Rex"
# }
#*/#########################################################################
setMethodS3("fitXYCurve", "matrix", function(X, weights=NULL, typeOfWeights=c("datapoint"), method=c("loess", "lowess", "spline", "robustSpline"), bandwidth=NULL, satSignal=2^16-1, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'X'
  if (ncol(X) != 2) {
    stop("Curve-fit normalization requires two channels only: ", ncol(X));
  }
  if (nrow(X) < 3) {
    stop("Curve-fit normalization requires at least three observations: ",
                                                                   nrow(X));
  }

  # Argument: 'satSignal'
  if (satSignal < 0) {
    stop("Argument 'satSignal' is negative: ", satSignal);
  }

  # Argument: 'method'
  method <- match.arg(method);
  zeroOneWeightsOnly <- (method %in% c("lowess", "spline"));

  # Argument: 'typeOfWeights'
  typeOfWeights <- match.arg(typeOfWeights);

  # Argument: 'weights'
  datapointWeights <- NULL;
  if (!is.null(weights)) {
    # If 'weights' is an object of a class with as.double(), cast it.
    weights <- as.double(weights);

    if (any(is.na(weights)))
      stop("Argument 'weights' must not contain NA values.");

    if (any(weights < 0 | weights > 1)) {
      stop("Argument 'weights' out of range [0,1]: ",
           paste(weights[weights < 0.0 | weights > 1.0], collapse=", "));
    }

    if (zeroOneWeightsOnly && any(weights > 0 & weights < 1)) {
      weights <- round(weights);
      warning("Weights were rounded to {0,1} since '", method, "' normalization supports only zero-one weights.");
    }

    weights <- as.vector(weights);
    if (length(weights) == 1) {
      weights <- rep(weights, length.out=nrow(X));
    } else if (length(weights) != nrow(X)) {
      stop("Argument 'weights' does not have the same length as the number of data points (rows) in the matrix: ", length(weights), " != ", nrow(X));
    }
    datapointWeights <- weights;
  } # if (!is.null(weights))


  # Argument: 'bandwidth'
  if (is.null(bandwidth)) {
    bandwidths <- c("loess"=0.75, "lowess"=0.3, "robustSpline"=0.75,
                    "spline"=0.75);
    bandwidth <- bandwidths[method];
  } else if (!is.numeric(bandwidth) || bandwidth <= 0 || bandwidth > 1) {
    stop("Argument 'bandwidth' must be in [0,1): ", bandwidth);
  } else if (length(bandwidth) != 1) {
    stop("Argument 'bandwidth' must be a scalar: ", paste(bandwidth, collapse=", "));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Prepare data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Use only positive non-saturated observations to estimate the
  # normalization function
  isValid <- (is.finite(X) & (X <= satSignal));
  isValid <- (isValid[,1] & isValid[,2]);
  X <- X[isValid,];
  x <- X[,1,drop=TRUE];
  y <- X[,2,drop=TRUE];

  if (!is.null(datapointWeights)) {
    datapointWeights <- datapointWeights[isValid];
  }

  # Not needed anymore
  X <- isValid <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Fit the curve
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "lowess") {
    keep <- if (!is.null(datapointWeights)) (datapointWeights > 0) else TRUE;
    x <- x[keep];
    y <- y[keep];
    # Not needed anymore
    keep <- NULL;
    fit <- lowess(x=x, y=y, f=bandwidth, ...);
    fit$predictY <- function(x) approx(fit, xout=x, ties=mean)$y;
  } else if (method == "loess") {
    fit <- loess(formula=y ~ x, weights=datapointWeights,
                 family="symmetric", degree=1, span=bandwidth,
                 control=loess.control(trace.hat="approximate",
                 iterations=5, surface="direct"), ...);

    fit$predictY <- function(x) predict(fit, newdata=x);
  } else if (method == "spline") {
    keep <- if (!is.null(datapointWeights)) (datapointWeights > 0) else TRUE;
    x <- x[keep];
    y <- y[keep];
    # Not needed anymore
    keep <- NULL;
    fit <- smooth.spline(x=x, y=y, spar=bandwidth, ...);
    fit$predictY <- function(x) predict(fit, x=x)$y;
  } else if (method == "robustSpline") {
    fit <- robustSmoothSpline(x=x, y=y, w=datapointWeights, spar=bandwidth, ...);
    fit$predictY <- function(x) predict(fit, x=x)$y;
  }

  class(fit) <- c("XYCurveFit", class(fit));

  fit;
}) # fitXYCurve()



############################################################################
# HISTORY:
# 2009-07-15
# o Created from normalizeCurveFit.R.
############################################################################
