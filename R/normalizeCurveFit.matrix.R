#########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizeCurveFit
# @alias "normalizeLoess.matrix"
# @alias "normalizeLowess.matrix"
# @alias "normalizeSpline.matrix"
# @alias "normalizeRobustSpline.matrix"
#
# \encoding{latin1}
#
# @title "Weighted curve-fit normalization between a pair of channels"
#
# \description{
#   @get "title".
#
#   This method will estimate a smooth function of the dependency
#   between the log-ratios and the log-intensity of the two channels and
#   then correct the log-ratios (only) in order to remove the dependency.
#   This is method is also known as \emph{intensity-dependent} or
#   \emph{lowess normalization}.
#
#   The curve-fit methods are by nature limited to paired-channel data.
#   There exist at least one method trying to overcome this limitation,
#   namely the cyclic-lowess [1], which applies the paired
#   curve-fit method iteratively over all pairs of channels/arrays.
#   Cyclic-lowess is not implented here.
#
#   We recommend that affine normalization [2] is used instead of curve-fit
#   normalization.
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
#   A Nx2 @matrix of the normalized two channels.
#   The fitted model is returned as attribute \code{modelFit}.
# }
#
# \details{
#  A smooth function \eqn{c(A)} is fitted throught data in \eqn{(A,M)},
#  where \eqn{M=log_2(y_2/y_1)} and \eqn{A=1/2*log_2(y_2*y_1)}. Data is
#  normalized by \eqn{M <- M - c(A)}.
#
#  Loess is by far the slowest method of the four, then lowess, and then
#  robust spline, which iteratively calls the spline method.
# }
#
# \section{Negative, non-positive, and saturated values}{
#  Non-positive values are set to not-a-number (@NaN).
#  Data points that are saturated in one or more channels are not used
#  to estimate the normalization function, but they are normalized.
# }
#
# \section{Missing values}{
#  The estimation of the normalization function will only be made
#  based on complete non-saturated observations, i.e. observations that
#  contains no @NA values nor saturated values as defined by \code{satSignal}.
# }
#
# \section{Weighted normalization}{
#  Each data point, that is, each row in \code{X}, which is a
#  vector of length 2, can be assigned a weight in [0,1] specifying how much
#  it should \emph{affect the fitting of the normalization function}.
#  Weights are given by argument \code{weights}, which should be a @numeric
#  @vector of length N. Regardless of weights, all data points are
#  \emph{normalized} based on the fitted normalization function.
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
# \references{
#   [1] M. \enc{Åstrand}{Astrand},
#       Contrast Normalization of Oligonucleotide Arrays,
#       Journal Computational Biology, 2003, 10, 95-102. \cr
#   [2] @include "../incl/BengtssonHossjer_2006.bib.Rdoc" \cr
# }
#
# \examples{
#  @include "../incl/normalizeCurveFit.matrix.Rex"
# }
#
# \seealso{
#   @seemethod "normalizeAffine".
# }
#*/#########################################################################
setMethodS3("normalizeCurveFit", "matrix", function(X, weights=NULL, typeOfWeights=c("datapoint"), method=c("loess", "lowess", "spline", "robustSpline"), bandwidth=NULL, satSignal=2^16-1, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'X'
  if (ncol(X) != 2)
    stop("Curve-fit normalization requires two channels only: ", ncol(X));
  if (nrow(X) < 3)
    stop("Curve-fit normalization requires at least three observations: ",
                                                                   nrow(X));

  # Argument: 'satSignal'
  if (satSignal < 0)
    stop("Argument 'satSignal' is negative: ", satSignal);

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
  # Convert non-positive signals to NaN. If not done here, the transform
  # (R,G) -> (A,M) -> (R,G) will no it.
  X[X <= 0] <- NaN;

  # Use only positive non-saturated observations to estimate the
  # normalization function
  isValid <- (is.finite(X) & (X <= satSignal));
  isValid <- (isValid[,1] & isValid[,2]);
  Y <- X[isValid,];

  if (!is.null(datapointWeights))
    datapointWeights <- datapointWeights[isValid];

  M <- log(Y[,1]/Y[,2], base=2);
  A <- log(Y[,1]*Y[,2], base=2)/2;
  rm(Y);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Estimate the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "lowess") {
    incl <- if (!is.null(datapointWeights)) (datapointWeights > 0) else TRUE;
    fit <- lowess(x=A[incl], y=M[incl], f=bandwidth, ...);
    fit$predictM <- function(newA) approx(fit, xout=newA, ties=mean)$y;
  } else if (method == "loess") {
    fit <- loess(formula=M ~ A, weights=datapointWeights,
                 family="symmetric", degree=1, span=bandwidth,
                 control=loess.control(trace.hat="approximate",
                 iterations=5, surface="direct"), ...);

    fit$predictM <- function(newA) predict(fit, newdata=newA);
  } else if (method == "spline") {
    incl <- if (!is.null(datapointWeights)) (datapointWeights > 0) else TRUE;
    fit <- smooth.spline(x=A[incl], y=M[incl], spar=bandwidth, ...);
    fit$predictM <- function(newA) predict(fit, x=newA)$y;
  } else if (method == "robustSpline") {
    fit <- robustSmoothSpline(x=A, y=M, w=datapointWeights, spar=bandwidth, ...);
    fit$predictM <- function(newA) predict(fit, x=newA)$y;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 4. Normalize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize all data
  M <- log(X[,1]/X[,2], base=2);
  A <- log(X[,1]*X[,2], base=2)/2;

  ok <- is.finite(A);
  M[ok] <- M[ok] - fit$predictM(A[ok]);
  rm(ok);
  X[,1] <- as.matrix(sqrt(2^(2*A+M)));
  X[,2] <- as.matrix(sqrt(2^(2*A-M)));

  rm(A,M);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 5. Return the normalized data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  attr(X, "modelFit") <- fit;
  X;
}) # normalizeCurveFit()


setMethodS3("normalizeLowess", "matrix", function(X, ...) {
  normalizeCurveFit(X, method="lowess", ...);
})

setMethodS3("normalizeLoess", "matrix", function(X, ...) {
  normalizeCurveFit(X, method="loess", ...);
})

setMethodS3("normalizeSpline", "matrix", function(X, ...) {
  normalizeCurveFit(X, method="spline", ...);
})

setMethodS3("normalizeRobustSpline", "matrix", function(X, ...) {
  normalizeCurveFit(X, method="robustSpline", ...);
})


############################################################################
# HISTORY:
# 2005-06-03
# o Added argument 'typeOfWeights' to make it similar to other normalization
#   methods, although only "datapoint" weights are allowed.
# o Removed argument '.fitOnly'.
# o renamed all "robust.spline" to "robustSpline".
# 2005-03-23
# o Updated normalizeCurveFit() so that approx() does not give warnings
#   about 'Collapsing to unique x values' when doing lowess normalization.
# 2005-02-02
# o Zero-one weights are now round off by round(w).
# o BUG FIX: Forgot to adjust weights vector in normalizeCurveFit() when
#   removing non-finite values from data.
# 2005-02-01
# o Added argument '.fitOnly'.
# 2005-01-24
# o Create an Rdoc example with MvsA and MvsM comparisons.
# 2005-01-23
# o Added aliases normalizeLowess() and normalizeLoess().
# o Created from normalizeCurveFit() in MAData().
############################################################################
