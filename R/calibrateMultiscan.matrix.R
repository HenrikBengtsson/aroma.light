#########################################################################/**
# @set "class=matrix"
# @RdocMethod calibrateMultiscan
#
# \encoding{latin1}
#
# @title "Weighted affine calibration of a multiple re-scanned channel"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An NxK @matrix (K>=2) where the columns represent the
#    multiple scans of one channel (a two-color array contains two
#    channels) to be calibrated.}
#  \item{weights}{If @NULL, non-weighted normalization is done.
#    If data-point weights are used, this should be a @vector of length 
#    N of data point weights used when estimating the normalization 
#    function.
#  }
#  \item{typeOfWeights}{A @character string specifying the type of 
#    weights given in argument \code{weights}.
#  }
#  \item{method}{A @character string specifying how the estimates are 
#    robustified.  See @seemethod "iwpca" for all accepted values.}
#  \item{constraint}{Constraint making the bias parameters identifiable.
#    See @seemethod "fitIWPCA" for more details.}
#  \item{satSignal}{Signals equal to or above this threshold is considered
#    saturated signals.}
#  \item{...}{Other arguments passed to @seemethod "fitIWPCA" and in
#   turn @seemethod "iwpca", e.g. \code{center} (see below).}
#  \item{average}{A @function to calculate the average signals between calibrated scans.}
#  \item{deviance}{A @function to calculate the deviance of the signals between calibrated scans.}
#  \item{project}{If @TRUE, the calibrated data points projected onto the
#    diagonal line, otherwise not. Moreover, if @TRUE, argument
#    \code{average} is ignored.}
#  \item{.fitOnly}{If @TRUE, the data will not be back-transform.}
# }
#
# \value{
#   If \code{average} is specified or \code{project} is @TRUE, 
#   an Nx1 @matrix is returned, otherwise an NxK @matrix is returned.
#   If \code{deviance} is specified, a deviance Nx1 @matrix is returned
#   as attribute \code{deviance}.
#   In addition, the fitted model is returned as attribute \code{modelFit}.
# }
#
# \section{Negative, non-positive, and saturated values}{
#   Affine multiscan calibration applies also to negative values, which are
#   therefor also calibrated, if they exist. 
#
#   Saturated signals in any scan are set to @NA. Thus, they will not be
#   used to estimate the calibration function, nor will they affect an
#   optional projection.
# }
#
# \section{Missing values}{
#   Only observations (rows) in \code{X} that contain all finite values are
#   used in the estimation of the alibration functions. Thus,
#   observations can be excluded by setting them to @NA.
# }
#
# \section{Weighted normalization}{
#  Each data point/observation, that is, each row in \code{X}, which is a
#  vector of length K, can be assigned a weight in [0,1] specifying how much
#  it should \emph{affect the fitting of the calibration function}. 
#  Weights are given by argument \code{weights},
#  which should be a @numeric @vector of length N. Regardless of weights, 
#  all data points are \emph{calibrated} based on the fitted calibration
#  function.
# }
# 
# \section{Robustness}{
#  By default, the model fit of multiscan calibration is done in \eqn{L_1}
#  (\code{method="L1"}). This way, outliers affect the parameter estimates
#  less than ordinary least-square methods.
#
#  When calculating the average calibrated signal from multiple scans,
#  by default the median is used, which further robustify against outliers.
#
#  For further robustness, downweight outliers such as saturated signals, 
#  if possible.
#
#  Tukey's biweight function is supported, but not used by default because 
#  then a "bandwidth" parameter has to selected. This can indeed be done
#  automatically by estimating the standard deviation, for instance using
#  MAD. However, since scanner signals have heteroscedastic noise 
#  (standard deviation is approximately proportional to the non-logged 
#  signal), Tukey's bandwidth parameter has to be a function of the 
#  signal too, cf. @see "stats::loess".  We have experimented with this 
#  too, but found that it does not significantly improve the robustness 
#  compared to \eqn{L_1}.
#  Moreover, using Tukey's biweight as is, that is, assuming homoscedastic
#  noise, seems to introduce a (scale dependent) bias in the estimates 
#  of the offset terms.
# }
#
# \section{Using a known/previously estimated offset}{
#  If the scanner offsets can be assumed to be known, for instance,
#  from prior multiscan analyses on the scanner, then it is possible
#  to fit the scanner model with no (zero) offset by specifying
#  argument \code{center=FALSE}.
#  Note that you cannot specify the offset.  Instead, subtract it
#  from all signals before calibrating, e.g.
#  \code{Xc <- calibrateMultiscan(X-e, center=FALSE)}
#  where \code{e} is the scanner offset (a scalar).
#  You can assert that the model is fitted without offset by
#  \code{stopifnot(all(attr(Xc, "modelFit")$adiag == 0))}.
# }
#
# \details{
#  Fitting is done by iterated re-weighted principal component analysis
#  (IWPCA).
# }
#
# @author
#
# \references{
#   [1] @include "../incl/BengtssonH_etal_2004.bib.Rdoc" \cr
# }
#
# \examples{\dontrun{# For an example, see help(normalizeAffine).}}
#
# \seealso{
#   @see "1. Calibration and Normalization".
#   @seemethod "normalizeAffine".
#   @seeclass
# }
#*/######################################################################### 
setMethodS3("calibrateMultiscan", "matrix", function(X, weights=NULL, typeOfWeights=c("datapoint"), method="L1", constraint="diagonal", satSignal=2^16-1, ..., average=median, deviance=NULL, project=FALSE, .fitOnly=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'X'
  if (ncol(X) < 2)
    stop("Multiscan calibratation requires at least two scans: ", ncol(X));
  if (nrow(X) < 3)
    stop("Multiscan calibratation requires at least three observations: ", nrow(X));
  
  # Argument: 'satSignal'
  if (satSignal < 0)
    stop("Argument 'satSignal' is negative: ", satSignal);

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

    weights <- as.vector(weights);

    if (length(weights) == 1) {
      weights <- rep(weights, length.out=nrow(X));
    } else if (length(weights) != nrow(X)) {
      stop("Argument 'weights' does not have the same length as the number of data points (rows) in the matrix: ", length(weights), " != ", nrow(X));
    }
    datapointWeights <- weights;
  }


  # Argument 'average':
  if (!is.null(average) && !is.function(average)) {
    throw("Argument 'average' must be a function or NULL: ", class(average)[1]);
  }

  # Argument 'deviance':
  if (!is.null(deviance) && !is.function(deviance)) {
    throw("Argument 'deviance' must be a function or NULL: ", class(deviance)[1]);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # 2. Prepare the data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Use non-saturated observations (non-finite values are taken care of by
  # the fitIWPCASpatial() function.
  X[(X >= satSignal)] <- NA;
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # 3. Fit the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  fit <- fitIWPCA(X, w=datapointWeights, method=method, constraint=constraint, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # 4. Backtransform
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.fitOnly == FALSE) {
    X <- backtransformAffine(X, a=fit, project=project);
    if (project == FALSE && !is.null(average)) {
      X <- apply(X, MARGIN=1, FUN=average, na.rm=TRUE);
      X <- as.matrix(X);
    }
    if (!is.null(deviance)) {
      deviance <- apply(X, MARGIN=1, FUN=deviance, na.rm=TRUE);
      attr(X, "deviance") <- as.matrix(deviance);
    }
  }
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # 5. Return the backtransformed data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  attr(X, "modelFit") <- fit;
  X;
}) # calibrateMultiscan()


############################################################################
# HISTORY:
# 2011-02-05
# o DOCUMENTATION: Added section on how to calibrate when scanner offsets
#   are supposed to be known/zero.
# o DOCUMENTATION: Fixed broken links to help for iwpca().
# 2005-06-03
# o Added argument 'typeOfWeights' to make it similar to other normalization
#   methods, although only "datapoint" weights are allowed.
# 2005-02-13
# o Made argument 'method="L1"' explicit and wrote a Rdoc comment about it
#   to document the fact that we have deliberately choosen not to use
#   "symmetric" Tukey's biweight.
# 2005-02-04
# o Put arguments 'average' and 'deviance' back again. It is much more
#   userfriendly. Averaging with median() is now the default.
# 2005-02-01
# o Added argument '.fitOnly'.
# 2005-01-24
# o Added argument 'weights' (instead of passing 'w' to fitIWPCA()).
# o Saturated values are not used to estimate the calibration function nor
#   are the used if data is projected.
# 2004-12-28
# o Added Rdoc comments on weights.
# 2004-06-28
# o BUG FIX: Missing braces in Rdoc comments.
# 2004-05-18
# o Removed averaging etc. That is now in its own function rowAverages().
# o The only difference between calibrateMultiscanSpatial() and
#   calibrateMultiscan() is how the parameters are fitted.
# 2004-05-14
# o Cleaned up. Making use of new backtransformAffine(), which makes the
#   code clearer. Explicit arguments that were just passed to iwpca() etc
#   are now passed as "..." to make the documentation simpler and less
#   confusing for the end user. Experts will follow "..." to iwpca().
############################################################################
