#########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizeAffine
#
# \encoding{latin1}
#
# @title "Weighted affine normalization between channels and arrays"
#
# \description{
#   @get "title".
#
#   This method will remove curvature in the M vs A plots that are
#   due to an affine transformation of the data. In other words, if there
#   are (small or large) biases in the different (red or green) channels,
#   biases that can be equal too, you will get curvature in the M vs A plots
#   and this type of curvature will be removed by this normalization method.
#
#   Moreover, if you normalize all slides at once, this method will also
#   bring the signals on the same scale such that the log-ratios for
#   different slides are comparable. Thus, do not normalize the scale of
#   the log-ratios between slides afterward.
#
#   It is recommended to normalize as many slides as possible in one run.
#   The result is that if creating log-ratios between any channels and any
#   slides, they will contain as little curvature as possible.
#
#   Furthermore, since the relative scale between any two channels on any
#   two slides will be one if one normalizes all slides (and channels) at
#   once it is possible to add or multiply with the \emph{same} constant
#   to all channels/arrays without introducing curvature. Thus, it is
#   easy to rescale the data afterwards as demonstrated in the example.
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An NxK @matrix (K>=2) where the columns represent the channels,
#    to be normalized.}
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
#  \item{satSignal}{Signals equal to or above this threshold will not
#    be used in the fitting.}
#  \item{...}{Other arguments passed to @seemethod "fitIWPCA" and in
#   turn @seemethod "iwpca". For example, the weight argument
#   of @seemethod "iwpca".  See also below.}
#  \item{.fitOnly}{If @TRUE, the data will not be back-transform.}
# }
#
# \value{
#   A NxK @matrix of the normalized channels.
#   The fitted model is returned as attribute \code{modelFit}.
# }
#
# \section{Negative, non-positive, and saturated values}{
#  Affine normalization applies equally well to negative values. Thus,
#  contrary to normalization methods applied to log-ratios, such as curve-fit
#  normalization methods, affine normalization, will not set these to @NA.
#
#  Data points that are saturated in one or more channels are not used
#  to estimate the normalization function, but they are normalized.
# }
#
# \section{Missing values}{
#  The estimation of the affine normalization function will only be made
#  based on complete non-saturated observations, i.e. observations that
#  contains no @NA values nor saturated values as defined by \code{satSignal}.
# }
#
# \section{Weighted normalization}{
#  Each data point/observation, that is, each row in \code{X}, which is a
#  vector of length K, can be assigned a weight in [0,1] specifying how much
#  it should \emph{affect the fitting of the affine normalization function}.
#  Weights are given by argument \code{weights},
#  which should be a @numeric @vector of length N. Regardless of weights,
#  all data points are \emph{normalized} based on the fitted normalization
#  function.
# }
#
# \section{Robustness}{
#  By default, the model fit of affine normalization is done in \eqn{L_1}
#  (\code{method="L1"}). This way, outliers affect the parameter estimates
#  less than ordinary least-square methods.
#
#  For further robustness, downweight outliers such as saturated signals,
#  if possible.
#
#  We do not use Tukey's biweight function for reasons similar to those
#  outlined in @seemethod "calibrateMultiscan".
# }
#
# \section{Using known/previously estimated channel offsets}{
#  If the channel offsets can be assumed to be known, then it is
#  possible to fit the affine model with no (zero) offset, which
#  formally is a linear (proportional) model, by specifying
#  argument \code{center=FALSE}.
#  In order to do this, the channel offsets have to be subtracted
#  from the signals manually before normalizing, e.g.
#  \code{Xa <- t(t(X)-a)} where \code{e} is @vector of length
#  \code{ncol(X)}.  Then normalize by
#  \code{Xn <- normalizeAffine(Xa, center=FALSE)}.
#  You can assert that the model is fitted without offset by
#  \code{stopifnot(all(attr(Xn, "modelFit")$adiag == 0))}.
# }
#
# \details{
#  A line is fitted robustly throught the \eqn{(y_R,y_G)} observations
#  using an iterated re-weighted principal component analysis (IWPCA),
#  which minimized the residuals that are orthogonal to the fitted line.
#  Each observation is down-weighted by the inverse of the absolute
#  residuals, i.e. the fit is done in \eqn{L_1}.
# }
#
# @author
#
# \references{
#   [1] @include "../incl/BengtssonHossjer_2006.bib.Rdoc" \cr
# }
#
# @examples "../incl/normalizeCurveFit.matrix.Rex"
#
# \seealso{
#   @seemethod "calibrateMultiscan".
# }
#*/#########################################################################
setMethodS3("normalizeAffine", "matrix", function(X, weights=NULL, typeOfWeights=c("datapoint"), method="L1", constraint=0.05, satSignal=2^16-1, ..., .fitOnly=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'X'
  if (ncol(X) < 2)
    stop("Affine normalization requires at least two channels: ", ncol(X));
  if (nrow(X) < 3)
    stop("Affine normalization requires at least three observations: ", nrow(X));

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
      stop("Argument 'weights' does not have the same length as the number of data points (rows9 in the matrix: ", length(weights), " != ", nrow(X));
    }
    datapointWeights <- weights;
  }


  # Argument: 'method'
  # Validate by fitIWPCA() -> iwpca()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Prepare data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Use only non-saturated observations to estimate the normalization
  # function (non-finite values are taken care of by fitIWPCA()).
  isSaturated <- (is.finite(X) & X >= satSignal);
  Xsat <- X[isSaturated];
  X[isSaturated] <- NA;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Fit the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit <- fitIWPCA(X, w=datapointWeights, method=method, constraint=constraint, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 4. Backtransform
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  X[isSaturated] <- Xsat;
  # Not needed anymore
  isSaturated <- Xsat <- NULL;

  if (.fitOnly == FALSE) {
    X <- backtransformAffine(X, a=fit);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 5. Return the backtransformed data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  attr(X, "modelFit") <- fit;
  X;
}) # normalizeAffine()


############################################################################
# HISTORY:
# 2011-02-05
# o DOCUMENTATION: Added section on how to normalize when channel offsets
#   are supposed to be known/zero.
# o DOCUMENTATION: Fixed broken links to help for iwpca().
# 2005-06-03
# o Added argument 'typeOfWeights' to make it similar to other normalization
#   methods, although only "datapoint" weights are allowed.
# 2005-02-27
# o Passes argument 'methods' to fitIWPCA() now.
# 2005-02-02
# o BUG FIX: isSaturated could contain NA.
# 2005-02-01
# o Added argument '.fitOnly'.
# o Added validation of argument 'weights'.
# 2005-01-24
# o Added argument 'weights' (instead of passing 'w' to fitIWPCA()).
# o Now, saturated functions are normalized, but just not used when
#   estimating the normalization function.
# 2005-01-23
# o Updated the Rdoc comments and error messages.
# 2004-12-28
# o Added Rdoc comments on weights.
# 2004-06-28
# o BUG FIX: Missing braces in Rdoc comments.
############################################################################
