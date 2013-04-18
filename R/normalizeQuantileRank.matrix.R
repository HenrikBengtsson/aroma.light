###########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizeQuantileRank
# @aliasmethod normalizeQuantile
#
# @title "Weighted sample quantile normalization"
#
# @synopsis
#
# \description{
#   Normalizes channels so they all have the same average sample
#   distributions.
#
#   The average sample distribution is calculated either robustly or not
#   by utilizing either \code{weightedMedian()} or \code{weighted.mean()}.
#   A weighted method is used if any of the weights are different from one.
# }
#
# \arguments{
#   \item{X}{a numerical NxK @matrix with the K columns representing the
#     channels and the N rows representing the data points.}
#   \item{robust}{If @TRUE, the (weighted) median function is used for
#            calculating the average sample distribution, otherwise the
#            (weighted) mean function is used.}
#   \item{ties}{Should ties be specially treated or not?}
#   \item{weights}{If @NULL, non-weighted normalization is done.
#     If channel weights, this should be a @vector of length K specifying
#     the weights for each channel.
#     If signal weights, it should be an NxK @matrix specifying the
#     weights for each signal.
#   }
#   \item{typeOfWeights}{A @character string specifying the type of
#     weights given in argument \code{weights}.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an NxK @matrix.
# }
#
# \section{Missing values}{
#   Missing values are excluded when estimating the "common" (the baseline)
#   distribution. Values that are @NA before remain @NA. No new @NAs are
#   introduced.
# }
#
# \section{Weights}{
#   Currently only channel weights are support due to the way quantile
#   normalization is done.
#   If signal weights are given, channel weights are calculated from these
#   by taking the mean of the signal weights in each channel.
# }
#
# @examples "../incl/normalizeQuantileRank.matrix.Rex"
#
# \author{
#   Adopted from Gordon Smyth (\url{http://www.statsci.org/}) in 2002 \& 2006.
#   Original code by Ben Bolstad at Statistics Department, University of
#   California.
#   Support for calculating the average sample distribution using (weighted)
#   mean or median was added by Henrik Bengtsson.
# }
#
# \seealso{
#   @see "stats::median", @see "matrixStats::weightedMedian",
#   @see "base::mean" and @see "stats::weighted.mean".
#   @seemethod "normalizeQuantileSpline".
# }
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("normalizeQuantileRank", "matrix", function(X, ties=FALSE, robust=FALSE, weights=NULL, typeOfWeights=c("channel", "signal"), ...) {
  zeroOneWeightsOnly <- TRUE;  # Until supported otherwise.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfChannels <- ncol(X);
  if(nbrOfChannels == 1)
    return(X);

  nbrOfObservations <- nrow(X);
  if(nbrOfObservations == 1)
    return(X);

  # Argument 'typeOfWeights':
  typeOfWeights <- match.arg(typeOfWeights);

  # Argument 'weights':
  channelWeights <- NULL;
  signalWeights <- NULL;
  if (!is.null(weights)) {
    # If 'weights' is an object of a class with as.double(), cast it.
    dim <- dim(weights);
    weights <- as.double(weights);
    dim(weights) <- dim;

    if (any(is.na(weights)))
      stop("Argument 'weights' must not contain NA values.");

    if (any(weights < 0 | weights > 1)) {
      stop("Argument 'weights' out of range [0,1]: ",
           paste(weights[weights < 0.0 | weights > 1.0], collapse=", "));
    }

    if (typeOfWeights == "channel") {
      if (length(weights) == 1) {
        weights <- rep(weights, length.out=nbrOfObservations);
      } else if (length(weights) != nbrOfObservations) {
        stop("Argument 'weights' (channel weights) does not have the same length as the number of rows in the matrix: ", length(weights), " != ", nbrOfChannels);
      }

      channelWeights <- weights;
    } else if (typeOfWeights == "signal") {
      if (!identical(dim(weights), dim(X))) {
        stop("Dimension of argument 'weights' (signal weights) does not match dimension of argument 'X': (", paste(dim(weights), collapse=","), ") != (", paste(dim(X), collapse=","), ")");
      }

      # Calculate channel weights
      channelWeights <- apply(weights, MARGIN=2, FUN=mean);

      if (zeroOneWeightsOnly && any(weights > 0 & weights < 1)) {
        weights <- round(weights);
        warning("Weights were rounded to {0,1} since quantile normalization supports only zero-one weights.");
      }

      signalWeights <- weights;
    }
  } # if (!is.null(weights))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 0. Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  maxNbrOfObservations <- nbrOfObservations;

  # Create a list S to hold the sorted values for each channels
  naValue <- as.double(NA);
  S <- matrix(naValue, nrow=maxNbrOfObservations, ncol=nbrOfChannels);

  # Create a list O to hold the ordered indices for each channels
  O <- vector("list", nbrOfChannels);

  # A vector specifying the number of observations in each column
  nbrOfFiniteObservations <- rep(maxNbrOfObservations, times=nbrOfChannels);

  # Construct the sample quantiles
  quantiles <- (0:(maxNbrOfObservations-1))/(maxNbrOfObservations-1);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Get the sample quantile for all channels (columns)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (cc in 1:nbrOfChannels) {
    values <- X[,cc];

    if (!is.null(signalWeights)) {
      values[signalWeights[,cc] == 0] <- NA;
    }

    # Order and sort the values
    Scc <- sort(values, index.return=TRUE);

    # The number of non-NA observations
    nobs <- length(Scc$x);

    # Has NAs?
    if(nobs < nbrOfObservations) {
      nbrOfFiniteObservations[cc] <- nobs;
      isOk <- !is.na(values);

      # Get the sample quantiles for those values
      bins <- (0:(nobs-1))/(nobs-1);

      # Record the order position for these values.
      O[[cc]] <- ((1:nbrOfObservations)[isOk])[Scc$ix];

      # Interpolate to get the values at positions specified by
      # 'quantile' using data points given by 'bins' and 'Scc$x'.
      Scc <- approx(x=bins, y=Scc$x, xout=quantiles, ties="ordered")$y;
    } else {
      O[[cc]] <- Scc$ix;
      Scc <- Scc$x;
    }

    S[,cc] <- Scc;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Calculate the average sample distribution, of each quantile
  #    across all columns. This can be done robustly or not and
  #    with weights or not.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  useWeighted <- (!is.null(channelWeights) & any(channelWeights != 1));
  if (robust) {
    if (useWeighted) {
      xTarget <- apply(S, MARGIN=1, FUN=weightedMedian, w=channelWeights, na.rm=TRUE);
    } else {
      xTarget <- apply(S, MARGIN=1, FUN=median, na.rm=TRUE);
    }
  } else {
    if (useWeighted) {
      xTarget <- apply(S, MARGIN=1, FUN=weighted.mean, w=channelWeights, na.rm=TRUE);
    } else {
      xTarget <- rowMeans(S, na.rm=TRUE);
    }
  }

  # Assert that xTarget is of then same length as number of observations
  stopifnot(length(xTarget) == maxNbrOfObservations);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. For all columns, get for each sample quantile the value of
  #    average sample distribution at that quantile.
  #
  # Input: X[r,c], xTarget[r], O[[c]][r], nbrOfFiniteObservations[c].
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (cc in 1:nbrOfChannels) {
    # Get the number of non-NA observations
    nobs <- nbrOfFiniteObservations[cc];

    # Has NAs?
    if(nobs < nbrOfObservations) {
      # Get the NAs
      isOk <- !is.na(X[,cc]);

      # Get the sample quantiles for those values
      if (ties) {
        r <- rank(X[isOk,cc]);
        bins <- (r-1)/(nobs-1);
      } else {
        bins <- (0:(nobs-1))/(nobs-1);
      }

      # Interpolate to get the m's at positions specified by
      # 'quantile' using data points given by 'bins' and 'xTarget'.
      if (ties) {
        idx <- isOk;
      } else {
        idx <- O[[cc]];
      }
      X[idx,cc] <- approx(x=quantiles, y=xTarget, xout=bins, ties="ordered")$y;
    } else {
      if (ties) {
        r <- rank(X[,cc]);
        bins <- (r-1)/(nobs-1);
        X[,cc] <- approx(x=quantiles, y=xTarget, xout=bins, ties="ordered")$y;
      } else {
        X[O[[cc]],cc] <- xTarget;
      }
    }
  }

  X;
})





##############################################################################
# HISTORY:
# 2011-04-12
# o Now using as.double(NA) instead of NA.
# 2008-04-14
# o Renamed normalizeQuantile() to normalizeQuantileRank().  Keeping the old
#   name for backward compatibility.
# 2006-05-12
# o Updated according to Gordon Smyth's package.
# 2006-02-08
# o Rd bug fix: Fixed broken links to help for median() and weighted.mean().
# 2005-06-03
# o Replaced arguments 'channelWeights' and 'signalWeights' with 'weights'
#   and 'typeOfWeights'.
# o Added Rdoc help on weights.
# 2005-03-23
# o Updated normalizeQuantile() so that approx() does not give warnings
#   about 'Collapsing to unique x values' when doing lowess normalization.
# 2005-02-02
# o Zero-one weights are now round off by round(w).
# 2005-02-01
# o Added argument 'signalWeights'.
# o Added validation of argument 'channelWeights'.
# 2005-01-27
# o Renamed argument 'A' to 'X'.
# o Renamed argument 'weights' to 'channelWeights'.
# o Making use of setMethodS3(). Added some more Rdoc comments.
# o Moved from R.basic to aroma.
# 2003-04-13
# o Updated the Rdoc comments.
# 2002-10-24
# o Adapted from source code by Gordon Smyth's smagks package.
##############################################################################
