###########################################################################/**
# @set "class=matrix"
# @RdocMethod normalizeAverage
# @alias normalizeAverage.list
#
# @title "Rescales channel vectors to get the same average"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{A @numeric NxK @matrix (or @list of length K).}
#   \item{baseline}{An @integer in [1,K] specifying which channel should be
#      the baseline.}
#   \item{avg}{A @function for calculating the average of one channel.}
#   \item{targetAvg}{The average that each channel should have afterwards.
#      If @NULL, the baseline column sets the target average.}
#   \item{...}{Additional arguments passed to the \code{avg} @function.}
# }
#
# \value{
#  Returns a normalized @numeric NxK @matrix (or @list of length K).
# }
#
# @author "HB"
#*/###########################################################################
setMethodS3("normalizeAverage", "matrix", function(x, baseline=1, avg=median, targetAvg=2200, ...) {
  # Estimate the scale for each channel
  scale <- apply(x, MARGIN=2, FUN=avg, ...);

  # The scale of the baseline column
  scale1 <- scale[baseline];

  # Standardize so that the 'baseline' column is not rescaled (has scale one).
  scale <- scale / scale1;

  # Rescale to target averages?
  if (!is.null(targetAvg)) {
    rho <- (scale1 / targetAvg);
    scale <- rho * scale;
  }

  # Rescale so that all channels have the same scale
  for (cc in 1:ncol(x)) {
    x[,cc] <- x[,cc] / scale[cc];
  }

  x;
}, private=TRUE)


setMethodS3("normalizeAverage", "list", function(x, baseline=1, avg=median, targetAvg=2200, ...) {
  # Estimate the scale for each channel
  scale <- lapply(x, FUN=avg, ...);
  scale <- unlist(scale, use.names=FALSE);
  scale1 <- scale[baseline];

  # Standardize so that the 'baseline' channel has scale one.
  scale <- scale / scale1;

  # Rescale to target averages?
  if (!is.null(targetAvg)) {
    rho <- (scale1 / targetAvg);
    scale <- rho * scale;
  }

  # Rescale so that all channels have the same scale
  for (cc in 1:length(x)) {
    x[[cc]] <- x[[cc]] / scale[cc];
  }

  x;
}, private=TRUE)




############################################################################
# HISTORY:
# 2007-06-04
# o Corrected minor ineffective typo in code.
# 2007-03-29
# o Added Rdoc comments.
# 2006-05-08
# o Created.
############################################################################
