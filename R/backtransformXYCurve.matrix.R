setMethodS3("backtransformXYCurve", "matrix", function(X, fit, targetChannel=1, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'X'
  if (ncol(X) != 2) {
    stop("Curve-fit normalization requires two channels only: ", ncol(X));
  }

  # Argument 'targetChannel':
# targetChannel <- Arguments$getIndex(targetChannel, range=c(1,ncol(X)));

  # Allocate result
  XN <- X;

  # Predict using only finite covariates (otherwise an error)
  keep <- which(is.finite(X[,targetChannel]));

  # Nothing to do?
  if (length(keep) > 0) {
    X <- X[keep,,drop=FALSE];

    xN <- fit$predictY(X[,targetChannel]);
    delta <- xN - X[,targetChannel];

    # Not needed anymore
    xN <- NULL;

    XN[keep,-targetChannel] <- X[,-targetChannel] - delta;

    # Not needed anymore
    keep <- delta <- NULL;
  }

  XN;
})

############################################################################
# HISTORY:
# 2009-07-15
# o Created.
############################################################################
