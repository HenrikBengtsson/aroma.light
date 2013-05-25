#########################################################################/**
# @set "class=matrix"
# @RdocMethod iwpca
#
# @title "Fits an R-dimensional hyperplane using iterative re-weighted PCA"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{N-times-K @matrix where N is the number of observations and
#    K is the number of dimensions.}
#  \item{w}{An N @vector of weights for each row (observation) in
#    the data matrix. If @NULL, all observations get the same weight.}
#  \item{R}{Number of principal components to fit. By default a line
#    is fitted.}
#  \item{method}{
#    If \code{"symmetric"} (or \code{"bisquare"}), Tukey's biweight
#    is used. If \code{"tricube"}, the tricube weight is used.
#    If \code{"L1"}, the model is fitted in \eqn{L_1}.
#    If a @function, it is used to calculate weights for next iteration
#    based on the current iteration's residuals.}
#  \item{maxIter}{Maximum number of iterations.}
#  \item{acc}{The (Euclidean) distance between two subsequent parameters
#    fit for which the algorithm is considered to have converged.}
#  \item{reps}{Small value to be added to the residuals before the
#    the weights are calculated based on their inverse. This is to avoid
#    infinite weights.}
#  \item{fit0}{A @list containing elements \code{vt} and \code{pc}
#    specifying an initial fit.
#    If @NULL, the initial guess will be equal to the (weighted) PCA fit.}
#  \item{...}{Additional arguments accepted by @seemethod "wpca".}
# }
#
# \value{
#   Returns the fit (a @list) from the last call to @seemethod "wpca"
#   with the additional elements \code{nbrOfIterations} and
#   \code{converged}.
# }
#
# \details{
#   This method uses weighted principal component analysis (WPCA) to fit a
#   R-dimensional hyperplane through the data with initial internal
#   weights all equal.
#   At each iteration the internal weights are recalculated based on
#   the "residuals".
#   If \code{method=="L1"}, the internal weights are 1 / sum(abs(r) + reps).
#   This is the same as \code{method=function(r) 1/sum(abs(r)+reps)}.
#   The "residuals" are orthogonal Euclidean distance of the principal
#   components R,R+1,...,K.
#   In each iteration before doing WPCA, the internal weighted are
#   multiplied by the weights given by argument \code{w}, if specified.
# }
#
# @author
#
# @examples "../incl/iwpca.matrix.Rex"
#
# \seealso{
#   Internally @seemethod "wpca" is used for calculating the weighted PCA.
# }
#
# @keyword "algebra"
#*/#########################################################################
setMethodS3("iwpca", "matrix", function(X, w=NULL, R=1, method=c("symmetric", "bisquare", "tricube", "L1"), maxIter=30, acc=1e-4, reps=0.02, fit0=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'w'
  if (!is.null(w))
    w <- rep(w, length.out=nrow(X));
  w0 <- w;

  # Argument: 'method'
  if (is.function(method)) {
    dummy <- method(1:5);
    if (!is.numeric(dummy))
      stop("Argument 'method' (weight function) does not return numeric values.");
    if (!is.vector(dummy))
      stop("Argument 'method' (weight function) does not return a vector.");
    if (length(dummy) != 5)
      stop("Argument 'method' (weight function) does not return the correct number of values.");
  } else {
    method <- match.arg(method);
  }

  # Argument: 'fit0'
  if (!is.null(fit0)) {
    if (!all(c("vt", "pc") %in% names(fit0))) {
      throw("Argument 'fit0' is missing element 'vt' or 'pc': ",
                                        paste(names(fit0), collapse=", "));
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Fit the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Ulast <- 1/.Machine$double.eps; # A large number
  iter <- 0;
  isConverged <- FALSE;
  w <- rep(1, length=nrow(X));
  while (!isConverged && iter < maxIter) {
    if (iter > 0 || is.null(fit0)) {
      iter <- iter + 1;

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # "Re-weight the weights"
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!is.null(w0))
        w <- w0 * w;

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Fit N-dimensional weighted PCA.
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      fit <- wpca(X, w=w, scale=FALSE, ...);
    } else {
      fit <- fit0;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the fitted line L
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get fitted eigenvectors (u1,u2,...,uN)
    # with ui*uj = 0; i!=j and ui*ui = 1.
    U <- fit$vt;
    colnames(U) <- rownames(U) <- NULL;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Check for convergence
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    isConverged <- (sum(abs(U-Ulast))/length(U) < acc);
    Ulast <- U;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Finally, update the weights
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Residuals in the "tailing" dimensions.
    r <- fit$pc[,-c(1:R), drop=FALSE];

    # Residuals in orthogonal Euclidean distance
    if (any(is.na(r))) {
      # Sometimes some residuals become NAs.
      r <- sqrt(rowSums(r^2, na.rm=TRUE));
    } else {
      r <- sqrt(rowSums(r^2));
    }

    # Down-weight points that are "far" away...
    if (is.character(method)) {
      if (method == "L1") {
        # Add small number to residuals to avoid infinite weights
        r <- abs(r) + reps;
        w <- 1/r;
      } else if (method %in% c("symmetric", "bisquare")) {
        # Add small number to residuals to avoid infinite weights
        r <- abs(r) + reps;
        r <- r/6;
        # Zero weights introduce NA's (for unknown reasons), therefore
        # with use a number very close to zero instead.
        w <- rep(.Machine$double.eps, length(r));
        ii <- (r < 1);
        w[ii] <- (1-r[ii]^2)^2;
        # Not needed anymore
        ii <- NULL;
      } else if (method == "tricube") {
        # Add small number to residuals to avoid infinite weights
        r <- abs(r) + reps;
        r <- r/6;
        # Zero weights introduce NA's (for unknown reasons), therefore
        # with use a number very close to zero instead.
        w <- rep(.Machine$double.eps, length(r));
        ii <- (r < 1);
        w[ii] <- (1-r[ii]^3)^3;
        # Not needed anymore
        ii <- NULL;
      }
    } else if (is.function(method)) {
      # Pass also the "fitted values" to the weight function.
      attr(r, "x") <- fit$pc[,1:R, drop=FALSE];
      attr(r, "reps") <- reps;
      w <- method(r);

      # Combine w_i = ||w_{i,j}||_2 (Euclidean distance), if needed.
      if (is.matrix(w) && ncol(w) > 1) {
        w <- sqrt(rowSums(w^2)/ncol(w));
      } else {
        w <- as.vector(w);
      }
    }

    # Sometimes some weights become NAs (not sure why). Set them to zero.
    if (any(is.na(w))) {
      w[is.na(w)] <- 0;
    }

    # Not needed anymore
    r <- NULL;
  } # while(...)

  fit$w <- w;
  fit$nbrOfIterations <- iter;
  fit$converged <- isConverged;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Return the estimated parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit;
}) # iwpca()


############################################################################
# HISTORY:
# 2006-04-25
# o Updated the example to first plot data from all viewpoints, then just
#   the lines.  Faster since the lines are only fitted once and nicer.
# 2005-05-03
# o Now test of argument 'fit0' checks if it is NULL.
# 2005-03-28
# o Second try with initial guess; added argument 'fit0'.
# 2005-02-08
# o Added "symmetric" (now default) and "tricube" too.
# 2005-02-07
# o Argument 'method' is now how the weights are calculated from the
#   residuals, not how residuals are combined across dimensions.
#   Method "L2" is therefore removed, because it corresponds to wpca().
# o Added support for weight functions via argument 'method'. The function
#   must take a matrix of residuals as the first argument.
# 2004-05-14
# o Made into a method of class matrix instead of a stand-alone function.
# 2004-04-18
# o Added support for weighted IWPCA; simply by weighing the weights.
# 2004-03-09
# o Now making use of our own weighted PCA method instead of the
#   multidim::acp() method, which is not maintained anymore.
# 2003-12-29
# o Generalized the method to fit a R-dimensional hyperplane instead of
#   just a line, which is the default fit.
# 2003-12-28
# o Added to the R.basic package.
# o Created by extracted code from RGData.fitMultiIWPCA().
############################################################################
