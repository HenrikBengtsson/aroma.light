#########################################################################/**
# @set "class=matrix"
# @RdocMethod wpca
#
# @title "Light-weight Weighted Principal Component Analysis"
#
# @synopsis
#
# \description{
#   Calculates the (weighted) principal components of a matrix, that is,
#   finds a new coordinate system (not unique) for representing the given
#   multivariate data such that
#    i) all dimensions are orthogonal to each other, and
#   ii) all dimensions have maximal variances.
# }
#
# \arguments{
#   \item{x}{An NxK @matrix.}
#   \item{w}{An N @vector of weights for each row (observation) in
#     the data matrix. If @NULL, all observations get the same weight,
#     that is, standard PCA is used.}
#   \item{center}{If @TRUE, the (weighted) sample mean column @vector is
#     subtracted from each column in \code{mat}, first.
#     If data is not centered, the effect will be that a linear subspace
#     that goes through the origin is fitted.}
#   \item{scale}{If @TRUE, each column in \code{mat} is
#     divided by its (weighted) root-mean-square of the
#     centered column, first.}
#   \item{method}{If \code{"dgesdd"} LAPACK's divide-and-conquer
#     based SVD routine is used (faster [1]), if \code{"dgesvd"}, LAPACK's
#     QR-decomposition-based routine is used, and if \code{"dsvdc"},
#     LINPACK's DSVDC(?) routine is used. The latter is just for
#     pure backward compatibility with R v1.7.0.
#   }
#   \item{swapDirections}{If @TRUE, the signs of eigenvectors
#     that have more negative than positive components are inverted.
#     The signs of corresponding principal components are also inverted.
#     This is only of interest when for instance visualizing or comparing
#     with other PCA estimates from other methods, because the
#     PCA (SVD) decompostion of a matrix is not unique.
#   }
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list with elements:
#   \item{pc}{An NxK @matrix where the column @vectors are the
#             principal components (a.k.a. loading vectors,
#             spectral loadings or factors etc).}
#   \item{d}{An K @vector containing the eigenvalues of the
#             principal components.}
#   \item{vt}{An KxK @matrix containing the eigenvector of the
#             principal components.}
#   \item{xMean}{The center coordinate.}
#
#   It holds that \code{x == t(t(fit$pc \%*\% fit$vt) + fit$xMean)}.
# }
#
# \section{Method}{
#   A singular value decomposition (SVD) is carried out.
#   Let X=\code{mat}, then the SVD of the matrix is \eqn{X = U D V'}, where
#   \eqn{U} and \eqn{V} are othogonal, and \eqn{D} is a diagonal matrix
#   with singular values. The principal returned by this method are \eqn{U D}.
#
#   Internally \code{La.svd()} (or \code{svd()}) of the \pkg{base}
#   package is used.
#   For a popular and well written introduction to SVD see for instance [2].
# }
#
# \examples{
#   @include "../incl/wpca.matrix.Rex"
#
#   if (dev.cur() > 1) dev.off()
#
#   @include "../incl/wpca2.matrix.Rex"
# }
#
# @author
#
# \references{
#   [1] J. Demmel and  J. Dongarra, \emph{DOE2000 Progress Report}, 2004.
#       \url{http://www.cs.berkeley.edu/~demmel/DOE2000/Report0100.html} \cr
#   [2] Todd Will,
#       \emph{Introduction to the Singular Value Decomposition},
#       UW-La Crosse, 2004.
#       \url{http://www.uwlax.edu/faculty/will/svd/} \cr
# }
#
# \seealso{
#   For a iterative re-weighted PCA method, see @seemethod "iwpca".
#   For Singular Value Decomposition, see @see "base::svd".
#   For other implementations of Principal Component Analysis functions see
#   (if they are installed):
#   @see "stats::prcomp" in package \pkg{stats} and
#   \code{pca()} in package \pkg{pcurve}.
# }
#
# @keyword "algebra"
#*/#########################################################################
setMethodS3("wpca", "matrix", function(x, w=NULL, center=TRUE, scale=FALSE, method=c("dgesdd", "dgesvd", "dsvdc"), swapDirections=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'x'
  x <- as.matrix(x);

  # The dimensions of 'x'
  N <- nrow(x);
  K <- ncol(x);

  # Standardizes the weights to [0,1] such that they sum to 1.
  if (!is.null(w)) {
    w <- rep(w, length.out=N);
    w <- w/sum(w);
    if (any(is.na(w)))
      stop("Argument 'w' has missing values.");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Weighted or non-weighted centering and rescaling of the data
  #
  # Note: The following split of (center == TRUE) and (center == FALSE)
  # is to minimize memory usage. In other words, the codes is longer,
  # but more memory efficient.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (center || scale) {
    if (is.null(w)) {
      # Calculates the standard column means
      xMean <- colMeans(x, na.rm=TRUE);           # a K vector
    } else {
      # Calculates the weighted column means (recall that sum(w) == 1)
      xMean <- as.vector(w %*% x);                # a K vector
    }

    if (center) {
      # Centers the data directly by subtracting the column means
      for (kk in 1:ncol(x))
        x[,kk] <- x[,kk] - xMean[kk];
    } else {
      # ...or calculates the centered data for rescaling
      xc <- x;
      for (kk in 1:ncol(x))
        xc[,kk] <- x[,kk] - xMean[kk];
    }

    if (scale) {
      if (is.null(w)) {
        # Non-weighted root-mean-squares
        rms <- function(v) {         # v - column vector of length N
          v <- v[!is.na(v)];
          sqrt(sum(v^2)/max(1, length(v)-1));
        }
      } else {
        # Weighted root-mean-squares
        rms <- function(v) {         # v - column vector of length N
          ok <- !is.na(v);
          v <- w[ok]*v[ok];
          sqrt(sum(v^2)/max(1, length(v)-1));
        }
      }

      if (center) {
        xRMS <- apply(x, MARGIN=2, FUN=rms);
      } else {
        xRMS <- apply(xc, MARGIN=2, FUN=rms);
        # Not needed anymore
        xc <- NULL;
      }

      for (kk in 1:ncol(x))
        x[,kk] <- x[,kk] / xRMS[kk];

      # Not needed anymore
      xRMS <- rms <- NULL;
    }
  } else {
    xMean <- rep(0, length=K);
  }

  # Weight the observations?
  if (!is.null(w)) {
    x <- sqrt(w)*x;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Singular Value Decomposition, i.e. X = U D V'
  #
  # "We compare DGESVD, the original QR-based routines from
  #  LAPACK 1.0, with DGESDD, the new divide-and-conquer based
  #  routine from LAPACK 3.0. The table below shows the speeds
  #  on several machines. The new routine is 5.7 to 16.8 times
  #  faster than the old routine. Part of the speed results
  #  from doing many fewer operations, and the rest comes from
  #  doing them faster (a high Mflop rate)." [1]
  # [1] J. Demmel and  J. Dongarra, DOE2000 Progress Report,
  #     http://www.cs.berkeley.edu/~demmel/DOE2000/Report0100.html
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (method == "dgesdd" || method == "dgesvd") {
    duvt <- La.svd(x);
  } else if (method == "dsvdc") {
    # For backward compatibility with R < 1.7.0. See ?svd for more info.
    duvt <- svd(x, LINPACK=TRUE);
    duvt$vt <- t(duvt$v);
  } else {
    stop(sprintf("Unknown LAPACK or LINPACK routine to solve SVD: %s", method));
  }

  # 'duvt' is a list with the follwing components:
  #
  # u  - a NxK matrix whose columns contain the left singular
  #      vectors (eigenvectors) of 'x'.
  #      It holds that t(u) %*% u == I

  # d  - a K vector containing the singular value of
  #      each principal component on its diagonal.
  #      It holds that d[1] >= d[2] >= ... d[K] >= 0.

  # vt - a KxK transposed matrix whose columns contain the right
  #      singular vectors (eigenvector) of 'x'.
  #      It holds that t(v) %*% v == I

  # Not need anymore, in case a local copy has been created!
  x <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 4. The PCA principal components
  #
  # (a.k.a. loading vectors, spectral loadings or factors).
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  d <- duvt$d;
  vt <- duvt$vt;
  pc <- duvt$u;
  # Not need anymore
  duvt <- NULL;

  # Note: D == diag(duvt$d) is memory expensive since the dimensions of D
  # is the same as the dimensions of 'x'. Thus, it unwise to do:
  # pc <- duvt$u %*% diag(duvt$d);
  for (kk in seq(length=N))
    pc[kk,] <- pc[kk,] * d;

  if (!is.null(w)) {
    # Rescale the principal components
    pc <- pc / sqrt(w);
    # Not need anymore
    w <- NULL;
  }

  if (swapDirections) {
    swap <- apply(vt, MARGIN=1, FUN=function(z) sum(sign(z)) < 0);

    # Which eigenvectors should swap signs?
    swap <- which(swap);
    for (kk in swap) {
      vt[kk,] <- -vt[kk,];
      pc[,kk] <- -pc[,kk];
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 4. Return the parameter estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list(pc=pc, d=d, vt=vt, xMean=xMean);
  class(res) <- "WPCAFit";
  res;
}) # wpca()


############################################################################
# HISTORY:
# 2006-06-26
# o Function would not work in R v2.4.0 devel, because argument 'method' was
#   removed from La.svd().
# 2006-04-25
# o Updated the URL to Todd Will's webpage.
# 2005-02-20
# o Added '...' to please R CMD check.
# 2005-02-20
# o Now using setMethodS3() and added '...' to please R CMD check.
# 2005-01-24
# o Added a check for missing values of argument 'w'.
# 2004-05-14
# o Made into a method of class matrix instead of a stand-alone function.
# 2003-03-09
# o Created! Verified that it gives similar results as acp().
############################################################################
