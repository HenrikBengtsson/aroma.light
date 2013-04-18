#########################################################################/**
# @set "class=matrix"
# @RdocMethod medianPolish
#
# @title "Median polish"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{N-times-K @matrix}
#  \item{tol}{A @numeric value greater than zero used as a threshold
#     to identify when the algorithm has converged.}
#  \item{maxIter}{Maximum number of iterations.}
#  \item{na.rm}{If @TRUE (@FALSE), @NAs are exclude (not exclude).
#     If @NA, it is assumed that \code{X} contains no @NA values.}
#  \item{.addExtra}{If @TRUE, the name of argument \code{X} is returned
#     and the returned structure is assigned a class.  This will make the
#     result compatible what @see "stats::medpolish" returns.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @list structure with elements:
#    \item{overall}{The fitted constant term.}
#    \item{row}{The fitted row effect.}
#    \item{col}{The fitted column effect.}
#    \item{residuals}{The residuals.}
#    \item{converged}{If @TRUE, the algorithm converged, otherwise not.}
# }
#
# \details{
#   The implementation of this method give identical estimates as
#   @see "stats::medpolish", but is about 3-5 times more efficient when
#   there are no @NA values.
# }
#
# @author "HB"
#
# @examples "../incl/medianPolish.matrix.Rex"
#
# \seealso{
#   @see "stats::medpolish".
# }
#
# @keyword "algebra"
#*/#########################################################################
setMethodS3("medianPolish", "matrix", function(X, tol=0.01, maxIter=10, na.rm=NA, ..., .addExtra=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  .psortKM <- matrixStats:::.psortKM;

  dim <- dim(X);
  nrow <- dim[1];
  ncol <- dim[2];

  if (.addExtra) {
    name <- deparse(substitute(X));
  }

  # Overall effects
  t <- 0;

  # Row effects
  r <- vector("double", nrow);

  # Column effects
  c <- vector("double", ncol);

  hasNa <- (!is.na(na.rm) && any(is.na(X)));
  if (hasNa) {
    oldSum <- 0;
    for (ii in 1:maxIter) {
      # Fit the row effects
      rdelta <- apply(X, MARGIN=1, FUN=median, na.rm=na.rm);
      X <- X - rdelta;
      r <- r + rdelta

      # Fit the overall effects
      delta <- median(c, na.rm=na.rm)
      c <- c - delta
      t <- t + delta

      # Fit the column effects
      cdelta <- apply(X, MARGIN=2, FUN=median, na.rm=na.rm);
      X <- X - matrix(cdelta, nrow=nrow, ncol=ncol, byrow=TRUE)
      c <- c + cdelta

      # Fit the overall effects
      delta <- median(r, na.rm=na.rm)
      r <- r - delta
      t <- t + delta

      # Fit the overall effects
      newSum <- sum(abs(X), na.rm=na.rm);
      converged <- (newSum == 0 || abs(newSum - oldSum) < tol * newSum);
      if (converged)
          break;

      oldSum <- newSum;
    } # for (ii ...)
  } else {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Optimized code for the case where there are no NAs
    #
    # Compared to medpolish(..., na.rm=FALSE), this version is about
    # 3-4 times faster.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    rhalf <- (nrow+1)/2;
    if (nrow %% 2 == 1) {
      # Get x(rhalf), where x(k) is k:th sorted value in x.
      rMedian <- function(x) .psortKM(x, k=rhalf);
    } else {
      # Average x(rhalf) and x(rhalf+1).
      rMedian <- function(x) sum(.psortKM(x, k=rhalf+1L, m=2L))/2;
    }

    chalf <- (ncol+1)/2;
    if (ncol %% 2 == 1) {
      # Get x(chalf), where x(k) is k:th sorted value in x.
      cMedian <- function(x) .psortKM(x, k=chalf);
    } else {
      # Average x(chalf) and x(chalf+1).
      cMedian <- function(x) sum(.psortKM(x, k=chalf+1L, m=2L))/2;
    }

    oldSum <- 0;
    for (ii in 1:maxIter) {
      # Fit the row effects
      rdelta <- apply(X, MARGIN=1, FUN=cMedian);
      X <- X - rdelta;
      r <- r + rdelta;

      # Fit the overall effects
      delta <- cMedian(c);
      c <- c - delta;
      t <- t + delta;

      # Fit the column effects
      cdelta <- apply(X, MARGIN=2, FUN=rMedian);
      X <- X - matrix(cdelta, nrow=nrow, ncol=ncol, byrow=TRUE)
      c <- c + cdelta;

      # Fit the overall effects
      delta <- rMedian(r);
      r <- r - delta;
      t <- t + delta;

      # Fit the overall effects
      newSum <- sum(abs(X), na.rm=FALSE);
      converged <- (newSum == 0 || abs(newSum - oldSum) < tol * newSum);
      if (converged)
          break;

      oldSum <- newSum;
    } # for (ii ...)
  }

  res <- list(overall=t, row=r, col=c, residuals=X, converged=converged);
  if (.addExtra) {
    res$name <- name;
    class(res) <- c("medianPolish", "medpolish");
  }

  res;
}) # medianPolish()


############################################################################
# HISTORY:
# 2012-09-12
# o ROBUSTNESS: Replaced an .Internal(psort(...)) call with a call to
#   matrixStats:::.psortKM() in medianPolish().
# 2012-04-16
# o Added local function psortGet() to medianPolish().
# 2006-05-16
# o Created from stats::medpolish().
############################################################################
