#########################################################################/**
# @class matrix
# @RdocMethod plotMvsMPairs
#
# @title "Plot log-ratios vs log-ratios for all pairs of columns"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{Nx2K @matrix where N is the number of observations and
#    2K is an even number of channels.}
#  \item{xlab,ylab}{Labels on the x and y axes.}
#  \item{xlim,ylim}{Plot range on the x and y axes.}
#  \item{pch}{Plot symbol used.}
#  \item{...}{Additional arguments accepted by @see "graphics::points".}
#  \item{add}{If @TRUE, data points are plotted in the current plot,
#    otherwise a new plot is created.}
# }
#
# \details{
#  Log-ratio are calculated by over paired columns, e.g. column 1 and 2,
#  column 3 and 4, and so on.
# }
#
# \value{
#   Returns nothing.
# }
#
# @author "HB"
#*/#########################################################################
setMethodS3("plotMvsMPairs", "matrix", function(X, xlab="M", ylab="M", xlim=c(-1,1)*6, ylim=xlim, pch=".", ..., add=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'X':
  if (ncol(X)/2 != round(ncol(X)/2))
    throw("Argument 'X' must have an even number of columns: ", ncol(X));

  if (!add) {
    plot(NA, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim);
  }

  # Do not plot (or generate false) M vs A for  non-positive signals.
  X[X <= 0] <- NA;

  npairs <- ncol(X)/2;
  for (kk in npairs-1) {
    R <- X[,2*kk-1];
    G <- X[,2*kk];
    M1 <- log(R/G, base=2);
    R <- X[,2*(kk+1)-1];
    G <- X[,2*(kk+1)];
    M2 <- log(R/G, base=2);
    points(M1,M2, pch=pch, ...);
  }
})

############################################################################
# HISTORY:
# 2005-09-06
# o Now non-positive signals are excluded.
# 2005-06-11
# o BUG FIX: Used 'rg' instead of 'X'.
# 2005-06-03
# o Created from the normalizeQuantile.matrix.Rex example.
############################################################################
