###########################################################################/**
# @set "class=matrix"
# @RdocMethod plotXYCurve
#
# @title "Plot the relationship between two variables as a smooth curve"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{X, Y}{Two @numeric NxK @matrix.}
#   \item{col}{A @vector of colors to be used for each of columns.}
#   \item{lwd}{A @vector of line widths to be used for each of columns.}
#   \item{xlim, ylim}{The x and y plotting limits.}
#   \item{xlab, ylab}{The x and y labels.}
#   \item{...}{Additional arguments passed to @see "plotXYCurve.numeric".}
#   \item{add}{If @TRUE, the graph is added to the current plot, otherwise
#     a new plot is created.}
# }
#
# \value{
#   Returns (invisibly) the curve fit.
# }
#
# \section{Missing values}{
#   Data points (x,y) with non-finite values are excluded.
# }
#
# \seealso{
#   Internally @see "plotXYCurve.numeric" is used.
# }
#
# @author "HB"
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("plotXYCurve", "matrix", function(X, Y, col=1:nrow(X), lwd=NULL, xlim=NULL, ylim=xlim, xlab=NULL, ylab=NULL, ..., add=FALSE) {
  if (identical((X), dim(Y))) {
    throw("Argument 'X' and 'Y' have different dimensions.");
  }

  if (is.null(xlab))
    xlab <- substitute(X);

  if (is.null(ylab))
    ylab <- substitute(Y);

  ncol <- ncol(X);
  if (is.null(col)) {
    col <- seq(length=ncol);
  } else {
    col <- rep(col, length.out=ncol);
  }

  if (!is.null(lwd))
    lwd <- rep(lwd, length.out=ncol);

  if (is.null(xlim)) {
    xlim <- range(X, na.rm=TRUE);
  }
  if (is.null(ylim)) {
    ylim <- range(Y, na.rm=TRUE);
  }

  for (kk in seq(length=ncol)) {
    plotXYCurve(X[,kk], Y[,kk], col=col[kk], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ..., add=add || (kk > 1));
  }

  box();

  invisible();
}) # plotXYCurve.matrix()


############################################################################
# HISTORY:
# 2007-02-04
# o Created.
############################################################################
