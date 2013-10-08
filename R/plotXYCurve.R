###########################################################################/**
# @RdocGeneric plotXYCurve
# @alias plotXYCurve.numeric
# @alias plotXYCurve.matrix
#
# @title "Plot the relationship between two variables as a smooth curve"
#
# \usage{
# @usage plotXYCurve,numeric
# @usage plotXYCurve,matrix
# }
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{x, y, X, Y}{Two @numeric @vectors of length N for one curve (K=1),
#     or two @numeric NxK @matrix:es for K curves.}
#   \item{col}{The color of each curve.
#     Either a scalar specifying the same value of all curves,
#     or a @vector of K curve-specific values.}
#   \item{lwd}{The line width of each curve.
#     Either a scalar specifying the same value of all curves,
#     or a @vector of K curve-specific values.}
#   \item{dlwd}{The width of each density curve.}
#   \item{dcol}{The fill color of the interior of each density curve.}
#   \item{xlim, ylim}{The x and y plotting limits.}
#   \item{xlab, ylab}{The x and y labels.}
#   \item{curveFit}{The @function used to fit each curve.  The two first
#     arguments of the function must take \code{x} and \code{y}, and the
#     function must return a @list with fitted elements \code{x} and
#     \code{y}.}
#   \item{...}{Additional arguments passed to @see "graphics::lines"
#     used to draw each curve.}
#   \item{add}{If @TRUE, the graph is added to the current plot, otherwise
#     a new plot is created.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \section{Missing values}{
#   Data points (x,y) with non-finite values are excluded.
# }
#
# @author "HB"
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("plotXYCurve", "numeric", function(x, y, col=1L, lwd=2, dlwd=1, dcol=NA, xlim=NULL, ylim=xlim, xlab=NULL, ylab=NULL, curveFit=smooth.spline, ..., add=FALSE) {
  if (is.null(xlab))
    xlab <- substitute(X);

  if (is.null(ylab))
    ylab <- substitute(Y);

  # Exclude non-finite data points
  ok <- (is.finite(x) & is.finite(y));
  x <- x[ok];
  y <- y[ok];

  # Create empty plot?
  if (!add) {
    par(mar=c(5,4,4,5)+0.1);
    suppressWarnings({
      plot(NA, xlim=xlim, ylim=ylim, xlab="", ylab="", ..., axes=FALSE);
    })

    cex <- par("cex.lab")*par("cex");
    mtext(xlab, side=1, line=3, cex=cex, col=par("col.lab"), font=par("font.lab"));
    mtext(ylab, side=4, line=3, cex=cex, col=par("col.lab"), font=par("font.lab"));

    abline(a=0, b=1, col="gray", lty=2L);
  }

  # Estimate and draw smooth function
  suppressWarnings({
    args <- list(x=x, y=y, ...);
    keep <- intersect(names(args), names(formals(smooth.spline)));
    args <- args[keep];
    fit <- do.call(smooth.spline, args=args);
    lines(fit, col=col, lwd=lwd, ...);
  })

  usr <- par("usr");

  # Limit the density plot to the plot region and data range
  rx <- range(x);
  rx[1L] <- max(rx[1], usr[1L]);
  rx[2L] <- min(rx[2], usr[2L]);

  ry <- range(y);
  ry[1L] <- max(ry[1], usr[3L]);
  ry[2L] <- min(ry[2], usr[4L]);

  # Estimate density of x
  d <- density(x, from=rx[1L], to=rx[2L]);
  n <- length(d$y);
  d$y[c(1L,n)] <- 0;
  d$y <- d$y / max(d$y, na.rm=TRUE);
  dx <- d;
  d$y <- 1/10*(usr[4L]-usr[3L])*d$y;
  d$y <- usr[4L]+d$y;
  polygon(d, col=dcol, border=col, lwd=dlwd, xpd=TRUE);

  # Estimate density of y
  d <- density(y, from=ry[1L], to=ry[2L]);
  n <- length(d$y);
  d$y[c(1L,n)] <- 0;
  d$y <- d$y / max(d$y, na.rm=TRUE);
  dy <- d;
  t <- d$x; d$x <- d$y; d$y <- t;
  d$x <- usr[1L]-1/10*(usr[2L]-usr[1L])*d$x;
  polygon(d, col=dcol, border=col, lwd=dlwd, xpd=TRUE);

  d <- dx;
  t <- d$x; d$x <- d$y; d$y <- t;
  d$x <- usr[1L]-1/10*(usr[2L]-usr[1L])*d$x;
  lines(d, col="black", lwd=0.618*dlwd, lty=2, xpd=TRUE);

  if (!add) {
    axis(side=1L);
    axis(side=4L);
    box();
  }

  invisible(fit);
}) # plotXYCurve.numeric()


setMethodS3("plotXYCurve", "matrix", function(X, Y, col=seq(length=nrow(X)), lwd=2, dlwd=1, dcol=NA, xlim=NULL, ylim=xlim, xlab=NULL, ylab=NULL, curveFit=smooth.spline, ..., add=FALSE) {
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
    plotXYCurve(X[,kk], Y[,kk], col=col[kk], lwd=lwd, dlwd=dlwd, dcol=dcol, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, curveFit=curveFit, ..., add=add || (kk > 1L));
  }

  box();

  invisible();
}) # plotXYCurve.matrix()


############################################################################
# HISTORY:
# 2013-10-08
# o BUG FIX: Argument 'lwd' of plotXYCurve(X, ...) was ignored if 'X'
#   was a matrix.
# o DOCUMENTATION: Now there is one combined plotXYCurve() help pages
#   for all data types.
# 2008-04-14
# o Replaced a R.utils::doCall() with a "cleanup" do.call().
# 2007-02-04
# o Created.
############################################################################
