###########################################################################/**
# @set "class=numeric"
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
#   \item{x, y}{Two @numeric @vector of length N.}
#   \item{lwd}{The width of the curve.}
#   \item{col}{The color of the curve.}
#   \item{dlwd}{The width of the density curves.}
#   \item{dcol}{The fill color of the interior of the density curves.}
#   \item{xlim, ylim}{The x and y plotting limits.}
#   \item{xlab, ylab}{The x and y labels.}
#   \item{curveFit}{The @function used to fit the curve.  The two first
#     arguments of the function must take \code{x} and \code{y}, and the
#     function must return a @list with fitted elements \code{x} and
#     \code{y}.}
#   \item{...}{Additional arguments passed to @see "graphics::lines"
#     used to draw the curve.}
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
#   @see "plotXYCurve.matrix".
# }
#
# @author "HB"
#
# @keyword "nonparametric"
# @keyword "multivariate"
# @keyword "robust"
#*/###########################################################################
setMethodS3("plotXYCurve", "numeric", function(x, y, lwd=2, col=1, dlwd=1, dcol=NA, xlim=NULL, ylim=xlim, xlab=NULL, ylab=NULL, curveFit=smooth.spline, ..., add=FALSE) {
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

    abline(a=0, b=1, col="gray", lty=2);
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
  rx[1] <- max(rx[1], usr[1]);
  rx[2] <- min(rx[2], usr[2]);

  ry <- range(y);
  ry[1] <- max(ry[1], usr[3]);
  ry[2] <- min(ry[2], usr[4]);

  # Estimate density of x
  d <- density(x, from=rx[1], to=rx[2]);
  n <- length(d$y);
  d$y[c(1,n)] <- 0;
  d$y <- d$y / max(d$y, na.rm=TRUE);
  dx <- d;
  d$y <- 1/10*(usr[4]-usr[3])*d$y;
  d$y <- usr[4]+d$y;
  polygon(d, col=dcol, border=col, lwd=dlwd, xpd=TRUE);

  # Estimate density of y
  d <- density(y, from=ry[1], to=ry[2]);
  n <- length(d$y);
  d$y[c(1,n)] <- 0;
  d$y <- d$y / max(d$y, na.rm=TRUE);
  dy <- d;
  t <- d$x; d$x <- d$y; d$y <- t;
  d$x <- usr[1]-1/10*(usr[2]-usr[1])*d$x;
  polygon(d, col=dcol, border=col, lwd=dlwd, xpd=TRUE);

  d <- dx;
  t <- d$x; d$x <- d$y; d$y <- t;
  d$x <- usr[1]-1/10*(usr[2]-usr[1])*d$x;
  lines(d, col="black", lwd=0.618*dlwd, lty=2, xpd=TRUE);

  if (!add) {
    axis(side=1);
    axis(side=4);
    box();
  }

  invisible(fit);
}) # plotXYCurve.numeric()


############################################################################
# HISTORY:
# 2008-04-14
# o Replaced a R.utils::doCall() with a "cleanup" do.call().
# 2007-02-04
# o Created.
############################################################################
