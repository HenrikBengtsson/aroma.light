#########################################################################/**
# @set "class=list"
# @RdocMethod plotDensity
# @alias plotDensity.data.frame
# @alias plotDensity.matrix
# @alias plotDensity.numeric
#
# @title "Plots density distributions for a set of vectors"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{A single of @list of @numeric @vectors, a @numeric @matrix,
#     or a @numeric @data.frame.}
#  \item{xlim,ylim}{@character @vector of length 2. The x and y limits.}
#  \item{xlab,ylab}{@character string for labels on x and y axis.}
#  \item{col}{The color(s) of the curves.}
#  \item{lty}{The types of curves.}
#  \item{lwd}{The width of curves.}
#  \item{...}{Additional arguments passed to @see "stats::density",
#    @see "graphics::plot", and @see "graphics::lines".}
#  \item{add}{If @TRUE, the curves are plotted in the current plot,
#   otherwise a new is created.}
# }
#
# @author "HB"
#*/#########################################################################
setMethodS3("plotDensity", "list", function(X, xlim=NULL, ylim=NULL, xlab=NULL, ylab="density (integrates to one)", col=1:length(X), lty=NULL, lwd=NULL, ..., add=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'X':
  nbrOfSamples <- length(X);

  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);

  # Argument 'col':
  if (is.null(col)) {
    col <- seq(length=nbrOfSamples);
  } else {
    col <- rep(col, length.out=nbrOfSamples);
  }

  # Argument 'lty':
  if (!is.null(lty))
    lty <- rep(lty, length.out=nbrOfSamples);

  # Argument 'lwd':
  if (!is.null(lwd))
    lwd <- rep(lwd, length.out=nbrOfSamples);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate all densities first and figure out the plot limits.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- list();
  xlimDef <- c(NA,NA);
  ylimDef <- c(0,NA);
  for(kk in 1:nbrOfSamples) {
    x <- X[[kk]];
    x <- x[is.finite(x)];
    suppressWarnings({
      d <- density(x, ...);
    })
    ds[[kk]] <- d;
    xlimDef <- range(c(xlimDef, range(d$x, na.rm=TRUE)), na.rm=TRUE);
    ylimDef <- range(c(ylimDef, range(d$y, na.rm=TRUE)), na.rm=TRUE);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot the densities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(xlim))
    xlim <- xlimDef;
  if (is.null(ylim))
    ylim <- ylimDef;

  if (add == FALSE) {
    suppressWarnings({
      plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...);
    })
  }

  for(kk in 1:nbrOfSamples) {
    suppressWarnings({
      lines(ds[[kk]], col=col[kk], lty=lty[kk], lwd=lwd[kk], ...);
    })
  }

  invisible(ds);
}) # plotDensity()



setMethodS3("plotDensity", "data.frame", function(X, xlab=NULL, ...) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotDensity(as.list(X), xlab=xlab, ...);
})



setMethodS3("plotDensity", "matrix", function(X, xlab=NULL, ...) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotDensity(as.data.frame(X), xlab=xlab, ...);
})


setMethodS3("plotDensity", "numeric", function(X, xlab=NULL, ...) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotDensity(list(X), xlab=xlab, ...);
})


##############################################################################
# HISTORY:
# 2006-05-12
# o Created.
##############################################################################

