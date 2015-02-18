#########################################################################/**
# @RdocGeneric plotDensity
# @alias plotDensity.list
# @alias plotDensity.data.frame
# @alias plotDensity.matrix
# @alias plotDensity.numeric
# @alias plotDensity.density
#
# @title "Plots density distributions for a set of vectors"
#
# \description{
#  @get "title".
# }
#
# \usage{
# @usage plotDensity,data.frame
# @usage plotDensity,matrix
# @usage plotDensity,numeric
# @usage plotDensity,list
# }
#
# \arguments{
#  \item{X}{A single of @list of @numeric @vectors or @see "stats::density"
#     objects, a @numeric @matrix, or a @numeric @data.frame.}
#  \item{W}{(optional) weights of similar data types and dimensions as
#     \code{X}.}
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
# \seealso{
#   Internally, @see "stats::density" is used to estimate the
#   empirical density.
# }
#
# @author "HB"
#*/#########################################################################
setMethodS3("plotDensity", "list", function(X, W=NULL, xlim=NULL, ylim=NULL, xlab=NULL, ylab="density (integrates to one)", col=1:length(X), lty=NULL, lwd=NULL, ..., add=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'X':
  nbrOfSamples <- length(X);

  # Argument 'W':
  if (is.numeric(W)) {
    nX <- sapply(X, FUN=length);
    if (any(nX != nX[1L])) {
      throw("If argument 'W' is a numeric vector or matrix, then all vectors of 'X' must of identical lengths, which is not the case.");
    }
    nX <- nX[1L];
    if (is.vector(W)) {
      nW <- length(W);
      if (nW != nX) {
        throw("Length of argument 'W' and the length of the elements of 'X' does not match: ", nW, " != ", nX);
      }
      # Coerce into a list of weights of the same number of elements as 'X'
      W <- rep(list(W), times=nbrOfSamples);
    } else if (is.matrix(W)) {
      nW <- nrow(W);
      if (nW != nX) {
        throw("Number of rows of argument 'W' and the length of the elements of 'X' does not match: ", nW, " != ", nX);
      }
      # Coerce into a list of weights of the same number of elements as 'X'
      Wx <- vector("list", length=ncol(W));
      for (kk in 1:ncol(W)) Wx[[kk]] <- W[,kk,drop=TRUE];
      W <- Wx;
      Wx <- NULL; # Not needed anymore
    }
  } # if (is.numeric(W))

  if (is.list(W)) {
    if (length(W) != nbrOfSamples) {
      throw("The lists of argument 'W' and 'X' do not have the same number of elements: ", length(W), " != ", nbrOfSamples);
    }
    for (kk in 1:nbrOfSamples) {
      w <- W[[kk]];
      nW <- length(w);
      nX <- length(X[[kk]]);
      if (nW != nX) {
        throw(sprintf("Element #%d of arguments 'W' and 'X' are of different lengths: %d != %d", kk, nW, nX));
      }
      if (any(w < 0)) throw("Argument 'W' contains negative weights.");
      w <- nW <- nX <- NULL; # Not needed anymore
    }
  } else if (!is.null(W)) {
    throw("Argument 'W' must be a list, a numeric vector, or a numeric matrix: ", class(W)[1L]);
  }

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
  xlimDef <- c(NA_real_,NA_real_);
  ylimDef <- c(0,NA_real_);
  for(kk in 1:nbrOfSamples) {
    x <- X[[kk]];
    if (inherits(x, "density")) {
      d <- x;
    } else {
      w <- W[[kk]];
      if (is.null(w)) {
        keep <- is.finite(x);
        x <- x[keep];
        keep <- NULL; # Not needed anymore
        suppressWarnings({
          d <- density(x, ...);
        });
        x <- NULL; # Not needed anymore
      } else {
        keep <- is.finite(x) & is.finite(w);
        x <- x[keep];
        w <- w[keep];
        keep <- NULL; # Not needed anymore
        # Standardize to sum(w) == 1
        w <- w / sum(w);
        suppressWarnings({
          d <- density(x, weights=w, ...);
        });
        x <- w <- NULL; # Not needed anymore
      }
    }
    ds[[kk]] <- d;
    xlimDef <- range(c(xlimDef, d$x), na.rm=TRUE);
    ylimDef <- range(c(ylimDef, d$y), na.rm=TRUE);
  } # for (kk ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot the densities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(xlim)) {
    xlim <- xlimDef;
  } else {
    for (kk in 1:2) if (is.na(xlim[kk])) xlim[kk] <- xlimDef[kk]
  }
  if (is.null(ylim)) {
    ylim <- ylimDef;
  } else {
    for (kk in 1:2) if (is.na(ylim[kk])) ylim[kk] <- ylimDef[kk]
  }

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



setMethodS3("plotDensity", "data.frame", function(X, ..., xlab=NULL) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotDensity(as.list(X), ..., xlab=xlab);
})



setMethodS3("plotDensity", "matrix", function(X, ..., xlab=NULL) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotDensity(as.data.frame(X), ..., xlab=xlab);
})


setMethodS3("plotDensity", "numeric", function(X, ..., xlab=NULL) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotDensity(list(X), ..., xlab=xlab);
})


setMethodS3("plotDensity", "density", function(X, ..., xlab=NULL) {
  # Argument 'xlab':
  if (is.null(xlab))
    xlab <- substitute(X);
  plotDensity(list(X), ..., xlab=xlab);
})


##############################################################################
# HISTORY:
# 2014-03-25
# o Now plotDensity() supports weights via argument 'W'.
# o Now plotDensity() also supports density() objects.
# 2006-05-12
# o Created.
##############################################################################
