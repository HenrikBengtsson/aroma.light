#########################################################################/**
# @class matrix
# @RdocMethod plotMvsA
#
# @title "Plot log-ratios vs log-intensities"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{Nx2 @matrix with two channels and N observations.}
#  \item{Alab,Mlab}{Labels on the x and y axes.}
#  \item{Alim,Mlim}{Plot range on the A and M axes.}
#  \item{aspectRatio}{Aspect ratio between \code{Mlim} and \code{Alim}.}
#  \item{pch}{Plot symbol used.}
#  \item{...}{Additional arguments accepted by @see "graphics::points".}
#  \item{add}{If @TRUE, data points are plotted in the current plot,
#    otherwise a new plot is created.}
# }
#
# \details{
#  Red channel is assumed to be in column one and green in column two.
#  Log-ratio are calculated taking channel one over channel two.
# }
#
# \value{
#   Returns nothing.
# }
#
# @author "HB"
#*/#########################################################################
setMethodS3("plotMvsA", "matrix", function(X, Alab="A", Mlab="M", Alim=c(0,16), Mlim=c(-1,1)*diff(Alim)*aspectRatio, aspectRatio=1, pch=".", ..., add=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'X':
  if (ncol(X) != 2) {
    throw("Argument 'X' must have exactly two columns: ", ncol(X));
  }

  if (!add) {
    plot(NA, xlab=Alab, ylab=Mlab, xlim=Alim, ylim=Mlim);
  }

  # Do not plot (or generate false) M vs A for  non-positive signals.
  X[X <= 0] <- NA;

  R <- as.double(X[,1]);
  G <- as.double(X[,2]);
  M <- log(R/G, base=2);
  A <- log(R*G, base=2)/2;
  points(A,M, pch=pch, ...);
})


############################################################################
# HISTORY:
# 2011-06-26
# o Added argument 'aspectRatio' to plotMvsA().  It can be used to adjust
#   the range of the 'Mlim' argument relative to the 'Alim' argument.
# 2005-09-06
# o Coercing to doubles to avoid overflow when multiplying to integers.
# o Now non-positive signals are excluded.
# 2005-06-11
# o BUG FIX: Used 'rg' instead of 'X' in R <- rg[,1] and G <- rg[,2].
# 2005-06-03
# o Created from the normalizeQuantile.matrix.Rex example.
############################################################################
