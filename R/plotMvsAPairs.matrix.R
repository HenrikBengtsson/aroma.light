#########################################################################/**
# @class matrix
# @RdocMethod plotMvsAPairs
#
# @title "Plot log-ratios/log-intensities for all unique pairs of data vectors"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{NxK @matrix where N is the number of observations and
#    K is the number of channels.}
#  \item{Alab,Mlab}{Labels on the x and y axes.}
#  \item{Alim,Mlim}{Plot range on the A and M axes.}
#  \item{pch}{Plot symbol used.}
#  \item{...}{Additional arguments accepted by @see "graphics::points".}
#  \item{add}{If @TRUE, data points are plotted in the current plot,
#    otherwise a new plot is created.}
# }
#
# \details{
#  Log-ratios and log-intensities are calculated for each neighboring pair
#  of channels (columns) and plotted. Thus, in total there will be K-1
#  data set plotted.
#
#  The colors used for the plotted pairs are 1, 2, and so on. To change
#  the colors, use a different color palette.
# }
#
# \value{
#   Returns nothing.
# }
#
# @author "HB"
#*/#########################################################################
setMethodS3("plotMvsAPairs", "matrix", function(X, Alab="A", Mlab="M", Alim=c(0,16), Mlim=c(-1,1)*diff(Alim), pch=".", ..., add=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'X'
  if (ncol(X) < 2)
    throw("Argument 'X' must have at least two columns: ", ncol(X));

  if (!add) {
    plot(NA, xlab=Alab, ylab=Mlab, xlim=Alim, ylim=Mlim);
  }

  nbrOfChannels <- ncol(X);

  # Do not plot (or generate false) M vs A for  non-positive signals.
  X[X <= 0] <- NA;

  col <- 1;
  for (ii in 1:(nbrOfChannels-1)) {
    Xii <- as.double(X[,ii]);
    for (jj in (ii+1):nbrOfChannels) {
      Xjj <- as.double(X[,jj]);
      M <- log(Xii/Xjj, base=2);
      A <- log(Xii*Xjj, base=2)/2;
      points(A,M, pch=pch, col=col, ...);
      col <- col + 1;
    }
  }
}) # plotMvsAPairs()



############################################################################
# HISTORY:
# 2005-09-06
# o Coercing to doubles to avoid overflow when multiplying to integers.
# o Now non-positive signals are excluded.
# 2005-06-03
# o Now using arguments 'Alab' and 'Mlab' instead of 'xlab' and 'ylab'.
#   Same for 'Alim' and 'Mlim'.
# 2005-04-07
# o Created Rdoc comments.
#############################################################################
