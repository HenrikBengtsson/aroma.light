#########################################################################/**
# @class matrix
# @RdocMethod sampleCorrelations
#
# @title "Calculates the correlation for random pairs of observations"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An NxK @matrix where N >= 2 and K >= 2.}
#  \item{MARGIN}{The dimension (1 or 2) in which the observations are.
#    If \code{MARGIN==1} (\code{==2}), each row (column) is an observation.}
#  \item{pairs}{If a Lx2 @matrix, the L index pairs for which the
#    correlations are calculated.
#    If @NULL, pairs of observations are sampled.}
#  \item{npairs}{The number of correlations to calculate.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @double @vector of length \code{npairs}.
# }
#
# @author "HB"
#
# @examples "../incl/sampleCorrelations.matrix.Rex"
#
# \seealso{
#   @see "base::sample".
# }
#
# \references{
#  [1] A. Ploner, L. Miller, P. Hall, J. Bergh & Y. Pawitan.
#      \emph{Correlation test to assess low-level processing of high-density
#      oligonucleotide microarray data}. BMC Bioinformatics, 2005, vol 6.
# }
#
# @keyword utilities
#*/#########################################################################
setMethodS3("sampleCorrelations", "matrix", function(X, MARGIN=1, pairs=NULL, npairs=max(5000, nrow(X)), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  corFast <- function(x, y, ...) {
    ## .Internal() calls are no longer allowed. /HB 2012-04-16
    ## # 3 == "pairwise.complete.obs"
    ## .Internal(cor(x, y, as.integer(3), FALSE));
    cor(x=x, y=y, use="pairwise.complete.obs", method="pearson");
  } # corFast()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'X'
  if (!is.matrix(X))
    throw("Argument 'X' must be a matrix: ", mode(X));
  if (nrow(X) < 2)
    throw("Argument 'X' must have more than two rows.");
  if (ncol(X) < 2)
    throw("Argument 'X' must have more than two columns.");

  # Argument 'MARGIN'
  if (MARGIN < 1 || MARGIN > 2)
    throw("Argument 'MARGIN' is out of range [1,2]: ", MARGIN);

  # Argument 'npairs'
  if (npairs < 1)
    throw("Argument 'npairs' must be equal or greater than one: ", npairs);

  # Get row/column-index pairs to calculate correlations for.
  if (is.null(pairs)) {
    pairs <- sampleTuples(dim(X)[MARGIN], size=npairs, length=2);
  } else {
    npairs <- nrow(pairs);
  }

  # Are 'pairs' and 'npairs' consistent with each other?
  if (nrow(pairs) < npairs) {
    throw("The number of pairs in 'pairs' is smaller than 'npairs': ",
                                            nrow(pairs), " < ", npairs);
  }

  # Pre-create result vector to optimize speed (and memory)
  cors <- rep(as.double(NA), times=npairs);

  if (MARGIN == 1) {
    for (kk in 1:npairs) {
      pair <- pairs[kk,];
      x <- X[pair[1],];
      y <- X[pair[2],];
      cors[kk] <- corFast(x,y);
    }
  } else {
    for (kk in 1:npairs) {
      pair <- pairs[kk,];
      x <- X[,pair[1]];
      y <- X[,pair[2]];
      cors[kk] <- corFast(x,y);
    }
  }

  cors;
}) # sampleCorrelations()


############################################################################
# HISTORY:
# 2012-04-16
# o sampleCorrelations() no longer utilizes .Internal() calls.
# o Added internal corFast() to sampleCorrelations().
# 2011-04-12
# o Now using NAs of the correct storage type.
# 2005-07-25
# o Added Rdoc comments with a simple example.
# 2005-04-07
# o Created.
############################################################################
