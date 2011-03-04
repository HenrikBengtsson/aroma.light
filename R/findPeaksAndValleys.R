###########################################################################/**
# @set "class=density"
# @RdocMethod findPeaksAndValleys
# @alias findPeaksAndValleys
# @alias findPeaksAndValleys.numeric
#
# @title "Finds extreme points in the empirical density estimated from data"
#
# \description{
#   @get "title".
# }
# 
# \usage{
#  \method{findPeaksAndValleys}{density}(x, tol=0, ...)
#  \method{findPeaksAndValleys}{numeric}(x, ..., tol=0, na.rm=TRUE)
# }
#
# \arguments{
#  \item{x}{A @numeric @vector containing data points or 
#     a @see "stats::density" object.}
#  \item{...}{Arguments passed to @see "stats::density".
#     Ignored if \code{x} is a @see "stats::density" object.}
#  \item{tol}{A non-negative @numeric threshold specifying the minimum
#    density at the extreme point in order to accept it.}
#  \item{na.rm}{If @TRUE, missing values are dropped, otherwise not.}
# }
#
# \value{
#   Returns a @data.frame (of class 'PeaksAndValleys') containing 
#   of "peaks" and "valleys" filtered by \code{tol}.
# }
#
# @examples "..\incl\findPeaksAndValleys.Rex"
#
# @author
#
# @keyword internal
#*/########################################################################### 
setMethodS3("findPeaksAndValleys", "density", function(x, tol=0, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'x':
  d <- x;

  # Argument 'tol':
  tol <- as.double(tol);
  stopifnot(length(tol) == 1);
  stopifnot(tol >= 0);

  delta <- diff(d$y);
  n <- length(delta);

  isPeak <- (delta[-n] > 0 & delta[-1] < 0);
  isValley <- (delta[-n] < 0 & delta[-1] > 0);
  isPeakOrValley <- (isPeak | isValley);

  idxs <- which(isPeakOrValley);
  types <- c("valley", "peak")[isPeak[idxs]+1];
  names(idxs) <- types;

  x <- d$x[idxs];
  y <- d$y[idxs];
  res <- data.frame(type=types, x=x, density=y);
  class(res) <- c("PeaksAndValleys", class(res));

  # Filter valleys by density?
  if (tol > 0) {
    res <- subset(res, density >= tol);
  }

  res;
}) # findPeaksAndValleys()


setMethodS3("findPeaksAndValleys", "numeric", function(x, ..., tol=0, na.rm=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'na.rm':
  na.rm <- as.logical(na.rm);
  stopifnot(length(na.rm) == 1);

  # Argument 'tol':
  tol <- as.double(tol);
  stopifnot(length(tol) == 1);
  stopifnot(tol >= 0);

  d <- density(x, na.rm=na.rm, ...);
  findPeaksAndValleys(d, tol=tol);
}) # findPeaksAndValleys()




############################################################################
# HISTORY:
# 2011-03-03 [HB]
# o BUG FIX: findPeaksAndValleys(x, to) were 'x' is numeric would use
#   partial match and interpret 'to' as argument 'tol' and not part of
#   '...' passed to density().  This problem was solved by placing '...'
#   before argument 'tol'.  Thanks Oscar Rueda at the Cancer Reasearch UK
#   for reporting and identify this bug.
# 2010-10-08 [HB]
# o Now findPeaksAndValleys() returns a object of class PeaksAndValleys,
#   which extends data.frame.
# 2010-10-06 [HB]
# o Added findPeaksAndValleys() for the 'density' class, which then
#   findPeaksAndValleys() for 'numeric' utilizes.
# 2010-04-04 [HB]
# o Made findPeaksAndValleys() an internal function in Rd.
# o Updated could to validate arguments with using R.utils::Arguments.
# o Corrected a non-defined Rdoc tag.
# 2009-11-03 [HB]
# o Added Rdoc comments with an example().
# 2009-03-06 [HB]
# o Created for doing quick naive genotyping of some TCGA normal samples in
#   order to highlight the centers of the clouds in a tumor-normal fracB
#   scatter plots.
############################################################################
