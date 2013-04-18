#########################################################################/**
# @RdocDefault sampleTuples
#
# @title "Sample tuples of elements from a set"
#
# \description{
#   @get "title".
#   The elements within a sampled tuple are unique, i.e. no two elements
#   are the same.
# }
#
# @synopsis
#
# \arguments{
#  \item{x}{A set of elements to sample from.}
#  \item{size}{The number of tuples to sample.}
#  \item{length}{The length of each tuple.}
#  \item{...}{Additional arguments passed to @see "base::sample".}
# }
#
# \value{
#   Returns a NxK @matrix where N = \code{size} and K = \code{length}.
# }
#
# @author "HB"
#
# @examples "../incl/sampleTuples.Rex"
#
# \seealso{
#   @see "base::sample".
# }
#
# @keyword utilities
#*/#########################################################################
setMethodS3("sampleTuples", "default", function(x, size, length, ...) {
  # Argument 'x':
  if (length(x) < 1)
    throw("Argument 'x' must be a vector of length one or greater.");

  # Argument 'size':
  if (size < 0)
    throw("Argument 'size' must be a non-negative integer: ", size);

  # Argument 'length':
  if (length < 1)
    throw("Argument 'length' must be one or greater: ", length);

  # Sample tuples
  naValue <- NA;
  storage.mode(naValue) <- storage.mode(x);
  tuples <- matrix(naValue, nrow=size, ncol=length);
  for (kk in seq(length=size)) {
    tuples[kk,] <- sample(x, size=length, ...);
  }

  tuples;
})


############################################################################
# HISTORY:
# 2011-04-12
# o Now using NAs of the correct storage type.
# 2005-07-25
# o Created generic sampleTuples().
# 2005-04-07
# o Created.
############################################################################
