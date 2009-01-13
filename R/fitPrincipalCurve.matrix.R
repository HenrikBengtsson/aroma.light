#########################################################################/**
# @set "class=matrix"
# @RdocMethod fitPrincipalCurve
#
# \encoding{latin1}
#
# @title "Fit a principal curve in K dimensions"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An NxK @matrix (K>=2) where the columns represent the dimension.}
#  \item{...}{Other arguments passed to @see "princurve::principal.curve".}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a principal.curve object (which is a @list).
#   See @see "princurve::principal.curve" for more details.
# }
#
# \section{Missing values}{
#  The estimation of the affine normalization function will only be made
#  based on complete observations, i.e. observations that contains no @NA
#  values in any of the channels.
# }
#
# @author
#
# \references{
#   [1] Hastie, T. and Stuetzle, W, \emph{Principal Curves}, JASA, 1989.
# }
#
# @examples "../incl/fitPrincipalCurve.matrix.Rex"
#
# \seealso{
#   @seemethod "backtransformPrincipalCurve".
#   @see "princurve::principal.curve".
# }
#*/######################################################################### 
setMethodS3("fitPrincipalCurve", "matrix", function(X, ..., verbose=FALSE) {
  require("princurve") || throw("Package not loaded: princurve");

  # The current implementation contains bugs. /HB 2008-05-26
  principal.curve <- principal.curve.hb;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  n <- nrow(X);
  p <- ncol(X);

  # Argument 'verbose':
  if (inherits(verbose, "Verbose")) {
  } else if (is.numeric(verbose)) {
    require("R.utils") || throw("Package not available: R.utils");
    verbose <- Verbose(threshold=verbose);
  } else {
    verbose <- as.logical(verbose);
    if (verbose) {
      require("R.utils") || throw("Package not available: R.utils");
      verbose <- Verbose(threshold=-1);
    }
  }
  if (verbose && inherits(verbose, "Verbose")) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting principal curve");
  verbose && cat(verbose, "Data size: ", n, "x", p);

  verbose && enter(verbose, "Identifying missing values");
  # princurve::principal.curve() does not handle missing values.
  keep <- rep(TRUE, n);
  for (cc in seq(length=p))
    keep <- keep & is.finite(X[,cc]);
  anyMissing <- (!all(keep));
  if (anyMissing)
    X <- X[keep,];
  verbose && exit(verbose);

  verbose && cat(verbose, "Data size after removing non-finite data points: ", nrow(X), "x", p);


  verbose && enter(verbose, "Calling principal.curve()");
  fit <- principal.curve(X, ...);
  verbose && exit(verbose);

  if (anyMissing) {
    values <- matrix(as.double(NA), nrow=n, ncol=p);
    values[keep,] <- fit$s;
    fit$s <- values;
    values <- rep(as.double(NA), times=n);
    for (ff in c("tag", "lambda")) {
      values[keep] <- fit[[ff]];
      fit[[ff]] <- values;
    }
  }

  verbose && exit(verbose);

  class(fit) <- c("PrincipalCurve", class(fit));

  fit;
}) # fitPrincipalCurve()



###########################################################################
# HISTORY:
# 2009-01-12
# o Updated code such that R.utils::Verbose is optional.
# 2008-10-08
# o Removed argument 'fixDimension'.  That constrain is taken care of
#   by backtransformPrincipalCurve().
# o Now the fitted object is of class PrincipalCurve that extends the
#   princurve::principal.curve class.
# 2008-10-07
# o Added Rdoc comments and an example.
# o Removed implementation for data.frame:s.
# 2008-10-03
# o Added argument 'fixDimension'.
# 2008-05-27
# o Added fitPrincipalCurve().
# o Created.
###########################################################################
