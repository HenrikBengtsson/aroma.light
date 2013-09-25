#########################################################################/**
# @set "class=matrix"
# @RdocMethod fitPrincipalCurve
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
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a principal.curve object (which is a @list).
#   See @see "princurve::principal.curve" for more details.
# }
#
# \section{Missing values}{
#  The estimation of the normalization function will only be made
#  based on complete observations, i.e. observations that contains no @NA
#  values in any of the channels.
# }
#
# @author "HB"
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

  # princurve v1.1-9 and before contains bugs. /HB 2008-05-26
  ver <- packageDescription("princurve")$Version;
  if (compareVersion(ver, "1.1-10") < 0) {
    throw("princurve v1.1-10 or newer is required: ", ver);
  }

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
    cat <- R.utils::cat;
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting principal curve");
  verbose && cat(verbose, "Data size: ", n, "x", p);

  verbose && enter(verbose, "Identifying missing values");
  # princurve::principal.curve() does not handle missing values.
  keep <- rep(TRUE, times=n);
  for (cc in seq(length=p)) {
    keep <- keep & is.finite(X[,cc]);
  }
  anyMissing <- (!all(keep));
  if (anyMissing) {
    X <- X[keep,, drop=FALSE];
  }
  verbose && exit(verbose);

  verbose && cat(verbose, "Data size after removing non-finite data points: ", nrow(X), "x", p);


  verbose && enter(verbose, "Calling principal.curve()");
  trace <- as.logical(verbose);
  t <- system.time({
    fit <- principal.curve(X, ..., trace=trace);
  });
  attr(fit, "processingTime") <- t;
  verbose && printf(verbose, "Converged: %s\n", fit$converged);
  verbose && printf(verbose, "Number of iterations: %d\n", fit$nbrOfIterations);
  verbose && printf(verbose, "Processing time/iteration: %.1fs (%.1fs/iteration)\n", t[3], t[3]/fit$nbrOfIterations);
  verbose && exit(verbose);

  # Expand, iff missing values were dropped
  if (anyMissing) {
    values <- matrix(as.double(NA), nrow=n, ncol=p);
    values[keep,] <- fit$s;
    dimnames(values) <- dimnames(fit$s);
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
# 2013-04-18
# o BUG FIX: fitPrincipalCurve() would not preserve dimension names
#   if data contain missing values.
# 2011-04-12
# o CLEANUP: Removed internal patch of principal.curve().  If an older
#   version than princurve v1.1-10 is used, an informative error is
#   thrown requesting an update.  The internal patch is part of the
#   offical princurve v1.1-10 release (since 2009-10-04).
# 2009-11-01
# o Now fitPrincipalCurve() bug-fixed princurve v1.1-10.  If earlier
#   version are available, it used the internal patch.
# 2009-07-15
# o Added attribute 'processingTime' to the fit object returned by
#   fitPrincipalCurve().
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
