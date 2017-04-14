############################################################################/**
# @RdocDefault robustSmoothSpline
#
# @title "Robust fit of a Smoothing Spline"
#
# @synopsis
#
# \description{
#   Fits a smoothing spline robustly using the \eqn{L_1} norm. Currently, the
#   algorithm is an \emph{iterative reweighted smooth spline} algorithm which
#   calls \code{smooth.spline(x,y,w,...)} at each iteration with the weights
#   \code{w} equal to the inverse of the absolute value of the residuals for
#   the last iteration step.
# }
#
# \arguments{
#   \item{x}{a @vector giving the values of the predictor variable, or a
#            @list or a two-column @matrix specifying \code{x} and \code{y}.
#            If \code{x} is of class \code{smooth.spline} then \code{x$x} is used
#            as the \code{x} values and \code{x$yin} are used as the \code{y}
#            values.}
#   \item{y}{responses. If \code{y} is missing, the responses are assumed to be
#            specified by \code{x}.}
#   \item{w}{a @vector of weights the same length as \code{x} giving the weights
#            to use for each element of \code{x}. Default value is equal weight
#            to all values.}
#   \item{...}{Other arguments passed to @see "stats::smooth.spline".}
#   \item{minIter}{the minimum number of iterations used to fit the smoothing
#            spline robustly. Default value is 3.}
#   \item{maxIter}{the maximum number of iterations used to fit the smoothing
#            spline robustly. Default value is 25.}
#   \item{sdCriteria}{Convergence criteria, which the difference between the
#            standard deviation of the residuals between two consecutive iteration
#            steps. Default value is 2e-4.}
#   \item{reps}{Small positive number added to residuals to avoid division by
#            zero when calculating new weights for next iteration.}
#   \item{tol}{Passed to @see "stats::smooth.spline" (R >= 2.14.0).}
#   \item{plotCurves}{If @TRUE, the fitted splines are added to the current
#         plot, otherwise not.}
# }
#
# \value{
#   Returns an object of class \code{smooth.spline}.
# }
#
# @examples "../incl/robustSmoothSpline.Rex"
#
# \seealso{
#   This implementation of this function was adopted from
#   @see "stats::smooth.spline" of the \pkg{stats} package.
#   Because of this, this function is also licensed under GPL v2.
# }
#
# @author
#
# @keyword "smooth"
# @keyword "robust"
#*/############################################################################
setMethodS3("robustSmoothSpline", "default", function(x, y=NULL, w=NULL, ..., minIter=3, maxIter=max(minIter, 50), sdCriteria=2e-4, reps=1e-15, tol=1e-6*IQR(x), plotCurves=FALSE) {
  requireNamespace("stats") || throw("Package not loaded: stats");  # smooth.spline()

  stats.smooth.spline <- smooth.spline;
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Verify arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'w'
  if (is.numeric(w)) {
    w <- as.double(w);
    if (anyMissing(w)) {
      stop("Weights with value NA are not allowed.");
    }
    if (any(w < 0 | w > 1)) {
      stop("Weights out of range [0,1]: ", paste(w[w < 0.0 | w > 1.0], collapse=", "));
    }
  } else if (!is.null(w)) {
    stop("Argument 'w' is of an unsupported datatype/class: ", class(weights)[1]);
  }

  # Argument: 'reps'
  if (!is.numeric(reps) || length(reps) != 1 || reps <= 0)
    throw("Argument 'reps' must be a single postive number.");

  # smooth.spline() next will only operate on unique x-values. For this reason,
  # we have to remove corresponding weights too. There is a small problem here;
  # if different weights are used for data points (x,y,w) with same x-value, which
  # data points (x,y,w) should be used? Here we use the first one only. /HB 2005-01-24
  uIdxs <- .whichUnique(x, tol = tol);
  nu <- length(uIdxs);
  w0 <- w[uIdxs];

  # WORKAROUND
  # We need to make sure that 'g$x == x' below. /HB 2011-10-10
  x <- x[uIdxs];
  y <- y[uIdxs];
  w <- w[uIdxs];
  uIdxs <- seq_along(x);

  if (inherits(x, "smooth.spline")) {
    g <- x;
  } else if (missing(w) || is.null(w)) {
    x <- as.vector(x);
    y <- as.vector(y);
    g <- stats.smooth.spline(x, y, ..., tol=tol);

    # Sanity check /HB 2011-10-10
    stopifnot(length(g$x) == nu);

    # Not needed anymore
    x <- y <- NULL;
  } else {
    x <- as.vector(x);
    y <- as.vector(y);
    w <- as.vector(w);
    g <- stats.smooth.spline(x, y, w=w, ..., tol=tol);

    # Sanity check /HB 2011-10-10
    stopifnot(length(g$x) == nu);

    # Not needed anymore
    x <- y <- w <- NULL;
  }

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Step 0. Initiation
  #
  # This will generate an object of class smooth.spline
  # containing the fields
  #  x   : the distinct `x' values in increasing order.
  #  y   : the fitted values corresponding to `x'.
  #  yin : the y values used at the unique `y' values.
  # From these the residuals can be calculated as
  #  r <- yin - y
  # The important is that we use these (x,yin) as our
  # (x,y) in the rest of the algorithm.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  sdR0 <- as.double(NA);
  col <- 0;
  ready <- FALSE;
  iter <- 0;
  while (!ready & iter < maxIter) {
    iter <- iter + 1;
    # Calculate the residuals and the weights
    r <- (g$yin-g$y);
    w <- 1/(abs(r)+reps); # Add a small constant for stability.

    # If the user specified weights initially, the weights
    # calculated from the inverse of the residuals are themselve
    # weighted by the user initial weights.
    if (!is.null(w0)) {
      w <- w0*w;
    }

    sdR <- sd(r);
    # Not needed anymore
    r <- NULL;

    if (iter > minIter) {
      if (!is.na(sdR0)) {
        dSd <- abs(sdR0-sdR);
        if (dSd < sdCriteria)
          break;
      }
    }

    # Remove "bad" weights. For instance, too large w's gives:
    #   Error in smooth.spline(g$x, g$yin, w = w, ...) :
    #     NA/NaN/Inf in foreign function call (arg 4)
    ok.weights <- (w != 0 & is.finite(w));
    if (!all(ok.weights))
      w[!ok.weights] <- 0;
    # Not needed anymore
    ok.weights <- NULL;

    g <- stats.smooth.spline(g$x, g$yin, w=w, ..., tol=tol);

    # Not needed anymore
    w <- NULL;

    if (plotCurves == TRUE)
      lines(g, col=(col<-col+1));

    sdR0 <- sdR;
  } # while ( ... )

  g$iter <- iter
  
  g;
}) # robustSmoothSpline()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.whichUnique <- function(x, ..., tol) {
  # We need to make sure that 'g$x == x' below. /HB 2011-10-10
  xx <- x;
  keep <- rep(TRUE, times=length(x));
  while (TRUE) {
    idxs <- which(keep);
    xx <- round((x[idxs] - mean(x[idxs]))/tol);  # de-mean to avoid possible overflow
    dups <- duplicated(xx);
    if (!any(dups)) {
      break;
    }
    keep[idxs[dups]] <- FALSE;
  } # while()
  nd <- keep;

  # Sanity check
  stopifnot(length(nd) == length(x));

  which(nd);
} # .whichUnique()


######################################################################
# HISTORY
# 2015-01-06
# o Using requestNamespace() instead of request().
# 2014-03-25
# o CLEANUP: The internal .Fortran() calls no longer pass DUP=FALSE,
#   which "may be disabled in future versions of R.".
# 2013-09-26
# o Now utilizing anyMissing().
# 2012-08-30
# o BUG FIX: Now local getNativeSplineFitFunction() sets up the
#   function such that it is called via a FortranRoutine object,
#   rather than by name.
# 2012-08-19
# o Added local getNativeSplineFitFunction() function to
#   robustSmoothSpline() which returns a wrapper to a proper
#   native and internal spline fit function of R.
# o Make it clear that robustSmoothSpline() is under GPL (>= 2),
#   because it is adapted from smooth.spline() of R by R Core Team.
#   Added a GPL source code header.
# 2011-10-10
# o Updated robustSmoothSpline() such that it works with the new
#   "uniqueness" scheme of smooth.spline() in R v2.14.0 and newer.
#   It is tricky, because robustSmoothSpline() is a reiterative
#   algorithm which requires that the choosen "unique" x:s does
#   not change in each iteration.  Previously, 'signif(x, 6)' was
#   used to identify unique x:s, which gives the same set of values
#   when called twice, whereas this is not true for the new choice
#   with 'round((x - mean(x))/tol)'.
# 2011-04-12
# o Now using as.double(NA) instead of NA, which is logical.
# o Interestingly, stats::smooth.spline() of R v2.14.0 now does
#   very similar speedups as robustSmoothSpline() has done
#   internally in its smooth.spline.fit() since 2002.  Great.
# o CLEANUP: Now robustSmoothSpline() utilizes stats:::n.knots()
#   internally, if running on R v2.14.0 or newer.
# 2008-07-20
# o MEMORY OPTIMIZATION: Removing more variables when done etc.
#   Helping the garbage collector by doing x <- as.vector(x) before
#   calling a function rather than having as.vector(x) as an argument.
# 2007-06-08
# o Added declaration 'nx <- 0' in robustSmoothSpline.matrix() in
#   order to please R CMD check R v2.6.0.
# 2007-01-01
# o Removed any code to make method backward compatibility with
#   R < 1.9.0, which was before 'modreg' was merged into 'stats'.
# 2005-06-03
# o Now making use of setMethodS3().
# o Renamed to robustSmoothSpline().
# o Copied from R.basic. At the same time, we speedup functions were made
#   into local functions.
# 2005-01-24
# o Added support for weights.
# 2002-04-21
# o Updated due to modreg is merged into stats from R v1.9.0.
# 2002-03-02
# o SPEED UP: Since robust. smooth. spline() now makes use of
#   the "home-made" smooth.spline.prepare() and smooth.spline0() it
#   speed up about three times on my test data; 32secs -> 9secs.
# o Splitted smooth.spline() into the two functions
#   smooth.spline.prepare() and smooth.spline.fit(). The purpose of
#   this is to speed up robust.spline(), especially when there are
#   duplicate x values!
# 2002-02-19
# o The idea of using w.org is not simple since the data points are
#   reorder by smooth.spline.
# o Made w <- as.vector(w).
# 2002-02-18
# o Created the Rd comments with an example adapted from
#   smooth.spline.
# o Made it possible to specify weights even in the robust estimation.
# o Added a assertion that the weights are non-illegal and not to
#   big.
# o Renamed to robust. smooth. spline() and made analogue to
#   smooth.spline().
# 2002-02-15
# o Created. It seems like the robust spline alorithm gives pretty
#   much the same result as lowess. If not, the differences are
#   quite small compared to the noise level of cDNA microarray data.
######################################################################
