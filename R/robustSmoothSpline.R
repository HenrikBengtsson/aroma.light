#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

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
#   Because of this, this function is also licensed under GPL~v2.
# }
#
# @author
#
# @keyword "smooth"
# @keyword "robust"
#*/############################################################################
setMethodS3("robustSmoothSpline", "default", function(x, y=NULL, w=NULL, ..., minIter=3, maxIter=max(minIter, 50), sdCriteria=2e-4, reps=1e-15, tol=1e-6*IQR(x), plotCurves=FALSE) {
  require(stats) || throw("Package not loaded: stats");  # smooth.spline()

  # To please RMD CMD check for R v2.6.0
  nx <- 0;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getNativeSplineFitFunction <- function() {
    # Locate all native Fortran routines
    pkgName <- "stats";
    fcns <- getDLLRegisteredRoutines(pkgName)$.Fortran;

    # Starting with R v2.15.1 patched (rev 60026)
    key <- "rbart";
    if (is.element(key, names(fcns))) {
      nparams <- fcns[[key]]$numParameters;
      if (nparams == 20) {
        expr <- parse(text="stats:::C_rbart");
        routine <- eval(expr);
        # Sanity check
        stopifnot(inherits(routine, "FortranRoutine"));
        fcn <- function(prep, ybar, wbar, yssw, nx, nk, ...) {
          .Fortran(routine,
             as.double(prep$penalty), as.double(prep$dofoff),
             x=as.double(prep$xbar), y=as.double(ybar), w=as.double(wbar),
             ssw=as.double(yssw), as.integer(nx), as.double(prep$knot),
             as.integer(prep$nk), coef=double(nk), ty=double(nx),
             lev=double(nx), crit=double(1), iparms=prep$iparms,
             spar=prep$spar, parms=unlist(prep$contr.sp[1:4]),
             scratch=double(17L * nk + 1L),
             ld4=4L, ldnk=1L, ier=integer(1L),
             DUP=TRUE, PACKAGE=pkgName);
        } # fcn()
        return(fcn);
      }

      throw(sprintf("Non-supported number of parameters for internal spline function %s(): %d", key, nparams));
    }

    # Prior to R v2.15.1 patched (rev 60026)
    key <- "qsbart";
    if (is.element(key, names(fcns))) {
      nparams <- fcns[[key]]$numParameters;
      if (nparams == 21) {
        fcn <- function(prep, ybar, wbar, yssw, nx, nk, ...) {
          .Fortran(key, as.double(prep$penalty), as.double(prep$dofoff),
             x=as.double(prep$xbar), y=as.double(ybar), w=as.double(wbar),
             ssw=as.double(yssw), as.integer(nx), as.double(prep$knot),
             as.integer(prep$nk), coef=double(nk), ty=double(nx),
             lev=double(nx), crit=double(1), iparms=prep$iparms,
             spar=prep$spar, parms=unlist(prep$contr.sp[1:4]),
             isetup=as.integer(0),
             scrtch=double((17 + nk) * nk),
             ld4=as.integer(4), ldnk=as.integer(1), ier=integer(1),
             DUP=FALSE, PACKAGE=pkgName);
        } # fcn()
        return(fcn);
      }

      throw(sprintf("Non-supported number of parameters for internal spline function %s(): %d", key, nparams));
    }

    # Failed to locate native function
    pd <- packageDescription("aroma.light");
    throw(sprintf("INTERNAL ERROR of robustSmoothSline(): Failed to locate an internal spline function for %s. Please report this to the package maintainer (%s) of %s.", R.version$version.string, pd$Maintainer, pd$Package));
  } # getNativeSplineFitFunction()


  ## Former internal n.kn() is now available as n.knots() in stats v2.14.0.
  rVer <- getRversion();
  if (rVer >= "2.14.0") {
    ## Cannot use n.knots <- stats:::n.knots because then
    ## R CMD check will complain with R 2.13.x and before.
    n.knots <- getAnywhere("n.knots")$obj[[1]];
    # Sanity check
    stopifnot(is.function(n.knots));

    whichUnique <- function(x, ...) {
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
    } # whichUnique()

    stats.smooth.spline <- smooth.spline;
  } else {
    n.knots <- function(n) {
        ## Number of inner knots
        if (n < 50L) n
        else trunc({
            a1 <- log2( 50)
            a2 <- log2(100)
            a3 <- log2(140)
            a4 <- log2(200)
            if	(n < 200L) 2^(a1+(a2-a1)*(n-50)/150)
            else if (n < 800L) 2^(a2+(a3-a2)*(n-200)/600)
            else if (n < 3200L) 2^(a3+(a4-a3)*(n-800)/2400)
            else 200 + (n-3200)^0.2
        })
    } # n.knots()

    whichUnique <- function(x, ...) {
      tx <- signif(x, 6);
      utx <- unique(sort(tx));
      otx <- match(utx, tx);
      otx;
    } # whichUnique()

    stats.smooth.spline <- function(..., tol) {
      smooth.spline(...);
    }
  } # if (rVer ...)


  smooth.spline.prepare <- function(x, w=NULL, df=5, spar=NULL, cv=FALSE, all.knots=FALSE, df.offset=0, penalty=1, control.spar=list(), tol=1e-6*IQR(x)) {
    sknotl <- function(x) {
      nk <- n.knots(n <- length(x))
      c(rep(x[1], 3), x[seq(1, n, len = nk)], rep(x[n], 3))
    }

    contr.sp <- list(low = -1.5, # low = 0.      was default till R 1.3.x
                     high = 1.5,
                     tol = 1e-4, # tol = 0.001   was default till R 1.3.x
                     eps = 2e-8, # eps = 0.00244 was default till R 1.3.x
                     maxit = 500, trace = getOption("verbose"));

    contr.sp[names(control.spar)] <- control.spar

    if (!all(sapply(contr.sp[1:4], is.numeric)) || contr.sp$tol < 0 ||
             contr.sp$eps <= 0 || contr.sp$maxit <= 0)
      stop("invalid `control.spar'")
    # ------------ Differences from smooth.spline BEGIN -----------
    n <- length(x)
    w <- if (is.null(w))
      rep(1, n)
    else {
      if (n != length(w))
        stop("lengths of 'x' and 'w' must match")
      if (any(w < 0))
        stop("all weights should be non-negative")
      if (all(w == 0))
        stop("some weights should be positive")
      (w * sum(w > 0))/sum(w)
    }

    uIdxs <- whichUnique(x);
    ux <- sort(x[uIdxs]);
    nx <- length(ux);
    if (nx == n) {
      ox <- 1:n;
    } else {
      ox <- match(x, ux);
    }

    if (nx <= 3)
        stop("need at least four unique `x' values")
    if (cv && nx < n)
        warning("crossvalidation with non-unique `x' seems doubtful")
    # ------------ Differences from smooth.spline BEGIN -----------
    # ...
    # ------------ Differences from smooth.spline END -----------
    r.ux <- ux[nx] - ux[1]
    xbar <- (ux - ux[1])/r.ux
    if (all.knots) {
      knot <- c(rep(xbar[1], 3), xbar, rep(xbar[nx], 3))
      nk <- nx + 2
    } else {
      knot <- sknotl(xbar)
      nk <- length(knot) - 4
    }
    ispar <- if (is.null(spar) || missing(spar)) {
      if (contr.sp$trace) -1 else 0
    } else
      1
    spar <- if (ispar == 1) as.double(spar) else double(1)
    icrit <- if (cv) 2 else 1
    dofoff <- df.offset
    if (!missing(df)) {
      if (df > 1 && df <= nx) {
        icrit <- 3
        dofoff <- df
      } else
        warning(paste("you must supply 1 < df <= n,  n = #{unique x} =", nx))
    }
    iparms <- as.integer(c(icrit, ispar, contr.sp$maxit))
    names(iparms) <- c("icrit", "ispar", "iter")

    object <- list(penalty=penalty, dofoff=dofoff, xbar=as.double(xbar), nx=nx, knot=knot, nk=nk, iparms=iparms, spar=spar, contr.sp=contr.sp, ox=ox, n=n, df.offset=df.offset, w=w, ux=ux, r.ux=r.ux);
    class(object) <- "smooth.spline.prepare";
    object;
  } # smooth.spline.prepare()



  smooth.spline.fit <- function(prep, y=NULL) {
    nx <- prep$nx
    n <- prep$n
    if (nx == n) {
      # Don't call tapply if not necessary. / HB 2002-03-02
      wbar <- prep$w
      ybar <- y
      yssw <- rep(0, n)
    } else {
      w <- prep$w
      ox <- prep$ox
      # The tapply is expensive! / HB 2002-03-02
##      tmp <- matrix(unlist(tapply(seq(along=y), ox, function(i, y, w) {
##           c(sum(w[i]), sum(w[i]*y[i]), sum(w[i]*y[i]^2))
##         }, y=y, w=w)), ncol=3, byrow = TRUE)
      tmp <- tapply(seq(along=y), INDEX=ox, FUN=function(i, y, w) {
           c(sum(w[i]), sum(w[i]*y[i]), sum(w[i]*y[i]^2))
         }, y=y, w=w);
      # Not needed anymore
      w <- y <- NULL;
      tmp <- unlist(tmp, use.names=FALSE);
      tmp <- matrix(tmp, ncol=3, byrow=TRUE);
      wbar <- tmp[, 1]
      ybar <- tmp[, 2]/ifelse(wbar > 0, wbar, 1)
      yssw <- sum(tmp[, 3] - wbar * ybar^2)
      # Not needed anymore
      tmp <- NULL;
    }

    nk <- prep$nk

    fit <- fitFcn(prep=prep, ybar=ybar, wbar=wbar, yssw=yssw, nx=nx, nk=nk);

    # Not needed anymore
    prep <- yssw <- nk <- NULL;

    fields <- c("coef", "ty", "lev", "spar", "parms", "crit", "iparms", "ier");
    fit <- fit[fields];
    fit$wbar <- wbar;
    fit$ybar <- ybar;
    # Not needed anymore
    wbar <- ybar <- NULL;

    fit;
  } # smooth.spline.fit()


  smooth.spline0 <- function(x, y=NULL, w=NULL, df=5, spar=NULL, cv=FALSE, all.knots=FALSE, df.offset=0, penalty=1, control.spar=list()) {
    if (inherits(x, "smooth.spline.prepare")) {
      prep <- x;
    } else {
      xy <- xy.coords(x,y);
      prep <- smooth.spline.prepare(x=xy$x, w=w, df=df, spar=spar, cv=cv, all.knots=all.knots, df.offset=df.offset, penalty=penalty, control.spar=control.spar);
      y <- xy$y;
      # Not needed anymore
      xy <- NULL;
    }

    fit <- smooth.spline.fit(prep, y=y);

    lev <- fit$lev;
    wbar <- fit$wbar;
    coef <- fit$coef;
    spar <- fit$spar;
    iparms <- fit$iparms;
    crit <- fit$crit;
    lambda <- unname(fit$parms["low"]);
    ybar <- fit$ybar;
    ty <- fit$ty;
    ier <-fit$ier;
    # Not needed anymore
    fit <- NULL;

    df <- sum(lev)
    if (is.na(df))
      stop("NA lev[]; probably smoothing parameter `spar' way too large!")
    if (ier > 0) {
      sml <- (spar < 0.5)
      wtxt <- paste("smoothing parameter value too", if (sml) "small" else "large")
      if (sml) {
        stop(wtxt)
      } else {
        ty <- rep(mean(y), nx)
        df <- 1
        warning(paste(wtxt, "setting df = 1  __use with care!__", sep="\n"))
      }
    }

    ox <- prep$ox;
    n <- prep$n;

    cv.crit <- if (cv) {
      ww <- wbar
      ww[!(ww > 0)] <- 1
      weighted.mean(((y - ty[ox])/(1 - (lev[ox] * w)/ww[ox]))^2, w)
      # Not needed anymore
      ww <- NULL;
    } else {
      weighted.mean((y - ty[ox])^2, w)/(1 - (df.offset + penalty * df)/n)^2
    }

    pen.crit <- sum(wbar * (ybar - ty) * ybar)

    knot <- prep$knot;
    nk <- prep$nk;
    ux <- prep$ux;
    r.ux <- prep$r.ux;
    # Not needed anymore
    prep <- NULL;

    fit.object <- list(knot=knot, nk=nk, min=ux[1], range=r.ux, coef=coef)
    class(fit.object) <- "smooth.spline.fit"
    # Not needed anymore
    knot <- nk <- r.ux <- NULL;
    object <- list(x=ux, y=ty, w=wbar, yin=ybar, lev=lev, cv.crit=cv.crit,
               pen.crit=pen.crit, crit=crit, df=df, spar=spar,
               lambda=lambda, iparms=iparms, fit=fit.object,
               call=match.call())
    # Not needed anymore
    ux <- ty <- wbar <- ybar <-  lev <- cv.crit <- pen.crit <- crit <- df <- spar <- lambda <- fit.object <- NULL;
    class(object) <- "smooth.spline"

    object
  } # smooth.spline0()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Verify arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'w'
  if (is.numeric(w)) {
    w <- as.double(w);
    if (any(is.na(w))) {
      stop("Weights with value NA are not allowed.");
    }
    if (any(w < 0 | w > 1)) {
      stop("Weights out of range [0,1]: ",
           paste(w[w < 0.0 | w > 1.0], collapse=", "));
    }
  } else if (!is.null(w)) {
    stop("Argument 'w' is of an unsupported datatype/class: ",
                                                         class(weights)[1]);
  }

  # Argument: 'reps'
  if (!is.numeric(reps) || length(reps) != 1 || reps <= 0)
    throw("Argument 'reps' must be a single postive number.");

  # smooth.spline() next will only operate on unique x-values. For this reason,
  # we have to remove corresponding weights too. There is a small problem here;
  # if different weights are used for data points (x,y,w) with same x-value, which
  # data points (x,y,w) should be used? Here we use the first one only. /HB 2005-01-24


  uIdxs <- whichUnique(x);
  nu <- length(uIdxs);
  w0 <- w[uIdxs];

  # WORKAROUND
  if (rVer >= "2.14.0") {
    # We need to make sure that 'g$x == x' below. /HB 2011-10-10
    x <- x[uIdxs];
    y <- y[uIdxs];
    w <- w[uIdxs];
    uIdxs <- seq(along=x);
  }

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

  # Precalculate a lot of thing for speeding up subsequent calls
  # to smooth.spline()
  spline.prep <- smooth.spline.prepare(x=g$x, w=g$w,...);

  # Sanity check
  stopifnot(with(spline.prep, {length(w) == length(ux)}));


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
  # Get a wrapper to the internal .Fortran() function
  # which is not part of the public API of R. /HB 2012-08-19
  fitFcn <- getNativeSplineFitFunction();

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

    g <- smooth.spline0(spline.prep, g$yin, w=w, ...);
    # Not needed anymore
    w <- NULL;

    if (plotCurves == TRUE)
      lines(g, col=(col<-col+1));

    sdR0 <- sdR;
  } # while ( ... )

  g;
}) # robustSmoothSpline()


######################################################################
# HISTORY
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
