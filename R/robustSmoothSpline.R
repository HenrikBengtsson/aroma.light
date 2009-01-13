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
#   @see "stats::smooth.spline".
# }
#
# @author
#
# @keyword "smooth"
# @keyword "robust"
#*/############################################################################
setMethodS3("robustSmoothSpline", "default", function(x, y=NULL, w=NULL, ..., minIter=3, maxIter=max(minIter, 50), sdCriteria=2e-4, reps=1e-15, plotCurves=FALSE) {
  require(stats) || throw("Package not loaded: stats");  # smooth.spline()

  # To please RMD CMD check for R v2.6.0
  nx <- 0;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  smooth.spline.prepare <- function(x, w=NULL, df=5, spar=NULL, cv=FALSE, all.knots=FALSE, df.offset=0, penalty=1, control.spar=list()) {
    sknotl <- function(x) {
      n.kn <- function(n) {
        if (n < 50) 
          n
        else trunc({
          a1 <- log(50, 2)
          a2 <- log(100, 2)
          a3 <- log(140, 2)
          a4 <- log(200, 2)
          if (n < 200) 2^(a1 + (a2 - a1) * (n - 50)/150) else if (n < 
            800) 2^(a2 + (a3 - a2) * (n - 200)/600) else if (n < 
            3200) 2^(a3 + (a4 - a3) * (n - 800)/2400) else 200 + 
            (n - 3200)^0.2
        })
      }
      nk <- n.kn(n <- length(x))
      c(rep(x[1], 3), x[seq(1, n, len = nk)], rep(x[n], 3))
    }
    contr.sp <- list(low=-1.5, high=1.5, tol=1e-04, eps=2e-08, maxit=500,
                     trace=getOption("verbose"))
    contr.sp[(namc <- names(control.spar))] <- control.spar
    if (!all(sapply(contr.sp[1:4], is.double)) || contr.sp$tol < 0 ||
             contr.sp$eps <= 0 || contr.sp$maxit <= 0) 
      stop("invalid `control.spar'")
    # ------------ Differences from smooth.spline BEGIN -----------
    n <- length(x)
    w <- if (is.null(w)) 
      rep(1, n)
    else {
      if (n != length(w)) 
        stop("lengths of x and w must match")
      if (any(w < 0)) 
        stop("all weights should be non-negative")
      if (all(w == 0)) 
        stop("some weights should be positive")
      (w * sum(w > 0))/sum(w)
    }
    x <- signif(x, 6)
    ux <- unique(sort(x))
    nx <- length(ux)
    if (nx <= 3) 
        stop("need at least for unique `x' values")
    if (cv && nx < n) 
        warning("crossvalidation with non-unique `x' seems doubtful")
    # ------------ Differences from smooth.spline BEGIN -----------
    if (nx == n)
      ox <- 1:n
    else
      ox <- match(x, ux)
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
      rm(w); rm(y);
      tmp <- unlist(tmp, use.names=FALSE);
      tmp <- matrix(tmp, ncol=3, byrow=TRUE);
      wbar <- tmp[, 1]
      ybar <- tmp[, 2]/ifelse(wbar > 0, wbar, 1)
      yssw <- sum(tmp[, 3] - wbar * ybar^2)
      rm(tmp);  # Cleanup / HB 2008-07-20
    }
  
    nk <- prep$nk
  
    fit <- .Fortran("qsbart", as.double(prep$penalty), as.double(prep$dofoff), 
      x=as.double(prep$xbar), y=as.double(ybar), w=as.double(wbar), 
      ssw=as.double(yssw), as.integer(nx), as.double(prep$knot), 
      as.integer(prep$nk), coef=double(nk), ty=double(nx), lev=double(nx), 
      crit=double(1), iparms=prep$iparms, spar=prep$spar, parms=unlist(prep$contr.sp[1:4]), 
      isetup=as.integer(0), scrtch=double((17 + nk) * nk), ld4=as.integer(4),
      ldnk=as.integer(1), ier=integer(1), DUP=FALSE,
      PACKAGE="stats");
    # Clean up. /HB 2008-07-20
    rm(prep, yssw, nk);

    fields <- c("coef", "ty", "lev", "spar", "parms", "crit", "iparms", "ier");
    fit <- fit[fields];
    fit$wbar <- wbar;
    fit$ybar <- ybar;
    rm(wbar, ybar);

    fit;
  } # smooth.spline.fit()


  smooth.spline0 <- function(x, y=NULL, w=NULL, df=5, spar=NULL, cv=FALSE, all.knots=FALSE, df.offset=0, penalty=1, control.spar=list()) {
    if (inherits(x, "smooth.spline.prepare")) {
      prep <- x;
    } else {
      xy <- xy.coords(x,y);
      prep <- smooth.spline.prepare(x=xy$x, w=w, df=df, spar=spar, cv=cv, all.knots=all.knots, df.offset=df.offset, penalty=penalty, control.spar=control.spar);
      y <- xy$y;
      rm(xy);
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
    rm(fit); # /HB 2008-07-20

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
      rm(ww)
    } else {
      weighted.mean((y - ty[ox])^2, w)/(1 - (df.offset + penalty * df)/n)^2
    }
  
    pen.crit <- sum(wbar * (ybar - ty) * ybar)

    knot <- prep$knot;
    nk <- prep$nk;
    ux <- prep$ux;
    r.ux <- prep$r.ux;
    rm(prep);

    fit.object <- list(knot=knot, nk=nk, min=ux[1], range=r.ux, coef=coef)
    class(fit.object) <- "smooth.spline.fit"
    rm(knot, nk, r.ux);
    object <- list(x=ux, y=ty, w=wbar, yin=ybar, lev=lev, cv.crit=cv.crit,
               pen.crit=pen.crit, crit=crit, df=df, spar=spar,
               lambda=lambda, iparms=iparms, fit=fit.object,
               call=match.call())
    rm(ux, ty, wbar, ybar,  lev, cv.crit, pen.crit, crit, df, spar, 
       lambda, fit.object); # /HB 2007-08-20
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
  tx <- signif(x, 6);
  utx <- unique(sort(tx));
  otx <- match(utx, tx);
  w0 <- w[otx];
  rm(tx, utx, otx); # /HB 2008-07-20

  if (inherits(x, "smooth.spline")) {
    g <- x;
  } else if (missing(w) || is.null(w)) {
    x <- as.vector(x);
    y <- as.vector(y);
    g <- smooth.spline(x, y, ...);
    rm(x,y); # HB /2008-07-20
  } else {
#    warning("Robust *weighted* smoothing splines are not implemented yet!")
    x <- as.vector(x);
    y <- as.vector(y);
    w <- as.vector(w);
    g <- smooth.spline(x, y, w=w, ...);
    rm(x,y,w); # HB /2008-07-20
  }
  
  # Precalculate a lot of thing for speeding up subsequent calls to smooth.spline()
  spline.prep <- smooth.spline.prepare(x=g$x, w=g$w,...);

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
  sdR0 <- NA;
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
    rm(r);  # HB /2008-07-20

    if (iter > minIter) {
      if (!is.na(sdR0)) {
        dSd <- abs(sdR0-sdR);
#       cat(", ", formatC(dSd, digits=3), sep="");
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
    rm(ok.weights); # HB /2008-07-20

    g <- smooth.spline0(spline.prep, g$yin, w=w, ...);
    rm(w); # HB /2008-07-20

    if (plotCurves == TRUE)
      lines(g, col=(col<-col+1));

    sdR0 <- sdR;
  }

  g;
}) # robustSmoothSpline()


######################################################################
# HISTORY
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
