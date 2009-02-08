principal.curve.hb <- function(x, start=NULL, thresh=0.001, plot.true=FALSE, maxit=10, stretch=2, smoother="smooth.spline", trace=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments:
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'smoother':
  if (is.function(smoother)) {
    smootherFcn <- smoother;
  } else {
    smooth.list <- c("smooth.spline", "lowess", "periodic.lowess");
    smoother <- match.arg(smoother, smooth.list);
    smootherFcn <- NULL;
  }

  # Argument 'stretch':
  stretches <- c(2, 2, 0);
  if (is.function(smoother)) {
    if (is.null(stretch))
      stop("Argument 'stretch' must be given if 'smoother' is a function.");
  } else {
    if(missing(stretch) || is.null(stretch)) {
      stretch <- stretches[match(smoother, smooth.list)];
    }
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(smootherFcn)) {
    # Setup the smoother function smootherFcn(lambda, xj, ...) which must
    # return fitted {y}:s.
    smootherFcn <- switch(smoother,
      lowess = function(lambda, xj, ...) {
        lowess(lambda, xj, ...)$y;
      },
  
      smooth.spline = function(lambda, xj, ..., df=5) {
        o <- order(lambda);
        lambda <- lambda[o];
        xj <- xj[o];
        fit <- smooth.spline(lambda, xj, ..., df=df, keep.data=FALSE);
        predict(fit, x=lambda)$y;
      },
  
      periodic.lowess = function(lambda, xj, ...) {
        periodic.lowess(lambda, xj, ...)$y;
      }
    ) # smootherFcn()

    # Should the fitted curve be bias corrected (in each iteration)?
    biasCorrectCurve <- (smoother == "periodic.lowess");
  } else {
    biasCorrectCurve <- FALSE;
  }



  this.call <- match.call()
  dist.old <- sum(diag(var(x)))
  d <- dim(x)
  n <- d[1]
  p <- d[2] 

  # You can give starting values for the curve
  if (missing(start) || is.null(start)) {
    # use largest principal component
    if (is.character(smoother) && smoother == "periodic.lowess") {
      start <- startCircle(x)
    } else {
      xbar <- colMeans(x)
      xstar <- scale(x, xbar, FALSE)
      svd.xstar <- svd(xstar)
      dd <- svd.xstar$d
      lambda <- svd.xstar$u[,1] * dd[1]
      tag <- order(lambda)
      s <- scale(outer(lambda, svd.xstar$v[,1]),  - xbar, FALSE)
      dist <- sum((dd^2)[-1]) * n
      start <- list(s=s, tag=tag, lambda=lambda, dist=dist)
    }
  } else if (!inherits(start, "principal.curve")) {
    # use given starting curve 
    if (is.matrix(start)) {
      start <- get.lam(x, start, stretch=stretch)
    } else {
      stop("Invalid starting curve: should be a matrix or principal.curve")
    }
  }

  pcurve <- start
  if (plot.true) {
    plot(x[,1:2], xlim=adjust.range(x[,1], 1.3999999999999999),
	 ylim=adjust.range(x[,2], 1.3999999999999999))
    lines(pcurve$s[pcurve$tag, 1:2])
  }

  it <- 0
  if (trace) {
    cat("Starting curve---distance^2: ", pcurve$dist, "\n", sep="");
  }


  # Pre-allocate nxp matrix 's'
  s <- matrix(NA, nrow=n, ncol=p);

  hasConverged <- (abs((dist.old - pcurve$dist)/dist.old) <= thresh);
  while (!hasConverged && it < maxit) {
    it <- it + 1;

    for(jj in 1:p) {
      s[,jj] <- smootherFcn(pcurve$lambda, x[,jj], ...);
    }

    dist.old <- pcurve$dist;


    # Finds the "projection index" for a matrix of points 'x',
    # when projected onto a curve 's'.  The projection index,
    # \lambda_f(x) [Eqn (3) in Hastie & Stuetzle (1989), is
    # the value of \lambda for which f(\lambda) is closest 
    # to x.
    pcurve <- get.lam(x, s, stretch=stretch);

    # Bias correct?
    if (biasCorrectCurve)
      pcurve <- bias.correct.curve(x, pcurve, ...)

    # Converged?
    hasConverged <- (abs((dist.old - pcurve$dist)/dist.old) <= thresh);

    if (plot.true) {
      plot(x[,1:2], xlim=adjust.range(x[,1], 1.3999999999999999), ylim=adjust.range(x[,2], 1.3999999999999999))
      lines(pcurve$s[pcurve$tag, 1:2])
    }

    if (trace) {
      cat("Iteration ", it, "---distance^2: ", pcurve$dist, "\n", sep="");
    }
  } # while()

  # Return fit
  structure(list(
    s = pcurve$s, 
    tag = pcurve$tag, 
    lambda = pcurve$lambda, 
    dist = pcurve$dist,
    converged = hasConverged,         # Added by HB
    nbrOfIterations = as.integer(it), # Added by HB
    call = this.call
  ), class="principal.curve");
} # principal.curve.hb()



###########################################################################
# HISTORY:
# 2009-02-08
# o BUG FIX: An error was thrown if 'smoother' was a function.
# o Cleaned up source code (removed comments).
# 2008-05-27
# o Benchmarking: For larger data sets, most of the time is spent in
#   get.lam().
# o BUG FIX: smooth.spline(x,y) will only use *and* return values for
#   "unique" {x}:s. This means that the fitted {y}:s maybe be fewer than
#   the input vector. In order to control for this, we use predict().
# o Now 'smoother' can also be a function taking arguments 'lambda', 'xj'
#   and '...' and return 'y' of the same length as 'lambda' and 'xj'.
# o Now arguments 'start' and 'stretch' can be NULL, which behaves the 
#   same as if they are "missing" [which is hard to emulate with for 
#   instance do.call()].
# o Added 'converged' and 'nbrOfIterations' to return structure.
# o SPEED UP/MEMORY OPTIMIZATION: Now the nxp matrix 's' is allocated only 
#   once. Before it was built up using cbind() once per iteration.
# o SPEED UP: Now the smoother function is identified/created before 
#   starting the algorithm, and not once per dimension and iteration.
###########################################################################
