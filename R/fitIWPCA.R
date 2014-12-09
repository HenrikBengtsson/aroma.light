########################################################################/**
# @RdocGeneric fitIWPCA
# @alias fitIWPCA.matrix
#
# @title "Robust fit of linear subspace through multidimensional data"
#
# \description{
#   @get "title".
# }
#
# \usage{
# @usage fitIWPCA,matrix
# }
#
# \arguments{
#  \item{X}{NxK @matrix where N is the number of observations and
#           K is the number of dimensions (channels).
#  }
#
#  \item{constraint}{A @character string or a @numeric value.
#   If @character it specifies which additional contraint to be used
#   to specify the offset parameters along the fitted line;
#
#   If \code{"diagonal"}, the offset vector will be a point on the line
#   that is closest to the diagonal line (1,...,1).
#   With this constraint, all bias parameters are identifiable.
#
#   If \code{"baseline"} (requires argument \code{baselineChannel}), the
#   estimates are such that of the bias and scale parameters of the
#   baseline channel is 0 and 1, respectively.
#   With this constraint, all bias parameters are identifiable.
#
#   If \code{"max"}, the offset vector will the point on the line that is
#   as "great" as possible, but still such that each of its components is
#   less than the corresponding minimal signal. This will guarantee that
#   no negative signals are created in the backward transformation.
#   If @numeric value, the offset vector will the point on the line
#   such that after applying the backward transformation there are
#   \code{constraint*N}. Note that \code{constraint==0} corresponds
#   approximately to \code{constraint=="max"}.
#   With the latter two constraints, the bias parameters are only
#   identifiable modulo the fitted line.
#  }
#
#  \item{baselineChannel}{Index of channel toward which all other
#    channels are conform.
#    This argument is required if \code{constraint=="baseline"}.
#    This argument is optional if \code{constraint=="diagonal"} and
#    then the scale factor of the baseline channel will be one. The
#    estimate of the bias parameters is not affected in this case.
#    Defaults to one, if missing.
#  }
#
#  \item{...}{Additional arguments accepted by @see "iwpca".
#   For instance, a N @vector of weights for each observation may be
#    given, otherwise they get the same weight.
#  }
#
#  \item{aShift, Xmin}{For internal use only.}
# }
#
# \value{
#   Returns a @list that contains estimated parameters and algorithm
#   details;
#
#   \item{a}{A @double @vector \eqn{(a[1],...,a[K])}with offset
#       parameter estimates.
#       It is made identifiable according to argument \code{constraint}.
#   }
#   \item{b}{A @double @vector \eqn{(b[1],...,b[K])}with scale
#       parameter estimates.  It is made identifiable by constraining
#       \code{b[baselineChannel] == 1}.
#       These estimates are idependent of argument \code{constraint}.
#   }
#   \item{adiag}{If identifiability constraint \code{"diagonal"},
#       a @double @vector \eqn{(adiag[1],...,adiag[K])}, where
#       \eqn{adiag[1] = adiag[2] = ... adiag[K]}, specifying the point
#       on the diagonal line that is closest to the fitted line,
#       otherwise the zero vector.
#   }
#   \item{eigen}{A KxK @matrix with columns of eigenvectors.
#   }
#   \item{converged}{@TRUE if the algorithm converged, otherwise @FALSE.
#   }
#   \item{nbrOfIterations}{The number of iterations for the algorithm
#                           to converge, or zero if it did not converge.
#   }
#
#   \item{t0}{Internal parameter estimates, which contains no more
#                            information than the above listed elements.
#   }
#   \item{t}{Always @NULL.}
# }
#
# \details{
#   This method uses re-weighted principal component analysis (IWPCA)
#   to fit a the nodel \eqn{y_n = a + bx_n + eps_n} where \eqn{y_n},
#   \eqn{a}, \eqn{b}, and \eqn{eps_n} are vector of the K and \eqn{x_n}
#   is a scalar.
#
#   The algorithm is:
#    For iteration i:
#    1) Fit a line \eqn{L} through the data close using weighted PCA
#       with weights \eqn{\{w_n\}}. Let
#         \eqn{r_n = \{r_{n,1},...,r_{n,K}\}}
#       be the \eqn{K} principal components.
#    2) Update the weights as
#         \eqn{w_n <- 1 / \sum_{2}^{K} (r_{n,k} + \epsilon_r)}
#       where we have used the residuals of all but the first principal
#       component.
#    3) Find the point a on \eqn{L} that is closest to the
#       line \eqn{D=(1,1,...,1)}. Similarily, denote the point on D that is
#       closest to \eqn{L} by \eqn{t=a*(1,1,...,1)}.
# }
#
# @author
#
# %examples "fitMultiIWPCA.matrix.Rex"
#
# \seealso{
#   This is an internal method used by the @see "calibrateMultiscan"
#   and @see "normalizeAffine" methods.
#   Internally the function @see "iwpca" is used to fit a line
#   through the data cloud and the function @see "distanceBetweenLines" to
#   find the closest point to the diagonal (1,1,...,1).
# }
#
# @keyword "algebra"
#*/########################################################################
setMethodS3("fitIWPCA", "matrix", function(X, constraint=c("diagonal", "baseline", "max"), baselineChannel=NULL, ..., aShift=rep(0, times=ncol(X)), Xmin=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 0. Define local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  statistic <- function(X, ..., constraint="diagonal", Xmin=NULL,
                               baselineChannel=1, aShift=rep(0, times=ncol(X))) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit an K-dimensional line through the data using iterative
    # re-weighted PCA.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fit <- iwpca(X, ...);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the fitted line L
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get the center of the fitted line...
    ax <- fit$xMean;
    names(ax) <- NULL;

    # ...and the fitted eigenvectors (u1,u2,...,uK)
    # with ui*uj = 0; i!=j and ui*ui = 1.
    U <- t(fit$vt);
    colnames(U) <- rownames(U) <- NULL;

    # The fitted scale parameters b=(b[1],b[2],...,b[K]) where
    # the elements are rescaled such that b[1] == 1.
    # [ min(b[i]) == 1. Before it was such that b[1] == 1, but this
    #   is probably better.                                               ]
    # [ Indeed not; this is not good if one do more than one estimate per
    #   array, e.g. printtip etc. /HB 2004-04-26                          ]
    # [ With the introduction of 'baselineChannel' it is possible to
    #   specify which channel should get scale one. / HB 2004-06-30       ]
    U1 <- U[,1];
    bx <- as.vector(U1/U1[baselineChannel]);

    # Shift the data.
    # [ This is for instance useful if fitting towards the diagonal line
    #   and resampling under H0: y_i = alpha + z_i and /HB 2004-01-02     ]
    ax <- ax + aShift;

    if (identical(constraint, "diagonal")) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the point t on the fitted line that is closest to the
      # points s on the "diagonal" line (1,1,...,1) in K-space.
      # [This works also for lines in two dimension.]
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # x(s) is the fitted line (the first IWPCA component)
      # y(t) is the diagonal line
      ay <- rep(0, times=length(ax));                         # (0,0,...,0)
      by <- rep(1, times=length(ay));                         # (1,1,...,1)

      dbl <- distanceBetweenLines(ax=ax,bx=bx, ay=ay,by=by);
      a     <- as.vector(dbl$xs);
      adiag <- as.vector(dbl$yt);
    } else if (identical(constraint, "baseline")) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the point t on the fitted line for which the bias parameter of
      # the baseline channel is zero, i.e. for which a[baselineChannel]==0.
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # The scale parameters are already such that b[baselineChannel]==1.
      #   y[c,i] = a[c] + b[c]*x[c,i] ; c = 1,...,C
      #   y[b,i] = a[b] +      x[b,i] ; b - baseline channel
      # Similar to the constraint=="max" reasoning:
      # For channel b, find t such that
      #   ax[b] + bx[b]*t == 0 <=> { bx[b]==1 } <=> t = -ax[b]
      # => a[c] <- ax[c] + bx[c]*t  <=> a[c] <- ax[c] - bx[c]*ax[b]
      a <- ax - bx*ax[baselineChannel];
      adiag <- rep(0, times=length(ax));
    } else if (identical(constraint, "max")) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the "greatest" point t on the fitted line that is within
      # the cube C whose upper limits are defined by the minimum value
      # in each channel.
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the minimal value of each X component.
      if (is.null(Xmin))
        Xmin <- colMins(X, na.rm=TRUE);
      # For each component k, find the value t such that
      #   ax[k] + bx[k]*t[k] == Xmin[k] <=> t[k] == (Xmin[k] - ax[k])/bx[k]
      t <- (Xmin-ax)/bx;
      # Choose minimum t[k]
      # Now, amax <<= Xmin if amax <- ax - bx[k]*min(t) where <<= (\prec)
      # means componentswise less or equal than.
      a <- ax + bx*min(t);
      adiag <- rep(0, times=length(ax));
    } else if (is.numeric(constraint)) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the "greatest" point t on the fitted line that is within
      # the cube C whose upper limits are defined by the alpha quantile
      # value in each channel.
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Find the alpha quantile value of each X component.
      if (is.null(Xmin)) {
        Xmin <- colQuantiles(X, probs=constraint, na.rm=TRUE);
      }
      # For each component k, find the value t such that
      #   ax[k] + bx[k]*t[k] == Xmin[k] <=> t[k] == (Xmin[k] - ax[k])/bx[k]
      t <- (Xmin-ax)/bx;
      # Choose minimum t[k]
      # Now, amax <<= Xmin if amax <- ax - bx[k]*min(t) where <<= (\prec)
      # means componentswise less or equal than.
      a <- ax + bx*min(t);
      adiag <- rep(0, times=length(ax));
    }

    # Return the statistic
    t <- c(a=a);
    t <- c(t, b=bx);
    t <- c(t, adiag=adiag);
    t <- c(t, U=as.vector(U));
    t <- c(t, niter=fit$nbrOfIterations * (fit$converged*2-1));
    t;
  } # statistic()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Verify the arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'X'
  if (!is.matrix(X))
    stop("Argument 'X' must be a matrix:", mode(X));
  N <- nrow(X);
  K <- ncol(X);
  if (K == 1) {
    stop("Argument 'X' must have two or more columns:", K);
  }
  if (N < K) {
    stop("Argument 'X' must have at least as many rows as columns:",
                                                                 N, "<", K);
  }

  # Argument: 'constraint'
  if (is.numeric(constraint)) {
    if (length(constraint) != 1)
      stop("Argument 'constraint' can not be a numerical vector.");
    if (constraint < 0 || constraint > 1)
      stop("Invalid value of argument 'constraint':", constraint);
  } else {
    constraint <- match.arg(constraint);
    if (identical(constraint, "baseline")) {
      if (is.null(baselineChannel)) {
        stop("Argument 'baselineChannel' must be given if 'constraint' is ",
                                                           "\"baseline\".");
      }
    }
  }

  # Argument: 'baselineChannel'
  if (!is.null(baselineChannel)) {
    if (!is.numeric(baselineChannel) || length(baselineChannel) != 1) {
      stop("Argument 'baselineChannel' must be a single numeric: ",
                                                           baselineChannel);
    }

    if (baselineChannel < 1 || baselineChannel > ncol(X)) {
      stop("Argument 'baselineChannel' is out of range [1,", ncol(X),"]: ",
                                                           baselineChannel);
    }

    if (!(constraint %in% c("baseline", "diagonal"))) {
      stop("Argument 'baselineChannel' must not be specified if ",
                 "argument 'constraint' is \"baseline\" or \"diagonal\": ",
                                                               constraint);
    }

    if (!is.null(Xmin)) {
      stop("Argument 'Xmin' must not be specified if 'baselineChannel' is",
                               " specified: ", paste(Xmin, collapse=", "));
    }
  } else {
    baselineChannel <- 1;
  }

  # Argument: 'aShift'
  if (is.null(aShift)) {
    aShift <- rep(0, times=ncol(X));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Prepare the data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(constraint, "max")) {
     Xmin <- colMins(X, na.rm=TRUE);
  } else {
     Xmin <- NULL;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Fit the model
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Use only finite observations
  isFinite <- apply(X, MARGIN=1L, FUN=function(r) all(is.finite(r)));

  # Number of finite observations
  N  <- sum(isFinite);

  # Validate the number of finite observations
  if (N < K) {
    stop("Argument 'X' must have at least as many non-finite ",
                         "observations ", "(rows) as columns:", N, "<", K);
  }

  t0 <- statistic(X[isFinite,], constraint=constraint,
          Xmin=Xmin, baselineChannel=baselineChannel, aShift=aShift, ...);
  t  <- NULL;

  # Extract the parameter estimates from the internal estimation vector.
  a <- t0[regexpr("^a[0-9]*$", names(t0)) != -1];

  b <- t0[regexpr("^b[0-9]*$", names(t0)) != -1];

  adiag <- t0[regexpr("^adiag[0-9]*$", names(t0)) != -1];

  U <- t0[regexpr("^U[0-9]*$", names(t0)) != -1];
  U <- matrix(U, nrow=sqrt(length(U)));

  niter <- as.integer(abs(t0["niter"]));
  converged <- (niter > 0);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 5. Return the parameter estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  list(a=a, b=b, adiag=adiag, eigen=U,
                   converged=converged, nbrOfIterations=niter, t0=t0, t=t);
}) # fitIWPCA()


###########################################################################
# HISTORY:
# 2013-09-26
# o Now utilizing colMins() and colQuantiles() of 'matrixStats'.
# 2011-02-05
# o DOCUMENTATION: Fixed broken links to help for iwpca().
# 2006-01-22
# o Added Rdoc help on the returned parameters.
# o If missing, 'baselineChannel' is now set to one before calling the
#   internal function.
# o Now fitIWPCA() does not return the data matrix. This is to save memory.
#   The calling algorithm can equally well add the data if it is needed.
# 2004-06-30
# o Added argument 'baselineChannel' with 'constraint' "baseline". This
#   is useful for instance when normalizing toward a common reference. In
#   such cases, the common-reference channel is unaffected and only the
#   other channel(s) is affinely transformed.
# 2004-06-28
# o BUG FIX: Forgot to exclude non-finite observations when fitting the
#   iwpca(). This was done before, but somehow it disappeared while
#   re-organizing the code.
###########################################################################
