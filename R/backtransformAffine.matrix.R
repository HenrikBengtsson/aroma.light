#########################################################################/**
# @set "class=matrix"
# @RdocMethod backtransformAffine
#
# @title "Reverse affine transformation"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{X}{An NxK @matrix containing data to be backtransformed.}
#  \item{a}{A scalar, @vector, a @matrix, or a @list.
#    First, if a @list, it is assumed to contained the elements \code{a}
#    and \code{b}, which are the used as if they were passed as seperate
#    arguments.
#    If a @vector, a matrix of size NxK is created which is then filled
#    \emph{row by row} with the values in the vector. Commonly, the
#    vector is of length K, which means that the matrix will consist of
#    copies of this vector stacked on top of each other.
#    If a @matrix, a matrix of size NxK is created which is then filled
#    \emph{column by column} with the values in the matrix (collected
#    column by column. Commonly, the matrix is of size NxK, or NxL with
#    L < K and then the resulting matrix consists of copies sitting
#    next to each other.
#    The resulting NxK matrix is subtracted from the NxK matrix \code{X}.
#  }
#  \item{b}{A scalar, @vector, a @matrix.
#    A NxK matrix is created from this argument. For details see
#    argument \code{a}.
#    The NxK matrix \code{X-a} is divided by the resulting NxK matrix.
#  }
#  \item{project}{
#    returned (K values per data point are returned).
#    If @TRUE, the backtransformed values "\code{(X-a)/b}" are projected
#    onto the line L(a,b) so that all columns
#    will be identical.
#  }
#  \item{...}{Not used.}
# }
#
# \value{
#   The "\code{(X-a)/b}" backtransformed NxK @matrix is returned.
#   If \code{project} is @TRUE, an Nx1 @matrix is returned, because
#   all columns are identical anyway.
# }
#
# \section{Missing values}{
#   Missing values remain missing values. If projected, data points that
#   contain missing values are projected without these.
# }
#
# @examples "../incl/backtransformAffine.matrix.Rex"
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("backtransformAffine", "matrix", function(X, a=NULL, b=NULL, project=FALSE, ...) {

  # Dimensions of 'X'
  nobs  <- nrow(X);
  ndims <- ncol(X);
  if (ndims == 1) {
    stop("Can not fit affine multiscan model. Matrix must contain at least two columns (scans): ", ndims);
  }

  # If argument 'a' is a list assume it contains the elements 'a' and 'b'.
  if (is.list(a)) {
    b <- a$b;
    a <- a$a;
  }

  # If 'a' and/or 'b' are vector convert them to row matrices.
  if (is.vector(a)) {
    # Create a full matrix and filled row by row with 'a'
    a <- matrix(a, nrow=nobs, ncol=ndims, byrow=TRUE);
  } else if (is.matrix(a)) {
    # Create a full matrix and filled column by column by the columns in 'a'
    t <- a;
    naValue <- as.double(NA);
    a <- matrix(naValue, nrow=nobs, ncol=ndims);
    for (cc in 1:ndims) {
      # Loop over the columns in a0 too.
      col <- ((cc-1) %% ncol(t)) + 1;
      value <- rep(t[,col], length.out=nobs);
      a[,cc] <- value;
    }
    # Not needed anymore
    t <- NULL;
  } else if (!is.null(a)) {
    stop(paste("Unknown data type of argument 'a':", class(a)[1]));
  }

  if (!project) {
    if (is.vector(b)) {
      # Create a full matrix and filled row by row with 'b'
      b <- matrix(b, nrow=nobs, ncol=ndims, byrow=TRUE);
    } else if (is.matrix(b)) {
      # Create a full matrix and filled column by column by the columns in 'b'
      t <- b;
      naValue <- as.double(NA);
      b <- matrix(naValue, nrow=nobs, ncol=ndims);
      for (cc in 1:ndims) {
        # Loop over the columns in a0 too.
        col <- ((cc-1) %% ncol(t)) + 1;
        value <- rep(t[,col], length.out=nobs);
        b[,cc] <- value;
      }
    } else if (!is.null(b)) {
      stop(paste("Unknown data type of argument 'b':", class(b)[1]));
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Subtract the bias and rescale
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(a))
    X <- X - a;

  ########################################################################
  # Alternative 2
  #
  #   i) Translate the fitted line L to L' such that L' goes through
  #      (a,a,...,a) and project the data y onto L' to obtain ytilde
  #  ii) from which we can calculate xtilde = (ytilde - a) / b
  # iii) Since all xtilde are the same take the first component to be our
  #      estimate xhat in y = a + b*xhat.
  ########################################################################
  if (project) {
    # In theory:
    #  ytilde <- projectUontoV(y-a,b) + a;
    #  xtilde <- (ytilde-a)/b;
    # In practice:
    X <- t(X);

    # 'b' standardized to a unit vector such that <v,v> == 1.
    v <- b / sqrt(sum(b^2));

    # projectUontoV(): U should be an NxK matrix and v an N vector.
    X <- projectUontoV(X,v, na.rm=TRUE);

    X <- X[1,] / b[1];                     # Note that here 'b' is a vector!
  } else {
    if (!is.null(b))
      X <- X / b;                          # Note that here 'b' is a matrix!
  }

  as.matrix(X);
}) # backtransformAffine()

############################################################################
# HISTORY:
# 2011-04-12
# o Now using as.double(NA) instead of NA.
# 2006-06-03
# o Minor to merge two different threads of this code.
# o Method passes the tests in the example code (again).
# 2005-01-24
# o Now missing values are excluded before projection.
# 2005-01-08
# o Now, if project is TRUE, only one column is returned.
# o Now the method is guaranteed to return a matrix by calling as.matrix().
# o The average of projected data is now no the same scale as the average on
#   non-projected data. Before the data was max(b) times too large.
# o Added argument 'project' to the Rdoc comments.
# o Added test in example to assert that the same matrix is returned if
#   projection on an identity transformation is applied.
# 2004-06-28
# o BUG FIX: Applied projection when project was FALSE and vice versa.
#   Projection did not work because 'b' was expanded into a full matrix.
#   Now this is only done if project == FALSE.
# 2004-05-14
# o Created. Extracted and generalize code from calibrateMultiscan(),
#   normalizeAffine() and calibrateMultiscanSpatial().
############################################################################

