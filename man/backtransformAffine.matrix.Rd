%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  backtransformAffine.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{backtransformAffine.matrix}
\alias{backtransformAffine.matrix}

\title{Reverse affine transformation}

\description{
  Reverse affine transformation.
}

\usage{
\method{backtransformAffine}{matrix}(X, a=NULL, b=NULL, project=FALSE, ...)
}

\arguments{
 \item{X}{An NxK \code{\link[base]{matrix}} containing data to be backtransformed.}
 \item{a}{A scalar, \code{\link[base]{vector}}, a \code{\link[base]{matrix}}, or a \code{\link[base]{list}}.
   First, if a \code{\link[base]{list}}, it is assumed to contained the elements \code{a}
   and \code{b}, which are the used as if they were passed as seperate
   arguments.
   If a \code{\link[base]{vector}}, a matrix of size NxK is created which is then filled
   \emph{row by row} with the values in the vector. Commonly, the
   vector is of length K, which means that the matrix will consist of
   copies of this vector stacked on top of each other.
   If a \code{\link[base]{matrix}}, a matrix of size NxK is created which is then filled
   \emph{column by column} with the values in the matrix (collected
   column by column. Commonly, the matrix is of size NxK, or NxL with
   L < K and then the resulting matrix consists of copies sitting
   next to each other.
   The resulting NxK matrix is subtracted from the NxK matrix \code{X}.
 }
 \item{b}{A scalar, \code{\link[base]{vector}}, a \code{\link[base]{matrix}}.
   A NxK matrix is created from this argument. For details see
   argument \code{a}.
   The NxK matrix \code{X-a} is divided by the resulting NxK matrix.
 }
 \item{project}{
   returned (K values per data point are returned).
   If \code{\link[base:logical]{TRUE}}, the backtransformed values "\code{(X-a)/b}" are projected
   onto the line L(a,b) so that all columns
   will be identical.
 }
 \item{...}{Not used.}
}

\value{
  The "\code{(X-a)/b}" backtransformed NxK \code{\link[base]{matrix}} is returned.
  If \code{project} is \code{\link[base:logical]{TRUE}}, an Nx1 \code{\link[base]{matrix}} is returned, because
  all columns are identical anyway.
}

\section{Missing values}{
  Missing values remain missing values. If projected, data points that
  contain missing values are projected without these.
}

\examples{
X <- matrix(1:8, nrow=4, ncol=2)
X[2,2] <- NA

print(X)

# Returns a 4x2 matrix
print(backtransformAffine(X, a=c(1,5)))

# Returns a 4x2 matrix
print(backtransformAffine(X, b=c(1,1/2)))

# Returns a 4x2 matrix
print(backtransformAffine(X, a=matrix(1:4,ncol=1)))

# Returns a 4x2 matrix
print(backtransformAffine(X, a=matrix(1:3,ncol=1)))

# Returns a 4x2 matrix
print(backtransformAffine(X, a=matrix(1:2,ncol=1), b=c(1,2)))

# Returns a 4x1 matrix
print(backtransformAffine(X, b=c(1,1/2), project=TRUE))

# If the columns of X are identical, and a identity
# backtransformation is applied and projected, the
# same matrix is returned.
X <- matrix(1:4, nrow=4, ncol=3)
Y <- backtransformAffine(X, b=c(1,1,1), project=TRUE)
print(X)
print(Y)
stopifnot(sum(X[,1]-Y) <= .Machine$double.eps)


# If the columns of X are identical, and a identity
# backtransformation is applied and projected, the
# same matrix is returned.
X <- matrix(1:4, nrow=4, ncol=3)
X[,2] <- X[,2]*2; X[,3] <- X[,3]*3;
print(X)
Y <- backtransformAffine(X, b=c(1,2,3))
print(Y)
Y <- backtransformAffine(X, b=c(1,2,3), project=TRUE)
print(Y)
stopifnot(sum(X[,1]-Y) <= .Machine$double.eps)
}

\seealso{
  For more information see \code{\link{matrix}}.
}
\keyword{methods}
