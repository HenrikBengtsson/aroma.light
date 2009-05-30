%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  backtransformPrincipalCurve.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{backtransformPrincipalCurve.matrix}
\alias{backtransformPrincipalCurve.matrix}
\alias{backtransformPrincipalCurve.matrix}

\alias{backtransformPrincipalCurve.numeric}

\title{Reverse transformation of principal-curve fit}

\description{
  Reverse transformation of principal-curve fit.
}

\usage{\method{backtransformPrincipalCurve}{matrix}(X, fit, dimensions=NULL, targetDimension=NULL, ...)}

\arguments{
 \item{X}{An NxK \code{\link[base]{matrix}} containing data to be backtransformed.}
 \item{fit}{An MxL principal-curve fit object of class
   \code{principal.curve} as returned by \code{\link[aroma.light:fitPrincipalCurve.matrix]{*fitPrincipalCurve}()}.
   Typically \eqn{L = K}, but not always.
 }
 \item{dimensions}{An (optional) subset of of D dimensions all in [1,L]
   to be returned (and backtransform).}
 \item{targetDimension}{An (optional) index specifying the dimension
   in [1,L] to be used as the target dimension of the \code{fit}.
   More details below.}
 \item{...}{Passed internally to \code{\link[stats]{smooth.spline}}.}
}

\value{
  The backtransformed NxK (or NxD) \code{\link[base]{matrix}}.
}

\section{Target dimension}{
  By default, the backtransform is such that afterward the signals are
  approximately proportional to the (first) principal curve as fitted
  by \code{\link[aroma.light:fitPrincipalCurve.matrix]{*fitPrincipalCurve}()}.  This scale and origin of this
  principal curve is not uniquely defined.
  If \code{targetDimension} is specified, then the backtransformed signals
  are approximately proportional to the signals of the target dimension,
  and the signals in the target dimension are unchanged.
}

\section{Subsetting dimensions}{
  Argument \code{dimensions} can be used to backtransform a subset of
  dimensions (K) based on a subset of the fitted dimensions (L).
  If \eqn{K = L}, then both \code{X} and \code{fit} is subsetted.
  If \eqn{K <> L}, then it is assumed that \code{X} is already
  subsetted/expanded and only \code{fit} is subsetted.
}

\examples{

# Consider the case where K=4 measurements have been done
# for the same underlying signals 'x'.  The different measurements
# have different systematic variation
#
#   y_k = f(x_k) + eps_k; k = 1,...,K.
#
# In this example, we assume non-linear measurement functions
#
#   f(x) = a + b*x + x^c + eps(b*x)
#
# where 'a' is an offset, 'b' a scale factor, and 'c' an exponential.
# We also assume heteroscedastic zero-mean noise with standard
# deviation proportional to the rescaled underlying signal 'x'.
#
# Furthermore, we assume that measurements k=2 and k=3 undergo the
# same transformation, which may illustrate that the come from
# the same batch. However, when *fitting* the model below we
# will assume they are independent.

# Transforms
a <- c(2, 15, 15,   3)
b <- c(2,  3,  3,   4)
c <- c(1,  2,  2, 1/2)
K <- length(a)

# The true signal
N <- 1000
x <- rexp(N)

# The noise
bX <- outer(b,x)
E <- apply(bX, MARGIN=2, FUN=function(x) rnorm(K, mean=0, sd=0.1*x))

# The transformed signals with noise
Xc <- t(sapply(c, FUN=function(c) x^c))
Y <- a + bX + Xc + E
Y <- t(Y)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit principal curve
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit principal curve through Y = (y_1, y_2, ..., y_K)
fit <- fitPrincipalCurve(Y)

# Flip direction of 'lambda'?
rho <- cor(fit$lambda, Y[,1], use="complete.obs")
flip <- (rho < 0)
if (flip) {
  fit$lambda <- max(fit$lambda, na.rm=TRUE)-fit$lambda
}

L <- ncol(fit$s)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Backtransform data according to model fit
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Backtransform toward the principal curve (the "common scale")
YN1 <- backtransformPrincipalCurve(Y, fit=fit)
stopifnot(ncol(YN1) == K)


# Backtransform toward the first dimension
YN2 <- backtransformPrincipalCurve(Y, fit=fit, targetDimension=1)
stopifnot(ncol(YN2) == K)


# Backtransform toward the last (fitted) dimension
YN3 <- backtransformPrincipalCurve(Y, fit=fit, targetDimension=L)
stopifnot(ncol(YN3) == K)


# Backtransform toward the third dimension (dimension by dimension)
# Note, this assumes that K == L.
YN4 <- Y
for (cc in 1:L) {
  YN4[,cc] <- backtransformPrincipalCurve(Y, fit=fit,
                                  targetDimension=1, dimensions=cc)
}
stopifnot(identical(YN4, YN2))


# Backtransform a subset toward the first dimension
# Note, this assumes that K == L.
YN5 <- backtransformPrincipalCurve(Y, fit=fit,
                               targetDimension=1, dimensions=2:3)
stopifnot(identical(YN5, YN2[,2:3]))
stopifnot(ncol(YN5) == 2)


# Extract signals from measurement #2 and backtransform according
# its model fit.  Signals are standardized to target dimension 1.
y6 <- Y[,2,drop=FALSE]
yN6 <- backtransformPrincipalCurve(y6, fit=fit, dimensions=2,
                                               targetDimension=1)
stopifnot(identical(yN6, YN2[,2,drop=FALSE]))
stopifnot(ncol(yN6) == 1)


# Extract signals from measurement #2 and backtransform according
# the the model fit of measurement #3 (because we believe these
# two have undergone very similar transformations.
# Signals are standardized to target dimension 1.
y7 <- Y[,2,drop=FALSE]
yN7 <- backtransformPrincipalCurve(y7, fit=fit, dimensions=3,
                                               targetDimension=1)
stopifnot(ncol(yN7) == 1)
stopifnot(cor(yN7, yN6) > 0.9999)
}

\seealso{
  \code{\link[aroma.light:fitPrincipalCurve.matrix]{*fitPrincipalCurve}()}
}
\keyword{methods}
