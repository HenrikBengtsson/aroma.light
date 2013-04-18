%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  wpca.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{wpca.matrix}
\alias{wpca.matrix}

\title{Light-weight Weighted Principal Component Analysis}

\usage{
\method{wpca}{matrix}(x, w=NULL, center=TRUE, scale=FALSE, method=c("dgesdd", "dgesvd", "dsvdc"), swapDirections=FALSE, ...)
}

\description{
  Calculates the (weighted) principal components of a matrix, that is,
  finds a new coordinate system (not unique) for representing the given
  multivariate data such that
   i) all dimensions are orthogonal to each other, and
  ii) all dimensions have maximal variances.
}

\arguments{
  \item{x}{An NxK \code{\link[base]{matrix}}.}
  \item{w}{An N \code{\link[base]{vector}} of weights for each row (observation) in
    the data matrix. If \code{\link[base]{NULL}}, all observations get the same weight,
    that is, standard PCA is used.}
  \item{center}{If \code{\link[base:logical]{TRUE}}, the (weighted) sample mean column \code{\link[base]{vector}} is
    subtracted from each column in \code{mat}, first.
    If data is not centered, the effect will be that a linear subspace
    that goes through the origin is fitted.}
  \item{scale}{If \code{\link[base:logical]{TRUE}}, each column in \code{mat} is
    divided by its (weighted) root-mean-square of the
    centered column, first.}
  \item{method}{If \code{"dgesdd"} LAPACK's divide-and-conquer
    based SVD routine is used (faster [1]), if \code{"dgesvd"}, LAPACK's
    QR-decomposition-based routine is used, and if \code{"dsvdc"},
    LINPACK's DSVDC(?) routine is used. The latter is just for
    pure backward compatibility with R v1.7.0.
  }
  \item{swapDirections}{If \code{\link[base:logical]{TRUE}}, the signs of eigenvectors
    that have more negative than positive components are inverted.
    The signs of corresponding principal components are also inverted.
    This is only of interest when for instance visualizing or comparing
    with other PCA estimates from other methods, because the
    PCA (SVD) decompostion of a matrix is not unique.
  }
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{list}} with elements:
  \item{pc}{An NxK \code{\link[base]{matrix}} where the column \code{\link[base]{vector}}s are the
            principal components (a.k.a. loading vectors,
            spectral loadings or factors etc).}
  \item{d}{An K \code{\link[base]{vector}} containing the eigenvalues of the
            principal components.}
  \item{vt}{An KxK \code{\link[base]{matrix}} containing the eigenvector of the
            principal components.}
  \item{xMean}{The center coordinate.}

  It holds that \code{x == t(t(fit$pc \%*\% fit$vt) + fit$xMean)}.
}

\section{Method}{
  A singular value decomposition (SVD) is carried out.
  Let X=\code{mat}, then the SVD of the matrix is \eqn{X = U D V'}, where
  \eqn{U} and \eqn{V} are othogonal, and \eqn{D} is a diagonal matrix
  with singular values. The principal returned by this method are \eqn{U D}.

  Internally \code{La.svd()} (or \code{svd()}) of the \pkg{base}
  package is used.
  For a popular and well written introduction to SVD see for instance [2].
}

\examples{
  for (zzz in 0) {

# This example requires plot3d() in R.basic [http://www.braju.com/R/]
if (!require(pkgName <- "R.basic", character.only=TRUE)) break

# -------------------------------------------------------------
# A first example
# -------------------------------------------------------------
# Simulate data from the model y <- a + bx + eps(bx)
x <- rexp(1000)
a <- c(2,15,3)
b <- c(2,3,15)
bx <- outer(b,x)
eps <- apply(bx, MARGIN=2, FUN=function(x) rnorm(length(x), mean=0, sd=0.1*x))
y <- a + bx + eps
y <- t(y)

# Add some outliers by permuting the dimensions for 1/3 of the observations
idx <- sample(1:nrow(y), size=1/3*nrow(y))
y[idx,] <- y[idx,c(2,3,1)]

# Down-weight the outliers W times to demonstrate how weights are used
W <- 10

# Plot the data with fitted lines at four different view points
N <- 4
theta <- seq(0,180,length=N)
phi <- rep(30, length.out=N)

# Use a different color for each set of weights
col <- topo.colors(W)

opar <- par(mar=c(1,1,1,1)+0.1)
layout(matrix(1:N, nrow=2, byrow=TRUE))
for (kk in seq(theta)) {
  # Plot the data
  plot3d(y, theta=theta[kk], phi=phi[kk])
 
  # First, same weights for all observations
  w <- rep(1, length=nrow(y))

  for (ww in 1:W) {
    # Fit a line using IWPCA through data
    fit <- wpca(y, w=w, swapDirections=TRUE)
 
    # Get the first principal component
    ymid <- fit$xMean
    d0 <- apply(y, MARGIN=2, FUN=min) - ymid
    d1 <- apply(y, MARGIN=2, FUN=max) - ymid
    b <- fit$vt[1,]
    y0 <- -b * max(abs(d0))
    y1 <-  b * max(abs(d1))
    yline <- matrix(c(y0,y1), nrow=length(b), ncol=2)
    yline <- yline + ymid
 
    points3d(t(ymid), col=col)
    lines3d(t(yline), col=col)

    # Down-weight outliers only, because here we know which they are.
    w[idx] <- w[idx]/2
  }

  # Highlight the last one
  lines3d(t(yline), col="red", lwd=3)
}

par(opar)

} # for (zzz in 0)
rm(zzz)


  if (dev.cur() > 1) dev.off()

  # -------------------------------------------------------------
# A second example
# -------------------------------------------------------------
# Data
x <- c(1,2,3,4,5)
y <- c(2,4,3,3,6)

opar <- par(bty="L")
opalette <- palette(c("blue", "red", "black"))
xlim <- ylim <- c(0,6)

# Plot the data and the center mass
plot(x,y, pch=16, cex=1.5, xlim=xlim, ylim=ylim)
points(mean(x), mean(y), cex=2, lwd=2, col="blue")


# Linear regression y ~ x
fit <- lm(y ~ x)
abline(fit, lty=1, col=1)

# Linear regression y ~ x through without intercept
fit <- lm(y ~ x - 1)
abline(fit, lty=2, col=1)


# Linear regression x ~ y
fit <- lm(x ~ y)
c <- coefficients(fit)
b <- 1/c[2]
a <- -b*c[1]
abline(a=a, b=b, lty=1, col=2)

# Linear regression x ~ y through without intercept
fit <- lm(x ~ y - 1)
b <- 1/coefficients(fit)
abline(a=0, b=b, lty=2, col=2)


# Orthogonal linear "regression"
fit <- wpca(cbind(x,y))

b <- fit$vt[1,2]/fit$vt[1,1]
a <- fit$xMean[2]-b*fit$xMean[1]
abline(a=a, b=b, lwd=2, col=3)

# Orthogonal linear "regression" without intercept
fit <- wpca(cbind(x,y), center=FALSE)
b <- fit$vt[1,2]/fit$vt[1,1]
a <- fit$xMean[2]-b*fit$xMean[1]
abline(a=a, b=b, lty=2, lwd=2, col=3)

legend(xlim[1],ylim[2], legend=c("lm(y~x)", "lm(y~x-1)", "lm(x~y)",
          "lm(x~y-1)", "pca", "pca w/o intercept"), lty=rep(1:2,3),
                     lwd=rep(c(1,1,2),each=2), col=rep(1:3,each=2))

palette(opalette)
par(opar)

}

\author{Henrik Bengtsson}

\references{
  [1] J. Demmel and  J. Dongarra, \emph{DOE2000 Progress Report}, 2004.
      \url{http://www.cs.berkeley.edu/~demmel/DOE2000/Report0100.html} \cr
  [2] Todd Will,
      \emph{Introduction to the Singular Value Decomposition},
      UW-La Crosse, 2004.
      \url{http://www.uwlax.edu/faculty/will/svd/} \cr
}

\seealso{
  For a iterative re-weighted PCA method, see \code{\link[aroma.light:iwpca.matrix]{*iwpca}()}.
  For Singular Value Decomposition, see \code{\link[base]{svd}}().
  For other implementations of Principal Component Analysis functions see
  (if they are installed):
  \code{\link[stats]{prcomp}} in package \pkg{stats} and
  \code{pca()} in package \pkg{pcurve}.
}


\keyword{methods}
\keyword{algebra}
