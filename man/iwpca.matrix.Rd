%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  iwpca.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{iwpca.matrix}
\alias{iwpca.matrix}

\title{Fits an R-dimensional hyperplane using iterative re-weighted PCA}

\description{
  Fits an R-dimensional hyperplane using iterative re-weighted PCA.
}

\usage{
\method{iwpca}{matrix}(X, w=NULL, R=1, method=c("symmetric", "bisquare", "tricube", "L1"), maxIter=30, acc=1e-04, reps=0.02, fit0=NULL, ...)
}

\arguments{
 \item{X}{N-times-K \code{\link[base]{matrix}} where N is the number of observations and
   K is the number of dimensions.}
 \item{w}{An N \code{\link[base]{vector}} of weights for each row (observation) in
   the data matrix. If \code{\link[base]{NULL}}, all observations get the same weight.}
 \item{R}{Number of principal components to fit. By default a line
   is fitted.}
 \item{method}{
   If \code{"symmetric"} (or \code{"bisquare"}), Tukey's biweight
   is used. If \code{"tricube"}, the tricube weight is used.
   If \code{"L1"}, the model is fitted in \eqn{L_1}.
   If a \code{\link[base]{function}}, it is used to calculate weights for next iteration
   based on the current iteration's residuals.}
 \item{maxIter}{Maximum number of iterations.}
 \item{acc}{The (Euclidean) distance between two subsequent parameters
   fit for which the algorithm is considered to have converged.}
 \item{reps}{Small value to be added to the residuals before the
   the weights are calculated based on their inverse. This is to avoid
   infinite weights.}
 \item{fit0}{A \code{\link[base]{list}} containing elements \code{vt} and \code{pc}
   specifying an initial fit.
   If \code{\link[base]{NULL}}, the initial guess will be equal to the (weighted) PCA fit.}
 \item{...}{Additional arguments accepted by \code{\link[aroma.light:wpca.matrix]{*wpca}()}.}
}

\value{
  Returns the fit (a \code{\link[base]{list}}) from the last call to \code{\link[aroma.light:wpca.matrix]{*wpca}()}
  with the additional elements \code{nbrOfIterations} and
  \code{converged}.
}

\details{
  This method uses weighted principal component analysis (WPCA) to fit a
  R-dimensional hyperplane through the data with initial internal
  weights all equal.
  At each iteration the internal weights are recalculated based on
  the "residuals".
  If \code{method=="L1"}, the internal weights are 1 / sum(abs(r) + reps).
  This is the same as \code{method=function(r) 1/sum(abs(r)+reps)}.
  The "residuals" are orthogonal Euclidean distance of the principal
  components R,R+1,...,K.
  In each iteration before doing WPCA, the internal weighted are
  multiplied by the weights given by argument \code{w}, if specified.
}

\author{Henrik Bengtsson}

\examples{
for (zzz in 0) {

# This example requires plot3d() in R.basic [http://www.braju.com/R/]
if (!require(pkgName <- "R.basic", character.only=TRUE)) break

# Simulate data from the model y <- a + bx + eps(bx)
x <- rexp(1000)
a <- c(2,15,3)
b <- c(2,3,4)
bx <- outer(b,x)
eps <- apply(bx, MARGIN=2, FUN=function(x) rnorm(length(x), mean=0, sd=0.1*x))
y <- a + bx + eps
y <- t(y)

# Add some outliers by permuting the dimensions for 1/10 of the observations
idx <- sample(1:nrow(y), size=1/10*nrow(y))
y[idx,] <- y[idx,c(2,3,1)]

# Plot the data with fitted lines at four different view points
opar <- par(mar=c(1,1,1,1)+0.1)
N <- 4
layout(matrix(1:N, nrow=2, byrow=TRUE))
theta <- seq(0,270,length=N)
phi <- rep(20, length.out=N)
xlim <- ylim <- zlim <- c(0,45);
persp <- list();
for (kk in seq(theta)) {
  # Plot the data
  persp[[kk]] <- plot3d(y, theta=theta[kk], phi=phi[kk], xlim=xlim, ylim=ylim, zlim=zlim)
}

# Weights on the observations
# Example a: Equal weights
w <- NULL
# Example b: More weight on the outliers (uncomment to test)
w <- rep(1, length(x)); w[idx] <- 0.8

# ...and show all iterations too with different colors.
maxIter <- c(seq(1,20,length=10),Inf)
col <- topo.colors(length(maxIter))
# Show the fitted value for every iteration
for (ii in seq(along=maxIter)) {
  # Fit a line using IWPCA through data
  fit <- iwpca(y, w=w, maxIter=maxIter[ii], swapDirections=TRUE)

  ymid <- fit$xMean
  d0 <- apply(y, MARGIN=2, FUN=min) - ymid
  d1 <- apply(y, MARGIN=2, FUN=max) - ymid
  b <- fit$vt[1,]
  y0 <- -b * max(abs(d0))
  y1 <-  b * max(abs(d1))
  yline <- matrix(c(y0,y1), nrow=length(b), ncol=2)
  yline <- yline + ymid

  for (kk in seq(theta)) {
    # Set pane to draw in
    par(mfg=c((kk-1) \%/\% 2, (kk-1) \%\% 2) + 1);
    # Set the viewpoint of the pane
    options(persp.matrix=persp[[kk]]);

    # Get the first principal component
    points3d(t(ymid), col=col[ii])
    lines3d(t(yline), col=col[ii])

    # Highlight the last one
    if (ii == length(maxIter))
      lines3d(t(yline), col="red", lwd=3)
  }
}

par(opar)

} # for (zzz in 0)
rm(zzz)
}

\seealso{
  Internally \code{\link[aroma.light:wpca.matrix]{*wpca}()} is used for calculating the weighted PCA.
}


\keyword{methods}
\keyword{algebra}
