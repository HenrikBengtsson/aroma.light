%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitPrincipalCurve.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitPrincipalCurve.matrix}
\alias{fitPrincipalCurve.matrix}
\alias{fitPrincipalCurve.matrix}


\encoding{latin1}

\title{Fit a principal curve in K dimensions}

\description{
  Fit a principal curve in K dimensions.
}

\usage{\method{fitPrincipalCurve}{matrix}(X, ..., verbose=FALSE)}

\arguments{
 \item{X}{An NxK \code{\link[base]{matrix}} (K>=2) where the columns represent the dimension.}
 \item{...}{Other arguments passed to \code{\link[princurve]{principal.curve}}.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a principal.curve object (which is a \code{\link[base]{list}}).
  See \code{\link[princurve]{principal.curve}} for more details.
}

\section{Missing values}{
 The estimation of the affine normalization function will only be made
 based on complete observations, i.e. observations that contains no \code{\link[base]{NA}}
 values in any of the channels.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] Hastie, T. and Stuetzle, W, \emph{Principal Curves}, JASA, 1989.
}

\examples{

# Simulate data from the model y <- a + bx + x^c + eps(bx)
x <- rexp(1000)
a <- c(2,15,3)
b <- c(2,3,4)
c <- c(1,2,1/2)
bx <- outer(b,x)
xc <- t(sapply(c, FUN=function(c) x^c))
eps <- apply(bx, MARGIN=2, FUN=function(x) rnorm(length(x), mean=0, sd=0.1*x))
y <- a + bx + xc + eps
y <- t(y)

# Fit principal curve through (y_1, y_2, y_3)
fit <- fitPrincipalCurve(y)

# Flip direction of 'lambda'?
rho <- cor(fit$lambda, y[,1], use="complete.obs")
flip <- (rho < 0)
if (flip) {
  fit$lambda <- max(fit$lambda, na.rm=TRUE)-fit$lambda
}


# Backtransform (y_1, y_2, y_3) to be proportional to each other
yN <- backtransformPrincipalCurve(y, fit=fit)

# Same backtransformation dimension by dimension
yN2 <- y
for (cc in 1:ncol(y)) {
  yN2[,cc] <- backtransformPrincipalCurve(y, fit=fit, dimensions=cc)
}
stopifnot(identical(yN2, yN))


xlim <- c(0, 1.04*max(x))
ylim <- range(c(y,yN), na.rm=TRUE)


# Display raw and backtransform data
layout(matrix(1:4, nrow=2, byrow=TRUE))
par(mar=c(4,4,2,1)+0.1)
for (rr in 1:2) {
  ylab <- substitute(y[c], list=list(c=rr))
  for (cc in 2:3) {
    if (cc == rr) {
      plot.new()
      next
    }
    xlab <- substitute(y[c], list=list(c=cc))
    plot(NA, xlim=ylim, ylim=ylim, xlab=xlab, ylab=ylab)
    abline(a=0, b=1, lty=2)
    points(y[,c(cc,rr)])
    points(yN[,c(cc,rr)], col="tomato")
  }
}


layout(matrix(1:4, nrow=2, byrow=TRUE))
par(mar=c(4,4,2,1)+0.1)
for (cc in 1:3) {
  ylab <- substitute(y[c], list=list(c=cc))
  plot(NA, xlim=xlim, ylim=ylim, xlab="x", ylab=ylab)
  points(x, y[,cc])
  points(x, yN[,cc], col="tomato")
}
}

\seealso{
  \code{\link[aroma.light:backtransformPrincipalCurve.matrix]{*backtransformPrincipalCurve}()}.
  \code{\link[princurve]{principal.curve}}.
}
\keyword{methods}
