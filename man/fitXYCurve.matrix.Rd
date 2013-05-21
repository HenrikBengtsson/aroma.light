%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitXYCurve.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitXYCurve.matrix}
\alias{fitXYCurve.matrix}
\alias{fitXYCurve}

\title{Fitting a smooth curve through paired (x,y) data}

\description{
  Fitting a smooth curve through paired (x,y) data.
}

\usage{
\method{fitXYCurve}{matrix}(X, weights=NULL, typeOfWeights=c("datapoint"), method=c("loess", "lowess",
  "spline", "robustSpline"), bandwidth=NULL, satSignal=2^16 - 1, ...)
}

\arguments{
 \item{X}{An Nx2 \code{\link[base]{matrix}} where the columns represent the two channels
   to be normalized.}
 \item{weights}{If \code{\link[base]{NULL}}, non-weighted normalization is done.
   If data-point weights are used, this should be a \code{\link[base]{vector}} of length
   N of data point weights used when estimating the normalization
   function.
 }
 \item{typeOfWeights}{A \code{\link[base]{character}} string specifying the type of
   weights given in argument \code{weights}.
 }
 \item{method}{\code{\link[base]{character}} string specifying which method to use when
   fitting the intensity-dependent function.
   Supported methods:
    \code{"loess"} (better than lowess),
    \code{"lowess"} (classic; supports only zero-one weights),
    \code{"spline"} (more robust than lowess at lower and upper
                     intensities; supports only zero-one weights),
    \code{"robustSpline"} (better than spline).
 }
 \item{bandwidth}{A \code{\link[base]{double}} value specifying the bandwidth of the
  estimator used.
 }
 \item{satSignal}{Signals equal to or above this threshold will not
   be used in the fitting.
 }
 \item{...}{Not used.}
}

\value{
  A named \code{\link[base]{list}} structure of class \code{XYCurve}.
}

\section{Missing values}{
 The estimation of the function will only be made based on complete
 non-saturated observations, i.e. observations that contains no \code{\link[base]{NA}}
 values nor saturated values as defined by \code{satSignal}.
}

\section{Weighted normalization}{
 Each data point, that is, each row in \code{X}, which is a
 vector of length 2, can be assigned a weight in [0,1] specifying how much
 it should \emph{affect the fitting of the normalization function}.
 Weights are given by argument \code{weights}, which should be a \code{\link[base]{numeric}}
 \code{\link[base]{vector}} of length N.

 Note that the lowess and the spline method only support zero-one
 \{0,1\} weights.
 For such methods, all weights that are less than a half are set to zero.
}

\section{Details on loess}{
 For \code{\link[stats]{loess}}, the arguments \code{family="symmetric"},
 \code{degree=1}, \code{span=3/4},
 \code{control=loess.control(trace.hat="approximate"},
 \code{iterations=5}, \code{surface="direct")} are used.
}

\author{Henrik Bengtsson}

\examples{
 # Simulate data from the model y <- a + bx + x^c + eps(bx)
x <- rexp(1000)
a <- c(2,15)
b <- c(2,1)
c <- c(1,2)
bx <- outer(b,x)
xc <- t(sapply(c, FUN=function(c) x^c))
eps <- apply(bx, MARGIN=2, FUN=function(x) rnorm(length(x), mean=0, sd=0.1*x))
Y <- a + bx + xc + eps
Y <- t(Y)

lim <- c(0,70)
plot(Y, xlim=lim, ylim=lim)

# Fit principal curve through a subset of (y_1, y_2)
subset <- sample(nrow(Y), size=0.3*nrow(Y))
fit <- fitXYCurve(Y[subset,], bandwidth=0.2)

lines(fit, col="red", lwd=2)

# Backtransform (y_1, y_2) keeping y_1 unchanged
YN <- backtransformXYCurve(Y, fit=fit)
points(YN, col="blue")
abline(a=0, b=1, col="red", lwd=2)

}
\keyword{methods}
