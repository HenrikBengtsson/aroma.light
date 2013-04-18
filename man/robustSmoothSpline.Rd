%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  robustSmoothSpline.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{robustSmoothSpline}
\alias{robustSmoothSpline.default}
\alias{robustSmoothSpline}

\title{Robust fit of a Smoothing Spline}

\usage{
\method{robustSmoothSpline}{default}(x, y=NULL, w=NULL, ..., minIter=3, maxIter=max(minIter, 50), sdCriteria=2e-04, reps=1e-15, tol=1e-06 * IQR(x), plotCurves=FALSE)
}

\description{
  Fits a smoothing spline robustly using the \eqn{L_1} norm. Currently, the
  algorithm is an \emph{iterative reweighted smooth spline} algorithm which
  calls \code{smooth.spline(x,y,w,...)} at each iteration with the weights
  \code{w} equal to the inverse of the absolute value of the residuals for
  the last iteration step.
}

\arguments{
  \item{x}{a \code{\link[base]{vector}} giving the values of the predictor variable, or a
           \code{\link[base]{list}} or a two-column \code{\link[base]{matrix}} specifying \code{x} and \code{y}.
           If \code{x} is of class \code{smooth.spline} then \code{x$x} is used
           as the \code{x} values and \code{x$yin} are used as the \code{y}
           values.}
  \item{y}{responses. If \code{y} is missing, the responses are assumed to be
           specified by \code{x}.}
  \item{w}{a \code{\link[base]{vector}} of weights the same length as \code{x} giving the weights
           to use for each element of \code{x}. Default value is equal weight
           to all values.}
  \item{...}{Other arguments passed to \code{\link[stats]{smooth.spline}}.}
  \item{minIter}{the minimum number of iterations used to fit the smoothing
           spline robustly. Default value is 3.}
  \item{maxIter}{the maximum number of iterations used to fit the smoothing
           spline robustly. Default value is 25.}
  \item{sdCriteria}{Convergence criteria, which the difference between the
           standard deviation of the residuals between two consecutive iteration
           steps. Default value is 2e-4.}
  \item{reps}{Small positive number added to residuals to avoid division by
           zero when calculating new weights for next iteration.}
  \item{tol}{Passed to \code{\link[stats]{smooth.spline}} (R >= 2.14.0).}
  \item{plotCurves}{If \code{\link[base:logical]{TRUE}}, the fitted splines are added to the current
        plot, otherwise not.}
}

\value{
  Returns an object of class \code{smooth.spline}.
}

\examples{
data(cars)
attach(cars)
plot(speed, dist, main = "data(cars)  &  robust smoothing splines")

# Fit a smoothing spline using L_2 norm
cars.spl <- smooth.spline(speed, dist)
lines(cars.spl, col = "blue")

# Fit a smoothing spline using L_1 norm
cars.rspl <- robustSmoothSpline(speed, dist)
lines(cars.rspl, col = "red")

# Fit a smoothing spline using L_2 norm with 10 degrees of freedom
lines(smooth.spline(speed, dist, df=10), lty=2, col = "blue")

# Fit a smoothing spline using L_1 norm with 10 degrees of freedom
lines(robustSmoothSpline(speed, dist, df=10), lty=2, col = "red")

legend(5,120, c(
    paste("smooth.spline [C.V.] => df =",round(cars.spl$df,1)),
    paste("robustSmoothSpline [C.V.] => df =",round(cars.rspl$df,1)),
    "standard with s( * , df = 10)", "robust with s( * , df = 10)"
  ), col = c("blue","red","blue","red"), lty = c(1,1,2,2), bg='bisque')
}

\seealso{
  This implementation of this function was adopted from
  \code{\link[stats]{smooth.spline}} of the \pkg{stats} package.
  Because of this, this function is also licensed under GPL~v2.
}

\author{Henrik Bengtsson}



\keyword{smooth}
\keyword{robust}
