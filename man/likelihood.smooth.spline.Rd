%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  likelihood.smooth.spline.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{likelihood.smooth.spline}
\alias{likelihood.smooth.spline}

\title{Calculate the log likelihood of a smoothing spline given the data}

\usage{
\method{likelihood}{smooth.spline}(object, x=NULL, y=NULL, w=NULL, base=exp(1),
  rel.tol=.Machine$double.eps^(1/8), ...)
}

\arguments{
  \item{object}{The smooth.spline object.}
  \item{x, y}{The x and y values for which the (weighted) likelihood will
  be calculated. If \code{x} is of type \code{xy.coords} any value of
  argument \code{y} will be omitted. If \code{x==NULL}, the x and y values
  of the smoothing spline will be used.}
  \item{w}{The weights for which the (weighted) likelihood will be
  calculated. If \code{\link[base]{NULL}}, weights equal to one are assumed.}
  \item{base}{The base of the logarithm of the likelihood. If \code{\link[base]{NULL}},
    the non-logged likelihood is returned.}
  \item{rel.tol}{The relative tolerance used in the call to
    \code{integrate}.}
  \item{...}{Not used.}
}

\description{
 Calculate the (log) likelihood of a spline given the data used to fit
 the spline, \eqn{g}. The likelihood consists of two main parts:
 1) (weighted) residuals sum of squares, and 2) a penalty term. The
 penalty term consists of a \emph{smoothing parameter} \eqn{lambda}
 and a \emph{roughness measure} of the spline
 \eqn{J(g) = \int g''(t) dt}. Hence, the overall log likelihood is
   \deqn{\log L(g|x) = (y-g(x))'W(y-g(x)) + \lambda J(g)}
 In addition to the overall likelihood, all its seperate
 components are also returned.

 Note: when fitting a smooth spline with \eqn{(x,y)} values where the
 \eqn{x}'s are \emph{not} unique, \code{smooth.spline} will replace
 such \eqn{(x,y)}'s with a new pair \eqn{(x,y')} where \eqn{y'} is a
 reweighted average on the original \eqn{y}'s. It is important to
 be aware of this. In such cases, the resulting \code{smooth.spline}
 object does \emph{not} contain all \eqn{(x,y)}'s and therefore this
 function will not calculate the weighted residuals sum of square on
 the original data set, but on the data set with unique \eqn{x}'s.
 See examples below how to calculate the likelihood for the spline with
 the original data.
}

\value{
  Returns the overall (log) likelihood of class
  \code{SmoothSplineLikelihood}, a class with the following attributes:
   \item{wrss}{the (weighted) residual sum of square}
   \item{penalty}{the penalty which is equal to \code{-lambda*roughness}.}
   \item{lambda}{the smoothing parameter}
   \item{roughness}{the value of the roughness functional given the
      specific smoothing spline and the range of data}
}

\details{
  The roughness penalty for the smoothing spline, \eqn{g}, fitted
  from data in the interval \eqn{[a,b]} is defined as
      \deqn{J(g) = \int_a^b g''(t) dt}
  which is the same as
      \deqn{J(g) = g'(b) - g'(a)}
  The latter is calculated internally by using
  \code{\link[stats]{predict.smooth.spline}}.
}

\examples{
# Define f(x)
f <- expression(0.1*x^4 + 1*x^3 + 2*x^2 + x + 10*sin(2*x))

# Simulate data from this function in the range [a,b]
a <- -2; b <- 5
x <- seq(a, b, length=3000)
y <- eval(f)

# Add some noise to the data
y <- y + rnorm(length(y), 0, 10)

# Plot the function and its second derivative
plot(x,y, type="l", lwd=4)

# Fit a cubic smoothing spline and plot it
g <- smooth.spline(x,y, df=16)
lines(g, col="yellow", lwd=2, lty=2)

# Calculating the (log) likelihood of the fitted spline
l <- likelihood(g)

cat("Log likelihood with unique x values:\n")
print(l)

# Note that this is not the same as the log likelihood of the
# data on the fitted spline iff the x values are non-unique
x[1:5] <- x[1]  # Non-unique x values
g <- smooth.spline(x,y, df=16)
l <- likelihood(g)

cat("\nLog likelihood of the *spline* data set:\n");
print(l)

# In cases with non unique x values one has to proceed as
# below if one want to get the log likelihood for the original
# data.
l <- likelihood(g, x=x, y=y)
cat("\nLog likelihood of the *original* data set:\n");
print(l)





}

\seealso{
  \code{\link[stats]{smooth.spline}} and \code{\link{robustSmoothSpline}}().
}

\author{Henrik Bengtsson}


\keyword{methods}
\keyword{smooth}
