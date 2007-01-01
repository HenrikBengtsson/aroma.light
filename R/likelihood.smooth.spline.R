###########################################################################/**
# @class smooth.spline
# @RdocMethod likelihood
#
# @title "Calculate the log likelihood of a smoothing spline given the data"
#
# @synopsis
#
# \arguments{
#   \item{object}{The smooth.spline object.}
#   \item{x, y}{The x and y values for which the (weighted) likelihood will
#   be calculated. If \code{x} is of type \code{xy.coords} any value of
#   argument \code{y} will be omitted. If \code{x==NULL}, the x and y values
#   of the smoothing spline will be used.}
#   \item{w}{The weights for which the (weighted) likelihood will be
#   calculated. If @NULL, weights equal to one are assumed.}
#   \item{base}{The base of the logarithm of the likelihood. If @NULL,
#     the non-logged likelihood is returned.}
#   \item{rel.tol}{The relative tolerance used in the call to
#     \code{integrate}.}
#   \item{...}{Not used.}
# }
#
# \description{
#  Calculate the (log) likelihood of a spline given the data used to fit 
#  the spline, \eqn{g}. The likelihood consists of two main parts:
#  1) (weighted) residuals sum of squares, and 2) a penalty term. The
#  penalty term consists of a \emph{smoothing parameter} \eqn{lambda}
#  and a \emph{roughness measure} of the spline
#  \eqn{J(g) = \int g''(t) dt}. Hence, the overall log likelihood is
#    \deqn{\log L(g|x) = (y-g(x))'W(y-g(x)) + \lambda J(g)}
#  In addition to the overall likelihood, all its seperate
#  components are also returned.
#
#  Note: when fitting a smooth spline with \eqn{(x,y)} values where the
#  \eqn{x}'s are \emph{not} unique, \code{smooth.spline} will replace
#  such \eqn{(x,y)}'s with a new pair \eqn{(x,y')} where \eqn{y'} is a
#  reweighted average on the original \eqn{y}'s. It is important to
#  be aware of this. In such cases, the resulting \code{smooth.spline}
#  object does \emph{not} contain all \eqn{(x,y)}'s and therefore this
#  function will not calculate the weighted residuals sum of square on
#  the original data set, but on the data set with unique \eqn{x}'s.
#  See examples below how to calculate the likelihood for the spline with
#  the original data.
# }
#
# \value{
#   Returns the overall (log) likelihood of class
#   \code{SmoothSplineLikelihood}, a class with the following attributes:
#    \item{wrss}{the (weighted) residual sum of square}
#    \item{penalty}{the penalty which is equal to \code{-lambda*roughness}.}
#    \item{lambda}{the smoothing parameter}
#    \item{roughness}{the value of the roughness functional given the
#       specific smoothing spline and the range of data}
# }
#
# \details{
#   The roughness penalty for the smoothing spline, \eqn{g}, fitted
#   from data in the interval \eqn{[a,b]} is defined as
#       \deqn{J(g) = \int_a^b g''(t) dt}
#   which is the same as 
#       \deqn{J(g) = g'(b) - g'(a)}
#   The latter is calculated internally by using
#   @see "stats::predict.smooth.spline".
# }
#
# @examples "../incl/likelihood.smooth.spline.Rex"
#
# \seealso{
#   @see "stats::smooth.spline" and @see "robustSmoothSpline".
# }
#
# @author
#
# @keyword "smooth"
#*/###########################################################################
setMethodS3("likelihood", "smooth.spline", function(object, x=NULL, y=NULL, w=NULL, base=exp(1), rel.tol=.Machine$double.eps^(1/8), ...) {
  ## require(stats); # smooth.spline()  # Loaded by default

  g <- object;

  if (is.null(x)) {
    x  <- g$x;
    y  <- g$yin;
    w  <- g$w;
    yg <- g$y;
  } else {
    xy <- xy.coords(x, y);
    if (is.element("w", names(x)))
      w <- x$w;
    x  <- xy$x;
    y  <- xy$y;
    if (is.null(w))
      w <- rep(1, length(x));
    yg <- NULL;
    ok <- (!is.na(x) & !is.na(y) & !is.na(w));
    if (any(ok == FALSE)) {
      x <- x[ok];
      y <- y[ok];
      z <- z[ok];
    }
  }

  # The weighted residuals sum of square
  if (is.null(yg))
    yg <- predict(g, x)$y;
  wrss <- sum(w * (y-yg)^2);
  
  # The smoothing parameter
  lambda <- g$lambda

  # The roughness score J(g) = \int_a^b (g''(t))^2 dt
  gbis <- smooth.spline(predict(g, x, deriv=2));
  ab <- range(x, na.rm=TRUE);
  Jg <- integrate(function(x) predict(gbis, x=x)$y,
                  lower=ab[1], upper=ab[2],
                  rel.tol=rel.tol, stop.on.error=FALSE)$value
  
  # The penalty term
  penalty <- -lambda * Jg;
  
  # The log likelihood
  l <- -(wrss + penalty);

  # Return the correct logarithm (if any)
  if (missing(base) || (!is.null(base) && base == exp(1))) {
  } else if (is.null(base)) {
    l <- exp(l);
  } else {
    l <- l*log(exp(1), base=base);;
  }

  attr(l, "base") <- base;
  attr(l, "wrss") <- wrss;
  attr(l, "lambda") <- lambda;
  attr(l, "roughness") <- Jg;
  attr(l, "penalty") <- penalty;
  class(l) <- "SmoothSplineLikelihood";
  l;
})




############################################################################
# HISTORY:
# 2007-01-01
# o Removed any code to make method backward compatibility with 
#   R < 1.9.0, which was before 'modreg' was merged into 'stats'.
# 2005-06-03
# o now returns an object of class SmoothSplineLikelihood.
# 2005-02-20
# o Now using setMethodS3() and added '...' to please R CMD check.
# 2002-04-21
# o Updated due to modreg is merged into stats from R v1.9.0.
# 2002-03-04
# o Returns an object of class likelihood. smooth.spline instead of a list.
# o Added the option to explicitly specify x, y and w.
# o Rename from logLikelihood(...) to likelihood(..., base=exp(1)).
# o BUG FIX: Forgot to take the square in the integral of J(g).
# 2002-03-02
# o Added Rdoc details about case with non unique x values.
# 2002-03-01
# o Wrote the Rdoc comments
# o Created.
############################################################################



