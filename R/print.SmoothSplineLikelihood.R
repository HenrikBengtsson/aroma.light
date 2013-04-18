###########################################################################/**
# @class SmoothSplineLikelihood
# @RdocMethod print
#
# @title "Prints an SmoothSplineLikelihood object"
#
# \description{
#  @get "title".  A SmoothSplineLikelihood object is returned by
#  \code{\link{likelihood.smooth.spline}()}.
# }
#
# @synopsis
#
# \arguments{
#   \item{x}{Object to be printed.}
#   \item{digits}{Minimal number of significant digits to print.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author "HB"
#*/###########################################################################
setMethodS3("print", "SmoothSplineLikelihood", function(x, digits=getOption("digits"), ...) {
  # To please R CMD check...
  object <- x;

  s <- paste("Likelihood of smoothing spline:", format(object, digits=digits), "\n");
  base <- attr(object, "base");
  s <- paste(s, "Log base:", format(base, digits=digits), "\n")
  wrss <- attr(object, "wrss");
  s <- paste(s, "Weighted residuals sum of square:", format(wrss, digits=digits), "\n");
  penalty <- attr(object, "penalty");
  s <- paste(s, "Penalty:", format(penalty, digits=digits), "\n");
  lambda <- attr(object, "lambda");
  s <- paste(s, "Smoothing parameter lambda:", format(lambda, digits=digits), "\n");
  roughness <- attr(object, "roughness");
  s <- paste(s, "Roughness score:", format(roughness, digits=digits), "\n");

  cat(s);

  invisible(object);
})

############################################################################
# HISTORY:
# 2005-06-03
# o Added Rdoc comments.
# o Extracted from likelihood.smooth.spline.R.
############################################################################
