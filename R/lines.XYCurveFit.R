setMethodS3("lines", "XYCurveFit", function(x, xNew=NULL, ...) {
  # To please R CMD check
  fit <- x;

  if (is.null(xNew)) {
    xNew <- fit$x;
    xNew <- sort(xNew);
    xNew <- xNew[!duplicated(xNew)];
  }
  y <- fit$predictY(xNew);
  lines(x=xNew, y=y, ...);
}) # lines()


############################################################################
# HISTORY:
# 2009-07-15
# o Created.
############################################################################
