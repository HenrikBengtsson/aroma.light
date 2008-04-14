setMethodS3("predict", "lowess", function(object, newdata=NULL, ties=mean, ...) {
  approx(object, xout=newdata, ties=ties, ...)$y;
}, private=TRUE) # predict()


############################################################################
# HISTORY:
# 2006-11-28
# o Created.
############################################################################
