# The trace of matrix X
tr <- function(X, na.rm=FALSE) {
  sum(diag(X), na.rm=na.rm);
}

# u and v are column vectors.
scalarProduct <- function(u,v, na.rm=FALSE) {
  colSums(u*v, na.rm=na.rm);
}

# u are column vectors, v a single vector for now.
projectUontoV <- function(u,v, na.rm=FALSE) {
  vN <- v / sqrt(sum(v^2, na.rm=na.rm));
  sp <- scalarProduct(u,vN, na.rm=na.rm);
  sp <- rep(sp, each=length(vN));
  sp <- matrix(sp, nrow=length(vN));
  v * sp;
}

############################################################################
# HISTORY:
# 2005-01-24
# o Added arg. na.rm=FALSE to tr(), scalarProduct(), and projectUontoV().
# 2003-11-06
# o Created for readability.
############################################################################
