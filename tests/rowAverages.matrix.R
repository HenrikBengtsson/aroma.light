library("aroma.light")

X <- matrix(1:30, nrow=5L, ncol=6L)
mu <- rowMeans(X)
sd <- apply(X, MARGIN=1L, FUN=sd)

y <- rowAverages(X)
stopifnot(all(y == mu))
stopifnot(all(attr(y,"deviance") == sd))
stopifnot(all(attr(y,"df") == ncol(X)))
