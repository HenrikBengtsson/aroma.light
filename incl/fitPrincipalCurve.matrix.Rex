
# Simulate data from the model y <- a + bx + x^c + eps(bx)
J <- 1000
x <- rexp(J)
a <- c(2,15,3)
b <- c(2,3,4)
c <- c(1,2,1/2)
bx <- outer(b,x)
xc <- t(sapply(c, FUN=function(c) x^c))
eps <- apply(bx, MARGIN=2, FUN=function(x) rnorm(length(b), mean=0, sd=0.1*x))
y <- a + bx + xc + eps
y <- t(y)

# Fit principal curve through (y_1, y_2, y_3)
fit <- fitPrincipalCurve(y, verbose=TRUE)

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


# Pairwise signals vs x before and after transform
layout(matrix(1:4, nrow=2, byrow=TRUE))
par(mar=c(4,4,3,2)+0.1)
for (cc in 1:3) {
  ylab <- substitute(y[c], env=list(c=cc))
  plot(NA, xlim=xlim, ylim=ylim, xlab="x", ylab=ylab)
  abline(h=a[cc], lty=3)
  mtext(side=4, at=a[cc], sprintf("a=%g", a[cc]), 
        cex=0.8, las=2, line=0, adj=1.1, padj=-0.2)
  points(x, y[,cc])
  points(x, yN[,cc], col="tomato")
  legend("topleft", col=c("black", "tomato"), pch=19,
                    c("orignal", "transformed"), bty="n")
}
title(main="Pairwise signals vs x before and after transform", outer=TRUE, line=-2)


# Pairwise signals before and after transform
layout(matrix(1:4, nrow=2, byrow=TRUE))
par(mar=c(4,4,3,2)+0.1)
for (rr in 3:2) {
  ylab <- substitute(y[c], env=list(c=rr))
  for (cc in 1:2) {
    if (cc == rr) { 
      plot.new()
      next
    }
    xlab <- substitute(y[c], env=list(c=cc))
    plot(NA, xlim=ylim, ylim=ylim, xlab=xlab, ylab=ylab)
    abline(a=0, b=1, lty=2)
    points(y[,c(cc,rr)])
    points(yN[,c(cc,rr)], col="tomato")
    legend("topleft", col=c("black", "tomato"), pch=19,
                      c("orignal", "transformed"), bty="n")
  }
}
title(main="Pairwise signals before and after transform", outer=TRUE, line=-2)
