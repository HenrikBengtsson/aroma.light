%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  distanceBetweenLines.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{distanceBetweenLines}
\alias{distanceBetweenLines.default}
\alias{distanceBetweenLines}

\title{Finds the shortest distance between two lines}

\description{
  Finds the shortest distance between two lines.

  Consider the two lines

    \eqn{x(s) = a_x + b_x*s} and \eqn{y(t) = a_y + b_y*t}

  in an K-space where the offset and direction \code{\link[base]{vector}}s are \eqn{a_x}
  and \eqn{b_x} (in \eqn{R^K}) that define the line \eqn{x(s)}
 (\eqn{s} is a scalar). Similar for the line \eqn{y(t)}.
  This function finds the point \eqn{(s,t)} for which \eqn{|x(s)-x(t)|}
  is minimal.
}

\usage{
\method{distanceBetweenLines}{default}(ax, bx, ay, by, ...)
}

\arguments{
 \item{ax,bx}{Offset and direction \code{\link[base]{vector}} of length K for line \eqn{z_x}.}
 \item{ay,by}{Offset and direction \code{\link[base]{vector}} of length K for line \eqn{z_y}.}
 \item{...}{Not used.}
}

\value{
  Returns the a \code{\link[base]{list}} containing
  \item{ax,bx}{The given line \eqn{x(s)}.}
  \item{ay,by}{The given line \eqn{y(t)}.}
  \item{s,t}{The values of \eqn{s} and \eqn{t} such that
       \eqn{|x(s)-y(t)|} is minimal.}
  \item{xs,yt}{The values of \eqn{x(s)} and \eqn{y(t)}
       at the optimal point \eqn{(s,t)}.}
  \item{distance}{The distance between the lines, i.e. \eqn{|x(s)-y(t)|}
       at the optimal point \eqn{(s,t)}.}
}

\author{Henrik Bengtsson}

\examples{
for (zzz in 0) {

# This example requires plot3d() in R.basic [http://www.braju.com/R/]
if (!require(pkgName <- "R.basic", character.only=TRUE)) break

layout(matrix(1:4, nrow=2, ncol=2, byrow=TRUE))

############################################################
# Lines in two-dimensions
############################################################
x <- list(a=c(1,0), b=c(1,2))
y <- list(a=c(0,2), b=c(1,1))
fit <- distanceBetweenLines(ax=x$a, bx=x$b, ay=y$a, by=y$b)

xlim <- ylim <- c(-1,8)
plot(NA, xlab="", ylab="", xlim=ylim, ylim=ylim)

# Highlight the offset coordinates for both lines
points(t(x$a), pch="+", col="red")
text(t(x$a), label=expression(a[x]), adj=c(-1,0.5))
points(t(y$a), pch="+", col="blue")
text(t(y$a), label=expression(a[y]), adj=c(-1,0.5))

v <- c(-1,1)*10;
xv <- list(x=x$a[1]+x$b[1]*v, y=x$a[2]+x$b[2]*v)
yv <- list(x=y$a[1]+y$b[1]*v, y=y$a[2]+y$b[2]*v)

lines(xv, col="red")
lines(yv, col="blue")

points(t(fit$xs), cex=2.0, col="red")
text(t(fit$xs), label=expression(x(s)), adj=c(+2,0.5))
points(t(fit$yt), cex=1.5, col="blue")
text(t(fit$yt), label=expression(y(t)), adj=c(-1,0.5))
print(fit)


############################################################
# Lines in three-dimensions
############################################################
x <- list(a=c(0,0,0), b=c(1,1,1))  # The 'diagonal'
y <- list(a=c(2,1,2), b=c(2,1,3))  # A 'fitted' line
fit <- distanceBetweenLines(ax=x$a, bx=x$b, ay=y$a, by=y$b)

xlim <- ylim <- zlim <- c(-1,3)
dummy <- t(c(1,1,1))*100;

# Coordinates for the lines in 3d
v <- seq(-10,10, by=1);
xv <- list(x=x$a[1]+x$b[1]*v, y=x$a[2]+x$b[2]*v, z=x$a[3]+x$b[3]*v)
yv <- list(x=y$a[1]+y$b[1]*v, y=y$a[2]+y$b[2]*v, z=y$a[3]+y$b[3]*v)

for (theta in seq(30,140,length=3)) {
  plot3d(dummy, theta=theta, phi=30, xlab="", ylab="", zlab="",
                             xlim=ylim, ylim=ylim, zlim=zlim)

  # Highlight the offset coordinates for both lines
  points3d(t(x$a), pch="+", col="red")
  text3d(t(x$a), label=expression(a[x]), adj=c(-1,0.5))
  points3d(t(y$a), pch="+", col="blue")
  text3d(t(y$a), label=expression(a[y]), adj=c(-1,0.5))

  # Draw the lines
  lines3d(xv, col="red")
  lines3d(yv, col="blue")

  # Draw the two points that are closest to each other
  points3d(t(fit$xs), cex=2.0, col="red")
  text3d(t(fit$xs), label=expression(x(s)), adj=c(+2,0.5))
  points3d(t(fit$yt), cex=1.5, col="blue")
  text3d(t(fit$yt), label=expression(y(t)), adj=c(-1,0.5))

  # Draw the distance between the two points
  lines3d(rbind(fit$xs,fit$yt), col="purple", lwd=2)
}

print(fit)

} # for (zzz in 0)
rm(zzz)
}

\references{
 [1] M. Bard and D. Himel, \emph{The Minimum Distance Between Two
    Lines in n-Space}, September 2001, Advisor Dennis Merino.\cr
 [2] Dan Sunday, \emph{Distance between Lines and Segments with
    their Closest Point of Approach},
    \url{http://geometryalgorithms.com/Archive/algorithm_0106/}.\cr
}


\keyword{algebra}
