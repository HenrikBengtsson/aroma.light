#########################################################################/**
# @RdocDefault distanceBetweenLines
#
# @title "Finds the shortest distance between two lines"
#
# \description{
#   @get "title".
#
#   Consider the two lines
#
#     \eqn{x(s) = a_x + b_x*s} and \eqn{y(t) = a_y + b_y*t}
#
#   in an K-space where the offset and direction @vectors are \eqn{a_x}
#   and \eqn{b_x} (in \eqn{R^K}) that define the line \eqn{x(s)}
#  (\eqn{s} is a scalar). Similar for the line \eqn{y(t)}.
#   This function finds the point \eqn{(s,t)} for which \eqn{|x(s)-x(t)|}
#   is minimal.
# }
#
# @synopsis
#
# \arguments{
#  \item{ax,bx}{Offset and direction @vector of length K for line \eqn{z_x}.}
#  \item{ay,by}{Offset and direction @vector of length K for line \eqn{z_y}.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns the a @list containing
#   \item{ax,bx}{The given line \eqn{x(s)}.}
#   \item{ay,by}{The given line \eqn{y(t)}.}
#   \item{s,t}{The values of \eqn{s} and \eqn{t} such that
#        \eqn{|x(s)-y(t)|} is minimal.}
#   \item{xs,yt}{The values of \eqn{x(s)} and \eqn{y(t)}
#        at the optimal point \eqn{(s,t)}.}
#   \item{distance}{The distance between the lines, i.e. \eqn{|x(s)-y(t)|}
#        at the optimal point \eqn{(s,t)}.}
# }
#
# @author "HB"
#
# @examples "../incl/distanceBetweenLines.Rex"
#
# \references{
#  [1] M. Bard and D. Himel, \emph{The Minimum Distance Between Two
#     Lines in n-Space}, September 2001, Advisor Dennis Merino.\cr
#  [2] Dan Sunday, \emph{Distance between Lines and Segments with
#     their Closest Point of Approach},
#     \url{http://geometryalgorithms.com/Archive/algorithm_0106/}.\cr
# }
#
# @keyword "algebra"
#*/#########################################################################
setMethodS3("distanceBetweenLines", "default", function(ax, bx, ay, by, ...) {
  if (length(ax) != length(bx)) {
    stop(sprintf("The length of the offset vector 'ax' (%d) and direction vector 'bx' (%d) are not equal.", length(ax), length(bx)));
  }

  if (length(ay) != length(by)) {
    stop(sprintf("The length of the offset vector 'ay' (%d) and direction vector 'by' (%d) are not equal.", length(ay), length(by)));
  }

  if (length(ax) != length(ay)) {
    stop(sprintf("The line x(s) and y(t) are of different dimensions: %d vs %d", length(ax), length(ay)));
  }

  if (length(ax) <= 1)
    stop(sprintf("The lines must be in two or more dimensions: %d", length(ax)));

  ax <- as.vector(ax);
  bx <- as.vector(bx);
  ay <- as.vector(ay);
  by <- as.vector(by);

  if (length(ax) == 2) {
    # Find (s,t) such that x(s) == y(t) where
    #   x(s) = a_x + b_x*s
    #   y(t) = a_y + b_y*t
    e <- (ax-ay);
    f <- e/by;
    g <- bx/by;
    s <- (f[2]-f[1])/(g[1]-g[2]);
    t <- f[1] + g[1]*s;
    d <- 0;
  } else {
    # Consider the two lines in an K-space
    #   x(s) = a_x + b_x*t    (line 1)
    #   y(t) = a_y + b_y*s    (line 2)
    # where s and t are scalars and the other vectors in R^K.

    # Some auxillary calculations
    A <- sum(bx*bx);
    B <- 2*(sum(bx*ax)-sum(bx*ay));
    C <- 2*sum(bx*by);
    D <- 2*(sum(by*ay)-sum(by*ax));
    E <- sum(by*by);
    F <- sum(ax*ax) + sum(ay*ay);

    # Shortest distance between the two lines (points)
    G <- C^2-4*A*E;
    d2 <- (B*C*D+B^2*E+C^2*F+A*(D^2-4*E*F))/G;
    d <- sqrt(d2);

    # The points that are closest to each other.
    t <- (2*A*D+B*C)/G;   # t is on y(t)
    s <- (C*t-B)/(2*A);   # s is on x(s)
  }

  # Get the coordinates of the two points on x(s) and y(t) that
  # are closest to each other.
  xs <- ax + bx*s;
  yt <- ay + by*t;

  list(ax=ax, bx=bx, ay=ay, by=by, s=s, t=t, xs=xs, yt=yt, distance=d);
}) # distanceBetweenLines()


############################################################################
# HISTORY:
# 2005-06-03
# o Made into a default method.
# 2003-12-29
# o Added Rdoc comments.
# o Created by generalizing from formet RGData$fitIWPCA() in com.braju.smax.
############################################################################
