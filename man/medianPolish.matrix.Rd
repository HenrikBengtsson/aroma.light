%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  medianPolish.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{medianPolish.matrix}
\alias{medianPolish.matrix}

\title{Median polish}

\description{
  Median polish.
}

\usage{
\method{medianPolish}{matrix}(X, tol=0.01, maxIter=10, na.rm=NA, ..., .addExtra=TRUE)
}

\arguments{
 \item{X}{N-times-K \code{\link[base]{matrix}}}
 \item{tol}{A \code{\link[base]{numeric}} value greater than zero used as a threshold
    to identify when the algorithm has converged.}
 \item{maxIter}{Maximum number of iterations.}
 \item{na.rm}{If \code{\link[base:logical]{TRUE}} (\code{\link[base:logical]{FALSE}}), \code{\link[base]{NA}}s are exclude (not exclude).
    If \code{\link[base]{NA}}, it is assumed that \code{X} contains no \code{\link[base]{NA}} values.}
 \item{.addExtra}{If \code{\link[base:logical]{TRUE}}, the name of argument \code{X} is returned
    and the returned structure is assigned a class.  This will make the
    result compatible what \code{\link[stats]{medpolish}} returns.}
 \item{...}{Not used.}
}

\value{
  Returns a named \code{\link[base]{list}} structure with elements:
   \item{overall}{The fitted constant term.}
   \item{row}{The fitted row effect.}
   \item{col}{The fitted column effect.}
   \item{residuals}{The residuals.}
   \item{converged}{If \code{\link[base:logical]{TRUE}}, the algorithm converged, otherwise not.}
}

\details{
  The implementation of this method give identical estimates as
  \code{\link[stats]{medpolish}}, but is about 3-5 times more efficient when
  there are no \code{\link[base]{NA}} values.
}

\author{Henrik Bengtsson}

\examples{

# Deaths from sport parachuting;  from ABC of EDA, p.224:
deaths <- matrix(c(14,15,14, 7,4,7, 8,2,10, 15,9,10, 0,2,0), ncol=3, byrow=TRUE)
rownames(deaths) <- c("1-24", "25-74", "75-199", "200++", "NA")
colnames(deaths) <- 1973:1975

print(deaths)

mp <- medianPolish(deaths)
mp1 <- medpolish(deaths, trace=FALSE)
print(mp)

ff <- c("overall", "row", "col", "residuals")
stopifnot(all.equal(mp[ff], mp1[ff]))

# Validate decomposition:
stopifnot(all.equal(deaths, mp$overall+outer(mp$row,mp$col,"+")+mp$resid))
}

\seealso{
  \code{\link[stats]{medpolish}}.
}


\keyword{methods}
\keyword{algebra}
