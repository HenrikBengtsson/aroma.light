%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitIWPCA.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitIWPCA.matrix}
\alias{fitIWPCA.matrix}
\alias{fitIWPCA.matrix}


\title{Robust fit of linear subspace through multidimensional data}

\description{
  Robust fit of linear subspace through multidimensional data.
}

\usage{\method{fitIWPCA}{matrix}(X, constraint=c("diagonal", "baseline", "max"), baselineChannel=NULL, ..., aShift=rep(0, ncol(X)), Xmin=NULL)}

\arguments{
 \item{X}{NxK \code{\link[base]{matrix}} where N is the number of observations and
          K is the number of dimensions (channels).
 }

 \item{constraint}{A \code{\link[base]{character}} string or a \code{\link[base]{numeric}} value.
  If \code{\link[base]{character}} it specifies which additional contraint to be used
  to specify the offset parameters along the fitted line;

  If \code{"diagonal"}, the offset vector will be a point on the line
  that is closest to the diagonal line (1,...,1).
  With this constraint, all bias parameters are identifiable.

  If \code{"baseline"} (requires argument \code{baselineChannel}), the
  estimates are such that of the bias and scale parameters of the
  baseline channel is 0 and 1, respectively.
  With this constraint, all bias parameters are identifiable.

  If \code{"max"}, the offset vector will the point on the line that is
  as "great" as possible, but still such that each of its components is
  less than the corresponding minimal signal. This will guarantee that
  no negative signals are created in the backward transformation.
  If \code{\link[base]{numeric}} value, the offset vector will the point on the line
  such that after applying the backward transformation there are
  \code{constraint*N}. Note that \code{constraint==0} corresponds
  approximately to \code{constraint=="max"}.
  With the latter two constraints, the bias parameters are only
  identifiable modulo the fitted line.
 }

 \item{baselineChannel}{Index of channel toward which all other
   channels are conform.
   This argument is required if \code{constraint=="baseline"}.
   This argument is optional if \code{constraint=="diagonal"} and
   then the scale factor of the baseline channel will be one. The
   estimate of the bias parameters is not affected in this case.
   Defaults to one, if missing.
 }

 \item{...}{Additional arguments accepted by \code{\link[aroma.light:iwpca.matrix]{*iwpca}()}.
  For instance, a N \code{\link[base]{vector}} of weights for each observation may be
   given, otherwise they get the same weight.
 }

 \item{aShift, Xmin}{For internal use only.}
}

\value{
  Returns a \code{\link[base]{list}} that contains estimated parameters and algorithm
  details;

  \item{a}{A \code{\link[base]{double}} \code{\link[base]{vector}} \eqn{(a[1],...,a[K])}with offset
      parameter estimates.
      It is made identifiable according to argument \code{constraint}.
  }
  \item{b}{A \code{\link[base]{double}} \code{\link[base]{vector}} \eqn{(b[1],...,b[K])}with scale
      parameter estimates.  It is made identifiable by constraining
      \code{b[baselineChannel] == 1}.
      These estimates are idependent of argument \code{constraint}.
  }
  \item{adiag}{If identifiability constraint \code{"diagonal"},
      a \code{\link[base]{double}} \code{\link[base]{vector}} \eqn{(adiag[1],...,adiag[K])}, where
      \eqn{adiag[1] = adiag[2] = ... adiag[K]}, specifying the point
      on the diagonal line that is closest to the fitted line,
      otherwise the zero vector.
  }
  \item{eigen}{A KxK \code{\link[base]{matrix}} with columns of eigenvectors.
  }
  \item{converged}{\code{\link[base:logical]{TRUE}} if the algorithm converged, otherwise \code{\link[base:logical]{FALSE}}.
  }
  \item{nbrOfIterations}{The number of iterations for the algorithm
                          to converge, or zero if it did not converge.
  }

  \item{t0}{Internal parameter estimates, which contains no more
                           information than the above listed elements.
  }
  \item{t}{Always \code{\link[base]{NULL}}.}
}

\details{
  This method uses re-weighted principal component analysis (IWPCA)
  to fit a the nodel \eqn{y_n = a + bx_n + eps_n} where \eqn{y_n},
  \eqn{a}, \eqn{b}, and \eqn{eps_n} are vector of the K and \eqn{x_n}
  is a scalar.

  The algorithm is:
   For iteration i:
   1) Fit a line \eqn{L} through the data close using weighted PCA
      with weights \eqn{\{w_n\}}. Let
        \eqn{r_n = \{r_{n,1},...,r_{n,K}\}}
      be the \eqn{K} principal components.
   2) Update the weights as
        \eqn{w_n <- 1 / \sum_{2}^{K} (r_{n,k} + \epsilon_r)}
      where we have used the residuals of all but the first principal
      component.
   3) Find the point a on \eqn{L} that is closest to the
      line \eqn{D=(1,1,...,1)}. Similarily, denote the point on D that is
      closest to \eqn{L} by \eqn{t=a*(1,1,...,1)}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

%examples "fitMultiIWPCA.matrix.Rex"

\seealso{
  This is an internal method used by the \code{\link[aroma.light:calibrateMultiscan.matrix]{*calibrateMultiscan}()}
  and \code{\link[aroma.light:normalizeAffine.matrix]{*normalizeAffine}()} methods.
  Internally the function \code{\link[aroma.light:iwpca.matrix]{*iwpca}()} is used to fit a line
  through the data cloud and the function \code{\link{distanceBetweenLines}}() to
  find the closest point to the diagonal (1,1,...,1).
}


\keyword{methods}
\keyword{algebra}
