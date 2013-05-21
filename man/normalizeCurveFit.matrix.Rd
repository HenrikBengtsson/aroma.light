%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeCurveFit.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeCurveFit.matrix}
\alias{normalizeCurveFit.matrix}
\alias{normalizeLoess.matrix}
\alias{normalizeLowess.matrix}
\alias{normalizeSpline.matrix}
\alias{normalizeRobustSpline.matrix}

\encoding{latin1}

\title{Weighted curve-fit normalization between a pair of channels}

\description{
  Weighted curve-fit normalization between a pair of channels.

  This method will estimate a smooth function of the dependency
  between the log-ratios and the log-intensity of the two channels and
  then correct the log-ratios (only) in order to remove the dependency.
  This is method is also known as \emph{intensity-dependent} or
  \emph{lowess normalization}.

  The curve-fit methods are by nature limited to paired-channel data.
  There exist at least one method trying to overcome this limitation,
  namely the cyclic-lowess [1], which applies the paired
  curve-fit method iteratively over all pairs of channels/arrays.
  Cyclic-lowess is not implented here.

  We recommend that affine normalization [2] is used instead of curve-fit
  normalization.
}

\usage{
\method{normalizeCurveFit}{matrix}(X, weights=NULL, typeOfWeights=c("datapoint"),
  method=c("loess", "lowess", "spline", "robustSpline"), bandwidth=NULL,
  satSignal=2^16 - 1, ...)
}

\arguments{
 \item{X}{An Nx2 \code{\link[base]{matrix}} where the columns represent the two channels
   to be normalized.}
 \item{weights}{If \code{\link[base]{NULL}}, non-weighted normalization is done.
   If data-point weights are used, this should be a \code{\link[base]{vector}} of length
   N of data point weights used when estimating the normalization
   function.
 }
 \item{typeOfWeights}{A \code{\link[base]{character}} string specifying the type of
   weights given in argument \code{weights}.
 }
 \item{method}{\code{\link[base]{character}} string specifying which method to use when
   fitting the intensity-dependent function.
   Supported methods:
    \code{"loess"} (better than lowess),
    \code{"lowess"} (classic; supports only zero-one weights),
    \code{"spline"} (more robust than lowess at lower and upper
                     intensities; supports only zero-one weights),
    \code{"robustSpline"} (better than spline).
 }
 \item{bandwidth}{A \code{\link[base]{double}} value specifying the bandwidth of the
  estimator used.
 }
 \item{satSignal}{Signals equal to or above this threshold will not
   be used in the fitting.
 }
 \item{...}{Not used.}
}

\value{
  A Nx2 \code{\link[base]{matrix}} of the normalized two channels.
  The fitted model is returned as attribute \code{modelFit}.
}

\details{
 A smooth function \eqn{c(A)} is fitted throught data in \eqn{(A,M)},
 where \eqn{M=log_2(y_2/y_1)} and \eqn{A=1/2*log_2(y_2*y_1)}. Data is
 normalized by \eqn{M <- M - c(A)}.

 Loess is by far the slowest method of the four, then lowess, and then
 robust spline, which iteratively calls the spline method.
}

\section{Negative, non-positive, and saturated values}{
 Non-positive values are set to not-a-number (\code{\link[base:is.finite]{NaN}}).
 Data points that are saturated in one or more channels are not used
 to estimate the normalization function, but they are normalized.
}

\section{Missing values}{
 The estimation of the normalization function will only be made
 based on complete non-saturated observations, i.e. observations that
 contains no \code{\link[base]{NA}} values nor saturated values as defined by \code{satSignal}.
}

\section{Weighted normalization}{
 Each data point, that is, each row in \code{X}, which is a
 vector of length 2, can be assigned a weight in [0,1] specifying how much
 it should \emph{affect the fitting of the normalization function}.
 Weights are given by argument \code{weights}, which should be a \code{\link[base]{numeric}}
 \code{\link[base]{vector}} of length N. Regardless of weights, all data points are
 \emph{normalized} based on the fitted normalization function.

 Note that the lowess and the spline method only support zero-one
 \{0,1\} weights.
 For such methods, all weights that are less than a half are set to zero.
}

\section{Details on loess}{
 For \code{\link[stats]{loess}}, the arguments \code{family="symmetric"},
 \code{degree=1}, \code{span=3/4},
 \code{control=loess.control(trace.hat="approximate"},
 \code{iterations=5}, \code{surface="direct")} are used.
}

\author{Henrik Bengtsson}

\references{
  [1] M. \enc{Åstrand}{Astrand},
      Contrast Normalization of Oligonucleotide Arrays,
      Journal Computational Biology, 2003, 10, 95-102. \cr
  [2] Henrik Bengtsson and Ola \enc{Hössjer}{Hossjer}, \emph{Methodological Study of Affine Transformations of Gene Expression Data}, Methodological study of affine transformations of gene expression data with proposed robust non-parametric multi-dimensional normalization method, BMC Bioinformatics, 2006, 7:100.
 \cr
}

\examples{
 pathname <- system.file("data-ex", "PMT-RGData.dat", package="aroma.light")
rg <- read.table(pathname, header=TRUE, sep="\t")
nbrOfScans <- max(rg$slide)

rg <- as.list(rg)
for (field in c("R", "G"))
  rg[[field]] <- matrix(as.double(rg[[field]]), ncol=nbrOfScans)
rg$slide <- rg$spot <- NULL
rg <- as.matrix(as.data.frame(rg))
colnames(rg) <- rep(c("R", "G"), each=nbrOfScans)

layout(matrix(c(1,2,0,3,4,0,5,6,7), ncol=3, byrow=TRUE))

rgC <- rg
for (channel in c("R", "G")) {
  sidx <- which(colnames(rg) == channel)
  channelColor <- switch(channel, R="red", G="green");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The raw data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  plotMvsAPairs(rg[,sidx])
  title(main=paste("Observed", channel))
  box(col=channelColor)
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The calibrated data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rgC[,sidx] <- calibrateMultiscan(rg[,sidx], average=NULL)

  plotMvsAPairs(rgC[,sidx])
  title(main=paste("Calibrated", channel))
  box(col=channelColor)
} # for (channel ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The average calibrated data
#
# Note how the red signals are weaker than the green. The reason
# for this can be that the scale factor in the green channel is
# greater than in the red channel, but it can also be that there
# is a remaining relative difference in bias between the green
# and the red channel, a bias that precedes the scanning.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rgCA <- rg
for (channel in c("R", "G")) {
  sidx <- which(colnames(rg) == channel)
  rgCA[,sidx] <- calibrateMultiscan(rg[,sidx])
}

rgCAavg <- matrix(NA, nrow=nrow(rgCA), ncol=2)
colnames(rgCAavg) <- c("R", "G");
for (channel in c("R", "G")) {
  sidx <- which(colnames(rg) == channel)
  rgCAavg[,channel] <- apply(rgCA[,sidx], MARGIN=1, FUN=median, na.rm=TRUE);
}

# Add some "fake" outliers
outliers <- 1:600
rgCAavg[outliers,"G"] <- 50000;

plotMvsA(rgCAavg)
title(main="Average calibrated (AC)")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Weight-down outliers when normalizing
weights <- rep(1, nrow(rgCAavg));
weights[outliers] <- 0.001;

# Affine normalization of channels
rgCANa <- normalizeAffine(rgCAavg, weights=weights)
# It is always ok to rescale the affine normalized data if its
# done on (R,G); not on (A,M)! However, this is only needed for
# esthetic purposes.
rgCANa <- rgCANa *2^1.4
plotMvsA(rgCANa)
title(main="Normalized AC")

# Curve-fit (lowess) normalization
rgCANlw <- normalizeLowess(rgCAavg, weights=weights)
plotMvsA(rgCANlw, col="orange", add=TRUE)

# Curve-fit (loess) normalization
rgCANl <- normalizeLoess(rgCAavg, weights=weights)
plotMvsA(rgCANl, col="red", add=TRUE)

# Curve-fit (robust spline) normalization
rgCANrs <- normalizeRobustSpline(rgCAavg, weights=weights)
plotMvsA(rgCANrs, col="blue", add=TRUE)

legend(x=0,y=16, legend=c("affine", "lowess", "loess", "r. spline"), pch=19,
       col=c("black", "orange", "red", "blue"), ncol=2, x.intersp=0.3, bty="n")


plotMvsMPairs(cbind(rgCANa, rgCANlw), col="orange", xlab=expression(M[affine]))
title(main="Normalized AC")
plotMvsMPairs(cbind(rgCANa, rgCANl), col="red", add=TRUE)
plotMvsMPairs(cbind(rgCANa, rgCANrs), col="blue", add=TRUE)
abline(a=0, b=1, lty=2)
legend(x=-6,y=6, legend=c("lowess", "loess", "r. spline"), pch=19,
       col=c("orange", "red", "blue"), ncol=2, x.intersp=0.3, bty="n")


}

\seealso{
  \code{\link[aroma.light:normalizeAffine.matrix]{*normalizeAffine}()}.
}
\keyword{methods}
