%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeFragmentLength.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeFragmentLength}
\alias{normalizeFragmentLength.default}
\alias{normalizeFragmentLength}

\title{Normalizes signals for PCR fragment-length effects}

\description{
 Normalizes signals for PCR fragment-length effects. Some or all signals are used to estimated the
 normalization function.  All signals are normalized.
}

\usage{
\method{normalizeFragmentLength}{default}(y, fragmentLengths, targetFcns=NULL, subsetToFit=NULL,
  onMissing=c("ignore", "median"), .isLogged=TRUE, ..., .returnFit=FALSE)
}

\arguments{
  \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length K of signals to be normalized
    across E enzymes.}
  \item{fragmentLengths}{An \code{\link[base]{integer}} KxE \code{\link[base]{matrix}} of fragment lengths.}
  \item{targetFcns}{An optional \code{\link[base]{list}} of E \code{\link[base]{function}}s; one per enzyme.
    If \code{\link[base]{NULL}}, the data is normalized to have constant fragment-length
    effects (all equal to zero on the log-scale).}
  \item{subsetToFit}{The subset of data points used to fit the
     normalization function.
     If \code{\link[base]{NULL}}, all data points are considered.}
  \item{onMissing}{Specifies how data points for which there is no
     fragment length is normalized.
     If \code{"ignore"}, the values are not modified.
     If \code{"median"}, the values are updated to have the same
     robust average as the other data points.
  }
  \item{.isLogged}{A \code{\link[base]{logical}}.}
  \item{...}{Additional arguments passed to \code{\link[stats]{lowess}}.}
  \item{.returnFit}{A \code{\link[base]{logical}}.}
}

\value{
  Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of the normalized signals.
}

\section{Multi-enzyme normalization}{
 It is assumed that the fragment-length effects from multiple enzymes
 added (with equal weights) on the intensity scale.
 The fragment-length effects are fitted for each enzyme separately based
 on units that are exclusively for that enzyme.
 \emph{If there are no or very such units for an enzyme, the assumptions
 of the model are not met and the fit will fail with an error.}
 Then, from the above single-enzyme fits the average effect across
 enzymes is the calculated for each unit that is on multiple enzymes.
}

\section{Target functions}{
  It is possible to specify custom target function effects for each
  enzyme via argument \code{targetFcns}.  This argument has to be a
  \code{\link[base]{list}} containing one \code{\link[base]{function}} per enzyme and ordered in the same
  order as the enzyme are in the columns of argument
  \code{fragmentLengths}.
  For instance, if one wish to normalize the signals such that their
  mean signal as a function of fragment length effect is contantly
  equal to 2200 (or the intensity scale), the use
  \code{targetFcns=function(fl, ...) log2(2200)} which completely
  ignores fragment-length argument 'fl' and always returns a
  constant.
  If two enzymes are used, then use
  \code{targetFcns=rep(list(function(fl, ...) log2(2200)), 2)}.

  Note, if \code{targetFcns} is \code{\link[base]{NULL}}, this corresponds to
  \code{targetFcns=rep(list(function(fl, ...) 0), ncol(fragmentLengths))}.

  Alternatively, if one wants to only apply minimial corrections to
  the signals, then one can normalize toward target functions that
  correspond to the fragment-length effect of the average array.
}

\examples{
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example 1: Single-enzyme fragment-length normalization of 6 arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number samples
I <- 9;

# Number of loci
J <- 1000;

# Fragment lengths
fl <- seq(from=100, to=1000, length.out=J);

# Simulate data points with unknown fragment lengths
hasUnknownFL <- seq(from=1, to=J, by=50);
fl[hasUnknownFL] <- NA;

# Simulate data
y <- matrix(0, nrow=J, ncol=I);
maxY <- 12;
for (kk in 1:I) {
  k <- runif(n=1, min=3, max=5);
  mu <- function(fl) {
    mu <- rep(maxY, length(fl));
    ok <- !is.na(fl);
    mu[ok] <- mu[ok] - fl[ok]^{1/k};
    mu;
  }
  eps <- rnorm(J, mean=0, sd=1);
  y[,kk] <- mu(fl) + eps;
}

# Normalize data (to a zero baseline)
yN <- apply(y, MARGIN=2, FUN=function(y) {
  normalizeFragmentLength(y, fragmentLengths=fl, onMissing="median");
})

# The correction factors
rho <- y-yN;
print(summary(rho));
# The correction for units with unknown fragment lengths
# equals the median correction factor of all other units
print(summary(rho[hasUnknownFL,]));

# Plot raw data
layout(matrix(1:9, ncol=3));
xlim <- c(0,max(fl, na.rm=TRUE));
ylim <- c(0,max(y, na.rm=TRUE));
xlab <- "Fragment length";
ylab <- expression(log2(theta));
for (kk in 1:I) {
  plot(fl, y[,kk], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
  ok <- (is.finite(fl) & is.finite(y[,kk]));
  lines(lowess(fl[ok], y[ok,kk]), col="red", lwd=2);
}

# Plot normalized data
layout(matrix(1:9, ncol=3));
ylim <- c(-1,1)*max(y, na.rm=TRUE)/2;
for (kk in 1:I) {
  plot(fl, yN[,kk], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
  ok <- (is.finite(fl) & is.finite(y[,kk]));
  lines(lowess(fl[ok], yN[ok,kk]), col="blue", lwd=2);
}


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example 2: Two-enzyme fragment-length normalization of 6 arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xbeef);

# Number samples
I <- 5;

# Number of loci
J <- 3000;

# Fragment lengths (two enzymes)
fl <- matrix(0, nrow=J, ncol=2);
fl[,1] <- seq(from=100, to=1000, length.out=J);
fl[,2] <- seq(from=1000, to=100, length.out=J);

# Let 1/2 of the units be on both enzymes
fl[seq(from=1, to=J, by=4),1] <- NA;
fl[seq(from=2, to=J, by=4),2] <- NA;

# Let some have unknown fragment lengths
hasUnknownFL <- seq(from=1, to=J, by=15);
fl[hasUnknownFL,] <- NA;

# Sty/Nsp mixing proportions:
rho <- rep(1, I);
rho[1] <- 1/3;  # Less Sty in 1st sample
rho[3] <- 3/2;  # More Sty in 3rd sample


# Simulate data
z <- array(0, dim=c(J,2,I));
maxLog2Theta <- 12;
for (ii in 1:I) {
  # Common effect for both enzymes
  mu <- function(fl) {
    k <- runif(n=1, min=3, max=5);
    mu <- rep(maxLog2Theta, length(fl));
    ok <- is.finite(fl);
    mu[ok] <- mu[ok] - fl[ok]^{1/k};
    mu;
  }

  # Calculate the effect for each data point
  for (ee in 1:2) {
    z[,ee,ii] <- mu(fl[,ee]);
  }

  # Update the Sty/Nsp mixing proportions
  ee <- 2;
  z[,ee,ii] <- rho[ii]*z[,ee,ii];

  # Add random errors
  for (ee in 1:2) {
    eps <- rnorm(J, mean=0, sd=1/sqrt(2));
    z[,ee,ii] <- z[,ee,ii] + eps;
  }
}


hasFl <- is.finite(fl);

unitSets <- list(
  nsp  = which( hasFl[,1] & !hasFl[,2]),
  sty  = which(!hasFl[,1] &  hasFl[,2]),
  both = which( hasFl[,1] &  hasFl[,2]),
  none = which(!hasFl[,1] & !hasFl[,2])
)

# The observed data is a mix of two enzymes
theta <- matrix(NA, nrow=J, ncol=I);

# Single-enzyme units
for (ee in 1:2) {
  uu <- unitSets[[ee]];
  theta[uu,] <- 2^z[uu,ee,];
}

# Both-enzyme units (sum on intensity scale)
uu <- unitSets$both;
theta[uu,] <- (2^z[uu,1,]+2^z[uu,2,])/2;

# Missing units (sample from the others)
uu <- unitSets$none;
theta[uu,] <- apply(theta, MARGIN=2, sample, size=length(uu))

# Calculate target array
thetaT <- rowMeans(theta, na.rm=TRUE);
targetFcns <- list();
for (ee in 1:2) {
  uu <- unitSets[[ee]];
  fit <- lowess(fl[uu,ee], log2(thetaT[uu]));
  class(fit) <- "lowess";
  targetFcns[[ee]] <- function(fl, ...) {
    predict(fit, newdata=fl);
  }
}


# Fit model only to a subset of the data
subsetToFit <- setdiff(1:J, seq(from=1, to=J, by=10))

# Normalize data (to a target baseline)
thetaN <- matrix(NA, nrow=J, ncol=I);
fits <- vector("list", I);
for (ii in 1:I) {
  lthetaNi <- normalizeFragmentLength(log2(theta[,ii]), targetFcns=targetFcns,
                     fragmentLengths=fl, onMissing="median",
                     subsetToFit=subsetToFit, .returnFit=TRUE);
  fits[[ii]] <- attr(lthetaNi, "modelFit");
  thetaN[,ii] <- 2^lthetaNi;
}


# Plot raw data
xlim <- c(0, max(fl, na.rm=TRUE));
ylim <- c(0, max(log2(theta), na.rm=TRUE));
Mlim <- c(-1,1)*4;
xlab <- "Fragment length";
ylab <- expression(log2(theta));
Mlab <- expression(M==log[2](theta/theta[R]));

layout(matrix(1:(3*I), ncol=I, byrow=TRUE));
for (ii in 1:I) {
  plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main="raw");

  # Single-enzyme units
  for (ee in 1:2) {
    # The raw data
    uu <- unitSets[[ee]];
    points(fl[uu,ee], log2(theta[uu,ii]), col=ee+1);
  }

  # Both-enzyme units (use fragment-length for enzyme #1)
  uu <- unitSets$both;
  points(fl[uu,1], log2(theta[uu,ii]), col=3+1);

  for (ee in 1:2) {
    # The true effects
    uu <- unitSets[[ee]];
    lines(lowess(fl[uu,ee], log2(theta[uu,ii])), col="black", lwd=4, lty=3);

    # The estimated effects
    fit <- fits[[ii]][[ee]]$fit;
    lines(fit, col="orange", lwd=3);

    muT <- targetFcns[[ee]](fl[uu,ee]);
    lines(fl[uu,ee], muT, col="cyan", lwd=1);
  }
}

# Calculate log-ratios
thetaR <- rowMeans(thetaN, na.rm=TRUE);
M <- log2(thetaN/thetaR);

# Plot normalized data
for (ii in 1:I) {
  plot(NA, xlim=xlim, ylim=Mlim, xlab=xlab, ylab=Mlab, main="normalized");
  # Single-enzyme units
  for (ee in 1:2) {
    # The normalized data
    uu <- unitSets[[ee]];
    points(fl[uu,ee], M[uu,ii], col=ee+1);
  }
  # Both-enzyme units (use fragment-length for enzyme #1)
  uu <- unitSets$both;
  points(fl[uu,1], M[uu,ii], col=3+1);
}

ylim <- c(0,1.5);
for (ii in 1:I) {
  data <- list();
  for (ee in 1:2) {
    # The normalized data
    uu <- unitSets[[ee]];
    data[[ee]] <- M[uu,ii];
  }
  uu <- unitSets$both;
  if (length(uu) > 0)
    data[[3]] <- M[uu,ii];

  uu <- unitSets$none;
  if (length(uu) > 0)
    data[[4]] <- M[uu,ii];

  cols <- seq(along=data)+1;
  plotDensity(data, col=cols, xlim=Mlim, xlab=Mlab, main="normalized");

  abline(v=0, lty=2);
}


}

\author{Henrik Bengtsson}

\references{
  [1] H. Bengtsson, R. Irizarry, B. Carvalho, and T. Speed, \emph{Estimation and assessment of raw copy numbers at the single locus level}, Bioinformatics, 2008.
 \cr
}



\keyword{nonparametric}
\keyword{robust}
