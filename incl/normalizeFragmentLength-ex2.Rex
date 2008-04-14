# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example 2: Two-enzyme fragment-length normalization of 6 arrays
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set.seed(0xbeef);

# Number samples
I <- 5;

# Number of loci
J <- 2000;

# Fragment lengths (two enzymes)
fl <- matrix(0, nrow=J, ncol=2);
fl[,1] <- seq(from=100, to=1000, length.out=J);
fl[,2] <- seq(from=1000, to=100, length.out=J);

# Let 1/2 of the units be on both enzymes
fl[seq(from=1, to=J, by=4),1] <- NA;
fl[seq(from=2, to=J, by=4),2] <- NA;

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
    maxLog2Theta - fl^{1/k};
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
  both = which( hasFl[,1] &  hasFl[,2])
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


# Normalize data (to a target baseline)
thetaN <- matrix(NA, nrow=J, ncol=I);
fits <- vector("list", I);
for (ii in 1:I) {
  lthetaNi <- normalizeFragmentLength(log2(theta[,ii]), targetFcns=targetFcns, fragmentLengths=fl, .returnFit=TRUE);
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
    lines(fit, col="orange", lwd=2);

    muT <- targetFcns[[ee]](fl[uu,ee]);
    lines(fl[uu,ee], muT, col="purple", lwd=1);
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
  plotDensity(data, col=1:3+1, xlim=Mlim, xlab=Mlab, main="normalized");
  abline(v=0, lty=2);
}

