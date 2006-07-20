library(aroma)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plotMvsA <- function(rg, xlab="A", ylab="M", Alim=c(0,16), main="", ..., add=FALSE) {
  if (!add) {
    plot(NA, xlab=xlab, ylab=ylab, xlim=Alim, ylim=c(-1,1)*diff(Alim))
    title(main=main)
  }
  
  R <- rg[,"R"]
  G <- rg[,"G"]
  M <- log(R/G, base=2)
  A <- log(R*G, base=2)/2
  points(A,M, pch=".", ...)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Load data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# One array was scanned four times at four different PMT settings.
rg0 <- RGData$read("PMT-RGData.dat", path=aroma$dataPath);
setLayout(rg0, Layout(4,4,17,17));       # Not really necessary!
setSlideName(rg0, c("500V","600V","700V","800V"));
keepSlides(rg0, slide="600V")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# "Add" some outliers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
outliers <- 1:50
rg0$G[outliers,] <- runif(n=length(outliers), min=40000, max=55000)

Alim <- c(-2,15)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Raw data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
png("../figures/NormalizeAffine-raw.png", width=210, height=210)
opar <- par(cex=0.6)
plot(as.MAData(rg0), xlim=Alim)
dev.off()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Non-weighted affine normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rgNW <- clone(rg0)
normalizeAffine(rgNW)

png("../figures/NormalizeAffine-nonweighted.png", width=210, height=210)
opar <- par(cex=0.6)
plot(as.MAData(rgNW), xlim=Alim)
dev.off()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Weighted affine normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rgW <- clone(rg0);

# Down-weight outliers
weights <- rep(1, nbrOfSpots(rg0))
weights[outliers] <- 1/200
setProbeWeights(rgW, weights=weights)

normalizeAffine(rgW)

png("../figures/NormalizeAffine-weighted.png", width=210, height=210)
opar <- par(cex=0.6)
plot(as.MAData(rgW), xlim=Alim)
dev.off()

