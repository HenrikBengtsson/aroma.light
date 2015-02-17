library("aroma.light")

pathname <- system.file("data-ex", "PMT-RGData.dat", package="aroma.light")
rg <- read.table(pathname, header=TRUE, sep="\t")
nbrOfScans <- max(rg$slide)

rg <- as.list(rg)
for (field in c("R", "G"))
  rg[[field]] <- matrix(as.double(rg[[field]]), ncol=nbrOfScans)
rg$slide <- rg$spot <- NULL
rg <- as.matrix(as.data.frame(rg))
colnames(rg) <- rep(c("R", "G"), each=nbrOfScans)

rgC <- rg

layout(matrix(c(1,2,0,3,4,0,5,6,7), ncol=3, byrow=TRUE))

for (channel in c("R", "G")) {
  sidx <- which(colnames(rg) == channel)
  channelColor <- switch(channel, R="red", G="green")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The raw data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  plotMvsAPairs(rg, channel=channel)
  title(main=paste("Observed", channel))
  box(col=channelColor)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The calibrated data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rgC[,sidx] <- calibrateMultiscan(rg[,sidx], average=NULL)

  plotMvsAPairs(rgC, channel=channel)
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
rgCA <- matrix(NA, nrow=nrow(rg), ncol=2)
colnames(rgCA) <- c("R", "G");
for (channel in c("R", "G")) {
  sidx <- which(colnames(rg) == channel)
  rgCA[,channel] <- calibrateMultiscan(rg[,sidx])
}

plotMvsA(rgCA)
title(main="Average calibrated")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The affine normalized average calibrated data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create a matrix where the columns represent the channels
# to be normalized.
rgCAN <- rgCA
# Affine normalization of channels
rgCAN <- normalizeAffine(rgCAN)

plotMvsA(rgCAN)
title(main="Affine normalized A.C.")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# It is always ok to rescale the affine normalized data if its
# done on (R,G); not on (A,M)! However, this is only needed for
# esthetic purposes.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rgCAN <- rgCAN * 2^5
plotMvsA(rgCAN)
title(main="Rescaled normalized")


