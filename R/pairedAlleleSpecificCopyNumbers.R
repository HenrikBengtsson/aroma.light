###########################################################################/**
# @RdocGeneric pairedAlleleSpecificCopyNumbers
# @alias pairedAlleleSpecificCopyNumbers.numeric
#
# @title "Calculating tumor-normal paired allele-specific copy number stratified on genotypes"
#
# \description{
#  @get "title".
#  The method is a single-sample (single-pair) method.
#  It requires paired tumor-normal parent-specific copy number signals.
# }
#
# \usage{
#  @usage pairedAlleleSpecificCopyNumbers,numeric
# }
#
# \arguments{
#  \item{thetaT, betaT}{Theta and allele-B fraction signals for the tumor.}
#  \item{thetaN, betaN}{Total and allele-B fraction signals for the
#     matched normal.}
#  \item{muN}{An optional @vector of length J containing
#     normal genotypes calls in (0,1/2,1,@NA) for (AA,AB,BB).}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @data.frame with elements \code{CT}, \code{betaT} and \code{muN}.
# }
#
# \seealso{
#  This definition of calculating tumor-normal paired ASCN is related
#  to how the @see "normalizeTumorBoost" method calculates normalized
#  tumor BAFs.
# }
#
# @author "PN, HB"
#*/###########################################################################
setMethodS3("pairedAlleleSpecificCopyNumbers", "numeric", function(thetaT, betaT, thetaN, betaN, muN=callNaiveGenotypes(betaN), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument: 'thetaT' & 'betaT':
  thetaT <- as.numeric(thetaT);
  betaT <- as.numeric(betaT);
  J <- length(thetaT);
  if (length(betaT) != J) {
    stop("The length of arguments 'betaT' and 'thetaT' differ: ",
                                                   length(betaT), " != ", J);
  }

  # Argument: 'thetaN' & 'betaN':
  thetaN <- as.numeric(thetaN);
  betaN <- as.numeric(betaN);
  if (length(thetaN) != J) {
    stop("The length of arguments 'thetaN' and 'thetaT' differ: ",
                                                   length(thetaN), " != ", J);
  }
  if (length(betaN) != J) {
    stop("The length of arguments 'betaN' and 'thetaN' differ: ",
                                                   length(betaN), " != ", J);
  }

  # Argument: 'muN':
  if (length(muN) != J) {
    stop("The length of arguments 'muN' and 'betaN' differ: ",
                                                   length(muN), " != ", J);
  }
  knownGenotypes <- c(0, 1/2, 1, NA);
  unknown <- which(!is.element(muN, knownGenotypes));
  n <- length(unknown);
  if (n > 0L) {
    unknown <- unique(muN[unknown]);
    stop("Argument 'muN' contains unknown values: ", hpaste(unknown));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate (thetaA,thetaB) for tumor and normal (for SNP only)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # SNPs are identifies as those loci that have non-missing 'betaTN' & 'muN'
  isSnp <- (!is.na(betaT) & !is.na(muN));
  nbrOfSnps <- sum(isSnp);

  thetaTs <- thetaT[isSnp] * matrix(c(betaT[isSnp], 1-betaT[isSnp]), ncol=2L);
  stopifnot(nrow(thetaTs) == nbrOfSnps);

  thetaNs <- thetaN[isSnp] * matrix(c(betaN[isSnp], 1-betaN[isSnp]), ncol=2L);
  stopifnot(nrow(thetaNs) == nbrOfSnps);

  homAs <- which(muN[isSnp] == 0);
  homBs <- which(muN[isSnp] == 1);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate tumor (CA,CB) conditioned on genotype (for SNP only)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  CTs <- thetaTs / thetaNs;
  CTs[homAs,1L] <- 2*CTs[homAs,1L];
  CTs[homAs,2L] <- 0;
  CTs[homBs,1L] <- 0;
  CTs[homBs,2L] <- 2*CTs[homBs,2L];

  thetaTs <- thetaNs <- NULL; # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Translate (CA,CB) to (CT,betaT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Total CN ratios
  CT <- betaT <- rep(NA_real_, times=J);
  CT[isSnp] <- CTs[,1L] + CTs[,2L];
  CT[!isSnp] <- 2 * thetaT[!isSnp] / thetaN[!isSnp];

  # Tumor BAFs
  betaT[isSnp] <- CTs[,2L] / CT[isSnp];

  CTs <- NULL; # Not needed anymore

  stopifnot(length(CT) == J);
  stopifnot(length(betaT) == J);

  # Sanity check: All homozygous tumor BAFs should be {0, 1}
  stopifnot(all(na.omit(betaT[isSnp][homAs]) == 0));
  stopifnot(all(na.omit(betaT[isSnp][homBs]) == 1));
  isSnp <- homAs <- homBs <- NULL; # Not needed anymore

  # Return value
  data <- data.frame(CT=CT, betaT=betaT, muN=muN);

  data;
}) # pairedAlleleSpecificCopyNumbers()


##############################################################################
# HISTORY:
# 2014-03-30 [HB in Juvisy]
# o Created from PN's description just to be on the same page. PN has argued
#   for this way of calculating ASCN's PN since our TumorBoost days (~2009).
#   PN has a high-level implementation elsewhere, but HB decided to do this
#   from scratch to get a low-level API similar to that of TumorBoost.
##############################################################################
