###########################################################################/**
# @RdocDocumentation "Non-documented objects"
#
# % Calibration and normalization functions
# @alias fitIWPCA
# @alias calibrateMultiscan
# @alias normalizeAffine
# @alias normalizeAverage
# @alias backtransformAffine
# @alias normalizeQuantile
# @alias normalizeQuantile.default
# @alias normalizeQuantileRank
# @alias normalizeQuantileSpline
# @alias normalizeCurveFit
# @alias normalizeLoess
# @alias normalizeLowess
# @alias normalizeRobustSpline
# @alias normalizeSpline
# @alias averageQuantile
#
# % Plot functions
# @alias plotMvsA
# @alias plotMvsAPairs
# @alias plotMvsMPairs
# @alias plotDensity
# @alias plotXYCurve
# @alias lines.XYCurveFit
#
# % Matrix operations
# @alias rowAverages
# @alias rowAverages.matrix
# @alias sampleCorrelations
# @alias medianPolish
# @alias backtransformXYCurve
# @alias backtransformXYCurve.matrix
#
# % Simple linear-algebra
# @alias projectUontoV
# @alias scalarProduct
# @alias tr
#
# % Advanced linear-algebra
# @alias iwpca
# @alias wpca
# @alias fitPrincipalCurve
# @alias backtransformPrincipalCurve
#
# % Miscellaneous statistical functions
# @alias likelihood
# @alias predict.lowess
#
# \description{
#   This page contains aliases for all "non-documented" objects that
#   \code{R CMD check} detects in this package.
#
#   Almost all of them are \emph{generic} functions that have specific
#   document for the corresponding method coupled to a specific class.
#   Other functions are re-defined by \code{setMethodS3()} to
#   \emph{default} methods. Neither of these two classes are non-documented
#   in reality.
#   The rest are deprecated methods.
# }
#
# @keyword internal
#*/###########################################################################

############################################################################
# HISTORY:
# 2005-02-10
# o Created to please R CMD check.
############################################################################
