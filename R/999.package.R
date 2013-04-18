#########################################################################/**
# @RdocPackage aroma.light
#
# \encoding{latin1}
#
# \description{
#   @eval "getDescription(aroma.light)"
# }
#
# \section{Requirements}{
#   This package requires the \pkg{R.oo} package [1].
# }
#
# \section{Installation}{
#   To install this package, see \url{http://www.braju.com/R/}.
#   Required packages are installed in the same way.
# }
#
# \section{To get started}{
#   For scanner calibration:
#   \enumerate{
#     \item see @see "calibrateMultiscan.matrix" - scan the same array two or more times to calibrate for scanner effects and extended dynamical range.
#   }
#
#   To normalize multiple single-channel arrays all with the same number of probes/spots:
#   \enumerate{
#     \item @see "normalizeAffine.matrix" - normalizes, on the intensity scale,  for differences in offset and scale between channels.
#     \item @see "normalizeQuantileRank.matrix", @see "normalizeQuantileSpline.matrix" - normalizes, on the intensity scale,  for differences in empirical distribution between channels.
#   }
#
#   To normalize multiple single-channel arrays with varying number probes/spots:
#   \enumerate{
#     \item @see "normalizeQuantileRank.list", @see "normalizeQuantileSpline.list" - normalizes, on the intensity scale, for differences in empirical distribution between channels.
#   }
#
#   To normalize two-channel arrays:
#   \enumerate{
#     \item @see "normalizeAffine.matrix" - normalizes, on the intensity scale, for differences in offset and scale between channels.  This will also correct for intensity-dependent affects on the log scale.
#     \item @see "normalizeCurveFit.matrix" - Classical intensity-dependent normalization, on the log scale, e.g. lowess normalization.
#   }
#
#   To normalize three or more channels:
#   \enumerate{
#     \item @see "normalizeAffine.matrix" - normalizes, on the intensity scale, for differences in offset and scale between channels.  This will minimize the curvature on the log scale between any two channels.
#   }
# }
#
# \section{Further readings}{
#   Several of the normalization methods proposed in [3]-[6] are
#   available in this package.
# }
#
# \section{How to cite this package}{
#   Whenever using this package, please cite [2] as\cr
#
#   @howtocite "aroma.light"
# }
#
# \section{Wishlist}{
#  Here is a list of features that would be useful, but which I have
#  too little time to add myself. Contributions are appreciated.
#  \itemize{
#    \item At the moment, nothing.
#  }
#
#  If you consider to contribute, make sure it is not already
#  implemented by downloading the latest "devel" version!
# }
#
# @author "*"
#
# \section{License}{
#   The releases of this package is licensed under
#   GPL version 2 or newer.
#
#   NB: Except for the \code{robustSmoothSpline()} method,
#   it is alright to distribute the rest of the package under
#   LGPL version 2.1 or newer.
#
#   The development code of the packages is under a private licence
#   (where applicable) and patches sent to the author fall under the
#   latter license, but will be, if incorporated, released under the
#   "release" license above.
# }
#
# \references{
#  Some of the reference below can be found at
#  \url{http://www.maths.lth.se/bioinformatics/publications/}.\cr
#
# [1] @include "../incl/BengtssonH_2003.bib.Rdoc" \cr
#
# [2] H. Bengtsson, \emph{aroma - An R Object-oriented Microarray
#     Analysis environment}, Preprints in Mathematical Sciences (manuscript
#     in preparation), Mathematical Statistics, Centre for Mathematical
#     Sciences, Lund University, 2004.\cr
#
# [3] @include "../incl/BengtssonHossjer_2006.bib.Rdoc" \cr
#
# [4] @include "../incl/BengtssonH_etal_2004.bib.Rdoc" \cr
#
# [5] @include "../incl/BengtssonH_etal_2008.bib.Rdoc" \cr
#
# [6] H. Bengtsson, \emph{Identification and normalization of plate effects
#     in cDNA microarray data}, Preprints in Mathematical Sciences,
#     2002:28, Mathematical Statistics, Centre for Mathematical Sciences,
#     Lund University, 2002.\cr
# }
#*/#########################################################################
