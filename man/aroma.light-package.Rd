%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{aroma.light-package}
\alias{aroma.light-package}
\alias{aroma.light}
\docType{package}

\title{Package aroma.light}


\encoding{latin1}

\description{
  Methods for microarray analysis that take basic data types such as matrices and lists of vectors.  These methods can be used standalone, be utilized in other packages, or be wrapped up in higher-level classes.
}

\section{Installation}{
  To install this package, see
  \url{http://bioconductor.org/packages/release/bioc/html/aroma.light.html}.
}

\section{To get started}{
  For scanner calibration:
  \enumerate{
    \item see \code{\link{calibrateMultiscan}}() - scan the same array two or more times to calibrate for scanner effects and extended dynamical range.
  }

  To normalize multiple single-channel arrays all with the same number of probes/spots:
  \enumerate{
    \item \code{\link{normalizeAffine}}() - normalizes, on the intensity scale,  for differences in offset and scale between channels.
    \item \code{\link{normalizeQuantileRank}}(), \code{\link{normalizeQuantileSpline}}() - normalizes, on the intensity scale,  for differences in empirical distribution between channels.
  }

  To normalize multiple single-channel arrays with varying number probes/spots:
  \enumerate{
    \item \code{\link{normalizeQuantileRank}}(), \code{\link{normalizeQuantileSpline}}() - normalizes, on the intensity scale, for differences in empirical distribution between channels.
  }

  To normalize two-channel arrays:
  \enumerate{
    \item \code{\link{normalizeAffine}}() - normalizes, on the intensity scale, for differences in offset and scale between channels.  This will also correct for intensity-dependent affects on the log scale.
    \item \code{\link{normalizeCurveFit}}() - Classical intensity-dependent normalization, on the log scale, e.g. lowess normalization.
  }

  To normalize three or more channels:
  \enumerate{
    \item \code{\link{normalizeAffine}}() - normalizes, on the intensity scale, for differences in offset and scale between channels.  This will minimize the curvature on the log scale between any two channels.
  }
}

\section{Further readings}{
  Several of the normalization methods proposed in [1]-[7] are
  available in this package.
}

\section{How to cite this package}{
  Whenever using this package, please cite one or more of [1]-[7].
}

\section{Wishlist}{
 Here is a list of features that would be useful, but which I have
 too little time to add myself. Contributions are appreciated.
 \itemize{
   \item At the moment, nothing.
 }

 If you consider to contribute, make sure it is not already
 implemented by downloading the latest "devel" version!
}

\author{Henrik Bengtsson, Pierre Neuvial}

\section{License}{
  The releases of this package is licensed under
  GPL version 2 or newer.

  NB: Except for the \code{robustSmoothSpline()} method,
  it is alright to distribute the rest of the package under
  LGPL version 2.1 or newer.

  The development code of the packages is under a private licence
  (where applicable) and patches sent to the author fall under the
  latter license, but will be, if incorporated, released under the
  "release" license above.
}

\references{
 Some of the reference below can be found at
 \url{http://www.aroma-project.org/publications/}.\cr

[1] H. Bengtsson, \emph{Identification and normalization of plate effects
    in cDNA microarray data}, Preprints in Mathematical Sciences,
    2002:28, Mathematical Statistics, Centre for Mathematical Sciences,
    Lund University, 2002.\cr

[2] H. Bengtsson, \emph{The R.oo package - Object-Oriented Programming with References Using Standard R Code}, In Kurt Hornik, Friedrich Leisch and Achim Zeileis, editors, Proceedings of the 3rd International Workshop on Distributed Statistical Computing (DSC 2003), March 20-22, Vienna, Austria. \url{http://www.ci.tuwien.ac.at/Conferences/DSC-2003/Proceedings/}
 \cr

[3] H. Bengtsson, \emph{aroma - An R Object-oriented Microarray
    Analysis environment}, Preprints in Mathematical Sciences (manuscript
    in preparation), Mathematical Statistics, Centre for Mathematical
    Sciences, Lund University, 2004.\cr

[4] H. Bengtsson, J. Vallon-Christersson and G. \enc{Jönsson}{Jonsson}, \emph{Calibration and assessment of channel-specific biases in microarray data with extended dynamical range}, BMC Bioinformatics, 5:177, 2004.
 \cr

[5] Henrik Bengtsson and Ola \enc{Hössjer}{Hossjer}, \emph{Methodological Study of Affine Transformations of Gene Expression Data}, Methodological study of affine transformations of gene expression data with proposed robust non-parametric multi-dimensional normalization method, BMC Bioinformatics, 2006, 7:100.
 \cr

[6] H. Bengtsson, R. Irizarry, B. Carvalho, and T. Speed, \emph{Estimation and assessment of raw copy numbers at the single locus level}, Bioinformatics, 2008.
 \cr

[7] H. Bengtsson, A. Ray, P. Spellman and T.P. Speed, \emph{A single-sample method for normalizing and combining full-resolutioncopy numbers from multiple platforms, labs and analysis methods}, Bioinformatics, 2009.
 \cr

[8] H. Bengtsson, P. Neuvial and T.P. Speed, \emph{TumorBoost: Normalization of allele-specific tumor copy numbers from a single pair of tumor-normal genotyping microarrays}, BMC Bioinformatics, 2010, 11:245. [PMID 20462408]
 \cr
}
\keyword{package}
