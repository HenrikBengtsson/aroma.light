%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  calibrateMultiscan.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{calibrateMultiscan.matrix}
\alias{calibrateMultiscan.matrix}
\alias{calibrateMultiscan.matrix}


 \encoding{latin1}

 \title{Weighted affine calibration of a multiple re-scanned channel}

 \description{
   Weighted affine calibration of a multiple re-scanned channel.
 }

 \usage{\method{calibrateMultiscan}{matrix}(X, weights=NULL, typeOfWeights=c("datapoint"), method="L1", constraint="diagonal", satSignal=2^16 - 1, ..., average=median, deviance=NULL, project=FALSE, .fitOnly=FALSE)}

 \arguments{
  \item{X}{An NxK \code{\link[base]{matrix}} (K>=2) where the columns represent the
    multiple scans of one channel (a two-color array contains two
    channels) to be calibrated.}
  \item{weights}{If \code{\link[base]{NULL}}, non-weighted normalization is done.
    If data-point weights are used, this should be a \code{\link[base]{vector}} of length
    N of data point weights used when estimating the normalization
    function.
  }
  \item{typeOfWeights}{A \code{\link[base]{character}} string specifying the type of
    weights given in argument \code{weights}.
  }
  \item{method}{A \code{\link[base]{character}} string specifying how the estimates are
    robustified.  See \code{\link{iwpca}}() for all accepted values.}
  \item{constraint}{Constraint making the bias parameters identifiable.
    See \code{\link[aroma.light:fitIWPCA.matrix]{*fitIWPCA}()} for more details.}
  \item{satSignal}{Signals equal to or above this threshold is considered
    saturated signals.}
  \item{...}{Other arguments passed to \code{\link[aroma.light:fitIWPCA.matrix]{*fitIWPCA}()} and in
   turn \code{\link{iwpca}}().}
  \item{average}{A \code{\link[base]{function}} to calculate the average signals between calibrated scans.}
  \item{deviance}{A \code{\link[base]{function}} to calculate the deviance of the signals between calibrated scans.}
  \item{project}{If \code{\link[base:logical]{TRUE}}, the calibrated data points projected onto the
    diagonal line, otherwise not. Moreover, if \code{\link[base:logical]{TRUE}}, argument
    \code{average} is ignored.}
  \item{.fitOnly}{If \code{\link[base:logical]{TRUE}}, the data will not be back-transform.}
 }

 \value{
   If \code{average} is specified or \code{project} is \code{\link[base:logical]{TRUE}},
   an Nx1 \code{\link[base]{matrix}} is returned, otherwise an NxK \code{\link[base]{matrix}} is returned.
   If \code{deviance} is specified, a deviance Nx1 \code{\link[base]{matrix}} is returned
   as attribute \code{deviance}.
   In addition, the fitted model is returned as attribute \code{modelFit}.
 }

 \section{Negative, non-positive, and saturated values}{
   Affine multiscan calibration applies also to negative values, which are
   therefor also calibrated, if they exist.

   Saturated signals in any scan are set to \code{\link[base]{NA}}. Thus, they will not be
   used to estimate the calibration function, nor will they affect an
   optional projection.
 }

 \section{Missing values}{
   Only observations (rows) in \code{X} that contain all finite values are
   used in the estimation of the alibration functions. Thus,
   observations can be excluded by setting them to \code{\link[base]{NA}}.
 }

 \section{Weighted normalization}{
  Each data point/observation, that is, each row in \code{X}, which is a
  vector of length K, can be assigned a weight in [0,1] specifying how much
  it should \emph{affect the fitting of the calibration function}.
  Weights are given by argument \code{weights},
  which should be a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length N. Regardless of weights,
  all data points are \emph{calibrated} based on the fitted calibration
  function.
 }

 \section{Robustness}{
  By default, the model fit of multiscan calibration is done in \eqn{L_1}
  (\code{method="L1"}). This way, outliers affect the parameter estimates
  less than ordinary least-square methods.

  When calculating the average calibrated signal from multiple scans,
  by default the median is used, which further robustify against outliers.

  For further robustness, downweight outliers such as saturated signals,
  if possible.

  Tukey's biweight function is supported, but not used by default because
  then a "bandwidth" parameter has to selected. This can indeed be done
  automatically by estimating the standard deviation, for instance using
  MAD. However, since scanner signals have heteroscedastic noise
  (standard deviation is approximately proportional to the non-logged
  signal), Tukey's bandwidth parameter has to be a function of the
  signal too, cf. \code{\link[stats]{loess}}.  We have experimented with this
  too, but found that it does not significantly improve the robustness
  compared to \eqn{L_1}.
  Moreover, using Tukey's biweight as is, that is, assuming homoscedastic
  noise, seems to introduce a (scale dependent) bias in the estimates
  of the offset terms.
 }

 \details{
  Fitting is done by iterated re-weighted principal component analysis
  (IWPCA).
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   [1] H. Bengtsson, J. Vallon-Christersson and G. \enc{Jönsson}{Jonsson}, \emph{Calibration and assessment of channel-specific biases in microarray data with extended dynamical range}, BMC Bioinformatics, 5:177, 2004.
 \cr
 }

 \examples{\dontrun{# For an example, see help(normalizeAffine).}}

 \seealso{
   \code{\link{1. Calibration and Normalization}}.
   \code{\link[aroma.light:normalizeAffine.matrix]{*normalizeAffine}()}.
   For more information see \code{\link{matrix}}.
 }
\keyword{methods}
