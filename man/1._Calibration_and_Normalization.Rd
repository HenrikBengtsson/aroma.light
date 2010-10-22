%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  901.CalibrationAndNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{1. Calibration and Normalization}
\alias{1. Calibration and Normalization}
\title{1. Calibration and Normalization}


 \encoding{latin1}

 \description{
   In this section we give \emph{our} recommendation on how spotted
   two-color (or multi-color) microarray data is best calibrated and
   normalized.
 }

 \section{Classical background subtraction}{
   We do \emph{not} recommend background subtraction in classical
   means where background is estimated by various image analysis
   methods.  This means that we will only consider foreground signals
   in the analysis.

   We estimate "background" by other means. In what is explain below,
   only a global background, that is, a global bias, is estimated
   and removed.
 }

 \section{Multiscan calibration}{
   In Bengtsson et al (2004) we give evidence that microarray scanners
   can introduce a significant bias in data. This bias, which is
   about 15-25 out of 65535, \emph{will} introduce intensity dependency
   in the log-ratios, as explained in Bengtsson &
   \enc{Hössjer}{Hossjer} (2006).

   In Bengtsson et al (2004) we find that this bias is stable across
   arrays (and a couple of months), but further research is needed
   in order to tell if this is true over a longer time period.

   To calibrate signals for scanner biases, scan the same array at
   multiple PMT-settings at three or more (K >= 3) different
   PMT settings (preferably in decreasing order).
   While doing this, \emph{do not adjust the laser power settings}.
   Also, do the multiscan \emph{without} washing, cleaning or by other
   means changing the array between subsequent scans.
   Although not necessary, it is preferred that the array
   remains in the scanner between subsequent scans. This will simplify
   the image analysis since spot identification can be made once
   if images aligns perfectly.

   After image analysis, read all K scans for the same array into the
   two matrices, one for the red and one for the green channel, where
   the K columns corresponds to scans and the N rows to the spots.
   It is enough to use foreground signals.

   In order to multiscan calibrate the data, for each channel
   separately call \code{Xc <- calibrateMultiscan(X)} where \code{X}
   is the NxK matrix of signals for one channel across all scans.  The
   calibrated signals are returned in the Nx1 matrix \code{Xc}.

   Multiscan calibration may sometimes be skipped, especially if affine
   normalization is applied immediately after, but we do recommend that
   every lab check at least once if their scanner introduce bias.
 }

 \section{Affine normalization}{
   In Bengtsson & \enc{Hössjer}{Hossjer} (2006), we carry out a detailed
   study on how biases in each channel introduce so called
   intensity-dependent log-ratios among other systematic artifacts.
   Data with (additive) bias in each channel is said to be \emph{affinely}
   transformed. Data without such bias, is said to be \emph{linearly}
   (proportionally) transform. Ideally, observed signals (data) is a
   linear (proportional) function of true gene expression levels.
 
   We do \emph{not} assume proportional observations. The scanner bias
   is real evidence that assuming linearity is not correct.
   Affine normalization corrects for affine transformation in data.
   Without control spots it is not possible to estimate the bias in each
   of the channels but only the relative bias such that after
   normalization the effective bias are the same in all channels.
   This is why we call it normalization and not calibration.

   In its simplest form, affine normalization is done by
   \code{Xn <- normalizeAffine(X)} where \code{X} is a Nx2 matrix with
   the first column holds the foreground signals from the red channel and
   the second holds the signals from the green channel.  If three- or
   four-channel data is used these are added the same way.  The normalized
   data is returned as a Nx2 matrix \code{Xn}.

   To normalize all arrays and all channels at once, one may put all
   data into one big NxK matrix where the K columns hold the all channels
   from the first array, then all channels from the second array and so
   on.  Then \code{Xn <- normalizeAffine(X)} will return the across-array
   and across-channel normalized data in the NxK matrix \code{Xn} where
   the colunms are stored in the same order as in matrix \code{X}.

   Equal effective bias in all channels is much better. First of all,
   any intensity-dependent bias in the log-ratios is removed \emph{for
   all non-differentially expressed genes}. There is still an
   intensity-dependent bias in the log-ratios for differentially expressed
   genes, but this is now symmetric around log-ratio zero.

   Affine normalization will (by default and recommended) normalize
   \emph{all} arrays together and at once. This will guarantee that
   all arrays are "on the same scale". Thus, it \emph{not} recommended
   to apply a classical between-array scale normalization afterward.
   Moreover, the average log-ratio will be zero after an affine
   normalization.
 
   Note that an affine normalization will only remove curvature in the
   log-ratios at lower intensities.
   If a strong intensity-dependent bias at high intensities remains,
   this is most likely due to saturation effects, such as too high PMT
   settings or quenching.

   Note that for a perfect affine normalization you \emph{should}
   expect much higher noise levels in the \emph{log-ratios} at lower
   intensities than at higher. It should also be approximately
   symmetric around zero log-ratio.
   In other words, \emph{a strong fanning effect is a good sign}.

   Due to different noise levels in red and green channels, different
   PMT settings in different channels, plus the fact that the
   minimum signal is zero, "odd shapes" may be seen in the log-ratio
   vs log-intensity graphs at lower intensities. Typically, these
   show themselves as non-symmetric in positive and negative log-ratios.
   Note that you should not see this at higher intensities.

   If there is a strong intensity-dependent effect left after the
   affine normalization, we recommend, for now, that a subsequent
   curve-fit or quantile normalization is done.
   Which one, we do not know.

   Why negative signals?
   By default, 5\% of the normalized signals will have a non-positive
   signal in one or both channels. \emph{This is on purpose}, although
   the exact number 5\% is chosen by experience. The reason for
   introducing negative signals is that they are indeed expected.
   For instance, when measure a zero gene expression level, there is
   a chance that the observed value is (should be) negative due to
   measurement noise. (For this reason it is possible that the scanner
   manufacturers have introduced scanner bias on purpose to avoid
   negative signals, which then all would be truncated to zero.)
   To adjust the ratio (or number) of negative signals allowed, use
   for example \code{normalizeAffine(X, constraint=0.01)} for 1\%
   negative signals. If set to zero (or \code{"max"}) only as much
   bias is removed such that no negative signals exist afterward.
   Note that this is also true if there were negative signals on
   beforehand.

   Why not lowess normalization?
   Curve-fit normalization methods such as lowess normalization are
   basically designed based on linearity assumptions and will for this
   reason not correct for channel biases. Curve-fit normalization
   methods can by definition only be applied to one pair of channels
   at the time and do therefore require a subsequent between-array
   scale normalization, which is by the way very ad hoc.

   Why not quantile normalization?
   Affine normalization can be though of a special case of quantile
   normalization that is more robust than the latter.
   See Bengtsson & \enc{Hössjer}{Hossjer} (2006) for details.
   Quantile normalization is probably better to apply than curve-fit
   normalization methods, but less robust than affine normalization,
   especially at extreme (low and high) intensities.
   For this reason, we do recommend to use affine normalization first,
   and if this is not satisfactory, quantile normalization may be applied.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{documentation}
