###########################################################################/**
# @RdocDefault normalizeFragmentLength
#
# @title "Normalizes signals for PCR fragment-length effects"
#
# \description{
#  @get "title". Some or all signals are used to estimated the 
#  normalization function.  All signals are normalized.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{y}{A @numeric @vector of length K of signals to be normalized
#     across E enzymes.}
#   \item{fragmentLengths}{An @integer KxE @matrix of fragment lengths.}
#   \item{targetFcns}{An optional @list of E @functions; one per enzyme.
#     If @NULL, the data is normalized to have constant fragment-length
#     effects (all equal to zero on the log-scale).}
#   \item{subsetToFit}{The subset of data points used to fit the 
#      normalization function.
#      If @NULL, all data points are considered.}
#   \item{onMissing}{Specifies how data points for which there is no
#      fragment length is normalized. 
#      If \code{"ignore"}, the values are not modified.
#      If \code{"median"}, the values are updated to have the same
#      robust average as the other data points.
#   }
#   \item{.isLogged}{A @logical.}
#   \item{...}{Additional arguments passed to @see "stats::lowess".}
#   \item{.returnFit}{A @logical.}
# }
#
# \value{
#   Returns a @numeric @vector of the normalized signals.
# }
#
# \section{Multi-enzyme normalization}{
#  It is assumed that the fragment-length effects from multiple enzymes
#  added (with equal weights) on the intensity scale.
#  The fragment-length effects are fitted for each enzyme separately based
#  on units that are exclusively for that enzyme. 
#  \emph{If there are no or very such units for an enzyme, the assumptions 
#  of the model are not met and the fit will fail with an error.}
#  Then, from the above single-enzyme fits the average effect across
#  enzymes is the calculated for each unit that is on multiple enzymes.
# }
#
# \section{Target functions}{
#   It is possible to specify custom target function effects for each
#   enzyme via argument \code{targetFcns}.  This argument has to be a
#   @list containing one @function per enzyme and ordered in the same
#   order as the enzyme are in the columns of argument 
#   \code{fragmentLengths}.
#   For instance, if one wish to normalize the signals such that their
#   mean signal as a function of fragment length effect is contantly
#   equal to 2200 (or the intensity scale), the use 
#   \code{targetFcns=function(fl, ...) log2(2200)} which completely
#   ignores fragment-length argument 'fl' and always returns a 
#   constant.
#   If two enzymes are used, then use
#   \code{targetFcns=rep(list(function(fl, ...) log2(2200)), 2)}.
#
#   Note, if \code{targetFcns} is @NULL, this corresponds to 
#   \code{targetFcns=rep(list(function(fl, ...) 0), ncol(fragmentLengths))}.
#
#   Alternatively, if one wants to only apply minimial corrections to
#   the signals, then one can normalize toward target functions that
#   correspond to the fragment-length effect of the average array.
# }
#
# \examples{
#   @include "../incl/normalizeFragmentLength-ex1.Rex"
#
#   @include "../incl/normalizeFragmentLength-ex2.Rex"
# }
#
# @author
#
# \references{
#   [1] @include "../incl/BengtssonH_etal_2008.bib.Rdoc" \cr
# } 
#
# @keyword "nonparametric"
# @keyword "robust" 
#*/###########################################################################
setMethodS3("normalizeFragmentLength", "default", function(y, fragmentLengths, targetFcns=NULL, subsetToFit=NULL, onMissing=c("median", "ignore"), .isLogged=TRUE, ..., .returnFit=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  y <- as.double(y);
  nbrOfDataPoints <- length(y);

  # Argument 'fragmentLengths':
  if (!is.matrix(fragmentLengths)) {
    if (is.vector(fragmentLengths)) {
      fragmentLengths <- as.matrix(fragmentLengths);
    } else {
      throw("Argument 'fragmentLengths' must be a matrix: ", 
                                                class(fragmentLengths)[[1]]);
    }
  }
  nbrOfEnzymes <- ncol(fragmentLengths);
  allEnzymes <- seq(length=nbrOfEnzymes);
  for (ee in allEnzymes) {
    fragmentLengths[,ee] <- as.double(fragmentLengths[,ee]);
  }

  # Argument 'targetFcns':
  if (!is.null(targetFcns)) {
    if (!is.list(targetFcns)) {
      if (nbrOfEnzymes == 1) {
        targetFcns <- list(targetFcns);
      } else {
        throw("Argument 'targetFcns' is not a list: ", class(targetFcns)[1]);
      }
    }
    if (length(targetFcns) != nbrOfEnzymes) {
      throw("Number of elements in 'targetFcns' does not match the number of columns in 'fragmentLengths': ", length(targetFcns), " != ", nbrOfEnzymes);
    }

    # Validate each element
    for (ee in allEnzymes) {
      if (!is.function(targetFcns[[ee]])) {
        throw("One element in 'targetFcns' is not a function: ", class(targetFcns[[ee]])[1]);
      }
    }
  }
  
  # Argument 'subsetToFit':
  if (!is.null(subsetToFit)) {
    subsetToFit <- as.integer(subsetToFit);
    if (length(subsetToFit) > nbrOfDataPoints) {
      throw("The length of argument 'subsetToFit' does not match the number of data points: ", length(subsetToFit), " != ", nbrOfDataPoints);
    }
  }

  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate normalization function and predict the signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit smooth curve to each enzyme separately
  hasFL <- is.finite(fragmentLengths);
  # Count the number of enzymes per units
  countFL <- rep(0, nbrOfDataPoints);
  for (ee in allEnzymes)
    countFL <- countFL + as.integer(hasFL[,ee]);
  isSingleEnzymed <- (countFL == 1);

  okY <- is.finite(y);

  # KxE matrix for sample (and target predictions)
  naValue <- as.double(NA);
  mu <- matrix(naValue, nrow=nbrOfDataPoints, ncol=nbrOfEnzymes);
  if (!is.null(targetFcns))
    muT <- matrix(naValue, nrow=nbrOfDataPoints, ncol=nbrOfEnzymes);

  if (.returnFit)
    fits <- vector("list", nbrOfEnzymes);

  for (ee in allEnzymes) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (a) Fit normalization function
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit only to units with known length and non-missing data points.
    ok <- (hasFL[,ee] & isSingleEnzymed & okY);

    # Sanity check
    if (sum(ok) == 0) {
      throw("Cannot fit normalization function to enzyme, because there are no (finite) data points that are unique to this enzyme: ", ee);
    }

    if (!is.null(subsetToFit)) {
      ok[-subsetToFit] <- FALSE;

      # Sanity check
      if (sum(ok) == 0) {
        throw("Cannot fit normalization function to enzyme, because there are no (finite) data points that are unique to this enzyme for the subset requested: ", ee);
      }
    }


    # All fragment lengths for current enzyme
    fl <- fragmentLengths[,ee];

    # Fit finite {(lambda, log2theta)_j} to data points j on current enzyme
    suppressWarnings({
      fit <- lowess(fl[ok], y[ok], ...);
      class(fit) <- "lowess";
    })
    rm(ok);

    if (.returnFit)
      fits[[ee]] <- fit;

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (b) Calculate correction factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate the correction factor for every data point on this enzyme
    ok <- (hasFL[,ee] & okY);
    mu[ok,ee] <- predict(fit, newdata=fl[ok]);

    if (.returnFit)
      fits[[ee]] <- list(fit=fit, mu=mu[,ee]);

    # Normalize toward a target function?
    if (!is.null(targetFcns)) {
      muT[ok,ee] <- targetFcns[[ee]](fl[ok]);
      if (.returnFit)
        fits[[ee]]$muT <- muT[,ee];
    }

    rm(fit, fl);

  } # for (ee ...)
  rm(hasFL, isSingleEnzymed);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate the *average* predicted signal across enzymes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sum on the non-log scale.
  if (.isLogged) {
    mu <- 2^mu;
    if (!is.null(targetFcns))
      muT <- 2^muT;
  }

  mu <- rowSums(mu, na.rm=TRUE);
#  mu <- mu / countFL; # Averaging (needed?!?)

  if (!is.null(targetFcns)) {
    muT <- rowSums(muT, na.rm=TRUE);
#    muT <- muT / countFL; # Averaging (needed?!?)
  }


  # Special case: Units with unknown fragment lengths
  if (onMissing != "ignore") {
    isMissing <- (countFL == 0);
    if (any(isMissing)) {
      if (onMissing == "median") {
        # Let the predicted value for these units be the robust average
        # of all other units (based on the assumption that the missing
        # fragment lengths are distributed as the known ones).
        mu[isMissing] <- median(mu[!isMissing], na.rm=TRUE);
        if (!is.null(targetFcns))
          muT[isMissing] <- median(muT[!isMissing], na.rm=TRUE);
      }
    }
  }

  if (.isLogged) {
    mu <- log2(mu);
    if (!is.null(targetFcns))
      muT <- log2(muT);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate the correction ("normalization") factor
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate correction factors
  if (is.null(targetFcns)) {
    dy <- mu;
  } else {
    dy <- (mu - muT);
  }
  rm(mu);
  if (!is.null(targetFcns))
    rm(muT);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Transform signals
  ok <- is.finite(dy) & okY;
  rm(okY);
  y[ok] <- y[ok] - dy[ok];
  
  if (.returnFit)
    attr(y, "modelFit") <- fits;

  y;
}, private=TRUE)


############################################################################
# HISTORY:
# 2008-09-10
# o Added argument 'onMissing' to normalizeFragmentLength() for specifying
#   how to normalize (if at all) data points for which the fragment lengths 
#   are unknown.
# 2008-05-10
# o BUG FIX: If the 'subsetToFit' was shorter than the number of data 
#   points, an exception was thrown.  The test was supposed to be assert
#   that the subset was not greater than the number of data points.
# 2008-04-14
# o Removed any usage of R.utils::Arguments.
# 2007-11-29
# o BUG FIX: The implemented multi-enzyme model was not the one in mind;
#   The correction for the multi-enzyme data points was not right.
#   Have now created an updated example that displays the normalized 
#   log-ratios (as a function of fragment length as well as they densities).
#   The example does also test the case for non-aliquot mixing proportions
#   between enzymes. This is actually automagically corrected for by the
#   way the model was set up, i.e. there is no need to estimate the 
#   mixing proportions.
# 2007-11-19
# o Added Rdoc examples. From these simulation examples, it looks like the
#   multi-enzyme normalization method works.
# o Updated normalizeFragmentLength() to handle multiple enzymes.
# 2006-11-28
# o Created.
############################################################################
