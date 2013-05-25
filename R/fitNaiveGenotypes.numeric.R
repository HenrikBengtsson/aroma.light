###########################################################################/**
# @set "class=numeric"
# @RdocMethod fitNaiveGenotypes
# @alias fitNaiveGenotypes
#
# @title "Fit naive genotype model from a normal sample"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{y}{A @numeric @vector of length J containing allele B fractions
#    for a normal sample.}
#  \item{cn}{An optional @numeric @vector of length J specifying the true
#    total copy number in \eqn{\{0,1,2,NA\}} at each locus.  This can be
#    used to specify which loci are diploid and which are not, e.g.
#    autosomal and sex chromosome copy numbers.}
#  \item{subsetToFit}{An optional @integer or @logical @vector specifying
#    which loci should be used for estimating the model.
#    If @NULL, all loci are used.}
#  \item{flavor}{A @character string specifying the type of algorithm used.}
#  \item{adjust}{A postive @double specifying the amount smoothing for
#    the empirical density estimator.}
#  \item{...}{Additional arguments passed to @see "findPeaksAndValleys".}
#  \item{censorAt}{A @double @vector of length two specifying the range
#    for which values are considered finite.  Values below (above) this
#    range are treated as -@Inf (+@Inf).}
#  \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a @list of @lists.
# }
#
# @author
#
# \seealso{
#   To call genotypes see @seemethod "callNaiveGenotypes".
#   Internally @see "findPeaksAndValleys" is used to identify the thresholds.
# }
#*/###########################################################################
setMethodS3("fitNaiveGenotypes", "numeric", function(y, cn=rep(2L, length(y)), subsetToFit=NULL, flavor=c("density", "fixed"), adjust=1.5, ..., censorAt=c(-0.1,1.1), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'y':
  J <- length(y);
  y <- as.double(y);

  # Argument 'cn':
  cn <- as.integer(cn);
  if (length(cn) == 1) {
    cn <- rep(cn, J);
  } else if (length(cn) != J) {
    stop("The length of argument 'cn' does not match 'y': ",
                                            length(cn), " != ", J);
  }
  uniqueCNs <- sort(unique(cn));
  unknown <- which(!is.element(uniqueCNs, c(0,1,2,NA)));
  if (length(unknown) > 0) {
    unknown <- paste(uniqueCNs[unknown], collapse=", ");
    stop("Argument 'cn' contains unknown CN levels: ", unknown);
  }

  # Argument 'subsetToFit':
  if (!is.null(subsetToFit)) {
    if (is.logical(subsetToFit)) {
      if (length(subsetToFit) != J) {
        stop("The length of argument 'subsetToFit' does not match 'y': ",
                                         length(subsetToFit), " != ", J);
      }
      subsetToFit <- which(subsetToFit);
    } else {
      subsetToFit <- as.integer(subsetToFit);
      subsetToFit <- sort(unique(subsetToFit));
      if (!all(1 <= subsetToFit & subsetToFit <= J)) {
        stop(sprintf("Some elements of argument 'subsetToFit' is out of range [1,%d].", J));
      }
    }
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'adjust':
  adjust <- as.double(adjust);
  if (length(adjust) != 1) {
    stop("Argument 'adjust' must be single value: ", adjust);
  }
  if (adjust <= 0) {
    stop("Argument 'adjust' must be positive: ", adjust);
  }

##   # Argument 'tol':
##   tol <- as.double(tol);
##   if (length(tol) != 1) {
##     stop("Argument 'tol' must be single value: ", tol);
##   }
##   if (tol <= 0) {
##     stop("Argument 'tol' must be positive: ", tol);
##   }

  # Argument 'censorAt':
  censorAt <- as.double(censorAt);
  stopifnot(length(censorAt) == 2);
  stopifnot(censorAt[1] <= censorAt[2]);

  # Argument 'verbose':
  if (inherits(verbose, "Verbose")) {
  } else if (is.numeric(verbose)) {
    require("R.utils") || throw("Package not available: R.utils");
    verbose <- Verbose(threshold=verbose);
  } else {
    verbose <- as.logical(verbose);
    if (verbose) {
      require("R.utils") || throw("Package not available: R.utils");
      verbose <- Verbose(threshold=-1);
    }
  }
  if (verbose && inherits(verbose, "Verbose")) {
    cat <- R.utils::cat;
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Fitting naive genotype model from normal allele B fractions (BAFs)");
  verbose && cat(verbose, "Flavor: ", flavor);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Adjust signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Censoring BAFs");
  verbose && cat(verbose, "Before:");
  verbose && summary(verbose, y);
  verbose && print(verbose, sum(is.finite(y)));
  # Censor values
  y[y < censorAt[1]] <- -Inf;
  y[y > censorAt[2]] <- +Inf;
  verbose && cat(verbose, "After:");
  verbose && summary(verbose, y);
  verbose && print(verbose, sum(is.finite(y)));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subsetting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(subsetToFit)) {
    verbose && enter(verbose, "Subsetting");
    verbose && cat(verbose, "Number of data points before: ", length(y));
    verbose && cat(verbose, "Number of true copy-number levels before: ", length(uniqueCNs));
    y <- y[subsetToFit];
    cn <- cn[subsetToFit];
    uniqueCNs <- sort(unique(cn));
    verbose && cat(verbose, "Number of data points afterward: ", length(y));
    verbose && cat(verbose, "Number of true copy-number levels afterward: ", length(uniqueCNs));
    verbose && exit(verbose);
  }

  # To please R CMD check
  type <- NULL; rm(list="type");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call genotypes
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitList <- list();
  for (kk in seq(along=uniqueCNs)) {
    cnKK <- uniqueCNs[kk];
    verbose && enter(verbose, sprintf("Copy number level #%d (C=%g) of %d", kk, cnKK, length(uniqueCNs)));

    keep <- which(cn == cnKK);
    yKK <- y[keep];

    # Exclude missing and non-finited values when fitting the density
    yT <- yKK[is.finite(yKK)];
    n <- length(yT);

    if (flavor == "density") {
      fit <- findPeaksAndValleys(yT, adjust=adjust, ...);
      verbose && cat(verbose, "Identified extreme points in density of BAF:");
      verbose && print(verbose, fit);

      fitValleys <- subset(fit, type == "valley");
      nbrOfGenotypeGroups <- nrow(fitValleys) + 1L;
      verbose && cat(verbose, "Local minimas (\"valleys\") in BAF:");
      verbose && print(verbose, fitValleys);
      tau <- fitValleys$x;
    } else if (flavor == "fixed") {
      args <- list(...);
      tau <- args$tau;
      if (is.null(tau)) {
        tau <- seq(length=cnKK) / (cnKK + 1L);
      }
      nbrOfGenotypeGroups <- length(tau) + 1L;
    }

    # Sanity check
    stopifnot(length(tau) == nbrOfGenotypeGroups - 1L);

    # Store
    fitKK <- list(
      flavor = flavor,
      cn=cnKK,
      nbrOfGenotypeGroups=nbrOfGenotypeGroups,  # Not really used
      tau=tau,
      n=n
    );

    if (flavor == "density") {
      fitKK$fit <- fit;
      fitKK$fitValleys <- fitValleys;
    }

    fitList[[kk]] <- fitKK;

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);

  class(fitList) <- c("NaiveGenotypeModelFit", class(fitList));

  fitList;
}) # fitNaiveGenotypes()


###########################################################################
# HISTORY:
# 2012-04-16
# o Added support for fitNaiveGenotypes(..., flavor="fixed").
# o GENERALIZATION: Now fitNaiveGenotypes() returns also 'flavor' and
#   'tau'.  The latter are the genotype threshholds used by the caller.
# 2010-10-14
# o TYPO FIX: Used name 'fitPeaks' instead of 'fitValleys'.
# 2010-10-12
# o New default of argument 'censorAt' for fitNaiveGenotypes().
# o BUG FIX: fitNaiveGenotypes(..., subsetToFit=<logical>) would throw
#   an exception reporting "Some elements of argument 'subsetToFit' is
#   out of range ...".
# 2010-10-07
# o Created from callNaiveGenotypes.R.
###########################################################################
