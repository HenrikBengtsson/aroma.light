# Version (development version)

## Documentation

 * Update redirecting and broken URLs.

 * Fix `R CMD check` notes on "Escaped LaTeX specials: \&".
 

# Version 3.32.0 [2023-04-25]

## Notes

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.18 for R (>= 4.4.0).


# Version 3.31.0 [2023-04-25]

## Notes

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.17 for R (>= 4.3.0).


# Version 3.30.0 [2023-04-25]

## Notes

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.17 for R (>= 4.3.0).


# Version 3.29.0 [2022-11-01]


# Version 3.28.0 [2022-11-01]


# Version 3.27.0 [2022-04-26]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.16 for R-devel.


# Version 3.26.0 [2022-04-26]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.15 for R (>= 4.2.0).


# Version 3.25.1 [2022-04-21]

## Bug Fixes

 * Fixed partial argument name used in `iwpca()`.


# Version 3.24.0 [2021-10-27]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.14 for R (>= 4.1.1).


# Version 3.23.1 [2021-08-19]

## Documentation

 * Update several citation URLs that were either broken or redirects
   elsewhere.


# Version 3.22.0 [2021-05-19]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.13 for R (>= 4.0.3).


# Version 3.20.0 [2020-10-27]

 * The version number was bumped for the Bioconductor release version, which is
   now Bioconductor 3.12 for R (>= 4.0.0).


# Version 3.18.0 [2020-04-27]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.11 for R (>= 3.6.1).


# Version 3.17.1 [2019-12-09]

## Code Refactoring

 * Now importing `throw()` from **R.oo** instead of **R.methodsS3**.


# Version 3.17.0 [2019-10-30]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.11 for R (>= 3.6.1).


# Version 3.16.0 [2019-10-30]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.10 for R (>= 3.6.1).


# Version 3.15.1 [2019-08-28]

## Bug Fixes

 * `wpca()` for matrices had a `length > 1 in coercion to logical` bug.


# Version 3.15.0 [2019-05-02]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.10 for R (>= 3.6.0).


# Version 3.14.0 [2019-05-02]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.9 for R (>= 3.6.0).


# Version 3.13.0 [2018-10-30]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.9 for R (>= 3.6.0).


# Version 3.12.0 [2018-10-30]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.8 for R (>= 3.5.0).


# Version 3.11.2 [2018-09-04]

## Code Refactoring

 * `fitPrincipalCurve()` now requires **princurve** (>= 2.1.2) and was
   updated to make use of new `principcal_curve` class instead of
   deprecated `principcal.curve` class.  This update "should not"
   affect the results, but see
   <https://github.com/rcannood/princurve/issues/8> for information of
   what has changed in the **princurve** package in this respect.


# Version 3.11.1 [2018-08-28]

 * Updated installation instructions in README.md.


# Version 3.11.0 [2018-04-30]

* The version number was bumped for the Bioconductor devel version,
  which is now Bioconductor 3.8 for R (>= 3.6.0).


# Version 3.10.0 [2018-04-30]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.7 for R (>= 3.5.0).
 

# Version 3.9.1 [2017-12-19]

## New Features

 * `robustSmoothSpline()` now supports using Tukey's biweight (in
   addition to already exising L1) estimators.  See argument `method`.
   Thanks to Aaron Lun at the Cancer Research UK Cambridge Institute
   for adding this feature.


# Version 3.9.0 [2017-10-30]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.7 for R (>= 3.5.0).


# Version 3.8.0 [2017-10-30]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.6 for R (>= 3.4.0).


# Version 3.7.0 [2017-04-23]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.6 for R (>= 3.4.0).


# Version 3.6.0 [2017-04-23]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.5 for R (>= 3.4.0).


# Version 3.5.1 [2017-04-14]

## Significant Changes

 * `robustSmoothSpline()` uses a re-weighted re-iterative algorithm
   that fits a smooth spline using `stats::smooth.spline()`,
   calculates the residuals and which are used to fit a re-weighted
   smooth spline and so on until converence. Due to updates to
   `stats::smooth.spline()` in R (>= 3.4.0) it is no longer feasible
   to maintain a highly optimized version of the algorithm, because it
   was based on internal `stats::smooth.spline()` code that has no
   completely changed.  Instead the re-iterative algorithm calls
   `stats::smooth.spline()` as is, which slows it down.  More
   importantly, it will now give slightly different estimates.

## Software Quality

 * In addition to continous integration (CI) tests and nightly
   Bioconductor tests, the package is now also tested regularly
   against all reverse package depencies available on CRAN and
   Bioconductor.


# Version 3.5.0 [2016-10-18]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.5 for R (>= 3.4.0).
 

# Version 3.4.0 [2016-10-18]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.4 for R (>= 3.3.1).
 

# Version: 3.3.0 [2016-05-03]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.4 for R (>= 3.3.0).
 

# Version 3.3.2 [2016-09-16]

## Code Refactoring

 * Using `NA_real_` (not just `NA`) everywhere applicable.

## Bug Fixes

 * `robustSmoothSpline()` gave an error since R-devel (>= 3.4.0
   r70682).


# Version 3.3.1 [2016-08-10]

## Code Refactoring

 * CLEANUP: Using `seq_len()` and `seq_along()` everywhere (Issue #8)


# Version 3.3.0 [2016-05-03]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.4 for R (>= 3.3.0).
 
 
# Version 3.2.0 [2016-05-03]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.3 for R (>= 3.3.0).
 
 

# Version 3.1.1 [2016-01-06]

 * Package requires R (>= 2.15.2).

 * CLEANUP: robustSmoothSpline() no longer generates messages that
   ".nknots.smspl() is now exported; use it instead of n.knots()" for
   R (>= 3.1.1).


# Version 3.1.0 [2015-10-23]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.3 for R (>= 3.3.0).


# Version 3.0.0 [2015-10-13]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.2 for R (>= 3.2.2).


# Version 2.99.0 [2015-10-06]

 * No changes.


# Version 2.9.0 [2015-09-17]

## Significant Changes

 * Argument `preserveScale` for `normalizeTumorBoost()` now defaults
   to FALSE. Since 1.33.3 (2014-04-30) it had no default and prior to
   that it was TRUE.


# Version 2.5.3 [2015-09-13]

## Software Quality

 * ROBUSTNESS: Explicitly importing core R functions.

## Bug Fixes

 * `rowAverages()` and `normalizeAverages()` would give an error if
   some of the argument default functions are overridden by
   non-functions of the same name in the calling environment.


# Version 2.5.2 [2015-06-16]

## Software Quality

 * Relaxed package test for `backtransformPrincipalCurve()`.


# Version 2.5.1 [2015-05-24]

 * Bumped package dependencies.

## Deprecated and Defunct

 * Removed obsolete `wpca(..., method = "dsvdc")`; was only needed for
   backward compatibility with R (< 1.7.0).


# Version 2.5.0 [2015-04-16]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.2 for R (>= 3.3.0).
 

# Version 2.4.0 [2015-04-16]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.1 for R (>= 3.2.0).
 

# Version 2.3.3 [2015-02-18]

## New Features

 * If a value of argument `xlim` or `ylim` for `plotDensity()` is NA,
   then it defaults to the corresponding extreme value of the data,
   e.g. `plotDensity(x, xlim = c(0, NA))`.


# Version 2.3.2 [2015-02-17]

## Software Quality

 * ROBUSTNESS: Added package tests. Code coverage is 76%.

## Code Refactoring

 * CLEANUP: Using `requestNamespace()` instead of `request()`.


# Version 2.3.1 [2014-12-08]

 * Same updates as in 2.2.1.
 
 
# Version 2.3.0 [2014-10-13]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 3.1 for R (>= 3.2.0).
 

# Version 2.2.1 [2014-12-08]

## Code Refactoring

 * Minor code cleanup.


# Version 2.2.0 [2014-10-13]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 3.0 for R (>= 3.1.1).
 

# Version 2.1.2 [2014-09-23]

 * Minor tweaks due to the move to GitHub.
 

# Version 2.1.1 [2014-09-16]

## Software Quality

 * Fixed some new `R CMD check` NOTEs.

## Code Refactoring

 * CLEANUP: Now importing **R.utils** (instead of only suggesting it).

 * IMPORTANT/CLEANUP: The **matrixStats** package is no longer
   attached with this package.  In other words, you now might have to
   add `library(matrixStats)` to your scripts.


# Version 2.1.0 [2014-04-11]

 * The version number was bumped for the Bioconductor devel version,
   which is now Bioconductor 2.15 for R (>= 3.1.0).
 

# Version 2.0.0 [2014-04-11]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 2.14 for R (>= 3.1.0).
 

# Version 1.99.3 [2014-03-31]

 * Bumped the version such that the next release will be 2.0.0.
 

# Version 1.33.3 [2014-03-30]

## Significant Changes

 * Argument `preserveScale` for `normalizeTumorBoost()` is now
   required. The goal with this is to in a future version migrate to
   use `preserveScale = FALSE` as the default (was `preserveScale =
   TRUE`) in order to avoid introducing a a global bias in the tumor
   allele B fraction of heterozygous SNPs.  The examples use
   `preserveScale = FALSE` now.

## New Features

 * Added `pairedAlleleSpecificCopyNumbers()`.


# Version 1.33.2 [2014-03-25]

## New Features

 * Now `plotDensity()` supports weights via argument `W`.


# Version 1.33.1 [2014-03-25]

## New Features

 * Now `plotDensity()` also supports `density()` objects.

## Code Refactoring

 * CLEANUP: `robustSmoothSpline()` no longer uses `DUP = FALSE` in an
   internal `.Fortran()` call.

 * Bumped up package dependencies.


# Version 1.33.0 [2013-10-14]

 * The version number was bumped for the Bioconductor devel version.
 
 
# Version 1.32.0 [2012-10-14]

 * The version number was bumped for the Bioconductor release version,
   which is now Bioconductor 2.13 for R (>= 3.0.0).


# Version 1.31.10 [2013-10-08]

## New Features

 * Added `averageQuantile()` for matrices in addition to lists.

## Performance

 * SPEEDUP: Now `normalizeQuantileSpline(..., sortTarget = TRUE)`
   sorts the target only once for lists of vectors just as done for
   matrices.

## Documentation

 * Merged the documentation for `normalizeQuantileSpline()` for all
   data types into one help page.  Same for `plotXYCurve()`.

## Code Refactoring

 * Bumped up package dependencies.

## Bug Fixes

 * Argument `lwd` of `plotXYCurve(X, ...)` was ignored if `X` was a
   matrix.


# Version 1.31.9 [2013-10-07]

## New Features

 * Now `library(aroma.light, quietly = TRUE)` attaches the package
   completely silently without any messages.

 * Now the `aroma.light` `Package` object is also available when the
   package is only loaded (but not attached).

## Documentation

 * Merged the documentation for `normalizeQuantileRank()` for numeric
   vectors and lists.

 * Now documention S3 methods under their corresponding generic
   function.


# Version 1.31.8 [2013-10-02]

## Documentation

 * More generic functions are now "aliased" under relevant
   corresponding methods.


# Version 1.31.7 [2013-09-27]

## Performance

 * SPEEDUP: Now all package functions utilizes **matrixStats**
   functions where possible, e.g. `anyMissing()`, `colMins()`, and
   `rowWeightedMedians()`.

## Code Refactoring

 * Bumped up package dependencies.


# Version 1.31.6 [2013-09-25]

## Code Refactoring

 * CLEANUP: Package no longer use a fallback attachment of the
   **R.oo** package upon attachment.


# Version 1.31.5 [2013-09-23]

## Software Quality

 * ROBUSTNESS: Now properly declaring all S3 methods in the NAMESPACE
   file.

## Performance

 * SPEEDUP/CLEANUP: `normalizeTumorBoost()` now uses `which()` instead
   of `whichVector()` of **R.utils**.  Before R (< 2.11.0), `which()`
   used to be 10x slower than `whichVector()`, but now it's 3x faster.

## Code Refactoring

 * CLEANUP: Now only using `Authors@R` in DESCRIPTION, which is
   possible since R (>= 2.14.0).  Hence the new requirement on the
   version of R.

 * Bumped up package dependencies.


# Version 1.31.4 [2013-09-10]

## Software Quality

 * CLEANUP: Now package explicitly imports what it needs from
   **matrixStats**.

## Code Refactoring

 * Bumped up package dependencies.


# Version 1.31.3 [2013-05-25]

## Performance

 * SPEEDUP: Removed all remaining `gc()` calls, which were in
   `normalizeQuantileSpline()`.

 * SPEEDUP: Replaced all `rm()` calls with NULL assignments.

## Code Refactoring

 * Updated the package dependencies.


# Version 1.31.2 [2013-05-20]

 * Same updates as in v1.30.2.
 

# Version 1.31.1 [2011-04-18]

 * Same updates as in v1.30.1.


# Version 1.31.0 [2013-04-03]

 * The version number was bumped for the Bioconductor devel version.


# Version 1.30.5 [2013-09-25]

## Software Quality

 * Backport from v1.31.5: Declaring all S3 methods in NAMESPACE.

## Code Refactoring

 * Backport from v1.31.5: `normalizeTumorBoost()` now uses `which()`,
   which also removes one dependency on **R.utils**.


# Version 1.30.4 [2013-09-25]

## Software Quality

 * Backport from v1.31.4: Now package explicitly imports what it needs
   from **matrixStats**.


# Version 1.30.3 [2013-09-25]

## Performance

 * Backport from v1.31.3: Removal of all `gc()` calls and removal of
   variables is now faster.

## Bug Fixes

 * Removed one stray `str()` debug output in `robustSmoothSpline()`.


# Version 1.30.2 [2013-05-20]

# Version 1.30.1 [2013-04-18]

## New Features

 * Now `backtransformPrincipalCurve()` preserves dimension names.

## Bug Fixes

 * `backtransformPrincipalCurve()` gave an error if the pricipal curve
   was fitted using data with missing values.

 * `fitPrincipalCurve()` would not preserve dimension names if data
   contain missing values.


# Version 1.30.0 [2012-04-03]

 * The version number was bumped for the Bioconductor release version,
   which now is Bioconductor 2.12 for R (>= 3.0.0).
 
 
# Version 1.29.0 [2012-10-01]

 * The version number was bumped for the Bioconductor devel version.
 
 
# Version 1.28.0 [2012-10-01]

 * The version number was bumped for the Bioconductor release version, which
   now is Bioconductor 2.11 for R (>= 2.15.1).
 
 
# Version 1.27.1 [2012-09-12]

## Software Quality

 * ROBUSTNESS: Replaced an `.Internal(psort(...))` call in
   `medianPolish()` with a call to **matrixStats**:::`.psortKM()`.


# Version 1.27.0 [2012-08-30]

## Code Refactoring

 * CLEANUP: Removed `weightedMedian()`, which has been moved to the
   **matrixStats** package.

 * BACKWARD COMPATIBILITY: Now package depends on the **matrixStats**
   (>= 0.5.2) package, so that `weightedMedian()` is still available
   when loading this package.  In future releases, **matrixStats**
   will be downgraded to only be a suggested package.


# Version 1.26.1 [2012-08-30]

## Code Refactoring

 * Updated the package dependencies.

## Bug Fixes

 * `robustSmoothSpline()` would not work with most recent R devel
   versions.


# Version 1.26.0 [2012-08-19]

## Significant Changes

 * Changed the license of aroma.light to GPL (>= 2) from LGPL (>= 2),
   because some of the implementation was adopted from GPL (>= 2)
   code, i.e. `robustSmoothSpline()` uses code from
   `stats::smooth.spline()`.

## Software Quality

 * `R CMD check` no longer warns about some examples depending on the
   **R.basic** package.


# Version 1.25.4 [2012-08-19]

## Software Quality

 * WORKAROUND: Now `robustSmoothSpline()` robustly locates the proper
   native R fit function for smooth splines, which vary with different
   releases of R.


# Version 1.25.3 [2012-04-16]

## Code Refactoring

 * Package no longer depends on **R.methodsS3**, only imports.


# Version 1.25.2 [2012-04-16]

## Software Quality

 * `R CMD check` no longer complaints about `.Internal()` calls.


# Version 1.25.1 [2012-04-16]

## New Features

 * Added support for `fitNaiveGenotypes(..., flavor = "fixed")`.

 * GENERALIZATION: Now `fitNaiveGenotypes()` returns also `flavor` and
   `tau`. The latter are the genotype thresholds used by the caller.


## Code Refactoring

 * CLEANUP: Dropped argument `flavor` of `callNaiveGenotypes()`; it is
   instead passed to `fitNaiveGenotypes()` via `...`.


# Version 1.25.0 [2012-03-30]

 * The version number was bumped for the Bioconductor devel version.
 

# Version 1.24.0 [2012-03-30]

 * The version number was bumped for the Bioconductor release version, which
   now is Bioconductor 2.10 for R (>= 2.15.0).


# Version 1.23.0 [2011-10-31]

 * The version number was bumped for the Bioconductor devel version.
 

# Version 1.22.0 [2011-10-31]

 * The version number was bumped for the Bioconductor release version,
   which now is Bioconductor 2.9 for R (>= 2.14.0).
 

# Version 1.21.2 [2011-10-10]

## Code Refactoring

 * Updated `robustSmoothSpline()` such that it works with the new
   "uniqueness" scheme of `smooth.spline()` in R v2.14.0 and newer.
   It is tricky, because `robustSmoothSpline()` is a reiterative
   algorithm which requires that the choosen "unique" x:s does not
   change in each iteration.  Previously, `signif(x, 6)` was used to
   identify unique x:s, which gives the same set of values when called
   twice, whereas this is not true for the new choice with `round((x -
   mean(x))/tol)`.


# Version 1.21.1 [2011-06-26]

## New Features

 * Added argument `aspectRatio` to `plotMvsA()`.  It can be used to
   adjust the range of the `Mlim` argument relative to the `Alim`
   argument.


# Version 1.21.0 [2011-04-13]

 * The version number was bumped for the Bioconductor devel version.
 

# Version 1.20.0 [2010-04-13]

 * The version number was bumped for the Bioconductor release version,
   which now is Bioconductor 2.8 for R (>= 2.13.0).
 

# Version 1.19.6 [2011-04-12]

## Code Refactoring

 * CLEANUP: Removed internal patch of `principal.curve()`.  If an
   older version than **princurve** v1.1-10 is used, an informative
   error is thrown requesting an update.  The internal patch is part
   of the offical **princurve** v1.1-10 release (since 2009-10-04).

 * Now all methods allocate objects with NAs of the appropriate mode.

## Known Issues

 * Recent updates to `smooth.spline()` in R v2.14.0 causes
   `robustSmoothSpline()` to break in some cases.


# Version 1.19.5 [2011-04-08]

## New Features

 * Now `weightedMedian()` returns NA:s of the same mode as argument
   `x`.


# Version 1.19.4 [2011-03-03]

# Version 1.19.3 [2011-02-05]

# Version 1.19.2 [2010-10-22]

# Version 1.19.1 [2010-10-18]

# Version 1.19.0 [2010-10-18]

 * The version number was bumped for the Bioconductor devel version.
 

# Version 1.18.4 [2011-03-03]

## Bug Fixes

 * `findPeaksAndValleys(x, to)` were `x` is numeric would use partial
   match and interpret `to` as argument `tol` and not part of `...`
   passed to `density()`. This problem was solved by placing `...`
   before argument `tol`.  Thanks Oscar Rueda at the Cancer Reasearch
   UK for reporting on and identifying this bug.


# Version 1.18.3 [2011-02-05]

## Documentation

 * Added paragraphs on how to do affine normalization when channel
   offsets are known/zero.  Same for multiscan calibration when
   scanner offsets are known/zero.

 * Fixed broken links to help for `iwpca()`.


# Version 1.18.2 [2010-10-22]

## Documentation

 * Minor clarifications to the help page on "1. Calibration and
   Normalization". This page is now also referenced in
   `help("calibrateMultiscan")`.


# Version 1.18.1 [2010-10-18]

## New Features

 * Argument `censorAt` for `fitNaiveGenotypes()` has new default.

 * These updates were supposed to be in v1.17.7, but we forgot to
   commit them to the Bioconductor repository before the new
   Bioconductor release.

## Bug Fixes

 * `fitNaiveGenotypes(..., subsetToFit = <logical>)` would throw an
   exception reporting `Some elements of argument `subsetToFit` is out
   of range ...`.


# Version 1.18.0 [2010-10-18]

 * The version number was bumped for the Bioconductor release version,
   which now is Bioconductor 2.7 for R (>= 2.12.0).
 

# Version 1.17.6 [2010-10-08]

## New Features

 * Now `findPeaksAndValleys()` returns a object of class
   `PeaksAndValleys`, which extends data.frame.


# Version 1.17.5 [2010-10-07]

## New Features

 * Added optional argument `fit` to `callNaiveGenotypes()` for passing
   a model fit returned by `fitNaiveGenotypes()`.  If not specified,
   `callNaiveGenotypes()` will call `fitNaiveGenotypes()` internally.

 * Added `fitNaiveGenotypes()`, which previously was only internal of
   `callNaiveGenotypes()`.


# Version 1.17.4 [2010-10-06]

## New Features

 * Added `findPeaksAndValleys()` for the `density` class, which then
   `findPeaksAndValleys()` for `numeric` utilizes.


# Version 1.17.3 [2010-09-18]

## Software Quality

 * ROBUSTNESS: Now `normalizeFragmentLength()` asserts that arguments
   `fragmentLengths` and `y` contain at least some finite values and
   specifies the same number of units.  In addition, the method also
   gives more informative error messages in case it cannot fit the
   normalization function due to non-finite values.


# Version 1.17.2 [2010-08-04]

## New Features

 * Added argument `preserveScale` to `normalizeTumorBoost()` to
   rescale the calibrated allele B fractions for heterozygous SNPs
   such that the compression relative to the homozgygotes is
   preserved.


# Version 1.17.1 [2010-07-23]

# Version 1.17.0 [2010-04-22]

 * The version number was bumped for the Bioconductor devel version.
 

# Version 1.16.1 [2010-07-23]

## New Features

 * Now `callNaiveGenotypes()` returns the model estimates as attribute
   `modelFit`.  This feature was supposed to be in v1.16.0.


# Version 1.16.0 [2010-04-22]

 * The version number was bumped for the Bioconductor release version,
   which now is Bioconductor 2.6 for R (>= 2.11.0).
 

# Version 1.15.4 [2010-04-08]

## Software Quality

 * R devel assumes ASCII encoding unless specified. Added explicit
   Latin-1 encoding to the DESCRIPTION file to `R CMD check` to pass.


# Version 1.15.3 [2010-04-04]

## New Features

 * Added `normalizeDifferencesToAverage()`, `normalizeTumorBoost()`,
   `callNaiveGenotypes()` and internal `findPeaksAndValleys()`, which
   all were moved from the **aroma.cn** package.


# Version 1.15.2 [2010-03-12]

## Bug Fixes

 * The example of `fitPrincipalCurve()` used a non-existing argument
   name in the calls to `substitute()`.  Thanks to Brian Ripley at
   University of Oxford for reporting this.


# Version 1.15.1 [2009-11-01]

## Code Refactoring

 * Now `fitPrincipalCurve()` only uses the internal bug-fix patch if a
   version earlier than **princurve** v1.1-10 is installed.


# Version 1.15.0 [2009-10-27]

 * The version number was bumped for the Bioconductor devel version.
 
 
# Version 1.14.0 [2009-10-27]

 * The version number was bumped for the Bioconductor release version,
   which now is Bioconductor 2.5 for R (>= 2.10.0).
 

# Version 1.13.6 [2009-10-20]

## Documentation

 * FIX: CITATION file reverted to that of v1.13.4.


# Version 1.13.5 [2009-10-08]

## Documentation

 * CITATION file (incorrectly) updated by Bioconductor maintainers.


# Version 1.13.4 [2009-09-23]

## Documentation

 * Fixed a few broken Rd links.


# Version 1.13.3 [2009-07-15]

## New Features

 * ADDED: `fitXYCurve()` and `backtransformXYCurve()`.

 * Added attribute `processingTime` to the fit object returned by
   `fitPrincipalCurve()`.


# Version 1.13.2 [2009-05-29]

# Version 1.13.1 [2009-05-13]

# Version 1.13.0 [2009-04-20]

# Version 1.12.2 [2009-05-29]

## Code Refactoring

 * Replacing old HOWTOCITE with a standard CITATION file.

## Bug Fixes

 * Previous bug fix in `backtransformPrincipalCurve()` regarding
   argument `dimension` broke the initial purpose of this
   argument. Since both use cases are still of interest, how the
   subsetting is done is now based on whether the number of dimensions
   of the input data and the model fit match or not. See help and
   example code for `backtransformPrincipalCurve.matrix()`.


# Version 1.12.1 [2009-05-13]

## Bug Fixes

 * `backtransformPrincipalCurve(..., dimensions)` did not subset the
   `X` matrix. Also, the method now returns a matrix of the same
   number of columns requested.  The Rd example now illustrates this.
   Thanks to Pierre Neuvial, UC Berkeley for the troublshooting and
   fix.


# Version 1.12.0 [2009-04-20]

 * The version number was bumped for the Bioconductor release version.
 
 
# Version 1.11.2 [2009-02-08]

## Bug Fixes

 * An error was thrown in `backtransformPrincipalCurve()` when
   argument `dimensions` was specified.


# Version 1.11.1 [2009-01-12]

## New Features

 * Added `fitPrincipalCurve()` and `backtransformPrincipalCurve()`.


# Version 1.11.0 [2008-10-21]

 * The version number was bumped for the Bioconductor devel version.
 
 
# Version 1.10.0 [2008-10-21]

 * The version number was bumped for the Bioconductor release version.
 
 
# Version 1.9.2 [2008-09-11]

## New Features

 * Added argument `onMissing` to `normalizeFragmentLength()` for
   specifying how to normalize, if at all, data points for which the
   fragment lengths are unknown.  For backward compatibility, we start
   off by having it `"ignore"` by default.

## Code Refactoring

 * MEMORY OPTIMIZATION: `robustSmoothSpline()` is now cleaning out
   more variables when no longer needed.


# Version 1.9.1 [2008-05-10]


# Version 1.9.0 [2008-04-29]

 * The version number was bumped for the Bioconductor devel version.
 
 
# Version 1.8.1 [2008-05-10]

## Bug Fixes

 * If the `subsetToFit` of `normalizeFragmentLength()` was shorter
   than the number of data points, an exception was thrown.  The test
   was supposed to assert that the subset was not greater than the
   number of data points.


# Version 1.8.0 [2008-04-29]

 * The version number was bumped for the Bioconductor release version.
 
 
# Version 1.7.2 [2008-04-14]

## New Features

 * Added `normalizeFragmentLength()`.

 * Added `normalizeQuantileSpline()`.

 * Renamed `normalizeQuantile()` to `normalizeQuantileRank()`.

 * Added `plotXYCurve()`.

 * Added `predict()` for the `lowess` class.


# Version 1.7.1 [2007-11-28]

## New Features

 * The startup message when loading the package is generated with
   `packageStartupMessage()` so that it can be suppressed.

## Documentation

 * TYPO: Corrected a spelling error in the help pages.

## Software Quality

 * Package passes `R CMD check` for R v2.6.1.

## Code Refactoring

 * Package now only suggest the **R.oo** package, and instead depends
   on the new **R.methodsS3**.


# Version 1.7.0 [2007-10-08]

 * The version number was bumped for the Bioconductor devel version.
 
 
# Version 1.6.0 [2007-10-08]

 * The version number was bumped for the Bioconductor release version.
 

# Version 1.5.2 [2007-08-10]

## Software Quality

 * Package pass `R CMD check` for R v2.6.0.


# Version 1.5.1 [NA]

## New Features

 * Added `normalizeAverage()`.


## Software Quality

 * Package pass `R CMD check` for R v2.6.0 with Rd encoding errors.


# Version 1.5.0 [2007-05-09]

 * The version number was bumped for the Bioconductor devel version.
 

# Version 1.4.0 [2007-05-09]

 * The version number was bumped up with the Bioconductor release.
 

# Version 1.3.1 [2007-01-15]

## Code Refactoring

 * Removed code to use **modreg** for backward compatibility with R <
   1.9.0.

 * Added **R.utils** to Suggests field of DESCRIPTION.


# Version 1.3.0 [2006-10-03]

 * The devel version number was bumped up with the Bioconductor
   release.
 

# Version 1.2.0 [2006-10-03]

 * The version number was bumped up with the Bioconductor release.


# Version 1.1.0 [2006-07-20]

## Significant Changes

 * Added to Bioconductor 1.9.

## New Features

 * Added some trial RSP pages. Try `browseRsp()` in the **R.rsp**
   package.


# Version 0.1.7 [2006-06-27]

## Code Refactoring

 * Made the package truely standalone except from **R.oo**.
   Previously package **R.basic** was used in some of the examples.


# Version 0.1.6 [2006-05-22]

## New Features

 * Added `medianPolish()` which is much faster than
   `stats::medpolish()` when there are no NA.

 * Added `plotDensity()` for list of vectors as well as for matrices.

 * Added `normalizeQuantile()` for lists of vectors as well as for a
   single vector of numerics.  To calculate the target quantile there
   is a new function `averageQuantile()`, which is also for lists of
   vectors.  It latter does not support robust estimatation of the
   average, because it safes memory.

 * Updated `normalizeQuantile()` for matrices according to the updates
   in the **limma** package.

## Code Refactoring

 * Added a namespace for the package.

 * Added `biocViews` since the package will eventually be added to the
   Bioconductor project.


# Version 0.1.5 [2006-04-21]

## Performance

 * Minor speedup to `weightedMedian()`, e.g. negative weights do no
   longer give and error, but are treated as zero instead.  This
   removes some overhead of the function.  Also, if it is known that
   there are no NAs that can be specified by `na.rm = NA`, which will
   skip NA checks.


# Version 0.1.4 [2006-03-28]

## Documentation

 * Updated broken Rd links.

 * Updated the references to publications.


# Version 0.1.3 [2006-01-22]

## New Features

 * Now `fitIWPCA()` does not return the data matrix. This is to save
   memory. The calling algorithm can equally well add the data if it
   is needed.

## Documentation

 * Added help on the returned parameters of `fitIWPCA()`.


# Version 0.1.2 [2005-09-06]

## New Features

 * All plot methods displaying log-ratios now assures that no fake
   log-ratios are calculated due to two negative raw signals.
   Similarily, methods display log-intensities now assures that the
   log-intensities are calculated as doubles to avoid possible
   overflow warnings for too large integers.


# Version 0.1.1 [2005-07-26]

## New Features

 * Added `sampleCorrelations()` and `sampleTuples()`.

 * Now argument `interpolate` of `weightedMedian()` defaults to TRUE
   only if `ties` is NULL.


# Version 0.1.0 [2005-06-03]

## Significant Changes

 * Created. Most of the matrix methods were copied from the
   **R.basic** and the **aroma** packages. The purpose of this package
   is to provide a standalone package, which does not require any of
   the aroma classes. This will allow the methods to be used by end
   users as is, or be utilized in other packages.
