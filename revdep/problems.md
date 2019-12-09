# apmsWAPP

<details>

* Version: 1.0
* Source code: https://github.com/cran/apmsWAPP
* Date/Publication: 2014-04-22 14:39:54
* Number of recursive dependencies: 51

Run `revdep_details(,"apmsWAPP")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    ...
    TSPM: no visible global function definition for ‘p.adjust’
    norm.inttable: no visible global function definition for ‘median’
    norm.inttable : <anonymous>: no visible global function definition for
      ‘quantile’
    saint_permF: no visible global function definition for ‘read.table’
    saint_permF: no visible global function definition for ‘write.table’
    saint_permF: no visible global function definition for ‘write.csv2’
    tspm_apms: no visible global function definition for ‘read.table’
    varFilter : <anonymous>: no visible global function definition for
      ‘median’
    varFilter: no visible global function definition for ‘median’
    varFilter: no visible global function definition for ‘quantile’
    Undefined global functions or variables:
      deviance glm hatvalues median p.adjust pchisq pf poisson qchisq qf
      quantile read.table residuals write.csv2 write.table
    Consider adding
      importFrom("stats", "deviance", "glm", "hatvalues", "median",
                 "p.adjust", "pchisq", "pf", "poisson", "qchisq", "qf",
                 "quantile", "residuals")
      importFrom("utils", "read.table", "write.csv2", "write.table")
    to your NAMESPACE file.
    ```

# aroma.affymetrix

<details>

* Version: 3.2.0
* Source code: https://github.com/cran/aroma.affymetrix
* URL: https://www.aroma-project.org/, https://github.com/HenrikBengtsson/aroma.affymetrix
* BugReports: https://github.com/HenrikBengtsson/aroma.affymetrix/issues
* Date/Publication: 2019-06-23 06:00:14 UTC
* Number of recursive dependencies: 80

Run `revdep_details(,"aroma.affymetrix")` for more info

</details>

## In both

*   checking installed package size ... NOTE
    ```
      installed size is  5.8Mb
      sub-directories of 1Mb or more:
        R             2.3Mb
        help          1.1Mb
        testScripts   1.3Mb
    ```

*   checking dependencies in R code ... NOTE
    ```
    Error in setGeneric("getX", function(object, type) standardGeneric("getX")) : 
      could not find function "setGeneric"
    ```

# aroma.core

<details>

* Version: 3.2.0
* Source code: https://github.com/cran/aroma.core
* URL: https://github.com/HenrikBengtsson/aroma.core, https://www.aroma-project.org/
* BugReports: https://github.com/HenrikBengtsson/aroma.core/issues
* Date/Publication: 2019-06-17 18:20:03 UTC
* Number of recursive dependencies: 46

Run `revdep_details(,"aroma.core")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'sfit', 'expectile', 'HaarSeg', 'mpcbs'
    ```

# EDASeq

<details>

* Version: 2.20.0
* Source code: https://github.com/cran/EDASeq
* URL: https://github.com/drisso/EDASeq
* BugReports: https://github.com/drisso/EDASeq/issues
* Date/Publication: 2019-10-29
* Number of recursive dependencies: 121

Run `revdep_details(,"EDASeq")` for more info

</details>

## In both

*   checking whether the namespace can be loaded with stated dependencies ... WARNING
    ```
    Error in setGeneric("Makesense", function(expr, lib, ...) standardGeneric("Makesense")) : 
      could not find function "setGeneric"
    Error: unable to load R code in package ‘geneplotter’
    Execution halted
    
    A namespace must be able to be loaded with just the base namespace
    loaded: otherwise if the namespace gets loaded by a saved object, the
    session will be unable to start.
    
    Probably some imports need to be declared in the NAMESPACE file.
    ```

*   checking examples ... WARNING
    ```
    Found the following significant warnings:
    
      Warning: 'GenomicRangesList' is deprecated.
    Deprecated functions may be defunct as soon as of the next release of
    R.
    See ?Deprecated.
    ```

*   checking R code for possible problems ... NOTE
    ```
    Error in setGeneric("Makesense", function(expr, lib, ...) standardGeneric("Makesense")) : 
      could not find function "setGeneric"
    Error: unable to load R code in package ‘geneplotter’
    Execution halted
    ```

# scone

<details>

* Version: 1.10.0
* Source code: https://github.com/cran/scone
* BugReports: https://github.com/YosefLab/scone/issues
* Date/Publication: 2019-10-29
* Number of recursive dependencies: 218

Run `revdep_details(,"scone")` for more info

</details>

## In both

*   checking whether the namespace can be loaded with stated dependencies ... WARNING
    ```
    Error in setGeneric("Makesense", function(expr, lib, ...) standardGeneric("Makesense")) : 
      could not find function "setGeneric"
    Error: unable to load R code in package ‘geneplotter’
    Execution halted
    
    A namespace must be able to be loaded with just the base namespace
    loaded: otherwise if the namespace gets loaded by a saved object, the
    session will be unable to start.
    
    Probably some imports need to be declared in the NAMESPACE file.
    ```

*   checking R code for possible problems ... NOTE
    ```
    Error in setGeneric("Makesense", function(expr, lib, ...) standardGeneric("Makesense")) : 
      could not find function "setGeneric"
    Error: unable to load R code in package ‘geneplotter’
    Execution halted
    ```

# scran

<details>

* Version: 1.14.5
* Source code: https://github.com/cran/scran
* Date/Publication: 2019-11-19
* Number of recursive dependencies: 188

Run `revdep_details(,"scran")` for more info

</details>

## In both

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/testthat.R’ failed.
    Last 13 lines of output:
      ══ testthat results  ═══════════════════════════════════════════════════════════
      [ OK: 5478 | SKIPPED: 0 | WARNINGS: 843 | FAILED: 9 ]
      1. Error: findMarkers works correctly with subsetting and spikes (@test-markers.R#43) 
      2. Error: pairwiseTTests works as expected with a design matrix (@test-pairwise-t.R#367) 
      3. Error: pairwiseTTests with linear models works across multiple cores (@test-pairwise-t.R#388) 
      4. Error: pairwiseTTests with linear models responds to non-standard level ordering (@test-pairwise-t.R#401) 
      5. Error: pairwiseTTests with design matrices responds to restriction (@test-pairwise-t.R#428) 
      6. Error: pairwiseTTests behaves as expected with subsetting (@test-pairwise-t.R#450) 
      7. Error: pairwiseTTests behaves as expected with log-transformation (@test-pairwise-t.R#497) 
      8. Error: pairwiseTTests behaves with standardization of the log-fold changes (@test-pairwise-t.R#526) 
      9. Error: pairwiseTTests fails gracefully with silly inputs (@test-pairwise-t.R#553) 
      
      Error: testthat unit tests failed
      Execution halted
      Error while shutting down parallel: unable to terminate some child processes
    ```

*   checking examples ... WARNING
    ```
    ...
      Warning: 'normalizeSCE' is deprecated.
      Warning: 'centreSizeFactors' is deprecated.
      Warning: handling of spike-ins via 'isSpike()' is deprecated.
      Warning: 'overlapExprs' is deprecated.
      Warning: 'parallelPCA' is deprecated.
      Warning: 'technicalCV2' is deprecated.
      Warning: 'isSpike<-' is deprecated.
      Warning: 'spikeNames' is deprecated.
      Warning: 'isSpike' is deprecated.
      Warning: 'technicalCV2' is deprecated.
      Warning: 'testVar' is deprecated.
      Warning: 'testVar' is deprecated.
      Warning: 'testVar' is deprecated.
      Warning: handling of spike-ins via 'isSpike()' is deprecated.
      Warning: 'normalizeSCE' is deprecated.
      Warning: 'centreSizeFactors' is deprecated.
      Warning: 'trendVar' is deprecated.
      Warning: 'trendVar' is deprecated.
    Deprecated functions may be defunct as soon as of the next release of
    R.
    See ?Deprecated.
    ```

*   checking installed package size ... NOTE
    ```
      installed size is 25.2Mb
      sub-directories of 1Mb or more:
        doc    2.1Mb
        libs  22.0Mb
    ```

# TIN

<details>

* Version: 1.18.0
* Source code: https://github.com/cran/TIN
* Date/Publication: 2019-10-29
* Number of recursive dependencies: 117

Run `revdep_details(,"TIN")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    ...
    scatterPlot: no visible global function definition for ‘pdf’
    scatterPlot: no visible global function definition for ‘bmp’
    scatterPlot: no visible global function definition for ‘plot’
    scatterPlot: no visible global function definition for ‘ave’
    scatterPlot: no visible global function definition for ‘axis’
    scatterPlot: no visible global function definition for ‘text’
    scatterPlot: no visible global function definition for ‘mtext’
    scatterPlot: no visible global function definition for ‘points’
    scatterPlot: no visible global function definition for ‘dev.off’
    Undefined global functions or variables:
      ave axis bmp colorRampPalette data dev.off dist hclust hist jpeg
      median mtext par pdf plot png points postscript quantile read.table
      text
    Consider adding
      importFrom("grDevices", "bmp", "colorRampPalette", "dev.off", "jpeg",
                 "pdf", "png", "postscript")
      importFrom("graphics", "axis", "hist", "mtext", "par", "plot",
                 "points", "text")
      importFrom("stats", "ave", "dist", "hclust", "median", "quantile")
      importFrom("utils", "data", "read.table")
    to your NAMESPACE file.
    ```

