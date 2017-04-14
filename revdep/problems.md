# Setup

## Platform

|setting  |value                        |
|:--------|:----------------------------|
|version  |R version 3.3.3 (2017-03-06) |
|system   |x86_64, linux-gnu            |
|ui       |X11                          |
|language |en                           |
|collate  |en_US.UTF-8                  |
|tz       |America/Los_Angeles          |
|date     |2017-04-14                   |

## Packages

|package     |*  |version |date       |source         |
|:-----------|:--|:-------|:----------|:--------------|
|aroma.light |   |3.4.0   |2017-03-10 |cran (@3.4.0)  |
|matrixStats |   |0.52.2  |2017-04-14 |cran (@0.52.2) |
|princurve   |   |1.1-12  |2013-04-25 |cran (@1.1-12) |
|R.methodsS3 |   |1.7.1   |2016-02-16 |cran (@1.7.1)  |
|R.oo        |   |1.21.0  |2016-11-01 |cran (@1.21.0) |
|R.utils     |   |2.5.0   |2016-11-07 |cran (@2.5.0)  |

# Check results

2 packages with problems

|package       |version | errors| warnings| notes|
|:-------------|:-------|------:|--------:|-----:|
|metaseqR      |1.14.0  |      1|        1|     4|
|oneChannelGUI |1.40.0  |      1|        0|     0|

## metaseqR (1.14.0)
Maintainer: Panagiotis Moulos <moulos@fleming.gr>

1 error  | 1 warning  | 4 notes

```
checking tests ... ERROR
  Running ‘runTests.R’ [12s/34s]
Running the tests in ‘tests/runTests.R’ failed.
Last 13 lines of output:
  ERROR in test_estimate_aufc_weights: Error in .check_ncores(cores) : 4 simultaneous processes spawned
  ERROR in test_metaseqr: Error in .check_ncores(cores) : 3 simultaneous processes spawned
  
  Test files with failing tests
  
     test_estimate_aufc_weights.R 
       test_estimate_aufc_weights 
  
     test_metaseqr.R 
       test_metaseqr 
  
  
  Error in BiocGenerics:::testPackage("metaseqR") : 
    unit tests failed for package metaseqR
  Execution halted

checking re-building of vignette outputs ... WARNING
Error in re-building vignettes:
  ...

The following objects are masked from 'package:ShortRead':

    left, right

Loading required package: lattice
    Welcome to 'DESeq'. For improved performance, usability and
... 8 lines ...
    plotMA

The following object is masked from 'package:BiocGenerics':

    plotMA

Loading required package: qvalue
Quitting from lines 119-159 (metaseqr-pdf.Rnw) 
Error: processing vignette 'metaseqr-pdf.Rnw' failed with diagnostics:
4 simultaneous processes spawned
Execution halted

checking package dependencies ... NOTE
Packages which this enhances but not available for checking: ‘TCC’ ‘RMySQL’

checking DESCRIPTION meta-information ... NOTE
Malformed Title field: should not end in a period.

checking dependencies in R code ... NOTE
'library' or 'require' calls in package code:
  ‘BSgenome’ ‘BiocInstaller’ ‘GenomicRanges’ ‘RMySQL’ ‘RSQLite’
  ‘Rsamtools’ ‘TCC’ ‘VennDiagram’ ‘parallel’ ‘rtracklayer’ ‘survcomp’
  ‘zoo’
  Please use :: or requireNamespace() instead.
  See section 'Suggested packages' in the 'Writing R Extensions' manual.

checking R code for possible problems ... NOTE
biasPlotToJSON: no visible binding for global variable ‘nams’
cddat: no visible global function definition for ‘assayData’
cddat: no visible global function definition for ‘ks.test’
cddat: no visible global function definition for ‘p.adjust’
cdplot: no visible global function definition for ‘plot’
cdplot: no visible global function definition for ‘lines’
countsBioToJSON: no visible binding for global variable ‘nams’
diagplot.avg.ftd : <anonymous>: no visible binding for global variable
  ‘sd’
... 246 lines ...
             "dev.off", "jpeg", "pdf", "png", "postscript", "tiff")
  importFrom("graphics", "abline", "arrows", "axis", "grid", "lines",
             "mtext", "par", "plot", "plot.new", "plot.window", "points",
             "text", "title")
  importFrom("methods", "as", "new")
  importFrom("stats", "as.dist", "cmdscale", "cor", "end", "ks.test",
             "mad", "median", "model.matrix", "na.exclude", "optimize",
             "p.adjust", "p.adjust.methods", "pchisq", "quantile",
             "rexp", "rnbinom", "runif", "sd", "start", "var")
to your NAMESPACE file (and ensure that your DESCRIPTION Imports field
contains 'methods').
```

## oneChannelGUI (1.40.0)
Maintainer: Raffaele A Calogero <raffaele.calogero@unito.it>

1 error  | 0 warnings | 0 notes

```
checking package dependencies ... ERROR
Packages required but not available: ‘affylmGUI’ ‘tkrplot’

Depends: includes the non-default packages:
  ‘Biobase’ ‘affylmGUI’ ‘tkrplot’ ‘tkWidgets’ ‘IRanges’ ‘Rsamtools’
  ‘Biostrings’ ‘siggenes’ ‘chimera’
Adding so many packages to the search path is excessive and importing
selectively is preferable.

See section ‘The DESCRIPTION file’ in the ‘Writing R Extensions’
manual.
```

