# apmsWAPP

Version: 1.0

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

# aroma.core

Version: 3.1.3

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      ‘sfit’ ‘expectile’ ‘HaarSeg’ ‘mpcbs’
    ```

# EDASeq

Version: 2.14.1

## In both

*   checking examples ... ERROR
    ```
    Running examples in ‘EDASeq-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: getGeneLengthAndGCContent
    > ### Title: Get gene length and GC-content
    > ### Aliases: getGeneLengthAndGCContent
    > 
    > ### ** Examples
    > 
    > getGeneLengthAndGCContent("ENSG00000012048", "hsa")
    Connecting to BioMart ...
    Request to BioMart web service failed.
    The BioMart web service you're accessing may be down.
    Check the following URL and see if this website is available:
    http://www.ensembl.org:80/biomart/martservice?type=registry&requestid=biomaRt
    Error in if (!grepl(x = registry, pattern = "^\n*<MartRegistry>")) { : 
      argument is of length zero
    Calls: getGeneLengthAndGCContent -> useMart -> listMarts
    Execution halted
    ```

*   checking tests ...
    ```
     ERROR
    Running the tests in ‘tests/getLengthAndGC.R’ failed.
    Last 13 lines of output:
      
          aperm, apply
      
      > library(yeastRNASeq)
      > 
      > getGeneLengthAndGCContent(id=c("ENSG00000012048", "ENSG00000139618"), org="hsa", mode = "biomart")
      Connecting to BioMart ...
      Request to BioMart web service failed.
      The BioMart web service you're accessing may be down.
      Check the following URL and see if this website is available:
      http://www.ensembl.org:80/biomart/martservice?type=registry&requestid=biomaRt
      Error in if (!grepl(x = registry, pattern = "^\n*<MartRegistry>")) { : 
        argument is of length zero
      Calls: getGeneLengthAndGCContent -> useMart -> listMarts
      Execution halted
    ```

*   checking re-building of vignette outputs ... WARNING
    ```
    ...
    locfit 1.5-9.1 	 2013-03-22
    
    Attaching package: 'locfit'
    
    The following objects are masked from 'package:ShortRead':
    
        left, right
    
    Loading required package: lattice
        Welcome to 'DESeq'. For improved performance,
        usability and functionality, please consider
        migrating to 'DESeq2'.
    Connecting to BioMart ...
    Request to BioMart web service failed.
    The BioMart web service you're accessing may be down.
    Check the following URL and see if this website is available:
    http://www.ensembl.org:80/biomart/martservice?type=registry&requestid=biomaRt
    Quitting from lines 554-555 (EDASeq.Rnw) 
    Error: processing vignette 'EDASeq.Rnw' failed with diagnostics:
    argument is of length zero
    Execution halted
    ```

*   checking dependencies in R code ... NOTE
    ```
    Namespace in Imports field not imported from: ‘methods’
      All declared Imports should be used.
    ```

*   checking R code for possible problems ... NOTE
    ```
    ...
    plotQuality,FastqFileList : <anonymous>: no visible global function
      definition for ‘as’
    plotQuality,FastqFileList: no visible global function definition for
      ‘lines’
    plotRLE,matrix: no visible global function definition for ‘abline’
    updateObject,SeqExpressionSet: no visible global function definition
      for ‘callNextMethod’
    Undefined global functions or variables:
      abline as available.packages biocLite callNextMethod countBam
      elementMetadata legend lines loess lowess matplot median narrow new
      pairs points predict quality quantile rainbow smoothScatter text
      validObject
    Consider adding
      importFrom("grDevices", "rainbow")
      importFrom("graphics", "abline", "legend", "lines", "matplot", "pairs",
                 "points", "smoothScatter", "text")
      importFrom("methods", "as", "callNextMethod", "new", "validObject")
      importFrom("stats", "loess", "lowess", "median", "predict", "quantile")
      importFrom("utils", "available.packages")
    to your NAMESPACE file (and ensure that your DESCRIPTION Imports field
    contains 'methods').
    ```

# TIN

Version: 1.12.0

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

