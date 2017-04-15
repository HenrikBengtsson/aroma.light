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

1 packages with problems

|package       |version | errors| warnings| notes|
|:-------------|:-------|------:|--------:|-----:|
|oneChannelGUI |1.40.0  |      0|        2|     6|

## oneChannelGUI (1.40.0)
Maintainer: Raffaele A Calogero <raffaele.calogero@unito.it>

0 errors | 2 warnings | 6 notes

```
checking whether package ‘oneChannelGUI’ can be installed ... WARNING
Found the following significant warnings:
  Warning: no DISPLAY variable so Tk is not available
  Warning: loading Rplot failed
See ‘/home/hb/repositories/aroma.light/revdep/checks/oneChannelGUI.Rcheck/00install.out’ for details.

checking sizes of PDF files under ‘inst/doc’ ... WARNING
  ‘gs+qpdf’ made some significant size reductions:
     compacted ‘Exon-level.analysis.pdf’ from 1395Kb to 598Kb
     compacted ‘RNAseq.pdf’ from 1979Kb to 385Kb
  consider running tools::compactPDF(gs_quality = "ebook") on these files

checking package dependencies ... NOTE
Depends: includes the non-default packages:
  ‘Biobase’ ‘affylmGUI’ ‘tkrplot’ ‘tkWidgets’ ‘IRanges’ ‘Rsamtools’
  ‘Biostrings’ ‘siggenes’ ‘chimera’
Adding so many packages to the search path is excessive and importing
selectively is preferable.

checking installed package size ... NOTE
  installed size is  6.6Mb
  sub-directories of 1Mb or more:
    doc   5.0Mb

checking DESCRIPTION meta-information ... NOTE
Malformed Description field: should contain one or more complete sentences.

checking dependencies in R code ... NOTE
'library' or 'require' calls in package code:
  ‘BSgenome.Hsapiens.UCSC.hg19’ ‘BSgenome.Mmusculus.UCSC.mm9’
  ‘BSgenome.Rnorvegicus.UCSC.rn4’ ‘Genominator’ ‘affy’ ‘affyPLM’
  ‘chipseq’ ‘maSigPro’
  Please use :: or requireNamespace() instead.
  See section 'Suggested packages' in the 'Writing R Extensions' manual.

checking R code for possible problems ... NOTE
.consistentFilter: warning in get("midas.p.Available", env =
  affylmGUIenvironment): partial argument match of 'env' to 'envir'
.consistentFilter: warning in get("AltSplRP.e.Available", env =
  affylmGUIenvironment): partial argument match of 'env' to 'envir'
.consistentFilter: warning in get("AltSplRP.e.p", env =
  affylmGUIenvironment): partial argument match of 'env' to 'envir'
.consistentFilter: warning in get("midas.p", env =
  affylmGUIenvironment): partial argument match of 'env' to 'envir'
.consistentFilter: warning in get("spliceIndexData", env =
... 2499 lines ...
  data("chrLength", package = "oneChannelGUI")
  data("chrLength", package = "oneChannelGUI")
  data("chrLength", package = "oneChannelGUI")
  data("chrLength", package = "oneChannelGUI")
File ‘oneChannelGUI/R/generaltoolsmenu.R’:
  data(HuExExonProbesetLocation)
  data(MoExExonProbesetLocation)
  data(RaExExonProbesetLocation)
File ‘oneChannelGUI/R/standalonefunctions.R’:
  data("chrLength", package = "oneChannelGUI")
See section ‘Good practice’ in ‘?data’.

checking Rd line widths ... NOTE
Rd file 'oneChannelGUI.Rd':
  \usage lines wider than 90 characters:
     wrapperMirnaCounts(working.dir, out.dir, org = "hsa", threads = 1, cutadapt.path = "/usr/local/bin/cutadapt", parallel = FALSE, ...)

Rd file 'standAloneBuildingLocalAnnotation.Rd':
  \usage lines wider than 90 characters:
        standAloneBuildingLocalAnnotation(libDirLocation = getwd(), netaffxUser = "myemail@somewhere.org", netaffxUserPw = "yourpassword", w ... [TRUNCATED]

These lines will be truncated in the PDF manual.
```

