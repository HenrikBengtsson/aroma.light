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

25 packages

|package           |version | errors| warnings| notes|
|:-----------------|:-------|------:|--------:|-----:|
|ACNE              |0.8.1   |      0|        0|     0|
|aroma.affymetrix  |3.1.0   |      0|        0|     0|
|aroma.cn          |1.6.1   |      0|        0|     0|
|aroma.core        |3.1.0   |      0|        0|     0|
|calmate           |0.12.1  |      0|        0|     0|
|EDASeq            |2.8.0   |      0|        0|     2|
|EnrichmentBrowser |2.4.6   |      0|        0|     3|
|HTSCluster        |2.0.8   |      0|        0|     0|
|HTSFilter         |1.14.1  |      0|        0|     0|
|metaRNASeq        |1.0.2   |      0|        0|     1|
|metaseqR          |1.14.0  |      1|        1|     4|
|MoonlightR        |1.0.0   |      0|        0|     2|
|MPAgenomics       |1.1.2   |      0|        0|     2|
|NSA               |0.0.32  |      0|        0|     6|
|oneChannelGUI     |1.40.0  |      1|        0|     0|
|PBNPA             |0.0.1   |      0|        0|     0|
|PECA              |1.10.0  |      0|        0|     1|
|PSCBS             |0.62.0  |      0|        0|     0|
|PureCN            |1.2.3   |      0|        0|     1|
|Repitools         |1.20.0  |      0|        0|     2|
|RUVSeq            |1.8.0   |      0|        0|     2|
|SpidermiR         |1.4.8   |      0|        0|     1|
|StarBioTrek       |1.0.3   |      0|        0|     2|
|TCGAbiolinks      |2.2.10  |      0|        0|     3|
|TIN               |1.6.0   |      0|        0|     2|

## ACNE (0.8.1)
Maintainer: Henrik Bengtsson <henrikb@braju.com>  
Bug reports: https://github.com/HenrikBengtsson/ACNE/issues

0 errors | 0 warnings | 0 notes

## aroma.affymetrix (3.1.0)
Maintainer: Henrik Bengtsson <henrikb@braju.com>  
Bug reports: https://github.com/HenrikBengtsson/aroma.affymetrix/issues

0 errors | 0 warnings | 0 notes

## aroma.cn (1.6.1)
Maintainer: Henrik Bengtsson <henrikb@braju.com>  
Bug reports: https://github.com/HenrikBengtsson/aroma.cn/issues

0 errors | 0 warnings | 0 notes

## aroma.core (3.1.0)
Maintainer: Henrik Bengtsson <henrikb@braju.com>  
Bug reports: https://github.com/HenrikBengtsson/aroma.core/issues

0 errors | 0 warnings | 0 notes

## calmate (0.12.1)
Maintainer: Henrik Bengtsson <henrikb@braju.com>  
Bug reports: https://github.com/HenrikBengtsson/calmate/issues

0 errors | 0 warnings | 0 notes

## EDASeq (2.8.0)
Maintainer: Davide Risso <risso.davide@gmail.com>  
Bug reports: https://github.com/drisso/EDASeq/issues

0 errors | 0 warnings | 2 notes

```
checking R code for possible problems ... NOTE
.availableOrgPkgs: no visible global function definition for
  ‘available.packages’
.gcLoess : ff: no visible global function definition for ‘quantile’
.gcLoess : ff: no visible global function definition for ‘loess’
.gcLoess : ff: no visible global function definition for ‘predict’
.gcLoess : ff: no visible global function definition for ‘median’
.gcQuant: no visible global function definition for ‘quantile’
.gcQuant : f : <anonymous>: no visible global function definition for
  ‘median’
... 114 lines ...
  validObject var
Consider adding
  importFrom("grDevices", "rainbow")
  importFrom("graphics", "abline", "legend", "lines", "matplot", "pairs",
             "points", "smoothScatter", "text")
  importFrom("methods", "as", "callNextMethod", "new", "validObject")
  importFrom("stats", "loess", "lowess", "median", "predict", "quantile",
             "var")
  importFrom("utils", "available.packages")
to your NAMESPACE file (and ensure that your DESCRIPTION Imports field
contains 'methods').

checking Rd line widths ... NOTE
Rd file 'newSeqExpressionSet.Rd':
  \usage lines wider than 90 characters:
                         normalizedCounts = matrix(data=NA, nrow=nrow(counts), ncol=ncol(counts), dimnames=dimnames(counts)),
                         offset = matrix(data=0, nrow=nrow(counts), ncol=ncol(counts), dimnames=dimnames(counts)),

Rd file 'withinLaneNormalization-methods.Rd':
  \usage lines wider than 90 characters:
     withinLaneNormalization(x, y, which=c("loess","median","upper","full"), offset=FALSE, num.bins=10, round=TRUE)

These lines will be truncated in the PDF manual.
```

## EnrichmentBrowser (2.4.6)
Maintainer: Ludwig Geistlinger <Ludwig.Geistlinger@bio.ifi.lmu.de>

0 errors | 0 warnings | 3 notes

```
checking dependencies in R code ... NOTE
Unexported object imported by a ':::' call: 'pathview:::parseKGML2Graph2'
  See the note in ?`:::` about the use of this operator.

checking R code for possible problems ... NOTE
Found the following assignments to the global environment:
File ‘EnrichmentBrowser/R/sbea.R’:
  assign("eset", eset, envir = .GlobalEnv)
  assign("local.de.ana", local.de.ana, envir = .GlobalEnv)
  assign(global.func, get(global.func), envir = .GlobalEnv)

Found the following calls to data() loading into the global environment:
File ‘EnrichmentBrowser/R/pathview2.R’:
  data("bods", package = "pathview")
  data("gene.idtype.list", package = "pathview")
  data("KEGGEdgeSubtype", package = "pathview")
File ‘EnrichmentBrowser/R/probeEset2geneEset.R’:
  data(korg, package = "pathview")
See section ‘Good practice’ in ‘?data’.

checking Rd line widths ... NOTE
Rd file 'de.ana.Rd':
  \usage lines wider than 90 characters:
             de.method = c("limma", "edgeR", "DESeq"), padj.method = "BH", stat.only=FALSE, min.cpm=2 )

These lines will be truncated in the PDF manual.
```

## HTSCluster (2.0.8)
Maintainer: Andrea Rau <andrea.rau@jouy.inra.fr>

0 errors | 0 warnings | 0 notes

## HTSFilter (1.14.1)
Maintainer: Andrea Rau <andrea.rau@jouy.inra.fr>

0 errors | 0 warnings | 0 notes

## metaRNASeq (1.0.2)
Maintainer: Guillemette Marot <guillemette.marot@inria.fr>

0 errors | 0 warnings | 1 note 

```
checking R code for possible problems ... NOTE
fishercomb: no visible global function definition for ‘pchisq’
fishercomb: no visible global function definition for ‘p.adjust’
invnorm : <anonymous>: no visible global function definition for
  ‘qnorm’
invnorm: no visible global function definition for ‘pnorm’
invnorm: no visible global function definition for ‘p.adjust’
sim.function : <anonymous>: no visible global function definition for
  ‘rnorm’
sim.function : <anonymous>: no visible global function definition for
  ‘rnbinom’
Undefined global functions or variables:
  p.adjust pchisq pnorm qnorm rnbinom rnorm
Consider adding
  importFrom("stats", "p.adjust", "pchisq", "pnorm", "qnorm", "rnbinom",
             "rnorm")
to your NAMESPACE file.
```

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

## MoonlightR (1.0.0)
Maintainer: Antonio Colaprico <antonio.colaprico@ulb.ac.be>, Catharina
 Olsen <colsen@ulb.ac.be>  
Bug reports: https://github.com/torongs82/Moonlight/issues

0 errors | 0 warnings | 2 notes

```
checking installed package size ... NOTE
  installed size is  5.9Mb
  sub-directories of 1Mb or more:
    data   3.1Mb
    doc    2.6Mb

checking Rd line widths ... NOTE
Rd file 'plotURA.Rd':
  \examples lines wider than 100 characters:
     plotURA(dataURA = dataURA[c(names(dataDual$TSG), names(dataDual$OCG)),],additionalFilename = "_example")

These lines will be truncated in the PDF manual.
```

## MPAgenomics (1.1.2)
Maintainer: Samuel Blanck <samuel.blanck@inria.fr>

0 errors | 0 warnings | 2 notes

```
checking dependencies in R code ... NOTE
'library' or 'require' calls in package code:
  ‘R.devices’ ‘R.filesets’ ‘R.methodsS3’ ‘R.oo’ ‘aroma.affymetrix’
  ‘aroma.cn’ ‘aroma.core’ ‘aroma.light’ ‘matrixStats’ ‘snowfall’
  Please use :: or requireNamespace() instead.
  See section 'Suggested packages' in the 'Writing R Extensions' manual.
Unexported object imported by a ':::' call: ‘cghseg:::segmeanCO’
  See the note in ?`:::` about the use of this operator.

checking R code for possible problems ... NOTE
.varregtimescount: no visible global function definition for ‘var’
CGHSEGaroma: no visible global function definition for ‘read.csv’
CGHSEGaroma : <anonymous>: no visible global function definition for
  ‘points’
CGHSEGaroma : <anonymous>: no visible global function definition for
  ‘lines’
CGHSEGaroma : <anonymous>: no visible global function definition for
  ‘write.table’
CGHcall: no visible global function definition for ‘mad’
... 43 lines ...
tumorboostPlot: no visible global function definition for ‘par’
tumorboostPlot: no visible global function definition for ‘axis’
tumorboostPlot: no visible global function definition for ‘points’
Undefined global functions or variables:
  axis head lines lm mad median optim par points read.csv sd var
  write.table
Consider adding
  importFrom("graphics", "axis", "lines", "par", "points")
  importFrom("stats", "lm", "mad", "median", "optim", "sd", "var")
  importFrom("utils", "head", "read.csv", "write.table")
to your NAMESPACE file.
```

## NSA (0.0.32)
Maintainer: Maria Ortiz-Estevez <mortizest@gmail.com>

0 errors | 0 warnings | 6 notes

```
checking package dependencies ... NOTE
Depends: includes the non-default packages:
  ‘R.methodsS3’ ‘MASS’ ‘matrixStats’ ‘R.oo’ ‘R.utils’ ‘aroma.core’
  ‘aroma.affymetrix’ ‘DNAcopy’
Adding so many packages to the search path is excessive and importing
selectively is preferable.

checking top-level files ... NOTE
Non-standard file/directory found at top level:
  ‘incl’

checking dependencies in R code ... NOTE
Packages in Depends field not imported from:
  ‘DNAcopy’ ‘MASS’ ‘R.methodsS3’ ‘R.oo’ ‘aroma.affymetrix’ ‘aroma.core’
  ‘matrixStats’
  These packages need to be imported from (in the NAMESPACE file)
  for when this namespace is loaded but not attached.

checking S3 generic/method consistency ... NOTE
Found the following apparent S3 methods exported but not registered:
  NSAByTotalAndFracB.matrix allocateOutputDataSets.NSANormalization
  allocateOutputDataSets.SNPsNormalization
  allocateOutputDataSets.SampleNormalization
  as.character.NSANormalization as.character.SNPsNormalization
  as.character.SampleNormalization findArraysTodo.NSANormalization
  findArraysTodo.SampleNormalization findUnitsTodo.SNPsNormalization
  fitNSA.matrix fitNSAcnPs.matrix getDataSets.NSANormalization
  getDataSets.SNPsNormalization getDataSets.SampleNormalization
  getName.NSANormalization getName.SNPsNormalization
  getName.SampleNormalization getOutputDataSets.NSANormalization
  getOutputDataSets.SNPsNormalization
  getOutputDataSets.SampleNormalization getPath.NSANormalization
  getPath.SNPsNormalization getPath.SampleNormalization
  getRootPath.NSANormalization getRootPath.SNPsNormalization
  getRootPath.SampleNormalization process.NSANormalization
  process.SNPsNormalization process.SampleNormalization
  sampleNByTotalAndFracB.numeric snpsNByTotalAndFracB.matrix
See section ‘Registering S3 methods’ in the ‘Writing R Extensions’
manual.

checking R code for possible problems ... NOTE
NB: .First.lib is obsolete and will not be used in R >= 3.0.0

.First.lib: no visible global function definition for
  ‘packageDescription’
NSAByTotalAndFracB.matrix: no visible global function definition for
  ‘throw’
NSAByTotalAndFracB.matrix: no visible global function definition for
  ‘str’
NSANormalization: no visible global function definition for ‘throw’
... 279 lines ...
  extractMatrix findUnitsTodo getAsteriskTags getChipType getFile
  getFullName getFullNames getGenomeInformation getName getNames
  getPath getPathname getPathnames getPositions getRam getRootPath
  getTags getUnitsOnChromosome hist median nbrOfFiles newInstance
  packageDescription rowAlls rowMedians segment setTags str throw trim
  verbose
Consider adding
  importFrom("graphics", "hist")
  importFrom("stats", "approxfun", "median")
  importFrom("utils", "packageDescription", "str")
to your NAMESPACE file.

checking Rd line widths ... NOTE
Rd file 'NSANormalization.Rd':
  \examples lines wider than 100 characters:
     by <- 50e3; # 50kb bins; you may want to try with other amounts of smoothing xOut <- seq(from=xRange[1], to=xRange[2], by=by);
     plot(getSignals(cnCNPS), getSignals(cnSNPS), xlim=Clim, ylim=Clim); abline(a=0, b=1, col="red", lwd=2);

These lines will be truncated in the PDF manual.
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

## PBNPA (0.0.1)
Maintainer: Gaoxiang Jia <GJia@SMU.edu>

0 errors | 0 warnings | 0 notes

## PECA (1.10.0)
Maintainer: Tomi Suomi <tomi.suomi@utu.fi>

0 errors | 0 warnings | 1 note 

```
checking Rd line widths ... NOTE
Rd file 'PECA.Rd':
  \usage lines wider than 90 characters:
     PECA_AffyBatch(affy=NULL, normalize=FALSE, test="t", type="median", paired=FALSE, progress=FALSE)

These lines will be truncated in the PDF manual.
```

## PSCBS (0.62.0)
Maintainer: Henrik Bengtsson <henrikb@braju.com>  
Bug reports: https://github.com/HenrikBengtsson/PSCBS/issues

0 errors | 0 warnings | 0 notes

## PureCN (1.2.3)
Maintainer: Markus Riester <markus.riester@novartis.com>

0 errors | 0 warnings | 1 note 

```
checking R code for possible problems ... NOTE
Found the following calls to data() loading into the global environment:
File ‘PureCN/R/bootstrapResults.R’:
  data(purecn.example.output)
File ‘PureCN/R/callAlterations.R’:
  data(purecn.example.output)
  data(purecn.example.output)
File ‘PureCN/R/callLOH.R’:
  data(purecn.example.output)
File ‘PureCN/R/createCurationFile.R’:
  data(purecn.example.output)
File ‘PureCN/R/curateResults.R’:
  data(purecn.example.output)
  data(purecn.example.output)
File ‘PureCN/R/plotAbs.R’:
  data(purecn.example.output)
File ‘PureCN/R/predictSomatic.R’:
  data(purecn.example.output)
File ‘PureCN/R/readCurationFile.R’:
  data(purecn.example.output)
See section ‘Good practice’ in ‘?data’.
```

## Repitools (1.20.0)
Maintainer: Mark Robinson <mark.robinson@imls.uzh.ch>

0 errors | 0 warnings | 2 notes

```
checking R code for possible problems ... NOTE
Found an obsolete/platform-specific call in the following function:
  ‘maskOut’
Found the platform-specific device:
  ‘windows’
dev.new() is the preferred way to open a new device, in the unlikely
event one is needed.
.cpgBoxplots: no visible global function definition for ‘pdf’
.cpgBoxplots: no visible global function definition for ‘par’
.cpgBoxplots: no visible global function definition for ‘dev.off’
... 291 lines ...
  rainbow read.table rect str t.test text title verbose
Consider adding
  importFrom("grDevices", "dev.off", "pdf", "rainbow")
  importFrom("graphics", "abline", "axis", "barplot", "bxp", "grid",
             "layout", "legend", "lines", "matlines", "matplot", "mtext",
             "par", "persp", "plot", "plot.new", "plot.window", "points",
             "polygon", "rect", "text", "title")
  importFrom("stats", "dbeta", "embed", "filter", "kmeans", "lm",
             "lowess", "p.adjust", "predict", "pt", "qnorm", "t.test")
  importFrom("utils", "read.table", "str")
to your NAMESPACE file.

checking Rd line widths ... NOTE
Rd file 'ChromaBlocks.Rd':
  \usage lines wider than 90 characters:
     ChromaBlocks(rs.ip, rs.input, organism, chrs, ipWidth=100, inputWidth=500, preset=NULL, blockWidth=NULL, minBlocks=NULL, extend=NULL, c ... [TRUNCATED]

Rd file 'GCbiasPlots.Rd':
  \usage lines wider than 90 characters:
                 cex = 0.2, pch.col = "black", line.col = "red", lty = 1, lwd = 2, verbose = TRUE)

Rd file 'absoluteCN.Rd':
... 57 lines ...

Rd file 'regionStats.Rd':
  \usage lines wider than 90 characters:
     regionStats(x, design = NULL, maxFDR=0.05, n.perm=5, window=600, mean.trim=.1, min.probes=10, max.gap=500, two.sides=TRUE, ndf, return. ... [TRUNCATED]
     regionStats(x, design = NULL, maxFDR=0.05, n.perm=5, window=600, mean.trim=.1, min.probes=10, max.gap=500, two.sides=TRUE, ind=NULL, re ... [TRUNCATED]

Rd file 'writeWig.Rd':
  \usage lines wider than 90 characters:
     writeWig(rs, seq.len = NULL, design=NULL, sample=20, drop.zero=TRUE, normalize=TRUE, verbose=TRUE)

These lines will be truncated in the PDF manual.
```

## RUVSeq (1.8.0)
Maintainer: Davide Risso <risso.davide@gmail.com>  
Bug reports: https://github.com/drisso/RUVSeq/issues

0 errors | 0 warnings | 2 notes

```
checking R code for possible problems ... NOTE
residuals.DGEGLM : <anonymous>: no visible global function definition
  for ‘poisson’
Undefined global functions or variables:
  poisson
Consider adding
  importFrom("stats", "poisson")
to your NAMESPACE file.

checking Rd line widths ... NOTE
Rd file 'RUVr.Rd':
  \usage lines wider than 90 characters:
     RUVr(x, cIdx, k, residuals, center=TRUE, round=TRUE, epsilon=1, tolerance=1e-8, isLog=FALSE)

These lines will be truncated in the PDF manual.
```

## SpidermiR (1.4.8)
Maintainer: Claudia Cava <claudia.cava@ibfm.cnr.it>  
Bug reports: https://github.com/claudiacava/SpidermiR/issues

0 errors | 0 warnings | 1 note 

```
checking R code for possible problems ... NOTE
.SpidermiRvisualize_gene: possible error in simpleNetwork(NetworkData,
  linkColour = "gray", textColour = "black", zoom = TRUE): unused
  argument (textColour = "black")
SpidermiRvisualize_plot_target: no visible binding for global variable
  ‘miRNAs’
SpidermiRvisualize_plot_target: no visible binding for global variable
  ‘mRNA_target’
Undefined global functions or variables:
  mRNA_target miRNAs
```

## StarBioTrek (1.0.3)
Maintainer: Claudia Cava <claudia.cava@ibfm.cnr.it>  
Bug reports: https://github.com/claudiacava/StarBioTrek/issues

0 errors | 0 warnings | 2 notes

```
checking installed package size ... NOTE
  installed size is 10.9Mb
  sub-directories of 1Mb or more:
    data   8.8Mb
    doc    1.9Mb

checking R code for possible problems ... NOTE
GE_matrix: no visible binding for global variable ‘path’
getKEGGdata: no visible binding for global variable ‘Carbohydrate’
getKEGGdata: no visible binding for global variable ‘Energy’
getKEGGdata: no visible binding for global variable ‘Lipid’
getKEGGdata: no visible binding for global variable ‘Aminoacid’
getKEGGdata: no visible binding for global variable ‘Glybio_met’
getKEGGdata: no visible binding for global variable ‘Cofa_vita_met’
getKEGGdata: no visible binding for global variable ‘Transcription’
getKEGGdata: no visible binding for global variable ‘Translation’
... 23 lines ...
matrix_plot: no visible binding for global variable ‘path’
plotting_cross_talk: no visible binding for global variable ‘path’
svm_classification: no visible binding for global variable ‘Target’
Undefined global functions or variables:
  Aminoacid Carbohydrate Cell_growth_and_death Cellular_community
  Circulatory_system Cofa_vita_met Digestive_system Endocrine_system
  Energy Excretory_system Folding_sorting_and_degradation Glybio_met
  Immune_system Lipid Nervous_system Replication_and_repair
  Sensory_system Signal_transduction
  Signaling_molecules_and_interaction Target Transcription Translation
  Transport_and_catabolism path
```

## TCGAbiolinks (2.2.10)
Maintainer: Antonio Colaprico <antonio.colaprico@ulb.ac.be>,
 Tiago Chedraoui Silva <tiagochst@usp.br>  
Bug reports: https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues

0 errors | 0 warnings | 3 notes

```
checking installed package size ... NOTE
  installed size is 65.7Mb
  sub-directories of 1Mb or more:
    R      1.1Mb
    data   6.4Mb
    doc   57.9Mb

checking R code for possible problems ... NOTE
GDCquery_clinic: no visible binding for global variable ‘portions’
TCGAvisualize_oncoprint: no visible binding for global variable ‘value’
calculate.pvalues : <anonymous>: no visible binding for global variable
  ‘aux’
Undefined global functions or variables:
  aux portions value

checking for unstated dependencies in vignettes ... NOTE
'library' or 'require' call not declared from: ‘DT’
```

## TIN (1.6.0)
Maintainer: Bjarne Johannessen <bjajoh@rr-research.no>

0 errors | 0 warnings | 2 notes

```
checking top-level files ... NOTE
Non-standard file/directory found at top level:
  ‘doc’

checking R code for possible problems ... NOTE
aberrantExonUsage: no visible global function definition for ‘quantile’
aberrantExonUsage: no visible global function definition for ‘ave’
clusterPlot: no visible global function definition for ‘dist’
clusterPlot: no visible global function definition for ‘hclust’
clusterPlot: no visible global function definition for
  ‘colorRampPalette’
clusterPlot: no visible global function definition for ‘par’
clusterPlot: no visible global function definition for ‘png’
clusterPlot: no visible global function definition for ‘jpeg’
... 50 lines ...
  importFrom("stats", "ave", "dist", "hclust", "median", "quantile")
  importFrom("utils", "data", "read.table")
to your NAMESPACE file.

Found the following assignments to the global environment:
File ‘TIN/R/aberrantExonUsage.R’:
  assign("quantiles", quantiles, envir = .GlobalEnv)
  assign("aberrantExons", aberrantExons, envir = .GlobalEnv)
File ‘TIN/R/correlationPlot.R’:
  assign("randomGeneSetsDist", B, envir = .GlobalEnv)
  assign("traPermutationsDist", L, envir = .GlobalEnv)
```

