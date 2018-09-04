# aroma.light: Light-Weight Methods for Normalization and Visualization of Microarray Data using Only Basic R Data Types


## Installation
R package aroma.light is available on [Bioconductor](https://www.bioconductor.org/packages/devel/bioc/html/aroma.light.html) and can be installed in R as:

```r
install.packages("BiocManager")
BiocManager::install("aroma.light")
```

### Pre-release version

To install the pre-release version that is available in Git branch `develop` on GitHub, use:
```r
remotes::install_github("HenrikBengtsson/aroma.light@develop")
```
This will install the package from source.  



## Contributions

This Git repository uses the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/aroma.light/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/aroma.light) branch contains the code of the latest release, which is exactly what is currently on [Bioconductor](https://www.bioconductor.org/packages/devel/bioc/html/aroma.light.html).

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [aroma.light repository](https://github.com/HenrikBengtsson/aroma.light).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/aroma.light">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/aroma-light">AppVeyor CI</a> when the PR is submitted.


## Software status

| Resource:     | Bioconductor        | Travis CI       | Appveyor         |
| ------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   | <a href="https://bioconductor.org/checkResults/release/bioc-LATEST/aroma.light/"><img border="0" src="https://bioconductor.org/shields/build/release/bioc/aroma.light.svg" alt="Build status"></a> (release)</br><a href="https://bioconductor.org/checkResults/devel/bioc-LATEST/aroma.light/"><img border="0" src="https://bioconductor.org/shields/build/devel/bioc/aroma.light.svg" alt="Build status"></a> (devel) | <a href="https://travis-ci.org/HenrikBengtsson/aroma.light"><img src="https://travis-ci.org/HenrikBengtsson/aroma.light.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/aroma-light"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/aroma.light?svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/gh/HenrikBengtsson/aroma.light"><img src="https://codecov.io/gh/HenrikBengtsson/aroma.light/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |
