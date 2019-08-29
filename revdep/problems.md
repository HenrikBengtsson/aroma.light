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
    Error in setClass("ProgressBarText", representation(steps = "integer",  : 
      could not find function "setClass"
    Error in setClass("ProgressBarText", representation(steps = "integer",  : 
      could not find function "setClass"
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
* Number of recursive dependencies: 45

Run `revdep_details(,"aroma.core")` for more info

</details>

## In both

*   checking package dependencies ... NOTE
    ```
    Packages suggested but not available for checking:
      'sfit', 'expectile', 'HaarSeg', 'mpcbs'
    ```

*   checking dependencies in R code ... NOTE
    ```
    Error in setGeneric("image", function(x, ...) standardGeneric("image")) : 
      could not find function "setGeneric"
    ```

# scone

<details>

* Version: 1.8.0
* Source code: https://github.com/cran/scone
* BugReports: https://github.com/YosefLab/scone/issues
* Date/Publication: 2019-05-02
* Number of recursive dependencies: 197

Run `revdep_details(,"scone")` for more info

</details>

## In both

*   checking R code for possible problems ... NOTE
    ```
    ...
    sconeReport : server: no visible global function definition for ‘theme’
    sconeReport : server: no visible global function definition for
      ‘element_blank’
    sconeReport : server: no visible global function definition for
      ‘ggplotly’
    sconeReport : server: no visible global function definition for
      ‘geom_violin’
    sconeReport : server: no visible global function definition for
      ‘coord_cartesian’
    sconeReport : server: no visible global function definition for
      ‘scale_fill_manual’
    sconeReport : server: no visible global function definition for
      ‘geom_point’
    sconeReport : server: no visible global function definition for
      ‘guides’
    Undefined global functions or variables:
      %>% aes coord_cartesian element_blank geom_bar geom_point geom_violin
      ggplot ggplotly guides labs plot_ly plotlyOutput renderVisNetwork
      scale_fill_manual theme visEdges visGroups visHierarchicalLayout
      visLegend visNetwork visNetworkOutput visNetworkProxy visOptions
      visSelectNodes ylim
    ```

# scran

<details>

* Version: 1.12.1
* Source code: https://github.com/cran/scran
* Date/Publication: 2019-05-27
* Number of recursive dependencies: 162

Run `revdep_details(,"scran")` for more info

</details>

## In both

*   checking whether the namespace can be loaded with stated dependencies ... WARNING
    ```
    Error in setClass("SCESet", contains = "ExpressionSet", slots = c(logExprsOffset = "numeric",  : 
      could not find function "setClass"
    Error: unable to load R code in package ‘scater’
    Execution halted
    
    A namespace must be able to be loaded with just the base namespace
    loaded: otherwise if the namespace gets loaded by a saved object, the
    session will be unable to start.
    
    Probably some imports need to be declared in the NAMESPACE file.
    ```

*   checking whether the namespace can be unloaded cleanly ... WARNING
    ```
    Error in classVersion("ExpressionSet") : 
      could not find function "classVersion"
    Error: unable to load R code in package ‘scater’
    Execution halted
    ```

*   checking installed package size ... NOTE
    ```
      installed size is 26.4Mb
      sub-directories of 1Mb or more:
        doc    1.5Mb
        libs  23.9Mb
    ```

# TIN

<details>

* Version: 1.16.0
* Source code: https://github.com/cran/TIN
* Date/Publication: 2019-05-02
* Number of recursive dependencies: 115

Run `revdep_details(,"TIN")` for more info

</details>

## In both

*   checking whether the package can be unloaded cleanly ... WARNING
    ```
    Error in globalVariables("naresid.omit") : 
      could not find function "globalVariables"
    Error: package or namespace load failed for ‘TIN’:
     unable to load R code in package ‘Hmisc’
    Execution halted
    ```

*   checking whether the namespace can be loaded with stated dependencies ... WARNING
    ```
    Error in globalVariables("naresid.omit") : 
      could not find function "globalVariables"
    Error: unable to load R code in package ‘Hmisc’
    Execution halted
    
    A namespace must be able to be loaded with just the base namespace
    loaded: otherwise if the namespace gets loaded by a saved object, the
    session will be unable to start.
    
    Probably some imports need to be declared in the NAMESPACE file.
    ```

*   checking dependencies in R code ... NOTE
    ```
    ...
    Call sequence:
    6: stop(msg, call. = FALSE, domain = NA)
    5: value[[3L]](cond)
    4: tryCatchOne(expr, names, parentenv, handlers[[1L]])
    3: tryCatchList(expr, classes, parentenv, handlers)
    2: tryCatch({
           attr(package, "LibPath") <- which.lib.loc
           ns <- loadNamespace(package, lib.loc)
           env <- attachNamespace(ns, pos = pos, deps, exclude, include.only)
       }, error = function(e) {
           P <- if (!is.null(cc <- conditionCall(e))) 
               paste(" in", deparse(cc)[1L])
           else ""
           msg <- gettextf("package or namespace load failed for %s%s:\n %s", 
               sQuote(package), P, conditionMessage(e))
           if (logical.return) 
               message(paste("Error:", msg), domain = NA)
           else stop(msg, call. = FALSE, domain = NA)
       })
    1: library(package, lib.loc = lib.loc, character.only = TRUE, verbose = FALSE)
    Execution halted
    ```

*   checking R code for possible problems ... NOTE
    ```
    Error in globalVariables("naresid.omit") : 
      could not find function "globalVariables"
    Error: unable to load R code in package ‘Hmisc’
    Execution halted
    ```

