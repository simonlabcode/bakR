## R CMD check results

There were no ERRORs or WARNINGs.

There were 5 NOTEs:

* This is a new release.

* checking to-level files ... NOTE
    Installed size 10.5Mb. 
    Sub-directories of 1 Mb or more: 
        libs 9.7Mb

  bakR uses Rcpp to compile Stan (a probabilistic progamming language) models and the compiled
  model objects are stored in Libs.

* checking compilation flags used ... NOTE
    Compilation used the following non-portable flag(s):
        '-mmmx', '-msse', '-msse3', '-msse4.1', '-msse4.2', '-mssse3'

  These flags only apply to my particular system setup and thus do not need to be portable.
  This NOTE does not show up when running checks on other platforms with the rhub package.

* Checking for GNU extensions in Makefiles ... NOTE
    GNU make is a SystemRequirements

  GNU make is needed for handling the Makevars files that rstantools uses to
  compile the Stan models in this package
	
* checking dependencies in R code ... NOTE
  Namespaces in Imports file not imported from: 'RcppParallel' 'rstantools'

  RccpParallel and rstantools are build-time dependencies

## Downstream dependencies
There are currently no downstream dependencies for this package
