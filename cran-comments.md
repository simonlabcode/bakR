## R CMD check results

There were no ERRORs or WARNINGs.

There were 3 NOTEs:

* checking to-level files ... NOTE
    Installed size 6.9Mb. 
    Sub-directories of 1 Mb or more: 
        libs 6.1Mb

  bakR uses Rcpp to compile Stan (a probabilistic progamming language) models and the compiled
  model objects are stored in Libs.

* Checking for GNU extensions in Makefiles ... NOTE
    GNU make is a SystemRequirements

  GNU make is needed for handling the Makevars files that rstantools uses to
  compile the Stan models in this package
	
* checking dependencies in R code ... NOTE
  Namespaces in Imports file not imported from: 'RcppParallel' 'rstantools'

  RccpParallel and rstantools are build-time dependencies

## Downstream dependencies
There are currently no downstream dependencies for this package
