## Comments for first submission (Date: 2022-10-07 03:53:58 UTC)

* Please do not start the description with "This package", package name, title or similar

  Changed DESCRIPTION accordingly.

* Please always write package names, software names and API names in single quotes in title and description

  Changed DESCRIPTION accordingly.

* If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year) <arXiv:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
auto-linking.

  Added citation and link to preprint describing methods.

* The Title field should be in title case

  Changed title so that strict title case did not conflict with common casing standards in the field of RNA sequencing.

* Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation.

  Added details of function return values to all exported functions, and converted unnecessarily exported functions to imported functions.

* \dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user.

  Previous uses of \dontrun{} replaced with runnable examples and wrapped in donttest{} due to their > 5 second execution times.

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
