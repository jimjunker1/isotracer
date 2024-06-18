## Reason for submitting an updated version of the package (1.1.5 -> 1.1.6)

- The main reason for this new release is to fix several bugs affecting
  calculations of compartment sizes and flows when considering steady states in
  a network model.
- Several minor improvements were also made (see summary below).


## Summary of update

Main updates:

- Fixed bug with how split compartments were handled when calculating steady
  state sizes. The issue was in the (internal) function
  `calculate_steady_state_one_row()` and affected the (exported) functions
  `tidy_steady_states()`, `tidy_flows()`, and `calculate_steady_state()`.
- Fixed bug in `tidy_flows()` which was present when using `steady_state =
  TRUE` and pulses were present in the network model. Now running
  `tidy_flows(..., steady_state = TRUE)` automatically ignores any pulse that
  might be present in the network model.
- Package version was bumped to 1.1.6.

Minor updates:

- Removed dependency on the tidyverse package in the "Suggests" field of
  `DESCRIPTION`.
- Added tests for `tidy_steady_states()` and `tidy_flows()`.
- Added a function `prop2delta()`.
- Improved argument handling for `delta2prop()`, which now provides helpful
  error messages.


## Test environments (on 2024-05-14)

- Debian 6.1.90-1, R 4.2.2 Patched (2022-11-10 r83330) -- run locally
- Windows Server 2022 x64 (build 20348), R 4.4.0 (2024-04-24 ucrt) -- run on
  https://win-builder.r-project.org/
- Windows Server 2022 x64 (build 20348), R Under development (unstable)
  (2024-05-12 r86534 ucrt) -- run on https://win-builder.r-project.org/
- macOS Ventura 13.3.1, R version 4.4.0 (2024-04-24) -- run on
  https://mac.r-project.org/


## Summary of R CMD check results

There were no ERRORs or WARNINGs.

There were 4 NOTES:


* checking CRAN incoming feasibility ... NOTE
Found the following (possibly) invalid URLs:
  URL: https://www.journals.uchicago.edu/doi/10.1086/708546
    From: inst/doc/tutorial-010-quick-start.html
    Status: 403
    Message: Forbidden

This URL is valid and points to the article describing the original method, on
the journal website.


*  checking installed package size ... NOTE
  installed size is  8.6Mb
  sub-directories of 1Mb or more:
    data   2.0Mb
    doc    2.3Mb
    libs   3.2Mb 

The size of the locally built tarball is 4031112 bytes, respecting the 5 MB
limit for CRAN packages.


*  checking dependencies in R code ... NOTE
Namespace in Imports field not imported from: ‘rstantools’
  All declared Imports should be used.

'rstantools' is present in the Imports field so that the package
installation/configuration is delegated to rstantools for compatibility with
future releases of rstan/StanHeaders. While rstantools is not used by our
package code in the './R/' folder, it is used in the files './configure' and
'./configure.win'.


*  checking for GNU extensions in Makefiles ... NOTE
GNU make is a SystemRequirements. 

The GNU make requirement is introduced by the dependency on the rstan package.


## Downstream dependencies

There are currently no downstream dependencies for this package.
