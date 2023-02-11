# bakR 0.4.0
* Added QC_check(), a function to perform quality control analysis with bakRFit objects. Looks for any problems in your data that will impair bakR's performance, generates a number of diagnostic visualizations, and makes suggestions about what to do next.
* Fixed plot coloring bug in plotMA() and plotVolcano().
* Fixed bug that led to problems when the number of -s4U replicates > +s4U replicates in one or more Exp_IDs
* Implemented improved U-content adjustment for MCMC implementation. Also impacts accuracy of StanRateEst = TRUE mutation rate estimation strategy.
* Increased default number of features to use for StanRateEst mutation rate estimation strategy.
* Improved scaling of NSSHeat() output matrix columns.
* Created a new function NSSHeat2() that implements a different mechanism scoring function than NSSHeat().

# bakR 0.3.0

* Optimized data preprocessing with data.table
* Increased simulation flexibility
* Added edge case error catching

# bakR 0.2.5

* Added checks to bakRData validator
* Fixed bug in NSSHeatmap that prevented adjusting padj cutoff

# bakR 0.2.4

* Addressed NOTEs to prepare for CRAN submission
* Removed previously deprecated function sim_bakRData()

# bakR 0.2.3

* Expanded discussion of NSS analysis in vignette

# bakR 0.2.2

* Fixed non-steady-state (NSS) uncertainty quantification
* Fixed uncertainty quantification for pulse-chase analysis
* Fixed fraction new estimation for pulse-chase analysis

# bakR 0.2.1

* Added pulse-chase analysis option

# bakR 0.2.0

* Improved default mutation rate estimation strategy
* Updated vignettes and added back non-steady state analysis strategy discussion
* Better commenting of fast_analysis (instead of opting for refactorization)

# bakR 0.1.1

* Changed ordering of vignettes on website
* Got rid of unnecessary stan models and functions
* Fixed TL_stan output documentation
* Corrected missing namespace issues

# bakR 0.1.0

* Made bakR public
