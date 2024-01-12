# bakR 1.0.1

# bakR 1.0.0
* Functions for visualizing (`VisualizeDropout`), quantifying (`QuantifyDropout`), and correcting (`CorrectDropout`) metabolic label-induced dropout of RNA during library preparation have been added. 
* New simulation function (`simulate_relative_bakRData`) which better captures the relative nature of RNA-seq and can accurately simulate dropout.
* New experimental function (`DissectMechanism`) for determining how likely that any observed differential expression is driven by transcriptional or post-transcriptional regulation. `DissectMechanism` is a rewrite and extension of the previously developed `NSSHeat2` function, which itself was an improvement of the now deprecated `NSSHeat`.
* Can now provide fraction new estimates (e.g., from a tool like GRAND-SLAM) as input to bakR. GRAND-SLAM input functionality is further supported by the new `GSprocess` function that will facilitate converting from GRAND-SLAM output to bakR input.
* `FnPCA` has been deprecated in favor of `FnPCA2` which accepts input differently and fixes some bugs.
* Read count filtering now includes two filters. One read count that all samples must pass, and one that only all replicates in a single Exp_ID need to pass. This facilitates identifying large increases or decreaes in expression.
* Several new vignettes to discuss much of the new functionality discussed above.
* Several small bug fixes

# bakR 0.4.4
* Small edit to configuration files that address compilation issues that can arise on some systems. Deals with "file too big" errors" during package installation from source.

# bakR 0.4.3
* Implemented long read sequencing data analysis strategy. Run bakRFit() with Long = TRUE to use k-means clustering (k = 2) for mutation rate and fraction new estimation. Need to have Ckmeans.1d.dp package installed to do this (not installed during bakR installation).

# bakR 0.4.2
* Fixed bug in reliableFeatures. high_p was supposed to be the maximum allowable mutation rate (# of mutations/# of Ts) in reads from -s4U controls, but was instead the maximum allowable average number of T-to-C mutations in reads from these controls.

# bakR 0.4.1
* Fixed bug in cBprocess that didn't properly check that features of interest provided by the FOI argument were valid.

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
