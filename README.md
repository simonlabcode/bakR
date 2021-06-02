# DynamicSeq
R package for analyzing nucleotide recoding high-throughput sequencing data
Kinetic parameter estimation and statistical testing is compatible with any enrichment free metabolic labeling mutational data.
# Installation
Currently, the only way to install the package is to:
  1. Clone the repository
  2. Navigate to the directory containing the cloned repository (not the repo directory itself)
  3. Call (from R or R_studio) devtools::install("DynamicSeq") then library(DynamicSeq)
# Bug Catching
Contact Isaac (isaac.vock@yale.edu or DM on Slack) with descriptions of bugs encountered or suggested improvements/functions. Thank you!
# Current Functions Implemented
All functions can be found in R folder, which each script including several similar functions. Concsise documentation is included as
comments. The existing functions available for beta testing are:

cB_to_Stan(cB, samp_list, c_list, type_list, mut_list, rep_list, tl, nreps, keep_input)
  * This function constructs list necessary to analyze data with Stan
  * cB is the unprocessed cB.rds file from the pipeline
  * samp_list is a vector with 1 entry per sample; each entry should be the sample name as it appears in the cB
  * c_list is a vector with 1 entry per control samples; each entry should be name of a no-s4U control samples as it appears in the cB
  * type_list is a vector with 1 entry per sample; 0 = no s4U, 1 = s4U fed; ith entry of type_list should describe ith entry of samp_list
  * mut_list is a vector with 1 entry per sample; indexes the experimental condition of each sample in samp_list (e.g., 1 = WT, > 1 = mutants)
  * rep_list is a vector with 1 entry per sample that indexes replicates; 1 = 1st replicate of ith sample in samp_list, 2 = 2nd replicate, etc.
  * tl is a single numerical value; label time for s4U feeds
  * nreps is a single numerical value; nubmer of replicates (assumes same number of replicates for all s4U fed samples)
  * keep_input is a two element vector; 1st element is highest mut rate accepted in control samples, 2nd element is read cutoff
  * returns list that can be passed to Stan
  * Note, this is only compatible with the new replicate Stan model (TL_poisson_replicate_model.stan)

cB_to_fast(same input as cB_to_Stan)
  * This function constructs dataframe necessary to analyze data with efficient analysis
  * Input is identical as that for cB_to_Stan but output is a dataframe
  * Returns dataframe that can be passed to fast_analysis()

fast_analysis(df, boot_iter)
  * This function estimates fraction news in each replicate along with the uncertainty in the estimate without Stan
  * df is a dataframe with the same form as the output of cB_to_fast
  * boot_iter Number of times to resample for bootstrapping; default is 50
  * Output is dataframe with fraction new and logit(fraction new) estimate for each replicate as well as the uncertainty for the logit(fraction new)
  
reliableFeatures(cB, c_list, high_p, totcut)
  * This function identifies features with low control sample mutation rates and high read count totals
  * cB is the cB.rds file from the pipeline with at least sample, XF, TC, and n columns
  * c_list is a vector with 1 entry per control sample; each entry should be the name of a no-s4U control sample as it appears in the cB
  * high_p is the mutation rate cutoff for the control sample
  * totcut is the total read count cutoff
