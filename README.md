# DynamicSeq
R package for analyzing nucleotide recoding high-throughput sequencing data
Kinetic parameter estimation and statistical testing is compatible with any enrichment free metabolic labeling mutational data.
# Installation
Currently, the only way to install the package is to:
  1. Clone the repository
  2. Navigate to the directory containing the cloned repository (not the repo directory itself)
  3. Call (from R or R_studio) devtools::install("DynamicSeq") then library(DynamicSeq)
  
If trying to run on the cluster, installation should be done through an interactive R session. The steps are as follows:
  1. Clone the repository to your ruddle account (or use Globus to transfer the locally cloned copy) and navigate to directory containing repo directory (not the repo directory      itself)
  2. Make a new directory to house installed packages (this prevents errors during installation due to permissions):
    
    
    mkdir ~/R/x86_64-pc-linux-gnu-library/4.0
    
    
  3. Start an interactive job and navigate to R:
    
    
    srun --pty -C oldest -p interactive bash
    module load R
    R
    
    
  4. Load devtools:
    
    
    library(devtools)
    
    
  5. Install DynamicSeq:
    
    
    devtools::install("DynamicSeq")
    
    
  6. Exit interactive session; now you can include library(DynamicSeq) in any of your batch jobs to use this package.
# Bug Catching
Contact Isaac (isaac.vock@yale.edu or DM on Slack) with descriptions of bugs encountered or suggested improvements/functions. Thank you!

