# Bayesian analysis of the kinetics of RNA (bakR)
R package for analyzing nucleotide recoding high-throughput sequencing data. 
Kinetic parameter estimation and statistical testing is compatible with any enrichment free metabolic labeling mutational data.
# Installation
There are two important dependencies that need to be installed prior to installing this package: Stan and rstan. Go to ``http://mc-stan.org`` and follow the installation instructions for your platform. The biggest challenge is getting a C++ compiler configured to work with your installation of R. The instructions at ``https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started`` are hopefully helpful on that front. Once that is complete, you can install bakR. Currently, the only way to install the package is to:
  1. Clone the repository
  2. Navigate to the directory containing the cloned repository (not the repo directory itself)
  3. Call (from R or R_studio) devtools::install("bakR") 
  4. Can load package with library(bakR) like normal

Small updates will constantly be made, so make sure to pull changes from Github and reinstall frequently.

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
    
    
  5. Install bakR:
    
    

    devtools::install("bakR")

    
    
  6. Exit interactive session; now you can include library(bakR) in any of your batch jobs to use this package.

# Documentation
There are currently two published vignettes to help get you up to speed with bakR:

  1. An introductory vignette (link [here](https://rpubs.com/isaacvock/923586)) that walks you through the basic bakR workflow with simulated data.
  2. A second vignette (link [here](https://rpubs.com/isaacvock/923576)) that dives into some specific follow-up analyses one can do after running bakR on a dataset. Currently, the only analysis discussed in this vignette is differential synthesis analysis (i.e., identifiying differences in RNA synthesis rates by combining bakR with differential expression analysis), but more analyses will be added to this vignette as time goes by.
# Bug Catching
Post descriptions of bugs and a simple reproducible example (if possible) in the Issues section of this repo. In fact, you should go to the Issues section with any question you have about bakR, and there are even helpful labels that you can append to your posts to make the nature of your request clear. If you email me (Isaac Vock) with a question/concern/suggestion, I will direct you to the Issues section. If you have basic use questions, I would suggest going through the vignettes linked above/ If these do not answer your question, then post your question to Issues.

