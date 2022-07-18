# Brief Description of bakR
bakR (Bayesian analysis of the kinetics of RNA) is an R package for performing differential kinetic analysis with nucleotide recoding high-throughput RNA sequencing (NR-seq) data. 
Kinetic parameter estimation and statistical testing is compatible with mutational data from any enrichment free NR-seq method (e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.).

# Why use bakR?
Differential expression analysis of RNA sequencing (RNA-seq) data can identify changes in cellular RNA levels, but cannot determine the kinetic mechanism underlying such changes. Previously, our lab and others addressed this shortcoming by developing nucleotide-recoding RNA-seq methods (NR-seq; e.g., TimeLapse-seq) to quantify changes in RNA synthesis and degradation kinetics. While advanced statistical models implemented in user-friendly software (e.g., DESeq2) have ensured the statistical rigor of differential expression analyses, no such tools that facilitate differential kinetic analysis with NR-seq exist. To address this need, we developed bakR, an R package that analyzes and compares NR-seq datasets. Differential kinetic analysis with bakR relies on a new statistical model of NR-seq data that shares data across transcripts in a stastistically principled manner using hierarchical modeling. Look out for our paper describing the model and its extensive validation soon!

# Installation
To install bakR from Github, you need to have a C++ compiler configured to rstan's (the R interface to the probabilistic programming language that bakR uses on the backend) liking. The best way to do this is to follow the Stan team's [helpful documentation](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) on installing rstan for your operating system. Once that is complete, you can install bakR as follows:

    install.packages("devtools") # if you haven't installed devtools already
    devtools::install_github("simonlabcode/bakR")

# Documentation
There are currently two published vignettes to help get you up to speed with using bakR:

  1. An introductory vignette (linked [here](https://rpubs.com/isaacvock/923586)) that walks you through the basic bakR workflow with simulated data.
  2. A second vignette (linked [here](https://rpubs.com/isaacvock/923576)) that dives into some specific follow-up analyses one can perform after running bakR on a dataset. Currently, the only analysis discussed in this vignette is differential synthesis analysis (i.e., identifiying differences in RNA synthesis rates by combining bakR with differential expression analysis), but more analyses will be added to this vignette as time goes by.
  
# Bug Catching and Further Questions
Post descriptions of bugs and a simple reproducible example (if possible) in the Issues section of this repo. In fact, you should go to the Issues section with any question you have about bakR, and there are even helpful labels that you can append to your posts to make the nature of your request clear. If you email me (Isaac Vock) with a question/concern/suggestion, I will direct you to the Issues section. If you have basic use questions, I would suggest going through the vignettes linked above. If these do not answer your question, then post your question to Issues.

