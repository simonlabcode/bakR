# Brief Description of bakR
bakR (Bayesian analysis of the kinetics of RNA) is an R package for performing differential kinetics analysis with nucleotide recoding high-throughput RNA sequencing (NR-seq) data. 
Kinetic parameter estimation and statistical testing is compatible with mutational data from any enrichment free NR-seq method (e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.).

# Why use bakR?
Differential expression analysis of RNA sequencing (RNA-seq) data can identify changes in cellular RNA levels, but cannot determine the kinetic mechanism underlying such changes. Previously, [our lab](https://simonlab.yale.edu/research/transcriptome-dynamics/timelapse-chemistry/) and others addressed this shortcoming by developing nucleotide-recoding RNA-seq methods (NR-seq; e.g., TimeLapse-seq) to quantify changes in RNA synthesis and degradation kinetics. While advanced statistical models implemented in user-friendly software (e.g., DESeq2) have ensured the statistical rigor of differential expression analyses, no such tools that facilitate differential kinetic analysis with NR-seq exist. To address this need, we developed bakR, an R package that analyzes and compares NR-seq datasets. Differential kinetics analysis with bakR relies on a Bayesian hierarchical model of NR-seq data to increase statistical power by sharing information across transcripts. bakR outperforms attempts to use single sample analysis tools (e.g., pulseR and GRAND-SLAM) for differential kinetics analyis. Check out [our preprint](https://www.biorxiv.org/content/10.1101/2022.09.02.505697v1) to learn more about the model and its extensive validation!

# Installation
To install bakR from Github, you need to have a C++ compiler configured to rstan's (the R interface to the probabilistic programming language [Stan](https://mc-stan.org/) that bakR uses on the backend) liking. The best way to do this is to follow the Stan team's [helpful documentation](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) on installing rstan for your operating system. Once that is complete, you can install bakR as follows:

    install.packages("devtools") # if you haven't installed devtools already
    devtools::install_github("simonlabcode/bakR")

# Documentation
There are currently two vignettes to help get you up to speed with using bakR:

  1. An introductory vignette (title: Differential Kinetic Analysis with bakR) that walks you through the basic bakR workflow with simulated data.
  2. A second vignette (title: Going beyond the defaults with bakR) that dives into some specific follow-up analyses one can perform after running bakR on a dataset. Currently, this vignette discusses differential synthesis analysis (i.e., identifiying differences in RNA synthesis rates by combining bakR with differential expression analysis) and performing analyses without assuming steady-state. More analyses will be added to this vignette as time goes by.
  
All vignettes are available on the [bakR website](https://simonlabcode.github.io/bakR/index.html) under the Articles section. [Here](https://github.com/simonlabcode/bakR) is the link to the bakR github as well if you need help getting back to the github from the website.

# Bug Catching and Further Questions
Post descriptions of bugs and a simple reproducible example (if possible) in the Issues section of this repo. In fact, you should go to the Issues section with any question you have about bakR, and there are even helpful labels that you can append to your posts to make the nature of your request clear. If you email me (Isaac Vock) with a question/concern/suggestion, I will direct you to the Issues section. If you have basic use questions, I would suggest going through the vignettes linked above. If these do not answer your question, then post your question to Issues.

