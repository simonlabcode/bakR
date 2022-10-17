# Brief Description of bakR
bakR (Bayesian analysis of the kinetics of RNA) is an R package for performing differential kinetics analysis with nucleotide recoding high-throughput RNA sequencing (NR-seq) data. 
Kinetic parameter estimation and statistical testing is compatible with mutational data from any enrichment free NR-seq method (e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.).

# Why use bakR?
Differential expression analysis of RNA sequencing (RNA-seq) data can identify changes in cellular RNA levels, but cannot determine the kinetic mechanism underlying such changes. Previously, [our lab](https://simonlab.yale.edu/research/transcriptome-dynamics/timelapse-chemistry/) and others addressed this shortcoming by developing nucleotide-recoding RNA-seq methods (NR-seq; e.g., TimeLapse-seq) to quantify changes in RNA synthesis and degradation kinetics. While advanced statistical models implemented in user-friendly software (e.g., DESeq2) have ensured the statistical rigor of differential expression analyses, no such tools that facilitate differential kinetic analysis with NR-seq exist. To address this need, we developed bakR, an R package that analyzes and compares NR-seq datasets. Differential kinetics analysis with bakR relies on a Bayesian hierarchical model of NR-seq data to increase statistical power by sharing information across transcripts. bakR outperforms attempts to use single sample analysis tools (e.g., pulseR and GRAND-SLAM) for differential kinetics analysis. Check out [our preprint](https://www.biorxiv.org/content/10.1101/2022.09.02.505697v1) to learn more about the model and its extensive validation!

# Installation
bakR is now available on CRAN! If you are using a Windows or Mac OS and R version 4.2.1 (or the newest develop version of R), then that means you don't need to configure a C++ compiler to install and use bakR. Eventually, this will also hold for Windows and Mac OS users with any version of R, CRAN just needs time to build the binaries for older R versions. If you are hoping to use bakR on any other operating system or are working with an older version of R, you need to first configure the C++ compiler (see the next paragraph for details and links). In either case, bakR can be installed as follows:

    install.packages("bakR") 

To install the newest version of bakR from Github, you need to have a C++ compiler configured to rstan's (the R interface to the probabilistic programming language [Stan](https://mc-stan.org/) that bakR uses on the backend) liking. The best way to do this is to follow the Stan team's [helpful documentation](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) on installing rstan for your operating system. Once that is complete, you can install bakR as follows:

    install.packages("devtools") # if you haven't installed devtools already
    devtools::install_github("simonlabcode/bakR")

# Documentation
There are currently two vignettes to help get you up to speed with using bakR:

  1. An introductory vignette (title: Differential Kinetic Analysis with bakR) that walks you through the basic bakR workflow with simulated data.
  2. A second vignette (title: Going beyond the defaults with bakR) that dives into some specific follow-up analyses one can perform after running bakR on a dataset. Currently, this vignette discusses differential synthesis analysis (i.e., identifying differences in RNA synthesis rates by combining bakR with differential expression analysis) and performing analyses without assuming steady-state. More analyses will be added to this vignette as time goes by.
  
All vignettes are available on the [bakR website](https://simonlabcode.github.io/bakR/index.html) under the Articles section. [Here](https://github.com/simonlabcode/bakR) is the link to the bakR github as well if you need help getting back to the github from the website.

# Obtaining the Necessary Input
As discussed in the introductory vignette, bakR requires data in the form of a so-called "cB", or counts binomial data frame. Each row of the cB data frame corresponds to a group of reads with identical mutational data, and the columns denote the sample from which the reads came, the feature the reads aligned to, the number of mutations of interest in the reads (e.g., T-to-C mutations), the number of mutable positions (e.g. Ts), and the number of such reads. It is reasonable to wonder "where am I supposed to get this information?" While there are a couple possibilities, perhaps the easiest and most widely applicable is [bam2bakR](https://github.com/simonlabcode/bam2bakR), a Snakemake implementation of the [TimeLapse pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by the Simon lab. bam2bakR takes as input aligned bam files and produces, among other things, the cB file required by bakR. Extensive documentation describing how to get bam2bakR up and running is available on its GitHub repo. Snakemake greatly facilitates running this pipeline on almost any computational infrastructure and bam2bakR uses the conda/mamba package manager to make setting up the necessary dependencies a breeze. Also look out for fastq2bakR in the coming days, a similar tool that additionally takes care of any necessary fastq preprocessing and alignment!

# Bug Catching and Further Questions
Post descriptions of bugs and a simple reproducible example (if possible) in the Issues section of this repo. In fact, you should go to the Issues section with any question you have about bakR, and there are even helpful labels that you can append to your posts to make the nature of your request clear. If you email me (Isaac Vock) with a question/concern/suggestion, I will direct you to the Issues section. If you have basic use questions, I would suggest going through the vignettes linked above. If these do not answer your question, then post your question to Issues.

