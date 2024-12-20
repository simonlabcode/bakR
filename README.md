# Brief Description of bakR <img src='man/figures/logo.png' align="right" height="225" />
bakR (Bayesian analysis of the kinetics of RNA) is an R package for performing differential kinetics analysis with nucleotide recoding high-throughput RNA sequencing (NR-seq) data. 
Kinetic parameter estimation and statistical testing is compatible with mutational data from any enrichment free NR-seq method (e.g., TimeLapse-seq, SLAM-seq, TUC-seq, etc.).

# Update: EZbakR, a complete rewrite of bakR, is out now! (10/14/2024)

If you are currently using bakR, we highly suggest checking out [EZbakR](https://github.com/isaacvock/EZbakR), a far more flexible alternative. EZbakR can do everything bakR can, and a lot more! EZbakR is accompanied by a similar improvement/extension of bam2bakR, called [fastq2EZbakR](https://github.com/isaacvock/fastq2EZbakR). Check out [our preprint](https://www.biorxiv.org/content/10.1101/2024.10.14.617411v1) for more details. While we will continue to maintain bakR, all future development will occur on EZbakR.

# Version 1.0.0 is out now! (06/27/2023)
A lot of functionality has been added, and I highly suggest all users of bakR to update to this version. There are also many new vignettes to discuss these new features. bakR v1.0.0 is now available for installation on CRAN! It is also currently available for installation from Github, as described below. Two major new additions are:

1. Ability to use [GRAND-SLAM](https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM) output (or fraction new estimates more generally) as bakR input
2. Strategy for correcting metabolic label related biases in kinetic parameter estimates and read counts 

# Why use bakR?
Differential expression analysis of RNA sequencing (RNA-seq) data can identify changes in cellular RNA levels, but cannot determine the kinetic mechanism underlying such changes. Previously, [our lab](https://simonlab.yale.edu/research/transcriptome-dynamics/timelapse-chemistry/) and others addressed this shortcoming by developing nucleotide-recoding RNA-seq methods (NR-seq; e.g., TimeLapse-seq) to quantify changes in RNA synthesis and degradation kinetics. While advanced statistical models implemented in user-friendly software (e.g., DESeq2) have ensured the statistical rigor of differential expression analyses, no such tools that facilitate differential kinetic analysis with NR-seq exist. To address this need, we developed bakR, an R package that analyzes and compares NR-seq datasets. Differential kinetics analysis with bakR relies on a Bayesian hierarchical model of NR-seq data to increase statistical power by sharing information across transcripts. bakR outperforms attempts to use single sample analysis tools (e.g., pulseR and GRAND-SLAM) for differential kinetics analysis. Check out [our manuscript in RNA](https://rnajournal.cshlp.org/content/29/7/958.abstract) to learn more about the model and its extensive validation!

# Installation
bakR is now available on CRAN! If you are using a Mac or Windows OS then that means you don't need to configure a C++ compiler to install and use bakR. Those not on a Mac Windows OS will need to first properly configure a C++ compiler; see the next paragraph for details and links describing how to do that. In either case, once you (and your compiler if necessary) are ready, bakR can be installed as follows:

    install.packages("bakR") 

To install the newest version of bakR from Github, you need to have a C++ compiler configured to rstan's (the R interface to the probabilistic programming language [Stan](https://mc-stan.org/) that bakR uses on the backend) liking. The best way to do this is to follow the Stan team's [helpful documentation](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) on installing rstan for your operating system. Once that is complete, you can install bakR as follows:

    install.packages("devtools") # if you haven't installed devtools already
    devtools::install_github("simonlabcode/bakR")

# Documentation
There are currently seven vignettes to help get you up to speed with using bakR:

  1. [An introductory vignette](https://simonlabcode.github.io/bakR/articles/Getting-Started.html) (title: Differential Kinetic Analysis with bakR) that walks you through the basic bakR workflow with simulated data.
  2. [A more concise version of the introductory vignette](https://simonlabcode.github.io/bakR/articles/bakR-Quickstart.html) that will get you up and running with bakR quickly (title: bakR for people in a hurry). Particularly appropriate for those who are very comfortable with adopting new bioinformatic tools.
  3. [Combining bakR with differential expression analysis](https://simonlabcode.github.io/bakR/articles/Differential-Synth.html) to perform differential synthesis rate analysis (title: Differential synthesis analysis with bakR and DESeq2).
  4. [How to use fraction new estimates (e.g., from a tool like GRAND-SLAM) as input to bakR](https://simonlabcode.github.io/bakR/articles/bakR-Fn.html), a new feature introduced in version 1.0.0 (title: GRAND-SLAM output/fn estimates as bakR input).
  5. [Correcting for disproportionate loss of s<sup>4</sup>U containing RNA](https://simonlabcode.github.io/bakR/articles/Dropout.html) (title: Correcting for dropout). This phenomenon, termed dropout, is discussed in two recent preprints, one from [our lab](https://www.biorxiv.org/content/10.1101/2023.05.24.542133v1) and one from the [Erhard lab](https://www.biorxiv.org/content/10.1101/2023.04.21.537786v1.full).
  6. [How to identify and deal with problems](https://simonlabcode.github.io/bakR/articles/Troubleshooting.html) that can crop up when analyzing NR-seq data (title: Troubleshooting analyses of NR-seq data with bakR).
  7. [Distinguishing transcriptional and post-transcriptional regulation](https://simonlabcode.github.io/bakR/articles/NSS.html), even when the steady-state assumption is partially violated (title: Steady-state quasi-independent mechanistic investigations). Describes a new and somewhat experimental function in bakR, `DissectMechanism`.
  
All vignettes are available on the [bakR website](https://simonlabcode.github.io/bakR/index.html) under the Articles section. [Here](https://github.com/simonlabcode/bakR) is the link to the bakR github as well if you need help getting back to the github from the website.

# Obtaining the Necessary Input
As discussed in the introductory vignette, bakR requires data in the form of a so-called "cB", or counts binomial data frame. Each row of the cB data frame corresponds to a group of reads with identical mutational data, and the columns denote the sample from which the reads came, the feature the reads aligned to, the number of mutations of interest in the reads (e.g., T-to-C mutations), the number of mutable positions (e.g. Ts), and the number of such reads. It is reasonable to wonder "where am I supposed to get this information?" While there are a couple possibilities, perhaps the easiest and most widely applicable is [bam2bakR](https://github.com/simonlabcode/bam2bakR), a Snakemake implementation of the [TimeLapse pipeline](https://bitbucket.org/mattsimon9/timelapse_pipeline/src/master/) developed by the Simon lab. bam2bakR takes as input aligned bam files and produces, among other things, the cB file required by bakR. Extensive documentation describing how to get bam2bakR up and running is available on its GitHub repo. Snakemake greatly facilitates running this pipeline on almost any computational infrastructure and bam2bakR uses the conda/mamba package manager to make setting up the necessary dependencies a breeze.

As of version 1.0.0, bakR can also take as input fraction new (sometimes referred to as new-to-total ratio, or NTR) estimates. These are obtainable via tools like [GRAND-SLAM](https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM), or perhaps a custom analysis pipeline that you developed while working with NR-seq datasets!

# Bug Catching and Further Questions
Post descriptions of bugs and a simple reproducible example (if possible) in the Issues section of this repo. In fact, you should go to the Issues section with any question you have about bakR, and there are even helpful labels that you can append to your posts to make the nature of your request clear. If you email me (Isaac Vock) with a question/concern/suggestion, I will direct you to the Issues section. If you have basic use questions, I would suggest going through the vignettes linked above. If these do not answer your question, then post your question to Issues.

