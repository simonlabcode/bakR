
#' Simulating nucleotide recoding data
#'
#' \code{sim_bakRData} simulates a `bakRData` object. It's output also includes the simulated
#' values of all kinetic parameters of interest. Only the number of genes (\code{ngene}) has to be set by the
#' user, but an extensive list of additional parameters can be adjusted.
#'
#' \code{sim_bakRData} simulates a `bakRData` object using an extensive list of
#' adjustable parameters. Average RNA kinetic parameters are drawn from biologically inspired
#' distributions. Replicate variability is simulated by drawing each replicate and feature's
#' fraction labeled (aka fraction new) from a Logit-Normal distribution with a heteroskedastic
#' variance term with average magnitude given by a read count vs. variance relationship, and each
#' replicate and feature's ksyn from a homoskedastic lognormal distribution. Read counts
#' can either be set to the same value for all simulated features or can be simulated according to
#' a heterodisperse negative binomial distribution.
#'
#' The number of Us in each sequencing read is drawn from a binomial distribution with number of trials
#' equal to the read length and probability of each nucleotide being a U of 0.25. Each read is assigned to the
#' labeled or unlabeled population according to a Bernoulli distribution with p = fraction new. The number of
#' mutations in each read are then drawn from one of two binomial distributions; if the read is assigned to the
#' labeled population, the number of mutations are drawn from a binomial distribution with number of trials equal
#' to the number of Us and probability of mutation = \code{p_new}; if the read is assigned to the unlabeled population,
#' the number of mutations is instead drawn from a binomial distribution with the same number of trials but with the probability
#' of mutation = \code{p_old}.
#'
#' Simulated read counts should be treated as if they are spike-in and RPKM normalized, so the same scale factor can be applied
#' to each sample when comparing the sequencing reads (like if you are performing differential expression analysis).
#'
#' Function to simulate a bakRData object according to a heteroskedastic beta-binomial model
#' of the fraction new.
#' @param ngene Number of genes to simulate data for
#' @param num_conds Number of experimental conditions (including the reference condition) to simulate
#' @param nreps Number of replicates to simulate
#' @param eff_sd Effect size; more specifically, the standard deviation of the normal distribution from which non-zero
#' changes in logit(fraction new) are pulled from.
#' @param eff_mean Effect size mean; mean of normal distribution from which non-zero changes in logit(fraction new) are pulled from.
#' Note, setting this to 0 does not mean that some of the significant effect sizes will be 0, as any exact integer is impossible
#' to draw from a continuous random number generator. Setting this to 0 just means that there is symmetric stabilization and destabilization
#' @param fn_mean Mean of fraction news of simulated transcripts in reference condition. The fraction of RNA from each transcript that is
#' s4U labeled (new) is drawn from a normal distribution with this mean
#' @param fn_sd Standard deviation of fraction news of simulated transcripts in reference condition. The fraction of RNA from each transcript that
#' is s4U labeled (new) is drawn from a normal distribution with this sd
#' @param kslog_c Synthesis rate constants will be drawn from a lognormal distribution with meanlog = \code{kslog_c} - mean(log(kd_mean)) where kd_mean
#' is determined from the fraction new simulated for each gene as well as the label time (\code{tl}).
#' @param kslog_sd Synthesis rate lognormal standard deviation; see kslog_c documentation for details
#' @param tl s4U label feed time
#' @param p_new s4U induced mutation rate. Can be a vector of length num_conds
#' @param p_old background mutation rate
#' @param read_lengths Total read length for each sequencing read (e.g., 200 means PE100 reads)
#' @param p_do Rate at which s4U containing reads are lost due to dropout; must be between 0 and 1
#' @param noise_deg_a Slope of trend relating log10(standardized read counts) to log(replicate variability)
#' @param noise_deg_b Intercept of trend relating log10(standardized read counts) to log(replicate variability)
#' @param noise_synth Homoskedastic variability of L2FC(ksyn)
#' @param sd_rep Variance of lognormal distribution from which replicate variability is drawn
#' @param low_L2FC_ks Most negative L2FC(ksyn) that can be simulated
#' @param high_L2FC_ks Most positive L2FC(ksyn) that can be simulated
#' @param num_kd_DE Vector where each element represents the number of genes that show a significant change in stability relative
#' to the reference. 1st entry must be 0 by definition (since relative to the reference the reference sample is unchanged)
#' @param num_ks_DE Same as num_kd_DE but for significant changes in synthesis rates.
#' @param scale_factor Factor relating RNA concentration (in arbitrary units) to average number of read counts
#' @param sim_read_counts Logical; if TRUE, read counts are simulated as coming from a heterodisperse negative binomial distribution
#' @param a1 Heterodispersion 1/reads dependence parameter
#' @param a0 High read depth limit of negative binomial dispersion parameter
#' @param nreads Number of reads simulated if sim_read_counts is FALSE
#' @param alpha shape1 parameter of the beta distribution from which U-contents (probability that a nucleotide in a read from a transcript is a U) are
#' drawn for each gene.
#' @param beta shape2 parameter of the beta distribution from which U-contents (probability that a nucleotide in a read from a transcript is a U) are
#' drawn for each gene.
#' @param STL logical; if TRUE, simulation is of STL-seq rather than a standard TL-seq experiment. The two big changes are that a short read lenght is required
#' (< 60 nt) and that every read for a particular feature will have the same number of Us. Only one read length is simulated for simplicity.
#' @param STL_len Average length of simulated STL-seq length. Since Pol II typically pauses about 20-60 bases
#' from the promoter, this should be around 40
#' @export
#' @return A list containing a simulated `bakRData` object as well as a list of simulated kinetic parameters of interest.
#' The contents of the latter list are:
#' \itemize{
#'  \item Effect_sim; Dataframe meant to mimic formatting of Effect_df that are part of \code{TL_stan} and \code{fast_analysis} output.
#'  \item Fn_mean_sim; Dataframe meant to mimic formatting of Regularized_ests that is part of \code{fast_analysis} output. Contains information
#'  about the true fraction new simulated in each condition (the mean of the normal distribution from which replicate fraction news are simulated)
#'  \item Fn_rep_sim; Dataframe meant to mimic formatting of Fn_Estimates that is part of \code{fast_analysis} output. Contains information
#'  about the fraction new simulated for each feature in each replicate of each condition.
#'  \item L2FC_ks_mean; The true L2FC(ksyn) for each feature in each experimental condition. The ith column corresponds to the L2FC(ksyn) when comparing
#'  the ith condition to the reference condition (defined as the 1st condition) so the 1st column is always all 0s
#'  \item RNA_conc; The average number of normalized read counts expected for each feature in each sample.
#' }
#'
sim_bakRData <- function(ngene, num_conds = 2L, nreps = 3L, eff_sd = 0.75, eff_mean = 0,
                               fn_mean = 0, fn_sd = 1, kslog_c = 0.8, kslog_sd = 0.95,
                               tl = 60, p_new = 0.05, p_old = 0.001, read_lengths = 200L,
                               p_do = 0, noise_deg_a = -0.3, noise_deg_b = -1.5, noise_synth = 0.1, sd_rep = 0.05,
                               low_L2FC_ks = -1, high_L2FC_ks = 1,
                               num_kd_DE = c(0L, as.integer(rep(round(as.integer(ngene)/2), times = as.integer(num_conds)-1))),
                               num_ks_DE = rep(0L, times = as.integer(num_conds)),
                               scale_factor = 150,
                               sim_read_counts = TRUE, a1 = 5, a0 = 0.01,
                               nreads = 50L, alpha = 25, beta = 75,
                               STL = FALSE, STL_len = 40){


  .Deprecated("Simulate_bakRData")

  ### Catch non-sensical values in all inputs

  # STL
  if(!is.logical(STL)){
    "STL must be logical (TRUE or FALSE)"
  }

  # STL_len
  if(STL){
    if(STL_len < 25){
      stop("STL_len is less than 25. Must be between 25 and 55.")
    }else if(STL_len > 55){
      stop("STL_len is greater than 55. Must be between 25 and 55.")
    }
  }

  # alpha parameter of beta distribution
  if(!is.numeric(alpha)){
    stop("alpha must be numeric!")
  }else if(alpha <= 1){
    stop("alpha must be > 1")
  }else if(!is.numeric(beta)){
    stop("beta must be numeric!")
  }else if(beta <=1){
    stop("beta must be > 1")
  }else if(alpha/(alpha + beta) < 0.1){
    warning("alpha and beta are such that the average U-content is less than 0.1. That is unreasonably low and may make many simulated transcripts unanalyzable.")
  }else if(alpha/(alpha + beta) > 0.75){
    warning("alpha and beta are such that the average U-content is > 0.75. I wish U-contents were this high, your simulation may not reflect real data though.")
  }

  # Number of genes (features) to simulate
  if(!is.numeric(ngene)){
    stop("ngene must be numeric")
  }else if(!is.integer(ngene)){
    ngene <- as.integer(ngene)
  }

  if(ngene < 1){
    stop("ngene must be > 0; it represents the number of features to simulate")
  }

  # Number of experimental conditions


  if(!is.numeric(num_conds)){
    stop("num_conds must be an numeric")
  }else if(!is.integer(num_conds)){
    num_conds <- as.integer(num_conds)
  }

  if(num_conds < 1){
    stop("num_conds must be > 0; it represents the number of experimental conditions")
  }else if(num_conds == 1){
    warning("You are only simluating a reference condition with no experimental conditions.")
  }

  # Number of replicates
  if(!is.numeric(nreps)){
    stop("nreps must be numeric")
  }else if(!is.integer(nreps)){
    nreps <- as.integer(nreps)
  }

  if(nreps < 1){
    stop("nreps must be > 0; it represents the number of replicates")
  }else if(nreps == 1){
    warning("You are only simulating a single replicate; All statistical models implemented in bakR require > 1 replicate.")
  }

  # Effect size distribution standard deviation
  if(!is.numeric(eff_sd)){
    stop("eff_sd must be numeric")
  }else if(eff_sd <= 0){
    stop("eff_sd must be > 0; it will be the sd parameter of a call to rnorm when simulating effect sizes")
  }else if (eff_sd > 4){
    warning("You are simulating an unusually large eff_sd (> 4)")
  }


  # Fraction new distribution sd
  if(!is.numeric(fn_sd)){
    stop("fn_sd must be numeric")
  }else if(fn_sd <= 0){
    stop("fn_sd must be > 0; it will be the sd parameter of a call to rnorm when simulating reference fraction news")
  }else if (fn_sd > 2){
    warning("You are simulating an unusually large fn_sd (> 2). This will lead to lots of extreme fraction news (close to 0 and 1).")
  }

  # Label time
  if(!is.numeric(tl)){
    stop("tl must be numeric")
  }else if(tl <= 0){
    stop("tl must be > 0; it represents the label time in minutes.")
  }

  # s4U mutation rate
  if(!all(is.numeric(p_new))){
    stop("p_new must be numeric")
  }else if(!all(p_new > 0)){
    stop("p_new must be > 0; it represents the mutation rate of s4U labeled transcripts")
  }else if(!all(p_new <= 1)){
    stop("p_new must be <= 1; it represents the mutation rate of s4U labeled transcripts")
  }else if(!all(p_new < 0.2)){
    warning("You have simulated an unusually large s4U induced mutation rate (>= 0.2)")
  }else if(!all((p_new - p_old) > 0)){
    stop("All p_new must be > p_old, since the background mutation rate (p_old) should always be less than the labeled mutation rate (p_new)")
  }

  # Background mutation rate
  if(!all(is.numeric(p_old))){
    stop("p_old must be numeric")
  }else if(!all(p_old > 0)){
    stop("p_old must be > 0; it represents the background mutation rate")
  }else if(!all(p_old <= 1)){
    stop("p_old must be <= 1; it represents the background mutation rate")
  }else if(!all(p_old < 0.01)){
    warning("You have simulated an unusually large background mutation rate (>= 0.01)")
  }

  # Read length
  if(!all(is.numeric(read_lengths))){
    stop("read_lengths must be numeric")
  }else if(!all(is.integer(read_lengths))){
    read_lengths <- as.integer(read_lengths)
  }

  if(!all(read_lengths > 0)){
    stop("read_lengths must be > 0; it represents the total length of sequencing reads")
  }else if(!all(read_lengths >= 20)){
    warning("You have simulated an unusually short read length (< 20 nucleotides)")
  }

  # Dropout probability
  if(!all(is.numeric(p_do))){
    stop("p_do must be numeric")
  }else if(!all(p_do <= 1)){
    stop("p_do must be <= 1; it represents the percentage of s4U containing RNA lost during library prep")
  }else if(!all(p_do >= 0)){
    stop("p_do must be >= 0; it represents the percentage of s4U containing RNA lost during library prep")
  }else if(!all(p_do == 0)){
    warning("You are simulating dropout; statistical models implemented by bakR do not correct for dropout and thus will provide biased estimates")
  }

  # Heteroskedastic Slope
  if(!is.numeric(noise_deg_a)){
    stop("noise_deg_a must be numeric")
  }else if(noise_deg_a > 0){
    stop("noise_deg_a must be < 0; it represents the slope of the log10(read depth) vs. log(replicate variability) trend")
  }else if(noise_deg_a == 0){
    warning("You are simulating fraction new homoskedasticity, which is not reflective of actual nucleotide recoding data")
  }

  # Synthesis variability
  if(!is.numeric(noise_synth)){
    stop("noise_synth must be numeric")
  }else if(noise_synth < 0){
    stop("noise_synth must be >= 0; it represents the homoskeastic variability in the synthesis rate")
  }

  # Replicate variability variability
  if(!is.numeric(sd_rep)){
    stop("sd_rep must be numeric")
  }else if(sd_rep < 0){
    stop("sd_rep must be >= 0; it represents the sdlog parameter of the log-normal distribution from which replicate variabilites are simulated")
  }else if(sd_rep == 0){
    warning("You are simulating no variability in the replicate variability. This can cause minor convergence issues in the completely hierarchical
            models implemented by bakR.")
  }

  # L2FC(ksyn) Bounds
  if(!all(is.numeric(c(high_L2FC_ks, low_L2FC_ks)))){
    stop("high_L2FC_ks and low_L2FC_ks must be numeric")
  }else if(high_L2FC_ks < low_L2FC_ks){
    stop("high_L2FC_ks must be greater than low_L2FC_ks; they represent the upper and lower bound of a uniform distribution respectively")
  }

  # Number of differentially stabilized features
  if(!all(is.numeric(num_kd_DE))){
    stop("All elements of num_kd_DE must be numeric")
  }else if(!all(is.integer(num_kd_DE))){
    num_kd_DE <- as.integer(num_kd_DE)
  }

  if (length(num_kd_DE) > num_conds){
    stop("num_kd_DE has too many elements; it should be a vector of length == num_conds")
  }else if (length(num_kd_DE) < num_conds){
    stop("num_kd_DE has too few elements; it should be a vector of length == num_conds")
  }else if(num_kd_DE[1] != 0){
    stop("The 1st element of num_kd_DE must equal 0 by definition, as it represents the number of features differentially stabilized in the reference
         condition relative to the reference condition")
  }else if(!all(num_kd_DE <= ngene)){
    stop("Not all elements of num_kd_DE are less than or equal to the total number of simulated features")
  }else if(!all(num_kd_DE >= 0)){
    stop("Not all elements of num_kd_DE are > 0")
  }

  # Number of differentially synthesized features
  if(!all(is.numeric(num_ks_DE))){
    stop("All elements of num_ks_DE must be numeric")
  }else if(!all(is.integer(num_ks_DE))){
    num_ks_DE <- as.integer(num_ks_DE)
  }

  if(length(num_ks_DE) > num_conds){
    stop("num_ks_DE has too many elements; it should be a vector of length == num_conds")
  }else if (length(num_ks_DE) < num_conds){
    stop("num_ks_DE has too few elements; it should be a vector of length == num_conds")
  }else if(num_ks_DE[1] != 0){
    stop("The 1st element of num_ks_DE must equal 0 by definition, as it represents the number of features differentially stabilized in the reference
         condition relative to the reference condition")
  }else if(!all(num_ks_DE <= ngene)){
    stop("Not all elements of num_ks_DE are less than or equal to the total number of simulated features")
  }else if(!all(num_ks_DE >= 0)){
    stop("Not all elements of num_ks_DE are > 0")
  }

  # Scale factor
  if(!is.numeric(scale_factor)){
    stop("scale_factor must be numeric")
  }else if(scale_factor <=0){
    stop("scale_factor must be > 0; it represents the factor by which the RNA concentration is multiplied to yield the average number of sequencing reads")
  }else if(scale_factor < 10){
    warning("You are simulating an unusually low scale_factor (< 10); small scale factors will lead to low read counts.")
  }

  # Simulate read count Boolean
  if(!is.logical(sim_read_counts)){
    stop("sim_read_counts must be a logical (TRUE or FALSE)")
  }else if(sim_read_counts){
    if(!is.integer(nreads)){
      stop("nreads must be an integer")
    }else if(nreads <= 0){
      stop("nreads must be > 0; it represents the number of sequencing reads to be simulated for all features")
    }
  }

  # Heterodispersion parameters

  if(!all(is.numeric(c(a1, a0)))){
    stop("a1 and a0 must be numeric")
  }else if(a1 < 0){
    stop("a1 must be >= 0; it relates the average read count to the negative binomial dispersion parameter")
  }else if(a0 <= 0){
    stop("a0 must be > 0; it represents the high read count limit of the negative binomial dispersion parameter")
  }

  # Parameters without numeric bounds
  if(!is.numeric(noise_deg_b)){
    stop("noise_deg_b must be numeric")
  }

  if(!is.numeric(eff_mean)){
    stop("eff_mean must be numeric")
  }



  if(length(p_new) == 1){
    p_new <- rep(p_new, times=num_conds)
  }else if(length(p_new) != num_conds){
    stop("p_new must be of length 1 or length == num_conds")
  }

  if(length(p_old) == 1){
    p_old <- rep(p_old, times=num_conds)
  }else if(length(p_old) != num_conds){
    stop("p_old must be of length 1 or length == num_conds")
  }

  if(length(read_lengths) == 1){
    read_lengths <- rep(read_lengths, times=num_conds)
  }else if(length(read_lengths) != num_conds){
    stop("read_lengths must be of length 1 or length == num_conds")
  }

  if(length(p_do) == 1){
    p_do <- rep(p_do, times=num_conds)
  }else if(length(p_do) != num_conds){
    stop("p_do must be of length 1 or length == num_conds")
  }


  U_cont <- stats::rbeta(ngene, alpha, beta)

  # Define helper functions:
  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))

  #Initialize matrices
  fn <- rep(0, times=ngene*nreps*num_conds)
  dim(fn) <- c(ngene, num_conds, nreps)

  Counts <- fn


  kd <- fn
  ks <- fn


  #Initialize vectors of mean values for each gene and condition
  fn_mean <- inv_logit(stats::rnorm(n=ngene, mean=fn_mean, sd=fn_sd))
  kd_mean <- -log(1-fn_mean)/tl
  ks_mean <- stats::rlnorm(n=ngene, meanlog = kslog_c + log(kd_mean), sdlog = kslog_sd)

  #ks_mean <- rexp(n = ngene, rate = (1/kd_mean)*5)

  effect_mean <- rep(0, times = ngene*num_conds)
  dim(effect_mean) <- c(ngene, num_conds)
  L2FC_ks_mean <- effect_mean
  L2FC_kd_mean <- effect_mean

  for(i in 1:ngene){

    #Make sure the user didn't input the wrong
    #number of significant genes
    if (i == 1 ){
      if (length(num_kd_DE) > num_conds){
        stop("num_kd_DE has too many elements")
      } else if (length(num_kd_DE) < num_conds){
        stop("num_kd_DE has too few elements")
      }

      if (length(num_ks_DE) > num_conds){
        stop("num_ks_DE has too many elements")
      } else if (length(num_ks_DE) < num_conds){
        stop("num_ks_DE has too few elements")
      }
    }

    for(j in 1:num_conds){
      if(j == 1){
        effect_mean[i,1] <- 0
        L2FC_ks_mean[i,1] <- 0
      }else{
        if(i < (ngene-num_kd_DE[j] + 1)){
          effect_mean[i,j] <- 0
        }else{
          effect_mean[i,j] <- stats::rnorm(n=1, mean=eff_mean, sd=eff_sd)
        }
        if(i < (ngene-num_ks_DE[j] + 1)){
          L2FC_ks_mean[i,j] <- 0
        }else{
          if (stats::runif(1) < 0.5){
            L2FC_ks_mean[i,j] <- stats::runif(n=1, min=low_L2FC_ks, max=high_L2FC_ks)
          }else{
            L2FC_ks_mean[i,j] <- stats::runif(n=1, min=-high_L2FC_ks, max=-low_L2FC_ks)
          }

        }
      }
    }
  }

  L2FC_kd_mean <- log2(log(1 - inv_logit(fn_mean + effect_mean))/log(1- inv_logit(fn_mean)))


  #Simulate read counts
  if (sim_read_counts == TRUE){
    L2FC_tot_mean <- L2FC_ks_mean - L2FC_kd_mean
    RNA_conc <- (ks_mean*2^(L2FC_ks_mean))/(kd_mean*2^(L2FC_kd_mean))
    a1 <- 5
    a0 <- 0.01

    for(i in 1:ngene){
      for(j in 1:num_conds){
        for(k in 1:nreps){
          Counts[i, j, k] <- stats::rnbinom(n=1, size=1/((a1/(scale_factor*RNA_conc[i,j])) + a0), mu = scale_factor*RNA_conc[i,j])
          #Counts[i, j, k] <- rpois(n=1, lambda = scale_factor*RNA_conc[i,j,k])

          if(Counts[i, j, k] < 5){
            Counts[i, j, k] <- Counts[i, j, k] + stats::rpois(n=1, lambda = 2) + 1
          }
        }
      }
    }
  } else{
    if(sim_from_data == FALSE){
      Counts <- rep(nreads, times= ngene*num_conds*nreps)
      dim(Counts) <- c(ngene, num_conds, nreps)
    }

  }

  standard_RNA <- matrix(0, nrow = ngene, ncol = num_conds)
  mean_RNA <- rep(0, times = num_conds)
  sd_RNA <- rep(0, times = num_conds)

  #SIMULATE L2FC OF DEG AND SYNTH RATE CONSTANTS; REPLICATE VARIABILITY SIMULATED
  for(j in 1:num_conds){
    mean_RNA[j] <- mean(log10(RNA_conc[,j]*scale_factor))
    sd_RNA[j] <- stats::sd(log10(RNA_conc[,j]*scale_factor))
    for(i in 1:ngene){
      standard_RNA[i,j] <- (log10(RNA_conc[i,j]*scale_factor) - mean_RNA[j])/sd_RNA[j]
      for(k in 1:nreps){
        fn[i, j, k] <- inv_logit(stats::rnorm(1, mean=(logit(fn_mean[i]) + effect_mean[i,j]), sd = stats::rlnorm(1, noise_deg_a*standard_RNA[i,j] + noise_deg_b, sd_rep )))
        ks[i,j,k] <- exp(stats::rnorm(1, mean=log((2^L2FC_ks_mean[i,j])*ks_mean[i]), sd=noise_synth))
      }
    }
  }


  kd <- -log(1 - fn)/tl



  l <- ngene

  p_do <- matrix(rep(0, times = num_conds*nreps), nrow = nreps, ncol = num_conds)
  fn_real <- fn

  for(j in 1:num_conds){
    for(k in 1:nreps){
      fn_real[,j,k] <- (fn[,j,k]*(1-p_do[k,j]))/(1 - p_do[k,j]*(1 - fn[,j,k]))
      Counts[,j,k] <- Counts[,j,k] - Counts[,j,k]*fn[,j,k]*p_do[k,j]
    }
  }



  # start_1 <- Sys.time()

  # This is one very huge function
  # It simulates TL-seq data, recording the number of TC mutations in each read, which is informed
  # by whatever the fraction new for the particular transcript is
  simulateData <- function(nmir = l,   # num of genes
                           fn_s4U = fn_real,   # fraction of s4U reads made after label introduction in non-heatshocked sample
                           #fn_s4U2 = fn_hs, # fraction of s4U reads made after label introduction in heatshocked sample
                           p_new_real_tc = p_new,                   # TC mutation rate in fed cells
                           p_old_real_tc = p_old,                  # TC mutation rate in unfed cells
                           read_length = read_lengths,
                           nreads = Counts, # per transript per sample
                           nsamp = (nreps*num_conds) + num_conds,
                           ctl = c(rep(1, times=nreps*num_conds), rep(0, times=num_conds)),   # cntl = 0 is no feed, cntl = 1 is feed
                           mt = c(rep(1:num_conds, each=nreps),seq(from=1,to=num_conds,by=1)),
                           replicate = c(rep(seq(from=1, to=nreps), times=num_conds), rep(1, times=num_conds))
                           #Could just generalize this, which is what next line does, repeating 1 for all but the last sample
                           #ctl = c(rep(1, times = nsamp-1),0) #assumes nsamp is odd, think it has to be
  ){



    # Start generating a vector with data
    sample_data <- vector('list', length = nsamp)
    mir_data <- vector('list', length = nmir)

    for (s in 1:nsamp){
      for (mir in 1:nmir){ #mir is feature number index, should change to gene or something

        r <- replicate[s] #Replicate number index
        MT <- mt[s]       #Experimental sample index
        readsize = read_length[MT]

        mir_pold_tc <- p_old[MT]
        mir_pnew_tc <- p_new[MT]

        #Simulate which reads are labeled
        newreads_tc <- purrr::rbernoulli(as.numeric(nreads[mir, MT, r]), p = as.numeric(fn_s4U[mir, MT, r]))# vector of reads, T/F is s4U labeled

        #Simulate the nubmer of Us in each read
        if(STL){
          nu <- abs(round(U_cont[mir]*STL_len) +  sign(stats::runif(nreads[mir, MT, r], min = -0.1, max = 0.1))*stats::rpois(nreads[mir, MT, r], 0.5))
        }else{
          nu <- stats::rbinom(n = nreads[mir,MT,r], size = readsize, prob = U_cont[mir])
        }

        #Number of reads that are new
        newreads_tc <- sum(newreads_tc)

        #Generate number of mutations for new and old reads
        if (!ctl[s]){ #If no s4U added, only old
          nmut_tc <- stats::rbinom(n = nreads[mir,MT,r], size=nu, prob = mir_pold_tc)
        }else {
          nmut_tc_new <- stats::rbinom(n=newreads_tc, size=nu[1:newreads_tc], prob=mir_pnew_tc)
          nmut_tc_old <- stats::rbinom(n=(nreads[mir, MT, r]-newreads_tc), size = nu[(newreads_tc+1):nreads[mir, MT, r]], prob=mir_pold_tc)
          nmut_tc <- c(nmut_tc_new, nmut_tc_old)
        }


        #Now make it look kind of like a cB file
        # use mirMut for each gene, cntl is if cntl or not, x is # of new reads
        df <- dplyr::tibble(S = rep(s, times = nreads[mir, MT, r]),  #starting to generate something that looks like a cB file
                     TP = rep(ctl[s], times = nreads[mir, MT, r]),
                     R = rep(r, times=nreads[mir, MT, r]),
                     MIR = rep(mir, times = nreads[mir, MT, r]), # same as XF or fnum, so just feature number
                     TC = nmut_tc,
                     MT = rep(MT, times=nreads[mir, MT, r]),
                     num_us = nu)
        #rep_data[[r]] <- df

        #rep_data <- bind_rows(rep_data)
        mir_data[[mir]] <- df
      }
      #mir_data <- dplyr::bind_rows(mir_data)
      sample_data[[s]] <- dplyr::bind_rows(mir_data)
    }
    sample_data <- dplyr::bind_rows(sample_data)
    sim_df <- list(nmir = nmir,
                   fn_s4U = fn_s4U,
                   p_new_real_tc = p_new_real_tc,
                   p_old_real_tc = p_old_real_tc,
                   nreads = nreads, # per miR per sample
                   nsamp = nsamp,
                   ctl = ctl,
                   mir_pnew_tc = mir_pnew_tc,
                   mir_pold_tc = mir_pold_tc,
                   sample_data = sample_data
    )
    return(sim_df)   # return a list of all the things we should know.
  }

  sim_df_1 <- simulateData()

  # end_1 <- Sys.time()


  # start_2 <- Sys.time()

  # Extract simulated cB and summarise data
  cB_sim_1 <- sim_df_1$sample_data

  rm(sim_df_1)

  cB_sim_1 <- cB_sim_1 %>% dplyr::group_by(S, TP, R, MIR, TC, MT, num_us) %>%
    dplyr::tally()

  # Define XF column
  cB_sim_1$XF <- cB_sim_1$MIR


  samp_list <- unique(cB_sim_1$S)

  type_list <- rep(0, times=length(samp_list))
  mut_list <- rep(0, times = length(samp_list))
  rep_list <- rep(0, times = length(samp_list))
  count <- 1
  for(i in samp_list){
    type_list[count] <- unique(cB_sim_1$TP[cB_sim_1$S == i])
    rep_list[count] <- unique(cB_sim_1$R[cB_sim_1$S == i])
    mut_list[count] <- unique(cB_sim_1$MT[cB_sim_1$S == i])
    count <- count + 1
  }


  colnames(cB_sim_1) <- c("sample", "TP", "R", "MIR", "TC", "MT", "nT", "n", "XF")

  metadf <- data.frame(tl = type_list*tl, Exp_ID = as.integer(mut_list))

  rownames(metadf) <- unique(cB_sim_1$sample)

  cB_sim_1$sample <- as.character(cB_sim_1$sample)

  bakRData <- bakR::bakRData(cB_sim_1, metadf)


  # end_2 <- Sys.time()


  # start_3 <- Sys.time()

  ## Create dataframe for Effect sizes and fn

  fn_vect <-  c()
  L2FC_kd_vect <- c()
  effect_vect <- c()
  fn_mean_vect <- c()


  for(j in 1:num_conds){
    if(j > 1){
      L2FC_kd_vect <- c(L2FC_kd_vect, L2FC_kd_mean[,j])
      effect_vect <- c(effect_vect, effect_mean[,j])
    }

    fn_mean_vect <- c(fn_mean_vect, logit(fn_mean) + effect_mean[,j])

    for(i in 1:ngene){
      fn_vect <- c(fn_vect, fn_real[i,j,])


    }
  }

  ## Make dataframes that are similar to Fit outputs

  Fn_rep_sim <- data.frame(Feature_ID = rep(rep(1:ngene, each = nreps), times = num_conds),
                           Replicate = rep(1:nreps, times = ngene*num_conds),
                           Exp_ID = rep(1:num_conds, each = ngene*nreps),
                           Logit_fn = logit(fn_vect),
                           fn = fn_vect)

  Fn_mean_sim <- data.frame(Feature_ID = rep(1:ngene, times = num_conds),
                            Exp_ID = rep(1:num_conds, each = ngene),
                            Avg_logit_fn = fn_mean_vect,
                            Avg_fn = inv_logit(fn_mean_vect))

  if(num_conds > 1){
    Effect_sim <- data.frame(Feature_ID = rep(1:ngene, times = (num_conds-1)),
                             Exp_ID = rep(2:num_conds, each = ngene),
                             L2FC_kdeg = L2FC_kd_vect,
                             effect = effect_vect)

    Effect_sim <- Effect_sim[order(Effect_sim$Feature_ID, Effect_sim$Exp_ID),]

    ## Order dataframes as they are in fit output
    Fn_rep_sim <- Fn_rep_sim[order(Fn_rep_sim$Feature_ID, Fn_rep_sim$Exp_ID, Fn_rep_sim$Replicate),]
    Fn_mean_sim <- Fn_mean_sim[order(Fn_mean_sim$Feature_ID, Fn_mean_sim$Exp_ID),]

    sim_data <- list(bakRData = bakRData,
                     sim_list = list(Effect_sim = Effect_sim,
                                     Fn_mean_sim = Fn_mean_sim,
                                     Fn_rep_sim = Fn_rep_sim,
                                     L2FC_ks_mean = L2FC_ks_mean,
                                     RNA_conc = RNA_conc*scale_factor,
                                     Counts = Counts) )

  }else{
    ## Order dataframes as they are in fit output
    Fn_rep_sim <- Fn_rep_sim[order(Fn_rep_sim$Feature_ID, Fn_rep_sim$Exp_ID, Fn_rep_sim$Replicate),]
    Fn_mean_sim <- Fn_mean_sim[order(Fn_mean_sim$Feature_ID, Fn_mean_sim$Exp_ID),]


    sim_data <- list(bakRData = bakRData,
                     sim_list = list(Fn_mean_sim = Fn_mean_sim,
                                     Fn_rep_sim = Fn_rep_sim,
                                     L2FC_ks_mean = L2FC_ks_mean,
                                     RNA_conc = RNA_conc*scale_factor,
                                     Counts = Counts) )
  }


  # end_3 <- Sys.time()
  #
  # time_3 <- round(end_3 - start_3, 3)
  # time_2 <- round(end_2 - start_2, 3)
  # time_1 <- round(end_1 - start_1, 3)
  #
  # sim_data$times <- c(time_1, time_2, time_3)

  return(sim_data)


}



#' Simulating nucleotide recoding data
#'
#' \code{Simulate_bakRData} simulates a `bakRData` object. It's output also includes the simulated
#' values of all kinetic parameters of interest. Only the number of genes (\code{ngene}) has to be set by the
#' user, but an extensive list of additional parameters can be adjusted.
#'
#' \code{Simulate_bakRData} simulates a `bakRData` object using a realistic generative model with many
#' adjustable parameters. Average RNA kinetic parameters are drawn from biologically inspired
#' distributions. Replicate variability is simulated by drawing a feature's
#' fraction new in a given replicate from a logit-Normal distribution with a heteroskedastic
#' variance term with average magnitude given by the chosen read count vs. variance relationship.
#' For each replicate, a feature's ksyn is drawn from a homoskedastic lognormal distribution. Read counts
#' can either be set to the same value for all simulated features or can be simulated according to
#' a heterodisperse negative binomial distribution. The latter is the default
#'
#' The number of Us in each sequencing read is drawn from a binomial distribution with number of trials
#' equal to the read length and probability of each nucleotide being a U drawn from a beta distribution. Each read is assigned to the
#' new or old population according to a Bernoulli distribution with p = fraction new. The number of
#' mutations in each read are then drawn from one of two binomial distributions; if the read is assigned to the
#' population of new RNA, the number of mutations are drawn from a binomial distribution with number of trials equal
#' to the number of Us and probability of mutation = \code{p_new}; if the read is assigned to the population of old RNA,
#' the number of mutations is instead drawn from a binomial distribution with the same number of trials but with the probability
#' of mutation = \code{p_old}. \code{p_new} must be greater than \code{p_old} because mutations in new RNA
#' arise from both background mutations that occur with probability \code{p_old} as well as metabolic label induced mutations
#'
#' Simulated read counts should be treated as if they are spike-in and RPKM normalized, so the same scale factor can be applied
#' to each sample when comparing the sequencing reads (e.g., if you are performing differential expression analysis).
#'
#' Function to simulate a bakRData object according to a realistic generative model
#' @param ngene Number of genes to simulate data for
#' @param num_conds Number of experimental conditions (including the reference condition) to simulate
#' @param nreps Number of replicates to simulate
#' @param eff_sd Effect size; more specifically, the standard deviation of the normal distribution from which non-zero
#' changes in logit(fraction new) are pulled from.
#' @param eff_mean Effect size mean; mean of normal distribution from which non-zero changes in logit(fraction new) are pulled from.
#' Note, setting this to 0 does not mean that some of the significant effect sizes will be 0, as any exact integer is impossible
#' to draw from a continuous random number generator. Setting this to 0 just means that there is symmetric stabilization and destabilization
#' @param fn_mean Mean of fraction news of simulated transcripts in reference condition. The logit(fraction) of RNA from each transcript that is
#' metabolically labeled (new) is drawn from a normal distribution with this mean
#' @param fn_sd Standard deviation of fraction news of simulated transcripts in reference condition. The logit(fraction) of RNA
#' from each transcript that is metabolically labeled (new) is drawn from a normal distribution with this sd
#' @param kslog_c Synthesis rate constants will be drawn from a lognormal distribution with meanlog = \code{kslog_c} - mean(log(kd_mean)) where kd_mean
#' is determined from the fraction new simulated for each gene as well as the label time (\code{tl}).
#' @param kslog_sd Synthesis rate lognormal standard deviation; see kslog_c documentation for details
#' @param tl metabolic label feed time
#' @param p_new metabolic label (e.g., s4U) induced mutation rate. Can be a vector of length num_conds
#' @param p_old background mutation rate
#' @param read_lengths Total read length for each sequencing read (e.g., PE100 reads correspond to read_lengths = 200)
#' @param p_do Rate at which metabolic label containing reads are lost due to dropout; must be between 0 and 1
#' @param noise_deg_a Slope of trend relating log10(standardized read counts) to log(replicate variability)
#' @param noise_deg_b Intercept of trend relating log10(standardized read counts) to log(replicate variability)
#' @param noise_synth Homoskedastic variability of L2FC(ksyn)
#' @param sd_rep Variance of lognormal distribution from which replicate variability is drawn
#' @param low_L2FC_ks Most negative L2FC(ksyn) that can be simulated
#' @param high_L2FC_ks Most positive L2FC(ksyn) that can be simulated
#' @param num_kd_DE Vector where each element represents the number of genes that show a significant change in stability relative
#' to the reference. 1st entry must be 0 by definition (since relative to the reference the reference sample is unchanged)
#' @param num_ks_DE Same as num_kd_DE but for significant changes in synthesis rates.
#' @param scale_factor Factor relating RNA concentration (in arbitrary units) to average number of read counts
#' @param sim_read_counts Logical; if TRUE, read counts are simulated as coming from a heterodisperse negative binomial distribution
#' @param a1 Heterodispersion 1/reads dependence parameter
#' @param a0 High read depth limit of negative binomial dispersion parameter
#' @param nreads Number of reads simulated if sim_read_counts is FALSE
#' @param alpha shape1 parameter of the beta distribution from which U-contents (probability that a nucleotide in a read from a transcript is a U) are
#' drawn for each gene.
#' @param beta shape2 parameter of the beta distribution from which U-contents (probability that a nucleotide in a read from a transcript is a U) are
#' drawn for each gene.
#' @param STL logical; if TRUE, simulation is of STL-seq rather than a standard TL-seq experiment. The two big changes are that a short read lenght is required
#' (< 60 nt) and that every read for a particular feature will have the same number of Us. Only one read length is simulated for simplicity.
#' @param STL_len Average length of simulated STL-seq length. Since Pol II typically pauses about 20-60 bases
#' from the promoter, this should be around 40
#' @importFrom magrittr %>%
#' @import data.table
#' @export
#' @return A list containing a simulated `bakRData` object as well as a list of simulated kinetic parameters of interest.
#' The contents of the latter list are:
#' \itemize{
#'  \item Effect_sim; Dataframe meant to mimic formatting of Effect_df that are part of \code{bakRFit(StanFit = TRUE)}, \code{bakRFit(HybridFit = TRUE)} and \code{bakRFit(bakRData object)} output.
#'  \item Fn_mean_sim; Dataframe meant to mimic formatting of Regularized_ests that is part of \code{bakRFit(bakRData object)} output. Contains information
#'  about the true fraction new simulated in each condition (the mean of the normal distribution from which replicate fraction news are simulated)
#'  \item Fn_rep_sim; Dataframe meant to mimic formatting of Fn_Estimates that is part of \\code{bakRFit(bakRData object)} output. Contains information
#'  about the fraction new simulated for each feature in each replicate of each condition.
#'  \item L2FC_ks_mean; The true L2FC(ksyn) for each feature in each experimental condition. The i-th column corresponds to the L2FC(ksyn) when comparing
#'  the i-th condition to the reference condition (defined as the 1st condition) so the 1st column is always all 0s
#'  \item RNA_conc; The average number of normalized read counts expected for each feature in each sample.
#' }
#'
Simulate_bakRData <- function(ngene, num_conds = 2L, nreps = 3L, eff_sd = 0.75, eff_mean = 0,
                         fn_mean = 0, fn_sd = 1, kslog_c = 0.8, kslog_sd = 0.95,
                         tl = 60, p_new = 0.05, p_old = 0.001, read_lengths = 200L,
                         p_do = 0, noise_deg_a = -0.3, noise_deg_b = -1.5, noise_synth = 0.1, sd_rep = 0.05,
                         low_L2FC_ks = -1, high_L2FC_ks = 1,
                         num_kd_DE = c(0L, as.integer(rep(round(as.integer(ngene)/2), times = as.integer(num_conds)-1))),
                         num_ks_DE = rep(0L, times = as.integer(num_conds)),
                         scale_factor = 150,
                         sim_read_counts = TRUE, a1 = 5, a0 = 0.01,
                         nreads = 50L, alpha = 25, beta = 75,
                         STL = FALSE, STL_len = 40){

  #### Parameter checks
  # STL
  if(!is.logical(STL)){
    "STL must be logical (TRUE or FALSE)"
  }

  # STL_len
  if(STL){
    if(STL_len < 25){
      stop("STL_len is less than 25. Must be between 25 and 55.")
    }else if(STL_len > 55){
      stop("STL_len is greater than 55. Must be between 25 and 55.")
    }
  }

  # alpha parameter of beta distribution
  if(!is.numeric(alpha)){
    stop("alpha must be numeric!")
  }else if(alpha <= 1){
    stop("alpha must be > 1")
  }else if(!is.numeric(beta)){
    stop("beta must be numeric!")
  }else if(beta <=1){
    stop("beta must be > 1")
  }else if(alpha/(alpha + beta) < 0.1){
    warning("alpha and beta are such that the average U-content is less than 0.1. That is unreasonably low and may make many simulated transcripts unanalyzable.")
  }else if(alpha/(alpha + beta) > 0.75){
    warning("alpha and beta are such that the average U-content is > 0.75. I wish U-contents were this high, your simulation may not reflect real data though.")
  }

  # Number of genes (features) to simulate
  if(!is.numeric(ngene)){
    stop("ngene must be numeric")
  }else if(!is.integer(ngene)){
    ngene <- as.integer(ngene)
  }

  if(ngene < 1){
    stop("ngene must be > 0; it represents the number of features to simulate")
  }

  # Number of experimental conditions


  if(!is.numeric(num_conds)){
    stop("num_conds must be an numeric")
  }else if(!is.integer(num_conds)){
    num_conds <- as.integer(num_conds)
  }

  if(num_conds < 1){
    stop("num_conds must be > 0; it represents the number of experimental conditions")
  }else if(num_conds == 1){
    warning("You are only simluating a reference condition with no experimental conditions.")
  }

  # Number of replicates
  if(!is.numeric(nreps)){
    stop("nreps must be numeric")
  }else if(!is.integer(nreps)){
    nreps <- as.integer(nreps)
  }

  if(nreps < 1){
    stop("nreps must be > 0; it represents the number of replicates")
  }else if(nreps == 1){
    warning("You are only simulating a single replicate; All statistical models implemented in bakR require > 1 replicate.")
  }

  # Effect size distribution standard deviation
  if(!is.numeric(eff_sd)){
    stop("eff_sd must be numeric")
  }else if(eff_sd <= 0){
    stop("eff_sd must be > 0; it will be the sd parameter of a call to rnorm when simulating effect sizes")
  }else if (eff_sd > 4){
    warning("You are simulating an unusually large eff_sd (> 4)")
  }


  # Fraction new distribution sd
  if(!is.numeric(fn_sd)){
    stop("fn_sd must be numeric")
  }else if(fn_sd <= 0){
    stop("fn_sd must be > 0; it will be the sd parameter of a call to rnorm when simulating reference fraction news")
  }else if (fn_sd > 2){
    warning("You are simulating an unusually large fn_sd (> 2). This will lead to lots of extreme fraction news (close to 0 and 1).")
  }

  # Label time
  if(!is.numeric(tl)){
    stop("tl must be numeric")
  }else if(tl <= 0){
    stop("tl must be > 0; it represents the label time in minutes.")
  }

  # s4U mutation rate
  if(!all(is.numeric(p_new))){
    stop("p_new must be numeric")
  }else if(!all(p_new > 0)){
    stop("p_new must be > 0; it represents the mutation rate of s4U labeled transcripts")
  }else if(!all(p_new <= 1)){
    stop("p_new must be <= 1; it represents the mutation rate of s4U labeled transcripts")
  }else if(!all(p_new < 0.2)){
    warning("You have simulated an unusually large s4U induced mutation rate (>= 0.2)")
  }else if(!all((p_new - p_old) > 0)){
    stop("All p_new must be > p_old, since the background mutation rate (p_old) should always be less than the labeled mutation rate (p_new)")
  }

  # Background mutation rate
  if(!all(is.numeric(p_old))){
    stop("p_old must be numeric")
  }else if(!all(p_old > 0)){
    stop("p_old must be > 0; it represents the background mutation rate")
  }else if(!all(p_old <= 1)){
    stop("p_old must be <= 1; it represents the background mutation rate")
  }else if(!all(p_old < 0.01)){
    warning("You have simulated an unusually large background mutation rate (>= 0.01)")
  }

  # Read length
  if(!all(is.numeric(read_lengths))){
    stop("read_lengths must be numeric")
  }else if(!all(is.integer(read_lengths))){
    read_lengths <- as.integer(read_lengths)
  }

  if(!all(read_lengths > 0)){
    stop("read_lengths must be > 0; it represents the total length of sequencing reads")
  }else if(!all(read_lengths >= 20)){
    warning("You have simulated an unusually short read length (< 20 nucleotides)")
  }

  # Dropout probability
  if(!all(is.numeric(p_do))){
    stop("p_do must be numeric")
  }else if(!all(p_do <= 1)){
    stop("p_do must be <= 1; it represents the percentage of s4U containing RNA lost during library prep")
  }else if(!all(p_do >= 0)){
    stop("p_do must be >= 0; it represents the percentage of s4U containing RNA lost during library prep")
  }else if(!all(p_do == 0)){
    warning("You are simulating dropout; statistical models implemented by bakR do not correct for dropout and thus will provide biased estimates")
  }

  # Heteroskedastic Slope
  if(!is.numeric(noise_deg_a)){
    stop("noise_deg_a must be numeric")
  }else if(noise_deg_a > 0){
    stop("noise_deg_a must be < 0; it represents the slope of the log10(read depth) vs. log(replicate variability) trend")
  }else if(noise_deg_a == 0){
    warning("You are simulating fraction new homoskedasticity, which is not reflective of actual nucleotide recoding data")
  }

  # Synthesis variability
  if(!is.numeric(noise_synth)){
    stop("noise_synth must be numeric")
  }else if(noise_synth < 0){
    stop("noise_synth must be >= 0; it represents the homoskeastic variability in the synthesis rate")
  }

  # Replicate variability variability
  if(!is.numeric(sd_rep)){
    stop("sd_rep must be numeric")
  }else if(sd_rep < 0){
    stop("sd_rep must be >= 0; it represents the sdlog parameter of the log-normal distribution from which replicate variabilites are simulated")
  }else if(sd_rep == 0){
    warning("You are simulating no variability in the replicate variability. This can cause minor convergence issues in the completely hierarchical
            models implemented by bakR.")
  }

  # L2FC(ksyn) Bounds
  if(!all(is.numeric(c(high_L2FC_ks, low_L2FC_ks)))){
    stop("high_L2FC_ks and low_L2FC_ks must be numeric")
  }else if(high_L2FC_ks < low_L2FC_ks){
    stop("high_L2FC_ks must be greater than low_L2FC_ks; they represent the upper and lower bound of a uniform distribution respectively")
  }

  # Number of differentially stabilized features
  if(!all(is.numeric(num_kd_DE))){
    stop("All elements of num_kd_DE must be numeric")
  }else if(!all(is.integer(num_kd_DE))){
    num_kd_DE <- as.integer(num_kd_DE)
  }

  if (length(num_kd_DE) > num_conds){
    stop("num_kd_DE has too many elements; it should be a vector of length == num_conds")
  }else if (length(num_kd_DE) < num_conds){
    stop("num_kd_DE has too few elements; it should be a vector of length == num_conds")
  }else if(num_kd_DE[1] != 0){
    stop("The 1st element of num_kd_DE must equal 0 by definition, as it represents the number of features differentially stabilized in the reference
         condition relative to the reference condition")
  }else if(!all(num_kd_DE <= ngene)){
    stop("Not all elements of num_kd_DE are less than or equal to the total number of simulated features")
  }else if(!all(num_kd_DE >= 0)){
    stop("Not all elements of num_kd_DE are > 0")
  }

  # Number of differentially synthesized features
  if(!all(is.numeric(num_ks_DE))){
    stop("All elements of num_ks_DE must be numeric")
  }else if(!all(is.integer(num_ks_DE))){
    num_ks_DE <- as.integer(num_ks_DE)
  }

  if(length(num_ks_DE) > num_conds){
    stop("num_ks_DE has too many elements; it should be a vector of length == num_conds")
  }else if (length(num_ks_DE) < num_conds){
    stop("num_ks_DE has too few elements; it should be a vector of length == num_conds")
  }else if(num_ks_DE[1] != 0){
    stop("The 1st element of num_ks_DE must equal 0 by definition, as it represents the number of features differentially stabilized in the reference
         condition relative to the reference condition")
  }else if(!all(num_ks_DE <= ngene)){
    stop("Not all elements of num_ks_DE are less than or equal to the total number of simulated features")
  }else if(!all(num_ks_DE >= 0)){
    stop("Not all elements of num_ks_DE are > 0")
  }

  # Scale factor
  if(!is.numeric(scale_factor)){
    stop("scale_factor must be numeric")
  }else if(scale_factor <=0){
    stop("scale_factor must be > 0; it represents the factor by which the RNA concentration is multiplied to yield the average number of sequencing reads")
  }else if(scale_factor < 10){
    warning("You are simulating an unusually low scale_factor (< 10); small scale factors will lead to low read counts.")
  }

  # Simulate read count Boolean
  if(!is.logical(sim_read_counts)){
    stop("sim_read_counts must be a logical (TRUE or FALSE)")
  }else if(sim_read_counts){
    if(!is.integer(nreads)){
      stop("nreads must be an integer")
    }else if(nreads <= 0){
      stop("nreads must be > 0; it represents the number of sequencing reads to be simulated for all features")
    }
  }

  # Heterodispersion parameters

  if(!all(is.numeric(c(a1, a0)))){
    stop("a1 and a0 must be numeric")
  }else if(a1 < 0){
    stop("a1 must be >= 0; it relates the average read count to the negative binomial dispersion parameter")
  }else if(a0 <= 0){
    stop("a0 must be > 0; it represents the high read count limit of the negative binomial dispersion parameter")
  }

  # Parameters without numeric bounds
  if(!is.numeric(noise_deg_b)){
    stop("noise_deg_b must be numeric")
  }

  if(!is.numeric(eff_mean)){
    stop("eff_mean must be numeric")
  }



  if(length(p_new) == 1){
    p_new <- rep(p_new, times=num_conds)
  }else if(length(p_new) != num_conds){
    stop("p_new must be of length 1 or length == num_conds")
  }

  if(length(p_old) == 1){
    p_old <- rep(p_old, times=num_conds)
  }else if(length(p_old) != num_conds){
    stop("p_old must be of length 1 or length == num_conds")
  }

  if(length(read_lengths) == 1){
    read_lengths <- rep(read_lengths, times=num_conds)
  }else if(length(read_lengths) != num_conds){
    stop("read_lengths must be of length 1 or length == num_conds")
  }

  if(length(p_do) == 1){
    p_do <- rep(p_do, times=num_conds)
  }else if(length(p_do) != num_conds){
    stop("p_do must be of length 1 or length == num_conds")
  }


  #### Simulation function from bakR

  # Average number of Us in reads from each simulated feature
  U_cont <- stats::rbeta(ngene, alpha, beta)

  # Define helper functions:
  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))

  #Initialize matrices
  fn <- rep(0, times=ngene*nreps*num_conds)
  dim(fn) <- c(ngene, num_conds, nreps)

  Counts <- fn


  kd <- fn
  ks <- fn


  #Initialize vectors of mean values for each gene and condition
  fn_mean <- inv_logit(stats::rnorm(n=ngene, mean=fn_mean, sd=fn_sd))
  kd_mean <- -log(1-fn_mean)/tl
  ks_mean <- stats::rlnorm(n=ngene, meanlog = kslog_c + log(kd_mean), sdlog = kslog_sd)


  effect_mean <- rep(0, times = ngene*num_conds)
  dim(effect_mean) <- c(ngene, num_conds)
  L2FC_ks_mean <- effect_mean
  L2FC_kd_mean <- effect_mean

  # Simualte true kinetic parameter values
  for(i in 1:ngene){

    #Make sure the user didn't input the wrong
    #number of significant genes
    if (i == 1 ){
      if (length(num_kd_DE) > num_conds){
        stop("num_kd_DE has too many elements")
      } else if (length(num_kd_DE) < num_conds){
        stop("num_kd_DE has too few elements")
      }

      if (length(num_ks_DE) > num_conds){
        stop("num_ks_DE has too many elements")
      } else if (length(num_ks_DE) < num_conds){
        stop("num_ks_DE has too few elements")
      }
    }

    # Simulate effect sizes (differences in kinetic parameters)
      # effect_mean = difference in logit(fraction new)s
      # L2FC_ks_mean = difference in L2FC(ksyn)
    for(j in 1:num_conds){
      if(j == 1){
        effect_mean[i,1] <- 0
        L2FC_ks_mean[i,1] <- 0
      }else{
        if(i < (ngene-num_kd_DE[j] + 1)){
          effect_mean[i,j] <- 0
        }else{
          effect_mean[i,j] <- stats::rnorm(n=1, mean=eff_mean, sd=eff_sd)
        }
        if(i < (ngene-num_ks_DE[j] + 1)){
          L2FC_ks_mean[i,j] <- 0
        }else{
          if (stats::runif(1) < 0.5){
            L2FC_ks_mean[i,j] <- stats::runif(n=1, min=low_L2FC_ks, max=high_L2FC_ks)
          }else{
            L2FC_ks_mean[i,j] <- stats::runif(n=1, min=-high_L2FC_ks, max=-low_L2FC_ks)
          }

        }
      }
    }
  }

  # Difference in L2FC(kdeg)
  L2FC_kd_mean <- log2(log(1 - inv_logit(fn_mean + effect_mean))/log(1- inv_logit(fn_mean)))


  #Simulate read counts
  if (sim_read_counts == TRUE){
    L2FC_tot_mean <- L2FC_ks_mean - L2FC_kd_mean
    RNA_conc <- (ks_mean*2^(L2FC_ks_mean))/(kd_mean*2^(L2FC_kd_mean))
    a1 <- 5
    a0 <- 0.01

    for(i in 1:ngene){
      for(j in 1:num_conds){
        for(k in 1:nreps){
          # DESeq2 model of heterodisperse read counts
          Counts[i, j, k] <- stats::rnbinom(n=1, size=1/((a1/(scale_factor*RNA_conc[i,j])) + a0), mu = scale_factor*RNA_conc[i,j])

          if(Counts[i, j, k] < 5){
            Counts[i, j, k] <- Counts[i, j, k] + stats::rpois(n=1, lambda = 2) + 1
          }
        }
      }
    }
  } else{
    if(sim_from_data == FALSE){
      Counts <- rep(nreads, times= ngene*num_conds*nreps)
      dim(Counts) <- c(ngene, num_conds, nreps)
    }

  }

  standard_RNA <- matrix(0, nrow = ngene, ncol = num_conds)
  mean_RNA <- rep(0, times = num_conds)
  sd_RNA <- rep(0, times = num_conds)

  # Simulate kdeg and ksyn for each feature and replicate.
    # replicate variability is simulated here
  for(j in 1:num_conds){
    mean_RNA[j] <- mean(log10(RNA_conc[,j]*scale_factor))
    sd_RNA[j] <- stats::sd(log10(RNA_conc[,j]*scale_factor))
    for(i in 1:ngene){
      standard_RNA[i,j] <- (log10(RNA_conc[i,j]*scale_factor) - mean_RNA[j])/sd_RNA[j]
      for(k in 1:nreps){
        fn[i, j, k] <- inv_logit(stats::rnorm(1, mean=(logit(fn_mean[i]) + effect_mean[i,j]), sd = stats::rlnorm(1, noise_deg_a*standard_RNA[i,j] + noise_deg_b, sd_rep )))
        ks[i,j,k] <- exp(stats::rnorm(1, mean=log((2^L2FC_ks_mean[i,j])*ks_mean[i]), sd=noise_synth))
      }
    }
  }

  # kdeg
  kd <- -log(1 - fn)/tl



  l <- ngene

  p_do <- matrix(rep(0, times = num_conds*nreps), nrow = nreps, ncol = num_conds)
  fn_real <- fn

  # Simulate dropout
  for(j in 1:num_conds){
    for(k in 1:nreps){
      fn_real[,j,k] <- (fn[,j,k]*(1-p_do[k,j]))/(1 - p_do[k,j]*(1 - fn[,j,k]))
      Counts[,j,k] <- Counts[,j,k] - Counts[,j,k]*fn[,j,k]*p_do[k,j]
    }
  }



  ## Parameter names from old function

  nmir <- l
  fn_s4U <- fn_real
  p_new_real_tc <- p_new
  p_old_real_tc <- p_old
  read_length <- read_lengths
  nreads <- Counts
  nsamp <- (nreps*num_conds) + num_conds
  ctl <- c(rep(1, times=nreps*num_conds), rep(0, times=num_conds))
  mt <- c(rep(1:num_conds, each=nreps),seq(from=1,to=num_conds,by=1))
  replicate <- c(rep(seq(from=1, to=nreps), times=num_conds), rep(1, times=num_conds))

  ### Simulate for each sequencing read, whether it is new or old

  # Read counts
  Reads_vect <- as.vector(Counts)

  # True average fraction news
  fn_vect <- as.vector(fn_real)

  # Total number of reads
  tot_reads <- sum(Reads_vect)

  # Numerical IDs for feature, experimental condition, replicate, and sample
  Gene_ID <- rep(1:l, times = num_conds*nreps)
  Exp_ID <- rep(rep(1:num_conds, each = l), times = nreps)
  Rep_ID <- rep(1:nreps, each = l*num_conds)
  Samp_ID <- rep(rep(seq(from = 1, to = num_conds*nreps, by = nreps), times = nreps) + rep(0:(nreps-1), each = num_conds), each = l)

  # Avg. number of Us in sequencing reads from each feature
  U_contents <- U_cont[Gene_ID]


  # For simulating STL-seq, assume a single pause site. This means that the U-content
  # is identical for all reads of a given length. To easily simulate this, use U-content
  # to calculate average number of Us, and add a bit of Poisson noise to the number of Us
  # to simualte variation about the average due to the polymerase pausing at slightly
  # different exact positions (i.e., it may go bit past or end up a bit short of the average
  # pause location, and thus incorporate a couple extra or a couple fewer Us.
  if(STL){
    nU <- abs(rep(round(U_contents*STL_len), times = Reads_vect) + sign(unlist(purrr::map(Reads_vect, stats::runif, min = -0.1, max = 0.1)))*unlist(purrr::map(Reads_vect, stats::rpois, lambda = 0.5)))
  }else{
    nU <- unlist(purrr::pmap(list(n = Reads_vect, size = read_length[Exp_ID], prob = U_contents),
                             stats::rbinom))
  }

  # 1 = new read; 0 = old read
  newness <- unlist(purrr::pmap(list(n = Reads_vect, p = fn_vect), purrr::rbernoulli))


  Gene_ID <- rep(Gene_ID, times = Reads_vect)
  Exp_ID <- rep(Exp_ID, times = Reads_vect)
  Rep_ID <- rep(Rep_ID, times = Reads_vect)
  Samp_ID <- rep(Samp_ID, times = Reads_vect)

  # Simulate number of mutations
  nmut <- stats::rbinom(n = length(nU), size = nU, prob = newness*p_new_real_tc[Exp_ID] + (1-newness)*p_old_real_tc[Exp_ID])



  # Create simualted cB file
  cB <- dplyr::tibble(S = Samp_ID,
               R = Rep_ID,
               MT = Exp_ID,
               FN = Gene_ID,
               TC = nmut,
               nT = nU,
               TP = rep(1, times = length(Samp_ID)))


  ### Simulate -s4U control data

  Gene_ID <- rep(1:l, times = num_conds*nreps)
  Exp_ID <- rep(rep(1:num_conds, each = l), times = nreps)
  Rep_ID <- rep(1:nreps, each = l*num_conds)
  Samp_ID <- rep(rep(seq(from = 1, to = num_conds*nreps, by = nreps), times = nreps) + rep(0:(nreps-1), each = num_conds), each = l)
  U_contents <- U_cont[Gene_ID]

  ctl_data <- which(Rep_ID == 1)

  Reads_ctl <- Reads_vect[ctl_data]
  Gene_ctl <- Gene_ID[ctl_data]
  Exp_ctl <- Exp_ID[ctl_data]
  Rep_ctl <- Rep_ID[ctl_data]
  Samp_ctl <- rep((max(Samp_ID) + 1):(max(Samp_ID) + num_conds), each = l )
  U_ctl <- U_contents[ctl_data]

  nU <- unlist(purrr::pmap(list(n = Reads_ctl, size = read_length[Exp_ctl], prob = U_ctl),
                    stats::rbinom))


  Gene_ctl <- rep(Gene_ctl, times = Reads_ctl)
  Exp_ctl <- rep(Exp_ctl, times = Reads_ctl)
  Rep_ctl <- rep(Rep_ctl, times = Reads_ctl)
  Samp_ctl <- rep(Samp_ctl, times = Reads_ctl)

  nmut <- stats::rbinom(n = length(nU), size = nU, prob = p_old_real_tc[Exp_ctl])

  cB_ctl <- data.frame(S = Samp_ctl,
                   R = Rep_ctl,
                   MT = Exp_ctl,
                   FN = Gene_ctl,
                   TC = nmut,
                   nT = nU,
                   TP = rep(0, times = length(Samp_ctl)))

  cB_final <- dplyr::bind_rows(list(cB, cB_ctl))


  # Extract simulated cB and summarise data
  cB_sim_1 <- data.table::setDT(cB_final)

  cols <- colnames(cB_sim_1)
  cB_sim_1 <- cB_sim_1[,.N, by = .(S, R, MT, FN, TC, nT, TP)]
  colnames(cB_sim_1) <- c(cols, "n")

  cB_sim_1 <- cB_sim_1[order(cB_sim_1$S),]


  # Define XF column
  cB_sim_1$XF <- cB_sim_1$FN


  # Map sample list to experimental characteristics
    # type = whether s4U labeled (1) or not (0)
    # mut = experimental condition ID
    # rep = replicate ID
  samp_list <- unique(cB_sim_1$S)

  type_list <- rep(0, times=length(samp_list))
  mut_list <- rep(0, times = length(samp_list))
  rep_list <- rep(0, times = length(samp_list))
  count <- 1
  for(i in samp_list){
    type_list[count] <- unique(cB_sim_1$TP[cB_sim_1$S == i])
    rep_list[count] <- unique(cB_sim_1$R[cB_sim_1$S == i])
    mut_list[count] <- unique(cB_sim_1$MT[cB_sim_1$S == i])
    count <- count + 1
  }


  colnames(cB_sim_1) <- c("sample", "R", "MT", "FN", "TC", "nT", "TP", "n", "XF")

  # metadata data frame for bakR
  metadf <- data.frame(tl = type_list*tl, Exp_ID = as.integer(mut_list))

  rownames(metadf) <- unique(cB_sim_1$sample)

  cB_sim_1$sample <- as.character(cB_sim_1$sample)

  # Make bakRData object
  bakRData <- bakR::bakRData(cB_sim_1, metadf)


  ## Create data frames storing simulated truths

  fn_vect <-  c()
  L2FC_kd_vect <- c()
  effect_vect <- c()
  fn_mean_vect <- c()


  for(j in 1:num_conds){
    if(j > 1){
      L2FC_kd_vect <- c(L2FC_kd_vect, L2FC_kd_mean[,j])
      effect_vect <- c(effect_vect, effect_mean[,j])
    }

    fn_mean_vect <- c(fn_mean_vect, logit(fn_mean) + effect_mean[,j])

    for(i in 1:ngene){
      fn_vect <- c(fn_vect, fn_real[i,j,])


    }
  }

  # Make data frames similar to bakRFit outputs in terms of ordering

  Fn_rep_sim <- data.frame(Feature_ID = rep(rep(1:ngene, each = nreps), times = num_conds),
                           Replicate = rep(1:nreps, times = ngene*num_conds),
                           Exp_ID = rep(1:num_conds, each = ngene*nreps),
                           Logit_fn = logit(fn_vect),
                           fn = fn_vect)

  Fn_mean_sim <- data.frame(Feature_ID = rep(1:ngene, times = num_conds),
                            Exp_ID = rep(1:num_conds, each = ngene),
                            Avg_logit_fn = fn_mean_vect,
                            Avg_fn = inv_logit(fn_mean_vect))

  if(num_conds > 1){
    Effect_sim <- data.frame(Feature_ID = rep(1:ngene, times = (num_conds-1)),
                             Exp_ID = rep(2:num_conds, each = ngene),
                             L2FC_kdeg = L2FC_kd_vect,
                             effect = effect_vect)

    Effect_sim <- Effect_sim[order(Effect_sim$Feature_ID, Effect_sim$Exp_ID),]

    ## Order dataframes as they are in fit output
    Fn_rep_sim <- Fn_rep_sim[order(Fn_rep_sim$Feature_ID, Fn_rep_sim$Exp_ID, Fn_rep_sim$Replicate),]
    Fn_mean_sim <- Fn_mean_sim[order(Fn_mean_sim$Feature_ID, Fn_mean_sim$Exp_ID),]

    sim_data <- list(bakRData = bakRData,
                     sim_list = list(Effect_sim = Effect_sim,
                                     Fn_mean_sim = Fn_mean_sim,
                                     Fn_rep_sim = Fn_rep_sim,
                                     L2FC_ks_mean = L2FC_ks_mean,
                                     RNA_conc = RNA_conc*scale_factor,
                                     Counts = Counts) )

  }else{
    ## Order dataframes as they are in fit output
    Fn_rep_sim <- Fn_rep_sim[order(Fn_rep_sim$Feature_ID, Fn_rep_sim$Exp_ID, Fn_rep_sim$Replicate),]
    Fn_mean_sim <- Fn_mean_sim[order(Fn_mean_sim$Feature_ID, Fn_mean_sim$Exp_ID),]


    sim_data <- list(bakRData = bakRData,
                     sim_list = list(Fn_mean_sim = Fn_mean_sim,
                                     Fn_rep_sim = Fn_rep_sim,
                                     L2FC_ks_mean = L2FC_ks_mean,
                                     RNA_conc = RNA_conc*scale_factor,
                                     Counts = Counts) )
  }

  return(sim_data)

}
