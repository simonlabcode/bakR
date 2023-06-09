#' Efficiently analyze nucleotide recoding data
#'
#' \code{fast_analysis} analyzes nucleotide recoding data maximum likelihood estimation with the L-BFGS-B algorithm
#' implemented by \code{stats::optim} combined with analytic solutions to simple Bayesian models to perform
#' approximate partial pooling. Output includes kinetic parameter estimates in each replicate, kinetic parameter estimates
#' averaged across replicates, and log-2 fold changes in the degradation rate constant (L2FC(kdeg)).
#' Averaging takes into account uncertainties estimated using the Fisher Information and estimates
#' are regularized using analytic solutions of fully Bayesian models. The result is that kdegs are
#' shrunk towards population means and that uncertainties are shrunk towards a mean-variance trend estimated as part of the analysis.
#'
#' Unless the user supplies estimates for pnew and pold, the first step of \code{fast_analysis} is to estimate the background
#' and metabolic label (will refer to as s4U for simplicity, though bakR is compatible with other metabolic labels such as s6G)
#' induced mutation rates. The former is best performed with a -s4U control sample, that is, a normal RNA-seq sample
#' that lacks a -s4U feed or TimeLapse chemistry conversion of s4U to a C analog. If this sample is missing, both background and
#' s4U induced mutation rates are estimated from the s4U fed samples. For the s4U mutation rate, features with sufficient read depth,
#' as defined by the \code{read_cut} parameter, and the highest mutation rates are assumed to be completely labeled. Thus, the
#' average mutation rates in these features is taken as the estimate of the s4U induced mutation rate in that sample. s4U induced mutation
#' rates are estimated on a per-sample basis as there is often much more variability in these mutation rates than in the background
#' mutation rates.
#'
#' If a -s4U control is included, the background mutation rate is estimated using all features in the control sample(s) with read depths
#' greater than \code{read_cut}. The average mutation rate among these features is taken as the estimated background mutation rate,
#' and that background is assumed to be constant for all samples. If a -s4U control is missing, then a strategy similar to that used
#' to estimate s4U induced mutation rates is used. In this case, the lowest mutation rate features with sufficient read depths are used,
#' and there average mutation rate is the background mutation rate estimate, as these features are assumed to be almost entirely unlabeled.
#' Another slightly more computationally intensive but more accurate strategy to estimate mutation rates is to set \code{StanRate} = TRUE.
#' This will fit a non-hierarchical mixture model to a small subset of transcripts using 'Stan'. The default in \code{bakRFit} is to use
#' 25 transcripts. If \code{StanRate} is TRUE, then a data list must be passed to \code{Stan_data} of the form that appears in the
#' bakRFit object's Data_list$Stan_data entry.
#'
#' Once mutation rates are estimated, fraction news for each feature in each sample are estimated. The approach utilized is MLE
#' using the L-BFGS-B algorithm implemented in \code{stats::optim}. The assumed likelihood function is derived from a Poisson mixture
#' model with rates adjusted according to each feature's empirical U-content (the average number of Us present in sequencing reads mapping
#' to that feature in a particular sample). Fraction new estimates are then converted to degradation rate constant estimates using
#' a solution to a simple ordinary differential equation model of RNA metabolism.
#'
#' Once fraction new and kdegs are estimated, the uncertainty in these parameters is estimated using the Fisher Information. In the limit of
#' large datasets, the variance of the MLE is inversely proportional to the Fisher Information evaluated at the MLE. Mixture models are
#' typically singular, meaning that the Fisher information matrix is not positive definite and asymptotic results for the variance
#' do not necessarily hold. As the mutation rates are estimated a priori and fixed to be > 0, these problems are eliminated. In addition, when assessing
#' the uncertainty of replicate fraction new estimates, the size of the dataset is the raw number of sequencing reads that map to a
#' particular feature. This number is often large (>100) which increases the validity of invoking asymptotics.
#'
#' With kdegs and their uncertainties estimated, replicate estimates are pooled and regularized. There are two key steps in this
#' downstream analysis. 1st, the uncertainty for each feature is used to fit a linear ln(uncertainty) vs. log10(read depth) trend,
#' and uncertainties for individual features are shrunk towards the regression line. The uncertainty for each feature is a combination of the
#' Fisher Information asymptotic uncertainty as well as the amount of variability seen between estimates. Regularization of uncertainty
#' estimates is performed using the analytic results of a Normal distribution likelihood with known mean and unknown variance and conjugate
#' priors. The prior parameters are estimated from the regression and amount of variability about the regression line. The strength of
#' regularization can be tuned by adjusting the \code{prior_weight} parameter, with larger numbers yielding stronger shrinkage towards
#' the regression line. The 2nd step is to regularize the average kdeg estimates. This is done using the analytic results of a
#' Normal distribution likelihood model with unknown mean and known variance and conjugate priors. The prior parameters are estimated from the
#' population wide kdeg distribution (using its mean and standard deviation as the mean and standard deviation of the normal prior).
#' In the 1st step, the known mean is assumed to be the average kdeg, averaged across replicates and weighted by the number of reads
#' mapping to the feature in each replicate. In the 2nd step, the known variance is assumed to be that obtained following regularization
#' of the uncertainty estimates.
#'
#' Effect sizes (changes in kdeg) are obtained as the difference in log(kdeg) means between the reference and experimental
#' sample(s), and the log(kdeg)s are assumed to be independent so that the variance of the effect size is the sum of the
#' log(kdeg) variances. P-values assessing the significance of the effect size are obtained using a moderated t-test with number
#' of degrees of freedom determined from the uncertainty regression hyperparameters and are adjusted for multiple testing using the Benjamini-
#' Hochberg procedure to control false discovery rates (FDRs).
#'
#' In some cases, the assumed ODE model of RNA metabolism will not accurately model the dynamics of a biological system being analyzed.
#' In these cases, it is best to compare logit(fraction new)s directly rather than converting fraction new to log(kdeg).
#' This analysis strategy is implemented when \code{NSS} is set to TRUE. Comparing logit(fraction new) is only valid
#' If a single metabolic label time has been used for all samples. For example, if a label time of 1 hour was used for NR-seq
#' data from WT cells and a 2 hour label time was used in KO cells, this comparison is no longer valid as differences in
#' logit(fraction new) could stem from differences in kinetics or label times.
#'
#'
#' @param df Dataframe in form provided by cB_to_Fast
#' @param pnew Labeled read mutation rate; default of 0 means that model estimates rate from s4U fed data. If pnew is provided by user, must be  a vector
#' of length == number of s4U fed samples. The 1st element corresponds to the s4U induced mutation rate estimate for the 1st replicate of the 1st
#' experimental condition; the 2nd element corresponds to the s4U induced mutation rate estimate for the 2nd replicate of the 1st experimental condition,
#' etc.
#' @param pold Unlabeled read mutation rate; default of 0 means that model estimates rate from no-s4U fed data
#' @param no_ctl Logical; if TRUE, then -s4U control is not used for background mutation rate estimation
#' @param read_cut Minimum number of reads for a given feature-sample combo to be used for mut rate estimates
#' @param features_cut Number of features to estimate sample specific mutation rate with
#' @param nbin Number of bins for mean-variance relationship estimation. If NULL, max of 10 or (number of logit(fn) estimates)/100 is used
#' @param prior_weight Determines extent to which logit(fn) variance is regularized to the mean-variance regression line
#' @param MLE Logical; if TRUE then replicate logit(fn) is estimated using maximum likelihood; if FALSE more conservative Bayesian hypothesis testing is used
#' @param ztest TRUE; if TRUE, then a z-test is used for p-value calculation rather than the more conservative moderated t-test. 
#' @param lower Lower bound for MLE with L-BFGS-B algorithm
#' @param upper Upper bound for MLE with L-BFGS-B algorithm
#' @param se_max Uncertainty given to those transcripts with estimates at the upper or lower bound sets. This prevents downstream errors due to
#' abnormally high standard errors due to transcripts with extreme kinetics
#' @param mut_reg If MLE has instabilities, empirical mut rate will be used to estimate fn, multiplying pnew by 1+mut_reg and pold by 1-mut_reg to regularize fn
#' @param p_mean Mean of normal distribution used as prior penalty in MLE of logit(fn)
#' @param p_sd Standard deviation of normal distribution used as prior penalty in MLE of logit(fn)
#' @param StanRate Logical; if TRUE, a simple 'Stan' model is used to estimate mutation rates for fast_analysis; this may add a couple minutes
#' to the runtime of the analysis.
#' @param Stan_data List; if StanRate is TRUE, then this is the data passed to the 'Stan' model to estimate mutation rates. If using the \code{bakRFit}
#' wrapper of \code{fast_analysis}, then this is created automatically.
#' @param null_cutoff bakR will test the null hypothesis of |effect size| < |null_cutoff|
#' @param NSS Logical; if TRUE, logit(fn)s are compared rather than log(kdeg) so as to avoid steady-state assumption.
#' @param Chase Logical; Set to TRUE if analyzing a pulse-chase experiment. If TRUE, kdeg = -ln(fn)/tl where fn is the fraction of
#' reads that are s4U (more properly referred to as the fraction old in the context of a pulse-chase experiment)
#' @param BDA_model Logical; if TRUE, variance is regularized with scaled inverse chi-squared model. Otherwise a log-normal
#' model is used.
#' @param multi_pold Logical; if TRUE, pold is estimated for each sample rather than use a global pold estimate.
#' @param Long Logical; if TRUE, long read optimized fraction new estimation strategy is used.
#' @param kmeans Logical; if TRUE, kmeans clustering on read-specific mutation rates is used to estimate pnews and pold.
#' @return List with dataframes providing information about replicate-specific and pooled analysis results. The output includes:
#' \itemize{
#'  \item Fn_Estimates; dataframe with estimates for the fraction new and fraction new uncertainty for each feature in each replicate.
#'  The columns of this dataframe are:
#'  \itemize{
#'   \item Feature_ID; Numerical ID of feature
#'   \item Exp_ID; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item Replicate; Numerical ID for replicate
#'   \item logit_fn; logit(fraction new) estimate, unregularized
#'   \item logit_fn_se; logit(fraction new) uncertainty, unregularized and obtained from Fisher Information
#'   \item nreads; Number of reads mapping to the feature in the sample for which the estimates were obtained
#'   \item log_kdeg; log of degradation rate constant (kdeg) estimate, unregularized
#'   \item kdeg; degradation rate constant (kdeg) estimate
#'   \item log_kd_se; log(kdeg) uncertainty, unregularized and obtained from Fisher Information
#'   \item sample; Sample name
#'   \item XF; Original feature name
#'  }
#'  \item Regularized_ests; dataframe with average fraction new and kdeg estimates, averaged across the replicates and regularized
#'  using priors informed by the entire dataset. The columns of this dataframe are:
#'  \itemize{
#'   \item Feature_ID; Numerical ID of feature
#'   \item Exp_ID; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item avg_log_kdeg; Weighted average of log(kdeg) from each replicate, weighted by sample and feature-specific read depth
#'   \item sd_log_kdeg; Standard deviation of the log(kdeg) estimates
#'   \item nreads; Total number of reads mapping to the feature in that condition
#'   \item sdp; Prior standard deviation for fraction new estimate regularization
#'   \item theta_o; Prior mean for fraction new estimate regularization
#'   \item sd_post; Posterior uncertainty
#'   \item log_kdeg_post; Posterior mean for log(kdeg) estimate
#'   \item kdeg; exp(log_kdeg_post)
#'   \item kdeg_sd; kdeg uncertainty
#'   \item XF; Original feature name
#'  }
#'  \item Effects_df; dataframe with estimates of the effect size (change in logit(fn)) comparing each experimental condition to the
#'  reference sample for each feature. This dataframe also includes p-values obtained from a moderated t-test. The columns of this
#'  dataframe are:
#'  \itemize{
#'   \item Feature_ID; Numerical ID of feature
#'   \item Exp_ID; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item L2FC(kdeg); Log2 fold change (L2FC) kdeg estimate or change in logit(fn) if NSS TRUE
#'   \item effect; LFC(kdeg)
#'   \item se; Uncertainty in L2FC_kdeg
#'   \item pval; P-value obtained using effect_size, se, and a z-test
#'   \item padj; pval adjusted for multiple testing using Benjamini-Hochberg procedure
#'   \item XF; Original feature name
#'  }
#'  \item Mut_rates; list of two elements. The 1st element is a dataframe of s4U induced mutation rate estimates, where the mut column
#'  represents the experimental ID and the rep column represents the replicate ID. The 2nd element is the single background mutation
#'  rate estimate used
#'  \item Hyper_Parameters; vector of two elements, named a and b. These are the hyperparameters estimated from the uncertainties for each
#'  feature, and represent the two parameters of a Scaled Inverse Chi-Square distribution. Importantly, a is the number of additional
#'  degrees of freedom provided by the sharing of uncertainty information across the dataset, to be used in the moderated t-test.
#'  \item Mean_Variance_lms; linear model objects obtained from the uncertainty vs. read count regression model. One model is run for each Exp_ID
#' }
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#'
#' # Simulate small dataset
#' sim <- Simulate_bakRData(300, nreps = 2)
#'
#' # Fit fast model to get fast_df
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Fit fast model with fast_analysis
#' Fast_Fit <- fast_analysis(Fit$Data_lists$Fast_df)
#' }
#' @export
fast_analysis <- function(df, pnew = NULL, pold = NULL, no_ctl = FALSE,
                          read_cut = 50, features_cut = 50,
                          nbin = NULL, prior_weight = 2,
                          MLE = TRUE, ztest = FALSE,
                          lower = -7, upper = 7,
                          se_max = 2.5,
                          mut_reg = 0.1,
                          p_mean = 0,
                          p_sd = 1,
                          StanRate = FALSE,
                          Stan_data = NULL,
                          null_cutoff = 0,
                          NSS = FALSE,
                          Chase = FALSE,
                          BDA_model = FALSE,
                          multi_pold = FALSE,
                          Long = FALSE,
                          kmeans = FALSE){


  # Bind variables locally to resolve devtools::check() Notes
  nT <- n <- TC <- mut <- reps <- fnum <- New_prob <- Old_prob <- NULL
  News <- Fn_rep_est <- lam_n <- lam_o <- logit_fn_rep <- totTC <- NULL
  totU <- tot_mut <- totUs <- avg_mut <- logit_fn_rep.x <- logit_fn_rep.y <- NULL
  kd_rep_est <- U_cont <- Exp_l_fn <- Fisher_lkd_num <- Fisher_lkd_den <- NULL
  Inv_Fisher_Logit_1 <- Inv_Fisher_Logit_2 <- Inv_Fisher_Logit_3 <- NULL
  tot_n <- Fisher_kdeg <- Fisher_Logit <- Logit_fn_se <- bin_ID <- kd_sd_log <- NULL
  intercept <- slope <- log_kd_rep_est <- avg_log_kd <- sd_log_kd <- NULL
  sd_post <- sdp <- theta_o <- log_kd_post <- effect_size <- effect_std_error <- NULL
  n_new <- nreads <- logit_fn_se <- log_kd_se <- NULL
  


  # Check input validity -------------------------------------------------------

  ## Check df
  if (!all(c("sample", "XF", "TC", "nT", "n", "fnum", "type", "mut", "reps", "tl") %in% names(df))) {
    stop("`df` must contain `sample`, `XF`, `TC`, `nT`, `n`, `fnum`, `type`, `mut`, `reps`, and `tl` columns")
  }

  ## Check pold
  if(length(pold) > 1 & !multi_pold){
    stop("pold must be NULL or length 1 if multi_pold is TRUE")
  }

  if(!is.null(pold)){
    if(!is.numeric(pold)){
      stop("pold must be numeric")
    }else if(pold < 0){
      stop("pold must be >= 0")
    }else if(pold > 1){
      stop("pold must be <= 1")
    }
  }

  ## Extract info about # of replicates in each condition
  nMT <- max(df$mut)

  # nreps must be vector where each element i is number of replicates
  # in Exp_ID = i
  nreps <- rep(0, times = nMT)
  for(i in 1:nMT){
    nreps[i] <- max(df$reps[df$mut == i & df$type == 1])
  }


  ## Check pnew
  if(!is.null(pnew)){
    if(!all(is.numeric(pnew))){
      stop("All elements of pnew must be numeric")
    }else if((length(pnew) != sum(nreps)) & (length(pnew) != 1) ){
      stop("pnew must be a vector of length == number of s4U fed samples in dataset, or length 1.")
    }else if(!all(pnew > 0)){
      stop("All elements of pnew must be > 0")
    }else if(!all(pnew <=1)){
      stop("All elements of pnew must be <= 1")
    }
  }
  
  if(!is.null(pold)){
    if(!all(is.numeric(pold))){
      stop("All elements of pold must be numeric")
    }else if(((length(pold) != sum(nreps)) & (length(pold) != 1)) & multi_pold ){
      stop("pold must be a vector of length == number of s4U fed samples in dataset, or length 1.")
    }else if(!all(pold >= 0)){
      stop("All elements of pold must be >= 0")
    }else if(!all(pold <1)){
      stop("All elements of pold must be < 1")
    }
  }

  ## Check no_ctl
  if(!is.logical(no_ctl)){
    stop("no_ctl must be logical (TRUE or FALSE)")
  }
  

  ## Check read_cut
  if(!is.numeric(read_cut)){
    stop("read_cut must be numeric")
  }else if(read_cut < 0){
    stop("read_cut must be >= 0")
  }

  ## Check features_cut
  if(!is.numeric(features_cut)){
    stop("features_cut must be numeric")
  }else if(!is.integer(features_cut)){
    features_cut <- as.integer(features_cut)
  }

  if(features_cut <= 0){
    stop("features_cut must be > 0")
  }

  ## Check nbin
  if(!is.null(nbin)){
    if(!is.numeric(nbin)){
      stop("nbin must be numeric")
    }else if(!is.integer(nbin)){
      nbin <- as.integer(nbin)
    }

    if(nbin <= 0){
      stop("nbin must be > 0 and is preferably greater than or equal to 10")
    }
  }
  
  ### Check ztest
  if(!is.logical(ztest)){
    stop("ztest must be logical (TRUE or FALSE)")
  }


  ## Check prior_weight
  if(!is.numeric(prior_weight)){
    stop("prior_weight must be numeric")
  }else if(prior_weight < 1){
    stop("prior_weight must be >= 1. The larger the value, the more each feature's fraction new uncertainty estimate is shrunk towards the
         log10(read counts) vs. log(uncertainty) regression line.")
  }


  ## Check MLE
  if(!is.logical(MLE)){
    stop("MLE must be logical (TRUE or FALSE)")
  }

  if(!MLE){
    stop("Bayesian hypothesis estimation strategy has been removed, set MLE to TRUE.")
  }

  ## Check lower and upper
  if(!all(is.numeric(c(lower, upper)))){
    stop("lower and upper must be numeric")
  }else if(upper < lower){
    stop("upper must be > lower. Upper and lower represent the upper and lower bounds used by stats::optim for MLE")
  }

  ## Check mut_reg
  if(!is.numeric(mut_reg)){
    stop("mut_reg must be numeric")
  }else if(mut_reg < 0){
    stop("mut_reg must be greater than 0")
  }else if(mut_reg > 0.5){
    stop("mut_reg must be less than 0.5")
  }

  ## Check p_mean
  if(!is.numeric(p_mean)){
    stop("p_mean must be numeric")
  }else if(abs(p_mean) > 2){
    warning("p_mean is far from 0. This might bias estimates considerably; tread lightly.")
  }

  ## Check p_sd
  if(!is.numeric(p_sd)){
    stop("p_sd must be numeric")
  }else if(p_sd <= 0){
    stop("p_sd must be greater than 0")
  }else if(p_sd < 0.5){
    warning("p_sd is pretty small. This might bias estimates considerably; tread lightly.")
  }

  ## Check StanRate
  if(!is.logical(StanRate)){
    stop("StanRate must be logical (TRUE or FALSE)")
  }

  ## Check null_cutoff
  if(!is.numeric(null_cutoff)){
    stop("null_cutoff must be numeric")
  }else if(null_cutoff < 0){
    stop("null_cutoff must be 0 or positive")
  }else if(null_cutoff > 2){
    warning("You are testing against a null hypothesis |L2FC(kdeg)| greater than 2; this might be too conservative")
  }


  ## Check Chase
  if(!is.logical(Chase)){
    stop("Chase must be logical (TRUE or FALSE)")
  }


  # ESTIMATE MUTATION RATES ----------------------------------------------------

  # Helper functions that I will use on multiple occasions
  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))


  if(Long | kmeans){
    
    message("Estimating pnew with kmeans clustering")
    
    # Make sure package for clustering is installed
    if (!requireNamespace("Ckmeans.1d.dp", quietly = TRUE))
      stop("To use kmeans estimation strategy requires 'Ckmeans.1d.dp' package which cannot be found. Please install 'Ckmeans.1d.dp' using 'install.packages('Ckmeans.1d.dp')'.")

    # Filter out -s4U data, wont' use that here
    df <- df[df$type == 1,]

    # Find all combos of mut and reps
    id_dict <- df %>%
      dplyr::select(mut, reps) %>%
      dplyr::distinct()

    New_data_estimate <- data.frame(pnew = 0, pold = 0, mut = id_dict$mut,
                                    reps = id_dict$reps)
    # Loop throw rows of id_dict
    for(i in 1:nrow(id_dict)){
      rows <- which(df$mut == id_dict$mut[i] & df$reps == id_dict$reps[i])
      rates <- rep(df$TC[rows]/df$nT[rows], times = df$n[rows])

      means <- Ckmeans.1d.dp::Ckmeans.1d.dp(rates, k = 2)$centers


      New_data_estimate$pnew[i] <- max(means)
      New_data_estimate$pold[i] <- min(means)

    }


    
    
    if(multi_pold){
      New_data_estimate <- New_data_estimate[,c("pnew", "pold", "mut", "reps")]
      
      Old_data_estimate <- New_data_estimate[,c("pold", "mut", "reps")]
      New_data_estimate <- New_data_estimate[,c("pnew", "mut", "reps")]

    }else{
      
      if((sum(df$type == 0) > 0) & !no_ctl & !multi_pold){
        message("Estimating unlabeled mutation rate with -s4U data")
        
        #Old mutation rate estimation
        Mut_data <- df
        
        Old_data <- Mut_data[Mut_data$type == 0, ]
        
        Old_data <- Old_data %>% dplyr::ungroup() %>%
          dplyr::summarise(mutrate = sum(n*TC)/sum(n*nT))
        
        pold <- Old_data$mutrate
      }else{
        # average out pold
        pold <- mean(New_data_estimate$pold)
        New_data_estimate <- New_data_estimate[,c("pnew", "mut", "reps")]
        
      }

      
    }


  }else{
    ### Old mutation rate estimation
    if((is.null(pnew) | is.null(pold)) & StanRate ){ # use Stan


      mut_fit <- rstan::sampling(stanmodels$Mutrate_est, data = Stan_data, chains = 1)
    }

    ## Make data frame of pnew estimates
    if(is.null(pnew)){

      # Get estimates from Stan fit summary
      if(StanRate){
        
        message("Estimating pnew with Stan output")

        # Even if replicates are imbalanced (more replicates of a given experimental condition
        # than another), the number of mutation rates Stan will estimate = max(nreps)*nMT
        # nreps = vector containing number of replicates for ith experimental condition
        # nMT = number of experimental conditions
        # Therefore, have to remove imputed mutation rates

        nrep_mut <- max(nreps)

        U_df <- df[(df$type == 1) & (df$XF %in% unique(Stan_data$sdf$XF)),] %>% dplyr::group_by(mut, reps) %>%
          dplyr::summarise(avg_T = sum(nT*n)/sum(n))

        U_df <- U_df[order(U_df$mut, U_df$reps),]

        # In theory this may have too many entries; some could be imputed
        # So going to need to associate a replicate/experimental ID with
        # each row and remove those that are not real
        # Especially important since I am dividing by U_df$avg_T
        pnew <- exp(as.data.frame(rstan::summary(mut_fit, pars = "log_lambda_n")$summary)$mean)

        rep_theory <- rep(seq(from = 1, to = nrep_mut), times = nMT)
        mut_theory <- rep(seq(from = 1, to = nMT), each = nrep_mut)

        rep_actual <- unlist(lapply(nreps, function(x) seq(1, x)))
        mut_actual <- rep(1:nMT, times = nreps)

        pnewdf <- data.frame(pnew = pnew,
                             R = rep_theory,
                             E = mut_theory)

        truedf <- data.frame(R = rep_actual,
                             E = mut_actual)


        ## Get final pnew estimate vector
        pnewdf <- dplyr::right_join(pnewdf, truedf, by = c("R", "E"))

        pnew <- pnewdf$pnew/U_df$avg_T


        if(!is.null(pold)){
          rm(mut_fit)
        }


        New_data_estimate <- data.frame(mut_actual, rep_actual, pnew)
        colnames(New_data_estimate) <- c("mut", "reps", "pnew")



      }else{ # Use binomial mixture model to estimate new read estimate mutation rate

        message("Estimating pnew with likelihood maximization")

        # Binomial mixture likelihood
        mixture_lik <- function(param, TC, nT, n){

          logl <- sum(n*log(inv_logit(param[3])*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[1])^TC)*((1 -inv_logit(param[1]))^(nT-TC)) +  (1-inv_logit(param[3]))*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[2])^TC)*((1 - inv_logit(param[2]))^(nT-TC)) ) )

          return(-logl)

        }

        # Remove unlabeled controls
        df_pnew <- df[df$type == 1,]


        # Summarize out genes
        df_pnew <- df_pnew %>%
          dplyr::group_by(TC, nT, mut, reps) %>%
          dplyr::summarise(n = sum(n))

        # Check to make sure likelihood is evalutable
        df_check <- df_pnew %>%
          dplyr::mutate(fn = mixture_lik(param = c(-7, -2, 0),
                                         TC = TC,
                                         nT = nT,
                                         n = n))

        if(sum(!is.finite(df_check$fn)) > 0){
          stop("The binomial log-likelihood is not computable for some of your data, meaning the default mutation rate estimation strategy won't work. Rerun bakRFit with StanRateEst set to TRUE to remedy this problem.")
        }

        rm(df_check)


        low_ps <- c(-9, -9, -9)
        high_ps <- c(0, 0, 9)


        # Fit mixture model to each sample
        df_pnew <- df_pnew %>%
          dplyr::group_by(mut, reps) %>%
          dplyr::summarise(pnew = inv_logit(max(stats::optim(par=c(-7, -2, 0), mixture_lik, TC = TC, nT = nT, n = n, method = "L-BFGS-B", lower = low_ps, upper = high_ps)$par[1:2])) )

        New_data_estimate <- df_pnew


      }


    }else{ # Construct pmut data frame from User input

      message("Using provided pnew estimates")
      
      if(length(pnew) == 1){ # replicate the single pnew provided

        pnew_vect <- rep(pnew, times = sum(nreps) )

        ## Compatible with unbalanced replicates
        rep_vect <- unlist(lapply(nreps, function(x) seq(1, x)))

        mut_vect <- rep(1:nMT, times = nreps)

        New_data_estimate <- data.frame(mut_vect, rep_vect, pnew_vect)
        colnames(New_data_estimate) <- c("mut", "reps", "pnew")

      } else{ # use vector of pnews provided

        ## Compatible with unbalanced replicates
        rep_vect <- unlist(lapply(nreps, function(x) seq(1, x)))

        mut_vect <- rep(1:nMT, times = nreps)
        New_data_estimate <- data.frame(mut_vect, rep_vect, pnew)
        colnames(New_data_estimate) <- c("mut", "reps", "pnew")
      }
    }

    ### Estimate mutation rate in old reads
    if(is.null(pold)){
      
      if((sum(df$type == 0) > 0) & !no_ctl & !multi_pold){
        message("Estimating unlabeled mutation rate with -s4U data")
        
        #Old mutation rate estimation
        Mut_data <- df
        
        Old_data <- Mut_data[Mut_data$type == 0, ]
        
        Old_data <- Old_data %>% dplyr::ungroup() %>%
          dplyr::summarise(mutrate = sum(n*TC)/sum(n*nT))
        
        pold <- Old_data$mutrate
        
      }else if(StanRate){ # Use Stan to estimate rates

        message("Estimating unlabeled mutation rate with Stan output")
        
        ## Extract estimate of mutation rate in old reads from Stan fit
        nrep_mut <- max(df$reps)

        # Calculate U-content
        U_df <- df[(df$type == 1) & (df$XF %in% unique(Stan_data$sdf$XF)),] %>% dplyr::group_by(mut, reps) %>%
          dplyr::summarise(avg_T = sum(nT*n)/sum(n))

        U_df <- U_df[order(U_df$mut, U_df$reps),]

        # Vector of pold estimates for each sample
        pold <- exp(as.data.frame(rstan::summary(mut_fit, pars = "log_lambda_o")$summary)$mean)

        # Balanced replicate vectors
        # For matching Stan estimates to a replicate and experimental condition ID
        rep_theory <- rep(seq(from = 1, to = nrep_mut), times = nMT)
        mut_theory <- rep(seq(from = 1, to = nMT), each = nrep_mut)

        # Actual replicate vectors
        # Will be same as rep_theory and mut_theory if there are the same
        # number of replicates in each experimental condition
        rep_actual <- unlist(lapply(nreps, function(x) seq(1, x)))
        mut_actual <- rep(1:nMT, times = nreps)

        polddf <- data.frame(pold = pold,
                             R = rep_theory,
                             E = mut_theory)

        truedf <- data.frame(R = rep_actual,
                             E = mut_actual)


        # Filter out imputed data
        polddf <- dplyr::right_join(polddf, truedf, by = c("R", "E"))
        
        pold <- polddf$pold/U_df$avg_T
      
        
        rm(mut_fit)
        rm(U_df)
        
        Old_data_estimate <- data.frame(mut_actual, rep_actual, pold)
        colnames(Old_data_estimate) <- c("mut", "reps", "pold")
        
        

        if(!multi_pold){
          # Final estimate of mutation rate in old reads
          pold <- mean(Old_data_estimate$pold)
          

        }
       
      }else{

          message("Estimating unlabeled mutation rate with likelihood maximization")

          # Binomial mixture likelihood
          mixture_lik <- function(param, TC, nT, n){

            logl <- sum(n*log(inv_logit(param[3])*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[1])^TC)*((1 -inv_logit(param[1]))^(nT-TC)) +  (1-inv_logit(param[3]))*(factorial(nT)/(factorial(nT-TC)*factorial(TC)))*(inv_logit(param[2])^TC)*((1 - inv_logit(param[2]))^(nT-TC)) ) )

            return(-logl)

          }

          # Remove unlabeled controls
          df_pold <- df[df$type == 1,]


          # Summarize out genes
          df_pold <- df_pold %>%
            dplyr::group_by(TC, nT, mut, reps) %>%
            dplyr::summarise(n = sum(n))


          # Check to make sure likelihood is evalutable
          df_check <- df_pold %>%
            dplyr::mutate(fn = mixture_lik(param = c(-7, -2, 0),
                                           TC = TC,
                                           nT = nT,
                                           n = n))

          if(sum(!is.finite(df_check$fn)) > 0){
            stop("The binomial log-likelihood is not computable for some of your data, meaning the default mutation rate estimation strategy won't work. Rerun bakRFit with StanRateEst set to TRUE to remedy this problem.")
          }

          rm(df_check)


          low_ps <- c(-9, -9, -9)
          high_ps <- c(0, 0, 9)


          if(multi_pold){
            # Fit mixture model to each sample
            df_pold <- df_pold %>%
              dplyr::group_by(mut,reps) %>%
              dplyr::summarise(pold = inv_logit(min(stats::optim(par=c(-7, -2, 0), mixture_lik, TC = TC, nT = nT, n = n, method = "L-BFGS-B", lower = low_ps, upper = high_ps)$par[1:2])) )
            
            
            Old_data_estimate <- df_pold
            
          }else{
            # Fit mixture model to each sample
            df_pold <- df_pold %>%
              dplyr::group_by(mut,reps) %>%
              dplyr::summarise(pold = inv_logit(min(stats::optim(par=c(-7, -2, 0), mixture_lik, TC = TC, nT = nT, n = n, method = "L-BFGS-B", lower = low_ps, upper = high_ps)$par[1:2])) ) %>%
              dplyr::ungroup() %>%
              dplyr::summarise(pold = mean(pold))
            
            
            
            pold <- df_pold$pold

          }



        }
    }else{
      
      if(multi_pold){
        message("Using provided pold estimates")
        
        if(length(pold) == 1){ # replicate the single pnew provided
          
          pold_vect <- rep(pold, times = sum(nreps) )
          
          ## Compatible with unbalanced replicates
          rep_vect <- unlist(lapply(nreps, function(x) seq(1, x)))
          
          mut_vect <- rep(1:nMT, times = nreps)
          
          Old_data_estimate <- data.frame(mut_vect, rep_vect, pold_vect)
          colnames(Old_data_estimate) <- c("mut", "reps", "pold")
          
        } else{ # use vector of pnews provided
          
          ## Compatible with unbalanced replicates
          rep_vect <- unlist(lapply(nreps, function(x) seq(1, x)))
          
          mut_vect <- rep(1:nMT, times = nreps)
          Old_data_estimate <- data.frame(mut_vect, rep_vect, pold)
          colnames(Old_data_estimate) <- c("mut", "reps", "pold")
        }
      }else{
        message("Using provided pold estimate")
      
      }


    }
  }

  if(!multi_pold){
    New_data_estimate$pold <- pold
    message(paste0(c("Estimated pnews and polds for each sample are:", utils::capture.output(New_data_estimate)), collapse = "\n"))

  }else{
    New_data_estimate <- dplyr::inner_join(New_data_estimate, Old_data_estimate, by = c("reps", "mut"))
    message(paste0(c("Estimated pnews and polds for each sample are:", utils::capture.output(New_data_estimate)), collapse = "\n"))
    
  }


  if(!all(New_data_estimate$pnew - New_data_estimate$pold > 0)){
    stop("All pnew must be > pold; did you input an unusually large pold?")
  }

  # Mutation rate estimates
  pmuts_list <- New_data_estimate



  # PREP DATA FOR ANALYSES -----------------------------------------------------
    # i) Filter out any -s4U control samples
    # ii) Make lookup table mapping label times to experimental conditions
    # iii) Make lookup table mapping sample names to experimental conditions and replicate IDs
    # iv) Make lookup table mapping feature number to feature name

  Mut_data <- df

  # Filter out -s4U control samples
  Mut_data <- Mut_data[Mut_data$type == 1,]

  tl_df <- Mut_data %>% dplyr::select(mut, tl) %>%
    dplyr::distinct()

  # Label time lookup table
  tl_df <- tl_df[order(tl_df$mut),]

  tl <- tl_df$tl

  ngene <- max(Mut_data$fnum)
  num_conds <- max(Mut_data$mut)

  nreps <- rep(0, times = num_conds)
  for(i in 1:num_conds){
    nreps[i] <- max(Mut_data$reps[Mut_data$mut == i & Mut_data$type == 1])
  }

  # Sample characteristics lookup table
  sample_lookup <- Mut_data[, c("sample", "mut", "reps")] %>% dplyr::distinct()

  # Feature lookup table
  feature_lookup <- Mut_data[,c("fnum", "XF")] %>% dplyr::distinct()



  # ESTIMATE FRACTION NEW ------------------------------------------------------


  ## Estimate fraction new in each replicate using U-content adjusted Poisson model
  message("Estimating fraction labeled")

  Mut_data <- dplyr::left_join(Mut_data, New_data_estimate, by = c("mut", "reps"))

  if(Long){

    # Calculate variance with delta approximation using beta distribution
    var_calc <- function(alpha, beta){

      EX <- alpha/(alpha + beta)
      VX <- (alpha*beta)/(((alpha + beta)^2)*(alpha + beta + 1))


      var1 <- (((1/EX) + 1/(1 - EX))^2)*VX
      var2 <- -(-1/(4*(EX^2)) + 1/(4*((1-EX)^2)))^2
      var3 <- VX^2

      totvar <- var1 - var2*var3

      return(totvar)
    }

    # convert logit(fn) se -> log(kdeg) se
    logit_to_log <- function(lfn, varx){

      ## composition strategy
      #se <- sqrt(((((1 + exp(lfn))^-2)*varx)/(-log(1 - inv_logit(lfn))))*varx)

      ## full monte
      part1 <- 1/(log( (exp(-lfn))/(1 + exp(-lfn)))*(1 + exp(-lfn)))
      se <- sqrt((part1^2)*varx)

      return(se)

    }


    # 1) Add new read identifier
    cutoff <- function(pnew, pold){
      return(((logit(pnew) + logit(pold))/2 + logit((pnew + pold)/2))/2 )
    }


    Mut_data <- Mut_data %>%
      dplyr::mutate(new = ifelse(logit(TC/nT) > cutoff(pnew, pold), TRUE, FALSE))

    # 2) Calculate number of new reads and total number of reads in each replicate for each feature
    Mut_data <- Mut_data %>%
      dplyr::group_by(fnum, mut, reps) %>%
      dplyr::summarise(nreads = sum(n),
                       n_new = sum(new))


    # 3) Calculate fraction new and uncertainty
    Mut_data_est <- Mut_data %>%
      dplyr::mutate(Fn_rep_est = (n_new + 1)/(nreads + 1),
                    logit_fn_rep = logit(Fn_rep_est),
                    logit_fn_se = sqrt(var_calc(n_new + 1, nreads + 1)),
                    kd_rep_est = -log(1 - Fn_rep_est)/tl[mut],
                    log_kd_rep_est = log(kd_rep_est),
                    log_kd_se = logit_to_log(logit_fn_rep, logit_fn_se^2))



  }else{ # MLE

    # Likelihood function for mixture model
    mixed_lik <- function(lam_n, lam_o, TC, n, logit_fn){
      logl <- sum(n*log(inv_logit(logit_fn)*(lam_n^TC)*exp(-lam_n) + (1-inv_logit(logit_fn))*(lam_o^TC)*exp(-lam_o) )) + log(stats::dnorm(logit_fn, mean = p_mean, sd = p_sd))
      return(-logl)
    }

    # Calculate logit(fn) MLE
    Mut_data_est <- Mut_data %>% dplyr::ungroup() %>% dplyr::mutate(lam_n = pnew*nT, lam_o = pold*nT) %>%
      dplyr::group_by(fnum, mut, reps, TC) %>%
      dplyr::summarise(lam_n = sum(lam_n*n)/sum(n), lam_o = sum(lam_o*n)/sum(n),
                       n = sum(n), .groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fnum, mut, reps) %>%
      dplyr::summarise(logit_fn_rep = stats::optim(0, mixed_lik, TC = TC, n = n, lam_n = sum(lam_n*n)/sum(n), lam_o = sum(lam_o*n)/sum(n), method = "L-BFGS-B", lower = lower, upper = upper)$par, nreads =sum(n), .groups = "keep") %>%
      dplyr::mutate(logit_fn_rep = ifelse(logit_fn_rep == lower, stats::runif(1, lower-0.2, lower), ifelse(logit_fn_rep == upper, stats::runif(1, upper, upper+0.2), logit_fn_rep))) %>%
      dplyr::ungroup()

    ## Look for numerical instabilities
    instab_df <- Mut_data_est %>% dplyr::filter(abs(logit_fn_rep) > upper )

    fnum_instab <- instab_df$fnum
    mut_instab <- instab_df$mut
    reps_instab <- instab_df$reps

    # Replace unstable estimates with estimate + jitter so as to not underestimate replicate variability
    instab_est <- dplyr::right_join(Mut_data, instab_df, by = c("fnum", "mut", "reps")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(totTC = TC*n, totU = nT*n) %>%
      dplyr::group_by(fnum, mut, reps) %>%
      dplyr::summarise(tot_mut = sum(totTC), totUs = sum(totU), pnew = mean(pnew), pold = mean(pold)) %>% dplyr::ungroup() %>%
      dplyr::mutate(avg_mut = tot_mut/totUs,
                    Fn_rep_est = (avg_mut - (1-mut_reg)*pold)/((1+mut_reg)*pnew - (1-mut_reg)*pold)) %>%
      dplyr::mutate(Fn_rep_est = ifelse(Fn_rep_est > 1, stats::runif(1, min = inv_logit(upper-0.1), max = 1) , ifelse(Fn_rep_est < 0, stats::runif(1, min = 0, max = inv_logit(lower+0.1)), Fn_rep_est ) )) %>%
      dplyr::mutate(logit_fn_rep = logit(Fn_rep_est)) %>%
      dplyr::select(fnum, mut, reps, logit_fn_rep)

    # Generate estimates on other useful scales

    if(Chase){
      Mut_data_est <- dplyr::left_join(Mut_data_est, instab_est, by = c("fnum", "mut", "reps")) %>%
        dplyr::mutate(logit_fn_rep = ifelse(abs(logit_fn_rep.x) > upper, -logit_fn_rep.y, -logit_fn_rep.x)) %>%
        dplyr::select(fnum, mut, reps, nreads, logit_fn_rep) %>%
        dplyr::mutate(Fn_rep_est = inv_logit(logit_fn_rep)) %>%
        dplyr::mutate(kd_rep_est = -log(1 - Fn_rep_est)/tl[mut])
    }else{
      Mut_data_est <- dplyr::left_join(Mut_data_est, instab_est, by = c("fnum", "mut", "reps")) %>%
        dplyr::mutate(logit_fn_rep = ifelse(abs(logit_fn_rep.x) > upper, logit_fn_rep.y, logit_fn_rep.x)) %>%
        dplyr::select(fnum, mut, reps, nreads, logit_fn_rep) %>%
        dplyr::mutate(Fn_rep_est = inv_logit(logit_fn_rep)) %>%
        dplyr::mutate(kd_rep_est = -log(1 - Fn_rep_est)/tl[mut])
    }

    message("Estimating per replicate uncertainties")

    Mut_data <- dplyr::left_join(Mut_data, Mut_data_est[, c("kd_rep_est" ,"logit_fn_rep", "fnum", "mut", "reps")], by = c("fnum", "mut", "reps"))

    ## Estimate Fisher Info and uncertainties
    if(Chase){
      Mut_data <- Mut_data %>% dplyr::ungroup() %>%
        dplyr::group_by(fnum, mut, reps, TC, pnew, pold, logit_fn_rep, kd_rep_est) %>%
        dplyr::summarise(U_cont = sum(nT*n)/sum(n), n = sum(n), .groups = "keep") %>%
        dplyr::mutate(Exp_l_fn = exp(-logit_fn_rep)) %>%
        dplyr::mutate(Inv_Fisher_Logit_3 = 1/(((pnew/pold)^TC)*exp(-(U_cont)*(pnew - pold)) - 1 )) %>%
        dplyr::mutate(Inv_Fisher_Logit_1 = 1 + Exp_l_fn ) %>%
        dplyr::mutate(Inv_Fisher_Logit_2 = ((1 + Exp_l_fn)^2)/Exp_l_fn) %>%
        dplyr::mutate(Fisher_lkd_num = tl[mut]*kd_rep_est*( ((pnew*U_cont)^TC)*exp(-pnew*U_cont) - ( ((pold*U_cont)^TC)*exp(-pold*U_cont) ) ) ) %>%
        dplyr::mutate(Fisher_lkd_den = (exp(kd_rep_est*tl[mut]) - 1)*((pnew*U_cont)^TC)*exp(-pnew*U_cont) + ((pold*U_cont)^TC)*exp(-pold*U_cont) ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(fnum, mut, reps) %>%
        dplyr::summarise(Fisher_kdeg = sum(n*((Fisher_lkd_num/Fisher_lkd_den)^2))/sum(n),
                         Fisher_Logit = sum(n/((Inv_Fisher_Logit_1 + Inv_Fisher_Logit_2*Inv_Fisher_Logit_3)^2))/sum(n),
                         tot_n = sum(n)) %>% #,
        #Fisher_fn = sum(n*((Fisher_fn_num/Fisher_fn_den)^2)), tot_n = sum(n)) %>%
        dplyr::mutate(log_kd_se = 1/sqrt(tot_n*Fisher_kdeg),
                      Logit_fn_se = 1/sqrt(tot_n*Fisher_Logit)) %>%
        dplyr::mutate(log_kd_se = ifelse(log_kd_se > se_max, se_max, log_kd_se),
                      Logit_fn_se = ifelse(Logit_fn_se > se_max, se_max, Logit_fn_se))#, Fn_se = 1/sqrt(tot_n*Fisher_fn))

    }else{
      Mut_data <- Mut_data %>% dplyr::ungroup() %>%
        dplyr::group_by(fnum, mut, reps, TC, pnew, pold, logit_fn_rep, kd_rep_est) %>%
        dplyr::summarise(U_cont = sum(nT*n)/sum(n), n = sum(n), .groups = "keep") %>%
        dplyr::mutate(Exp_l_fn = exp(logit_fn_rep)) %>%
        dplyr::mutate(Inv_Fisher_Logit_3 = 1/(((pnew/pold)^TC)*exp(-(U_cont)*(pnew - pold)) - 1 )) %>%
        dplyr::mutate(Inv_Fisher_Logit_1 = 1 + Exp_l_fn ) %>%
        dplyr::mutate(Inv_Fisher_Logit_2 = ((1 + Exp_l_fn)^2)/Exp_l_fn) %>%
        dplyr::mutate(Fisher_lkd_num = tl[mut]*kd_rep_est*( ((pnew*U_cont)^TC)*exp(-pnew*U_cont) - ( ((pold*U_cont)^TC)*exp(-pold*U_cont) ) ) ) %>%
        dplyr::mutate(Fisher_lkd_den = (exp(kd_rep_est*tl[mut]) - 1)*((pnew*U_cont)^TC)*exp(-pnew*U_cont) + ((pold*U_cont)^TC)*exp(-pold*U_cont) ) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(fnum, mut, reps) %>%
        dplyr::summarise(Fisher_kdeg = sum(n*((Fisher_lkd_num/Fisher_lkd_den)^2))/sum(n),
                         Fisher_Logit = sum(n/((Inv_Fisher_Logit_1 + Inv_Fisher_Logit_2*Inv_Fisher_Logit_3)^2))/sum(n),
                         tot_n = sum(n)) %>% #,
        #Fisher_fn = sum(n*((Fisher_fn_num/Fisher_fn_den)^2)), tot_n = sum(n)) %>%
        dplyr::mutate(log_kd_se = 1/sqrt(tot_n*Fisher_kdeg),
                      Logit_fn_se = 1/sqrt(tot_n*Fisher_Logit)) %>%
        dplyr::mutate(log_kd_se = ifelse(log_kd_se > se_max, se_max, log_kd_se),
                      Logit_fn_se = ifelse(Logit_fn_se > se_max, se_max, Logit_fn_se))#, Fn_se = 1/sqrt(tot_n*Fisher_fn))

    }



    Mut_data_est$logit_fn_se <- Mut_data$Logit_fn_se
    Mut_data_est$log_kd_se <- Mut_data$log_kd_se

    Mut_data_est <- Mut_data_est %>% dplyr::mutate(log_kd_rep_est = log(kd_rep_est))


  }
  
  logit_fn <- Mut_data_est$logit_fn_rep

  # Report numerical instabilities from maximum likelihood estimation
  if(!(all(logit_fn < upper) & all(logit_fn > lower))){
    num_unstable <- sum(logit_fn >= upper) + sum(logit_fn <= lower)
    tot_ests <- length(logit_fn)
    prcnt_unstable <- round((num_unstable/tot_ests)*100,3)
    warning(paste0(num_unstable, " out of " , tot_ests, " (", prcnt_unstable,"%)", " logit(fn) estimates are at the upper or lower bounds set. These likely represent features with limited data or extremely stable/unstable features. If the number of boundary estimates is concerning, try increasing the magnitude of upper and lower."))
  }


  # ESTIMATE VARIANCE VS. READ COUNT TREND -------------------------------------

  
  out <- avg_and_regularize(Mut_data_est, nreps, sample_lookup, feature_lookup,
                            nbin = nbin, NSS = NSS, 
                            BDA_model = BDA_model, null_cutoff = null_cutoff,
                            Mutrates = pmuts_list, ztest = ztest)
  
  return(out)

}

#' Efficiently average replicates of nucleotide recoding data and regularize
#'
#' \code{avg_and_regularize} pools and regularizes replicate estimates of kinetic parameters. There are two key steps in this
#' downstream analysis. 1st, the uncertainty for each feature is used to fit a linear ln(uncertainty) vs. log10(read depth) trend,
#' and uncertainties for individual features are shrunk towards the regression line. The uncertainty for each feature is a combination of the
#' Fisher Information asymptotic uncertainty as well as the amount of variability seen between estimates. Regularization of uncertainty
#' estimates is performed using the analytic results of a Normal distribution likelihood with known mean and unknown variance and conjugate
#' priors. The prior parameters are estimated from the regression and amount of variability about the regression line. The strength of
#' regularization can be tuned by adjusting the \code{prior_weight} parameter, with larger numbers yielding stronger shrinkage towards
#' the regression line. The 2nd step is to regularize the average kdeg estimates. This is done using the analytic results of a
#' Normal distribution likelihood model with unknown mean and known variance and conjugate priors. The prior parameters are estimated from the
#' population wide kdeg distribution (using its mean and standard deviation as the mean and standard deviation of the normal prior).
#' In the 1st step, the known mean is assumed to be the average kdeg, averaged across replicates and weighted by the number of reads
#' mapping to the feature in each replicate. In the 2nd step, the known variance is assumed to be that obtained following regularization
#' of the uncertainty estimates.
#'
#' Effect sizes (changes in kdeg) are obtained as the difference in log(kdeg) means between the reference and experimental
#' sample(s), and the log(kdeg)s are assumed to be independent so that the variance of the effect size is the sum of the
#' log(kdeg) variances. P-values assessing the significance of the effect size are obtained using a moderated t-test with number
#' of degrees of freedom determined from the uncertainty regression hyperparameters and are adjusted for multiple testing using the Benjamini-
#' Hochberg procedure to control false discovery rates (FDRs).
#'
#' In some cases, the assumed ODE model of RNA metabolism will not accurately model the dynamics of a biological system being analyzed.
#' In these cases, it is best to compare logit(fraction new)s directly rather than converting fraction new to log(kdeg).
#' This analysis strategy is implemented when \code{NSS} is set to TRUE. Comparing logit(fraction new) is only valid
#' If a single metabolic label time has been used for all samples. For example, if a label time of 1 hour was used for NR-seq
#' data from WT cells and a 2 hour label time was used in KO cells, this comparison is no longer valid as differences in
#' logit(fraction new) could stem from differences in kinetics or label times.
#'
#'
#' @param Mut_data_est Dataframe with fraction new estimation information. Required columns are: 
#' \itemize{
#'  \item fnum; numerical ID of feature
#'  \item reps; numerical ID of replicate
#'  \item mut; numerical ID of experimental condition (Exp_ID)
#'  \item logit_fn_rep; logit(fn) estimate
#'  \item kd_rep_est; kdeg estimate
#'  \item log_kd_rep_est; log(kdeg) estimate
#'  \item logit_fn_se; logit(fn) estimate uncertainty
#'  \item log_kd_se; log(kdeg) estimate uncertainty
#' }
#' @param nbin Number of bins for mean-variance relationship estimation. If NULL, max of 10 or (number of logit(fn) estimates)/100 is used
#' @param nreps Vector of number of replicates in each experimental condition
#' @param sample_lookup Dictionary mapping sample names to various experimental details
#' @param feature_lookup Dictionary mapping feature IDs to original feature names
#' @param null_cutoff bakR will test the null hypothesis of |effect size| < |null_cutoff|
#' @param NSS Logical; if TRUE, logit(fn)s are compared rather than log(kdeg) so as to avoid steady-state assumption.
#' @param Chase Logical; Set to TRUE if analyzing a pulse-chase experiment. If TRUE, kdeg = -ln(fn)/tl where fn is the fraction of
#' reads that are s4U (more properly referred to as the fraction old in the context of a pulse-chase experiment)
#' @param BDA_model Logical; if TRUE, variance is regularized with scaled inverse chi-squared model. Otherwise a log-normal
#' model is used.
#' @param Mutrates List containing new and old mutation rate estimates
#' @param ztest TRUE; if TRUE, then a z-test is used for p-value calculation rather than the more conservative moderated t-test. 
#' @return List with dataframes providing information about replicate-specific and pooled analysis results. The output includes:
#' \itemize{
#'  \item Fn_Estimates; dataframe with estimates for the fraction new and fraction new uncertainty for each feature in each replicate.
#'  The columns of this dataframe are:
#'  \itemize{
#'   \item Feature_ID; Numerical ID of feature
#'   \item Exp_ID; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item Replicate; Numerical ID for replicate
#'   \item logit_fn; logit(fraction new) estimate, unregularized
#'   \item logit_fn_se; logit(fraction new) uncertainty, unregularized and obtained from Fisher Information
#'   \item nreads; Number of reads mapping to the feature in the sample for which the estimates were obtained
#'   \item log_kdeg; log of degradation rate constant (kdeg) estimate, unregularized
#'   \item kdeg; degradation rate constant (kdeg) estimate
#'   \item log_kd_se; log(kdeg) uncertainty, unregularized and obtained from Fisher Information
#'   \item sample; Sample name
#'   \item XF; Original feature name
#'  }
#'  \item Regularized_ests; dataframe with average fraction new and kdeg estimates, averaged across the replicates and regularized
#'  using priors informed by the entire dataset. The columns of this dataframe are:
#'  \itemize{
#'   \item Feature_ID; Numerical ID of feature
#'   \item Exp_ID; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item avg_log_kdeg; Weighted average of log(kdeg) from each replicate, weighted by sample and feature-specific read depth
#'   \item sd_log_kdeg; Standard deviation of the log(kdeg) estimates
#'   \item nreads; Total number of reads mapping to the feature in that condition
#'   \item sdp; Prior standard deviation for fraction new estimate regularization
#'   \item theta_o; Prior mean for fraction new estimate regularization
#'   \item sd_post; Posterior uncertainty
#'   \item log_kdeg_post; Posterior mean for log(kdeg) estimate
#'   \item kdeg; exp(log_kdeg_post)
#'   \item kdeg_sd; kdeg uncertainty
#'   \item XF; Original feature name
#'  }
#'  \item Effects_df; dataframe with estimates of the effect size (change in logit(fn)) comparing each experimental condition to the
#'  reference sample for each feature. This dataframe also includes p-values obtained from a moderated t-test. The columns of this
#'  dataframe are:
#'  \itemize{
#'   \item Feature_ID; Numerical ID of feature
#'   \item Exp_ID; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item L2FC(kdeg); Log2 fold change (L2FC) kdeg estimate or change in logit(fn) if NSS TRUE
#'   \item effect; LFC(kdeg)
#'   \item se; Uncertainty in L2FC_kdeg
#'   \item pval; P-value obtained using effect_size, se, and a z-test
#'   \item padj; pval adjusted for multiple testing using Benjamini-Hochberg procedure
#'   \item XF; Original feature name
#'  }
#'  \item Mut_rates; list of two elements. The 1st element is a dataframe of s4U induced mutation rate estimates, where the mut column
#'  represents the experimental ID and the rep column represents the replicate ID. The 2nd element is the single background mutation
#'  rate estimate used
#'  \item Hyper_Parameters; vector of two elements, named a and b. These are the hyperparameters estimated from the uncertainties for each
#'  feature, and represent the two parameters of a Scaled Inverse Chi-Square distribution. Importantly, a is the number of additional
#'  degrees of freedom provided by the sharing of uncertainty information across the dataset, to be used in the moderated t-test.
#'  \item Mean_Variance_lms; linear model objects obtained from the uncertainty vs. read count regression model. One model is run for each Exp_ID
#' }
#' @importFrom magrittr %>%
avg_and_regularize <- function(Mut_data_est, nreps, sample_lookup, feature_lookup,
                               nbin = NULL, NSS = FALSE, Chase = FALSE,
                               BDA_model = FALSE, null_cutoff = 0,
                               Mutrates = NULL, ztest = FALSE){
  

  ### Check Mut_data_est
  expected_cols <- c("nreads", "fnum", "reps",
                              "mut", "logit_fn_rep", "kd_rep_est", "log_kd_rep_est", 
                              "logit_fn_se", "log_kd_se")
  
  if(!is.data.frame(Mut_data_est)){
    stop("Mut_data_est must be a data frame!")
  }
  
  if(!all(expected_cols %in% colnames(Mut_data_est))){
    stop("Mut_data_est is lacking some necessary columns. See documentation (?avg_and_regularize())
         for list of necessary columns")
  }
  
  ### Check nreps
  if(!is.numeric(nreps)){
    stop("nreps must be a numeric vector of the number of replicates in each Exp_ID!")
  }

  ### Check sample_lookup
  expected_cols <- c("sample", "mut", "reps")
  if(!all(expected_cols %in% colnames(sample_lookup))){
    stop("sample_lookup is lacking necessary columns. It should contain three columns: 1) sample, which is the name of the same. 2) mut, which is the Exp_ID, and 3) reps which is the numerical replicate ID ")
  }
  
  if(!is.data.frame(sample_lookup)){
    stop("sample_lookup must be a data frame!")
  }
  
  ### Check feature_lookup
  expected_cols <- c("XF", "fnum")
  if(!all(expected_cols %in% colnames(feature_lookup))){
    stop("feature_lookup is lacking necessary columns. It should contain at least two columns: 1) XF, the name of each feature analyzed, and 2) fnum, the numerical ID for the corresponding feature")
  }
  
  if(!is.data.frame(feature_lookup)){
    stop("feature_lookup must be a data frame!")
  }
  
  ## Check nbin
  if(!is.null(nbin)){
    if(!is.numeric(nbin)){
      stop("nbin must be numeric")
    }else if(!is.integer(nbin)){
      nbin <- as.integer(nbin)
    }
    
    if(nbin <= 0){
      stop("nbin must be > 0 and is preferably greater than or equal to 10")
    }
  }
  
  ## Check NSS
  if(!is.logical(NSS)){
    stop("NSS must be logical (TRUE or FALSE)")
  }
  
  ### Check Chase
  if(!is.logical(Chase)){
    stop("Chase must be logical (TRUE or FALSE)")
  }
  
  ### Check BDA_model
  if(!is.logical(BDA_model)){
    stop("BDA_model must be logical (TRUE or FALSE)")
  }
  
  ### Check ztest
  if(!is.logical(ztest)){
    stop("ztest must be logical (TRUE or FALSE)")
  }
  
  ## Check null_cutoff
  if(!is.numeric(null_cutoff)){
    stop("null_cutoff must be numeric")
  }else if(null_cutoff < 0){
    stop("null_cutoff must be 0 or positive")
  }else if(null_cutoff > 2){
    warning("You are testing against a null hypothesis |L2FC(kdeg)| greater than 2; this might be too conservative")
  }
  
  
  # Bind variables locally to resolve devtools::check() Notes
  fnum <- mut <- logit_fn_rep <- bin_ID <- kd_sd_log <- intercept <- NULL
  slope <- log_kd_rep_est <- avg_log_kd <- sd_log_kd <- sd_post <- NULL
  sdp <- theta_o <- log_kd_post <- effect_size <- effect_size_std_error <- NULL
  
  
  # Helper functions that I will use on multiple occasions
  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))
  
  
  # ESTIMATE VARIANCE VS. READ COUNT TREND -------------------------------------
  
  ngene <- max(Mut_data_est$fnum)
  nMT <- max(Mut_data_est$mut)
  
  ## Now affiliate each fnum, mut with a bin Id based on read counts,
  ## bin data by bin_ID and average log10(reads) and log(sd(logit_fn))
  
  if(is.null(nbin)){
    nbin <- max(c(round(ngene*sum(nreps)/100), 10))
  }
  
  message("Estimating read count-variance relationship")
  
  
  if(NSS){
    Binned_data <- Mut_data_est %>% dplyr::group_by(fnum, mut) %>%
      dplyr::summarise(nreads = sum(nreads), kd_sd_log = log(sqrt(1/sum(1/((stats::sd(logit_fn_rep)^2) +logit_fn_se^2 ) ) ) )) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(bin_ID = as.numeric(Hmisc::cut2(nreads, g = nbin))) %>% dplyr::group_by(bin_ID, mut) %>%
      dplyr::summarise(avg_reads = mean(log10(nreads)), avg_sd = mean(kd_sd_log))
    
    ## Estimate read count vs. replicate variability trend
    
    lm_list <- vector("list", length = nMT)
    lm_var <- lm_list
    
    # One linear model for each experimental condition
    for(i in 1:nMT){
      heterosked_lm <- stats::lm(avg_sd ~ avg_reads, data = Binned_data[Binned_data$mut == i,] )
      h_int <- summary(heterosked_lm)$coefficients[1,1]
      h_slope <- summary(heterosked_lm)$coefficients[2,1]
      
      if(h_slope > 0){
        h_slope <- 0
        h_int <- mean(Binned_data$avg_sd[Binned_data$mut == i])
      }
      
      lm_list[[i]] <- c(h_int, h_slope)
      
      lm_var[[i]] <- stats::var(stats::residuals(heterosked_lm))
    }
    
    # Put linear model fit extrapolation into convenient data frame
    true_vars <-  Mut_data_est %>% dplyr::group_by(fnum, mut) %>%
      dplyr::summarise(nreads = sum(nreads), kd_sd_log = log(sqrt(1/sum(1/((stats::sd(logit_fn_rep)^2) + logit_fn_se^2 ) ) ) )) %>%
      dplyr::ungroup() %>% dplyr::group_by(fnum, mut) %>% dplyr::mutate(slope = lm_list[[mut]][2], intercept = lm_list[[mut]][1]) %>%
      dplyr::group_by(mut) %>%
      dplyr::summarise(true_var = stats::var(kd_sd_log - (intercept + slope*log10(nreads) ) ))
    
    log_kd <- as.vector(Mut_data_est$log_kd_rep_est)
    logit_fn <- as.vector(Mut_data_est$logit_fn_rep)
    
  }else{
    Binned_data <- Mut_data_est %>% dplyr::group_by(fnum, mut) %>%
      dplyr::summarise(nreads = sum(nreads), kd_sd_log = log(sqrt(1/sum(1/((stats::sd(log_kd_rep_est)^2) + log_kd_se^2 ) ) ) )) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(bin_ID = as.numeric(Hmisc::cut2(nreads, g = nbin))) %>% dplyr::group_by(bin_ID, mut) %>%
      dplyr::summarise(avg_reads = mean(log10(nreads)), avg_sd = mean(kd_sd_log))
    
    ## Estimate read count vs. replicate variability trend
    
    lm_list <- vector("list", length = nMT)
    lm_var <- lm_list
    
    # One linear model for each experimental condition
    for(i in 1:nMT){
      heterosked_lm <- stats::lm(avg_sd ~ avg_reads, data = Binned_data[Binned_data$mut == i,] )
      h_int <- summary(heterosked_lm)$coefficients[1,1]
      h_slope <- summary(heterosked_lm)$coefficients[2,1]
      
      if(h_slope > 0){
        h_slope <- 0
        h_int <- mean(Binned_data$avg_sd[Binned_data$mut == i])
      }
      
      lm_list[[i]] <- c(h_int, h_slope)
      
      lm_var[[i]] <- stats::var(stats::residuals(heterosked_lm))
    }
    
    # Put linear model fit extrapolation into convenient data frame
    true_vars <-  Mut_data_est %>% dplyr::group_by(fnum, mut) %>%
      dplyr::summarise(nreads = sum(nreads), kd_sd_log = log(sqrt(1/sum(1/((stats::sd(log_kd_rep_est)^2) + log_kd_se^2 ) ) ) )) %>%
      dplyr::ungroup() %>% dplyr::group_by(fnum, mut) %>% dplyr::mutate(slope = lm_list[[mut]][2], intercept = lm_list[[mut]][1]) %>%
      dplyr::group_by(mut) %>%
      dplyr::summarise(true_var = stats::var(kd_sd_log - (intercept + slope*log10(nreads) ) ),
                       true_nat_var = stats::var(exp(kd_sd_log)^2 - exp((intercept + slope*log10(nreads) ))^2 ))
    
    log_kd <- as.vector(Mut_data_est$log_kd_rep_est)
    logit_fn <- as.vector(Mut_data_est$logit_fn_rep)
    
  }
  
  
  
  # PREP DATA FOR REGULARIZATION -----------------------------------------------
  
  # Vectors of use
  kd_estimate <- as.vector(Mut_data_est$kd_rep_est)
  log_kd_se <- as.vector(Mut_data_est$log_kd_se)
  logit_fn_se <- as.vector(Mut_data_est$logit_fn_se)
  Replicate <- as.vector(Mut_data_est$reps)
  Condition <- as.vector(Mut_data_est$mut)
  Gene_ID <- as.vector(Mut_data_est$fnum)
  nreads <- as.vector(Mut_data_est$nreads)
  
  rm(Mut_data_est)
  
  # Create new data frame with kinetic parameter estimates and uncertainties
  df_fn <- data.frame(logit_fn, logit_fn_se, Replicate, Condition, Gene_ID, nreads, log_kd, kd_estimate, log_kd_se)
  
  
  # Remove vectors no longer of use
  rm(logit_fn)
  rm(logit_fn_se)
  rm(Replicate)
  rm(Condition)
  rm(Gene_ID)
  rm(nreads)
  
  df_fn <- df_fn[order(df_fn$Gene_ID, df_fn$Condition, df_fn$Replicate),]
  
  
  # Relabel logit(fn) as log(kdeg) for steady-state independent analysis
  if(NSS){
    colnames(df_fn) <- c("log_kd", "log_kd_se", "Replicate", "Condition", "Gene_ID", "nreads", "logit_fn", "kd_estimate", "logit_fn_se")
    
    df_fn$kd_estimate <- inv_logit(df_fn$log_kd)
  }
  
  
  
  # AVERAGE REPLICATE DATA AND REGULARIZE WITH INFORMATIVE PRIORS --------------
  
  message("Averaging replicate data and regularizing estimates")
  
  #Average over replicates and estimate hyperparameters
  avg_df_fn_bayes <- df_fn %>% dplyr::group_by(Gene_ID, Condition) %>%
    dplyr::summarize(avg_log_kd = stats::weighted.mean(log_kd, 1/log_kd_se),
                     sd_log_kd = sqrt(1/sum(1/((stats::sd(log_kd)^2) + log_kd_se^2 ) ) ),
                     nreads = sum(nreads)) %>% dplyr::ungroup() %>%
    dplyr::group_by(Condition) %>%
    dplyr::mutate(sdp = stats::sd(avg_log_kd)) %>%
    dplyr::mutate(theta_o = mean(avg_log_kd)) %>%
    dplyr::ungroup()
  

  
  #Calcualte population averages
  # sdp <- sd(avg_df_fn$avg_logit_fn) # Will be prior sd in regularization of mean
  # theta_o <- mean(avg_df_fn$avg_logit_fn) # Will be prior mean in regularization of mean
  
  var_pop <- mean(avg_df_fn_bayes$sd_log_kd^2) # Will be prior mean in regularization of sd
  var_of_var <- mean(true_vars$true_nat_var) # Will be prior variance in regularization of sd
  
  
  ## Regularize standard deviation estimate
  
  # Estimate hyperpriors with method of moments
  
  a_hyper <- 2*(var_pop^2)/var_of_var + 4
  # that serves as prior degrees of freedom
  
  b_hyper <- (var_pop*(a_hyper - 2))/a_hyper
  
  
  #browser()
  if(BDA_model){ # Not yet working well; inverse chi-squared model from BDA3
    avg_df_fn_bayes <- avg_df_fn_bayes %>% dplyr::group_by(Gene_ID, Condition) %>%
      dplyr::mutate(sd_post = (sd_log_kd*nreps[Condition] + a_hyper*exp(lm_list[[Condition]][1] + lm_list[[Condition]][2]*log10(nreads)))/(a_hyper + nreps[Condition] - 2) ) %>%
      dplyr::mutate(log_kd_post = (avg_log_kd*(nreps[Condition]*(1/(sd_post^2))))/(nreps[Condition]/(sd_post^2) + (1/sdp^2)) + (theta_o*(1/sdp^2))/(nreps[Condition]/(sd_post^2) + (1/sdp^2))) %>%
      dplyr::mutate(kdeg = exp(log_kd_post) ) %>%
      dplyr::mutate(kdeg_sd = sd_post*exp(log_kd_post) ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Gene_ID) %>%
      dplyr::mutate(effect_size = log_kd_post - log_kd_post[Condition == 1]) %>%
      dplyr::mutate(effect_std_error = ifelse(Condition == 1, sd_post/sqrt(nreps[1]), sqrt((sd_post[Condition == 1]^2)/nreps[1] + (sd_post^2)/nreps[Condition] ) )) %>%
      dplyr::mutate(L2FC_kdeg = effect_size*log2(exp(1))) %>%
      dplyr::mutate(pval = pmin(1, 2*stats::pnorm((abs(effect_size) - null_cutoff)/effect_std_error, lower.tail = FALSE))) %>%
      dplyr::ungroup()
    
  }else{
    # Regularize estimates with Bayesian models and empirically informed priors
    avg_df_fn_bayes <- avg_df_fn_bayes %>% dplyr::group_by(Gene_ID, Condition) %>%
      dplyr::mutate(sd_post = exp( (log(sd_log_kd)*nreps[Condition]*(1/true_vars$true_var[Condition]) + (1/lm_var[[Condition]])*(lm_list[[Condition]][1] + lm_list[[Condition]][2]*log10(nreads)))/(nreps[Condition]*(1/true_vars$true_var[Condition]) + (1/lm_var[[Condition]])) )) %>%
      dplyr::mutate(log_kd_post = (avg_log_kd*(nreps[Condition]*(1/(sd_post^2))))/(nreps[Condition]/(sd_post^2) + (1/sdp^2)) + (theta_o*(1/sdp^2))/(nreps[Condition]/(sd_post^2) + (1/sdp^2))) %>%
      dplyr::mutate(kdeg = exp(log_kd_post) ) %>%
      dplyr::mutate(kdeg_sd = sd_post*exp(log_kd_post) ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(Gene_ID) %>%
      dplyr::mutate(effect_size = log_kd_post - log_kd_post[Condition == 1]) %>%
      dplyr::mutate(effect_std_error = ifelse(Condition == 1, sd_post, sqrt(sd_post[Condition == 1]^2 + sd_post^2))) %>%
      dplyr::mutate(L2FC_kdeg = effect_size*log2(exp(1)))
    
    if(ztest){
      avg_df_fn_bayes <- avg_df_fn_bayes %>%
        dplyr::mutate(pval = pmin(1, 2*stats::pnorm((abs(effect_size) - null_cutoff)/effect_std_error, lower.tail = FALSE))) %>%
        dplyr::ungroup()
      
    }else{
      
      avg_df_fn_bayes <- avg_df_fn_bayes %>%
        dplyr::mutate(pval = pmin(1, 2*stats::pt((abs(effect_size) - null_cutoff)/effect_std_error, df = 2*(nreps[Condition]-1) + 2*a_hyper, lower.tail = FALSE))) %>%
        dplyr::ungroup()
        
      
    }
    
    
  }
  
  
  
  # STATISTICAL TESTING --------------------------------------------------------
  
  
  message("Assessing statistical significance")
  
  ## Populate various data frames with important fit information
  
  effects <- avg_df_fn_bayes$effect_size[avg_df_fn_bayes$Condition > 1]
  ses <- avg_df_fn_bayes$effect_std_error[avg_df_fn_bayes$Condition > 1]
  
  pval <- avg_df_fn_bayes$pval[avg_df_fn_bayes$Condition > 1]
  padj <- stats::p.adjust(pval, method = "BH")
  
  
  # ORGANIZE OUTPUT ------------------------------------------------------------
  
  Genes_effects <- avg_df_fn_bayes$Gene_ID[avg_df_fn_bayes$Condition > 1]
  Condition_effects <- avg_df_fn_bayes$Condition[avg_df_fn_bayes$Condition > 1]
  
  L2FC_kdegs <- avg_df_fn_bayes$L2FC_kdeg[avg_df_fn_bayes$Condition > 1]
  
  Effect_sizes_df <- data.frame(Genes_effects, Condition_effects, L2FC_kdegs, effects, ses, pval, padj)
  
  colnames(Effect_sizes_df) <- c("Feature_ID", "Exp_ID", "L2FC_kdeg", "effect", "se", "pval", "padj")
  
  # Add sample information to output
  df_fn <- merge(df_fn, sample_lookup, by.x = c("Condition", "Replicate"), by.y = c("mut", "reps"))
  
  # Add feature name information to output
  df_fn <- merge(df_fn, feature_lookup, by.x = "Gene_ID", by.y = "fnum")
  avg_df_fn_bayes <- merge(avg_df_fn_bayes, feature_lookup, by.x = "Gene_ID", by.y = "fnum")
  Effect_sizes_df <- merge(Effect_sizes_df, feature_lookup, by.x = "Feature_ID", by.y = "fnum")
  
  # Order output
  df_fn <- df_fn[order(df_fn$Gene_ID, df_fn$Condition, df_fn$Replicate),]
  avg_df_fn_bayes <- avg_df_fn_bayes[order(avg_df_fn_bayes$Gene_ID, avg_df_fn_bayes$Condition),]
  Effect_sizes_df <- Effect_sizes_df[order(Effect_sizes_df$Feature_ID, Effect_sizes_df$Exp_ID),]
  
  
  avg_df_fn_bayes <- avg_df_fn_bayes[,c("Gene_ID", "Condition", "avg_log_kd", "sd_log_kd", "nreads", "sdp", "theta_o", "sd_post",
                                        "log_kd_post", "kdeg", "kdeg_sd", "XF")]
  
  colnames(avg_df_fn_bayes) <- c("Feature_ID", "Exp_ID", "avg_log_kdeg", "sd_log_kdeg", "nreads", "sdp", "theta_o", "sd_post",
                                 "log_kdeg_post", "kdeg", "kdeg_sd", "XF")
  
  if(NSS){
    colnames(df_fn) <- c("Feature_ID", "Exp_ID", "Replicate", "logit_fn", "logit_fn_se", "nreads", "log_kdeg", "fn", "log_kd_se", "sample", "XF")
  }else{
    colnames(df_fn) <- c("Feature_ID", "Exp_ID", "Replicate", "logit_fn", "logit_fn_se", "nreads", "log_kdeg", "kdeg", "log_kd_se", "sample", "XF")
  }
  
  if(is.null(Mutrates)){
    # Convert to tibbles because I like tibbles better
    fast_list <- list(dplyr::as_tibble(df_fn), dplyr::as_tibble(avg_df_fn_bayes), dplyr::as_tibble(Effect_sizes_df), c(a = a_hyper, b = b_hyper), lm_list)
    
    names(fast_list) <- c("Fn_Estimates", "Regularized_ests", "Effects_df", "Hyper_Parameters", "Mean_Variance_lms")
    
    class(fast_list) <- "FastFit"
    
    message("All done! Run QC_checks() on your bakRFit object to assess the quality of your data and get recommendations for next steps.")
    
  }else{
    # Convert to tibbles because I like tibbles better
    fast_list <- list(dplyr::as_tibble(df_fn), dplyr::as_tibble(avg_df_fn_bayes), dplyr::as_tibble(Effect_sizes_df), Mutrates, c(a = a_hyper, b = b_hyper), lm_list)
    
    names(fast_list) <- c("Fn_Estimates", "Regularized_ests", "Effects_df", "Mut_rates", "Hyper_Parameters", "Mean_Variance_lms")
    
    class(fast_list) <- "FastFit"
    
    message("All done! Run QC_checks() on your bakRFit object to assess the quality of your data and get recommendations for next steps.")
    
  }
  
  return(fast_list)
  
}
