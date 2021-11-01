
#' Extract data for efficient analysis from cB
#'
#' This function processes the cB into a form necessary for the parameter estimation function
#' that does not call Stan and is thus much more efficient and scalable.
#' @param cB_raw cB file generated from TL-seq pipeline
#' @param samp_list vector of names of control samples
#' @param type_list vector with 1 entry per sample; 0 = no s4U, 1 = s4U fed
#' @param mut_list vector with 1 entry per sample; 1 = reference, > 1 = different experimental conditions (e.g., KO of protein X)
#' @param rep_list vector with 1 entry per sample that indexes replicate; 1 = 1st replicate, 2 = 2nd replicate, etc.
#' @param tl single numerical value; s4U label time used in s4U fed samples
#' @param keep_input two element vector; 1st element is highest mut rate accepted in control samples, 2nd element is read count cutoff
#' @param FOI Features of interest; character vector containing names of features to analyze
#' @param concat Boolean; If TRUE, FOI is concatenated with output of reliableFeatures
#' @importFrom magrittr %>%
#' @return returns dataframe that can be passed to fast analysis
#' @export
cBtofast <- function(cB_raw,
                     samp_list,
                     type_list,
                     mut_list,
                     rep_list,
                     tl,
                     keep_input=c(0.2, 50),
                     FOI = c(),
                     concat = TRUE){

  cB <- cB_raw %>%
    dplyr::select(sample, XF, TC, n, nT)

  c_list <- samp_list[type_list == 0]

  names(type_list) <- samp_list
  names(mut_list) <- samp_list
  names(rep_list) <- samp_list

  nreps <- max(rep_list)

  # Helper function:
  getType <- function(s) type_list[paste(s)]
  getMut <- function(s) mut_list[paste(s)]
  getRep <- function(s) rep_list[paste(s)]

  # Get reliable features:
  if(concat == TRUE | is.null(FOI)){
    reliables <- DynamicSeq::reliableFeatures(cB = cB, c_list = c_list, high_p = keep_input[1], totcut = keep_input[2])
    keep <- c(FOI, reliables[!(reliables %in% FOI)])
  }else{
    keep <- FOI
  }
  # Get only the desired features:

  ranked_features_df  <- cB %>%
    dplyr::ungroup() %>%
    dplyr::filter(XF %in% keep) %>%
    dplyr::group_by(XF) %>%
    dplyr::summarize(n = sum(n)) %>%
    dplyr::mutate(fnum = order(-n)) %>%
    dplyr::arrange(fnum) %>%
    dplyr::select(XF, fnum)


  sdf <- cB %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample, XF, TC, nT) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::right_join(ranked_features_df, by = 'XF')

  d = sdf
  slist = samp_list
  tlist = type_list
  mlist = mut_list
  rlist = rep_list
  kp = keep


  df <- d %>%
    dplyr::ungroup() %>%
    dplyr::filter(sample %in% slist)

  df$type <- paste(df$sample) %>% purrr::map_dbl(function(x) getType(x))
  df$type <- as.integer(df$type)

  df$mut <- paste(df$sample) %>% purrr::map_dbl(function(x) getMut(x))
  df$mut <- as.integer(df$mut)

  df$reps <- paste(df$sample) %>% purrr::map_dbl(function(x) getRep(x))
  df$reps <- as.integer(df$reps)


  return(df)
}

#' Efficiently analyze nucleotide recoding data
#'
#' \code{fast_analysis} analyzes nucleotide recoding data using either maximum likelihood estimation with the L-BFGS-B algorithm
#' implemented by \code{stats::optim} or a more conservative Bayesian hypothesis testing strategy. Output includes fraction new
#' estimates for individual estimates as well as fraction news and effect sizes (L2FC(kdeg)s and changes in logit(fraction new))
#' averaged across replicate data. Averaging takes into account uncertainties estimated using the Fisher Information and estimates
#' are regularized using analytical results of fully Bayesian models. The result is that fraction news are shrunk towards population means
#' and that uncertainties are shrunk towards a mean-variance trend estimated as part of the analysis.
#'
#' Unless the user supplies estimates for pnew and pold, the first step of \code{fast_analysis} is to estimate the background
#' and s4U induced mutation rates. The former is best performed with a -s4U control sample, that is, a normal RNA-seq sample
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
#'
#' Once mutation rates are estimated, fraction news for each feature in each sample are estimated. The default approach utilized is MLE
#' using the L-BFGS-B algorithm implemented in \code{stats::optim}. The assumed likelihood function is derived from a Poisson mixture
#' model with rates adjusted according to each feature's empirical U-content (the average number of Us present in sequencing reads mapping
#' to that feature in a particular sample). An alternative model is a more conservative Bayesian hypothesis testing model with an uninformative
#' labeled vs. unlabeled prior. The Bayesian hypothesis testing model estimates the posterior probability that a set of reads with equivalent mutation
#' and U content are labeled using a binomial likelihood and a prior of 1/2 labeled and 1/2 unlabeled. For example, if the background
#' mutation rate is 0.001 and the s4U induced mutation rate is 0.1, then reads with 1 mutation in 25 Us have a posterior probability of being
#' labeled of 0.891. If there are 15 such reads, then 15 x 0.891 (or 13.365) of those reads are called labeled and 15 x (1 - 0.891)
#' (or 1.635) are called unlabeled. The estimated fraction new is then the number of reads called labeled divided by the total number
#' of reads. This approach biases fraction new estimates towards 0.5 and is thus a more conservative estimate of the fraction new
#' in each sample, as it is very skeptical of extreme fraction new values. This approach avoids some of the numerical instabilities
#' of MLE strategies and can thus be favorable if running into problems performing the MLE fit.
#'
#' Once fraction news are estimated, the uncertainty in the fraction new is estimated using the Fisher Information. In the limit of
#' large datasets, the variance of the MLE is inversely related to the Fisher Information evaluated at the MLE. Mixture models are
#' typically singular, meaning that the Fisher information matrix is not positive definite and asymptotic results for the variance
#' do not necessarily hold. This can be understood as arising due to nothing stopping the model from estimating a background mutation
#' rate of 0. As the mutation rates are estimated a priori and fixed to be > 0, these problems are eliminated. In addition, when assessing
#' the uncertainty of replicate fraction new estimates, the size of the dataset is the raw number of sequencing reads that map to a
#' particular feature. This number is often large (>100) which increases the validity of invoking asymptotics.
#'
#' With fraction news and their uncertainties estimated, replicate estimates are pooled and regularized. There are two key steps in this
#' downstream analysis. 1st, the uncertainty for each feature is used to extrapolate a linear log(uncertainty) vs. log10(read depth) trend,
#' and uncertainties for individual features are shrunk towards the regression line. The uncertainty for each feature is a combination of the
#' Fisher Information asymptotic uncertainty as well as the amount of variability seen between estimates. Regularization of uncertainty
#' estimates is performed using the analytic results of a Normal distribution likelihood with known mean and unknown variance and conjugate
#' priors. The prior parameters are estimated from the regression and amount of variability about the regression line. The strength of
#' regularization can be tuned by adjusting the \code{prior_weight} parameter, with larger numbers yielding stronger shrinkage towards
#' the regression line. The 2nd step is to regularize the average fraction new estimates. This is done using the analytic results of a
#' Normal distribution likelihood model with unknown mean and known variance and conjugate priors. The prior parameters are estimated from the
#' population wide fraction new distribution (using its mean and standard deviation as the mean and standard deviation of the normal prior).
#' In the 1st step, the known mean is assumed to be the average fraction new, averaged across replicates and weighted by the number of reads
#' mapping to the feature in each replicate. In the 2nd step, the known variance is assumed to be that obtained following regularization
#' of the uncertainty estimates
#'
#' Effect sizes (changes in fraction new) are obtained as the difference in logit(fraction new) means between the reference and experimental
#' sample(s), and the logit(fraction new)s are assumed to be independent so that the variance of the effect size is the sum of the
#' logit(fraction new) variances. P-values assessing the significance of the effect size are obtained using a moderated t-test with number
#' of degrees of freedom determined from the uncertainty regression hyperparameters and are adjusted for multiple testing using the Benjamini-
#' Hochberg procedure to control false discovery rates (FDRs).
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
#' @param lower Lower bound for MLE with L-BFGS-B algorithm
#' @param upper Upper bound for MLE with L-BFGS-B algorithm
#' @param mut_reg If MLE has instabilities, empircal mut rate will be used to estimate fn, multiplying pnew by 1+mut_reg and pold by 1-mut_reg to regularize fn
#' @param p_mean Mean of normal distribution used as prior penalty in MLE of logit(fn)
#' @param p_sd Standard deviation of normal distribution used as prior peanlty in MLE of logit(fn)
#' @return List with dataframes providing information about replicate-specific and pooled analysis results. The output includes:
#' \itemize{
#'  \item Fn_Estimates; dataframe with estimates for the fraction new and fraction new uncertainty for each feature in each replicate.
#'  The columns of this dataframe are:
#'  \itemize{
#'   \item Gene_ID; Numerical ID of feature
#'   \item Condition; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item Replicate; Numerical ID for replicate
#'   \item logit_fn; logit(fraction new) estimate, unregularized
#'   \item logit_fn_se; logit(fraction new) uncertainty, unregularized and obtained from Fisher Information
#'   \item fn_estimate; fraction new estimate (inverse logit of logit_fn)
#'   \item nreads; Number of reads mapping to the feature in the sample for which the estimates were obtained
#'   \item sample; Sample name
#'   \item XF; Original feature name
#'  }
#'  \item Regularized_ests; dataframe with average fraction new and kdeg estimates, averaged across the replicates and regularized
#'  using priors informed by the entire dataset. The columns of this dataframe are:
#'  \itemize{
#'   \item Gene_ID; Numerical ID of feature
#'   \item Condition; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item avg_logit_fn; Weighted average of logit(fn) from each replicate, weighted by sample and feature-specific read depth
#'   \item sd_logit_fn; Variability of the logit(fn) estimates, estimated using minimum MSE variance estimator
#'   \item nreads; Total number of reads mapping to the feature in that condition
#'   \item sdp; Prior standard deviation for fraction new estimate regularization
#'   \item theta_o; Prior mean for fraction new estimate regularization
#'   \item sd_post; Posterior uncertainty
#'   \item logit_fn_post; Posterior logit(fraction new)
#'   \item kdeg; average kdeg estimate
#'   \item kdeg_sd; kdeg uncertainty
#'   \item XF; Original feature name
#'  }
#'  \item Effects_df; dataframe with estimates of the effect size (change in logit(fn)) comparing each experimental condition to the
#'  reference sample for each feature. This dataframe also includes p-values obtained from a moderated t-test. The columns of this
#'  dataframe are:
#'  \itemize{
#'   \item Genes_effects; Numerical ID of feature
#'   \item Condition_effects; Numerical ID for experimental condition (Exp_ID from metadf)
#'   \item L2FC_kdeg; L2FC(kdeg) estimate
#'   \item effects; Change in logit(fraction) new comparing reference and experimental sample(s)
#'   \item ses; Uncertainty in effect size
#'   \item pval; P-value obtained using effect_size, effect_std_error, and a moderated t-test
#'   \item padj; pval adjusted for multiple testing using Benjamini-Hochberg procedure
#'   \item XF; Original feature name
#'  }
#'  \item Mut_rates; list of two elements. The 1st element is a dataframe of s4U induced mutation rate estimates, where the mut column
#'  represents the experimental ID and the rep column represents the replicate ID. The 2nd element is the single background mutation
#'  rate estimate used
#'  \item Hyper_Parameters; vector of two elements, named a and b. These are the hyperparameters estimated from the uncertainties for each
#'  feature, and represent the two parameters of a Scaled Inverse Chi-Square distribution. Importantly, a is the number of additional
#'  degrees of freedom provided by the sharing of uncertainty information across the dataset, to be used in the moderated t-test.
#'  \item Mean_Variance_lm; a linear model object obtained from the uncertainty vs. read count regression model.
#' }
#' @importFrom magrittr %>%
#' @export
fast_analysis <- function(df, pnew = NULL, pold = NULL, no_ctl = FALSE,
                          read_cut = 50, features_cut = 10,
                          nbin = NULL, prior_weight = 2,
                          MLE = TRUE,
                          lower = -7, upper = 7,
                          se_max = 2.5,
                          mut_reg = 0.1,
                          p_mean = 0,
                          p_sd = 1){

  ## Check pold
  if(length(pold) > 1){
    stop("pold must be NULL or length 1")
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

  nMT <- max(df$mut)
  nreps <- max(df$reps)

  ## Check pnew
  if(!is.null(pnew)){
    if(!all(is.numeric(pnew))){
      stop("All elements of pnew must be numeric")
    }else if(length(pnew) != nMT*nreps){
      stop("pnew must be a vector of length == number of s4U fed samples in dataset")
    }else if(!all(pnew > 0)){
      stop("All elements of pnew must be > 0")
    }else if(!all(pnew <=1)){
      stop("All elements of pnew must be <= 1")
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

  ## Check read_cut
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

  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))

  #Old mutation rate estimation

  #Trim df and name columns
  if(is.null(pnew)){
    message("Estimating labeled mutation rate")
    Mut_data <- df

    ##New Mutation rate Estimation
    # Extract only s4U labeled data to estimate s4U mut rate
    New_data <- Mut_data[Mut_data$type == 1, ]

    # Calculate avg. mut rate at each row
    New_data$avg_mut <- New_data$TC/New_data$nT

    # Remove rows with NAs
    New_data <- New_data[!is.na(New_data$avg_mut),]

    # calculate total number of mutations
    # which is the avg. for that row of dataframe times n
    # the number of reads that had the identical average
    New_data$weight_mut <- New_data$avg_mut*New_data$n

    # This is to estimate the total mutation rate for each gene in
    # each replicate and each experimental condition
    New_data_summary <- New_data %>%
      dplyr::group_by(reps, mut, fnum) %>% # group by gene, replicate ID, and experiment ID
      dplyr::summarise(avg_mut = sum(weight_mut)/sum(n), n = sum(n))

    # Order datalist so that it's ordered by sample and then avg mutation rate
    # Goal is to use the highest avg. mutation rates to estimate s4U mutation rate,
    # assuming that highest mutation rates are from fast turnover, compeltely
    # labeled transcripts
    New_data_ordered <- New_data_summary[order(New_data_summary$mut, New_data_summary$reps, New_data_summary$avg_mut, decreasing=TRUE), ]

    ## This part has some magic numbers I should get rid of
    ## or move to user input

    New_data_cutoff <- New_data_ordered[New_data_ordered$n > read_cut,]

    # Check to make sure that the number of features that made it past the
    # read count filter is still more than the total number of features required for
    # mutation rate estimate
    check <- New_data_cutoff %>% dplyr::ungroup() %>%
      dplyr::count(mut, reps, sort = TRUE)

    if(sum(check$n < features_cut) > 0){
      stop("Not enough features made it past the read cutoff filter in one sample; try decreasing read_cut or features_cut")
    }else{
      New_data_estimate <- New_data_cutoff %>% dplyr::group_by(mut, reps) %>%
        dplyr::summarise(pnew = mean(avg_mut[1:features_cut]))
      message(paste0(c("Estimated pnews for each sample are:", capture.output(New_data_estimate)), collapse = "\n"))
    }


  }else{ # Need to construct pmut dataframe from User input
    nMT <- max(df$mut)
    nreps <- max(df$reps)
    if(length(pnew) == 1){

      pnew_vect <- rep(pnew, times = (nMT*nreps))

      rep_vect <- rep(seq(from = 1, to = nreps), times = nMT)

      mut_vect <- rep(seq(from = 1, to = nMT), each = nreps)

      New_data_estimate <- data.frame(mut_vect, rep_vect, pnew_vect)
      colnames(New_data_estimate) <- c("mut", "reps", "pnew")

    } else if( length(pnew) != (nMT*nreps) ){
      stop("User inputted pnew is not of length 1 or of length equal to number of samples")
    } else{
      rep_vect <- rep(seq(from = 1, to = nreps), times = nMT)

      mut_vect <- rep(seq(from = 1, to = nMT), each = nreps)
      New_data_estimate <- data.frame(mut_vect, rep_vect, pnew)
      colnames(New_data_estimate) <- c("mut", "reps", "pnew")
    }
  }

  if(is.null(pold)){

    if((sum(df$type == 0) == 0) | (no_ctl)){ # Estimate using low mutation rate features
      message("Estimating unlabeled mutation rate")

      New_data <- df

      # Calculate avg. mut rate at each row
      New_data$avg_mut <- New_data$TC/New_data$nT

      # Remove rows with NAs
      New_data <- New_data[!is.na(New_data$avg_mut),]

      # calculate total number of mutations
      # which is the avg. for that row of dataframe times n
      # the number of reads that had the identical average
      New_data$weight_mut <- New_data$avg_mut*New_data$n

      # This is to estimate the total mutation rate for each gene in
      # each replicate and each experimental condition
      New_data_summary <- New_data %>%
        dplyr::group_by(fnum) %>% # group by gene, replicate ID, and experiment ID
        dplyr::summarise(avg_mut = sum(weight_mut)/sum(n), n = sum(n), .groups = "keep")

      # Order datalist so that it's ordered by sample and then avg mutation rate
      # Goal is to use the lowest avg. mutation rates to estimate s4U mutation rate,
      # assuming that highest mutation rates are from fast turnover, completely
      # labeled transcripts
      New_data_ordered <- New_data_summary[order(New_data_summary$avg_mut, decreasing=FALSE), ]


      # Filter out for high read depth features
      New_data_cutoff <- New_data_ordered[New_data_ordered$n > read_cut,]

      # Check to make sure that the number of features that made it past the
      # read count filter is still more than the total number of features required for
      # mutation rate estimate
      check <- nrow(New_data_cutoff)

      if(check < features_cut){
        stop("Not enough features made it past the read cutoff filter in one sample; try decreasing read_cut or features_cut")
      }else{
        pold <- stats::weighted.mean(New_data_cutoff$avg_mut[1:features_cut], w = New_data_cutoff$n[1:features_cut])
        message(paste(c("Estimated pold is: ", pold), collapse = " "))
      }


    }else{ # Estimate using -s4U data
      message("Estimating unlabeled mutation rate")

      #Old mutation rate estimation
      Mut_data <- df

      Old_data <- Mut_data[Mut_data$type == 0, ]

      Old_data$avg_mut <- Old_data$TC/Old_data$nT

      # Remove rows with NAs
      Old_data <- Old_data[!is.na(Old_data$avg_mut),]

      Old_data$weight_mut <- Old_data$avg_mut*Old_data$n

      #Old_data$n <- rep(1, times=nrow(Old_data))

      Old_data_summary <- Old_data %>%
        dplyr::group_by(reps, mut, fnum) %>% # group by gene, replicate ID, and experiment ID
        dplyr::summarise(avg_mut = sum(weight_mut)/sum(n), n = sum(n))

      # Order data differently than for s4U mut rate estimation
      # Difference is that every mutation is a background mutation in these samples
      # So we just want the highest confidence estimation, meaning we should only
      # order by read counts
      Old_data_ordered <- Old_data_summary[order(Old_data_summary$n, decreasing=TRUE), ]

      ## This part has some magic numbers I should get rid of
      ## or move to user input

      Old_data_cutoff <- Old_data_ordered[Old_data_ordered$n > read_cut,]

      # Check to make sure that the number of features that made it past the
      # read count filter is still more than the total number of features required for
      # mutation rate estimate
      check <- nrow(Old_data_cutoff)

      if(check < features_cut){
        stop("Not enough features made it past the read cutoff filter in one sample; try decreasing read_cut or features_cut")
      }else{
        pold <- stats::weighted.mean(Old_data_cutoff$avg_mut[1:features_cut], w = Old_data_cutoff$n[1:features_cut])
        message(paste(c("Estimated pold is: ", pold), collapse = " "))
      }
    }



  }

  if(!all(New_data_estimate$pnew - pold > 0)){
    stop("All pnew must be > pold; did you input an unusually large pold?")
  }

  pmuts_list <- list(New_data_estimate, pold)

  Mut_data <- df

  Mut_data <- Mut_data[Mut_data$type == 1,]

  ngene <- max(Mut_data$fnum)
  num_conds <- max(Mut_data$mut)
  nreps <- max(Mut_data$reps)

  sample_lookup <- Mut_data[, c("sample", "mut", "reps")] %>% dplyr::distinct()
  feature_lookup <- Mut_data[,c("fnum", "XF")] %>% dplyr::distinct()

  # Estimate fraction new in each replicate using binomial model
  message("Estimating fraction labeled")

  Mut_data <- dplyr::left_join(Mut_data, New_data_estimate, by = c("mut", "reps"))


  if(!MLE){
    # Bayesian Hypothesis Testing Method
    Mut_data_est <- Mut_data %>% dplyr::group_by(fnum, mut, reps, TC, nT) %>%
      # mutate(avg_mut = TC/nT) %>%
      # #mutate(prior_new = ifelse(avg_mut >= (pnew_est - 0.01), 0.99, (avg_mut + 0.01)/pnew_est )) %>%
      # mutate(prior_new = 0.9)%>%
      dplyr::mutate(New_prob = stats::dbinom(TC, size=nT, prob=pnew)) %>%
      dplyr::mutate(Old_prob = stats::dbinom(TC, size = nT, prob = pold)) %>%
      dplyr::mutate(News = n*(New_prob/(New_prob + Old_prob))) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fnum, mut, reps) %>%
      dplyr::summarise(nreads = sum(n), Fn_rep_est = sum(News)/nreads) %>%
      dplyr::mutate(logit_fn_rep = ifelse(Fn_rep_est == 1, logit(0.999), ifelse(Fn_rep_est == 0, logit(0.001), logit(Fn_rep_est)))) %>%
      dplyr::ungroup()
  }else{
    # MLE
    mixed_lik <- function(lam_n, lam_o, TC, n, logit_fn){
      logl <- sum(n*log(inv_logit(logit_fn)*(lam_n^TC)*exp(-lam_n) + (1-inv_logit(logit_fn))*(lam_o^TC)*exp(-lam_o) )) + log(dnorm(logit_fn, mean = p_mean, sd = p_sd))
      return(-logl)
    }

    Mut_data_est <- Mut_data %>% dplyr::ungroup() %>% dplyr::mutate(lam_n = pnew*nT, lam_o = pold*nT) %>%
      dplyr::group_by(fnum, mut, reps, TC) %>%
      dplyr::summarise(lam_n = sum(lam_n*n)/sum(n), lam_o = sum(lam_o*n)/sum(n),
                       n = sum(n), .groups = "keep") %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fnum, mut, reps) %>%
      dplyr::summarise(logit_fn_rep = optim(0, mixed_lik, TC = TC, n = n, lam_n = sum(lam_n*n)/sum(n), lam_o = sum(lam_o*n)/sum(n), method = "L-BFGS-B", lower = lower, upper = upper)$par, nreads =sum(n), .groups = "keep") %>%
      dplyr::mutate(logit_fn_rep = ifelse(logit_fn_rep == lower, runif(1, lower-0.2, lower), ifelse(logit_fn_rep == upper, runif(1, upper, upper+0.2), logit_fn_rep))) %>%
      dplyr::ungroup()

    ## Look for numerical instabilities
    instab_df <- Mut_data_est %>% dplyr::filter(abs(logit_fn_rep) > upper )

    fnum_instab <- instab_df$fnum
    mut_instab <- instab_df$mut
    reps_instab <- instab_df$reps

    instab_est <- dplyr::right_join(Mut_data, instab_df, by = c("fnum", "mut", "reps")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(totTC = TC*n, totU = nT*n) %>%
      dplyr::group_by(fnum, mut, reps) %>%
      dplyr::summarise(tot_mut = sum(totTC), totUs = sum(totU), pnew = mean(pnew), pold = mean(pold)) %>% dplyr::ungroup() %>%
      dplyr::mutate(avg_mut = tot_mut/totUs,
                    Fn_rep_est = (avg_mut - (1-mut_reg)*pold)/((1+mut_reg)*pnew - (1-mut_reg)*pold)) %>%
      dplyr::mutate(Fn_rep_est = ifelse(Fn_rep_est > 1, runif(1, min = inv_logit(upper-0.1), max = 1) , ifelse(Fn_rep_est < 0, runif(1, min = 0, max = inv_logit(lower+0.1)), Fn_rep_est ) )) %>%
      dplyr::mutate(logit_fn_rep = logit(Fn_rep_est)) %>%
      dplyr::select(fnum, mut, reps, logit_fn_rep)

    Mut_data_est <- dplyr::left_join(Mut_data_est, instab_est, by = c("fnum", "mut", "reps")) %>%
      dplyr::mutate(logit_fn_rep = ifelse(abs(logit_fn_rep.x) > upper, logit_fn_rep.y, logit_fn_rep.x)) %>%
      dplyr::select(fnum, mut, reps, nreads, logit_fn_rep) %>%
      dplyr::mutate(Fn_rep_est = inv_logit(logit_fn_rep)) %>%
      dplyr::mutate(kd_rep_est = -log(1 - Fn_rep_est)/tl[mut])

  }

  message("Estimating per replicate uncertainties")

  Mut_data <- dplyr::left_join(Mut_data, Mut_data_est[, c("kd_rep_est" ,"logit_fn_rep", "fnum", "mut", "reps")], by = c("fnum", "mut", "reps"))

  ## Estimate Fisher Info and uncertainties
  ## Could make more efficient by summarizing over nT info and using U_content to adjust pnew*avg_U
  Mut_data <- Mut_data %>% dplyr::ungroup() %>%
    dplyr::group_by(fnum, mut, reps, TC, pnew, logit_fn_rep, kd_rep_est) %>%
    dplyr::summarise(U_cont = sum(nT*n)/sum(n), n = sum(n), .groups = "keep") %>%
    dplyr::mutate(Fisher_lkd_num = tl[mut]*kd_rep_est*( ((pnew*U_cont)^TC)*exp(-pnew*U_cont) - ( ((pold*U_cont)^TC)*exp(-pold*U_cont) ) ) ) %>%
    dplyr::mutate(Fisher_lkd_den = (exp(kd_rep_est*tl[mut]) - 1)*((pnew*U_cont)^TC)*exp(-pnew*U_cont) + ((pold*U_cont)^TC)*exp(-pold*U_cont) ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(fnum, mut, reps) %>%
    dplyr::summarise(Fisher_kdeg = sum(n*((Fisher_lkd_num/Fisher_lkd_den)^2))/sum(n),tot_n = sum(n)) %>% #,
    #Fisher_fn = sum(n*((Fisher_fn_num/Fisher_fn_den)^2)), tot_n = sum(n)) %>%
    dplyr::mutate(log_kd_se = 1/sqrt(tot_n*Fisher_kdeg)) %>%
    dplyr::mutate(log_kd_se = ifelse(log_kd_se > se_max, se_max, log_kd_se))#, Fn_se = 1/sqrt(tot_n*Fisher_fn))


  Mut_data_est$log_kd_se <- Mut_data$log_kd_se

  Mut_data_est <- Mut_data_est %>% dplyr::mutate(log_kd_rep_est = log(kd_rep_est))

  # Mut_data_est$fn_se = Mut_data$Fn_se

  ## Now affiliate each fnum, mut with a bin Id based on read counts,
  ## bin data by bin_ID and average log10(reads) and log(sd(logit_fn))

  if(is.null(nbin)){
    nbin <- max(c(round(ngene*num_conds*nreps/100), 10))
  }

  message("Estimating read count-variance relationship")

  Binned_data <- Mut_data_est %>% dplyr::group_by(fnum, mut) %>%
    dplyr::summarise(nreads = sum(nreads), kd_sd_log = log(sqrt(1/sum(1/((sd(log_kd_rep_est)^2) + log_kd_se^2 ) ) ) )) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bin_ID = as.numeric(Hmisc::cut2(nreads, g = nbin))) %>% dplyr::group_by(bin_ID, mut) %>%
    dplyr::summarise(avg_reads = mean(log10(nreads)), avg_sd = mean(kd_sd_log))

  ## Regress avg_reads vs. avg_sd
  lm_list <- vector("list", length = nMT)
  lm_var <- lm_list

  for(i in 1:nMT){
    heterosked_lm <- stats::lm(avg_sd ~ avg_reads, data = Binned_data[Binned_data$mut == i,] )
    h_int <- summary(heterosked_lm)$coefficients[1,1]
    h_slope <- summary(heterosked_lm)$coefficients[2,1]
    lm_list[[i]] <- c(h_int, h_slope)

    lm_var[[i]] <- var(stats::residuals(heterosked_lm))
  }


  true_vars <-  Mut_data_est %>% dplyr::group_by(fnum, mut) %>%
    dplyr::summarise(nreads = sum(nreads), kd_sd_log = log(sqrt(1/sum(1/((sd(log_kd_rep_est)^2) + log_kd_se^2 ) ) ) )) %>%
    dplyr::ungroup() %>% dplyr::group_by(fnum, mut) %>% dplyr::mutate(slope = lm_list[[mut]][2], intercept = lm_list[[mut]][1]) %>%
    dplyr::group_by(mut) %>%
    dplyr::summarise(true_var = var(kd_sd_log - (intercept + slope*log10(nreads) ) ))

  log_kd <- as.vector(Mut_data_est$log_kd_rep_est)
  logit_fn <- as.vector(Mut_data_est$logit_fn_rep)

  if(!(all(logit_fn < upper) & all(logit_fn > lower))){
    num_unstable <- sum(logit_fn >= upper) + sum(logit_fn <= lower)
    tot_ests <- length(logit_fn)
    prcnt_unstable <- round((num_unstable/tot_ests)*100,3)
    warning(paste0(num_unstable, " out of " , tot_ests, " (", prcnt_unstable,"%)", " logit(fn) estimates are at the upper or lower bounds set. These likely represent features with limited data or extremely stable/unstable features. If the number of boundary estimates is concerning, try increasing the magnitude of upper and lower."))
  }

  kd_estimate <- as.vector(Mut_data_est$kd_rep_est)
  # fn_se <- as.vector(Mut_data_est$fn_se)
  log_kd_se <- as.vector(Mut_data_est$log_kd_se)
  Replicate <- as.vector(Mut_data_est$reps)
  Condition <- as.vector(Mut_data_est$mut)
  Gene_ID <- as.vector(Mut_data_est$fnum)
  nreads <- as.vector(Mut_data_est$nreads)

  rm(Mut_data_est)

  df_fn <- data.frame(logit_fn, Replicate, Condition, Gene_ID, nreads, log_kd, kd_estimate, log_kd_se)

  # Remove vectors no longer of use
  rm(logit_fn)
  rm(Replicate)
  rm(Condition)
  rm(Gene_ID)
  rm(nreads)

  df_fn <- df_fn[order(df_fn$Gene_ID, df_fn$Condition, df_fn$Replicate),]

  nreps <- max(df_fn$Replicate)

  message("Averaging replicate data and regularizing estimates")

  #Average over replicates and estimate hyperparameters
  avg_df_fn_bayes <- df_fn %>% dplyr::group_by(Gene_ID, Condition) %>%
    dplyr::summarize(avg_log_kd = stats::weighted.mean(log_kd, 1/log_kd_se),
                     sd_log_kd = sqrt(1/sum(1/((sd(log_kd)^2) + log_kd_se^2 ) ) ),
                     nreads = sum(nreads)) %>% dplyr::ungroup() %>%
    dplyr::group_by(Condition) %>%
    dplyr::mutate(sdp = sd(avg_log_kd)) %>%
    dplyr::mutate(theta_o = mean(avg_log_kd)) %>%
    # mutate(var_pop = mean(sd_logit_fn^2)) %>%
    # mutate(var_of_var = var(sd_logit_fn^2)) %>%
    # mutate(two_params = 8*(var_pop^4)/var_of_var) %>%
    # mutate(roots = RConics::cubic(c(1, -(4 + 2*(var_pop^2)/var_of_var), two_params, two_params))) %>%
    # mutate(a_hyper = roots[(roots > 2) & (!is.complex(roots))]) %>%
    # mutate(b_hyper = (var_pop*(a_hyper - 2))/a_hyper) %>%
    dplyr::ungroup()



  #Calcualte population averages
  # What I need to do is calculate these parameters for each Condition
  #sdp <- sd(avg_df_fn$avg_logit_fn) # Will be prior sd in regularization of mean
  #theta_o <- mean(avg_df_fn$avg_logit_fn) # Will be prior mean in regularization of mean

  var_pop <- mean(avg_df_fn_bayes$sd_log_kd^2) # Will be prior mean in regularization of sd
  var_of_var <- stats::var(avg_df_fn_bayes$sd_log_kd^2) # Will be prior variance in regularization of sd


  ## Regularize standard deviation estimate
  # Estimate hyperpriors with method of moments
  two_params <- 8*(var_pop^4)/var_of_var
  b <- c(1, -(4 + 2*(var_pop^2)/var_of_var), two_params, two_params)

  roots <- RConics::cubic(b)

  a_hyper <- roots[(roots > 2) & (!is.complex(roots))]
  # Divide this by 2 to get inverse-gamma hyperprior
  # that serves as prior degrees of freedom

  b_hyper <- (var_pop*(a_hyper - 2))/a_hyper


  # Regularize estimates with Bayesian models and empirically informed priors
  avg_df_fn_bayes <- avg_df_fn_bayes %>% dplyr::group_by(Gene_ID, Condition) %>%
    #dplyr::mutate(sd_post = sqrt((a_hyper*b_hyper + nreps*sd_logit_fn)/(a_hyper + nreps - 2))) %>%
    dplyr::mutate(sd_post = exp( (log(sd_log_kd)*nreps*(1/true_vars$true_var[Condition]) + (1/lm_var[[Condition]])*(lm_list[[Condition]][1] + lm_list[[Condition]][2]*log10(nreads)))/(nreps*(1/true_vars$true_var[Condition]) + (1/lm_var[[Condition]])) )) %>%
    dplyr::mutate(log_kd_post = (avg_log_kd*(nreps*(1/(sd_post^2))))/(nreps/(sd_post^2) + (1/sdp^2)) + (theta_o*(1/sdp^2))/(nreps/(sd_post^2) + (1/sdp^2))) %>%
    #dplyr::mutate(sd_post = sqrt(1/(nreps/(sd_post^2) + (1/sdp^2)))) %>%
    dplyr::mutate(kdeg = exp(log_kd_post) ) %>%
    dplyr::mutate(kdeg_sd = sd_post*exp(log_kd_post) ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Gene_ID) %>%
    dplyr::mutate(effect_size = log_kd_post - log_kd_post[Condition == 1]) %>%
    dplyr::mutate(effect_std_error = ifelse(Condition == 1, sd_post, sqrt(sd_post[Condition == 1]^2 + sd_post^2))) %>%
    dplyr::mutate(L2FC_kdeg = effect_size) %>%
    dplyr::mutate(pval = 2*pt(-abs(effect_size/effect_std_error), df = 2*nreps - 2 + 2*a_hyper)) %>%
    dplyr::ungroup()

  # Calcuate lfsr using ashr package

  message("Assessing statistical significance")

  effects <- avg_df_fn_bayes$effect_size[avg_df_fn_bayes$Condition > 1]
  ses <- avg_df_fn_bayes$effect_std_error[avg_df_fn_bayes$Condition > 1]

  pval <- avg_df_fn_bayes$pval[avg_df_fn_bayes$Condition > 1]
  padj <- stats::p.adjust(pval, method = "BH")

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


  #hyperpars <- c(sdp, theta_o, var_pop, var_of_var, a_hyper, b_hyper)
  #names(hyperpars) <- c("Mean Prior sd", "Mean prior mean", "Variance prior mean", "Variance prior variance", "Variance hyperparam a", "Variance hyperparam b")

  avg_df_fn_bayes <- avg_df_fn_bayes[,c("Gene_ID", "Condition", "avg_log_kd", "sd_log_kd", "nreads", "sdp", "theta_o", "sd_post",
                                        "log_kd_post", "kdeg", "kdeg_sd", "XF")]

  colnames(avg_df_fn_bayes) <- c("Feature_ID", "Exp_ID", "avg_log_kdeg", "sd_log_kdeg", "nreads", "sdp", "theta_o", "sd_post",
                                 "log_kdeg_post", "kdeg", "kdeg_sd", "XF")

  colnames(df_fn) <- c("Feature_ID", "Exp_ID", "Replicate", "logit_fn", "nreads", "log_kdeg", "kdeg", "log_kd_se", "sample", "XF")

  #fast_list <- list(estimate_df, avg_df_fn_bayes, Effect_sizes_df, pmuts_list, hyperpars)
  fast_list <- list(df_fn, avg_df_fn_bayes, Effect_sizes_df, pmuts_list, c(a = a_hyper, b = b_hyper), lm_list)

  names(fast_list) <- c("Fn_Estimates", "Regularized_ests", "Effects_df", "Mut_rates", "Hyper_Parameters", "Mean_Variance_lms")

  class(fast_list) <- "FastFit"


}
