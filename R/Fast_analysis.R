
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

#' This function efficiently analyzes data to provide fraction new estimates and uncertainties
#'
#' @param df Dataframe in form provided by cB_to_Fast
#' @param pnew Labeled read mutation rate; default of 0 means that model estimates rate from s4U fed data
#' @param pold Unlabeled read mutation rate; default of 0 means that model estimates rate from no-s4U fed data
#' @param read_cut Minimum number of reads for a given feature-sample combo to be used for mut rate estimates
#' @param features_cut Number of features to estimate sample specific mutation rate with
#' @param prior_weight Determines extent to which logit(fn) variance is regularized to the mean-variance regression line
#' @return list with dataframe of replicate specific estimates as well as dataframe of pooled estimates
#' @importFrom magrittr %>%
#' @export
fast_analysis <- function(df, pnew = NULL, pold = NULL, read_cut = 50, features_cut = 10, nbin = NULL, prior_weight = 2){

  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))

  #Old mutation rate estimation

  #Trim df and name columns
  if(is.null(pnew)){
    message("Estimating labeled mutation rate")
    Mut_data <- df[,1:9]
    colnames(Mut_data) <- c("sample", "XF", "TC", "nT", "n", "fnum", "type", "mut", "reps")

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
      message(paste(c("Estimated pnews are: ", New_data_estimate$pnew), collapse = " "))
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
    message("Estimating unlabeled mutation rate")

    #Old mutation rate estimation
    Mut_data <- df[,1:9]
    colnames(Mut_data) <- c("sample", "XF", "TC", "nT", "n", "fnum", "type", "mut", "reps")

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
      pold <- mean(Old_data_cutoff$avg_mut[1:features_cut])
      message(paste(c("Estimated pold is: ", pold), collapse = " "))
    }

  }


  pmuts_list <- list(New_data_estimate, pold)

  Mut_data <- df

  Mut_data <- Mut_data[Mut_data$type == 1,]

  ngene <- max(Mut_data$fnum)
  num_conds <- max(Mut_data$mut)
  nreps <- max(Mut_data$reps)

  fn_rep_est <- rep(0, times=ngene*num_conds*nreps)
  dim(fn_rep_est) <- c(ngene, num_conds, nreps)

  R_ID <- fn_rep_est
  MT_ID <- R_ID
  FN_ID <- R_ID

  # Estimate fraction new in each replicate using binomial model
  message("Estimating fraction labeled")

  Mut_data <- merge(Mut_data, New_data_estimate, by = c("mut", "reps"))

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

  message("Estimating per replicate uncertainties")

  Mut_data <- merge(Mut_data, Mut_data_est[, c("logit_fn_rep", "fnum", "mut", "reps")], by = c("fnum", "mut", "reps"))

  ## Estimate Fisher Info and uncertainties
  ## Could make more efficient by summarizing over nT info and using U_content to adjust pnew*avg_U
  Mut_data <- Mut_data %>% ungroup() %>%
    dplyr::group_by(fnum, mut, reps, TC, pnew, logit_fn_rep) %>%
    dplyr::summarise(U_cont = sum(nT*n)/sum(n), n = sum(n), .groups = "keep") %>%
    dplyr::mutate(Exp_l_fn = exp(logit_fn_rep)) %>%
    # dplyr::mutate(fn = inv_logit(Mut_data_est$logit_fn_rep[(Mut_data_est$mut == mut) & (Mut_data_est$reps == reps) & (Mut_data_est$fnum == fnum)])) %>%
    # dplyr::mutate(Fisher_fn_num = ((pnew_est*nT)^TC)*exp(-pnew_est*nT) - ((pold*nT)^TC)*exp(-pold*nT) ) %>%
    # dplyr::mutate(Fisher_fn_den = fn*Fisher_fn_num + ((pold*nT)^TC)*exp(-pold*nT)) %>%
    dplyr::mutate(Inv_Fisher_Logit_3 = 1/(((pnew/pold)^TC)*exp(-(U_cont)*(pnew - pold)) - 1 )) %>%
    dplyr::mutate(Inv_Fisher_Logit_1 = 1 + Exp_l_fn ) %>%
    dplyr::mutate(Inv_Fisher_Logit_2 = ((1 + Exp_l_fn)^2)/Exp_l_fn) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(fnum, mut, reps) %>%
    dplyr::summarise(Fisher_Logit = sum(n/((Inv_Fisher_Logit_1 + Inv_Fisher_Logit_2*Inv_Fisher_Logit_3)^2))/sum(n), tot_n = sum(n)) %>% #,
    #Fisher_fn = sum(n*((Fisher_fn_num/Fisher_fn_den)^2)), tot_n = sum(n)) %>%
    dplyr::mutate(Logit_fn_se = 1/sqrt(tot_n*Fisher_Logit)) #, Fn_se = 1/sqrt(tot_n*Fisher_fn))

  Mut_data_est$logit_fn_se = Mut_data$Logit_fn_se
  # Mut_data_est$fn_se = Mut_data$Fn_se

  ## Now affiliate each fnum, mut with a bin Id based on read counts,
  ## bin data by bin_ID and average log10(reads) and log(sd(logit_fn))

  if(is.null(nbin)){
    nbin <- max(c(round(ngene*num_conds*nreps/100), 10))
  }

  message("Estimating mean-variance relationship")

  Binned_data <- Mut_data_est %>% dplyr::group_by(fnum, mut) %>%
    dplyr::summarise(nreads = sum(nreads),fn_sd_log = log(sqrt(1/sum(1/(sd(logit_fn_rep)^2 + logit_fn_se^2 ) ) ) )) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bin_ID = as.numeric(Hmisc::cut2(nreads, g = nbin))) %>% dplyr::group_by(bin_ID) %>%
    dplyr::summarise(avg_reads = mean(log10(nreads)), avg_sd = mean(fn_sd_log))

  ## Regress avg_reads vs. avg_sd
  heterosked_lm <- stats::lm(avg_sd ~ avg_reads, data = Binned_data )

  h_int <- summary(heterosked_lm)$coefficients[1,1]
  h_slope <- summary(heterosked_lm)$coefficients[2,1]


  logit_fn <- as.vector(Mut_data_est$logit_fn_rep)
  fn_estimate <- as.vector(Mut_data_est$Fn_rep_est)
  # fn_se <- as.vector(Mut_data_est$fn_se)
  logit_fn_se <- as.vector(Mut_data_est$logit_fn_se)
  Replicate <- as.vector(Mut_data_est$reps)
  Condition <- as.vector(Mut_data_est$mut)
  Gene_ID <- as.vector(Mut_data_est$fnum)
  nreads <- as.vector(Mut_data_est$nreads)

  rm(Mut_data_est)

  df_fn <- data.frame(logit_fn, logit_fn_se, fn_estimate, Replicate, Condition, Gene_ID, nreads)

  df_fn <- df_fn[order(df_fn$Gene_ID, df_fn$Condition, df_fn$Replicate),]

  nreps <- max(df_fn$Replicate)

  message("Averaging replicate data and regularizing estimates")

  #Average over replicates and estimate hyperparameters
  avg_df_fn_bayes <- df_fn %>% dplyr::group_by(Gene_ID, Condition) %>%
    dplyr::summarize(avg_logit_fn = stats::weighted.mean(logit_fn, 1/logit_fn_se),
              sd_logit_fn = sqrt(1/sum(1/(sd(logit_fn)^2 + logit_fn_se^2 ) ) ),
              nreads = sum(nreads)) %>% dplyr::ungroup() %>%
    dplyr::group_by(Condition) %>%
    dplyr::mutate(sdp = sd(avg_logit_fn)) %>%
    dplyr::mutate(theta_o = mean(avg_logit_fn)) %>%
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

  var_pop <- mean(avg_df_fn_bayes$sd_logit_fn^2) # Will be prior mean in regularization of sd
  var_of_var <- stats::var(avg_df_fn_bayes$sd_logit_fn^2) # Will be prior variance in regularization of sd


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
    dplyr::mutate(sd_post = exp( (log(sd_logit_fn)*nreps + prior_weight*(h_int + h_slope*log10(nreads)))/(nreps + prior_weight) )) %>%
    dplyr::mutate(logit_fn_post = (avg_logit_fn*(nreps*(1/(sd_post^2))))/(nreps/(sd_post^2) + (1/sdp^2)) + (theta_o*(1/sdp^2))/(nreps/(sd_post^2) + (1/sdp^2))) %>%
    dplyr::mutate(sd_post = sqrt(1/(nreps/(sd_post^2) + (1/sdp^2)))) %>%
    dplyr::mutate(kdeg = -log(1 - inv_logit(logit_fn_post))) %>%
    dplyr::mutate(kdeg_sd = sd_post*(exp(logit_fn_post)/((1 + exp(-logit_fn_post))^2))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Gene_ID) %>%
    dplyr::mutate(effect_size = logit_fn_post - logit_fn_post[Condition == 1]) %>%
    dplyr::mutate(effect_std_error = ifelse(Condition == 1, sd_post, sqrt(sd_post[Condition == 1]^2 + sd_post^2))) %>%
    dplyr::mutate(L2FC_kdeg = ifelse(Condition == 1, 0, log2(kdeg/kdeg[Condition == 1]))) %>%
    dplyr::mutate(pval = 2*pt(-abs(effect_size/effect_std_error), df = 2*nreps - 2 + a_hyper)) %>%
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

  fit <- ashr::ash(effects, ses, method = "fdr")
  lfsr <- fit$result$lfsr

  Effect_sizes_df <- data.frame(Genes_effects, Condition_effects, L2FC_kdegs, effects, ses, lfsr, pval, padj)


  #hyperpars <- c(sdp, theta_o, var_pop, var_of_var, a_hyper, b_hyper)
  #names(hyperpars) <- c("Mean Prior sd", "Mean prior mean", "Variance prior mean", "Variance prior variance", "Variance hyperparam a", "Variance hyperparam b")

  #fast_list <- list(estimate_df, avg_df_fn_bayes, Effect_sizes_df, pmuts_list, hyperpars)
  fast_list <- list(df_fn, avg_df_fn_bayes, Effect_sizes_df, pmuts_list, c(a = a_hyper, b = b_hyper), heterosked_lm)

  names(fast_list) <- c("Fn_Estimates", "Regularized_ests", "Effects_df", "Mut_rates", "Hyper_Parameters", "Mean_Variance_lm")

  return(fast_list)

}
