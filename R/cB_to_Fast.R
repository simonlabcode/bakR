
#' Extract data for efficient analysis from cB
#'
#' This function processes the cB into a form necessary for the parameter estimation function
#' that does not call Stan and is thus much more efficient and scalable.
#' @param cB cB file generated from TL-seq pipeline
#' @param samp_list vector of names of all samples
#' @param c_list vector of names of control (no s4U fed) samples
#' @param type_list vector with 1 entry per sample; 0 = no s4U, 1 = s4U fed
#' @param mut_list vector with 1 entry per sample; 1 = WT, > 1 = different experimental conditions (e.g., KO of protein X)
#' @param rep_list vector with 1 entry per sample that indexes replicate; 1 = 1st replicate, 2 = 2nd replicate, etc.
#' @param tl single numerical value; s4U label time used in s4U fed samples
#' @param nreps single numerical value; number of replicates (assumes same number of replicates for each mut_list index)
#' @param keep_input two element vector; 1st element is highest mut rate accepted in control samples, 2nd element is read count cutoff
#' @importFrom magrittr %>%
#' @return returns dataframe that can be passed to fast analysis
#' @export
cBtofast <- function(cB_raw,
                       c_list,
                       type_list,
                       mut_list,
                       rep_list,
                       tl,
                       nreps,
                       keep_input = c(0.2, 50)){

  cB <- cB_raw %>%
    dplyr::select(sample, XF, GF, TC, n, io, nT)

  samp_list <- unique(cB$sample)

  names(type_list) <- samp_list
  names(mut_list) <- samp_list
  names(rep_list) <- samp_list


  # Helper function:
  getType <- function(s) type_list[paste(s)]
  getMut <- function(s) mut_list[paste(s)]
  getRep <- function(s) rep_list[paste(s)]

  # Get reliable features:
  keep <- DynamicSeq::reliableFeatures(cB, c_list, high_p=keep_input[1], totcut = keep_input[2])

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
#' @param boot_iter Number of times to resample for bootstrapping; default is 50
#' @param pnew Labeled read mutation rate; default of 0 means that model estimates rate from s4U fed data
#' @param pold Unlabeled read mutation rate; default of 0 means that model estimates rate from no-s4U fed data
#' @return list with dataframe of replicate specific estimates as well as dataframe of pooled estimates
#' @importFrom magrittr %>%
#' @export
fast_analysis <- function(df, boot_iter = 50, pnew = 0, pold = 0){

  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))

  #Old mutation rate estimation

  #Trim df and name columns
  if(pnew == 0){
    Mut_data <- df[,1:9]
    colnames(Mut_data) <- c("sample", "XF", "TC", "nT", "n", "fnum", "type", "mut", "reps")

    #New Mutation rate Estimation
    New_data <- Mut_data[Mut_data$type == 1, ]

    New_data$avg_mut <- New_data$TC/New_data$nT
    New_data$avg_mut[is.na(New_data$avg_mut)] <- 0
    New_data$weight_mut <- New_data$avg_mut*New_data$n

    New_data_summary <- New_data %>%
      dplyr::group_by(reps, mut, fnum) %>%
      dplyr::do(purrr::invoke_map_dfc(list(purrr::map_df),
                                      list(list(dplyr::select(., weight_mut), sum),
                                           list(dplyr::select(., n), sum))
      )
      )


    New_data_summary$avg_mut <- New_data_summary$weight_mut/New_data_summary$n

    New_data_ordered <- New_data_summary[order(New_data_summary$avg_mut, decreasing=TRUE), ]

    New_data_cutoff <- New_data_ordered[New_data_ordered$n > 100,]

    pnew <- mean(New_data_cutoff$avg_mut[1:30])

  }

  if(pold == 0){
    #Old mutation rate estimation
    Mut_data <- df[,1:9]
    colnames(Mut_data) <- c("sample", "XF", "TC", "nT", "n", "fnum", "type", "mut", "reps")

    Old_data <- Mut_data[Mut_data$type == 0, ]

    Old_data$avg_mut <- Old_data$TC/Old_data$nT
    Old_data$avg_mut[is.na(Old_data$avg_mut)] <- 0
    Old_data$weight_mut <- Old_data$avg_mut*Old_data$n

    #Old_data$n <- rep(1, times=nrow(Old_data))

    Old_data_summary <- Old_data %>%
      dplyr::group_by(reps, mut, fnum) %>%
      dplyr::do(purrr::invoke_map_dfc(list(purrr::map_df),
                                      list(list(dplyr::select(., weight_mut), sum),
                                           list(dplyr::select(., n), sum))
      )
      )


    Old_data_summary$avg_mut <- Old_data_summary$weight_mut/Old_data_summary$n

    Old_data_ordered <- Old_data_summary[order(Old_data_summary$n, decreasing=TRUE), ]

    Old_data_cutoff <- Old_data_ordered[Old_data_ordered$n > 100,]

    pold <- mean(Old_data_ordered$avg_mut[1:30])

  }



  pmuts <- c(pold, pnew)

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

  resamps <- boot_iter
  fn_boot <- rep(0, times=resamps)
  std_dev <- fn_rep_est


  for(i in 1:ngene){
    for(j in 1:num_conds){
      for(k in 1:nreps){
        #Extract Relevant Data
        TCs <- Mut_data$TC[(Mut_data$reps == k) & (Mut_data$mut == j) & (Mut_data$fnum == i)]
        Us <- Mut_data$nT[(Mut_data$reps == k) & (Mut_data$mut == j) & (Mut_data$fnum == i)]

        #Calculate newness likelihood for each gene
        pnews <- stats::dbinom(TCs, size=Us, prob=pnew)
        polds <- stats::dbinom(TCs, size=Us, prob=pold)

        new_reads <- sum(pnews/(pnews + polds))

        L <- length(TCs)

        fn_rep_est[i, j, k] <- new_reads/L




        #Bootstrapping
        for(l in 1:resamps){
          #Simulate TimeLapse data with fn_est and estimate mut rates
          TCs <- rep(0, times=L)
          Us <- TCs
          Us <- stats::rbinom(n=L, size=200, prob=0.25)
          new_ID <- purrr::rbernoulli(n=L, p=fn_rep_est[i,j,k]) + 1

          TCs <- stats::rbinom(n=L, size=Us, prob=pmuts[new_ID])

          pnews <- stats::dbinom(TCs, size=Us, prob=pnew)
          polds <- stats::dbinom(TCs, size=Us, prob=pold)

          new_reads <- sum(pnews/(pnews + polds))

          fn_boot[l] <- new_reads/L

          if(fn_boot[l] == 1){
            fn_boot[l] <- 1 - 1/(length(TCs) + 1)
          }else if(fn_boot[l] == 0){
            fn_boot[l] <- 1/(length(TCs) + 1)
          }
        }
        #Standard deviations of bootstrapped estimates
        std_dev[i, j, k] <- stats::sd(logit(fn_boot))
        R_ID[i, j, k] <- k
        FN_ID[i, j, k] <- i
        MT_ID[i, j, k] <- j
      }
    }
  }


  logit_fn_rep <- logit(fn_rep_est)

  R <- rep(1:nreps, times=)

  logit_fn <- as.vector(logit_fn_rep)
  logit_fn_sd <- as.vector(std_dev)
  fn_estimate <- as.vector(fn_rep_est)
  Replicate <- as.vector(R_ID)
  Condition <- as.vector(MT_ID)
  Gene_ID <- as.vector(FN_ID)


  estimate_df <- data.frame(logit_fn, logit_fn_sd, fn_estimate, Replicate, Condition, Gene_ID)

  df_fn <- estimate_df[order(estimate_df$Gene_ID, estimate_df$Condition, estimate_df$Replicate),]

  nreps <- max(df_fn$Replicate)

  #Average over replicates
  avg_df_fn <- df_fn %>% dplyr::group_by(Gene_ID, Condition) %>%
    summarize(avg_logit_fn = mean(logit_fn),
              sd_boot = mean(logit_fn_sd),
              sd_logit_fn = sd(logit_fn))

  #Calcualte population averages
  sdp <- mean(avg_df_fn$sd_logit_fn)
  theta_o <- mean(avg_df_fn$avg_logit_fn)

  #Adjust average fn according to Bayesian normal model with known sd
  avg_df_fn_bayes <- avg_df_fn %>%
    mutate(avg_logit_fn = (avg_logit_fn*(nreps*(1/sd_boot)))/(nreps/sd_boot + sdp) + (theta_o*sdp)/(nreps/sd_boot + sdp) )


  fn_list <- list(estimate_df, avg_df_fn_bayes, pmuts)

  return(fn_list)

}
