#' Identify features (e.g., transcripts) with high quality data
#'
#' This function identifies all features (e.g., transcripts, exons, etc.) for which the mutation rate
#' is below a set threshold in the control (no s4U) sample and which have more reads than a set threshold
#' in all samples. If there is no -s4U sample, then only the read count cutoff is considered.
#'
#' @param obj Object of class DynamicSeqData
#' @param high_p highest mutation rate accepted in control samples
#' @param totcut readcount cutoff
#' @importFrom magrittr %>%
#' @return vector of gene names that passed reliability filter
#' @export
reliableFeatures <- function(obj,
                             high_p = 0.2,
                             totcut = 50){

  cB <- obj$cB
  nsamps <- length(unique(cB$sample))



  if(sum(obj$metadf$tl == 0) > 0){
    y <- obj$cB %>%
      dplyr::ungroup() %>%
      dplyr::filter(sample %in% unique(sample),
                    !grepl('__', XF)) %>%
      dplyr::mutate(totTC = TC*n*ifelse(obj$metadf[sample, "tl"]==0, 1, 0) ) %>%
      dplyr::group_by(sample, XF) %>%
      dplyr::summarize(tot_mut = sum(totTC),
                       totcounts = sum(n)) %>%
      dplyr::filter(totcounts >= totcut) %>%
      dplyr::filter(tot_mut/totcounts < high_p) %>%
      dplyr::ungroup( ) %>%
      dplyr::group_by(XF) %>%
      dplyr::summarize(counts = dplyr::n()) %>%
      dplyr::filter(counts == nsamps) %>%
      dplyr::select(XF) %>%
      unlist() %>%
      unique()
  }else{
    y <- obj$cB %>%
      dplyr::ungroup() %>%
      dplyr::filter(sample %in% unique(sample),
                    !grepl('__', XF)) %>%
      dplyr::group_by(sample, XF) %>%
      dplyr::summarize(totcounts = sum(n)) %>%
      dplyr::filter(totcounts >= totcut) %>%
      dplyr::ungroup( ) %>%
      dplyr::group_by(XF) %>%
      dplyr::summarize(counts = dplyr::n()) %>%
      dplyr::filter(counts == nsamps) %>%
      dplyr::select(XF) %>%
      unlist() %>%
      unique()
  }


  return(y)
}

#' Curate data in DynamicSeqData object for statistical modeling
#'
#' This function obtains the data structures necessary to analyze nucleotide recoding sequencing data with any of the
#' statistical models implemented by \code{DynamicSeqFit}.
#'
#' @param obj An object of class DynamicSeqData
#' @param keep_input two element vector; 1st element is highest mut rate accepted in control samples, 2nd element is read count cutoff
#' @param Stan Boolean; if TRUE, then data_list that can be passed to Stan is curated
#' @param Fast Boolean; if TRUE, then dataframe that can be passed to fast_analysis() is curated
#' @param FOI Features of interest; character vector containing names of features to analyze
#' @param concat Boolean; If TRUE, FOI is concatenated with output of reliableFeatures
#' @return returns list of objects that can be passed to \code{TL_stan} and/or \code{fast_analysis}. Those objects are:
#' \itemize{
#'  \item Stan_data; list that can be passed to \code{TL_stan} with Hybrid_Fit = FALSE. Consistents of metadata as well as data that
#'  Stan will analyze. Data to be analyzed consists of equal length vectors. The contents of Stan_data are:
#'  \itemize{
#'   \item NE; Number of datapoints for Stan to analyze (NE = Number of Elements)
#'   \item NF; Number of features in dataset
#'   \item TP; Numerical indicator of s4U feed (0 = no s4U feed, 1 = s4U fed)
#'   \item FE; Numerical indicator of feature
#'   \item num_mut; Number of U-to-C mutations observed in a particular set of reads
#'   \item MT; Numerical indicator of experimental condition (Exp_ID from metadf)
#'   \item nMT; Number of experimental conditions
#'   \item R; Numerical indicator of replicate
#'   \item nrep; Number of replicates (analysis requires same number of replicates of all conditions)
#'   \item num_obs; Number of reads with identical data (number of mutations, feature of origin, and sample of origin)
#'   \item tl; Vector of label times for each experimental condition
#'   \item U_cont; Log2-fold-difference in U-content for a feature in a sample relative to average U-content for that sample
#'   \item Avg_Reads; Standardized log10(average read counts) for a particular feature in a particular condition, averaged over
#'   replicates
#'   \item sdf; Dataframe that maps numerical feature ID to original feature name. Also has read depth information
#'   \item sample_lookup; Lookup table relating MT and R to the original sample name
#'  }
#'  \item Fast_df; A data frame that can be passed to \code{fast_analysis}. The contents of Fast_df are:
#'  \itemize{
#'   \item sample; Original sample name
#'   \item XF; Original feature name
#'   \item TC; Number of T to C mutations
#'   \item nT; Number of Ts in read
#'   \item n; Number of identical observations
#'   \item fnum; Numerical indicator of feature
#'   \item type; Numerical indicator of s4U feed (0 = no s4U feed, 1 = s4U fed)
#'   \item mut; Numerical indicator of experimental condition (Exp_ID from metadf)
#'   \item reps; Numerical indicator of replicate
#'  }
#'  \item Count_Matrix; A matrix with read count information. Each column represents a sample and each row represents a feature.
#'  Each entry is the raw number of read counts mapping to a particular feature in a particular sample. Column names are the corresponding
#'  sample names and row names are the corresponding feature names.
#' }
#' @importFrom magrittr %>%
#' @export
cBprocess <- function(obj,
                       keep_input=c(0.2, 50),
                       Stan = TRUE,
                       Fast = TRUE,
                       FOI = c(),
                       concat = TRUE){

  ## Check obj
  if(class(obj) != "DynamicSeqData"){
    stop("obj must be of class DynamicSeqData")
  }

  ## Check keep_input
  if(!all(is.numeric(keep_input))){
    stop("All elements of keep_input must be numeric")
  }else if(length(keep_input) != 2){
    stop("keep_input must be a vector of length 2. The 1st element should be a number between 0 and 1 representing the maximum acceptable mutation rate
         in the no s4U control sample. The 2nd element should be anumber > 0 representing the read count cutoff for filtered features")
  }else if(keep_input[1] < 0){
    stop("1st element of keep_input must be >= 0")
  }else if(keep_input[1] > 1){
    stop("1st element of keep_input must be <= 1")
  }else if(keep_input[2] < 0){
    stop("2nd element of keep_input must be >= 0")
  }

  ## Check Stan
  if(!is.logical(Stan)){
    stop("Stan must be logical (TRUE or FALSE)")
  }

  ## Check Fast_prep
  if(!is.logical(Fast)){
    stop("Fast must be logical (TRUE or FALSE")
  }

  ## Check FOI
  if(!is.null(FOI)){
    if(typeof(obj$cB$XF) != typeof(FOI)){
      warning("FOI should be the same data type as cB$XF in the DynamicSeqData object; if it is not none of the feature of interest will be found
            in the cB.")
    }
  }

  ## Check concat
  if(!is.logical(concat)){
    stop("concat must be logical (TRUE or FALSE)")
  }


  cB <- obj$cB
  metadf <- obj$metadf

  samp_list <- unique(cB$sample)

  c_list <- rownames(metadf[metadf$tl == 0,])

  s4U_list <- samp_list[!(samp_list %in% c_list)]

  type_list <- ifelse(metadf[samp_list, "tl"] == 0, 0, 1)
  mut_list <- metadf[samp_list, "Exp_ID"]
  rep_list <- metadf[samp_list,] %>% dplyr::group_by(tl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup() %>% dplyr::select(r_id)
  rep_list <- rep_list$r_id

  metadf <- metadf[samp_list, ] %>% dplyr::group_by(tl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup()

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
    message("Finding reliable Features")

    reliables <- DynamicSeq::reliableFeatures(obj, high_p = keep_input[1], totcut = keep_input[2])
    keep <- c(FOI, reliables[!(reliables %in% FOI)])
  }else{
    keep <- FOI
  }



  message("Filtering out unwanted or unreliable features")

  # This df is created in both cases
  ranked_features_df  <- cB %>%
    dplyr::ungroup() %>%
    dplyr::filter(XF %in% keep) %>%
    dplyr::group_by(XF) %>%
    dplyr::summarize(n = sum(n)) %>%
    dplyr::mutate(fnum = order(XF)) %>%
    dplyr::arrange(fnum) %>%
    dplyr::select(XF, fnum)


  message("Processing data...")

  # Create count dataframe
  Counts_df <- cB %>%
    dplyr::ungroup() %>%
    dplyr::filter(XF %in% keep) %>%
    dplyr::group_by(XF, sample) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::right_join(ranked_features_df, by = 'XF') %>% dplyr::ungroup()


  # Make count matrix
  Cnt_mat <- matrix(0, ncol = length(samp_list), nrow = length(unique(Counts_df$XF)))

  for(s in seq_along(samp_list)){
    Cnt_mat[,s] <- Counts_df$n[Counts_df$sample == samp_list[s]]
  }

  rownames(Cnt_mat) <- Counts_df$XF[Counts_df$sample == samp_list[1]]
  colnames(Cnt_mat) <- samp_list

  rm(Counts_df)

  ## U content estimation
  sdf_U <- cB %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample, XF, TC, nT) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::right_join(ranked_features_df, by = 'XF') %>%
    dplyr::ungroup()

  slist = samp_list
  tlist = type_list
  mlist = mut_list
  rlist = rep_list



  kp = keep

  df_U <- sdf_U %>%
    dplyr::ungroup() %>%
    dplyr::filter(sample %in% slist)

  df_U$type <- paste(df_U$sample) %>% purrr::map_dbl(function(x) getType(x))
  df_U$type <- as.integer(df_U$type)

  df_U$mut <- paste(df_U$sample) %>% purrr::map_dbl(function(x) getMut(x))
  df_U$mut <- as.integer(df_U$mut)

  df_U$reps <- paste(df_U$sample) %>% purrr::map_dbl(function(x) getRep(x))
  df_U$reps <- as.integer(df_U$reps)

  sample_lookup <- df_U[df_U$type == 1, c("sample", "mut", "reps")] %>% dplyr::distinct()


  df_global_U <- df_U[df_U$type == 1, ] %>% dplyr::group_by(reps, mut) %>%
    dplyr::summarise(tot_avg_Us = sum(nT*n)/sum(n)) %>% dplyr::ungroup()

  df_feature_U <- df_U[df_U$type == 1, ] %>% dplyr::group_by(reps, mut, fnum) %>%
    dplyr::summarise(feature_avg_Us = sum(nT*n)/sum(n)) %>% dplyr::ungroup()

  df_U_tot <- merge(df_global_U, df_feature_U, by = c("mut", "reps"))

  df_U_tot <- df_U_tot %>% dplyr::mutate(U_factor = log(feature_avg_Us/tot_avg_Us)) %>%
    dplyr::select(mut, reps, fnum, U_factor)

  if(Stan){

    sdf <- cB %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample, XF, TC) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::right_join(ranked_features_df, by = 'XF') %>% dplyr::ungroup()

    ##### Run Stan model ---------

    d = sdf


    df <- d %>%
      dplyr::ungroup() %>%
      dplyr::filter(sample %in% slist)

    df$type <- paste(df$sample) %>% purrr::map_dbl(function(x) getType(x))
    df$type <- as.integer(df$type)

    df$mut <- paste(df$sample) %>% purrr::map_dbl(function(x) getMut(x))
    df$mut <- as.integer(df$mut)

    df$reps <- paste(df$sample) %>% purrr::map_dbl(function(x) getRep(x))
    df$reps <- as.integer(df$reps)


    df <- df  %>%
      dplyr::group_by(XF, fnum, type, mut, TC, reps) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fnum) %>%
      dplyr::arrange(type, .by_group = TRUE)



    df <- merge(df, df_U_tot, by = c("fnum", "mut", "reps"))
    df <- df[order(df$fnum, df$mut, df$reps), ]

    FE = df$fnum
    NE <- dim(df)[1]
    NF <- length(kp)
    TP <- df$type
    MT <- df$mut
    nMT <- length(unique(MT))
    R <- df$reps

    num_mut <- df$TC # Number of mutations
    num_obs <- df$n # Number of times observed

    ## Calculate Avg. Read Counts
    Avg_Counts <- df %>% dplyr::ungroup() %>% dplyr::group_by(fnum, mut) %>%
      dplyr::summarise(Avg_Reads = sum(n)/nreps) %>% dplyr::ungroup()

    Avg_Reads <- matrix(0, ncol = nMT, nrow = NF)

    tls <-rep(0, times = nMT)

    for(f in 1:NF){
      for(i in 1:nMT){
        Avg_Reads[f,i] <- (mean(log10(Avg_Counts$Avg_Reads[(Avg_Counts$mut == i) & (Avg_Counts$fnum == f)])) - mean(log10(Avg_Counts$Avg_Reads[Avg_Counts$mut == i])))/sd(log10(Avg_Counts$Avg_Reads[Avg_Counts$mut == i]))

      }
    }

    for(m in 1:nMT){
      tls[m] <- unique(metadf$tl[(metadf$Exp_ID == m) & (metadf$tl != 0)])
    }



    # Will have to change to make Stan models compatible
    # with a vector or matrix of tls
    data_list <- list(
      NE = NE, #Number of reads
      NF = NF, # Number of gene
      TP = TP, # Flag the control sample EV_0
      FE = FE, # gene
      num_mut = num_mut,
      MT = MT,
      nMT = nMT,
      R = R,
      nrep = nreps,
      num_obs = num_obs,
      tl = tls,
      U_cont = df$U_factor,
      Avg_Reads = Avg_Reads,
      sdf = sdf,
      sample_lookup = sample_lookup
    )

  }


  if(Stan & Fast){
    out <- list(data_list, df_U, Cnt_mat)
    names(out) <- c("Stan_data", "Fast_df", "Count_Matrix")
  }else if(!Stan){
    out <- list(df_U, Cnt_mat)
    names(out) <- c("Fast_df", "Count_Matrix")
  }else{
    out <- list(data_list, Cnt_mat)
    names(out) <- c("Stan_data", "Count_Matrix")
  }

  return(out)
}
