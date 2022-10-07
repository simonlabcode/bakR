#' Identify features (e.g., transcripts) with high quality data
#'
#' This function identifies all features (e.g., transcripts, exons, etc.) for which the mutation rate
#' is below a set threshold in the control (no s4U) sample and which have more reads than a set threshold
#' in all samples. If there is no -s4U sample, then only the read count cutoff is considered. Additional
#' filtering options are only relevant if working with short RNA-seq read data. This includes filtering out
#' features with extremely low empirical U-content (i.e., the average number of Us in sequencing reads from
#' that feature) and those with very few reads having at least 3 Us in them.
#'
#' @param obj Object of class bakRData
#' @param high_p highest mutation rate accepted in control samples
#' @param totcut readcount cutoff
#' @param Ucut Must have a fraction of reads with 2 or less Us less than this cutoff in all samples
#' @param AvgU Must have an average number of Us greater than this
#' @importFrom magrittr %>%
#' @return vector of gene names that passed reliability filter
#' @examples
#' \donttest{
#'
#' # Load cB
#' data("cB_small")
#'
#' # Load metadf
#' data("metadf")
#'
#' # Create bakRData
#' bakRData <- bakRData(cB_small, metadf)
#'
#' # Find reliable features
#' features_to_keep <- reliableFeatures(obj = bakRData)
#' }
#' @export
reliableFeatures <- function(obj,
                             high_p = 0.2,
                             totcut = 50,
                             Ucut = 0.25,
                             AvgU = 4){

  # Bind variables locally to resolve devtools::check() Notes
  XF <- TC <- n <- totTC <- nT <- n2U <- nmore <- totcounts <- NULL
  tot_mut <- f2U <- avgU <- counts <- NULL

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
                       totcounts = sum(n),
                       avgU = sum(nT*n)/sum(n),
                       n2U = sum(n[nT <= 2]),
                       nmore = sum(n[nT > 2])) %>%
      dplyr::mutate(f2U = n2U/(nmore + n2U)) %>%
      dplyr::filter(totcounts >= totcut) %>%
      dplyr::filter(tot_mut/totcounts < high_p) %>%
      dplyr::filter(f2U < Ucut) %>%
      dplyr::filter(avgU > AvgU) %>%
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
      dplyr::summarize(totcounts = sum(n),
                       avgU = sum(nT*n)/sum(n),
                       n2U = sum(n[nT <= 2]),
                       nmore = sum(n[nT > 2])) %>%
      dplyr::mutate(f2U = n2U/(nmore + n2U)) %>%
      dplyr::filter(totcounts >= totcut,
                    f2U < Ucut) %>%
      dplyr::filter(avgU > AvgU) %>%
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

#' Curate data in bakRData object for statistical modeling
#'
#' \code{cBprocess} creates the data structures necessary to analyze nucleotide recoding RNA-seq data with any of the
#' statistical model implementations in \code{bakRFit}. The input to \code{cBprocess} must be an object of class
#' `bakRData`. The output can contain data passable to \code{fast_analysis} (if Fast == TRUE), \code{TL_stan} with StanFit = TRUE if
#' Stan == TRUE, or both.
#'
#' The 1st step executed by \code{cBprocess} is to find the names of features which are deemed "reliable". A reliable feature is one with
#' sufficient read coverage in every single sample (i.e., > totcut reads in all samples) and limited mutation content in all -s4U
#' control samples (i.e., < high_p mutation rate in all samples lacking s4U feeds). This is done with a call to \code{reliableFeatures}.
#' The 2nd step is to extract only reliableFeatures from the cB dataframe in the `bakRData` object. During this process, a numerical
#' ID is given to each reliableFeature, with the numerical ID corresponding to the order in which each feature is found in the original cB
#' (this might typically be alphabetical order).
#'
#' The 3rd step is to prepare a dataframe where each row corresponds to a set of n identical reads (that is they come from the same sample
#' and have the same number of mutations and Us). Part of this process involves assigning an arbitrary numerical ID to each replicate in each
#' experimental condition. The numerical ID will correspond to the order the sample appears in metadf. The outcome of this step is multiple
#' dataframes with variable information content. These include a dataframe with information about read counts in each sample, one which logs
#' the U-contents of each feature, one which is compatible with \code{fast_analysis} and thus groups reads by their number of mutations as
#' well as their number of Us, and one which is compatible with \code{TL_stan} with StanFit == TRUE and thus groups ready by only their number
#' of mutations. At the end of this step, two other smaller data structures are created, one which is an average count matrix (a count matrix
#' where the ith row and jth column corresponds to the average number of reads mappin to feature i in experimental condition j, averaged over
#' all replicates) and the other which is a sample lookup table that relates the numerical experimental and replicate IDs to the original
#' sample name.
#'
#' If FOI is non-null and concat == TRUE, the features listed in FOI will be included in the list of reliable features that make it past
#' filtering. If FOI is non-null and concat == FALSE, the features listed in FOI will be the only reliable features that make it past filtering.
#' If FOI is null and concat == FALSE, an error will be thrown.
#'
#'
#' @param obj An object of class bakRData
#' @param high_p Numeric; Any transcripts with a mutation rate (number of mutations / number of Ts in reads) higher than this in any no s4U control
#' samples are filtered out
#' @param totcut Numeric; Any transcripts with less than this number of sequencing reads in any sample are filtered out
#' @param Ucut Numeric; All transcripts must have a fraction of reads with 2 or less Us less than this cutoff in all samples
#' @param AvgU Numeric; All transcripts must have an average number of Us greater than this cutoff in all samples
#' @param Stan Boolean; if TRUE, then data_list that can be passed to 'Stan' is curated
#' @param Fast Boolean; if TRUE, then dataframe that can be passed to fast_analysis() is curated
#' @param FOI Features of interest; character vector containing names of features to analyze
#' @param concat Boolean; If TRUE, FOI is concatenated with output of reliableFeatures
#' @return returns list of objects that can be passed to \code{TL_stan} and/or \code{fast_analysis}. Those objects are:
#' \itemize{
#'  \item Stan_data; list that can be passed to \code{TL_stan} with Hybrid_Fit = FALSE. Consists of metadata as well as data that
#'  'Stan' will analyze. Data to be analyzed consists of equal length vectors. The contents of Stan_data are:
#'  \itemize{
#'   \item NE; Number of datapoints for 'Stan' to analyze (NE = Number of Elements)
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
#'   \item Avg_Reads_natural; Unstandardized average read counts for a particular feature in a particular condition, averaged over
#'   replicates. Used for \code{plotMA}
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
#' @examples
#' \donttest{
#'
#' # Load cB
#' data("cB_small")
#'
#' # Load metadf
#' data("metadf")
#'
#' # Create bakRData
#' bakRData <- bakRData(cB_small, metadf)
#'
#' # Preprocess data
#' data_for_bakR <- cBprocess(obj = bakRData)
#' }
#' @export
cBprocess <- function(obj,
                       high_p = 0.2,
                       totcut = 50,
                       Ucut = 0.25,
                       AvgU = 4,
                       Stan = TRUE,
                       Fast = TRUE,
                       FOI = c(),
                       concat = TRUE){


  # Bind variables locally to resolve devtools::check() Notes
  tl <- ctl <- Exp_ID <- r_id <- XF <- n <- fnum <- TC <- nT <- reps <- NULL
  mut <- feature_avg_Us <- tot_avg_Us <- U_factor <- type <- NULL


  ## Check obj
  if(!inherits(obj, "bakRData")){
    stop("obj must be of class bakRData")
  }

  ## Check high_p
  if(!is.numeric(high_p)){
    stop("high_p must be numeric")
  }else if( (high_p < 0) | (high_p > 1) ){
    stop("high_p must be between 0 and 1")
  }else if (high_p < 0.01){
    warning("high_p is abnormally low (< 0.01); many features will by pure chance have a higher mutation rate than this in a -s4U control and thus get filtered out")
  }

  ## Check totcut
  if(!is.numeric(totcut)){
    stop("totcut must be numeric")
  }else if( totcut < 0 ){
    stop("totcut must be greater than 0")
  }else if(totcut > 5000){
    warning("totcut is abnormally high (> 5000); many features will not have this much coverage in every sample and thus get filtered out.")
  }

  ## Check Ucut
  if(!is.numeric(Ucut)){
    stop("Ucut must be numeric")
  }else if( Ucut < 0 ){
    stop("Ucut must be greater than 0")
  }else if(Ucut > 0.5 ){
    warning("Ucut is abnormally high; you are allowing > 50% of reads to have 2 or less Us.")
  }

  ## Check AvgU
  if(!is.numeric(AvgU)){
    stop("AvgU must be numeric")
  }else if(AvgU < 0){
    stop("AvgU must be greater than or equal to 0")
  }else if (AvgU > 50){
    warning("AvgU is abnormally high; you are requiring an average number of Us greater than 50")
  }else if(AvgU < 4){
    warning("AvgU is abnormally low; you are allowing an average of less than 4 Us per read, which may model convergence issues.")
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
      warning("FOI should be the same data type as cB$XF in the bakRData object; if it is not none of the feature of interest will be found
            in the cB.")
    }
  }else{
    if(concat == FALSE){
      stop("concat cannot be FALSE if FOI is null; this would cause no features to make it past filtering")
    }
  }

  ## Check concat
  if(!is.logical(concat)){
    stop("concat must be logical (TRUE or FALSE)")
  }

  ## Create vectors that map sample to important characteristics:
    # samp_list = vector of sample names
    # c_list = vector of -s4U control sample names
    # s4U_list = vector of samples that are treated with s4U
    # type_list = 1 if sample i in samp_list is s4U treated; 0 if not
    # mut_list = numerical ID for each experimental condition as in metadf
      # 1 = reference
      # > 1 = experimental conditions
    # rep_list = vector of replicate IDs


  cB <- obj$cB
  metadf <- obj$metadf

  samp_list <- unique(cB$sample)

  c_list <- rownames(metadf[metadf$tl == 0,])

  s4U_list <- samp_list[!(samp_list %in% c_list)]

  type_list <- ifelse(metadf[samp_list, "tl"] == 0, 0, 1)
  mut_list <- metadf[samp_list, "Exp_ID"]

  rep_list <- metadf[samp_list,] %>% dplyr::mutate(ctl = ifelse(tl == 0, 0, 1)) %>%
    dplyr::group_by(ctl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup() %>% dplyr::select(r_id)
  rep_list <- rep_list$r_id

  # Add replicate ID and s4U treatment status to metadf
  metadf <- metadf[samp_list, ] %>% dplyr::mutate(ctl = ifelse(tl == 0, 0, 1)) %>%
    dplyr::group_by(ctl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup()

  names(type_list) <- samp_list
  names(mut_list) <- samp_list
  names(rep_list) <- samp_list

  # Make vector of number of replicates of each condition
  nreps <- rep(0, times = max(mut_list))
  for(i in 1:max(mut_list)){
    nreps[i] <- max(rep_list[mut_list == i])
  }


  # Helper function:
  getType <- function(s) type_list[paste(s)]
  getMut <- function(s) mut_list[paste(s)]
  getRep <- function(s) rep_list[paste(s)]

  # Get reliable features:
  if(concat == TRUE | is.null(FOI)){
    message("Finding reliable Features")

    reliables <- bakR::reliableFeatures(obj, high_p = high_p, totcut = totcut, Ucut = Ucut, AvgU = AvgU)
    keep <- c(FOI, reliables[!(reliables %in% FOI)])
  }else{
    keep <- FOI
  }

  if((length(keep) == 0) | (is.null(keep))){
    stop("No features made it past filtering.Try increasing the read count or -s4U background mutation rate cutoffs.")
  }


  message("Filtering out unwanted or unreliable features")

  # Map each reliable feature to a numerical feature ID (fnum)
  ranked_features_df  <- cB %>%
    dplyr::ungroup() %>%
    dplyr::filter(XF %in% keep) %>%
    dplyr::group_by(XF) %>%
    dplyr::summarize(n = sum(n)) %>%
    dplyr::mutate(fnum = order(XF)) %>%
    dplyr::arrange(fnum) %>%
    dplyr::select(XF, fnum)


  message("Processing data...")

  # Make data frame with read count information
  Counts_df <- cB %>%
    dplyr::ungroup() %>%
    dplyr::filter(XF %in% keep) %>%
    dplyr::group_by(XF, sample) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::right_join(ranked_features_df, by = 'XF') %>% dplyr::ungroup()


  # Make count matrix that is DESeq2 compatible
  Cnt_mat <- matrix(0, ncol = length(samp_list), nrow = length(unique(Counts_df$XF)))

  for(s in seq_along(samp_list)){
    Cnt_mat[,s] <- Counts_df$n[Counts_df$sample == samp_list[s]]
  }

  rownames(Cnt_mat) <- Counts_df$XF[Counts_df$sample == samp_list[1]]
  colnames(Cnt_mat) <- samp_list

  rm(Counts_df)

  # Empirical U-content calculations
    # = average number of Us in sequencing reads originating from each feature
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

  ## Add sample characteristic details to U-content data frame

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


  ## Calculate global average U-content so that feature-specific difference
  ## from average can be calculated

  df_global_U <- df_U[df_U$type == 1, ] %>% dplyr::group_by(reps, mut) %>%
    dplyr::summarise(tot_avg_Us = sum(nT*n)/sum(n)) %>% dplyr::ungroup()

  df_feature_U <- df_U[df_U$type == 1, ] %>% dplyr::group_by(reps, mut, fnum) %>%
    dplyr::summarise(feature_avg_Us = sum(nT*n)/sum(n)) %>% dplyr::ungroup()

  df_U_tot <- dplyr::left_join(df_global_U, df_feature_U, by = c("mut", "reps"))

  # U_factor is log-fold difference in feature specific U-content from global average
    # Used in Stan model to properly adjust population average Poisson mutation rates
  df_U_tot <- df_U_tot %>% dplyr::mutate(U_factor = log(feature_avg_Us/tot_avg_Us)) %>%
    dplyr::select(mut, reps, fnum, U_factor)

  if(Stan){

    # Filter out unreliable features and assign feature ID to each feature
    sdf <- cB %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample, XF, TC) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::right_join(ranked_features_df, by = 'XF') %>% dplyr::ungroup()


    ### Curate data for 'Stan' models

    d = sdf


    ## Add sample characteristic info to data frame

    df <- d %>%
      dplyr::ungroup() %>%
      dplyr::filter(sample %in% slist)

    df$type <- paste(df$sample) %>% purrr::map_dbl(function(x) getType(x))
    df$type <- as.integer(df$type)

    df$mut <- paste(df$sample) %>% purrr::map_dbl(function(x) getMut(x))
    df$mut <- as.integer(df$mut)

    df$reps <- paste(df$sample) %>% purrr::map_dbl(function(x) getRep(x))
    df$reps <- as.integer(df$reps)

    ## Remove any unnecessary columns
    df <- df  %>%
      dplyr::group_by(XF, fnum, type, mut, TC, reps) %>%
      dplyr::summarise(n = sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(fnum) %>%
      dplyr::arrange(type, .by_group = TRUE)



    # Add U-content information
    df <- dplyr::left_join(df, df_U_tot, by = c("fnum", "mut", "reps"))

    df <- df[order(df$fnum, df$mut, df$reps), ]

    # Feature ID
    FE = df$fnum

    # Number of data points
    NE <- dim(df)[1]

    # Number of features analyzed
    NF <- length(kp)

    # s4U ID (1 = s4U treated, 0 = not)
    TP <- df$type

    # Experimental condition ID
    MT <- df$mut

    # Number of experimental conditions
    nMT <- length(unique(MT))

    # Replicate ID
    R <- df$reps

    # Number of mutations
    num_mut <- df$TC

    # Number of identical observations
    num_obs <- df$n

    ## Calculate Avg. Read Counts
    Avg_Counts <- df %>% dplyr::ungroup() %>% dplyr::group_by(fnum, mut) %>%
      dplyr::summarise(Avg_Reads = sum(n)/nreps) %>% dplyr::ungroup()

    Avg_Reads <- matrix(0, ncol = nMT, nrow = NF)
    Avg_Reads_natural <- Avg_Reads

    tls <-rep(0, times = nMT)

    # Calculate average read counts on log10 and natural scales
      # log10 scale read counts used in 'Stan' model
      # natural scale read counts used in plotting function (plotMA())
    for(f in 1:NF){
      for(i in 1:nMT){
        Avg_Reads[f,i] <- (mean(log10(Avg_Counts$Avg_Reads[(Avg_Counts$mut == i) & (Avg_Counts$fnum == f)])) - mean(log10(Avg_Counts$Avg_Reads[Avg_Counts$mut == i])))/stats::sd(log10(Avg_Counts$Avg_Reads[Avg_Counts$mut == i]))
        Avg_Reads_natural[f,i] <- mean(Avg_Counts$Avg_Reads[(Avg_Counts$mut == i) & (Avg_Counts$fnum == f)])
      }
    }

    # s4U label time in each experimental condition
    for(m in 1:nMT){
      tls[m] <- unique(metadf$tl[(metadf$Exp_ID == m) & (metadf$tl != 0)])
    }

    # Add tl info to fast_df
    tl_df <- data.frame(tl = tls,
                        mut = 1:length(tls))

    df_U <- dplyr::left_join(df_U, tl_df, by = "mut")


    # data passed to 'Stan' model (MCMC implementation)
    data_list <- list(
      NE = NE,
      NF = NF,
      TP = TP,
      FE = FE,
      num_mut = num_mut,
      MT = MT,
      nMT = nMT,
      R = R,
      nrep_vect = nreps,
      nrep = max(nreps),
      num_obs = num_obs,
      tl = tls,
      U_cont = df$U_factor,
      Avg_Reads = Avg_Reads,
      Avg_Reads_natural = Avg_Reads_natural,
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
