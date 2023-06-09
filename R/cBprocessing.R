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
  tot_mut <- f2U <- avgU <- counts <- type <- totTs <- NULL
  `.` <- list

  cBdt <- as.data.frame(obj$cB)
  metadf <- obj$metadf

  nsamps <- length(unique(cBdt$sample))


  samp_list <- unique(cBdt$sample)

  c_list <- rownames(metadf[metadf$tl == 0,])

  s4U_list <- samp_list[!(samp_list %in% c_list)]

  type_list <- ifelse(metadf[samp_list, "tl"] == 0, 0, 1)

  # Create mut and reps dictionary
  ID_dict <- data.frame(sample = rownames(metadf),
                        type = type_list)


  cBdt <- dplyr::left_join(cBdt, ID_dict, by = "sample") %>%
    dplyr::ungroup()


  cBdt <- data.table::setDT(cBdt)


  if(sum(obj$metadf$tl == 0) > 0){

    cBdt <- cBdt[sample %in% unique(sample) & !grepl("__", XF)]


    cBdt <- cBdt[, `:=`(totTC = TC * n * abs(type - 1))]


    cBdt <- cBdt[, .(tot_mut = sum(totTC), totcounts = sum(n), totTs = sum(n*nT),
                     avgU = sum(nT * n)/sum(n), n2U = sum(n[nT <=2]),
                     nmore = sum(n[nT > 2])), keyby = .(sample,XF)]

    cBdt <- cBdt[, `:=`(f2U = n2U/(nmore + n2U))]

    cBdt <- cBdt[(totcounts >= totcut) & (tot_mut/totTs < high_p) &
                   (f2U < Ucut) & (avgU > AvgU)]

    cBdt <- cBdt[, .(counts = .N), keyby = .(XF)]

    cBdt <- dplyr::as_tibble(cBdt)

    cBdt <- cBdt[cBdt$counts == nsamps,]

    y <- unique(unlist(cBdt$XF))



  }else{

    cBdt <- cBdt[sample %in% unique(sample) & !grepl("__", XF)]


    cBdt <- cBdt[, .(totcounts = sum(n),
                     avgU = sum(nT * n)/sum(n), n2U = sum(n[nT <=2]),
                     nmore = sum(n[nT > 2])), keyby = .(sample,XF)]

    cBdt <- cBdt[, `:=`(f2U = n2U/(nmore + n2U))]

    cBdt <- cBdt[(totcounts >= totcut) &
                   (f2U < Ucut) & (avgU > AvgU)]

    cBdt <- cBdt[, .(counts = .N), keyby = .(XF)]

    cBdt <- dplyr::as_tibble(cBdt)

    cBdt <- cBdt[cBdt$counts == nsamps,]

    y <- unique(unlist(cBdt$XF))


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
  mut <- feature_avg_Us <- tot_avg_Us <- U_factor <- type <- Ucont <- NULL
  `.` <- list


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


  cB <- dplyr::as_tibble(obj$cB)
  metadf <- as.data.frame(obj$metadf)

  samp_list <- unique(cB$sample)

  c_list <- rownames(metadf[metadf$tl == 0,])

  s4U_list <- samp_list[!(samp_list %in% c_list)]

  type_list <- ifelse(metadf[samp_list, "tl"] == 0, 0, 1)
  mut_list <- metadf[samp_list, "Exp_ID"]

  # If number of -s4U controls > +s4U controls, wrap R_ID for -s4Us
    # So if there are 3 -s4U replicates and 2 + s4U, -s4U R_ID should be: 1, 2, 1
    # Let R_raw but R_ID for -s4U calculated as with +s4U samples (so 1, 2, 3 in above example)
    # Then R_ID = ((R_raw - 1) %% nreps_+s4U) + 1, where nreps_+s4U is number of replicates in +s4U sample
  rep_list <- metadf[samp_list,] %>% dplyr::mutate(ctl = ifelse(tl == 0, 0, 1)) %>%
    dplyr::group_by(ctl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup()

  ### Calculate -s4U adjusted r_id
  # 1) Determine number of +s4U replicates in each sample
  nrep_s4U <- rep_list %>% dplyr::filter(ctl == 1) %>%
    dplyr::group_by(Exp_ID) %>%
    dplyr::summarise(nreps = max(r_id)) %>% dplyr::ungroup()

  nrep_s4U <- nrep_s4U$nreps

  # 2) Adjust -s4U replicate ID accordingly
  rep_list <- rep_list %>%
    dplyr::group_by(ctl, Exp_ID) %>%
    dplyr::mutate(r_id = ifelse(ctl == 0, ((r_id - 1) %% nrep_s4U[Exp_ID]) + 1, r_id)) %>%
    dplyr::ungroup()


  rep_list <- rep_list$r_id

  # Create mut and reps dictionary
  ID_dict <- data.frame(sample = rownames(metadf),
                        reps = rep_list,
                        mut = mut_list,
                        type = type_list)


  # Add replicate ID and s4U treatment status to metadf
  metadf <- metadf[samp_list, ] %>% dplyr::mutate(ctl = ifelse(tl == 0, 0, 1)) %>%
    dplyr::group_by(ctl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup()

  names(type_list) <- samp_list
  names(mut_list) <- samp_list
  names(rep_list) <- samp_list

  # Make vector of number of replicates of each condition
  nreps <- rep(0, times = max(mut_list))
  for(i in 1:max(mut_list)){
    nreps[i] <- max(rep_list[mut_list == i & type_list == 1])
  }

  # Get reliable features:
  message("Finding reliable Features")
  reliables <- bakR::reliableFeatures(obj, high_p = high_p, totcut = totcut, Ucut = Ucut, AvgU = AvgU)

  if(is.null(FOI)){
    keep <- reliables
  }else{
    if(concat == TRUE){
      # FOI still need to pass filtering to make sure bakRFit doesn't break
      min_reliables <- bakR::reliableFeatures(obj, high_p = 1, totcut = 1, Ucut = Ucut, AvgU = AvgU)
      keep <- unique(c(intersect(FOI, min_reliables), reliables))
    }else{
      min_reliables <- bakR::reliableFeatures(obj, high_p = 1, totcut = 1, Ucut = Ucut, AvgU = AvgU)
      keep <- intersect(FOI, min_reliables)
    }
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
  cB <- data.table::setDT(cB[cB$XF %in% ranked_features_df$XF, ])

  sdf_U <- cB[, .(n = sum(n)), by = .(sample, XF, TC, nT)]

  cB <- dplyr::as_tibble(cB)


  slist = samp_list
  tlist = type_list
  mlist = mut_list
  rlist = rep_list


  colnames(sdf_U) <- c("sample", "XF", "TC", "nT", "n")

  sdf_U <- dplyr::as_tibble(sdf_U) %>%
    dplyr::right_join(ranked_features_df, by = 'XF') %>%
    dplyr::ungroup()


  kp = keep

  ## Add sample characteristic details to U-content data frame

  df_U <- sdf_U %>%
    dplyr::ungroup() %>%
    dplyr::filter(sample %in% slist)

  df_U <- dplyr::left_join(df_U, ID_dict, by = "sample")

  sample_lookup <- df_U[df_U$type == 1, c("sample", "mut", "reps")] %>% dplyr::distinct()


  ## Calculate global average U-content so that feature-specific difference
  ## from average can be calculated

  # df_global_U <- df_U[df_U$type == 1, ] %>% dplyr::group_by(reps, mut) %>%
  #   dplyr::summarise(tot_avg_Us = sum(nT*n)/sum(n)) %>% dplyr::ungroup()
  #
  # df_feature_U <- df_U[df_U$type == 1, ] %>% dplyr::group_by(reps, mut, fnum) %>%
  #   dplyr::summarise(feature_avg_Us = sum(nT*n)/sum(n)) %>% dplyr::ungroup()
  #
  # df_U_tot <- dplyr::left_join(df_global_U, df_feature_U, by = c("mut", "reps"))
  #
  # # U_factor is log-fold difference in feature specific U-content from global average
  #   # Used in Stan model to properly adjust population average Poisson mutation rates
  # df_U_tot <- df_U_tot %>% dplyr::mutate(U_factor = log(feature_avg_Us/tot_avg_Us)) %>%
  #   dplyr::select(mut, reps, fnum, U_factor)

  df_U_tot <- df_U[df_U$type == 1, ] %>% dplyr::group_by(reps, mut) %>%
    dplyr::summarise(tot_avg_Us = sum(nT*n)/sum(n)) %>% dplyr::ungroup()


  if(Stan){

    # Filter out unreliable features and assign feature ID to each feature
    sdf <- cB %>%
      dplyr::ungroup() %>%
      dplyr::group_by(sample, XF, TC) %>%
      dplyr::summarise(Ucont = sum(nT*n)/sum(n),
                       n = sum(n)) %>%
      dplyr::right_join(ranked_features_df, by = 'XF') %>% dplyr::ungroup()


    ### Curate data for 'Stan' models

    d = sdf


    ## Add sample characteristic info to data frame

    df <- d %>%
      dplyr::ungroup() %>%
      dplyr::filter(sample %in% slist)

    df <- dplyr::left_join(df, ID_dict, by = "sample")


    ## Remove any unnecessary columns
    df <- df  %>%
      dplyr::group_by(fnum) %>%
      dplyr::arrange(type, .by_group = TRUE)



    # Add U-content information
    df <- dplyr::left_join(df, df_U_tot, by = c("mut", "reps")) %>%
      dplyr::mutate(U_factor = log(Ucont/tot_avg_Us))

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
    Avg_Counts <- df %>% dplyr::ungroup() %>% dplyr::group_by(fnum, mut, reps) %>%
      dplyr::summarise(Avg_Reads = sum(n)) %>%
      dplyr::group_by(fnum, mut) %>%
      dplyr::summarise(Avg_Reads = mean(Avg_Reads)) %>%
      dplyr::ungroup()


    tls <-rep(0, times = nMT)

    # Calculate average read counts on log10 and natural scales
      # log10 scale read counts used in 'Stan' model
      # natural scale read counts used in plotting function (plotMA())
    Avg_Counts <- Avg_Counts[order(Avg_Counts$mut, Avg_Counts$fnum),]


    Avg_Reads <- matrix(log10(Avg_Counts$Avg_Reads), ncol = nMT, nrow = NF)
    Avg_Reads_natural <- matrix(Avg_Counts$Avg_Reads, ncol = nMT, nrow = NF)


    # standardize
    Avg_Reads <- scale(Avg_Reads)

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

#' Curate data in bakRFnData object for statistical modeling
#'
#' \code{fn_process} creates the data structures necessary to analyze nucleotide recoding RNA-seq data with the
#' MLE and Hybrid implementations in \code{bakRFit}. The input to \code{fn_process} must be an object of class
#' `bakRFnData`. 
#'
#' \code{fn_process} first filters out features with less than totcut reads in any sample. It then
#' creates the necessary data structures for analysis with \code{bakRFit} and some of the visualization
#' functions (namely \code{plotMA}).
#'
#' If FOI is non-null and concat == TRUE, then all the features making it past read count filtering will be included
#' in the output. This is the same behavior as if FOI is null. If FOI is non-null and concat == FALSE, the features 
#' listed in FOI will be the only reliable features that make it past filtering. NOTE: FOIs must be deemed 
#' "reliable" (in this case that means making it past read count filtering) to make it past filtering. This is
#' because bakRFit will break otherwise.
#' If FOI is null and concat == FALSE or TRUE, then only the features making it past read count filtering
#' will be kept.
#'
#'
#' @param obj An object of class bakRFnData
#' @param totcut Numeric; Any transcripts with less than this number of sequencing reads in any sample are filtered out
#' @param FOI Features of interest; character vector containing names of features to analyze. 
#' If FOI is non-null and concat == TRUE, the features listed in FOI will be included in the list of reliable features that make it past
#' filtering. If FOI is non-null and concat == FALSE, the features listed in FOI will be the only reliable features that make it past filtering.
#' If FOI is null and concat == FALSE or TRUE, then only the features making it past read count filtering
#' will be kept.
#' @param concat Boolean; If TRUE, FOI is concatenated with output of reliableFeatures
#' @param Chase Boolean; if TRUE, pulse-chase analysis strategy is implemented
#' @return returns list of objects that can be passed to \code{TL_stan} and/or \code{fast_analysis}. Those objects are:
#' \itemize{
#'  \item Stan_data; list that can be passed to \code{TL_stan} with Hybrid_Fit = TRUE. Consists of metadata as well as data that
#'  `Stan` will analyze. Data to be analyzed consists of equal length vectors. The contents of Stan_data are:
#'  \itemize{
#'   \item NE; Number of datapoints for 'Stan' to analyze (NE = Number of Elements)
#'   \item NF; Number of features in dataset
#'   \item TP; Numerical indicator of s4U feed (0 = no s4U feed, 1 = s4U fed)
#'   \item FE; Numerical indicator of feature
#'   \item num_mut; Number of U-to-C mutations observed in a particular set of reads
#'   \item MT; Numerical indicator of experimental condition (Exp_ID from metadf)
#'   \item nMT; Number of experimental conditions
#'   \item R; Numerical indicator of replicate
#'   \item nrep; Number of replicates (maximum across experimental conditions)
#'   \item nrep_vect; Vector of number of replicates in each experimental condition
#'   \item tl; Vector of label times for each experimental condition
#'   \item Avg_Reads; Standardized log10(average read counts) for a particular feature in a particular condition, averaged over
#'   replicates
#'   \item sdf; Dataframe that maps numerical feature ID to original feature name. Also has read depth information
#'   \item sample_lookup; Lookup table relating MT and R to the original sample name
#'  }
#'  \item Fn_est; A data frame containing fraction new estimates:
#'  \itemize{
#'   \item sample; Original sample name
#'   \item XF; Original feature name
#'   \item fn; Fraction new estimate
#'   \item n; Number of reads
#'   \item Feature_ID; Numerical ID for each feature
#'   \item Replicate; Numerical ID for each replicate
#'   \item Exp_ID; Numerical ID for each experimental condition
#'   \item tl; s4U label time
#'   \item logit_fn; logit of fraction new estimate
#'   \item kdeg; degradation rate constant estimate
#'   \item log_kdeg; log of degradation rate constant estimate
#'   \item logit_fn_se; Uncertainty of logit(fraction new) estimate
#'   \item log_kd_se; Uncertainty of log(kdeg) estimate
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
fn_process <- function(obj, totcut = 50, Chase = FALSE, FOI = c(), concat = TRUE){
  
  ## Check obj
  if(!inherits(obj, "bakRFnData")){
    stop("obj must be of class bakRFnData")
  }
  

  ## Check totcut
  if(!is.numeric(totcut)){
    stop("totcut must be numeric")
  }else if( totcut < 0 ){
    stop("totcut must be greater than 0")
  }else if(totcut > 5000){
    warning("totcut is abnormally high (> 5000); many features will not have this much coverage in every sample and thus get filtered out.")
  }
  
  ## Check FOI
  if(!is.null(FOI)){
    if(typeof(obj$fns$XF) != typeof(FOI)){
      warning("FOI should be the same data type as fns$XF in the bakRFnData object; if it is not none of the feature of interest will be found
            in the fn data frame.")
    }
  }
  
  # Define helper functions:
  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))
  
  
  # Bind to NULL
  tl <- ctl <- Exp_ID <- r_id <- XF <- n <- npass <- NULL
  Feature_ID <- fn <- var <- global_mean <- global_var <- alpha_p <- beta_p <- NULL

  metadf <- obj$metadf
  fns <- obj$fns
  
  message("Mapping sample name to sample characteristics")
  
  samp_list <- unique(fns$sample)
  
  c_list <- rownames(metadf[metadf$tl == 0,])
  
  s4U_list <- samp_list[!(samp_list %in% c_list)]
  
  type_list <- ifelse(metadf[samp_list, "tl"] == 0, 0, 1)
  mut_list <- metadf[samp_list, "Exp_ID"]
  
  # If number of -s4U controls > +s4U controls, wrap R_ID for -s4Us
  # So if there are 3 -s4U replicates and 2 + s4U, -s4U R_ID should be: 1, 2, 1
  # Let R_raw but R_ID for -s4U calculated as with +s4U samples (so 1, 2, 3 in above example)
  # Then R_ID = ((R_raw - 1) %% nreps_+s4U) + 1, where nreps_+s4U is number of replicates in +s4U sample
  rep_list <- metadf[samp_list,] %>% dplyr::mutate(ctl = ifelse(tl == 0, 0, 1)) %>%
    dplyr::group_by(ctl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup()
  
  ### Calculate -s4U adjusted r_id
  # 1) Determine number of +s4U replicates in each sample
  nrep_s4U <- rep_list %>% dplyr::filter(ctl == 1) %>%
    dplyr::group_by(Exp_ID) %>%
    dplyr::summarise(nreps = max(r_id)) %>% dplyr::ungroup()
  
  nrep_s4U <- nrep_s4U$nreps
  
  # 2) Adjust -s4U replicate ID accordingly
  rep_list <- rep_list %>%
    dplyr::group_by(ctl, Exp_ID) %>%
    dplyr::mutate(r_id = ifelse(ctl == 0, ((r_id - 1) %% nrep_s4U[Exp_ID]) + 1, r_id)) %>%
    dplyr::ungroup()
  
  
  rep_list <- rep_list$r_id
  
  # Create mut and reps dictionary
  ID_dict <- data.frame(sample = rownames(metadf),
                        Replicate = rep_list,
                        Exp_ID = mut_list,
                        Type = type_list)
  
  
  # Add replicate ID and s4U treatment status to metadf
  metadf <- metadf[samp_list, ] %>% dplyr::mutate(ctl = ifelse(tl == 0, 0, 1)) %>%
    dplyr::group_by(ctl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup()
  
  names(type_list) <- samp_list
  names(mut_list) <- samp_list
  names(rep_list) <- samp_list
  
  # Make vector of number of replicates of each condition
  nreps <- rep(0, times = max(mut_list))
  for(i in 1:max(mut_list)){
    nreps[i] <- max(rep_list[mut_list == i & type_list == 1])
  }
  
  
  # Find features with above the threshold in all samples
  nsamps <- nrow(metadf)
  
  message("Filtering out low coverage features")

  # Read count filter
  keep <- fns %>%
    dplyr::group_by(XF) %>%
    dplyr::summarise(npass = sum(n >= totcut)) %>%
    dplyr::filter(npass == nsamps) %>%
    dplyr::select(XF)
  
  if(is.null(FOI)){
  
    keep <- keep$XF
    
  }else{
    
    if(concat){
      
      # FOI can only be kept if it has at least 1 read in all samples
      min_keep <- fns %>%
        dplyr::group_by(XF) %>%
        dplyr::summarise(npass = sum(n >= 1)) %>%
        dplyr::filter(npass == nsamps) %>%
        dplyr::select(XF)

      keep <- unique(c(intersect(FOI, min_keep$XF), keep$XF))

    }else{
      
      # FOI can only be kept if it has at least 1 read in all samples
      min_keep <- fns %>%
        dplyr::group_by(XF) %>%
        dplyr::summarise(npass = sum(n >= 1)) %>%
        dplyr::filter(npass == nsamps) %>%
        dplyr::select(XF)
      
      keep <- intersect(FOI, min_keep$XF)
    }
  }

  
  fns <- fns[fns$XF %in% keep,]
  
  # Map each reliable feature to a numerical feature ID (fnum)
  ranked_features_df  <- fns %>%
    dplyr::ungroup() %>%
    dplyr::select(XF) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Feature_ID = order(XF)) %>%
    dplyr::arrange(Feature_ID) %>%
    dplyr::select(XF, Feature_ID)
  
  message("Processing data...")
  
  # Make count matrix that is DESeq2 compatible
  Cnt_mat <- matrix(0, ncol = length(samp_list), nrow = length(unique(keep)))
  
  # Order by XF to make sure XFs together in fns
  fns <- fns[order(fns$XF),]
  
  for(s in seq_along(samp_list)){
    Cnt_mat[,s] <- fns$n[fns$sample == samp_list[s]]
  }
  
  rownames(Cnt_mat) <- fns$XF[fns$sample == samp_list[1]]
  colnames(Cnt_mat) <- samp_list
  
  
  # Add feature_ID and sample characteristics
  fns <- dplyr::inner_join(fns, ranked_features_df, by = "XF")
  fns <- dplyr::inner_join(fns, ID_dict, by = "sample")
  
  # Add label time
  tl_df <- metadf[metadf$tl > 0,c("tl", "Exp_ID", "r_id")]
  colnames(tl_df) <- c("tl", "Exp_ID", "Replicate")
  
  fns <- dplyr::inner_join(fns, tl_df, by = c("Exp_ID", "Replicate"))
  
  # Remove -s4U data from fns
  fn_ctl_data <- fns[fns$Type == 0, colnames(fns) != "Type"]
  fn_ctl_data$tl <- 0
  
  fns <- fns[fns$Type == 1,]
  
  fns <- fns[,colnames(fns) != "Type"]
  
  # Convert fn to various reparameterizations of interest
  fns$logit_fn <- logit(fns$fn)
  fns$kdeg <- -log(1 - fns$fn)/fns$tl
  fns$log_kdeg <- log(fns$kdeg)
  
  
  ### Estimate uncertainty
  
  ## Procedure:
  # 1) Estimate beta distribution prior empirically
  # 2) Determine beta posterior
  # 3) Use beta sd as uncertainty
  
  ## Estimate beta prior
  fn_means <- fns %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(global_mean = mean(fn),
                     global_var = stats::var(fn))
  
  priors <- fn_means %>%
    dplyr::mutate(alpha_p = global_mean*(((global_mean*(1-global_mean))/global_var) - 1),
                  beta_p = alpha_p*(1 - global_mean)/global_mean)
  
  ## Add priors
  fns <- dplyr::inner_join(fns, priors, by = "sample")
  
  ## Estimate logit uncertainty
  lfn_calc <- function(alpha, beta){
    
    EX <- alpha/(alpha + beta)
    VX <- (alpha*beta)/(((alpha + beta)^2)*(alpha + beta + 1))
    
    
    totvar <- (((1/EX) + 1/(1 - EX))^2)*VX
    
    return(totvar)
  }
  
  ## Estimate log(kdeg) uncertainty
  lkdeg_calc <- function(alpha, beta){
    
    EX <- alpha/(alpha + beta)
    VX <- (alpha*beta)/(((alpha + beta)^2)*(alpha + beta + 1))
    
    
    totvar <- (( 1/(log(1-EX)*(1-EX)) )^2)*VX
    
    return(totvar)
  }
  
  fns <- fns %>%
    dplyr::mutate(logit_fn_se = sqrt(lfn_calc(alpha_p + n*fn, n + beta_p)),
                  log_kd_se = sqrt(lkdeg_calc(alpha_p + n*fn, n + beta_p)))
  
  fns <- fns[,!(colnames(fns) %in% c("alpha_p", "beta_p", "global_mean", "global_var"))]
  
  
  ### Average standardized reads
  
  # Dataset characteristics
  nMT <- max(metadf$Exp_ID)
  NF <- max(fns$Feature_ID)
  
  ## Calculate Avg. Read Counts
  Avg_Counts <- fns %>% dplyr::ungroup() %>%
    dplyr::group_by(Feature_ID, Exp_ID) %>%
    dplyr::summarise(Avg_Reads = mean(n), .dplyr.summarise.inform = FALSE) %>%
    dplyr::ungroup()
  
  # Calculate average read counts on log10 and natural scales
  # log10 scale read counts used in 'Stan' model
  # natural scale read counts used in plotting function (plotMA())
  Avg_Counts <- Avg_Counts[order(Avg_Counts$Exp_ID, Avg_Counts$Feature_ID),]
  
  
  Avg_Reads <- matrix(log10(Avg_Counts$Avg_Reads), ncol = nMT, nrow = NF)
  Avg_Reads_natural <- matrix(Avg_Counts$Avg_Reads, ncol = nMT, nrow = NF)
  
  
  # standardize
  Avg_Reads <- scale(Avg_Reads)
  
  
  ### Generate input necessary for Hybrid model
  
  # s4U label time in each experimental condition
  tls <- rep(0, times = nMT)
  for(m in 1:nMT){
    tls[m] <- unique(metadf$tl[(metadf$Exp_ID == m) & (metadf$tl != 0)])
  }
  
  # Add tl info to fast_df
  tl_df <- data.frame(tl = tls,
                      mut = 1:length(tls))
  
  # Sample lookup
  colnames(ID_dict) <- c("sample", "Replicate", "Exp_ID", "Type")
  sample_lookup <- ID_dict[ID_dict$Type == 1, c("sample", "Exp_ID", "Replicate")] %>% dplyr::distinct()
  
  # Feature number lookup
  sdf <- ranked_features_df
  colnames(sdf) <- c("XF", "fnum")
  
  data_list <- list(
    NE = nrow(fns),
    NF = max(fns$Feature_ID),
    MT = fns$Exp_ID,
    FE = fns$Feature_ID,
    tl = tl_df$tl,
    logit_fn_rep = fns$logit_fn,
    fn_se = fns$logit_fn_se,
    Avg_Reads = Avg_Reads,
    Avg_Reads_natural = Avg_Reads_natural,
    nMT = max(fns$Exp_ID),
    R = fns$Replicate,
    nrep = max(fns$Replicate),
    sample_lookup = sample_lookup, 
    sdf = sdf,
    mutrates = data.frame(),
    nrep_vect = nreps,
    Chase = as.integer(Chase) 
  )
  
  out <- list(Stan_data = data_list, Count_Matrix = Cnt_mat,
                               Fn_est = dplyr::as_tibble(fns),
                               Ctl_data = dplyr::as_tibble(fn_ctl_data))
  
  return(out)
}