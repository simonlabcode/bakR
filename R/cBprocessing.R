# I can add match.arg and a string input denoting if user wants output to be for
#' Identify features (e.g., transcripts) with high quality data
#'
#' This function identifies all features (e.g., transcripts, exons, etc.) for which the mutation rate
#' is below a set threshold in the control (no s4U) sample and which have more reads than a set threshold
#' in at least one control sample
#' @param obj Object of class DynamicSeqData
#' @param high_p highest mutation rate accepted in control samples
#' @param totcut readcount cutoff
#' @importFrom magrittr %>%
#' @return list of gene names that passed reliability filter
#' @export
reliableFeatures <- function(obj,
                             high_p = 0.2,
                             totcut = 1000){

  cB <- obj$cB
  nsamps <- length(unique(cB$sample))


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
  return(y)
}

# I should give them an option about the name of the column they want to sort by
# I could even give users an option to specify column names or order
#' Extract data for Stan analysis from cB
#'
#' This function obtains the data list necessary to analyze TL-seq data with Stan.
#' @param obj An object of class DynamicSeqData
#' @param keep_input two element vector; 1st element is highest mut rate accepted in control samples, 2nd element is read count cutoff
#' @param Stan Boolean; if TRUE, then data_list that can be passed to Stan is curated
#' @param Fast Boolean; if TRUE, then dataframe that can be passed to fast_analysis() is curated
#' @param FOI Features of interest; character vector containing names of features to analyze
#' @param concat Boolean; If TRUE, FOI is concatenated with output of reliableFeatures
#' @return returns list which can be passed to Stan analyses
#' @importFrom magrittr %>%
#' @export
cBprocess <- function(obj,
                       keep_input=c(0.2, 50),
                       Stan = TRUE,
                       Fast = FALSE,
                       FOI = c(),
                       concat = TRUE){
  cB <- obj$cB
  metadf <- obj$metadf

  samp_list <-unique(cB$sample)

  c_list <- rownames(metadf[metadf$tl == 0,])

  s4U_list <- samp_list[!(samp_list %in% c_list)]

  type_list <- ifelse(metadf[samp_list, "tl"] == 0, 0, 1)
  mut_list <- metadf[samp_list, "Exp_ID"]
  rep_list <- metadf[samp_list,] %>% dplyr::group_by(tl, Exp_ID) %>% dplyr::mutate(r_id = 1:length(tl)) %>% dplyr::ungroup() %>% dplyr::select(r_id)
  rep_list <- rep_list$r_id


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
    reliables <- DynamicSeq::reliableFeatures(obj, high_p = keep_input[1], totcut = keep_input[2])
    keep <- c(FOI, reliables[!(reliables %in% FOI)])
  }else{
    keep <- FOI
  }




  # This df is created in both cases
  ranked_features_df  <- cB %>%
    dplyr::ungroup() %>%
    dplyr::filter(XF %in% keep) %>%
    dplyr::group_by(XF) %>%
    dplyr::summarize(n = sum(n)) %>%
    dplyr::mutate(fnum = order(XF)) %>%
    dplyr::arrange(fnum) %>%
    dplyr::select(XF, fnum)


  # Create count dataframe
  Counts_df <- cB %>%
    dplyr::ungroup() %>%
    dplyr::filter(XF %in% keep) %>%
    dplyr::filter(sample %in% s4U_list) %>%
    dplyr::group_by(XF, sample) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::right_join(ranked_features_df, by = 'XF') %>% dplyr::ungroup()

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
      tl = unique(metadf$tl[metadf$tl > 0])[1],
      U_cont = df$U_factor,
      sdf = sdf
    )

  }


  if(Stan & Fast){
    out <- list(data_list, df_U, Counts_df)
    names(out) <- c("Stan_data", "Fast_df", "Counts_df")
  }else if(!Stan){
    out <- list(df_U, Counts_df)
    names(out) <- c("Fast_df", "Counts_df")
  }else{
    out <- list(data_list, Counts_df)
    names(out) <- c("Stan_data", "Counts_df")
  }

  return(out)
}
