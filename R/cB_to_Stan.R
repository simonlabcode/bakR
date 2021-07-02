# I can add match.arg and a string input denoting if user wants output to be for
#' Identify features (e.g., transcripts) with high quality data
#'
#' This function identifies all features (e.g., transcripts, exons, etc.) for which the mutation rate
#' is below a set threshold in the control (no s4U) sample and which have more reads than a set threshold
#' in at least one control sample
#' @param cB cB.rds file from pipeline
#' @param c_list vector of names of control (no s4U fed) samples
#' @param high_p highest mutation rate accepted in control samples
#' @param totcut readcount cutoff
#' @importFrom magrittr %>%
#' @return list of gene names that passed reliability filter
#' @export
reliableFeatures <- function(cB,
                             c_list,
                             high_p = 0.2,
                             totcut = 1000){

  y <- cB %>%
    dplyr::ungroup() %>%
    dplyr::filter(sample %in% c_list,
                  !grepl('__', XF)) %>%
    dplyr::mutate(totTC = TC*n) %>%
    dplyr::group_by(sample, XF) %>%
    dplyr::summarize(tot_mut = sum(totTC),
                     totcounts = sum(n)) %>%
    dplyr::filter(totcounts >= totcut) %>%
    dplyr::filter(tot_mut/totcounts < high_p) %>%
    dplyr::ungroup( ) %>%
    dplyr::group_by(XF) %>%
    dplyr::summarize(counts = dplyr::n()) %>%
    dplyr::filter(counts == length(c_list)) %>%
    dplyr::select(XF) %>%
    unlist() %>%
    unique()
  return(y)
}


#' Extract data for Stan analysis from cB
#'
#' This function obtains the data list necessary to analyze TL-seq data with Stan.
#' @param cB_raw cB file generated from TL-seq pipeline
#' @param samp_list vector of names of samples
#' @param type_list vector with 1 entry per sample; 0 = no s4U, 1 = s4U fed
#' @param mut_list vector with 1 entry per sample; 1 = reference, > 1 = different experimental conditions (e.g., KO of protein X)
#' @param rep_list vector with 1 entry per sample that indexes replicate; 1 = 1st replicate, 2 = 2nd replicate, etc.
#' @param tl single numerical value; s4U label time used in s4U fed samples
#' @param keep_input two element vector; 1st element is highest mut rate accepted in control samples, 2nd element is read count cutoff
#' @param FOI Features of interest; character vector containing names of features to analyze
#' @param concat Boolean; If TRUE, FOI is concatenated with output of reliableFeatures
#' @return returns list which can be passed to Stan analyses
#' @importFrom magrittr %>%
#' @export
cBtoStan <- function(cB_raw,
                       samp_list,
                       type_list,
                       mut_list,
                       rep_list,
                       tl,
                       keep_input=c(0.2, 50),
                       FOI = c(),
                       concat = TRUE){
  cB <- cB_raw %>%
    dplyr::select(sample, XF, GF, TC, n, io)

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
    dplyr::group_by(sample, XF, TC) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::right_join(ranked_features_df, by = 'XF')

  ##### Run Stan model ---------

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


  df <- df  %>%
    dplyr::group_by(XF, fnum, type, mut, TC, reps) %>%
    dplyr::summarise(n = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(fnum) %>%
    dplyr::arrange(type, .by_group = TRUE)


  FE = df$fnum
  NE <- dim(df)[1]
  NF <- length(kp)
  TP <- df$type
  MT <- df$mut
  nMT <- length(unique(MT))
  R <- df$reps

  num_mut <- df$TC # Number of mutations
  num_obs <- df$n # Number of times observed

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
    tl = tl
  )

  return(data_list)
}
