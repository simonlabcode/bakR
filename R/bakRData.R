#' bakRData object constructor for internal use
#'
#' This function efficiently creates an object of class bakRData
#' without performing rigorous checks
#' @param cB Dataframe with columns corresponding to feature ID, number of Ts, number of mutations, sample ID, and number of identical observations
#' @param metadf Dataframe detailing s4U label time and experimental ID of each sample
new_bakRData <- function(cB, metadf){
  stopifnot(is.data.frame(cB))
  stopifnot(is.data.frame(metadf))
  structure(list(cB = cB, metadf = metadf), class = "bakRData")
}

#' bakR Data object validator
#'
#' This functions ensures that input for bakRData object construction is valid
#' @importFrom magrittr %>%
#' @param obj An object of class bakRData
validate_bakRData <- function(obj){

  tl <- Exp_ID <- check <- NULL

  vals <- unclass(obj)
  cB <- vals$cB
  metadf <- vals$metadf

  ## Check if rownames of metadf are same as unique(cB$sample)
  if("sample" %in% colnames(cB)){
    if(!identical(rownames(metadf), unique(cB$sample))){
      stop(
        "Row names of metadf are not same as unique(cB$sample).
        Make sure the order in which samples appear in cB
        matches the rownames of metadf",
        call. = FALSE
      )
    }
  }else{
    if(!identical(rownames(metadf), unique(cB[,2]))){
      stop(
        "Row names of metadf are not same as unique(cB[,2]).
        Make sure 2nd column of cB contains sample names and
        that these match the rownames of metadf",
        call. = FALSE
      )
    }
  }


  ## Check if Numerical data in cB is all integer
  if(sum(purrr::map_dbl(unclass(cB), is.numeric)) < 3){
    stop(
      "There are less than 3 columns of cB containing numeric data.
      cB should have a column corresponding to numbers of mutations (TC),
      numbers of Ts in a read (nT), and number of identical observations (n),
      all of which should be integer data.",
      call. = FALSE
    )
  }else if(sum(purrr::map_dbl(unclass(cB), is.integer)) < 3){
    cB$TC <- as.integer(cB$TC)
    cB$nT <- as.integer(cB$nT)
    cB$n <- as.integer(cB$n)
  }


  ## Check that all data in metadf is numerical

  if(sum(!purrr::map_dbl(unclass(metadf), is.numeric)) > 0 ){
    stop(
      "Not all columns of metadf contain numeric data. Metadf should contain 2 columns,
      one corresponding to numerical IDs of experimental conditions for each sample
      (Exp_ID) and one corresponding to s4U label times (tl; 0 if an unlabeled control)",
      call. = FALSE
    )
  }

  ## Check that one column in metadf is integer

  if(sum(purrr::map_dbl(unclass(metadf), is.integer)) < 1){
    metadf$Exp_ID <- as.integer(metadf$Exp_ID)
  }

  ## Check to make sure n > 0 in cB

  if("n" %in% colnames(cB)){
    if(sum(cB$n < 1) > 0){
      stop(
        "n column of cB contains values < 1. This column should tally the number of identical
        observations for each combination of gene ID, mutational and T content, and sample ID."
      )
    }
  }else if(sum(cB[,5] < 1) > 0){
    stop(
      "5th column of cB contains values < 1. This column should tally the number of identical
        observations for each combination of gene ID, mutational and T content, and sample ID."
    )
  }

  ## Check to make sure TC >= 0 in cB

  if("TC" %in% colnames(cB)){
    if(sum(cB$TC < 0) > 0){
      stop(
        "TC column of cB contains values < 0. This column should tally the number of mutations
        in each read and thus should be >= 0."
      )
    }
  }else if(sum(cB[,3] < 0) > 0){
    stop(
      "3rd column of cB contains values < 0. This column should tally the number of mutations
        in each read and thus should be >= 0."
    )
  }

  ## Check to make sure that nT >= 0 in cB
  if("nT" %in% colnames(cB)){
    if(sum(cB$nT < 0) > 0){
      stop(
        "nT column of cB contains values < 0. This column should tally the number of Ts (Us in RNA)
        in each read and thus should be >= 0."
      )
    }
  }else if(sum(cB[,4] < 0) > 0){
    stop(
      "4th column of cB contains values < 0. This column should tally the number of Ts (Us in RNA)
        in each read and thus should be >= 0."
    )
  }

  ## Check to make sure that all Exp_IDs have at least 2 s4U samples
  nMT <- length(unique(metadf$Exp_ID))
  replicate_check <- metadf %>%
    dplyr::filter(tl > 0) %>%
    dplyr::mutate(check = 1) %>%
    dplyr::group_by(Exp_ID) %>%
    dplyr::summarise(check = sum(check))

  if(any(replicate_check$check < 2) | (nrow(replicate_check) != nMT)){
    stop("All experimental IDs must have at least 2 s4U fed replicates.")
  }

  ## Check that Exp_ID ranges from 1 to the max Exp_ID, hitting all integers in between
  Exp_min <- min(replicate_check$Exp_ID)
  Exp_max <- max(replicate_check$Exp_ID)
  Exp_ID_actual <- replicate_check$Exp_ID


  if(Exp_min != 1){
    stop("The minimum Exp_ID in metadf must be 1!")
  }else if(!all(Exp_min:Exp_max %in% Exp_ID_actual)){ # Hard to see how this could ever be the case
    stop("Exp_ID must contain all integers from 1 to the max Exp_ID.")
  }else if(!all(Exp_ID_actual %in% Exp_min:Exp_max)){
    stop("Exp_ID must contain all integers from 1 to the max Exp_ID. Did you skip an integer?")
  }

  obj

}


#' bakR Data object helper function for users
#'
#' This function creates an object of class bakRData
#' @param cB Dataframe with columns corresponding to feature ID, number of Ts, number of mutations, sample ID, and number of identical observations
#' @param metadf Dataframe detailing s4U label time and experimental ID of each sample
#' @return A bakRData object. This has two components: a data frame describing experimental
#' details (metadf) and a data frame containing the NR-seq data (cB).
#' @examples
#' # Load cB
#' data("cB_small")
#'
#' # Load metadf
#' data("metadf")
#'
#' # Create bakRData object
#' bakRData <- bakRData(cB_small, metadf)
#'
#' @export
bakRData <- function(cB, metadf){

  cB <- as.data.frame(cB)
  metadf <- as.data.frame(metadf)

  cB_cols <- c("XF", "sample", "TC", "nT", "n")
  meta_cols <- c("tl", "Exp_ID")

  ## Check number of columns of CB
  if(ncol(cB) < 5){
    stop(
      "There are less than 5 columns in cB. cB should contain at least 5 columns,
      one corresponding to gene names (XF), one to sample names (sample),
      one to number of T to C mutations (TC), one to number of Ts in read (nT),
      and one to the number of identical observations (n).",
      call. = FALSE
    )
  }

  if(ncol(metadf) < 2){
    stop(
      "There are less than 2 columns in metadf. Metadf should contain at least 2 columns,
      one corresponding to numerical IDs of experimental conditions for each sample
      (Exp_ID) and one corresponding to s4U label times (tl; 0 if an unlabeled control)",
      call. = FALSE
    )
  }

  ## Add column names to cB and metadf as necessary
  if(!all(cB_cols %in% colnames(cB))){
    colnames(cB)[1:5] <- cB_cols
    warning("Renamed first 5 columns of cB to XF, sample, TC, nT, and n.
    If this does not reflect the content of those columns, properly rearrange cB columns and rerun bakRData().")
  }

  if(!all(meta_cols %in% colnames(metadf))){
    colnames(metadf)[1:2] <- meta_cols
    warning("Renamed first 2 columns of metadf to tl and Exp_ID.
    If this does not reflect the content of those columns, properly rearrange metadf columns and rerun bakRData().")
  }

  if(max(metadf$Exp_ID) == 1){
    warning("You only have data from one experimental condition. bakR is designed for comparing
    two or more sets of samples, and thus some basic bakR use-cases will fail with your
    data. If you are hoping to just estimate kinetic parmaters for your one set of samples,
    then bakRFit(bakRData, StanRateEst = FALSE) will work. No other model implementaiton
    (i.e., with any of StanRateEst, StanFit, or HybridFit set to TRUE) will work with this kind
    of data.")
  }

  validate_bakRData(new_bakRData(cB, metadf))


}


#' bakRFnData object constructor for internal use
#'
#' This function efficiently creates an object of class bakRFnData
#' without performing rigorous checks
#' @param fns Dataframe with columns corresponding to sample names (sample), feature IDs (XF),
#' fraction new estimates (fn), and number of sequencing reads (nreads) 
#' @param metadf Dataframe detailing s4U label time and experimental ID of each sample
new_bakRFnData <- function(fns, metadf){
  stopifnot(is.data.frame(fns))
  stopifnot(is.data.frame(metadf))
  structure(list(fns = fns, metadf = metadf), class = "bakRFnData")
}

#' bakRFnData object validator
#'
#' This functions ensures that input for bakRFnData object construction is valid
#' @importFrom magrittr %>%
#' @param obj An object of class bakRFnData
validate_bakRFnData <- function(obj){
  
  tl <- Exp_ID <- check <- NULL
  
  vals <- unclass(obj)
  fns <- vals$fns
  metadf <- vals$metadf
  
  if("se" %in% colnames(fns)){
    se_provided <- TRUE
  }else{
    se_provided <- FALSE
  }
  
  ## Check if rownames of metadf are same as unique(cB$sample)
  if("sample" %in% colnames(fns)){
    if(!identical(rownames(metadf), unique(fns$sample))){
      stop(
        "Row names of metadf are not same as unique(fns$sample).
        Make sure the order in which samples appear in fns
        matches the rownames of metadf",
        call. = FALSE
      )
    }
  }else{
    if(!identical(rownames(metadf), unique(fns[,2]))){
      stop(
        "Row names of metadf are not same as unique(fns[,2]).
        Make sure 2nd column of fns contains sample names and
        that these match the rownames of metadf",
        call. = FALSE
      )
    }
  }
  
  
  ## Check if Numerical data in fns is all integer
  if(sum(purrr::map_dbl(unclass(fns), is.numeric)) < ifelse(se_provided, 3, 2) ){
    
    if(se_provided){
      stop(
        "There are less than 3 columns of fns containing numeric data.
      fns should have a column corresponding to the fraction new estimates (fn),
      fn estimate uncertainties (se), and number of sequencing reads (n).",
        call. = FALSE
      )
    }else{
      stop(
        "There are less than 2 columns of fns containing numeric data.
      fns should have a column corresponding to the fraction new estimates (fn),
      and number of sequencing reads (n).",
        call. = FALSE
      )
    }

  }
  
  
  ## Check that all data in metadf is numerical
  
  if(sum(!purrr::map_dbl(unclass(metadf), is.numeric)) > 0 ){
    stop(
      "Not all columns of metadf contain numeric data. Metadf should contain 2 columns,
      one corresponding to numerical IDs of experimental conditions for each sample
      (Exp_ID) and one corresponding to s4U label times (tl)",
      call. = FALSE
    )
  }
  
  ## Check that one column in metadf is integer
  
  if(sum(purrr::map_dbl(unclass(metadf), is.integer)) < 1){
    metadf$Exp_ID <- as.integer(metadf$Exp_ID)
  }
  
  ## Check to make sure n > 0 in fns
  
  if("n" %in% colnames(fns)){
    if(sum(fns$n < 0) > 0){
      stop(
        "n column of fns contains negative values. This column should tally the number of reads 
        mapping to the corresponding feature in the corresponding sample."
      )
    }
  }else if(sum(fns[,ifelse(se_provided, 5, 4)] < 0) > 0){
    
    if(se_provided){
      stop(
        "5th column of fns contains values < 0. This column should tally the number of reads 
        mapping to the corresponding feature in the corresponding sample."
      )
    }else{
      stop(
        "4th column of fns contains values < 0. This column should tally the number of reads 
        mapping to the corresponding feature in the corresponding sample."
      )
    }

  }
  
  
  ## Check to make sure that fn between 0 and 1 in fns
  if("fn" %in% colnames(fns)){
    if(sum(fns$fn < 0 | fns$fn > 1) > 0){
      stop(
        "fn column of fns contains values < 0 or > 1. This column should represent
        the fraction new estimate for the corresponding feature in the corresponding
        sample, and thus should be betweeen 0 and 1."
      )
    }
  }else if(sum(fns[,3] < 0 | fns[,3] > 1) > 0){
    stop(
      "3rd column of fns contains values < 0 or > 1. This column should represent
        the fraction new estimate for the corresponding feature in the corresponding
        sample, and thus should be betweeen 0 and 1."
    )
  }
  
  ## Check to make sure that fn se is between 0 and 1
  if("se" %in% colnames(fns)){
    if(sum(fns$se < 0 | fns$se > 1) > 0){
      stop(
        "se column of fns contains values < 0 or > 1. This column should represent
        the fraction new estimate uncertainty for the corresponding feature in the corresponding
        sample, and thus should be betweeen 0 and 1."
      )
    }
  }else if(sum(fns[,4] < 0 | fns[,4] > 1) > 0){
    stop(
      "4th column of fns contains values < 0 or > 1. This column should represent
        the fraction new estimate uncertainty for the corresponding feature in the corresponding
        sample, and thus should be betweeen 0 and 1."
    )
  }
  
  ## Check to make sure that there is no negative tl
  if(sum(metadf$tl < 0) > 1){
    stop(
      "tl is < 0 for some entries of metadf. tl represents the duration of s4U labeling
      and thus is strictly positive or 0 (for -s4U controls)."
    )
  }
  
  ## Check to make sure that all Exp_IDs have at least 2 s4U samples
  nMT <- length(unique(metadf$Exp_ID))
  replicate_check <- metadf %>%
    dplyr::filter(tl > 0) %>%
    dplyr::mutate(check = 1) %>%
    dplyr::group_by(Exp_ID) %>%
    dplyr::summarise(check = sum(check))
  
  if(any(replicate_check$check < 2) | (nrow(replicate_check) != nMT)){
    stop("All experimental IDs must have at least 2 s4U fed replicates.")
  }
  
  ## Check that Exp_ID ranges from 1 to the max Exp_ID, hitting all integers in between
  Exp_min <- min(replicate_check$Exp_ID)
  Exp_max <- max(replicate_check$Exp_ID)
  Exp_ID_actual <- replicate_check$Exp_ID
  
  
  if(Exp_min != 1){
    stop("The minimum Exp_ID in metadf must be 1!")
  }else if(!all(Exp_min:Exp_max %in% Exp_ID_actual)){ # Hard to see how this could ever be the case
    stop("Exp_ID must contain all integers from 1 to the max Exp_ID.")
  }else if(!all(Exp_ID_actual %in% Exp_min:Exp_max)){
    stop("Exp_ID must contain all integers from 1 to the max Exp_ID. Did you skip an integer?")
  }
  
  obj
  
}


#' bakRFnData object helper function for users
#'
#' This function creates an object of class bakRFnData
#' @param fns Dataframe with columns corresponding to sample names (sample), feature IDs (XF),
#' fraction new estimates (fn), and number of sequencing reads (nreads) 
#' @param metadf Dataframe detailing s4U label time and experimental ID of each sample
#' @return A bakRFnData object. This has two components: a data frame describing experimental
#' details (metadf) and a data frame containing the fraction new estimates (fns).
#' @examples
#' ### NEED TO ADD EXAMPLE DATA
#' # Load cB
#' data("cB_small")
#'
#' # Load metadf
#' data("metadf")
#'
#' # Create bakRData object
#' bakRData <- bakRData(cB_small, metadf)
#'
#' @export
bakRFnData <- function(fns, metadf){
  

  fns <- as.data.frame(fns)
  metadf <- as.data.frame(metadf)
  
  if("se" %in% colnames(fns)){
    fns_cols <- c("XF", "sample", "fn", "se", "n")
    se_provided <- TRUE
    
  }else{
    fns_cols <- c("XF", "sample", "fn", "n")
    se_provided <- FALSE
    
  }
  
  meta_cols <- c("tl", "Exp_ID")
  
  ## Check number of columns of fns
  if(ncol(fns) < 4){
    stop(
      "There are less than 4 columns in fns. fns should contain at least 4 columns,
      one corresponding to feature names (XF), one to sample names (sample),
      one to the fraction new estimate (fn), and one to the number of sequencing reads
      mapping to XF.",
      call. = FALSE
    )
  }
  
  if(ncol(metadf) < 2){
    stop(
      "There are less than 2 columns in metadf. Metadf should contain at least 2 columns,
      one corresponding to numerical IDs of experimental conditions for each sample
      (Exp_ID) and one corresponding to s4U label times (tl; 0 if an unlabeled control)",
      call. = FALSE
    )
  }
  
  ## Add column names to fns and metadf as necessary
  if(!all(fns_cols %in% colnames(fns))){
    
    if(se_provided){
      colnames(fns)[1:5] <- fns_cols
      warning("Renamed first 5 columns of fns to XF, sample, fn, se, and n.
      If this does not reflect the content of those columns, properly rearrange fns columns and rerun bakRFnData().")
    }else{
      colnames(fns)[1:4] <- fns_cols
      warning("Renamed first 4 columns of fns to XF, sample, fn, and n.
      If this does not reflect the content of those columns, properly rearrange fns columns and rerun bakRFnData().")
    }

  }
  
  if(!all(meta_cols %in% colnames(metadf))){
    colnames(metadf)[1:2] <- meta_cols
    warning("Renamed first 2 columns of metadf to tl and Exp_ID.
    If this does not reflect the content of those columns, properly rearrange metadf columns and rerun bakRFnData().")
  }
  
  if(max(metadf$Exp_ID) == 1){
    warning("You only have data from one experimental condition. bakR is designed for comparing
    two or more sets of samples, and thus some basic bakR use-cases will fail with your
    data. If you are hoping to just estimate kinetic parmaters for your one set of samples,
    then bakRFnFit(bakRFnData) will work. No other model implementaiton
    (e.g., with HybridFit set to TRUE) will work with this kind
    of data.")
  }
  
  validate_bakRFnData(new_bakRFnData(fns, metadf))
  
  
}
