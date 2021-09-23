#' DynamicSeqData object constructor for internal use
#'
#' This function efficiently creates an object of class DynamicSeqData
#' without performing rigorous checks
#' @param cB Dataframe with columns corresponding to feature ID, number of Ts, number of mutations, sample ID, and number of identical observations
#' @param metadata Dataframe detailing s4U label time and experimental ID of each sample
new_DynamicSeqData <- function(cB, metadf){
  stopifnot(is.data.frame(cB))
  stopifnot(is.data.frame(metadf))
  structure(list(cB = cB, metadf = metadf), class = "DynamicSeqData")
}

#' DynamicSeq Data object validator
#'
#' This functions ensures that input for DynamicSeqData object construction is valid
#'
#' @param obj An object of class DynamicseqData
#' @export
validate_DynamicSeqData <- function(obj){

  vals <- unclass(obj)
  cB <- vals$cB
  metadf <- vals$metadf

  ## Check if rownames of metadf are same as unique(cB$sample)
  if("sample" %in% colnames(cB)){
    if(!(all(rownames(metadf) %in% unique(cB$sample)) & all(unique(cB$sample) %in% rownames(metadf) ))){
      stop(
        "Row names of metadf are not same as unique(cB$sample).
        Make sure the order in which samples appear in cB
        matches the rownames of metadf",
        call. = FALSE
      )
    }
  }else{
    if(!(all(rownames(metadf) %in% unique(cB[,2])) & all(unique(cB[,2]) %in% rownames(metadf)))){
      stop(
        "Row names of metadf are not same as unique(cB[,2]).
        Make sure 2nd column of cB contains sample names and
        that these match the rownames of metadf",
        call. = FALSE
      )
    }
  }


  ## Check if Numerical data in cB is all integer
  if(sum(purrr::map_dbl(unclass(cB), is.integer)) < 3){
    stop(
      "There are less than 3 columns of cB containing integer data.
      cB should have a column corresponding to numbers of mutations (TC),
      numbers of Ts in a read (nT), and number of identical observations (n),
      all of which should be integer data.",
      call. = FALSE
    )
  }#else if(sum(purrr::map_dbl(unclass(cB), is.integer)) > 3){
  #   stop(
  #     "There are more than 3 columns of cB containing integer data. Remove
  #     unnecessary columns (i.e., those that do not correspond to the number of
  #     mutations (TC), numbers of Ts in a read (nT), and number of identical
  #     observations (n).",
  #     call. = FALSE
  #   )
  # }


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
    stop(
      "No columns in metadf contain strictly integer data. The column containing numerical IDs
      for experimental conditions (Exp_ID) should be strictly integer data. If Exp_ID is correctly
      defined, try casting the column to an int with as.integer().",
      call. = FALSE
    )
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

  obj

}


#' DynamicSeq Data object helper function for users
#'
#' This function creates an object of class DynamicSeqData
#' @param cB Dataframe with columns corresponding to feature ID, number of Ts, number of mutations, sample ID, and number of identical observations
#' @param metadata Dataframe detailing s4U label time and experimental ID of each sample
#' @export
DynamicSeqData <- function(cB, metadf){

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
    If this does not reflect the content of those columns, properly rearrange cB columns and rerun DynamicSeqData().")
  }

  if(!all(meta_cols %in% colnames(metadf))){
    colnames(metadf)[1:2] <- meta_cols
    warning("Renamed first 2 columns of metadf to tl and Exp_ID.
    If this does not reflect the content of those columns, properly rearrange metadf columns and rerun DynamicSeqData().")
  }

  validate_DynamicSeqData(new_DynamicSeqData(cB, metadf))


}
