#' Prep GRAND-SLAM output for \code{bakRFnData}
#'
#' This function creates a fraction new estimate data frame that can be passed to
#' \code{bakRFnData}, using main .tsv file output by GRAND-SLAM.
#' @param GS Table of read counts and NTR (fraction new) estimate parameters output by GRAND-SLAM. Corresponds
#' to the *run_name*.tsv file included in GRAND-SLAM output
#' @param use_symbol Logical; if TRUE, then Symbol column rather than Gene column is used
#' as feature column (XF) in output data frame.
#' @return A data frame that can be passed as the \code{fns} parameter to \code{bakRFnData}
#' @importFrom magrittr %>%
#' @examples
#' # Load GRAND-SLAM table
#' data("GS_table")
#'
#'
#' # Create bakRData object
#' fns <- GSprocessing(GS_table)
#'
#' @export
GSprocessing <- function(GS, use_symbol = FALSE){
  
  Gene <- Readcount <- alpha <- NULL
  
  GS <- dplyr::as_tibble(GS)
  
  ### Make sure that the relevant columns exist
  
  # Readcount
  if(sum(grepl("Readcount", colnames(GS))) == 0){
    stop("Readcount columns are not present! GRAND-SLAM outputs many tables, did you pass in the correct one?")
  }
  
  # alpha
  if(sum(grepl("alpha", colnames(GS))) == 0){
    stop("alpha columns are not present! GRAND-SLAM outputs many tables, did you pass in the correct one?")
  }
  
  # beta
  if(sum(grepl("beta", colnames(GS))) == 0){
    stop("beta columns are not present! GRAND-SLAM outputs many tables, did you pass in the correct one?")
  }
  
  if(!use_symbol){
    if(sum(grepl("Gene", colnames(GS))) == 0){
      stop("Gene column is not present! GRAND-SLAM outputs many tables, did you pass in the correct one?")
    }
    
    if(sum(grepl("Gene", colnames(GS))) > 1){
      stop("There is more than one Gene column present! GRAND-SLAM outputs many tables, did you pass in the correct one?")
    }
    
  }else{
    if(sum(grepl("Symbol", colnames(GS))) == 0){
      stop("Symbol column is not present! GRAND-SLAM outputs many tables, did you pass in the correct one?")
    }
    
    if(sum(grepl("Symbol", colnames(GS))) > 1){
      stop("There is more than one Symbol column present! GRAND-SLAM outputs many tables, did you pass in the correct one?")
    }
    
  }
  
  ### Make sure alpha and beta columns are present in equal abundance
  if(sum(grepl("alpha", colnames(GS))) != sum(grepl("beta", colnames(GS)))){
    stop("There are not the same number of alpha columns as beta columns! GRAND-SLAM outputs many tables, did you pass in the correct one?")
  }
  
  
  if(!use_symbol){
    # Filter out useless columns
    GS <- GS[,grepl("Readcount", colnames(GS)) |
               grepl("alpha", colnames(GS)) |
               grepl("beta", colnames(GS)) |
               grepl("Gene", colnames(GS))]
  }else{
    # Filter out useless columns
    GS <- GS[,grepl("Readcount", colnames(GS)) |
               grepl("alpha", colnames(GS)) |
               grepl("beta", colnames(GS)) |
               grepl("Symbol", colnames(GS))]
  }

  
  # Standard deviation of a beta distribution
  sd_beta <- function(alpha, beta){
    var <- (alpha*beta)/(((alpha + beta)^2)*(alpha + beta + 1))
    return(sqrt(var))
  }
  
  # Pivot longer to get one sample per row
  GS <- GS %>% 
    tidyr::pivot_longer(
      cols = !Gene,
      names_to = c("sample", ".value"), 
      names_sep = " "
    ) %>%
    dplyr::filter(Readcount > 0) %>%
    dplyr::mutate(fn = ifelse(is.na(alpha), 0, alpha/(alpha + beta))) %>%
    dplyr::mutate(se = ifelse(is.na(alpha), 0, sd_beta(alpha, beta)))
  
  if(!use_symbol){
    # Keep relevant columns and rename as necessary
    GS <- GS[,c("sample", "Gene", "fn", "se", "Readcount")]
    
  }else{
    # Keep relevant columns and rename as necessary
    GS <- GS[,c("sample", "Symbol", "fn", "se", "Readcount")]
    
  }

  colnames(GS) <- c("sample", "XF", "fn", "se", "n")
  
  return(GS)

  
  
}