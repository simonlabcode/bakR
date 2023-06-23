#' Prep GRAND-SLAM output for \code{bakRFnData}
#'
#' This function creates a fraction new estimate data frame that can be passed to
#' \code{bakRFnData}, using main .tsv file output by GRAND-SLAM.
#' @param GS Table of read counts and NTR (fraction new) estimate parameters output by GRAND-SLAM. Corresponds
#' to the *run name*.tsv file included in GRAND-SLAM output
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
  
  GS <- dplyr::as_tibble(GS)
  
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