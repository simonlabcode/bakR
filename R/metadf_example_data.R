#' Example meatdf data frame
#'
#' metadf dataframe describing the data present in the cB file that
#' can be loaded with \code{data(cB_example)}. The contents are discussed
#' in great detail in the Getting_started vignette.
#'
#' @docType data
#'
#' @usage data(metadf_example)
#'
#' @format A dataframe; row names are samples in the corresponding
#' cB, in the order they appear in the cB. Columns are tl (s4U label time)
#' and Exp_ID (numerical ID of reference and experimental conditions)
#'
#' @keywords metadf
#'
#'
#'
#' @examples
#' data(cB_example)
#' data(metadf_example)
#' dyndata <- DynamicSeqData(cB_small, metadf_small)
