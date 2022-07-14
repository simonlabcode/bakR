#' Example meatdf data frame
#'
#' metadf dataframe describing the data present in the cB file that
#' can be loaded with \code{data(cB_small)}. The contents are discussed
#' in great detail in the Getting_started vignette.
#'
#' @docType data
#'
#' @usage data(metadf)
#'
#' @format A dataframe with 6 rows and 2 variables: row names are samples in the corresponding
#' cB
#' \describe{
#' \item{tl}{time of s4U labeling, in hours}
#' \item{Exp_ID}{numerical ID of reference and experimental conditions; 1 is reference and 2 is the single experimental condition}
#'
#' }
#'
#' @keywords metadf
#'
#'
#'
#' @examples
#' data(cB_small)
#' data(metadf)
#' bakRdat <- bakRData(cB_small, metadf)
"metadf"
