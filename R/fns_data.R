#' Example fraction news (fns) data frame
#'
#' Subset of fraction new estimates for dataset published by Luo et al. (2020).
#' Fraction new estimates, uncertainties, and read counts are included for 300 
#' genes to keep the file size small.
#'
#' @docType data
#'
#' @usage data(fns)
#'
#' @format A dataframe with 1,800 rows and 5 variables. Input to \code{bakRFndata}
#' \describe{
#' \item{sample}{Sample name}
#' \item{XF}{Name of feature (e.g., ENSEMBL gene ID)}
#' \item{fn}{Estimate of fraction of reads from feature that were new}
#' \item{se}{Uncertainty in fraction new estimate (optional in bakRFnData)}
#' \item{n}{Number of sequencing reads}
#' }
#'
#' @keywords fns
#'
#' @references Luo et al. (2020) Biochemistry. 59(42), 4121-4142
#'
#'
#' @examples
#' data(fns)
#' data(metadf)
#' bakRFndataobj <- bakRFnData(fns, metadf)
"fns"