#' Example cB data frame
#'
#' Subset of a cB file from the DCP2 dataset published in Luo et al. 2020.
#' The original file is large (69 MB), so the example cB file has been
#' downsampled and contains only 10 genes (rather than 25012). The columns
#' are described in the Getting_Started vignette.
#'
#' @docType data
#'
#' @usage data(cB_small)
#'
#' @format A dataframe with 2614 rows and 5 variables; each row corresponds to a group of sequencing reads
#' \describe{
#' \item{sample}{Sample name}
#' \item{TC}{Number of T-to-C mutations}
#' \item{nT}{Number of Ts}
#' \item{XF}{Name of feature to which the group of reads map; usually a gene name}
#' \item{n}{Number of identical sequencing reads}
#' }
#'
#' @keywords cB_small
#'
#' @references Luo et al. (2020) Biochemistry. 59(42), 4121-4142
#'
#'
#' @examples
#' data(cB_small)
#' data(metadf)
#' dyndata <- bakRData(cB_small, metadf)
"cB_small"
