#' Example cB data frame
#'
#' Subset of a cB file from the DCP2 dataset published in Luo et al. 2020.
#' The original file is large (69 MB), so the example cB file has been
#' downsampled and contains only 10 genes (rather than 25012). The columns
#' are described in the Getting_Started vignette.
#'
#' @docType data
#'
#' @usage data(cB_example)
#'
#' @format A dataframe
#'
#' @keywords cB
#'
#' @references Luo et al. (2020) Biochemistry. 59(42), 4121-4142
#'
#'
#' @examples
#' data(cB_example)
#' data(metadf_example)
#' dyndata <- DynamicSeqData(cB_small, metadf)
