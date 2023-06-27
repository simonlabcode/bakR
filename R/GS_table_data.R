#' Example cB data frame
#'
#' Subset of a GRAND-SLAM main output table from anlaysis of a dataset published in Luo et al. 2020.
#' Data for 300 randomly selected genes is included to keep file size small.
#'
#' @docType data
#'
#' @usage data(GS_table)
#'
#' @format A dataframe with 300 rows and 63 variables; each row corresponds to GRAND-SLAM parameter
#' estimates for a single gene and 6 different samples (4 +s4U and 2 -s4U). Description of
#' all columns can be found on [GRAND-SLAM wiki](https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM)
#'
#' @keywords GS_table
#'
#' @references Luo et al. (2020) Biochemistry. 59(42), 4121-4142
#'
#'
#' @examples
#' data(GS_table)
#' data(metadf)
#' fns <- GSprocessing(GS_table)
#' bdfo <- bakRFnData(fns, metadf)
"GS_table"