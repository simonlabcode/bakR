#' Bayesian hierarchical mixture model with Stan for TL-seq data analysis
#'
#' @export
#' @param data_list list to pass to Stan of form given by cBtoStan
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
TL_stan <- function(data_list, ...) {
  out <- rstan::sampling(stanmodels$Replicates, data = data_list, ...)
  return(out)
}
