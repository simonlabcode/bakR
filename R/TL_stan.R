#' Bayesian hierarchical mixture model with Stan for TL-seq data analysis
#'
#' @export
#' @param data_list list to pass to Stan of form given by cBtoStan
#' @param keep_fit if TRUE, Stan fit object is included in output; typically large file so default FALSE.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
TL_stan <- function(data_list, keep_fit = FALSE, ...) {
  fit <- rstan::sampling(stanmodels$Replicates, data = data_list, ...)

  # Extract alpha to get number of genes
  genes_summary <- rstan::summary(fit, pars = "alpha", probs = c(0.5))$summary

  genes_summary <- genes_summary[, c("50%","mean", "sd")]

  genes_df <- as.data.frame(genes_summary)
  ngs <- nrow(genes_df)

  rm(genes_df)

  # Extract kdeg to get number of conditions
  MT_summary <- rstan::summary(fit, pars = "kd", probs = c(0.5))$summary

  MT_summary <- MT_summary[, c("50%","mean", "sd")]

  MT_df <- as.data.frame(MT_summary)
  nconds <- (nrow(MT_df)/ngs) - 1


  #Extract frac_new to get number of replicates
  R_summary <- rstan::summary(fit, pars = "frac_new", probs = c(0.5))$summary

  R_summary <- R_summary[, c("50%","mean", "sd")]

  R_df <- as.data.frame(R_summary)
  reps <- nrow(R_df)/((nconds+1)*ngs)

  rm(R_df)

  nreps <- rep(reps, times=nconds)

  #Extract effect sizes
  eff_summary <- rstan::summary(fit, pars = c("eff"), probs = c(0.5))$summary

  #Pull out mean and standard deviation of parameter estimate
  eff_gauss <- eff_summary[, c("50%","mean", "sd")]


  #Convert to data frame:
  eff_gauss <- as.data.frame(eff_gauss)


  #Extract L2FC_kdeg
  L2FC_summary <- rstan::summary(fit, pars = "L2FC_kd", probs = c(0.5))$summary

  #Pull out mean and standard deviation of parameter estimate
  L2FC_summary <- L2FC_summary[, c("50%","mean", "sd")]


  #Convert to data frame:
  L2FC_df <- as.data.frame(L2FC_summary)

  # Effects dataframe setup
  F_ID <- rep(seq(from=1, to=ngs), each=nconds) # Feature number vector
  Exp_ID <- rep(seq(from=1, to=nconds), times=ngs) # Experimental condition vector
  Effect <- eff_gauss$mean # Effect size vector
  Se <- eff_gauss$sd # Effect standard deviation vector
  L2FC_kdeg <- L2FC_df$mean # L2FC(kdeg) vector
  L2FC_kdeg_sd <- L2FC_df$sd # L2FC(kdeg) sd vector

  # Perform statistical significance analysis
  fit_stat <- ashr::ash(Effect, Se, method = "fdr")
  lfsr <- fit_stat$result$lfsr
  lfdr <- fit_stat$result$lfdr

  # Kdeg dataframe setup
  F_ID_kd <- rep(seq(from=1, to=ngs), each=(nconds+1))
  Exp_ID_kd <- rep(seq(from=1, to=(nconds+1)), times=ngs)
  kdeg <- MT_df$mean # kdeg vector
  kdeg_sd <- MT_df$sd # kdeg sd vector

  Effects_df <- data.frame(F_ID, Exp_ID, L2FC_kdeg, L2FC_kdeg_sd, Effect, Se, lfsr, lfdr)
  Kdeg_df <- data.frame(F_ID_kd, Exp_ID_kd, kdeg, kdeg_sd)

  colnames(Effects_df) <- c("Feature ID", "Exp. ID", "L2FC(kdeg)", "L2FC(kdeg) sd", "Effect", "Se", "lfsr", "lfdr")
  colnames(Kdeg_df) <- c("Feature ID", "Exp. ID", "kdeg", "kdeg_sd")

  if(keep_fit == FALSE){
    out <- list(Effects_df, Kdeg_df)
    names(out) <- c("Effects_df", "Kdeg_df")
  }else{
    out <- list(Effects_df, Kdeg_df, fit)
    names(out) <- c("Effects_df", "Kdeg_df", "Stan_fit")
  }

  return(out)
}
