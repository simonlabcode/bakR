#' Bayesian hierarchical mixture model with Stan for TL-seq data analysis
#'
#' \code{TL_stan} analyzes nucleotide recoding sequencing data with a fully
#' Bayesian hierarchical model implemented in the PPL Stan.
#'
#' Details of the model can be found in Vock et al. 2021. In short, mutations
#' are modeled as coming from a Poisson distribution with rate parameter
#' adjusted by the empirical U-content of each feature analyzed. Features
#' represent whatever the user defined them to be when constructing the
#' DynamicSeq data object. Typical feature categories are genes, exons, etc.
#' Multiple test adjusted significance values are also calculated for all
#' parameters quantifying changes in degradation and synthesis rate constants.
#' Local false-discovery rates and false sign rates (lfdr and lfsr,
#' respectively) are provided for each feature and each experimental condition.
#' Kinetic parameter changes are calcualted with respect to the user defined
#' reference sample, but additional comparisons can be made downstream (see
#' vignette for details).
#'
#' @export
#' @param data_list list to pass to Stan of form given by cBtoStan
#' @param Hybrid_Fit if TRUE, Hybrid Stan model that takes as data output of fast_analysis is run.
#' @param keep_fit if TRUE, Stan fit object is included in output; typically large file so default FALSE.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
TL_stan <- function(data_list, Hybrid_Fit = FALSE, Pooled = TRUE, keep_fit = FALSE, ...) {

  if(Pooled){

    if(Hybrid_Fit){
      fit <- rstan::sampling(stanmodels$Hybrid_Pooled, data = data_list, ...)

    }else{
      fit <- rstan::sampling(stanmodels$Heterosked_Pooled, data = data_list, ...)

    }

  }else{
    if(Hybrid_Fit){
      fit <- rstan::sampling(stanmodels$Hybrid, data = data_list, ...)
    }else{
      fit <- rstan::sampling(stanmodels$Heterosked, data = data_list, ...)
    }
  }


  fn_summary <- rstan::summary(fit, pars = "mu_rep_logit_fn", probs = c(0.5))$summary

  fn_summary <- fn_summary[, c("50%","mean", "sd")]

  fn_summary <- as.data.frame(fn_summary)

  colnames(fn_summary) <- c("median", "logit_fn", "sd")


  # Get number of features from data
  ngs <- data_list$NF

  # Extract kdeg to get number of conditions (could get from nMT but need kdeg df anyway)
  MT_summary <- rstan::summary(fit, pars = "kd", probs = c(0.5))$summary

  MT_summary <- MT_summary[, c("50%","mean", "sd")]

  MT_summary <- as.data.frame(MT_summary)

  nconds <- data_list$nMT - 1


  #Extract frac_new to get number of replicates
  reps <- data_list$nrep

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

  rm(eff_gauss)

  # Fn dataframe setup
  F_ID_fn <- rep(1:ngs, each = (nconds+1)*nreps)
  Exp_ID_fn <- rep(rep(1:(nconds+1), each = nreps), times = ngs)
  R_ID_fn <- rep(1:nreps, times = (nconds+1)*ngs)
  logit_fn <- fn_summary$logit_fn
  fn_se <- fn_summary$sd

  rm(fn_summary)

  # Kdeg dataframe setup
  F_ID_kd <- rep(seq(from=1, to=ngs), each=(nconds+1))
  Exp_ID_kd <- rep(seq(from=1, to=(nconds+1)), times=ngs)
  kdeg <- MT_summary$mean # kdeg vector
  kdeg_sd <- MT_summary$sd # kdeg sd vector

  rm(MT_summary)

  Effects_df <- data.frame(F_ID, Exp_ID, L2FC_kdeg, L2FC_kdeg_sd, Effect, Se)
  Kdeg_df <- data.frame(F_ID_kd, Exp_ID_kd, kdeg, kdeg_sd)
  Fn_df <- data.frame(F_ID_fn, Exp_ID_fn, R_ID_fn, logit_fn, fn_se)

  colnames(Effects_df) <- c("Genes_effects", "Condition_effects", "L2FC_kdegs", "L2FC_kd_sds", "effects", "ses")
  colnames(Kdeg_df) <- c("Feature_ID", "Exp_ID", "kdeg", "kdeg_sd")
  colnames(Fn_df) <- c("Gene_ID", "Condition", "Replicate", "logit_fn", "logit_fn_se")

  Fn_df <- merge(Fn_df, data_list$sample_lookup, by.x = c("Condition", "Replicate"), by.y = c("mut", "reps"))

  Fn_df <- Fn_df[order(Fn_df$Gene_ID, Fn_df$Condition, Fn_df$Replicate),]

  if(keep_fit == FALSE){
    fit_summary <- as.data.frame(rstan::summary(fit)$summary)
    out <- list(Effects_df, Kdeg_df, Fn_df, fit_summary)
    names(out) <- c("Effects_df", "Kdeg_df", "Fn_Estimates","Fit_Summary")
  }else{
    out <- list(Effects_df, Kdeg_df, fit)
    names(out) <- c("Effects_df", "Kdeg_df", "Fn_Estimates","Stan_fit")
  }


  return(out)
}
