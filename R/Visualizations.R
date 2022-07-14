#' Creating PCA plots with logit(fn) estimates
#'
#' This function creates a 2-component PCA plot using logit(fn) estimates.
#'
#' @param obj Object contained within output of \code{bakRFit}. So, either Fast_Fit (MLE implementation fit),
#' Stan_Fit (MCMC implementation fit), or Hybrid_Fit (Hybrid implementation fit)
#' @param log_kdeg Boolean; if TRUE, then log(kdeg) estimates used for PCA rather than logit(fn). Currently
#' only compatible with Fast_Fit
#' @export
FnPCA <- function(obj, log_kdeg = FALSE){
  ### Extract logit(fn)

  if(log_kdeg){
    if(class(obj) != "FastFit"){
      stop("log(kdeg) PCA is curently only compatible with Fast_Fit.")
    }

    logit_fn_df <- obj$Fn_Estimates[,c("log_kdeg", "Feature_ID", "Exp_ID", "Replicate", "sample")]
    colnames(logit_fn_df) <- c("logit_fn", "Feature_ID", "Exp_ID", "Replicate", "sample")
  }else{
    logit_fn_df <- obj$Fn_Estimates[,c("logit_fn", "Feature_ID", "Exp_ID", "Replicate", "sample")]

  }

  ### Create sample to [MT, R] lookup table

  sample_lookup <- logit_fn_df[,c("sample", "Exp_ID", "Replicate")] %>% dplyr::distinct()


  ### Create logit_fn matrix

  logit_fn_mat <- matrix(0, ncol = nrow(sample_lookup), nrow = max(logit_fn_df$Feature_ID))

  count <- 1
  for(i in sample_lookup$sample){
    logit_fn_mat[,count] <- logit_fn_df$logit_fn[logit_fn_df$sample == i]
    count <- count + 1
  }

  ## Cast to a dataframe with colnames = sample names

  logit_fn_mat <- as.data.frame(logit_fn_mat)
  colnames(logit_fn_mat) <- sample_lookup$sample

  ### Perform PCA

  fn_pcs <- stats::prcomp(t(logit_fn_mat), center = TRUE, scale. = FALSE)


  ### Extract loadings and PCs

  fn_eigenvect <- fn_pcs$x

  fn_PC1 <- fn_eigenvect[,c("PC1")]
  fn_PC2 <- fn_eigenvect[,c("PC2")]

  Exp_ID <- as.factor(rep(1:max(sample_lookup$Exp_ID), each = max(sample_lookup$Replicate)))


  fn_pca_df <- data.frame(PC1 = fn_PC1,
                          PC2 = fn_PC2,
                          Exp_ID = Exp_ID)

  fn_prop_var <- unclass(summary(fn_pcs))$importance[c("Proportion of Variance"),]

  ### Plot

  (g_pca <- ggplot2::ggplot(fn_pca_df, ggplot2::aes(x = PC1, y = PC2, color = Exp_ID)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::xlab(paste0("PC1 (", fn_prop_var[1]*100, "% of Var.)")) +
      ggplot2::ylab(paste0("PC2 (", fn_prop_var[2]*100, "% of Var.)")))

  return(g_pca)

}

#' Creating L2FC(kdeg) volcano plot from fit objects
#'
#' This function creates a L2FC(kdeg) volcano plot.
#' Plots are colored according to statistical significance and sign of L2FC(kdeg).
#'
#' @param obj Object contained within output of \code{bakRFit}. So, either Fast_Fit (MLE implementation fit),
#' Stan_Fit (MCMC implementation fit), or Hybrid_Fit (Hybrid implementation fit)
#' @param FDR False discovery rate to control at for significance assessment
#' @param Exps Vector of Experimental IDs to include in plot; must only contain elements within 2:(# of experimental IDs)
#' @param Exp_shape Logical indicating whether to use Expeirmental ID as factor determining point shape in volcano plot
#' @export
plotVolcano <- function(obj, FDR = 0.05, Exps = NULL, Exp_shape = FALSE){

  # Pretty plotting theme
  theme_mds <-    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                 panel.grid.minor = ggplot2::element_blank(),

                                 panel.background = ggplot2::element_blank(),

                                 axis.line.x = ggplot2::element_line(colour = "black"),

                                 axis.line.y = ggplot2::element_line(colour = "black"),

                                 axis.ticks = ggplot2::element_line(colour = "black"),

                                 axis.text = ggplot2::element_text(color="black", size = 12),

                                 axis.title = ggplot2::element_text(color="black", size = 18),

                                 strip.background = ggplot2::element_blank())

  ## Extract L2FC(kdeg) and padj
  L2FC_df <- obj$Effects_df[,c("L2FC_kdeg", "padj", "Exp_ID")]

  ## Add significance ID
  L2FC_df <- L2FC_df %>% dplyr::mutate(conclusion = ifelse(padj < FDR, ifelse(L2FC_kdeg < 0, "Stabilized", "Destabilized"), "Not Sig."))

  if(is.null(Exps)){
    Exps <- 2:max(as.integer(L2FC_df$Exp_ID))
  }

  if(Exp_shape){ # Plot different experimental conditions together using different shapes
    ggplot2::ggplot(L2FC_df[L2FC_df$Exp_ID %in% Exps, ], ggplot2::aes(x = L2FC_kdeg,y = -log10(padj), color = conclusion,  shape = as.factor(Exp_ID))) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::theme_classic() +
      ggplot2::ylab(bquote(log[10](p[adj]))) +
      ggplot2::xlab(bquote(L2FC(k[deg]))) +
      ggplot2::scale_color_manual(values = c("#FFC20A", "gray","#0C7BDC")) +
      theme_mds
  }else{ # Plot a subset of experimental conditions
    ggplot2::ggplot(L2FC_df[L2FC_df$Exp_ID %in% Exps, ], ggplot2::aes(x = L2FC_kdeg,y = -log10(padj), color = conclusion )) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::theme_classic() +
      ggplot2::ylab(bquote(log[10](p[adj]))) +
      ggplot2::xlab(bquote(L2FC(k[deg]))) +
      ggplot2::scale_color_manual(values = c("#FFC20A", "gray","#0C7BDC")) +
      theme_mds


  }


}


#' Creating L2FC(kdeg) MA plot from fit objects
#'
#' This function outputs a L2FC(kdeg) MA plot. Plots are colored according to statistical
#' significance and the sign of L2FC(kdeg)
#'
#' @param obj Object of class bakRFit outputted by \code{bakRFit} function
#' @param Model String identifying implementation for which you want to generate an MA plot
#' @param FDR False discovery rate to control at for significance assessment
#' @param Exps Vector of Experimental IDs to include in plot; must only contain elements within 2:(# of experimental IDs)
#' @param Exp_shape Logical indicating whether to use Expeirmental ID as factor determining point shape in volcano plot
#' @export
#'
plotMA <- function(obj, Model = c("MLE", "Hybrid", "MCMC"), FDR = 0.05, Exps = NULL, Exp_shape = FALSE){

  Model <- match.arg(Model)

  # Pretty plotting theme
  theme_mds <-    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                 panel.grid.minor = ggplot2::element_blank(),

                                 panel.background = ggplot2::element_blank(),

                                 axis.line.x = ggplot2::element_line(colour = "black"),

                                 axis.line.y = ggplot2::element_line(colour = "black"),

                                 axis.ticks = ggplot2::element_line(colour = "black"),

                                 axis.text = ggplot2::element_text(color="black", size = 12),

                                 axis.title = ggplot2::element_text(color="black", size = 18),

                                 strip.background = ggplot2::element_blank())

  ## Extract L2FC(kdeg) and padj
  if(Model == "MLE"){
    L2FC_df <- obj$Fast_Fit$Effects_df[,c("L2FC_kdeg", "padj", "Exp_ID", "Feature_ID")]

  }else if(Model == "Hybrid"){
    L2FC_df <- obj$Hybrid_Fit$Effects_df[,c("L2FC_kdeg", "padj", "Exp_ID", "Feature_ID")]

  }else if(Model == "MCMC"){
    L2FC_df <- obj$Stan_Fit$Effects_df[,c("L2FC_kdeg", "padj", "Exp_ID", "Feature_ID")]

  }

  Avg_reads_natural <- obj$Data_lists$Stan_data$Avg_Reads_natural

  ## Extract L2FC(kdeg) and padj

  ## Add significance ID
  L2FC_df <- L2FC_df %>% dplyr::mutate(conclusion = ifelse(padj < FDR, ifelse(L2FC_kdeg < 0, "Stabilized", "Destabilized"), "Not Sig."))

  if(is.null(Exps)){
    Exps <- 2:max(as.integer(L2FC_df$Exp_ID))
  }

  ## Calc avg reads per pairs of conditions
  Reads <- matrix(0, ncol = max(L2FC_df$Exp_ID), nrow = max(L2FC_df$Feature_ID))
  for(i in seq_along(unique(L2FC_df$Exp_ID))){
    Reads[,i] <- log10(rowMeans(Avg_reads_natural[,c(1,i+1)]))
  }

  ## Add read cnt info to L2FC_df

  L2FC_df <- L2FC_df %>% dplyr::group_by(Feature_ID, Exp_ID)  %>% dplyr::mutate(Read_ct = Reads[Feature_ID, Exp_ID-1]) %>% dplyr::ungroup()

  if(Exp_shape){
    ggplot2::ggplot(L2FC_df[L2FC_df$Exp_ID %in% Exps, ], ggplot2::aes(x = Read_ct,y = L2FC_kdeg, color = conclusion,  shape = as.factor(Exp_ID))) +
      ggplot2::geom_point(size = 1) +
      ggplot2::theme_classic() +
      ggplot2::xlab(expression(log[10](Avg.~read~count))) +
      ggplot2::ylab(bquote(L2FC(k[deg]))) +
      ggplot2::scale_color_manual(values = c("#FFC20A", "gray","#0C7BDC")) +
      theme_mds
  }else{
    ggplot2::ggplot(L2FC_df[L2FC_df$Exp_ID %in% Exps, ], ggplot2::aes(x = Read_ct,y = L2FC_kdeg, color = conclusion )) +
      ggplot2::geom_point(size = 1) +
      ggplot2::theme_classic() +
      ggplot2::xlab(expression(log[10](Avg.~read~count))) +
      ggplot2::ylab(bquote(L2FC(k[deg]))) +
      ggplot2::scale_color_manual(values = c("#FFC20A", "gray","#0C7BDC")) +
      theme_mds
  }


}

#' Creating a L2FC(kdeg) matrix that can be passed to heatmap functions
#'
#' \code{Heatmap_kdeg} creates a matrix where each column represents a pair of samples (reference and experimental) and each
#' row represents a feature. The entry in the ith row and jth column is the L2FC(kdeg) for feature i when comparing sample with
#' experimental ID j+1 to the reference sample
#'
#' @param obj Object of class FastFit or HMCFit outputted by respective analysis functions
#' @param zscore Logical; if TRUE, then each matrix entry is log-odds fold change in the fraction new (a.k.a the effect size) divided by
#' the uncertainty in the effect size
#' @param filter_sig Logical; if TRUE, then only features which have a statistically significant L2FC(kdeg) in at least one comparison
#' are kept
#' @param FDR Numeric; False discovery to control at if filter_sig is TRUE.
#' @export
Heatmap_kdeg <- function(obj, zscore = FALSE, filter_sig = FALSE, FDR = 0.05){
  ## Extract L2FC(kdeg) and padj
  if(zscore){


    L2FC_df <- obj$Effects_df[,c("effect", "se", "Exp_ID", "Feature_ID", "padj")]

    if(class(obj) == "FastFit"){
      df <- obj$Hyper_Parameters[1] + 2*max(obj$Fn_Estimates$Replicate) - 2

      L2FC_df <- L2FC_df %>% dplyr::mutate(zscore = ifelse(effect < 0, 1, -1)*stats::qt(padj, df = df) )
    }


    if(filter_sig){
      L2FC_df <- L2FC_df %>% dplyr::mutate(Sig = ifelse(padj < FDR, 1, 0))

      sig_df <- L2FC_df %>% dplyr::group_by(Feature_ID) %>% dplyr::summarise(Sig = sum(Sig), .groups = "keep") %>% dplyr::ungroup()

      fnum_keep <- sig_df$Feature_ID[sig_df$Sig > 0]

      L2FC_df <- L2FC_df[L2FC_df$Feature_ID %in% fnum_keep,]

      if(nrow(L2FC_df) == 0){
        stop("Filtering for significance removed all features; try increasing FDR or not filtering for significance")
      }

    }

    L2FC_df <- L2FC_df[order(L2FC_df$Feature_ID, L2FC_df$Exp_ID),]

    L2FC_mat <- matrix(0, nrow = length(unique(L2FC_df$Feature_ID)), ncol = max(L2FC_df$Exp_ID)-1)

    for(i in seq_along(unique(L2FC_df$Exp_ID))){
      if(class(obj) == "FastFit"){
        L2FC_mat[,i] <- L2FC_df$zscore[L2FC_df$Exp_ID == (i+1)]

      }else{
        L2FC_mat[,i] <- L2FC_df$effect[L2FC_df$Exp_ID == (i+1)]/L2FC_df$se[L2FC_df$Exp_ID == (i+1)]

      }
    }

  }else{
    L2FC_df <- obj$Effects_df[,c("L2FC_kdeg", "Exp_ID", "Feature_ID", "padj")]

    if(filter_sig){
      L2FC_df <- L2FC_df %>% dplyr::mutate(Sig = ifelse(padj < FDR, 1, 0))

      sig_df <- L2FC_df %>% dplyr::group_by(Feature_ID) %>% dplyr::summarise(Sig = sum(Sig), .groups = "keep") %>% dplyr::ungroup()

      fnum_keep <- sig_df$Feature_ID[sig_df$Sig > 0]

      L2FC_df <- L2FC_df[L2FC_df$Feature_ID %in% fnum_keep,]

      if(nrow(L2FC_df) == 0){
        stop("Filtering for significance removed all features; try increasing FDR or not filtering for significance")
      }

    }

    L2FC_df <- L2FC_df[order(L2FC_df$Feature_ID, L2FC_df$Exp_ID),]

    L2FC_mat <- matrix(0, nrow = length(unique(L2FC_df$Feature_ID)), ncol = max(L2FC_df$Exp_ID)-1)

    for(i in seq_along(unique(L2FC_df$Exp_ID))){
      L2FC_mat[,i] <- L2FC_df$L2FC_kdeg[L2FC_df$Exp_ID == (i+1)]
    }
  }

  return(L2FC_mat)





}
