#' Creating PCA plots with logit(fn) estimates
#'
#' This function creates a 2-component PCA plot with logit(fn) estimates.
#'
#' @param obj Object of class FastFit or HMCFit outputted by respective analysis functions
#' @param ... Further arguments passed to or from other methods
#' @export
FnPCA <- function(obj, ...){
  ### Extract logit(fn)

  logit_fn_df <- obj$Fn_Estimates[,c("logit_fn", "Gene_ID", "Condition", "Replicate", "sample")]

  ### Create sample to [MT, R] lookup table

  sample_lookup <- logit_fn_df[,c("sample", "Condition", "Replicate")] %>% dplyr::distinct()


  ### Create logit_fn matrix

  logit_fn_mat <- matrix(0, ncol = nrow(sample_lookup), nrow = max(logit_fn_df$Gene_ID))

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

  Exp_ID <- as.factor(rep(1:max(sample_lookup$Condition), each = max(sample_lookup$Replicate)))


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
#' This function creates L2FC(kdeg) volcano plot.
#' Plots are colored according to statistical significance and sign of L2FC(kdeg).
#'
#' @param obj Object of class FastFit or HMCFit outputted by respective analysis functions
#' @param FDR False discovery rate to control at for significance assessment
#' @param Exps Vector of Experimental IDs to include in plot; must only contain elements within 2:(# of experimental IDs)
#' @param Exp_shape Logical indicating whether to use Expeirmental ID as factor determining point shape in volcano plot
plotVolcano <- function(obj, FDR = 0.05, Exps = NULL, Exp_shape = FALSE, ...){
  ## Extract L2FC(kdeg) and padj
  L2FC_df <- obj$Effects_df[,c("L2FC_kdegs", "padj", "Condition_effects")]

  ## Add significance ID
  L2FC_df <- L2FC_df %>% dplyr::mutate(conclusion = ifelse(padj < FDR, ifelse(L2FC_kdegs < 0, "Stabilized", "Destabilized"), "Not Sig."))

  if(is.null(Exps)){
    Exps <- 2:max(as.integer(L2FC_df$Condition_effects))
  }

  if(Exp_shape){
    ggplot2::ggplot(L2FC_df[L2FC_df$Condition_effects %in% Exps, ], ggplot2::aes(x = L2FC_kdegs,y = -log10(padj), color = conclusion,  shape = as.factor(Condition_effects))) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::theme_classic() +
      ggplot2::xlab("L2FC(kdeg)") +
      ggplot2::ylab("-log10(padj)")
  }else{
    ggplot2::ggplot(L2FC_df[L2FC_df$Condition_effects %in% Exps, ], ggplot2::aes(x = L2FC_kdegs,y = -log10(padj), color = conclusion )) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::theme_classic() +
      ggplot2::xlab("L2FC(kdeg)") +
      ggplot2::ylab("-log10(padj)")
  }


}

