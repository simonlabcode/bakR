#' Creating PCA plots with logit(fn) estimates
#'
#' This function creates a 2-component PCA plot using logit(fn) estimates. 
#' \code{FnPCA} has been deprecated in favor of `FnPCA2`. The latter accepts a full
#' bakRFit as input and handles imbalanced replicates.
#'
#' @param obj Object contained within output of \code{bakRFit}. So, either Fast_Fit (MLE implementation fit),
#' Stan_Fit (MCMC implementation fit), or Hybrid_Fit (Hybrid implementation fit)
#' @param log_kdeg Boolean; if TRUE, then log(kdeg) estimates used for PCA rather than logit(fn). Currently
#' only compatible with Fast_Fit
#' @return A ggplot object.
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Simulate data for 500 genes and 2 replicates
#' sim <- Simulate_bakRData(500, nreps = 2)
#'
#' # Fit data with fast implementation
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Fn PCA
#' FnPCA(Fit$Fast_Fit)
#'
#' # log(kdeg) PCA
#' FnPCA(Fit$Fast_Fit, log_kdeg = TRUE)
#'
#' }
#' @export
FnPCA <- function(obj, log_kdeg = FALSE){

  .Deprecated("FnPCA2")  

  # Bind variables locally to resolve devtools::check() Notes
  PC1 <- PC2 <- NULL

  ### Extract logit(fn)

  if(log_kdeg){
    if(!inherits(obj, "FastFit")){
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

#' Creating PCA plots with logit(fn) estimates
#'
#' This function creates a 2-component PCA plot using logit(fn) or log(kdeg) estimates.
#'
#' @param obj bakRFit object
#' @param Model String identifying implementation for which you want to generate a PCA plot
#' @param log_kdeg Boolean; if TRUE, then log(kdeg) estimates used for PCA rather than logit(fn). Currently
#' only compatible with MLE implementation
#' @return A ggplot object.
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Simulate data for 500 genes and 2 replicates
#' sim <- Simulate_bakRData(500, nreps = 2)
#'
#' # Fit data with fast implementation
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Fn PCA
#' FnPCA2(Fit$Fast_Fit)
#'
#' # log(kdeg) PCA
#' FnPCA2(Fit$Fast_Fit, log_kdeg = TRUE)
#'
#' }
#' @export
FnPCA2 <- function(obj, Model = c("MLE", "Hybrid", "MCMC"), log_kdeg = FALSE){
  
  # Bind variables locally to resolve devtools::check() Notes
  PC1 <- PC2 <- NULL
  
  Model <- match.arg(Model)
  
  ### Extract logit(fn)
  
  if(log_kdeg){
    if(Model == "MCMC"){
      stop("log(kdeg) PCA plots are only compatible with MLE implementation")
    }else if(Model == "Hybrid"){
      stop("log(kdeg) PCA plots are only compatible with MLE implementation")
    }else{
      logit_fn_df <- obj$Fast_Fit$Fn_Estimates[,c("log_kdeg", "Feature_ID", "Exp_ID", "Replicate", "sample")]
      colnames(logit_fn_df) <- c("logit_fn", "Feature_ID", "Exp_ID", "Replicate", "sample")
    }
  }else{
    if(Model == "MCMC"){
      logit_fn_df <- obj$Stan_Fit$Fn_Estimates[,c("logit_fn", "Feature_ID", "Exp_ID", "Replicate", "sample")]
      
    }else if(Model == "Hybrid"){
      logit_fn_df <- obj$Hybrid_Fit$Fn_Estimates[,c("logit_fn", "Feature_ID", "Exp_ID", "Replicate", "sample")]
      
    }else{
      logit_fn_df <- obj$Fast_Fit$Fn_Estimates[,c("logit_fn", "Feature_ID", "Exp_ID", "Replicate", "sample")]
      
    }
  }
  
  if(log_kdeg){
    if(!inherits(obj, "FastFit")){
      stop("log(kdeg) PCA is curently only compatible with Fast_Fit.")
    }
    
    logit_fn_df <- obj$Fast_Fit$Fn_Estimates[,c("log_kdeg", "Feature_ID", "Exp_ID", "Replicate", "sample")]
    colnames(logit_fn_df) <- c("logit_fn", "Feature_ID", "Exp_ID", "Replicate", "sample")
  }else{
    logit_fn_df <- obj$Fast_Fit$Fn_Estimates[,c("logit_fn", "Feature_ID", "Exp_ID", "Replicate", "sample")]
    
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
  
  Exp_ID <- as.factor(rep(1:max(sample_lookup$Exp_ID), 
                          times = obj$Data_lists$Stan_data$nrep_vect))
  
  
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
#' @param Exp_shape Logical indicating whether to use Experimental ID as factor determining point shape in volcano plot
#' @return A ggplot object. Each point represents a transcript. The x-axis is the
#' log-2 fold change in the degradation rate constant and the y-axis is the log-10
#' transformed multiple test adjusted p value.
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Simulate data for 500 genes and 2 replicates
#' sim <- Simulate_bakRData(500, nreps = 2)
#'
#' # Fit data with fast implementation
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Volcano plot
#' plotVolcano(Fit$Fast_Fit)
#'
#' }
#' @export
plotVolcano <- function(obj, FDR = 0.05, Exps = NULL, Exp_shape = FALSE){

  # Bind variables locally to resolve devtools::check() Notes
  padj <- L2FC_kdeg <- conclusion <- Exp_ID <- NULL

  # Pretty plotting theme
  theme_mds <-    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                 panel.grid.minor = ggplot2::element_blank(),

                                 panel.background = ggplot2::element_blank(),

                                 axis.line.x = ggplot2::element_line(colour = "black"),

                                 axis.line.y = ggplot2::element_line(colour = "black"),

                                 axis.ticks = ggplot2::element_line(colour = "black"),

                                 axis.text = ggplot2::element_text(color="black", size = 10),

                                 axis.title = ggplot2::element_text(color="black", size = 14),

                                 strip.background = ggplot2::element_blank())

  ## Extract L2FC(kdeg) and padj
  L2FC_df <- obj$Effects_df[,c("L2FC_kdeg", "padj", "Exp_ID")]

  ## Add significance ID
  L2FC_df <- L2FC_df %>% dplyr::mutate(conclusion = ifelse(padj < FDR, ifelse(L2FC_kdeg < 0, "Stabilized", "Destabilized"), "Not Sig."))

  if(is.null(Exps)){
    Exps <- 2:max(as.integer(L2FC_df$Exp_ID))
  }

  ## Check which results are present and plan color scheme appropriately
  conclusions <- unique(L2FC_df$conclusion)

  if(length(conclusions) == 3){
    colors <- c("#FFC20A", "gray","#0C7BDC")
  }else if(length(conclusions) == 2){
    if(all(c("Stabilized", "Destabilized") %in% conclusions)){
      colors <- c("#FFC20A", "#0C7BDC")
    }else if(all(c("Stabilized", "Not Sig.") %in% conclusions)){
      colors <- c("gray","#0C7BDC")
    }else if(all(c("Destabilized", "Not Sig.") %in% conclusions)){
      colors <- c("#FFC20A", "gray")
    }
  }else{
    if(conclusions == "Stabilized"){
      colors <- c("#0C7BDC")
    }else if(conclusions == "Destabilized"){
      colors <- c("#FFC20A")
    }else{
      colors <- "gray"
    }

  }


  if(Exp_shape){ # Plot different experimental conditions together using different shapes
    ggplot2::ggplot(L2FC_df[L2FC_df$Exp_ID %in% Exps, ], ggplot2::aes(x = L2FC_kdeg,y = -log10(padj), color = conclusion,  shape = as.factor(Exp_ID))) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::theme_classic() +
      ggplot2::ylab(bquote(-log[10](p[adj]))) +
      ggplot2::xlab(bquote(L2FC(k[deg]))) +
      ggplot2::scale_color_manual(values = colors) +
      theme_mds
  }else{ # Plot a subset of experimental conditions
    ggplot2::ggplot(L2FC_df[L2FC_df$Exp_ID %in% Exps, ], ggplot2::aes(x = L2FC_kdeg,y = -log10(padj), color = conclusion )) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::theme_classic() +
      ggplot2::ylab(bquote(-log[10](p[adj]))) +
      ggplot2::xlab(bquote(L2FC(k[deg]))) +
      ggplot2::scale_color_manual(values = colors) +
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
#' @param Exp_shape Logical indicating whether to use Experimental ID as factor determining point shape in volcano plot
#' @return A ggplot object. Each point represents a transcript. The
#' x-axis is log-10 transformed replicate average read counts,
#' y-axis is the log-2 fold-change in the degradation rate constant.
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Simulate data for 500 genes and 2 replicates
#' sim <- Simulate_bakRData(500, nreps = 2)
#'
#' # Fit data with fast implementation
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Volcano plot
#' plotMA(Fit, Model = "MLE")
#'
#' }
#' @export
#'
plotMA <- function(obj, Model = c("MLE", "Hybrid", "MCMC"), FDR = 0.05, Exps = NULL, Exp_shape = FALSE){

  # Bind variables locally to resolve devtools::check() Notes
  padj <- L2FC_kdeg <- Feature_ID <- Exp_ID <- Read_ct <- conclusion <- NULL


  Model <- match.arg(Model)

  # Pretty plotting theme
  theme_mds <-    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                 panel.grid.minor = ggplot2::element_blank(),

                                 panel.background = ggplot2::element_blank(),

                                 axis.line.x = ggplot2::element_line(colour = "black"),

                                 axis.line.y = ggplot2::element_line(colour = "black"),

                                 axis.ticks = ggplot2::element_line(colour = "black"),

                                 axis.text = ggplot2::element_text(color="black", size = 10),

                                 axis.title = ggplot2::element_text(color="black", size = 14),

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

  ## Check which results are present and plan color scheme appropriately
  conclusions <- unique(L2FC_df$conclusion)

  if(length(conclusions) == 3){
    colors <- c("#FFC20A", "gray","#0C7BDC")
  }else if(length(conclusions) == 2){
    if(all(c("Stabilized", "Destabilized") %in% conclusions)){
      colors <- c("#FFC20A", "#0C7BDC")
    }else if(all(c("Stabilized", "Not Sig.") %in% conclusions)){
      colors <- c("gray","#0C7BDC")
    }else if(all(c("Destabilized", "Not Sig.") %in% conclusions)){
      colors <- c("#FFC20A", "gray")
    }
  }else{
    if(conclusions == "Stabilized"){
      colors <- c("#0C7BDC")
    }else if(conclusions == "Destabilized"){
      colors <- c("#FFC20A")

    }else{
      colors <- "gray"
    }

  }

  ## Add read cnt info to L2FC_df

  L2FC_df <- L2FC_df %>% dplyr::group_by(Feature_ID, Exp_ID)  %>% dplyr::mutate(Read_ct = Reads[Feature_ID, Exp_ID-1]) %>% dplyr::ungroup()

  if(Exp_shape){
    ggplot2::ggplot(L2FC_df[L2FC_df$Exp_ID %in% Exps, ], ggplot2::aes(x = Read_ct,y = L2FC_kdeg, color = conclusion,  shape = as.factor(Exp_ID))) +
      ggplot2::geom_point(size = 1) +
      ggplot2::theme_classic() +
      ggplot2::xlab(expression(log[10](Avg.~read~count))) +
      ggplot2::ylab(bquote(L2FC(k[deg]))) +
      ggplot2::scale_color_manual(values = colors) +
      theme_mds
  }else{
    ggplot2::ggplot(L2FC_df[L2FC_df$Exp_ID %in% Exps, ], ggplot2::aes(x = Read_ct,y = L2FC_kdeg, color = conclusion )) +
      ggplot2::geom_point(size = 1) +
      ggplot2::theme_classic() +
      ggplot2::xlab(expression(log[10](Avg.~read~count))) +
      ggplot2::ylab(bquote(L2FC(k[deg]))) +
      ggplot2::scale_color_manual(values = colors) +
      theme_mds
  }


}

#' Creating a L2FC(kdeg) matrix that can be passed to heatmap functions
#'
#' \code{Heatmap_kdeg} creates a matrix where each column represents a pair of samples (reference and experimental) and each
#' row represents a feature. The entry in the ith row and jth column is the L2FC(kdeg) for feature i when comparing sample with
#' experimental ID j+1 to the reference sample
#'
#' @param obj Object outputted by \code{bakRFit}
#' @param zscore Logical; if TRUE, then each matrix entry is log-odds fold change in the fraction new (a.k.a the effect size) divided by
#' the uncertainty in the effect size
#' @param filter_sig Logical; if TRUE, then only features which have a statistically significant L2FC(kdeg) in at least one comparison
#' are kept
#' @param FDR Numeric; False discovery to control at if filter_sig is TRUE.
#' @return A matrix. Rows represent transcripts which were differentially expressed
#' and columns represent (from left to right) differential kinetics z-score,
#' differential expression z-score, and a mechanism score where positive represents
#' synthesis driven and negative degradation driven changes in expression.
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Simulate data
#' sim <- Simulate_bakRData(1000)
#'
#' # Fit data with fast implementation
#' Fit <- bakRFit(sim$bakRData)
#'
#' # L2FC(kdeg) heatmap matrix
#' L2FC_kdeg_heat <- Heatmap_kdeg(Fit$Fast_Fit)
#'
#' }
#' @export
Heatmap_kdeg <- function(obj, zscore = FALSE, filter_sig = FALSE, FDR = 0.05){


  # Bind variables locally to resolve devtools::check() Notes
  effect <- padj <- Feature_ID <- Sig <- NULL

  ## Extract L2FC(kdeg) and padj
  if(zscore){


    L2FC_df <- obj$Effects_df[,c("effect", "se", "Exp_ID", "Feature_ID", "padj")]

    if(inherits(obj, "FastFit")){
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
      if(inherits(obj, "FastFit")){
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

#' Visualize dropout
#' 
#' \code{VisualizeDropout} fits dropout model with \code{QuantifyDropout},
#' reports the fit results, and then generates a ggplot object showing the
#' data used to infer the fit as well as the fitted nonlinear trend.
#' 
#' @param obj bakRFit or bakRFnFit object
#' @param keep_data Logical; if TRUE, will return data used to make plots along with
#' the plots themselves
#' @param no_message Logical; if TRUE, will not output message regarding estimated
#' rates of dropout in each sample
#' @return If keep_data is FALSE, then a list of `ggplot` objects are returned, one
#' for each +s4U sample. The plots show the relationship between a feature's fraction new
#' and the difference between its +s4U and -s4U read coverage. Nonlinear-least squares fit
#' is plotted on top of points as a blue line. If keep_data is TRUE, then the data used
#' to make the plots is returned in addition to the list of plots.
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Simulate data for 500 genes and 2 replicates with 40% dropout
#' sim <- Simulate_relative_bakRData(500, nreps = 2, p_do = 0.4)
#'
#' # Fit data with fast implementation
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Quantify dropout
#' DO_plots <- VisualizeDropout(Fit)
#'
#' }
#' @export
VisualizeDropout <- function(obj,
                             keep_data = FALSE,
                             no_message = FALSE){
  
  # Address "no visible binding" NOTEs
  type <- mut <- reps <- fn <- dropout <- pdo <- NULL
  
  if(inherits(obj, "bakRFit")){
    
    if(sum(obj$Data_lists$Fast_df$type == 0) == 0){
      stop("You do not have any -s4U control data!")
    }
    
    
    # Check that -s4U samples exist for all mut
    check <- obj$Data_lists$Fast_df %>%
      dplyr::filter(type == 0) %>%
      dplyr::select(mut) %>%
      dplyr::distinct()
    
    check <- check$mut
    check <- as.integer(check[order(check)])
  }else{
    
    if(is.null(obj$Data_lists$Ctl_data)){
      stop("You do not have any -s4U control data!")
    }
    
    # Check that -s4U samples exist for all mut
    check <- obj$Data_lists$Ctl_data %>%
      dplyr::select(Exp_ID) %>%
      dplyr::distinct()
    
    check <- check$Exp_ID
    check <- as.integer(check[order(check)])
  }
  
  
  if(!identical(check, 1:obj$Data_lists$Stan_data$nMT)){
    stop("You do not have at least one replicate of -s4U data for all experimental conditions!")
  }
  
  Dropout <- QuantifyDropout(obj, keep_data = TRUE, no_message = TRUE)
  
  Data_d <- Dropout$Input_data
  Fit_d <- Dropout$Dropout_df
  
  if(!no_message){
    message(paste0(c("Estimated rates of dropout are:", 
                     utils::capture.output(as.data.frame(Fit_d[,c("Exp_ID", "Replicate", "pdo")]))),
                   collapse = "\n"))
  }

  
  combined_df <- dplyr::inner_join(Fit_d, Data_d, by = dplyr::join_by(Exp_ID == mut, 
                                                                      Replicate == reps))
  
  ### Make plots
  exp_vect <- Fit_d$Exp_ID
  rep_vect <- Fit_d$Replicate
  
  dropout_plots <- vector(mode = "list", length = length(exp_vect))
  
  
  for(i in 1:nrow(Fit_d)){
    data_sub <- combined_df[combined_df$Exp_ID == exp_vect[i] & combined_df$Replicate == rep_vect[i],]
    
    ymax <- max(data_sub$dropout)
    npoints <- nrow(data_sub)
    k <- 0.75
    alpha <- exp(-(log10(npoints) - 1)*k)
    if(alpha > 1){
      alpha <- 1
    }
    
    dropout_plots[[i]] <- ggplot2::ggplot() + 
      ggplot2::geom_point(data = data_sub, ggplot2::aes(x = fn, y = dropout),alpha = alpha ) + 
      ggplot2::geom_line(data = data_sub, ggplot2::aes(x = fn, 
                                                       y = (-(scale*pdo)*fn)/((1-pdo) + fn*pdo) + scale),
                         color = "blue",
                         linewidth = 1) + 
      ggplot2::theme_classic() + 
      ggplot2::ylim(c(0, ymax + 0.5)) + 
      ggplot2::ggtitle('Blue line = nls fit')
    
  }
  
  names(dropout_plots) <- paste0("ExpID_", exp_vect, "_Rep_", rep_vect)
  

  if(keep_data){
    return(list(dropout_plots = dropout_plots,
                Input_data = Data_d,
                Model_fit = Fit_d))
    
  }else{
    return(dropout_plots)
    
  }
  
  
  
}
