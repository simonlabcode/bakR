#' Construct heatmap for non-steady state (NSS) analysis
#'
#' This uses the output of bakR and a differential expression analysis software to construct
#' a dataframe that can be passed to pheatmap::pheatmap(). This heatmap will display the
#' result of a steady-state quasi-independent analysis of NR-seq data.
#'
#' @param bakRFit bakRFit object
#' @param DE_df dataframe of required format with differential expression analysis results. See
#' Further-Analyses vignette for details on what this dataframe should look like
#' @param bakRModel Model fit from which bakR implementation should be used? Options are MLE, Hybrid,
#' or MCMC
#' @param DE_cutoff padj cutoff for calling a gene differentially expressed
#' @param bakR_cutoff padj cutoff for calling a fraction new significantly changed
#' @param Exp_ID Exp_ID of experimental sample whose comparison to the reference sample you want to use.
#' Only one reference vs. experimental sample comparison can be used at a time
#' @param lid Maximum absolute value for standardized score present in output. This is for improving
#' aesthetics of any heatmap generated with the output.
#' @importFrom magrittr %>%
#' @return returns dataframe that can be passed to pheatmap::pheatmap()
#' @examples
#' \dontrun{
#' # See vignettes for example of this function in action
#' NSSHeatmap <- NSSHeat(bakRFit, DESeq2_df, bakRModel = "MLE")
#' }
#' @export
NSSHeat <- function(bakRFit,
                    DE_df,
                    bakRModel = c("MLE", "Hybrid", "MCMC"),
                    DE_cutoff = 0.05,
                    bakR_cutoff = 0.05,
                    Exp_ID = 2,
                    lid = 4){

  bakRModel <- match.arg(bakRModel)


  # Bind variables locally to resolve devtools::check() Notes
  padj <- log2FoldChange <- effect <- se <- Sig <- Sig_bakR <- score_bakR <- NULL
  mech <- stat <- DE_score <- Mech_score <- bakR_score <- NULL

  ### Checks
  if(sum(c("XF", "log2FoldChange", "stat", "padj") %in% colnames(DE_df)) < 4){
    stop("You are missing necessary columns in DE_df. Columns named XF, log2FoldChange,
         stat, and padj must be included. See the Further-Analyses vignette for details")
  }

  if(sum(colnames(DE_df) %in% c("XF", "log2FoldChange", "stat", "padj")) > 4){
    stop("Looks like you have repeat columns of the same name in DE_df. Make sure column
         names are unique and that each column is correctly labeled and try again.")
  }

  if(!inherits(bakRFit, "bakRFit")){
    stop("Input for bakRFit is not an object of class bakRFit.")
  }


  if(bakRModel == "MLE"){
    NSS_eff <- bakRFit$Fast_Fit$Effects_df


    if(is.null(NSS_eff)){
      stop("Your bakRFit is missing the Fast_Fit object, yet bakRModel == MLE.")
    }

    NSS_eff$XF <- as.character(NSS_eff$XF)


  }else if(bakRModel == "Hybrid"){
    NSS_eff <- bakRFit$Hybrid_Fit$Effects_df

    if(is.null(NSS_eff)){
      stop("Your bakRFit is missing the Hybrid_Fit object, yet bakRModel == Hybrid.
           Did you mean to use MLE or MCMC fit?")
    }

    NSS_eff$XF <- as.character(NSS_eff$XF)


  }else if(bakRModel == "MCMC"){
    NSS_eff <- bakRFit$Stan_Fit$Effects_df

    if(is.null(NSS_eff)){
      stop("Your bakRFit is missing the Stan_Fit object, yet bakRModel == MCMC.
           Did you mean to use MLE or Hybrid fit?")
    }

    NSS_eff$XF <- as.character(NSS_eff$XF)

  }


  if(nrow(NSS_eff) != nrow(DE_df)){
    stop("DE_df and the effect size dataframe from bakRFit (Effects_df) have results
         for a different number of features. Make sure the same set of genes were anlayzed
         with bakR and the differential expression analysis tool.")
  }


  DE_df <- DE_df %>%
    dplyr::mutate(Sig = ifelse(padj < DE_cutoff, ifelse(log2FoldChange < 0, "Down", "Up"), "NS"))


  DE_reso <-  DE_df[DE_df$padj < 0.05 & !is.na(DE_df$padj),]

  DE_XF <- DE_reso$XF


  NSS_eff_DE <- NSS_eff %>%
    dplyr::mutate(Sig_bakR = ifelse(padj < 0.05, ifelse(effect < 0, "Down", "Up"), "NS"),
           score_bakR = effect/se)


  NSS_eff_DE <- NSS_eff_DE[, c("pval", "padj", "score_bakR", "XF", "Sig_bakR")]

  colnames(NSS_eff_DE) <- c("bakR_pval", "bakR_padj", "score_bakR", "XF", "Sig_bakR")


  NSS_eff_DE <- NSS_eff_DE[NSS_eff_DE$XF %in% DE_XF,]

  #browser()

  test <- dplyr::right_join(NSS_eff_DE, DE_reso, by = "XF")


  ## Assess mechanism of differential expression
  test <- test %>%
    dplyr::mutate(mech = ifelse( Sig == "NS", "NS",  ifelse(is.na(Sig_bakR), "ksyn",
                                                     ifelse(Sig == "Down", ifelse(Sig_bakR == "Up", "kdeg", "ksyn"  ),
                                                            ifelse(Sig == "Up", ifelse(Sig_bakR == "Down", "kdeg", "ksyn"), "ksyn") ) ) ))



  ## Calculate mechanism score
  test_stat <- test %>%
    dplyr::mutate(score_bakR = ifelse(is.na(score_bakR), 0, score_bakR)) %>% dplyr::rowwise() %>%
    dplyr::mutate(mech_stat = ifelse(mech == "ksyn", ifelse( sign(score_bakR) == sign(stat), mean( c(abs(score_bakR), abs(stat)) ), abs(stat)/(abs(score_bakR) + 2) ) ,
                              ifelse(mech == "kdeg", -mean( c(abs(score_bakR), abs(stat)) ) , 0) ))


  heatmap_df <- dplyr::tibble(DE_score = test_stat$stat,
                       Mech_score = test_stat$mech_stat,
                       XF = test_stat$XF,
                       bakR_score = test_stat$score_bakR)


  sd_DE <- stats::sd(heatmap_df$DE_score)
  sd_bakR <- stats::sd(heatmap_df$bakR_score)
  sd_mech <- stats::sd(heatmap_df$Mech_score)


  ## Standardize columns
  heatmap_df <- heatmap_df %>%
    dplyr::mutate(DE_score = DE_score/sd_DE,
           Mech_score = Mech_score/sd_mech,
           bakR_score = bakR_score/sd_bakR)


  ## Scale all columns equally
  max_DE <- max(heatmap_df$DE_score)
  min_DE <- min(heatmap_df$DE_score)

  max_bakR <- max(heatmap_df$bakR_score)
  min_bakR <- min(heatmap_df$bakR_score)

  max_mech <- max(heatmap_df$Mech_score)
  min_mech <- min(heatmap_df$Mech_score)

  heatmap_df <- heatmap_df %>%
    dplyr::mutate(Mech_score = ifelse(Mech_score < 0, Mech_score*(lid/-min_mech), Mech_score*(lid/max_mech))) %>%
    dplyr::mutate(DE_score = ifelse(DE_score < 0, DE_score*(lid/-min_DE), DE_score*(lid/max_DE)) ) %>%
    dplyr::mutate(bakR_score = ifelse(bakR_score < 0, bakR_score*(lid/-min_bakR), bakR_score*(lid/max_bakR)) )


  heatmap_df <- as.data.frame(heatmap_df)
  row.names(heatmap_df) <- heatmap_df$XF
  heatmap_df <- heatmap_df[,c("bakR_score", "DE_score", "Mech_score")]

  heatmap_df <- heatmap_df[order(heatmap_df$Mech_score),]

  return(heatmap_df)

}
