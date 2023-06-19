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
#' @return returns data frame that can be passed to pheatmap::pheatmap()
#' @examples
#' \donttest{
#' # Simulate small dataset
#' sim <- Simulate_bakRData(100, nreps = 2)
#'
#' # Analyze data with bakRFit
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Number of features that made it past filtering
#' NF <- nrow(Fit$Fast_Fit$Effects_df)
#'
#' # Simulate mock differential expression data frame
#' DE_df <- data.frame(XF = as.character(1:NF),
#'                        log2FoldChange = stats::rnorm(NF, 0, 2))
#'
#' DE_df$stat <- DE_df$log2FoldChange/0.5
#'
#' DE_df$padj <- 2*stats::dnorm(-abs(DE_df$stat))
#'
#' # make heatmap matrix
#' Heatmap <- NSSHeat(bakRFit = Fit,
#'                DE_df = DE_df,
#'                bakRModel = "MLE")
#'
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


  DE_reso <-  DE_df[DE_df$padj < DE_cutoff & !is.na(DE_df$padj),]

  DE_XF <- DE_reso$XF


  NSS_eff_DE <- NSS_eff %>%
    dplyr::mutate(Sig_bakR = ifelse(padj < bakR_cutoff, ifelse(effect < 0, "Down", "Up"), "NS"),
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

  abs_max_DE <- max(c(abs(c(max_DE, min_DE))))
  abs_max_bakR <- max(c(abs(c(max_bakR, min_bakR))))
  abs_max_mech <- max(c(abs(c(max_mech, min_mech))))

  # Changed to leave mechanism score completely unperturbed
  heatmap_df <- heatmap_df %>%
    #dplyr::mutate(Mech_score = ifelse(Mech_score < 0, Mech_score*(lid/-min_mech), Mech_score*(lid/max_mech))) %>%
    dplyr::mutate(DE_score = DE_score*(abs_max_mech/abs_max_DE)) %>%
    dplyr::mutate(bakR_score =  bakR_score*(abs_max_mech/abs_max_bakR))

  heatmap_df <- as.data.frame(heatmap_df)
  row.names(heatmap_df) <- heatmap_df$XF
  heatmap_df <- heatmap_df[,c("bakR_score", "DE_score", "Mech_score")]

  heatmap_df <- heatmap_df[order(heatmap_df$Mech_score),]

  return(heatmap_df)

}


#' Construct heatmap for non-steady state (NSS) analysis with improved mechanism score
#'
#' This uses the output of bakR and a differential expression analysis software to construct
#' a dataframe that can be passed to pheatmap::pheatmap(). This heatmap will display the
#' result of a steady-state quasi-independent analysis of NR-seq data.
#'
#' Unlike NSSHeat, NSSHeat2 uses a mechanism scoring function that is not discontinuous
#' as the "degradation driven" vs. "synthesis driven" boundary. Instead, the score
#' approaches 0 as the function approaches the boundary from either side.
#' 
#' In addition, NSSHeat2 now defines a null model for the purpose of p value calculation using
#' the mechanism score. Under the null hypothesis, the mechanism score is the product of two
#' normal distributions with unit variance, one which has a non-zero mean. Simulation is used
#' to estimate the integral of this distribution, and the number of draws (which determines the
#' precision of the p value estimate) is determined by the \code{sims} parameter.
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
#' @param sims Number of simulation draws from null distribution for mechanism p value calculation
#' @importFrom magrittr %>%
#' @return returns list of data frames: heatmap_df and NSS_stats. 
#' The heatmap_dfdata frame can be passed to pheatmap::pheatmap().
#' The NSS_stats data frame contains all of the information passed to NSS_stats as well
#' as the raw mechanism scores. It also has p values calculated from the mechanism z scores.
#' @examples
#' \donttest{
#' # Simulate small dataset
#' sim <- Simulate_bakRData(100, nreps = 2)
#'
#' # Analyze data with bakRFit
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Number of features that made it past filtering
#' NF <- nrow(Fit$Fast_Fit$Effects_df)
#'
#' # Simulate mock differential expression data frame
#' DE_df <- data.frame(XF = as.character(1:NF),
#'                        log2FoldChange = stats::rnorm(NF, 0, 2))
#'
#' DE_df$stat <- DE_df$log2FoldChange/0.5
#'
#' DE_df$padj <- 2*stats::dnorm(-abs(DE_df$stat))
#'
#' # make heatmap matrix
#' Heatmap <- NSSHeat2(bakRFit = Fit,
#'                DE_df = DE_df,
#'                bakRModel = "MLE")
#'
#' }
#' @export
NSSHeat2 <- function(bakRFit,
                     DE_df,
                     bakRModel = c("MLE", "Hybrid", "MCMC"),
                     DE_cutoff = 0.05,
                     bakR_cutoff = 0.05,
                     Exp_ID = 2,
                     sims = 10000000){

  bakRModel <- match.arg(bakRModel)


  # Bind variables locally to resolve devtools::check() Notes
  padj <- pval <- log2FoldChange <- effect <- se <- Sig <- Sig_bakR <- score_bakR <- NULL
  mech <- stat <- DE_score <- Mech_score <- bakR_score <- NULL
  mech_stat <- NULL

  ### Checks
  if(sum(c("XF", "log2FoldChange", "stat", "pval", "padj") %in% colnames(DE_df)) < 5){
    stop("You are missing necessary columns in DE_df. Columns named XF, log2FoldChange,
         stat, pval, and padj must be included. See the Further-Analyses vignette for details")
  }

  if(sum(colnames(DE_df) %in% c("XF", "log2FoldChange", "stat", "pval", "padj")) > 5){
    stop("Looks like you have repeat columns of the same name in DE_df. Make sure column
         names are unique and that each column is correctly labeled and try again.")
  }

  if(!inherits(bakRFit, "bakRFit")){
    stop("Input for bakRFit is not an object of class bakRFit.")
  }

  if(!is.numeric(sims)){
    stop("sims must be numeric!")
  }else if(sims <= 100){
    stop("sims must be strictly greater than 100.")
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


  DE_reso <-  DE_df[DE_df$padj < DE_cutoff & !is.na(DE_df$padj),]

  DE_XF <- DE_reso$XF


  NSS_eff_DE <- NSS_eff %>%
    dplyr::mutate(score_bakR = effect/se)


  NSS_eff_DE <- NSS_eff_DE[, c("XF", "pval", "padj", "score_bakR")]

  colnames(NSS_eff_DE) <- c( "XF", "bakR_pval", "bakR_padj", "score_bakR")

  XF_both <- intersect(NSS_eff_DE$XF, DE_XF)

  NSS_eff_DE <- NSS_eff_DE[NSS_eff_DE$XF %in% XF_both,]
  DE_reso <- DE_reso[DE_reso$XF %in% XF_both,]


  #browser()

  test <- dplyr::right_join(NSS_eff_DE, DE_reso, by = "XF")


  ## Calculate zf and zde
  if(sum(NSS_eff_DE$bakR_padj < bakR_cutoff) > 0){
    zfn <- min(abs(NSS_eff_DE$score_bakR[NSS_eff_DE$bakR_padj < bakR_cutoff]))

  }else{
    zfn <- max(abs(NSS_eff_DE$score_bakR))
  }


  test_stat <- test %>%
    dplyr::mutate(mech_stat = ifelse(stat > 0,
                                     (score_bakR + zfn)*(stat),
                                     (score_bakR - zfn)*(stat)))
  
  ## Null is product of two independent normal distributions, both with unit variance and one with
  ## non-zero mean.
  
  # Simulate from null for empirical p-value calc
  null_x <- stats::rnorm(sims, mean = zfn)
  null_y <- stats::rnorm(sims)
  null_xy <- null_x*null_y
  
  # Calculate p value and multiple-test adjust
  test_stat <- test_stat %>% dplyr::rowwise() %>%
    dplyr::mutate(mech_pval = 2*sum(null_xy > abs(mech_stat))/sims ) %>%
    dplyr::ungroup()
  
  rm(null_x)
  rm(null_y)
  rm(null_xy)
  
  test_stat$mech_padj <- stats::p.adjust(test_stat$mech_pval, method= "BH")
  
  ## Calculate meta analysis p value (p value that either expression or fraction new has changed)
  test_stat$meta_pval <- stats::pchisq(-2*(log(test_stat$bakR_pval) + log(test_stat$pval)), 
                                       df = 4,
                                       lower.tail = FALSE)
  
  test_stat$meta_padj <- stats::p.adjust(test_stat$meta_pval, method = "BH")
  
  
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

  abs_max_DE <- max(c(abs(c(max_DE, min_DE))))
  abs_max_bakR <- max(c(abs(c(max_bakR, min_bakR))))
  abs_max_mech <- max(c(abs(c(max_mech, min_mech))))

  # Changed to leave mechanism score completely unperturbed
  heatmap_df <- heatmap_df %>%
    #dplyr::mutate(Mech_score = ifelse(Mech_score < 0, Mech_score*(lid/-min_mech), Mech_score*(lid/max_mech))) %>%
    dplyr::mutate(DE_score = DE_score*(abs_max_mech/abs_max_DE)) %>%
    dplyr::mutate(bakR_score =  bakR_score*(abs_max_mech/abs_max_bakR))

  # Instead, Let me preserve the shape intelligently

  heatmap_df <- as.data.frame(heatmap_df)
  row.names(heatmap_df) <- heatmap_df$XF
  heatmap_df <- heatmap_df[,c("bakR_score", "DE_score", "Mech_score")]

  heatmap_df <- heatmap_df[order(heatmap_df$Mech_score),]

  # Compile output
  nss_list <- list(heatmap_df = heatmap_df,
                   NSS_stats = test_stat)
  
  return(nss_list)

}
