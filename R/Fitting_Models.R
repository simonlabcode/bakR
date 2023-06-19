#' Estimating kinetic parameters from nucleotide recoding RNA-seq data
#'
#' \code{bakRFit} analyzes nucleotide recoding RNA-seq data to estimate
#' kinetic parameters relating to RNA stability and changes in RNA
#' stability induced by experimental perturbations. Several statistical
#' models of varying efficiency and accuracy can be used to fit data.
#'
#' If \code{bakRFit} is run on a bakRData object, \code{cBprocess}
#' and then \code{fast_analysis} will always be called. The former will generate the processed
#' data that can be passed to the model fitting functions (\code{fast_analysis}
#' and \code{TL_Stan}). The call to \code{fast_analysis} will generate a list of dataframes
#' containing information regarding the \code{fast_analysis} fit. \code{fast_analysis} is always
#' called because its output is required for both \code{Hybrid_fit} and \code{TL_Stan}.
#'
#' If \code{bakRFit} is run on a bakRFit object, \code{cBprocess} will not be called again,
#' as the output of \code{cBprocess} will already be contained in the bakRFit object. Similarly,
#' \code{fast_analysis} will not be called again unless bakRFit is rerun on the bakRData object.
#' or if \code{FastRerun} is set to TRUE. If you want to generate model fits using different parameters for cBprocess,
#' you will have to rerun \code{bakRFit} on the bakRData object.
#' 
#' If \code{bakRFit} is run on a bakRFnData object, \code{fn_process} and then \code{avg_and_regularize}
#' will always be called. The former will generate the processed data that can be passed to the model
#' fitting functions (\code{avg_and_regularize} and \code{TL_Stan}, the latter only with HybridFit = TRUE).
#' 
#' If \code{bakRFit} is run on a bakRFnFit object. \code{fn_process} will not be called again, as
#' the output of \code{fn_process} will already be contained in the bakRFnFit object. Similary,
#' \code{avg_and_regularize} will not be called unless \code{bakRFit} is rerun on the bakRData object,
#' or if \code{FastRerun} is set to TRUE. If you want to generate model fits using different
#' parameters for \code{fn_process}, you will have to rerun \code{bakRFit} on the bakRData object.
#'
#' See the documentation for the individual fitting functions for details regarding how they analyze nucleotide
#' recoding data. What follows is a brief overview of how each works
#'
#' \code{fast_analysis} (referred to as the MLE implementation in the bakR paper)
#' either estimates mutation rates from + and (if available) - s4U samples or uses mutation rate estimates
#' provided by the user to perform maximum likelihood estimation (MLE) of the fraction of RNA that is labeled for each
#' replicate of nucleotide recoding data provided. Uncertainties for each replicate's estimate are approximated using
#' asymptotic results involving the Fisher Information and assuming known mutation rates. Replicate data
#' is pooled using an approximation to hierarchical modeling that relies on analytic solutions to simple Bayesian models.
#' Linear regression is used to estimate the relationship between read depths and replicate variability for uncertainty
#' estimation regularization, again performed using analytic solutions to Bayesian models.
#'
#' \code{TL_Stan} with Hybrid_Fit set to TRUE (referred to as the Hybrid implementation in the bakR paper)
#' takes as input estimates of the logit(fraction new) and uncertainty provided by \code{fast_analysis}.
#' It then uses 'Stan' on the backend to implement a hierarchical model that pools data across replicates and the dataset
#' to estimate effect sizes (L2FC(kdeg)) and uncertainties. Replicate variability information is pooled across each experimental
#' condition to regularize variance estimates using a hierarchical linear regression model.
#'
#' The default behavior of \code{TL_Stan} (referred to as the MCMC implementation in the bakR paper)
#' is to use 'Stan' on the back end to implement a U-content exposure adjusted Poisson mixture model
#' to estimate fraction news from the mutational data. Partial pooling of replicate variability estimates
#' is performed as with the Hybrid implementation.
#'
#'
#' @param obj bakRData object produced by \code{bakRData} or a bakRFit object produced by \code{bakRFit}
#' @param StanFit Logical; if TRUE, then the MCMC implementation is run. Will only be used if \code{obj}
#' is a \code{bakRFit} object
#' @param HybridFit Logical; if TRUE, then the Hybrid implementation is run. Will only be used if \code{obj}
#' is a \code{bakRFit} object
#' @param high_p Numeric; Any transcripts with a mutation rate (number of mutations / number of Ts in reads) higher than this in any -s4U control
#' samples (i.e., samples that were not treated with s4U) are filtered out
#' @param totcut Numeric; Any transcripts with less than this number of sequencing reads in any sample are filtered out
#' @param Ucut Numeric; All transcripts must have a fraction of reads with 2 or less Us less than this cutoff in all samples
#' @param AvgU Numeric; All transcripts must have an average number of Us greater than this cutoff in all samples
#' @param FastRerun Logical; only matters if a bakRFit object is passed to \code{bakRFit}. If TRUE, then the Stan-free
#' model implemented in \code{fast_analysis} is rerun on data, foregoing fitting of either of the 'Stan' models.
#' @param FOI Features of interest; character vector containing names of features to analyze
#' @param concat Logical; If TRUE, FOI is concatenated with output of reliableFeatures
#' @param StanRateEst Logical; if TRUE, a simple 'Stan' model is used to estimate mutation rates for fast_analysis; this may add a couple minutes
#' to the runtime of the analysis.
#' @param RateEst_size Numeric; if StanRateEst is TRUE, then data from RateEst_size genes are used for mutation rate estimation. This can be as low
#' as 1 and should be kept low to ensure maximum efficiency
#' @param low_reads Numeric; if StanRateEst is TRUE, then only features with more than low_reads reads in all samples will be used for mutation rate estimation
#' @param high_reads Numeric; if StanRateEst is TRUE, then only features with less than high_read reads in all samples will be used for mutation rate estimation.
#' A high read count cutoff is as important as a low read count cutoff in this case because you don't want the fraction labeled of chosen features to be
#' extreme (e.g., close to 0 or 1), and high read count features are likely low fraction new features.
#' @param chains Number of Markov chains to sample from. 1 should suffice since these are validated models. Running more chains is generally
#' preferable, but memory constraints can make this unfeasible.
#' @param NSS Logical; if TRUE, logit(fn)s are directly compared to avoid assuming steady-state
#' @param Chase Logical; Set to TRUE if analyzing a pulse-chase experiment. If TRUE, kdeg = -ln(fn)/tl where fn is the fraction of
#' reads that are s4U (more properly referred to as the fraction old in the context of a pulse-chase experiment).
#' @param BDA_model Logical; if TRUE, variance is regularized with scaled inverse chi-squared model. Otherwise a log-normal
#' model is used.
#' @param multi_pold Logical; if TRUE, pold is estimated for each sample rather than use a global pold estimate.
#' @param Long Logical; if TRUE, long read optimized fraction new estimation strategy is used.
#' @param kmeans Logical; if TRUE, kmeans clustering on read-specific mutation rates is used to estimate pnews and pold.
#' @param ztest Logical; if TRUE and the MLE implementation is being used, then a z-test will be used for p-value calculation
#' rather than the more conservative moderated t-test.
#' @param Fisher Logical; if TRUE, Fisher information is used to estimate logit(fn) uncertainty. Else, a less conservative binomial model is used, which
#' can be preferable in instances where the Fisher information strategy often drastically overestimates uncertainty
#' (i.e., low coverage or low pnew).
#' @param ... Arguments passed to either \code{fast_analysis} (if a bakRData object)
#' or \code{TL_Stan} and \code{Hybrid_fit} (if a bakRFit object)
#' @return bakRFit object with results from statistical modeling and data processing. Objects possibly included are:
#' \itemize{
#'  \item Fast_Fit; Always will be present. Output of \code{fast_analysis}
#'  \item Hybrid_Fit; Only present if HybridFit = TRUE. Output of \code{TL_stan}
#'  \item Stan_Fit; Only present if StanFit = TRUE. Output of \code{TL_stan}
#'  \item Data_lists; Always will be present. Output of \code{cBprocess} with Fast and Stan == TRUE
#' }
#' @importFrom magrittr %>%
#' @examples
#' \donttest{
#' # Simulate data for 1000 genes, 2 replicates, 2 conditions
#' simdata <- Simulate_bakRData(1000, nreps = 2)
#'
#' # You always must fit fast implementation before any others
#' Fit <- bakRFit(simdata$bakRData)
#'
#' }
#'
#' @export
bakRFit <- function(obj, StanFit = FALSE, HybridFit = FALSE,
                          high_p = 0.2,
                          totcut = 50,
                          Ucut = 0.25,
                          AvgU = 4,
                          FastRerun = FALSE,
                          FOI = c(),
                          concat = TRUE,
                          StanRateEst = FALSE,
                          RateEst_size = 30,
                          low_reads = 100,
                          high_reads = 500000,
                          chains = 1, NSS = FALSE,
                          Chase = FALSE, BDA_model = FALSE, multi_pold = FALSE,
                          Long = FALSE, kmeans = FALSE, ztest = FALSE, Fisher = TRUE,
                          ...){

  # Bind variables locally to resolve devtools::check() Notes
  old_fnum <- effect <- se <- pval <- NULL

  ## Check StanFit
  if(!is.logical(StanFit)){
    stop("StanFit must be logical (TRUE or FALSE)")
  }

  ## Check HybridFit
  if(!is.logical(HybridFit)){
    stop("HybridFit must be logical (TRUE or FALSE")
  }

  ## Check high_p
  if(!is.numeric(high_p)){
    stop("high_p must be numeric")
  }else if( (high_p < 0) | (high_p > 1) ){
    stop("high_p must be between 0 and 1")
  }else if (high_p < 0.01){
    warning("high_p is abnormally low (< 0.01); many features will by pure chance have a higher mutation rate than this in a -s4U control and thus get filtered out")
  }

  ## Check totcut
  if(!is.numeric(totcut)){
    stop("totcut must be numeric")
  }else if( totcut < 0 ){
    stop("totcut must be greater than 0")
  }else if(totcut > 5000){
    warning("totcut is abnormally high (> 5000); many features will not have this much coverage in every sample and thus get filtered out.")
  }

  ## Check Ucut
  if(!is.numeric(Ucut)){
    stop("Ucut must be numeric")
  }else if( Ucut < 0 ){
    stop("Ucut must be greater than or equal to 0")
  }else if(Ucut > 0.5 ){
    warning("Ucut is abnormally high; you are allowing > 50% of reads to have 2 or less Us.")
  }


  ## Check AvgU
  if(!is.numeric(AvgU)){
    stop("AvgU must be numeric")
  }else if(AvgU < 0){
    stop("AvgU must be greater than or equal to 0")
  }else if (AvgU > 50){
    warning("AvgU is abnormally high; you are requiring an average number of Us greater than 50")
  }else if(AvgU < 4){
    warning("AvgU is abnormally low; you are allowing an average of less than 4 Us per read, which may model convergence issues.")
  }


  ## Check concat
  if(!is.logical(concat)){
    stop("concat must be logical (TRUE or FALSE)")
  }

  ## Check StanRateEst
  if(!is.logical(StanRateEst)){
    stop("StanRateEst must be logical (TRUE or FALSE)")
  }

  ## Check Chase
  if(!is.logical(Chase)){
    stop("Chase must be logical (TRUE or FALSE)")
  }
  
  ## Check multi_pold
  if(!is.logical(multi_pold)){
    stop("multi_pold must be logical (TRUE or FALSE)")
  }

  ## Check RateEst_size
  if(!is.numeric(RateEst_size)){
    stop("RateEst_size must be numeric")
  }else{
    RateEst_size <- as.integer(RateEst_size)
    if(RateEst_size < 1){
      stop("RateEst_size must be greater than 0")
    }else if(RateEst_size > 100){
      warning("You have set RateEst_size to a number greater than 100. This will reduce efficiency without much benefit to estimate accuracy.")
    }else if(RateEst_size < 10){
      warning("You have set RateEst_size to a number less than 10. This may lead to inaccurate mutation rate estimates")
    }
  }

  ## Check low_reads
  if(!is.numeric(low_reads)){
    stop("low_reads must be numeric")
  }else if(low_reads < 0){
    stop("low_reads must be greater than 0")
  }

  if(low_reads < 50){
    warning("low_reads is less than 50. This may lead to inaccurate mutation rate estimation; tread lightly")
  }else if(low_reads > 1000){
    warning("low_reads is greater than 1000. There may not be enough features meeting this criterion.")
  }

  ## Check high_reads
  if(!is.numeric(high_reads)){
    stop("high_reads must be numeric")
  }else if(high_reads < 0){
    stop("high_reads must be greater than 0")
  }

  if(high_reads < 1000){
    warning("high_reads is less than 1000. This may lead to inaccurate mutation rate estimation; tread lightly")
  }else if(high_reads > 1000000){
    warning("high reads is greater than 1000000. This may mean that chosen features are almost completely unlabeled; tread lightly")
  }

  ## Check low_reads vs. high_reads
  if(high_reads < low_reads){
    stop("high_reads must be greater than low_reads")
  }else if(high_reads == low_reads){
    warning("low_reads and high_reads are equal. This is a very small window of allowable read depths, so no features might pass this condition.")
  }

  ## Check number of chains
  if(!is.numeric(chains)){
    stop("chains must be numerical")
  }else if(chains < 1){
    stop("chains must be greater than or equal to 1")
  }else if(!is.integer(chains)){
    chains <- as.integer(chains)
  }

  if(chains > 4){
    warning("You are running more than 4 chains. This may be more than necessary and may sacrifice computational performance.")
  }
  
  ### Check ztest
  if(!is.logical(ztest)){
    stop("ztest must be logical (TRUE or FALSE)")
  }
  
  ## Check Fisher
  if(!is.logical(Fisher)){
    stop("Fisher must be logical (TRUE or FALSE")
  }

  if(inherits(obj, "bakRData")){

    # Preprocess data
    data_list <- bakR::cBprocess(obj, high_p = high_p, totcut = totcut, Ucut = Ucut,
                                       AvgU = AvgU,
                                       Stan = TRUE,
                                       Fast = TRUE,
                                       FOI = FOI,
                                       concat = concat)

    # Use Stan to estimate mutation rates
    if(StanRateEst){

      ## find features to keep
      Cnt_Mat <- data_list$Count_Matrix

      Cnt_Mat <- (Cnt_Mat >= low_reads) & (Cnt_Mat <= high_reads)

      fnum_choose <- which(rowSums(Cnt_Mat) == ncol(Cnt_Mat))

      # Check that enough features made it past filtering
      if(length(fnum_choose) < RateEst_size){
        stop("Not enough features have read depths between low_reads and high_reads. Try increasing high_reads and/or decreasing low_reads.")
      }

      ## Create small cB for mutation rate estimation Stan model

      XF_choose <- sample(unique(data_list$Stan_data$sdf$XF[data_list$Stan_data$sdf$fnum %in% fnum_choose]), RateEst_size, replace = FALSE)

      cB_small <- obj$cB[obj$cB$XF %in% XF_choose,]

      bakRData2 <- new_bakRData(cB_small, obj$metadf)

      mutrate_list <- bakR::cBprocess(bakRData2, FOI = XF_choose, concat = FALSE)

      mutrate_list$Stan_data$nU <- sum(mutrate_list$Fast_df$nT*mutrate_list$Fast_df$n)/sum(mutrate_list$Fast_df$n)

      
      # Run MLE implementation
      fast_list <- bakR::fast_analysis(data_list$Fast_df, Stan_data = mutrate_list$Stan_data, StanRate = TRUE,
                                       BDA_model = BDA_model, Chase = Chase, multi_pold = multi_pold,
                                       NSS = NSS, Long = Long, kmeans = kmeans, ztest = ztest,
                                       Fisher = Fisher, ...)

    }else{
      # Run MLE implementation
      fast_list <- bakR::fast_analysis(data_list$Fast_df,
                                       BDA_model = BDA_model, Chase = Chase, multi_pold = multi_pold,
                                       NSS = NSS, Long = Long, kmeans = kmeans, ztest = ztest, 
                                       Fisher = Fisher, ...)
    }



    # Make final object to be passed to users
    Fit_lists <- list(Fast_Fit = fast_list,
                        Data_lists = data_list)

    class(Fit_lists) <- "bakRFit"
    return(Fit_lists)

  }else if(inherits(obj, "bakRFit")){

    if(FastRerun & (StanFit | HybridFit)){
      stop("Can only rerun MLE implementation or run Hybrid/MCMC implementation, not both. If you want to rerun the MLE implementation
           and then run either the MCMC or Hybrid implementations, use separate calls to bakRFit().")
    }

    if(FastRerun){
      if(StanRateEst){
        Cnt_Mat <- obj$Data_list$Count_Matrix

        Cnt_Mat <- (Cnt_Mat >= low_reads) & (Cnt_Mat <= high_reads)

        fnum_choose <- sample(which(rowSums(Cnt_Mat) == ncol(Cnt_Mat)), size = RateEst_size)

        if(length(fnum_choose) < RateEst_size){
          stop("Not enough features have read depths between low_reads and high_reads. Try increasing high_reads and/or decreasing low_reads.")
        }

        ## Construct mini-cB
        Stan_data <- obj$Data_list$Stan_data

        data_df <- data.frame(old_fnum = Stan_data$FE,
                              TP = Stan_data$TP,
                              MT = Stan_data$MT,
                              num_mut = Stan_data$num_mut,
                              R = Stan_data$R,
                              num_obs = Stan_data$num_obs,
                              U_cont = Stan_data$U_cont)


        ## Filter for features to be used in mutation rate estimation
        data_df <- data_df %>% dplyr::filter(old_fnum %in% fnum_choose)

        new_FE <- data.frame(old_fnum = unique(data_df$old_fnum[order(data_df$old_fnum)]),
                             FE = 1:RateEst_size)

        data_df <- dplyr::left_join(data_df, new_FE, by = "old_fnum")


        # List for mutation rate estimation Stan model
        mutrate_list <- list(FE = data_df$FE,
                          TP = data_df$TP,
                          MT = data_df$MT,
                          R = data_df$R,
                          num_mut = data_df$num_mut,
                          num_obs = data_df$num_obs,
                          U_cont = data_df$U_cont,
                          nMT = Stan_data$nMT,
                          nrep = Stan_data$nrep,
                          tl = Stan_data$tl,
                          NE = length(data_df$FE),
                          NF = max(data_df$FE),
                          sdf = Stan_data$sdf[Stan_data$sdf$fnum %in% fnum_choose,],
                          nU = sum(obj$Data_lists$Fast_df$nT*obj$Data_lists$Fast_df$n)/sum(obj$Data_lists$Fast_df$n))

        rm(Stan_data)

        # Run MLE implementation
        fast_list <- bakR::fast_analysis(obj$Data_lists$Fast_df, Stan_data = mutrate_list, StanRate = TRUE,
                                         NSS = NSS,
                                         BDA_model = BDA_model, multi_pold = multi_pold,
                                         Chase = Chase, Long = Long, kmeans = kmeans, ztest = ztest,
                                         Fisher = Fisher,
                                         ...)

      }else{
        # Run MLE implementation
        fast_list <- bakR::fast_analysis(obj$Data_lists$Fast_df,
                                         NSS = NSS,
                                         BDA_model = BDA_model, multi_pold = multi_pold,
                                         Chase = Chase, Long = Long, kmeans = kmeans, ztest = ztest,
                                         ...)

      }


      obj$Fast_Fit <- fast_list

      return(obj)
    }

    # Run MCMC implementation
    if(StanFit){

      obj$Data_lists$Stan_data$Chase <- as.integer(Chase)
      
      ## Can calculate # of Us
      obj$Data_lists$Stan_data$nU <- sum(obj$Data_lists$Fast_df$nT*obj$Data_lists$Fast_df$n)/sum(obj$Data_lists$Fast_df$n)

      Stan_list <- TL_stan(obj$Data_lists$Stan_data, NSS = NSS, chains = chains, ...)

      ## Calculate and adjust p-values
      Effects <- Stan_list$Effects_df

      Effects <- Effects %>% dplyr::mutate(pval = 2*stats::pnorm(-abs(effect/se))) %>%
        dplyr::mutate(padj = stats::p.adjust(pval, method = "BH"))

      Stan_list$Effects_df <- Effects

      rm(Effects)

      class(Stan_list) <- "HMCFit"

      obj$Stan_Fit <- Stan_list

    }

    # Run Hybrid implementation
    if(HybridFit){

      ## Get fn estimates from MLE implementation and curate data for Stan
      
      Rep_Fn <- obj$Fast_Fit$Fn_Estimates
      
      data_list <- list(
        NE = nrow(Rep_Fn),
        NF = max(Rep_Fn$Feature_ID),
        MT = Rep_Fn$Exp_ID,
        FE = Rep_Fn$Feature_ID,
        tl = obj$Data_lists$Stan_data$tl,
        logit_fn_rep = Rep_Fn$logit_fn,
        fn_se = Rep_Fn$logit_fn_se,
        Avg_Reads = obj$Data_lists$Stan_data$Avg_Reads,
        nMT = max(Rep_Fn$Exp_ID),
        R = Rep_Fn$Replicate,
        nrep = max(Rep_Fn$Replicate),
        sample_lookup = obj$Data_lists$Stan_data$sample_lookup,
        sdf = obj$Data_lists$Stan_data$sdf,
        mutrates = obj$Fast_Fit$Mut_rates,
        nrep_vect = obj$Data_lists$Stan_data$nrep_vect,
        Chase = as.integer(Chase)
      )
      
      
      rm(Rep_Fn)
      
      # Run Hybrid implementation
      Stan_list <- TL_stan(data_list, Hybrid_Fit = TRUE, NSS = NSS, chains = chains, ...)
      
      ## Calculate and adjust p-values
      
      Effects <- Stan_list$Effects_df
      
      Effects <- Effects %>% dplyr::mutate(pval = 2*stats::pnorm(-abs(effect/se))) %>%
        dplyr::mutate(padj = stats::p.adjust(pval, method = "BH"))
      
      Stan_list$Effects_df <- Effects
      
      rm(Effects)
      
      class(Stan_list) <- "HybridModelFit"
      
      obj$Hybrid_Fit <- Stan_list
      
    }

    if(!(StanFit | HybridFit)){
      stop("Either StanFit or HybridFit has to be true (or both).")
    }

    return(obj)

  }else if(inherits(obj, "bakRFnData")){
    
    data_list <- fn_process(obj, totcut = totcut, Chase = Chase, FOI = FOI, concat = concat)
    
    ## Necessary hacky preprocessing
    Mut_data_est <- data_list$Fn_est
    
    # There is a better way to do this...
    Mut_data_est <- Mut_data_est[,c("XF", "sample", "fn", "n",
                                    "Feature_ID", "Replicate", "Exp_ID", "tl",
                                    "logit_fn", "kdeg", "log_kdeg", "logit_fn_se", "log_kd_se")]
    
    colnames(Mut_data_est) <- c("XF", "sample", "fn", "nreads", "fnum", "reps",
                                "mut", "tl", "logit_fn_rep", "kd_rep_est", "log_kd_rep_est", 
                                "logit_fn_se", "log_kd_se")
    
    sample_lookup <- data_list$Stan_data$sample_lookup
    feature_lookup <- data_list$Stan_data$sdf
    
    colnames(sample_lookup) <- c("sample", "mut", "reps")
    
    ngene <- max(Mut_data_est$fnum)
    num_conds <- max(Mut_data_est$mut)
    nMT <- num_conds
    nreps <- rep(0, times = num_conds)
    for(i in 1:num_conds){
      nreps[i] <- max(Mut_data_est$reps[Mut_data_est$mut == i])
    }
    
    fast_list <- avg_and_regularize(Mut_data_est, nreps, sample_lookup, 
                                    feature_lookup, NSS = NSS, BDA_model = BDA_model,
                                    ...)
    
    Fit_lists <- list(Fast_Fit = fast_list,
                      Data_lists = data_list)
    
    class(Fit_lists) <- "bakRFnFit"
    return(Fit_lists)
    

    
  }else if(inherits(obj, "bakRFnFit")){
    
    if(FastRerun){
      data_list <- obj$Data_lists
      
      ## Necessary hacky preprocessing
      Mut_data_est <- data_list$Fn_est
      
      colnames(Mut_data_est) <- c("XF", "sample", "fn", "nreads", "fnum", "reps",
                                  "mut", "tl", "logit_fn_rep", "kd_rep_est", "log_kd_rep_est", 
                                  "logit_fn_se", "log_kd_se")
      
      sample_lookup <- data_list$Stan_data$sample_lookup
      feature_lookup <- data_list$Stan_data$sdf
      
      colnames(sample_lookup) <- c("sample", "mut", "reps")
      
      ngene <- max(Mut_data_est$fnum)
      num_conds <- max(Mut_data_est$mut)
      nMT <- num_conds
      nreps <- rep(0, times = num_conds)
      for(i in 1:num_conds){
        nreps[i] <- max(Mut_data_est$reps[Mut_data_est$mut == i])
      }
      
      fast_list <- avg_and_regularize(Mut_data_est, nreps, sample_lookup, 
                                      feature_lookup, NSS = NSS, BDA_model = BDA_model,
                                      ...)
      
      obj$Fast_Fit <- fast_list
      return(obj)
    }else{
      message("Running Hybrid implementation. Cannot run MCMC implementation with bakRFnFit input")
      
      Rep_Fn <- obj$Fast_Fit$Fn_Estimates
      
      data_list <- obj$Data_lists$Stan_data
      
      sample_lookup <- data_list$sample_lookup
      colnames(sample_lookup) <- c("sample", "mut", "reps")
      data_list$sample_lookup <- sample_lookup
      
      
      rm(Rep_Fn)
      
      # Run Hybrid implementation
      Stan_list <- TL_stan(data_list, Hybrid_Fit = TRUE, NSS = NSS, chains = chains, ...)
      
      ## Calculate and adjust p-values
      
      Effects <- Stan_list$Effects_df
      
      Effects <- Effects %>% dplyr::mutate(pval = 2*stats::pnorm(-abs(effect/se))) %>%
        dplyr::mutate(padj = stats::p.adjust(pval, method = "BH"))
      
      Stan_list$Effects_df <- Effects
      
      rm(Effects)
      
      class(Stan_list) <- "HybridModelFit"
      
      obj$Hybrid_Fit <- Stan_list
      
      return(obj)
      
    }
    
    
  }else{
    stop("obj is not of class bakRData, bakRFit, bakRFnData, or bakRFnFit")
  }



}
