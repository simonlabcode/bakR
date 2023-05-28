#' Correcting for metabolic labeling induced RNA dropout
#'
#' Dropout is the name given to a phenomenon originally identified by our lab and
#' further detailed in two independent publications (<a href="https://www.biorxiv.org/content/10.1101/2023.05.24.542133v1" >Zimmer et al. (2023)</a>,
#' and <a href="https://www.biorxiv.org/content/10.1101/2023.04.21.537786v1"> Berg et al. (2023)</a>). 
#' Dropout is the under-representation
#' of reads from RNA containing metabolic label (4-thiouridine or 6-thioguanidine most commonly).
#' Loss of 4-thiouridine (s4U) containing RNA on plastic surfaces and RT dropoff caused by 
#' modifications on s4U introduced by recoding chemistry have been attributed as the likely
#' causes of this phenmenon. While protocols can be altered in ways to drastically reduce this
#' source of dropout, you may still have datasets that you want to analyze with bakR collected
#' with suboptimal handling. That is where \code{CorrectDropout} comes in.
#'
#' \code{CorrectDropout} estimates the percentage of 4-thiouridine containing RNA
#' that was lost during library preparation (pdo). It then uses this estimate of pdo
#' to correct fraction new estimates and read counts. Both corrections are analytically
#' derived from a rigorous generative model of NR-seq data. Importantly, the read count
#' correction preserves the total library size to avoid artifically inflating read counts.
#' 
#' @importFrom magrittr %>%
#' @param obj bakRFit object
#' @param scale_init Numeric; initial estimate for -s4U/+s4U scale factor. This is the factor
#' difference in RPM normalized read counts for completely unlabeled transcripts (i.e., highly stable
#' transcript) between the +s4U and -s4U samples.
#' @param pdo_init Numeric; initial estimtae for the dropout rate. This is the probability
#' that an s4U labeled RNA molecule is lost during library prepartion. 
#' @param ... Additional (optional) parameters to be passed to \code{stats::nls()}
#' @export
CorrectDropout <- function(obj, scale_init = 1.05, pdo_init = 0.3,
                           ...){
  
  ### Checks
  # 1) Input must be a bakRFit object
  # 2) There must be -s4U controls for all experimental conditions
  
  # Check obj
  if(!(inherits(obj, "bakRFit") | inherits(obj, "bakRFnFit")) ){
    if(inherits(obj, "bakRData")){
      stop("You provided a bakRData object. You need to run bakRFit before running CorrectDropout")
    }else{
      stop("You did not provide a bakRFit object!")
    }
  }
  
  if(inherits(obj, "bakRFit")){
    # Check that -s4U samples exist for all mut
    check <- obj$Data_lists$Fast_df %>%
      dplyr::filter(type == 0) %>%
      dplyr::select(mut) %>%
      dplyr::distinct()
    
    check <- check$mut
  }else{
    # Check that -s4U samples exist for all mut
    check <- obj$Data_lists$Ctl_data %>%
      dplyr::select(Exp_ID) %>%
      dplyr::distinct()
    
    check <- check$Exp_ID
  }

  
  if(!identical(check, 1:obj$Data_lists$Stan_data$nMT)){
    stop("You do not have at least one replicate of -s4U data for all experimental conditions!")
  }
  
  # Check scale_init
  if(!is.numeric(scale_init)){
    stop("scale_init must be numeric!")
  }else if(scale_init < 1){
    stop("scale_init must be >= 1.")
  }else if(scale_init > 5){
    warning("scale_init is set to an unusually high number. The +s4U/-s4U scale
            factor is likely just over 1")
  }
  
  # Check pdo_init
  if(!is.numeric(pdo_init)){
    stop("pdo_init must be numeric!")
  }else if(pdo_init <= 0){
    stop("pdo_init must be > 0")
  }else if(pdo_init >= 1){
    stop("pdo_init must be < 1")
  }
  
  # Define helper functions:
  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))
  
  if(is.null(obj$Data_lists$Dropout_df)){

    do_df <- QuantifyDropout(obj, ...)
    
  }else{
    do_df <- obj$Data_lists$Dropout_df
  }
  
  
  
  # ggplot2::ggplot(model_df[model_df$reps == 2 & model_df$mut == 1,], ggplot2::aes(x = fn, y = dropout)) + 
  #   ggplot2::geom_point(alpha = 0.1) + 
  #   ggplot2::theme_minimal() + 
  #   ggplot2::xlab("Fraction new") + 
  #   ggplot2::ylab("Dropout") + 
  #   ggplot2::ylim(c(0, 3))
  
  
  # Correct fraction news --------------------------------------------------------
  
  correct_fn <- function(fndo, pdo){
    fndo/((1-pdo) + fndo*pdo)
  }
  
  # Collect biased fn estimates
  Fn_bias <- obj$Fast_Fit$Fn_Estimates
  
  Fn_bias <- Fn_bias[,c("XF", "nreads", "sample", "Exp_ID", "Replicate", "logit_fn")]
  
  Fn_bias$fndo <- inv_logit(Fn_bias$logit_fn)
  
  Fn_bias <- dplyr::inner_join(Fn_bias, do_df, by = c("Exp_ID", "Replicate"))
  
  # Correct
  Fn_bias$fn_corrected <- correct_fn(Fn_bias$fndo, Fn_bias$pdo)
  Fn_bias$logit_fn_corrected <- logit(Fn_bias$fn_corrected)

  
  # Correct read counts ----------------------------------------------------------
  
  # Calculate global fraction new
  FnG <- Fn_bias %>%
    dplyr::group_by(Exp_ID, Replicate) %>%
    dplyr::summarise(fnGdo = sum(nreads*fndo)/sum(nreads),
                     pdo = mean(pdo)) %>%
    dplyr::mutate(fnG = correct_fn(fnGdo, pdo))
  
  # Add global fraction new info to Fn_bias
  Fn_bias <- dplyr::inner_join(Fn_bias, FnG[,c("Exp_ID", "Replicate", "fnG")], by = c("Exp_ID", "Replicate"))
  
  
  # Correct read counts
  correct_reads <- function(reads, fn, fnG, pdo){
    reads_c <- reads*((1 - fnG*pdo)/(1 - fn*pdo))
  }
  
  Fn_bias <- Fn_bias %>%
    dplyr::mutate(reads_corrected = round(correct_reads(nreads, fn_corrected, fnG, pdo)))
  
  
  # Rerun MLE, going thru bakRFnData ---------------------------------------------
  
  
  if(inherits(obj, "bakRFit")){
    # Prep fns
    Fns <- Fn_bias[,c("XF", "sample", "reads_corrected", "fn_corrected")]
    colnames(Fns) <- c("XF", "sample", "n", "fn")
    
    
    # Make metadf
    metadf <- obj$Data_lists$Fast_df %>%
      dplyr::select(sample, mut, tl, type) %>%
      dplyr::mutate(tl = ifelse(type == 0, 0, tl)) %>%
      dplyr::select(sample, mut, tl) %>%
      dplyr::distinct()
    
    metadf <- data.frame(Exp_ID = metadf$mut,
                         tl = metadf$tl,
                         row.names = metadf$sample)
    
    
    # Add -s4U data to Fns
    ctl_samps <- rownames(metadf[metadf$tl == 0,])
    
    XF <- unique(Fns$XF)
    Fn_ctl <- data.frame(fn = 0,
                         XF = rep(XF, times = length(ctl_samps)),
                         sample = rep(ctl_samps, each = length(XF)))
    
    ctl_reads <- obj$Data_lists$Fast_df %>%
      dplyr::filter(type == 0) %>%
      dplyr::group_by(XF, sample) %>%
      dplyr::summarise(n = sum(n))
    
    Fn_ctl <- dplyr::inner_join(Fn_ctl, ctl_reads, by = c("XF", "sample"))
    Fns <- dplyr::bind_rows(list(Fns, Fn_ctl[,c("XF", "sample", "fn", "n")]))
    
  }else{
    # Prep fns
    Fns <- Fn_bias[,c("XF", "sample", "reads_corrected", "fn_corrected")]
    colnames(Fns) <- c("XF", "sample", "n", "fn")
    
    # Make metadf
    meta_s4U <- Fit$Data_lists$Fn_est %>%
      dplyr::select(Exp_ID, tl, sample) %>%
      dplyr::distinct()
    meta_ctl <- Fit$Data_lists$Ctl_data %>%
      dplyr::select(Exp_ID, tl, sample) %>%
      dplyr::distinct()
    
    meta_tibble <- dplyr::bind_rows(meta_s4U, meta_ctl)
    
    metadf <- data.frame(Exp_ID = meta_tibble$Exp_ID,
                         tl = meta_tibble$tl,
                         row.names = meta_tibble$sample)
    
    # Add -s4U data to Fns
    ctl_samps <- rownames(metadf[metadf$tl == 0,])
    
    XF <- unique(Fns$XF)
    Fn_ctl <- data.frame(fn = 0,
                         XF = rep(XF, times = length(ctl_samps)),
                         sample = rep(ctl_samps, each = length(XF)))
    
    ctl_reads <- Fit$Data_lists$Ctl_data %>%
      dplyr::group_by(XF, sample) %>%
      dplyr::summarise(n = sum(n))
    
    Fn_ctl <- dplyr::inner_join(Fn_ctl, ctl_reads, by = c("XF", "sample"))
    Fns <- dplyr::bind_rows(list(Fns, Fn_ctl[,c("XF", "sample", "fn", "n")]))
    
    
  }

  
  # Make bakRFnData
  bakRFnData <- bakRFnData(Fns, metadf)
  
  # Rerun fit with corrected reads and fns
  Fit_correct <- bakRFit(bakRFnData, FOI = unique(Fns$XF), concat = FALSE)
  
  
  # Make Final Fit object --------------------------------------------------------
  
  Fit_Final <- obj
  Fit_Final$Fast_Fit <- Fit_correct$Fast_Fit
  Fit_Final$Data_lists$Count_Matrix_corrected <- Fit_correct$Data_lists$Count_Matrix
  Fit_Final$Data_lists$Dropout_df <- do_df
  Stan_data <- obj$Data_lists$Stan_data
  
  Stan_data$Avg_Reads <- Fit_correct$Data_lists$Stan_data$Avg_Reads
  Stan_data$Avg_Reads_natural <- Fit_correct$Data_lists$Stan_data$Avg_Reads_natural
  Fit_Final$Data_lists$Stan_data <- Stan_data
  
  return(Fit_Final)
  
  
}

#' Fit dropout model to quantify dropout frequency
#'
#' \code{QuantifyDropout} estimates the percentage of 4-thiouridine containing RNA
#' that was lost during library preparation (pdo).
#' @importFrom magrittr %>%
#' @param obj bakRFit object
#' @param keep_data Logical; if TRUE, will return list with two elements. First element
#' is the regular return (data frame with dropout quantitied), and the second element
#' will be the data frame that was used for fitting the dropout model. This is useful
#' if wanting to visualize the fit.
#' @param scale_init Numeric; initial estimate for -s4U/+s4U scale factor. This is the factor
#' difference in RPM normalized read counts for completely unlabeled transcripts (i.e., highly stable
#' transcript) between the +s4U and -s4U samples.
#' @param pdo_init Numeric; initial estimtae for the dropout rate. This is the probability
#' that an s4U labeled RNA molecule is lost during library prepartion. 
#' @param no_message Logical; if TRUE, will not output message regarding estimated
#' rates of dropout in each sample
#' @param ... Additional (optional) parameters to be passed to \code{stats::nls()}
#' @export
QuantifyDropout <- function(obj, scale_init = 1.05, pdo_init = 0.3,
                            keep_data = FALSE, no_message = FALSE,
                            ...){
  ### Checks
  # 1) Input must be a bakRFit object
  # 2) There must be -s4U controls for all experimental conditions
  # 3) Make sure scale_init is > 0 and close to 1
  # 4) Make sure pdo_init is between 0 and 1
  
  # Check obj
  if(!(inherits(obj, "bakRFit") | inherits(obj, "bakRFnFit")) ){
    if(inherits(obj, "bakRData")){
      stop("You provided a bakRData object. You need to run bakRFit before running CorrectDropout")
    }else{
      stop("You did not provide a bakRFit object!")
    }
  }
  
  if(inherits(obj, "bakRFit")){
    # Check that -s4U samples exist for all mut
    check <- obj$Data_lists$Fast_df %>%
      dplyr::filter(type == 0) %>%
      dplyr::select(mut) %>%
      dplyr::distinct()
    
    check <- check$mut
  }else{
    # Check that -s4U samples exist for all mut
    check <- obj$Data_lists$Ctl_data %>%
      dplyr::select(Exp_ID) %>%
      dplyr::distinct()
    
    check <- check$Exp_ID
  }
  
  
  if(!identical(check, 1:obj$Data_lists$Stan_data$nMT)){
    stop("You do not have at least one replicate of -s4U data for all experimental conditions!")
  }
  
  # Check scale_init
  if(!is.numeric(scale_init)){
    stop("scale_init must be numeric!")
  }else if(scale_init < 1){
    stop("scale_init must be >= 1.")
  }else if(scale_init > 5){
    warning("scale_init is set to an unusually high number. The +s4U/-s4U scale
            factor is likely just over 1")
  }
  
  # Check pdo_init
  if(!is.numeric(pdo_init)){
    stop("pdo_init must be numeric!")
  }else if(pdo_init <= 0){
    stop("pdo_init must be > 0")
  }else if(pdo_init >= 1){
    stop("pdo_init must be < 1")
  }
  
  
  # Define helper functions:
  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))
  
  # Compile necessary data -------------------------------------------------------
  
  if(inherits(obj, "bakRFit")){
    # Calculate number of reads in each sample
    total_reads <- obj$Data_lists$Fast_df
    total_reads <- total_reads %>%
      dplyr::group_by(mut, reps, type) %>%
      dplyr::summarise(total_reads = sum(n))
    
    # Calculate number of reads for each feature in each sample
    RPMs <- obj$Data_lists$Fast_df
    RPMs <- RPMs %>%
      dplyr::group_by(fnum, mut, reps, type) %>%
      dplyr::summarise(reads = sum(n))
    
    # RPM normalize
    RPMs <- dplyr::inner_join(RPMs, total_reads, by = c("mut", "reps", "type"))
    RPMs$RPM <- RPMs$reads/(RPMs$total_reads/1000000)
    
    # Separate control from +s4U
    ctl_RPMs <- RPMs[RPMs$type == 0,]
    colnames(ctl_RPMs) <- c("fnum", "mut", "reps", "type", "reads", "total_reads", "ctl_RPM")
    ctl_RPMs <- ctl_RPMs[,c("fnum", "mut", "reps", "ctl_RPM")] %>%
      dplyr::group_by(fnum, mut) %>%
      dplyr::summarise(ctl_RPM = mean(ctl_RPM))
    
    s4U_RPMs <- RPMs[RPMs$type == 1,c("fnum", "mut", "reps", "RPM")]
    
    # Combine +s4U and -s4U RPMs
    s4U_RPMs <- dplyr::inner_join(s4U_RPMs, ctl_RPMs, by = c("fnum", "mut"))
    
    # Get fraction new estimates
    Fns <- obj$Fast_Fit$Fn_Estimates
    
    Fns <- Fns[,c("Feature_ID", "Exp_ID", "Replicate", "logit_fn", "XF", "sample")]
    
    colnames(Fns) <- c("fnum", "mut", "reps", "logit_fn", "XF", "sample")
    
    Fns$fn <- inv_logit(Fns$logit_fn)
    
    # Add fn info and dropout calc
    model_df <- dplyr::inner_join(s4U_RPMs, Fns, by = c("fnum", "mut", "reps"))
    model_df$dropout <- model_df$RPM/model_df$ctl_RPM
    
  }else{
    # Calculate number of reads in each -s4U sample
    total_reads_ctl <- Fit$Data_lists$Ctl_data %>%
      dplyr::group_by(Exp_ID, Replicate) %>%
      dplyr::summarise(total_reads = sum(n))
    total_reads_ctl$type <- 0
    
    # Calculate number of reads in each +s4U samle
    total_reads <- Fit$Data_lists$Fn_est %>%
      dplyr::group_by(Exp_ID, Replicate) %>%
      dplyr::summarise(total_reads = sum(n))
    total_reads$type <- 1
    
    total_reads <- dplyr::bind_rows(list(total_reads, total_reads_ctl))
    
    # Calculate number of reads for each feature in each sample
    RPMs <- Fit$Data_lists$Fn_est
    RPMs <- RPMs %>%
      dplyr::group_by(Feature_ID, Exp_ID, Replicate) %>%
      dplyr::summarise(reads = sum(n))
    RPMs$type <- 1
    
    # RPM normalize
    RPMs <- dplyr::inner_join(RPMs, total_reads, by = c("Exp_ID", "Replicate", "type"))
    RPMs$RPM <- RPMs$reads/(RPMs$total_reads/1000000)
    
    # Calculate number of reads for each feature in each -s4U sample
    ctl_RPMs <- Fit$Data_list$Ctl_data
    ctl_RPMs <- ctl_RPMs %>%
      dplyr::group_by(Feature_ID, Exp_ID, Replicate) %>%
      dplyr::summarise(reads = sum(n))
    ctl_RPMs$type <- 0
    
    
    # RPM normalize
    ctl_RPMs <- dplyr::inner_join(ctl_RPMs, total_reads, by = c("Exp_ID", "Replicate", "type"))
    ctl_RPMs$ctl_RPM <- ctl_RPMs$reads/(ctl_RPMs$total_reads/1000000)
    
    # Average over -s4U replicates
    ctl_RPMs <- ctl_RPMs %>%
      dplyr::group_by(Feature_ID, Exp_ID) %>%
      dplyr::summarise(ctl_RPM = mean(ctl_RPM))
    
    # Combine +s4U and -s4U RPMs
    s4U_RPMs <- dplyr::inner_join(RPMs[,c("Feature_ID", "Exp_ID", "Replicate", "RPM")], 
                                  ctl_RPMs, by = c("Feature_ID", "Exp_ID"))
    colnames(s4U_RPMs) <- c("fnum", "mut", "reps", "RPM", "ctl_RPM")
    
    # Get fraction new estimates
    Fns <- Fit$Fast_Fit$Fn_Estimates
    
    Fns <- Fns[,c("Feature_ID", "Exp_ID", "Replicate", "logit_fn", "XF", "sample")]
    
    colnames(Fns) <- c("fnum", "mut", "reps", "logit_fn", "XF", "sample")
    
    Fns$fn <- inv_logit(Fns$logit_fn)
    
    # Add fn info and dropout calc
    model_df <- dplyr::inner_join(s4U_RPMs, Fns, by = c("fnum", "mut", "reps"))
    model_df$dropout <- model_df$RPM/model_df$ctl_RPM
    
  }
  
  
  # Fit dropout model ------------------------------------------------------------
  
  # Loop over data
  nMT <- Fit$Data_lists$Stan_data$nMT
  nreps <- Fit$Data_lists$Stan_data$nrep_vect
  
  pdos <- rep(0, times = sum(nreps))
  scales <- pdos
  
  count <- 1
  
  # Fit model
  for(i in 1:nMT){
    for(j in 1:nreps[i]){
      fit <- stats::nls(dropout ~ (-(scale*pdo)*fn)/((1-pdo) + fn*pdo) + scale,
                        data = model_df[model_df$reps == j & model_df$mut == i,],
                        start = list(scale = scale_init, pdo = pdo_init),
                        ...)
      
      fit_summ <- summary(fit)
      
      ests <- fit_summ$coefficients
      
      pdos[count] <- ests[rownames(ests) == "pdo","Estimate"]
      scales[count] <- ests[rownames(ests) == "scale","Estimate"]
      
      if( pdos[count] < 0 ){
        pdos[count] <- 0
      }else if( pdos[count] >= 1){
        stop("Dropout was estimated to be almost 100%, in one of your samples, which is likely an estimation error.
              Is your data of unusually low coverage or metabolic label incorporation rate? This can lead to estimation problems.")
      }
      count <- count + 1
    }
  }
  
  # Add dropout info
  do_df <- data.frame(pdo = pdos,
                      scale = scales,
                      Exp_ID = rep(1:nMT, times = nreps) )
  
  do_df <- do_df %>%
    dplyr::group_by(Exp_ID) %>%
    dplyr::mutate(Replicate = as.integer(1:length(pdo)))
  
  if(!no_message){
    message(paste0(c("Estimated rates of dropout are:", 
                     utils::capture.output(as.data.frame(do_df[,c("Exp_ID", "Replicate", "pdo")]))),
                   collapse = "\n"))
  }
  
  if(keep_data){
    output <- list(Dropout_df = do_df,
                   Input_data = model_df)
  }else{
    return(do_df)
  }
  
  
  
}

