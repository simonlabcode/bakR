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
#' @param ... Additoinal (optional) parameters to be passed to \code{stats::nls()}
#' @export
CorrectDropout <- function(obj, ...){
  
  ### Checks
  # 1) Input must be a bakRFit object
  # 2) There must be -s4U controls for all experimental conditions
  
  # Check obj
  if(!inherits(obj, "bakRFit")){
    if(inherits(obj, "bakRData")){
      stop("You provided a bakRData object. You need to run bakRFit before running CorrectDropout")
    }else{
      stop("You did not provide a bakRFit object!")
    }
  }
  
  # Check that -s4U samples exist for all mut
  check <- obj$Data_lists$Fast_df %>%
    dplyr::filter(type == 0) %>%
    dplyr::select(mut) %>%
    dplyr::distinct()
  
  if(!identical(check$mut, 1:obj$Data_lists$Stan_data$nMT)){
    stop("You do not have at least one replicate of -s4U data for all experimental conditions!")
  }
  
  
  # Define helper functions:
  logit <- function(x) log(x/(1-x))
  inv_logit <- function(x) exp(x)/(1+exp(x))
  
  # Compile necessary data -------------------------------------------------------
  
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
  
  # Fit dropout model ------------------------------------------------------------
  
  # Loop over data
  nMT <- Fit$Data_lists$Stan_data$nMT
  nreps <- Fit$Data_lists$Stan_data$nrep_vect
  
  pdos <- rep(0, times = sum(nreps))
  
  count <- 1
  
  # Fit model
  for(i in 1:nMT){
    for(j in 1:nreps[i]){
      fit <- stats::nls(dropout ~ (-(scale*pdo)*fn)/((1-pdo) + fn*pdo) + scale,
                 data = model_df[model_df$reps == j & model_df$mut == i,],
                 start = list(scale = 1, pdo = 0.5),
                 ...)
      
      fit_summ <- summary(fit)
      
      ests <- fit_summ$coefficients
      
      pdos[count] <- ests[rownames(ests) == "pdo","Estimate"]
      if(pdos[count] < 0 ){
        pdos[count] <- 0
      }else if(pdos[count] >= 1{
        stop("Dropout was estimated to be almost 100%, in one of your samples, which is likely an estimation error.
              Is your data of unusually low coverage or metabolic label incorporation rate? This can lead to estimation problems.")
      })
      count <- count + 1
    }
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
  
  # Add dropout info
  do_df <- data.frame(pdo = pdos,
                      Exp_ID = rep(1:nMT, times = nreps) )
  
  do_df <- do_df %>%
    dplyr::group_by(Exp_ID) %>%
    dplyr::mutate(Replicate = as.integer(1:length(pdo)))
  
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
