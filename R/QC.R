#' Check data quality and make suggestions to user about what analyses to run.
#'
#' \code{QC_checks} takes as input a bakRFit object and uses the Fast_Fit object to assess
#' data quality and make suggestions about which implementation to run next. QC_checks
#' takes into account the mutation rates in all samples, the fraction new distributions, the reproducibility
#' of fraction new estimates, and the read lengths. It then outputs a number of
#' diagnostic plots that might alert users to problems in their data. It also
#' outputs messages informing users what implementation is best used next.
#'
#' @param obj bakRFit object
#' @return A list with 3 components:
#' \itemize{
#'  \item raw_mutrates. This is a plot of the raw T-to-C mutation rates in all samples
#'  analyzed by bakR. It includes horizontal lines as reference for what could be
#'  considered "too low" to be useful in s4U fed samples.
#'  \item conversion_rates. This is a plot of the estimated T-to-C mutation rates
#'  in new and old reads. Thus, each bar represents the probability that a U in
#'  a new/old read is mutated. It includes horizontal lines as reference for what could
#'  be considered good mutation rates.
#'  \item correlation_plots. This is a list of ggplot objects. Each is a scatter plot
#'  comparing estimates of the fraction new in one replicate to another replicate
#'  in the same experimental condition. A y=x guide line is included to reveal any
#'  estimation biases.
#' }
#' @examples
#' \donttest{
#' # Simulate data for 500 genes and 2 replicates
#' sim <- Simulate_bakRData(500, nreps = 2)
#'
#' # Fit data with fast implementation
#' Fit <- bakRFit(sim$bakRData)
#'
#' # Run QC
#' QC <- QC_checks(Fit)
#'
#' }
#' @importFrom magrittr %>%
#' @export
QC_checks <- function(obj){
  

  Exp_ID <- Replicate <- logit_fn <- fn_1 <- fn_2 <- NULL
  type <- mut <- reps <- pnew <- mutrate <- TC <- n <- nT <- ctl <- NULL
  
  ## Function for assessing fraction new correlation
  assess_fn_cor <- function(obj2, Bad_data, bakRFn){
    
    ### Assess fraction new distribution
    Fns <- obj2$Fn_Estimates
    
    # Calculate average fraction new in each sample
    avg_fns <- Fns %>%
      dplyr::group_by(Exp_ID, Replicate) %>%
      dplyr::summarise(avg_logit_fn = mean(logit_fn))
    
    message(paste0(c("Average logit(fraction news) for each sample are:", utils::capture.output(avg_fns)), collapse = "\n"))
    
    message("Reminder: a logit fraction new of 0 means a fraction new of 0.5, which would be ideal.")
    
    avg_fns <- avg_fns$avg_logit_fn
    
    if(all(dplyr::between(avg_fns, -2, 2))){
      message("The average logit(fraction news) in all samples are between -2 and 2, suggesting an appropriate label time!")
    }
    
    if(any(dplyr::between(avg_fns, -4, -2) )){
      warning("The average logit(fraction news) are relatively low (between -4 and -2) in one or more samples, suggesting your label time was a bit short. This will limit bakR's ability to identify kinetic differences")
      
      if(!bakRFn){
        message("Low fraction news impair bakR's default mutation rate estimation strategy. I suggest rerunning bakRFit with FastRerun and StanRateEst = TRUE, particularly if some of the estimated mutation rates are oddly low (< 0.01) in a subset of samples.")
      }
      
    }
    
    if(any(dplyr::between(avg_fns, 2, 4))){
      warning("The average logit(fraction news) are relatively high (between 2 and 4) in one or more samples, suggesting your label time was a big long. This will limit bakR's ability to identify kinetic differences")
    }
    
    if(any(avg_fns < -4)){
      warning("The average logit(fraction news) are extremely low (less than -4) in one or more samples, suggesting your label time was too short. It will be difficult for bakR to identify any kinetic differences.")
      if(!bakRFn){
        message("Low fraction news impair bakR's default mutation rate estimation strategy. I suggest rerunning bakRFit with FastRerun and StanRateEst = TRUE, particularly if some of the estimated mutation rates are oddly low (< 0.01) in a subset of samples.")
      }
      Bad_data <- TRUE
    }
    
    if(any(avg_fns > 4)){
      warning("The average logit(fraction news) are extremely high (greater than 4) in one or more samples, suggesting your label time was too long. It will be difficult for bakR to identify any kinetic differences.")
      Bad_data <- TRUE
    }
    
    
    ### Assess fraction new correlation
    
    # How many replicates in each Exp_ID?
    nreps <- Fns %>%
      dplyr::group_by(Exp_ID) %>%
      dplyr::summarise(nreps = max(Replicate)) %>%
      dplyr::select(nreps)
    
    nreps <- nreps$nreps
    
    # calculate correlations between each set of replicates
    ncalcs <- sum(choose(nreps, 2))
    
    Exps <- rep(0, times = ncalcs)
    Rep_ID1 <- Exps
    Rep_ID2 <- Exps
    fn_cors <- Exps
    cor_plots <- vector(mode = "list", length = ncalcs)
    
    count <- 1
    for(i in 1:length(nreps)){
      
      
      for(j in 1:(nreps[i]-1)){
        for(k in (j+1):nreps[i]){
          
          Exps[count] <- i
          
          
          fn_cors[count] <- stats::cor(Fns$logit_fn[Fns$Exp_ID == i & Fns$Replicate == j],
                                       Fns$logit_fn[Fns$Exp_ID == i & Fns$Replicate == k])
          
          cor_df <- data.frame(fn_1 = Fns$logit_fn[Fns$Exp_ID == i & Fns$Replicate == j],
                               fn_2 = Fns$logit_fn[Fns$Exp_ID == i & Fns$Replicate == k])
          
          coeff <- stats::cor(cor_df$fn_1,
                              cor_df$fn_2)
          
          npoints <- nrow(cor_df)
          k <- 0.75
          alpha <- exp(-(log10(npoints) - 1)*k)
          if(alpha > 1){
            alpha <- 1
          }
          
          cor_plots[[count]] <- ggplot2::ggplot(cor_df, ggplot2::aes(x = fn_1, y = fn_2)) +
            ggplot2::geom_point(alpha = alpha) +
            ggplot2::xlab(paste0("logit(fn) replicate " , j, ", condition ", i)) +
            ggplot2::ylab(paste0("logit(fn) replicate " , k, ", condition ", i)) +
            ggplot2::ggtitle(paste0("logit(fn) correlation (p = ", round(coeff, digits = 3), ")"),
                             subtitle = "Points ideally closely follow red y = x line") +
            ggplot2::geom_abline(slope = 1, intercept = 0, color = "red") +
            ggplot2::theme_classic()
          
          Rep_ID1[count] <- j
          Rep_ID2[count] <- k
          
          count <- count + 1
        }
      }
      
      
    }
    
    
    ## Make correlation matrix
    nr <- max(Fns$Feature_ID)
    nc <- length(unique(Fns$sample))
    
    fn_mat <- matrix(0, nrow = nr, ncol = nc)
    
    samps <- unique(Fns$sample)
    
    count <- 1
    for(i in samps){
      fn_mat[,count] <- Fns$logit_fn[Fns$sample == i]
      count <- count + 1
    }
    
    fn_cor_mat <- stats::cor(fn_mat)
    
    
    
    fn_cors <- data.frame(Exp_ID = Exps,
                          Rep_ID1 = Rep_ID1,
                          Rep_ID2 = Rep_ID2,
                          correlation = fn_cors)
    
    message(paste0(c("logit(fn) correlations for each pair of replicates are:", utils::capture.output(fn_cors)), collapse = "\n"))
    
    if(any(fn_cors$correlation < 0.7)){
      warning("logit(fraction new) correlation is low in one or more samples. Did you properly identify replicates in the metadf of your bakRData object?")
    }else{
      message("logit(fn) correlations are high, suggesting good reproducibility!")
    }
    
    out <- list(correlation_plots = cor_plots,
                correlation_matrix = fn_cor_mat,
                Bad_data = Bad_data)
    
  }
  
  ### Extract Fast_Fit to be used for diagnostic analyses
  Fit <- obj$Fast_Fit
  
  MCMC_next <- FALSE
  
  Bad_data <- FALSE
  
  if(inherits(obj, "bakRFnFit")){
    
    fn_assessment <- assess_fn_cor(Fit, bakRFn = TRUE, Bad_data = Bad_data)
    
    Bad_data <- fn_assessment$Bad_data
    
    ### Make suggestions
    if(Bad_data){
      message("Some aspects of your data may limit bakR's ability to detect differential kinetics. Check warning messages for details.")
    }else{
      message("I suggest running the Hybrid implementation next. This can be done with bakRFit(Fit, HybridFit = TRUE), where Fit is your bakRFit object.")
    }
    
    glist <- list(correlation_plots = fn_assessment$correlation_plots,
                  correlation_matrix = fn_assessment$correlation_matrix)
    
    return(glist)
    
  }else{
    
    ### Assess mutation rates
    Mutation_Rates <- Fit$Mut_rates
    
    # Mutation rates
    mutrates <- Mutation_Rates
    
    if(all(mutrates$pnew > 0.02)){
      message("Mutation rates in new reads looks good!")
    }
    
    if(any(dplyr::between(mutrates$pnew, 0.0099, 0.0201))){
      warning("Mutation rates in new reads is somewhat low in one or more samples.")
      
      MCMC_next <- TRUE
    }
    
    if(any(dplyr::between(mutrates$pnew, 0.0069, 0.01))){
      warning("Mutation rates in new reads are below 1% one or more samples. This significanlty reduces bakR's ability to identify differential kinetics.")
      
      MCMC_next <- TRUE
      
    }
    
    if(any(mutrates$pnew < 0.007)){
      warning("Mutation rates in new reads are below 0.7% in one or more samples. It is very difficult to identify kinetic differences with such low mutation rates.")
      
      Bad_data <- TRUE
    }
    
    # polds
    if(all(mutrates$pold < 0.004)){
      message("Background mutation rate looks good!")
    }else if(all(mutrates$pold < 0.01)){
      warning("Background mutation rate is a bit high. Did you account for SNPs when counting mutations?")
    }else{
      warning("Background mutation rate is high (>= 1%). Did you properly identify -s4U control samples in the metadf of your bakRData object?")
      
      Bad_data <- TRUE
    }
    
    fn_assessment <- assess_fn_cor(Fit, Bad_data = Bad_data, 
                  bakRFn = FALSE)
    
    Bad_data <- fn_assessment$Bad_data
    
    
    ### Make suggestions
    if(Bad_data){
      message("Some aspects of your data may limit bakR's ability to detect differential kinetics. Check warning messages for details.")
    }else if(MCMC_next){
      message("Given your low mutation rates, I suggest running the MCMC implementation next. This can be done with bakRFit(Fit, StanFit = TRUE), where Fit is your bakRFit object.")
    }else{
      message("I suggest running the Hybrid implementation next. This can be done with bakRFit(Fit, HybridFit = TRUE), where Fit is your bakRFit object.")
    }
    
    
    ### Create visualizations
    
    # U-to-C mutation rates
    
    muts <- Fit$Mut_rates
    
    
    pnews <- muts[,c("mut", "reps", "pnew")]
    
    fast_df <- obj$Data_lists$Fast_df
    
    samp_ID <- fast_df %>% dplyr::filter(type == 1) %>%
      dplyr::select(sample, mut, reps) %>%
      dplyr::distinct()
    
    pnews <- dplyr::left_join(pnews, samp_ID, by = c("mut", "reps"))
    
    pnews <- pnews %>%
      dplyr::mutate(mutrate = as.factor("new"))
    
    polds <- data.frame(mut = pnews$mut,
                        reps = pnews$reps,
                        pnew = muts$pold,
                        sample = pnews$sample,
                        mutrate = as.factor("old"))
    
    pnews <- dplyr::bind_rows(pnews, polds)
    
    pnews$sample <- as.factor(pnews$sample)
    
    # Pretty plotting theme
    theme_mds <-    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                   panel.grid.minor = ggplot2::element_blank(),
                                   
                                   panel.background = ggplot2::element_blank(),
                                   
                                   axis.line.x = ggplot2::element_line(colour = "black"),
                                   
                                   axis.line.y = ggplot2::element_line(colour = "black"),
                                   
                                   axis.ticks = ggplot2::element_line(colour = "black"),
                                   
                                   title = ggplot2::element_text(color = "black", size = 10),
                                   
                                   axis.text = ggplot2::element_text(color="black", size = 10),
                                   
                                   axis.title = ggplot2::element_text(color = "black", size = 12),
                                   
                                   strip.background = ggplot2::element_blank())
    
    
    
    g_conversion <- ggplot2::ggplot(pnews, ggplot2::aes(x = sample, y = pnew, fill = mutrate)) +
      theme_mds + 
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("Mutation rate") +
      ggplot2::ggtitle("New (red) and old (gray) read mutation rates",
                       subtitle = "Red bars ideally above blue line; Gray pold bar ideally below black line") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
      ggplot2::scale_fill_manual(values = c("#A1121B", "darkgray")) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = 0.02, linetype = "dotted", size = 1.5, color = "blue") +
      ggplot2::geom_hline(yintercept = 0.004, linetype = "dotted", size = 1.5, color = "black")
    
    # Raw mutation rates
    fast_df <- obj$Data_lists$Fast_df
    
    pnews <- fast_df %>%
      dplyr::group_by(sample, type) %>%
      dplyr::summarise(mutrate = sum(TC*n)/sum(nT*n))
    
    pnews <- pnews %>%
      dplyr::mutate(ctl = as.factor(ifelse(type == 1, "labeled", "unlabeled")))
    
    pnews$sample <- as.factor(pnews$sample)
    
    # Pretty plotting theme
    theme_mds <-    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                                   panel.grid.minor = ggplot2::element_blank(),
                                   
                                   panel.background = ggplot2::element_blank(),
                                   
                                   axis.line.x = ggplot2::element_line(colour = "black"),
                                   
                                   axis.line.y = ggplot2::element_line(colour = "black"),
                                   
                                   axis.ticks = ggplot2::element_line(colour = "black"),
                                   
                                   title = ggplot2::element_text(color = "black", size = 10),
                                   
                                   axis.text = ggplot2::element_text(color="black", size = 10),
                                   
                                   axis.title = ggplot2::element_text(color = "black", size = 12),
                                   
                                   strip.background = ggplot2::element_blank())
    
    
    
    
    g_raw <- ggplot2::ggplot(pnews, ggplot2::aes(x = sample, y = mutrate, fill = ctl)) +
      theme_mds + 
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("Raw mutation rate") +
      ggplot2::ggtitle("Raw mutation rates (gray = -s4U; red = +s4U)",
                       subtitle = "Gray bars ideally below black line. Red bars ideally well above black line.") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
      ggplot2::scale_fill_manual(values = c("#A1121B", "darkgray")) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::geom_hline(yintercept = 0.004, linetype = "dotted", size = 1.5, color = "black")
    
    
    # Fraction new correlations
    
    glist <- list(raw_mutrates = g_raw,
                  conversion_rates = g_conversion,
                  correlation_plots = fn_assessment$correlation_plots,
                  correlation_matrix = fn_assessment$correlation_matrix)
    
    return(glist)
  }
  
  

}
