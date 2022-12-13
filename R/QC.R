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
#' @importFrom magrittr %>%
#' @export
QC_checks <- function(obj){

  ### Extract Fast_Fit to be used for diagnostic analyses
  Fit <- obj$Fast_Fit

  ### Assess mutation rates
  Mutation_Rates <- Fit$Mut_rates

  # pnews
  pnews <- Mutation_Rates[[1]]

  if(all(pnews$pnew > 0.02)){
    message("Mutation rates in new reads looks good!")
  }

  if(any(dplyr::between(pnews$pnew, 0.0099, 0.0201))){
    warning("Mutation rates in new reads is somewhat low in one or more samples.")
  }

  if(any(dplyr::between(pnews$pnew, 0.0069, 0.01))){
    warning("Mutation rates in new reads are below 1% one or more samples. This significanlty reduces bakR's ability to identify differential kinetics.")
  }

  if(any(pnews$pnew < 0.007)){
    warning("Mutation rates in new reads are below 0.7% in one or more samples. It is nearly impossible to identify kinetic differences with such low mutation rates.")
  }

  # polds
  polds <- Mutation_Rates[[2]]

  if(polds < 0.004){
    message("Background mutation rate looks good!")
  }else if(polds < 0.01){
    warning("Background mutation rate is a bit high. Did you account for SNPs when counting mutations?")
  }else{
    warning("Background mutation rate is extremely high. Did you properly identify -s4U control samples in the metadf of your bakRData object?")
  }


  ### Assess fraction new distribution
  Fns <- Fit$Fn_Estimates

  # Calculate average fraction new in each sample
  avg_fns <- Fns %>%
    dplyr::group_by(Exp_ID, Replicate) %>%
    dplyr::summarise(avg_logit_fn = mean(logit_fn))

  message(paste0(c("Average logit(fraction news) for each sample are:", utils::capture.output(avg_fns)), collapse = "\n"))

  message("Reminder: a logit fraction new of 0 means a fraction new of 0.5, which would be ideal.")

  avg_fns <- avg_fns$avg_logit_fn

  if(all(dplyr::between(avg_fns, -2, 2))){
    message("Fraction news look good, suggesting an appropriate label time!")
  }

  if(any(dplyr::between(avg_fns, -4, -2) )){
    warning("Fraction news are low in one or more samples, suggesting your label time was a bit short. This will limit bakR's ability to identify kinetic differences")
  }

  if(any(dplyr::between(avg_fns, 2, 4))){
    warning("Fraction news are high in one or more samples, suggesting your label time was a big long. This will limit bakR's ability to identify kinetic differences")
  }

  if(any(avg_fns < -4)){
    warning("Fraction news are extremely low in one or more samples, suggesting your label time was too short. It will be difficult for bakR to identify any kinetic differences.")
  }

  if(any(avg_fns > 4)){
    warning("Fraction news are extremely high in one or more samples, suggesting your label time was too long. It will be difficult for bakR to identify any kinetic differences.")
  }


  ### Assess fraction new correlation

  # How many replicates in each Exp_ID?
  nreps <- Fns %>%
    dplyr::group_by(Exp_ID) %>%
    dplyr::summarise(nreps = max(Replicate)) %>%
    dplyr::select(nreps)

  nreps <- nreps$nreps

  # calculate correlations between each set of replicates
  Exps <- c()
  Rep_ID1 <- c()
  Rep_ID2 <- c()
  fn_cors <- c()

  count <- 1
  for(i in 1:length(nreps)){

    Exps[count] <- i

    for(j in 1:(nreps[i]-1)){
      for(k in (j+1):nreps[i]){
        fn_cors[count] <- cor(Fns$logit_fn[Fns$Exp_ID == i & Fns$Replicate == j],
                              Fns$logit_fn[Fns$Exp_ID == i & Fns$Replicate == k])

        Rep_ID1[count] <- j
        Rep_ID2[count] <- k

        count <- count + 1
      }
    }


  }

  fn_cors <- data.frame(Exp_ID = Exps,
                        Rep_ID1 = Rep_ID1,
                        Rep_ID2 = Rep_ID2,
                        correlation = fn_cors)

  message(paste0(c("logit(fn) correlations for each pair of replicates are:", utils::capture.output(fn_cors)), collapse = "\n"))

  if(any(fn_cors$correlation < 0.8)){
    warning("logit(fraction new) correlation is low in one or more samples. Did you properly identify replicates in the metadf of your bakRData object?")
  }else{
    message("logit(fn) correlations are high, suggesting good reproducibility!")
  }


  ### Create visualizations

}
