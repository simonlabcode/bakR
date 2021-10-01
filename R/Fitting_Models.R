#' Estimating kinetic parameters from nucleotide recoding data
#'
#' \code{DynamicFit} analyzes nucleotide recoding data to estimate
#' kinetic parameters relating to RNA stability and changes in RNA
#' stability induced by experimental perturbations. Several statistical
#' models of varying efficiency and performance can be used to fit data.
#'
#' If \code{DynamicFit} is run on a DynamicSeqData object, \code{cBprocess}
#' and then \code{fast_analysis} will always be called. The former will generate the processed
#' data that can be passed to the model fitting functions (\code{fast_analysis}, \code{Hybrid_fit},
#' and \code{TL_Stan}). The call to \code{fast_analysis} will generate a list of dataframes
#' containing information regarding the \code{fast_analysis} fit. \code{fast_analysis} is always
#' called because its output is required for both \code{Hybrid_fit} and \code{TL_Stan}. In both cases,
#' \code{fast_analysis} estimates a parameter that determines the number of degress of freedom of
#' the t-distribution to use when performing a moderated t-test for significance of L2FC(kdeg)s.
#' For \code{Hybrid_fit}, the logit(fraction new) replicate estimates and Fisher Information derived
#' approximate uncertainties are passed as data to a model that performs partial pooling approximated
#' by the procedure implemented in \code{fast_analysis}.
#'
#' If \code{DynamicFit} is run on a DynamicSeqFit object, \code{cBprocess} will not be called again,
#' as the output of \code{cBprocess} will already be contained in the DynamicSeqFit object. Similarly,
#' \code{fast_analysis} will not be called again unless DynamicSeqFit is rerun on the DynamicSeqData object. If you want to
#' generate model fits using different parameters for cBprocess, you will have to rerun \code{DynamicFit}
#' on the DynamicSeqData object.
#'
#' See the doumentation for the indvidual fitting functions for details regarding how they analyze nucleotide
#' recoding data. What follows is a brief overview of how each works
#'
#' \code{fast_analysis} either estimates mutation rates from + and - s4U samples or uses mutation rate estimates
#' provided by the user to perform maximum likelihood estimation (MLE) of the fraction of RNA that is labeled for each
#' replicate of nucleotide recoding data provided. Uncertainties for each replicate's estimate are approximated using
#' asymptotic results involving the Fisher Information and assuming known mutation rates. Replicate data
#' is pooled using an approximation to hierarchical modeling that relies on analytic solutions to simple Bayesian models.
#' Linear regression is used to estimate the relationship between read depths and replicate variability for uncertainty
#' estimation regularization, again performed using analytic solutions to Bayesian models.
#'
#' \code{Hybrid_fit} takes as input estimates of the logit(fraction new) and uncertainty provided by \code{fast_analysis}.
#' It then uses Stan on the backend to implement a hierarchical modeling that pools data across replicates and the dataset
#' to estimate effect sizes (L2FC(kdeg)) and uncertainties. Replicate variability information is pooled across each experimental
#' condition to regularize variance estimates using a hierarchical linear regression model.
#'
#' \code{TL_Stan} uses Stan on the back end to implement a model similar to \code{Hybrid_fit}. The difference is that while
#' \code{Hybrid_fit} requires fraction new estimates from \code{fast_analysis}, \code{TL_Stan} implements a U-content exposure
#' adjusted Poisson model to estimate fraction news while also using hierarchical modeling to partially pool information across
#' replicates and the entire dataset.
#'
#' Since \code{fast_analysis} relies on -s4U control data to estimate background mutation rates, user inputted mutation rate
#' estimates will have to be provided if such a control is not present. If the -s4U control is absent and no mutation rate
#' estimate is inputted, an error will be thrown.
#'
#' @param obj DynamicSeqData object produced by \code{DynamicSeqData} or a DynamicSeqFit object produced by \code{DynamicFit}
#' @param StanFit Logical; if TRUE, then Fully Bayesian Hierarchical model is implemented to estimate fraction news and pool information
#' across the dataset
#' @param HybridFit Logical; if TRUE, then a Fully Bayesian Hierarchical model is run using
#' estimates for the fraction new and fraction new uncertainty obtained via the efficient
#' statistical model.
#' @param keep_input Two element vector; 1st element is highest mut rate accepted in control samples, 2nd element is read count cutoff
#' @param Stan_prep Logical; if TRUE, then data_list that can be passed to Stan is curated
#' @param Fast_prep Logical; if TRUE, then dataframe that can be passed to fast_analysis() is curated
#' @param FOI Features of interest; character vector containing names of features to analyze
#' @param concat Logical; If TRUE, FOI is concatenated with output of reliableFeatures
#' @param ... Arguments passed to either \code{fast_analysis} (if a DynamicSeqData object)
#' or \code{TL_Stan} and \code{Hybrid_fit} (if a DynamicSeqFit object)
#' @return DynamicSeqFit object with results from statistical modeling and data processing
#' @export
DynamicSeqFit <- function(obj, StanFit = TRUE, HybridFit = FALSE,
                          keep_input = c(0.2, 50),
                          Stan_prep = TRUE,
                          Fast_prep = TRUE,
                          FOI = c(),
                          concat = TRUE,
                          ...){

  if(class(obj) == "DynamicSeqData"){

    data_list <- DynamicSeq::cBprocess(obj, keep_input = keep_input,
                                       Stan = Stan_prep,
                                       Fast = Fast_prep,
                                       FOI = FOI,
                                       concat = concat)

    fast_list <- DynamicSeq::fast_analysis(data_list$Fast_df, ...)


    Fit_lists <- list(Fast_Fit = fast_list,
                        Data_lists = data_list)

    class(Fit_lists) <- "DynamicSeqFit"
    return(Fit_lists)

  }else if(class(obj) == "DynamicSeqFit"){

    if(StanFit){

      Stan_list <- DynamicSeq::TL_stan(obj$Data_lists$Stan_data, ...)

      Effects <- Stan_list$Effects_df

      ## Calculate p-value using moderated t-test
      dfs <- 2*max(obj$Fast_Fit$Fn_Estimates$Replicate) - 2 + as.numeric(obj$Fast_Fit$Hyper_Parameters[1])

      Effects <- Effects %>% dplyr::mutate(padj = 2*stats::pt(-abs(effects/ses), df = dfs))

      Stan_list$Effects_df <- Effects

      rm(Effects)

      class(Stan_list) <- "HMCFit"

      obj$Stan_Fit <- Stan_list

    }
    if(HybridFit){

      Rep_Fn <- obj$Fast_Fit$Fn_Estimates

      data_list <- list(
        NE = nrow(Rep_Fn),
        NF = max(Rep_Fn$Gene_ID),
        MT = Rep_Fn$Condition,
        FE = Rep_Fn$Gene_ID,
        tl = obj$Data_lists$Stan_data$tl,
        logit_fn_rep = Rep_Fn$logit_fn,
        fn_se = Rep_Fn$logit_fn_se,
        Avg_Reads = obj$Data_lists$Stan_data$Avg_Reads,
        nMT = max(Rep_Fn$Condition),
        R = Rep_Fn$Replicate,
        nrep = max(Rep_Fn$Replicate),
        sample_lookup = obj$Data_lists$Stan_data$sample_lookup,
        sdf = obj$Data_lists$Stan_data$sdf
      )

      rm(Rep_Fn)

      Stan_list <- DynamicSeq::TL_stan(data_list, Hybrid_Fit = TRUE, ...)

      Effects <- Stan_list$Effects_df

      ## Calculate p-value using moderated t-test
      dfs <- 2*max(obj$Fast_Fit$Fn_Estimates$Replicate) - 2 + as.numeric(obj$Fast_Fit$Hyper_Parameters[1])

      Effects <- Effects %>% dplyr::mutate(padj = 2*stats::pt(-abs(effects/ses), df = dfs))

      Stan_list$Effects_df <- Effects

      rm(Effects)

      class(Stan_list) <- "HybridModelFit"

      obj$Hybrid_Fit <- Stan_list
    }

    return(obj)

  }else{
    stop("obj is not of class DynamicSeqData or DynamicSeqFit")
  }



}
