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
#' \code{fast_analysis} will not be called again unless \code{Fast} is set to \code{TRUE}. If you want to
#' generate model fits using different parameters for cBprocess, you will have to rerun \code{DynamicFit}
#' on the DynamicSeqData object.
#'
#' See the doumentation for the indvidual fitting functions for details regarding how they analyze nucleotide
#' recoding data. What follows is a brief overview of how each works
#'
#' \code{fast_analysis} either estimates mutation rates from + and - s4U samples or uses mutation rate estimates
#' provided by the user to perform maximum likelihood estimation (MLE) of the fraction of RNA that is labeled for each
#' replicate of nucleotide recoding data provided. Uncertainties for each replicate's estimate are approximated using
#' asymptotic results involving the Fisher Information. Replicate data is pooled using an approximation to hierarchical
#' modeling that relies on analytic solutions to simple Bayesian models. Linear regression is used to estimate the
#' relationship between read depths and replicate variabiltiy for uncertainty estimation regularization, again performed
#' using analytic solutions to models.
#'
#' \code{Hybrid_fit} takes as input estimates of the logit(fraction new) and uncertainty provided by \code{fast_analysis}.
#' It then uses Stan on the backend to implement a hierarchical modeling that pools data across replicates and the dataset
#' to estimate effect sizes (L2FC(kdeg)) and uncertainties. Replicate variability information is pooled across each experimental
#' condition to regularize variance estimates using a hierarchical linear regression model.
#'
#' \code{TL_Stan} uses Stan on the back end to impelement a model similar to \code{Hybrid_fit}. The difference is that while
#' \code{Hybrid_fit} requires fraction new estimates from \code{fast_analysis}, \code{TL_Stan} implements a U-content exposure
#' adjusted Poisson model to estimate fraction news while also using hierarchical modeling to partially pool information across
#' replicates and the entire dataset.
#'
#' Since \code{fast_analysis} relies on -s4U control data to estimate background mutation rates, user inputted mutation rate
#' estimates will have to be provided if such a control is not present. If the -s4U control is absent and no mutation rate
#' estimate is inputted, a non-hierarchical (if Stan = FALSE) or the fully hierarchical (if Stan = TRUE, implemented by \code{TL_Stan})
#' model will first be run on the data to estimate mutation rates, and these mutation rates as well as the fraction news and uncertainties
#' estimated in each replicate by these models will be passed to \code{fast_analysis}.
#'
#' @param obj DynamicSeqData object produced by \code{DynamicSeqData} or a DynamicSeqFit object produced by \code{DynamicFit}
#' @param StanFit Logical; if TRUE, then Fully Bayesian Hierarchical model is implemented
#' @param HybridFit Logical; if TRUE, then a Fully Bayesian Hierarchical model is run using
#' estimates for the fraction new and fraction new uncertainty obtained via the efficient
#' statistical model.
#' @param FastRefit Logical; if TRUE, then the efficient statistical model that does not
#' use Stan in the backend will be reran even if results from a previous run are present
#' in the DynamicSeqData object.
#' @param ... Arguments passed to \code{cBprocess}, \code{fast_analysis}, \code{TL_stan}, or \code{Hybrid_fit}.
#' See the documentation for these functions to learn what arguments are possible and their
#' default values
#' @return DynamicSeqFit object with results from statistical modeling and data processing
DynamicSeqFit <- function(obj, StanFit, HybridFit, FastRefit, ...){

  if(class(obj) == "DynamicSeqData"){
    data_list <- DynamicSeq::cBprocess(obj, ...)

    fast_list <- DynamicSeq::fast_analysis(data_list$Fast_df, ...)

    if(StanFit){

      Stan_list <- DynamicSeq::TL_stan(data_list$Stan_data)

    }

    if(HybridFit){

      Hybrid_list <- DynamicSeq::Hybrid_Fit(data_list, fast_list)

    }

    if(all(c(StanFit, HybridFit))){
      Fit_lists <- list(Fast_Fit = fast_list,
                        Stan_Fit = Stan_list,
                        Hybrid_Fit = Hybrid_list,
                        Data_lists = data_list)
    }else if(StanFit){
      Fit_lists <- list(Fast_Fit = fast_list,
                        Stan_Fit = Stan_list,
                        Data_lists = data_list)
    }else if (HybridFit){
      Fit_lists <- list(Fast_Fit = fast_list,
                        Hybrid_Fit = Hybrid_list,
                        Data_lists = data_list)
    }else{
      Fit_lists <- list(Fast_Fit = fast_list,
                        Data_lists = data_list)
    }


  }else if(class(obj) == "DynamicSeqFit"){


  }else{
    stop("obj is not of class DynamicSeqData or DynamicSeqFit")
  }


}
