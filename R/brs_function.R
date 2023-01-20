## Functions for running BRS

#' Import BRS from Python
#'
#' For testing purposes. This function imports Python code required for running BRS
.importBRS <- function(){
  #reticulate::py_install(system.file("python", "BOAmodel.py", package="BRS"))
  reticulate::import_from_path("BOAmodel", path=system.file("python", package="BRS")) #delay=T?
  reticulate::source_python(system.file("python", "bootstrap.py", package="BRS"))
  #return(BOA)
}

#' Run BRS
#'
#' This function runs BRS. It allows the user either to run BRS only once on the
#' original data or to run BRS on bootstrapped samples
#' @param df dataframe of binary features
#' @param Y vector of binary outcome
#' @param maxLen integer maximum length of rules
#' @param trainProp numeric for proportion of data to use as training data
#' @param numIter integer for number of iterations in simulated annealing
#' @param numChain integer for number of chains of simulated annealing
#' @param numMine integer for number of rules to mine
#' @param supp integer for percent minimum support
#' @param alpha_1 numeric for alpha_+ from the paper
#' @param alpha_2 alpha_-
#' @param beta_1 beta_+
#' @param beta_2 beta_-
#' @param prior_type string for the prior type. Either "beta" for BRS-BetaBinomial or "poisson" for BRS-Poisson
#' @param alpha_l vector of alpha_l for l=1...maxLen. If set to NULL and prior_type="beta", values will be automatically generated. Ignored if prior_type="poisson"
#' @param beta_l vector of beta_l for l=1...maxLen. If set to NULL and prior_type="beta", values will be automatically generated. Ignored if prior_type="poisson"
#' @param lambda numeric rate parameter for the prior on the number of rules
#' @param nu numeric rate parameter for the prior on the length of rules
#' @param print logical whether to print progress of algorithm
#' @param bootstrap logical for whether to bootstrap
#' @param reps integer for number of bootstrap reps
#' @param sampleSize integer for bootstrap sample size. If set to NULL, default is the number of observations in data
#' @return list of rule sets from all bootstrap replications
#' @return indices of bootstrap samples
#' @return accuracy, true positive rate, and false positive rate on test data for each rule set
#' @export
BRS <- function(df, Y, maxLen, trainProp=.7,
                numIter=500L, numChain=2L, numMine=5000L, supp=5L,
                alpha_1=50L, beta_1=1L, alpha_2=50L, beta_2=1L,
                prior_type="beta",
                alpha_l=NULL, beta_l=NULL,
                lambda=NULL, nu=NULL,
                print=FALSE,
                bootstrap=FALSE, reps=1L, sampleSize=NULL,
                seed=NULL){
  # set seed
  if(!is.null(seed)){
    reticulate::py_set_seed(as.integer(seed))
  }

  # default sample size is number of observations
  if(bootstrap & is.null(sampleSize)){
    sampleSize <- nrow(df)
  }

  # load required Python code
  reticulate::import_from_path("BOAmodel", path=system.file("python", package="brs")) #delay=T?
  reticulate::source_python(system.file("python", "bootstrap.py", package="brs"))

  # run BRS
  results <- bootstrapBOA(supp, maxLen, numMine,
                as.integer(alpha_1), as.integer(beta_1),
                as.integer(alpha_2), as.integer(beta_2),
                prior_type,
                alpha_l, beta_l,
                lambda, nu,
                numIter, numChain, trainProp,
                print,
                df, Y,
                bootstrap, reps, as.integer(sampleSize))
  allRuleSets <- results[1][[1]]
  allIndices <- results[2][[1]]
  allStats <- results[3][[1]]
  return(list("Rule Sets"=allRuleSets, "Indices"=allIndices, "Stats"=allStats))
}
