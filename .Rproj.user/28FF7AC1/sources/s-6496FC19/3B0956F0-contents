
#' Find aggregate rule set
#'
#' Find the aggregate rule set from a list of bootstrapped BRS rule sets
#'
#' @param fit the output from the BRS function. A list whose first element is a list of rule sets and whose second element is a list of bootstrap indices. The third element is ignored.
#' @param X data frame or matrix of the data, excluding the outcome
#' @param Y vector of outcomes
#' @param maxLen maximum length of a rule possible
#' @param split logical for whether to split the sample into a training set on which the aggregate rule set is found and a test set on which that rule set's performance is evaluated
#' @param train numeric for proportion of the data to use as training data. If split=F, this argument is ignored.
#' @param maxRules integer for the maximum number of rules in the aggregate rule set
#' @param stat the statistic on which to evaluate the aggregated rule sets. Currently only accuracy is supported
#' @param topRules integer for the number of high prevalence rules of each length to consider
#' @param minProp numeric for proportion of times a rule must appear in order to be considered
#' @param simplify logical for whether equivalent rules are combined for determining prevalence
#' @param oppmat a matrix with two columns and K rows, where K is the length of
#'        the list oppind. The kth row contains values v1 and v2 (i.e.,
#'        v1=oppmat[k,1] and v2=oppmat[k,2]) such that for any variable var in
#'        oppind[[k]], var_v1 and !var_v2 are equivalent. v1 should be the
#'        prefered return value.
#' @param oppind a list of vectors of variables. Each vector oppind[[k]] contains
#'        variables var such that var_v1 and !var_v2 are equivalent, where v1
#'        and v2 form the kth row of oppmat, v1=oppmat[k,1] and v2=oppmat[k,2]
#' @return the aggregate rule set, which has the highest stat out of all possible rule sets constructed from at most maxRules candidate rules
#' @export
agg_BRS <- function(fit, X, Y, maxLen, split= F, train=.7, maxRules=3, stat="acc", topRules=5, minProp=0, simplify=F, oppmat=NULL, oppind=NULL) {
  if (!split) {  # use whole sample to get which one to use
    X_test <- X
    Y_test <- Y
  } else {
    if (length(train)==1) { # user specifies train proportion instead of indices
      train <- sample(1:nrow(X), size=train*nrow(X))
    }
    X_test <- X[-train, ]
    Y_test <- Y[-train]
    X <- X[train, ]
    Y <- Y[train]
  }

  # out of bootstrapped-sample
  test_ind <- fit[["Indices"]]

  rules <- get_topRules(X=X, Y=Y, fit=fit, maxLen = maxLen,
                        topRules=topRules, minProp=minProp,
                        simplify=simplify, oppmat=oppmat, oppind=oppind)[["rules"]]

  # indices of rules; first column is index in list, second column is row number in rules[[x]]
  ind <- do.call(rbind, lapply(c(1:length(rules))[!sapply(rules, is.null)], function(x) cbind(x, 1:nrow(rules[[x]]))))

  combn <- lapply(1:(min(maxRules, nrow(ind))), function(x) gtools::combinations(n=nrow(ind), r=x))


  bestRSs <- list()
  bestStats <- list()
  for (i in 1:length(combn)) {
    maxStat <- 0
    for (j in 1:nrow(combn[[i]])) {
      this <- as.data.frame(matrix(ind[combn[[i]][j,],], nrow=i))
      rs <- lapply(1:nrow(this), function(x) unlist(rules[[this[x, 1]]][this[x, 2], ]))
      #this_stat <- get_stat(Yhat = get_Yhat(rs, df=X), Y = Y, stat = stat)
      this_stat <- mean(sapply(test_ind, function(x) get_stat(Yhat = get_Yhat(rs, df=X[-x, ]), Y = Y[-x], stat = stat))) # number of obs: length((1:nrow(X))[-x]))
      if (this_stat > maxStat) {
        bestRSs[[i]] <- rs
        maxStat <- this_stat
      }
    }
    bestStats[[i]] <- get_stat(Yhat = get_Yhat(bestRSs[[i]], df=X_test), Y = Y_test, stat = stat)
  }

  #return(list("Rule Sets"=bestRSs, "Performance"=bestStats, train))

  # only return best rule set
  return(bestRSs[[which.max(unlist(bestStats))]])
}
