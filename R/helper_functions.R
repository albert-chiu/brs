## Helper functions

#' Function that simplifies feature if possible
#' @param originalFeature feature to simplify
#' @param featureNames vector of all feature names, repeated according to multiplicity
#' @param featureValues vector of all values corresponding to features (same order as featureNames)
#' @return name of feature to use
#' Original feature if no _neg suffix
#' if _neg suffix: original feature and value if none or multiple alternative values, other feature value if exactly one alternative
simplifyFeature <- function(originalFeature, featureNames, featureValues){
  split_feat <- strsplit(originalFeature, "_")[[1]]
  name <- split_feat[1]
  value <- split_feat[2]
  # Check if _neg and if there is exactly one other possible value
  if(tail(split_feat, 1) == "neg" & (sum(featureNames==name) == 2)){
    # Return other feature and value
    return(paste(name, featureValues[featureNames==name & featureValues!=value], sep="_"))
  }
  else {
    return(originalFeature)
  }
}

#'
simplifyCondition <- function(cond, oppmat, oppind, feats){
  split <- strsplit(cond, "_")[[1]]
  if (tail(split, 1) != "neg") {  # if not negation, then return original
    return(cond)
  } else {
    fg <- split[1]
    oppmat <- matrix(oppmat[which(unlist(lapply(oppind, function(x) length(which(x==fg))>0))), ], ncol=2)
    val <- split[length(split)-1]
    ind <- as.matrix(which(oppmat==val, arr.ind = T))
    if (nrow(oppmat)>0) {
      if (ncol(ind) > 1) {
        for (i in 1:nrow(ind)) {
          row <- ind[i, 1]
          col <- ifelse(ind[i, 2]==1, yes=2, no=1)  # other column is oppmatosite value
          feat <- paste(fg, oppmat[row, col], sep="_")
          if (feat %in% feats) {  # check if this oppmatosite value is one of the possible features
            return(feat)
          }
        }
      } else {
        col <- ifelse(ind==1, yes=2, no=1)  # other column is oppmatosite value
        feat <- paste(fg, oppmat[1, col], sep="_")
        feat <- paste(fg, oppmat[row, col], sep="_")
        if (feat %in% feats) {  # check if this oppmatosite value is one of the possible features
          return(feat)
        }
      }
    }
  }
  return(cond)  # not one of the possible features
}


#' Get labels for features
#'
#' Get the labels of a features to use in graphs
#'
#' @param feat feature or vector of features to get the label of; should be of the form "feature_value"
#' @param labels_df dataframe of unique feature and corresponding labels; first column is the feature, second column is corresponding label
#' @param neg_label prefix to use for negative features
#' @return label to use for feat
getLabel <- function(feat, labels_df, neg_label){
  label <- c()
  for(i in 1:length(feat)){
    neg <- ""
    split_feat <- strsplit(feat[i], "_")[[1]]
    if ( tail(split_feat, 1) == "neg" ){
      neg <- neg_label
      feat[i] <- paste(split_feat[1:(length(split_feat)-1)], collapse="_")
    }
    label[i] <- paste(neg, labels_df[labels_df[,1]==feat[i],2], sep="")
  }
  return(label)
}



#' Get the number of cases that satisfy a rule
#'
#' @param rule rule to check: must be in the format feature_value or feature_value_neg
#' @param df data to check
#' @param Y if included, checks that case's Y=y_val
#' @param y_val value Y should take
#' @return number of cases in data that satisfies rule
numCases <- function(rule, df, Y=NULL, y_val=1){
  feats <- c() # Features to check
  values <- c() # Values of features
  for(i in 1:length(rule)){ # Loops through all features in rule to get feature (w/o _neg) and value
    feats[i] <- rule[i]
    values[i] <- 1
    split_feat <- strsplit(feats[i], "_")[[1]]

    if ( tail(split_feat, 1) == "neg" ) {
    #if(length(split_feat) == 3){ # Checks if _neg
      values[i] <- 0
      feats[i] <- paste(split_feat[1:(length(split_feat)-1)], collapse="_")
    }
  }
  if(is.null(Y)){
    val <- sum(apply(df[,feats] == matrix(rep(values, nrow(df)), ncol=length(values), byrow=T), 1, all)) # Number of cases that satisfy rule
  } else { # must satisfy Y
    val <- sum(apply(cbind(Y, df[,feats]) == matrix(rep(c(y_val,values), nrow(df)), ncol=(length(values)+1), byrow=T), 1, all)) # Number of cases that satisfy rule and Y
  }
  return(val)
}

#'
getFeat <- function(rule) {
  split_feat <- strsplit(rule, "_")[[1]]
  feat <- ifelse(tail(split_feat, 1) == "neg",
         yes=paste(split_feat[1:(length(split_feat)-1)], collapse="_"),
         no=paste(split_feat, collapse="_"))
  return(feat)
}

get_Yhat <- function(ruleset, df) {
  yhat <- rep(F, times=nrow(df))
  for(rule in ruleset) {
    split <- strsplit(rule, "_")
    vals <- unlist(lapply(split, function(x) x[length(x)]!="neg"))
    feats <- sapply(1:length(split), function(i) paste(split[[i]][1:(length(split[[i]])-!vals[i])], collapse="_"))
    yhat <- (yhat | apply(df[,feats] == matrix(rep(vals, nrow(df)), ncol=length(vals), byrow=T), MARGIN=1, all))
  }
  return(yhat)
}

get_stat <- function(Yhat, Y, stat="acc") {
  if (stat=="acc") {
    return(sum(Yhat==Y)/length(Y))
  }
}


get_topRules <- function(X, Y, fit, maxLen, topRules=5, minProp=0, simplify=F, oppmat=NULL, oppind=NULL) {
  # names of features (repeated according to multiplicity), and their values
  # needed for simplifyFeature (changes _neg to positive if only two possible values)
  allRuleSets <- fit[["Rule Sets"]]
  allFeatures <- colnames(X)

  feat_split = strsplit(allFeatures, "_")
  feat_names <- c()
  feat_values <- c()
  for(i in 1:length(feat_split)){
    feat_names <- append(feat_names, feat_split[[i]][1])
    feat_values <- append(feat_values, feat_split[[i]][2])
  }
  feat_values[is.na(feat_values)] <- 1

  # find most frequent rules of each length
  byLen <- vector(mode="list", length=maxLen) # rules grouped by length
  for(i in 1:length(allRuleSets)){ # loop through all rule sets
    ruleSet <- allRuleSets[[i]]
    for(j in 1:length(ruleSet)){ # loop through all rules in rule set
      rule <- sort(ruleSet[[j]]) # sort alphabetically
      if(simplify){
        for(k in 1:length(rule)){ # loop through each condition in the rule
          rule[k] <- simplifyCondition(rule[k], oppmat=oppmat, oppind=oppind, allFeatures) # simplify if possible
        }
      }
      len <- length(rule)
      byLen[[len]] <- rbind(byLen[[len]], rule) # append to group of rules of same length
    }
  }

  # get top topRules rules of each length
  rules <- vector(mode="list", length=maxLen) ## names of the most frequent rules
  freq <- vector(mode="list", length=maxLen) ## frequencies of the most frequent rules
  reps <- length(allRuleSets)
  for(len in 1:maxLen){
    if (!is.null(byLen[[len]])) {
      temp <- as.data.frame(byLen[[len]]) %>% dplyr::group_by_all() %>% dplyr::tally(sort=T) ## group together all columns and count number of times each row (rule) appears
      thisFreq <- c(temp[,len+1]/reps)[[1]] # frequency as a proporiton of number of bootstraps
      keep <- (thisFreq >= minProp)  # indices of rules that appear sufficiently many time
      if ( sum(keep) > 0 ) {
        rules[[len]] <- as.data.frame(temp[keep,1:len][1:min(topRules, sum(keep)),])
        freq[[len]] <- thisFreq[keep][1:min(topRules, sum(keep))]
      }
      colnames(rules[[len]]) <- NULL
    }
  }

  return(list("rules"=rules, "freq"=freq))
}
