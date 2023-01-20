## Functions for making bar plot

#' Make bar plot
#'
#' This function makes a bar plot of rule frequency and coverage
#' @importFrom magrittr %>%
#' @param df dataframe of binary features
#' @param Y vector of binary outcome
#' @param fit output from function \code{BRS}
#' @param featureLabels data frame whose first column is feature names (as they appear in df) and whose second column is the corresponding labels to be displayed on the graph
#' @param maxLen integer maximum length of rules
#' @param topRules the max number of rules to plot for each length
#' @param and string for how an "and" statement will be concatenated
#' @param neg string for how a "not" statement will be prefaced
#' @param minProp minimum proportion of rule sets in which a rule must appear before being plotted
#' @param simplify whether to simplify negative rules when possible. If true, then a negative rule is changed to it's positive counterpart if and only if there are only two possible values for the corresponding variable (e.g. "A_1_neg" is changed to "A_0" if and only if the only two possible values for "A" are 0 and 1)
#' @param heightBuffer amount by which to shift bottom of graph up to make room for label
#' @param textSize \code{size} parameter passed to \code{ggplot2::element_text}
#' @return bar plot with frequency and coverage of rules
#' @export
.plot_bar <- function(df, Y, fit, featureLabels,
                     maxLen, topRules=10,
                     and=" AND ", neg="NOT ",
                     minProp=0, simplify=T, oppmat=NULL, oppind=NULL,
                     heightBuffer=1, plotBuffer=0, textSize=16,
                     boot_rep=100){

  allRuleSets <- fit[["Rule Sets"]]
  allIndices <- fit[["Indices"]]
  if ( length(allIndices) == 0 ) {
    allIndices <- lapply(1:boot_rep, function(x) sample(nrow(df), nrow(df), replace=T))
  }
  allFeatures <- colnames(df)

  # names of features (repeated according to multiplicity), and their values
  # needed for simplifyFeature (changes _neg to positive if only two possible values)
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
    }
  }

  # calculate coverage statistics
  tp <- vector(mode="list", length=maxLen)
  fp <- vector(mode="list", length=maxLen)
  for(len in 1:maxLen){
    if(!is.null(rules[[len]]) && nrow(rules[[len]])>0){
      tp[[len]] <- .getTP(rules=rules[[len]], allIndices=allIndices, reps=reps, df=df, Y=Y)
      fp[[len]] <- .getFP(rules=rules[[len]], allIndices=allIndices, reps=reps, df=df, Y=Y)
    }
  }

  # limits for plot
  maxStat <- max(unlist(tp))
  minStat <- -max(unlist(fp))
  maxFreq <- max(unlist(freq))

  # format data for plotting
  p_data_freq <- vector(mode="list", length=maxLen)
  p_data_stats <- vector(mode="list", length=maxLen)
  for(len in 1:maxLen){
    if(!is.null(rules[[len]])){
      rules[[len]][] <- lapply(rules[[len]], function(x) getLabel(x, featureLabels, neg_label=neg)) # convert rules into readable format
      rules_cat <- apply(rules[[len]], 1, function(x) paste(x, collapse=and)) # concat rules into single "and" phrase
      rules_cat <- factor(rules_cat, rules_cat)
      p_data_freq[[len]] <- data.frame("rules"=rules_cat, "freq"=freq[[len]])
      p_data_stats[[len]] <- data.frame("rules"=rules_cat, tp[[len]])
    }
  }

  # make plots
  p <- list()
  makeLabel <- TRUE # add X axis label to last plot only
  heights <- unlist(lapply(freq, function(x) length(x)))
  for(len in maxLen:1){
    if ( !is.null(rules[[len]]) ){
      p_freq <- ggplot2::ggplot(p_data_freq[[len]])+
        ggplot2::geom_bar(ggplot2::aes(x=reorder(rules, freq), y=freq), stat="identity", width=.7, alpha=.5)+
        ggplot2::lims(y=c(maxFreq, 0))+
        ggplot2::coord_flip()+
        ggplot2::theme_bw()+
        ggplot2::xlab(paste0("Length ", len))
      p_stats <- ggplot2::ggplot(p_data_stats[[len]])+
        # tp
        ggplot2::geom_bar(ggplot2::aes(x=reorder(rules, p_data_freq[[len]]$freq), y=median), stat="identity", width=.7, alpha=.5)+
        ggplot2::geom_errorbar(ggplot2::aes(x=reorder(rules, p_data_freq[[len]]$freq), min=min, max=max), width=.5, alpha=.7, size=.25)+
        # fp
        ggplot2::geom_bar(ggplot2::aes(x=reorder(rules, p_data_freq[[len]]$freq), y=-fp[[len]]$median), fill="red", stat="identity", width=.7, alpha=.5)+
        ggplot2::geom_errorbar(ggplot2::aes(x=reorder(rules, p_data_freq[[len]]$freq), min=-fp[[len]]$min, max=-fp[[len]]$max), color="red", width=.5, alpha=.7, size=.25)+
        ggplot2::coord_flip()+
        ggplot2::theme_bw()+
        ggplot2::lims(y=c(minStat,maxStat))

      # add X-label if last plot; otherwise, make blank
      if(makeLabel){
        p_freq <- p_freq+
          ggplot2::theme(text=ggplot2::element_text(size=textSize))+
          ggplot2::ylab("Prevalence")
        p_stats <- p_stats+
          ggplot2::theme(text=ggplot2::element_text(size=textSize), axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank())+
          ggplot2::ylab("Coverage")
        heights[len] <- heights[len] + heightBuffer # account for the x-axis labels
        makeLabel <- FALSE
      } else {
        p_freq <- p_freq+
          ggplot2::theme(text=ggplot2::element_text(size=textSize),
                         axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank())
        p_stats <- p_stats+
          ggplot2::theme(text=ggplot2::element_text(size=textSize),
                         axis.title.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank(),
                         axis.title.y=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks.y=ggplot2::element_blank())
      }

      p[[len]] <- cbind(ggplot2::ggplotGrob(p_freq), ggplot2::ggplotGrob(p_stats))
    }
  }

  heights <- (heights+plotBuffer)/sum(heights)
  return(cowplot::plot_grid(plotlist=p,align = "v", nrow = maxLen, rel_heights = heights))
}


#' Make bar plot
#'
#' Make a bar plot of prevalence and coverage for rules produced by BRS w/
#' bootstrapping. Separate functions are used for voting and democracy because
#' the bars that are darkened (those corresponding to the aggregated rule set)
#' differ and are set manually within the function.
#'
#' @param df data, without outcome
#' @param Y outcome
#' @param fit the output of the BRS function with bootstrap=T, which is a list
#'            whose first entry is a list of rule sets and whose second entry is
#'            the indices of each bootstrap sample (the  list's third element is
#'            not used).
#' @param featureLabels the labels_df input for the getLabels function
#' @param maxLen the maximum length of a rule allowed
#' @param topRules the maximum number of rules of each length to plot
#' @param and a string for how the logical 'and' operator will appear in labels
#' @param neg a string for how the logical 'not' operator will appear in labels
#' @param minProp the minimum proportion of times a rule has to appear to be
#'                included
#' @param simplify a logical for whether rules from the same equivalence class
#'                 will be treated as the same (e.g., 'not low wealth'
#'simplifies to 'medium or high wealth')
#' @param oppmat the oppmat input for the simplifyCondition function. Only
#'               required if simplify=T
#' @param oppind the oppind input for the simplifyCondition function. Only
#'               required if simplify=T
#' @param heightBuffer a graphical parameter. Larger values make plots for
#'                     different rule lengths closer when there are unequal
#'                     numbers of rules of each length
#' @param plotBuffer a graphical parameter to create white space
#' @param titleSize a graphical parameter to specify text size for titles
#' @param rule_text_size a graphical parameter to specify text size for rule
#'                       name labels
#' @param number_size a graphical parameter for text size of numbers
#' @param boot_rep the number of bootstraps
#' @param neg_color color for negative coverage
#' @param darken which features to darken
#' @return a bar plot of the most prevalent rules of each length, sorted by
#'         prevalence, and their prevalence and (negative and positive)
#'         Darkens the bars for the rules chosen by our aggregation
#'         method for the democracy example.
### allows darkening; to be fixed ###
plot_bar <- function(df, Y, fit, featureLabels, maxLen, topRules=10,
                     and=" AND ", neg="NOT ", minProp=0, simplify=T,
                     oppmat=NULL, oppind=NULL, heightBuffer=1, plotBuffer=0,
                     titleSize=16, rule_text_size=16, number_size=16,
                     boot_rep=100, neg_color = "red",
                     darken = NULL){

    neg_alpha <- 1 # .5

    allRuleSets <- fit[["Rule Sets"]]
    allIndices <- fit[["Indices"]]
    if ( length(allIndices) == 0 ) {
      allIndices <- lapply(1:boot_rep, function(x) sample(nrow(df), nrow(df), replace=T))
    }
    allFeatures <- colnames(df)

    # names of features (repeated according to multiplicity), and their values
    # needed for simplifyFeature (changes _neg to positive if only two possible values)
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
      }
    }

    # calculate coverage statistics
    tp <- vector(mode="list", length=maxLen)
    fp <- vector(mode="list", length=maxLen)
    for(len in 1:maxLen){
      if(!is.null(rules[[len]]) && nrow(rules[[len]])>0){
        tp[[len]] <- .getTP(rules=rules[[len]], allIndices=allIndices, reps=reps, df=df, Y=Y)
        fp[[len]] <- .getFP(rules=rules[[len]], allIndices=allIndices, reps=reps, df=df, Y=Y)
      }
    }

    # limits for plot
    maxStat <- max(unlist(tp))
    minStat <- -max(unlist(fp))
    maxFreq <- max(unlist(freq))

    # format data for plotting
    p_data_freq <- vector(mode="list", length=maxLen)
    p_data_stats <- vector(mode="list", length=maxLen)
    for(len in 1:maxLen){
      if(!is.null(rules[[len]])){
        rules[[len]][] <- lapply(rules[[len]], function(x) getLabel(x, featureLabels, neg_label=neg)) # convert rules into readable format
        rules_cat <- apply(rules[[len]], 1, function(x) paste(x, collapse=and)) # concat rules into single "and" phrase
        rules_cat <- factor(rules_cat, rules_cat)
        p_data_freq[[len]] <- data.frame("rules"=rules_cat, "freq"=freq[[len]])
        p_data_stats[[len]] <- data.frame("rules"=rules_cat, tp[[len]])
      }
    }

    if (is.null(darken)) {
      darken_none <- T
    }

    # make plots
    p <- list()
    makeLabel <- TRUE # add X axis label to last plot only
    heights <- unlist(lapply(freq, function(x) length(x)))
    for(len in maxLen:1){

      if (darken_none == T) {
        darken <- rep(FALSE, times=nrow(p_data_freq[[len]]))
      } else {
# to add code to automatically darken aggregated rules
      }

      if ( !is.null(rules[[len]]) ){
        p_freq <- ggplot2::ggplot(p_data_freq[[len]])+
          ggplot2::geom_bar(ggplot2::aes(x=reorder(rules, freq), y=freq), stat="identity", width=.7, alpha=ifelse(darken, yes=1, no=.5))+
          ggplot2::lims(y=c(maxFreq, 0))+
          ggplot2::coord_flip()+
          ggplot2::theme_bw()+
          ggplot2::xlab(paste0("Length ", len))+
          ggplot2::theme(axis.text.y=ggplot2::element_text(size=rule_text_size))+
          ggplot2::scale_fill_grey()
        p_stats <- ggplot2::ggplot(p_data_stats[[len]])+
          # tp
          ggplot2::geom_bar(ggplot2::aes(x=reorder(rules, p_data_freq[[len]]$freq), y=median), stat="identity", width=.7, alpha=ifelse(darken, yes=1, no=.5))+
          ggplot2::geom_errorbar(ggplot2::aes(x=reorder(rules, p_data_freq[[len]]$freq), min=min, max=max), width=.5, alpha=ifelse(darken, yes=1, no=.7), size=.5)+
          # fp
          ggpattern::geom_bar_pattern(ggplot2::aes(x=reorder(rules, p_data_freq[[len]]$freq), y=-fp[[len]]$median),
                                      stat="identity",
                                      width=.7,
                                      fill = 'white',
                                      colour = ifelse(darken, yes=neg_color, no=prettyGraphs::add.alpha(neg_color, alpha=neg_alpha)),
                                      pattern_color = neg_color,
                                      pattern_alpha = ifelse(darken, yes=1, no=neg_alpha),
                                      pattern_fill = neg_color,
                                      pattern = 'stripe',
                                      pattern_density = .5,
                                      pattern_spacing = 0.1) +
          ggplot2::geom_errorbar(ggplot2::aes(x=reorder(rules, p_data_freq[[len]]$freq), min=-fp[[len]]$min, max=-fp[[len]]$max), color=neg_color, width=.5, alpha=ifelse(darken, yes=1, no=.7), size=.5)+
          ggplot2::coord_flip()+
          ggplot2::theme_bw()+
          ggplot2::lims(y=c(minStat,maxStat))+
          ggplot2::scale_fill_grey()

        # add X-label if last plot; otherwise, make blank
        if(makeLabel){
          p_freq <- p_freq+
            ggplot2::theme(text=ggplot2::element_text(size=titleSize),
                           axis.text.x=ggplot2::element_text(size=number_size))+
            ggplot2::ylab("Prevalence")
          p_stats <- p_stats+
            ggplot2::theme(text=ggplot2::element_text(size=titleSize),
                           axis.title.y=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks.y=ggplot2::element_blank())+
            ggplot2::ylab("Coverage")
          heights[len] <- heights[len] + heightBuffer # account for the x-axis labels
          makeLabel <- FALSE
        } else {
          p_freq <- p_freq+
            ggplot2::theme(text=ggplot2::element_text(size=titleSize),
                           axis.title.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank())
          p_stats <- p_stats+
            ggplot2::theme(text=ggplot2::element_text(size=titleSize),
                           axis.title.x=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.ticks.x=ggplot2::element_blank(),
                           axis.title.y=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks.y=ggplot2::element_blank())
        }

        p[[len]] <- cbind(ggplot2::ggplotGrob(p_freq), ggplot2::ggplotGrob(p_stats))
      }
    }

    heights <- (heights+plotBuffer)/sum(heights)

    return(cowplot::plot_grid(plotlist=p,align = "v", nrow = maxLen, rel_heights = heights))
}


#' Get proportion of true positives
.getTP <- function(rules, allIndices, reps, df, Y){
  stats <- c()
  for(i in 1:nrow(rules)){ # Loop through rules
    coverage <- c()
    for(j in 1:reps){ # Loop through bootstrap samples
      ind <- allIndices[[j]]
      coverage[j] <- numCases(rule=unlist(rules[i,]), df=df[ind,], Y=Y[ind]) # number of true positives in that bootstrap
    }
    stats <- rbind(stats,
                   c(quantile(coverage, .025), median(coverage), quantile(coverage, .975))/length(Y))
  }
  colnames(stats) <- c("min", "median", "max") # min=.025 quantile, max=.975 quantile
  return(data.frame(stats))
}

#' Get proportion of false positives
.getFP <- function(rules, allIndices, reps, df, Y){
  stats <- c()
  for(i in 1:nrow(rules)){ # Loop through rules
    coverage <- c()
    for(j in 1:reps){ # Loop through bootstrap samples
      ind <- allIndices[[j]]
      coverage[j] <- numCases(rule=unlist(rules[i,]), df=df[ind,], Y=Y[ind], y_val=0) # number of false positives in that bootstrap
    }
    stats <- rbind(stats,
                   c(quantile(coverage, .025), median(coverage), quantile(coverage, .975))/length(Y))
  }
  colnames(stats) <- c("min", "median", "max") # min=.025 quantile, max=.975 quantile
  return(data.frame(stats))
}

