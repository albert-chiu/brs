## Functions for making chord diagrams

#' Make Chord Diagram
#'
#' Makes a chord diagram representing the frequency of interactions between features
#' @param fit output from function \code{BRS}
#' @param featureNames names of features as they appeared in df
#' @param featureGroups dataframe whose first column is the features (only the prefixes, not the values) in \code{featureNames} and whose second column is the corresponding labels to display
#' @param linkColors vector of colors to use for chords
#' @param gridColors vector of colors to use for arcs
#' @param minProp minimum proportion of times a relationship must appear to be plotted
#' @param textSize \code{cex} parameter to be passed to \code{par}
#' @return chord diagram
#' @export
plot_chord <- function(fit, featureNames, featureGroups,
                       linkColors, gridColors, bgLinkColor="gray",
                       maxLen, lengths=c(1:maxLen), dims=c(1, maxLen),
                       minProp=0, textSize=1, line_arg=1, side_mar=0, top_mar=0){

  allRuleSets <- fit[["Rule Sets"]]

  df <- .get_df_chord(allRuleSets, featureGroups, maxLen=maxLen, minProp=minProp)

  print(df)
  # plot
  circlize::circos.clear()
  par(mfrow=dims, mar = c(0, side_mar, top_mar, side_mar), cex=textSize)
  for (len in lengths) {
    colors <- rep(bgLinkColor, times=nrow(df))
    ind <- which(df[, "deg"]==len)
    if (length(ind) > 0) {
      colors[ind] <- linkColors[as.numeric(df[ind, "color_ind"])]
    }
    circlize::chordDiagram(df[, 1:3],
                           link.sort=T,
                           grid.col = gridColors, col = colors,
                           annotationTrack = c("name", "grid"),
                           annotationTrackHeight = c(0.03, 0.05),
                           self.link = 1)
    title(paste0("Length ", len), line=line_arg)
  }
  circlize::circos.clear()

  # adjacency matrix -----------------
  if(F){
  amat <- .amat_old(allRuleSets, featureNames, featureGroups)
  amat[amat<minProp*length(allRuleSets)] <- 0

  circlize::circos.clear()
  par(cex = textSize, mar = c(0, 0, 0, 0), mfrow=c(1,1))
  circlize::chordDiagram(amat,
                         grid.col = gridColors, row.col = linkColors,
                         annotationTrack = c("name", "grid"),
                         self.link = 1)
  circlize::circos.clear()
  }

  if(F) {
  # get adjacency matrices
  mats <- .amat(allRuleSets, featureGroups, maxLen=maxLen, minProp=minProp)
  amats <- mats[[1]]
  ind_mats <- mats[[2]]

  amat <- amats[[1]]
  for (i in 2:maxLen) {
    amat <- amat + amats[[i]]
  }

  # plot
  #circlize::circos.par(gap.after = 1/rep(1,nrow(amat_group)))
  par(cex = textSize, mar = c(0, 0, 0, 0), mfrow=dims)
  for (i in 1:maxLen) {
    colorMat <- matrix(bgLinkColor, nrow=nrow(amat), ncol=ncol(amat))
    ind <- which(amats[[i]]!=0)
    if (length(ind) > 0) {
      colorMat[ind] <- linkColors[ind_mats[[i]][ind]]
    }
    circlize::chordDiagram(amat,
                           grid.col = gridColors, col = colorMat,
                           annotationTrack = c("name", "grid"),
                           self.link = 1)
  }
  circlize::circos.clear()


  circlize::chordDiagram(amat,
                         grid.col = gridColors, row.col = linkColors,
                         annotationTrack = c("name", "grid"),
                         self.link = 1)
  }
}


#' Make data frame for circilize
#'
#' @return data frame
.get_df_chord <- function(allRuleSets, featureGroups, maxLen, minProp){
  # create separate lists of rules for each length
  allRules <- vector(mode="list", length=maxLen)
  for (ruleSet in allRuleSets) {
    for (rule in ruleSet) {
      allRules[[length(rule)]] <- c(allRules[[length(rule)]],
                                    paste(sort(sapply(rule, function(x)
                                      featureGroups[featureGroups[,1]==.getFeatGroup(x), 2][1])), collapse="__"))
    }
  }

  # tabulate interactions
  counts <- vector(mode="list", length=maxLen)
  for (len in 1:maxLen) {
    #counts[[len]] <- as.data.frame( table(allRules[[len]]) )
    #counts[[len]] <- counts[[len]][counts[[len]]["Freq"] >= minProp*length(allRuleSets), ]
    counts[[len]] <- table(allRules[[len]])
    counts[[len]] <- counts[[len]][counts[[len]] >= minProp*length(allRuleSets)]
  }

  # make dfs
  df <- c()

  if (length(counts[[1]]) > 0) {
    for (i in 1:length(counts[[1]])) {
      count <- counts[[1]][i]
      df <- rbind.data.frame(df, c(names(count), names(count), count, 1, i))
    }
  }

  # interactions
  for (len in 2:maxLen) {
    if (length(counts[[len]]) > 0) {
      for (i in 1:length(counts[[len]])) {
        count <- counts[[len]][i]
        feats <- strsplit(names(count), "__")[[1]]
        for (j in 2:length(feats)) {
          df <- rbind.data.frame(df, c(feats[1], feats[j], count, len, max(0, as.numeric(df[df[, 4]==(len-1), 5]))+i))
        }
      }
    }
  }

  colnames(df) <- c("from", "to", "freq", "deg", "color_ind")

  return(df)
}


#' Make an Adjacency Matrix
#'
#' @return adjacency matrix
.amat <- function(allRuleSets, featureGroups, maxLen, minProp){
  # create separate lists of rules for each length
  allRules <- vector(mode="list", length=maxLen)
  for (ruleSet in allRuleSets) {
    for (rule in ruleSet) {
      allRules[[length(rule)]] <- c(allRules[[length(rule)]],
                                    paste(sort(sapply(rule, .getFeatGroup)), collapse="__"))
    }
  }

  # tabulate interactions
  counts <- vector(mode="list", length=maxLen)
  for (len in 1:maxLen) {
    #counts[[len]] <- as.data.frame( table(allRules[[len]]) )
    #counts[[len]] <- counts[[len]][counts[[len]]["Freq"] >= minProp*length(allRuleSets), ]
    counts[[len]] <- table(allRules[[len]])
    counts[[len]] <- counts[[len]][counts[[len]] >= minProp*length(allRuleSets)]
  }

  # make adjacency matrices
  allFeats <- unique(featureGroups[,1])
  empty_amat <- matrix(0, nrow=length(allFeats), ncol=length(allFeats)) # Empty adjacency matrix
  rownames(empty_amat) <- allFeats; colnames(empty_amat) <- allFeats

  amats <- lapply(1:maxLen, function(x) empty_amat)
  ind_mats <- amats

  # alone
  if (length(counts[[1]]) > 0) {
    for (i in 1:length(counts[[1]])) {
      count <- counts[[1]][i]
      amats[[1]][names(count), names(count)] <- count
      ind_mats[[1]][names(count), names(count)] <- i
    }
  }

  # interactions
  for (len in 2:maxLen) {
    if (length(counts[[len]]) > 0) {
      for (i in 1:length(counts[[len]])) {
        count <- counts[[len]][i]
        feats <- strsplit(names(count), "__")[[1]]
        for (j in 2:length(feats)) {
          amats[[len]][feats[1], feats[j]] <- count
          ind_mats[[len]][feats[1], feats[j]] <- max(ind_mats[[len-1]])+i  # keep track of which interaction it's a part of
        }
      }
    }
  }

  return(list(amats, ind_mats))
}

#'
.getFeatGroup <- function(rule) {
  return(strsplit(rule, "_")[[1]][1])
}



# old amat function
.amat_old <- function(allRuleSets, featureNames, featureGroups){
  feat_split = strsplit(featureNames, "_")
  feat_names <- c()
  feat_values <- c()
  for(i in 1:length(feat_split)){
    feat_names <- append(feat_names, feat_split[[i]][1])
    feat_values <- append(feat_values, feat_split[[i]][2])
  }

  feats_appearances <- NULL
  allFeatGroups <- NULL
  for(i in 1:length(allRuleSets)){ # for ith rule set
    ruleSet <- allRuleSets[[i]]
    for(j in 1:length(ruleSet)){# for jth rule in ith rule set
      rule <- ruleSet[[j]]
      for(k in 1:length(rule)){ # extract each feature
        feats_appearances <- append(feats_appearances, simplifyFeature(rule[k], feat_names, feat_values))
        allFeatGroups <- append(allFeatGroups, strsplit(rule[k], "_")[[1]][1])
      }
    }
  }
  freq_featGroup <- sort(table(allFeatGroups), decreasing = TRUE)
  amat_group <- matrix(0, nrow = length(names(freq_featGroup)), ncol = length(names(freq_featGroup)), dimnames = list(names(freq_featGroup),names(freq_featGroup))) # Empty adjacency matrix
  for(i in 1:length(allRuleSets)){ # for ith rule set
    ruleSet <- allRuleSets[[i]]
    for(j in 1:length(ruleSet)){# for jth rule
      rule <- ruleSet[[j]]
      # Check that there are more than 1 rule
      if(length(rule) > 1){
        # Connections between rules
        for(k in 1:(length(rule)-1)){ # for kth feature in rule, up to 2nd to last feature
          feat1 <- strsplit(rule[k], "_")[[1]][1]
          for(l in (k+1):length(rule)){ # for each rule after the kth rule
            feat2 <- strsplit(rule[l], "_")[[1]][1]
            # Order by frequency (link goes from most to least frequent; alphabetical if even frequency)
            if(freq_featGroup[feat1] > freq_featGroup[feat2]){
              feats <- c(feat1, feat2)
            } else if(freq_featGroup[feat2] > freq_featGroup[feat1]){
              feats <- c(feat2, feat1)
            } else{
              feats <- sort(c(feat1, feat2))
            }
            amat_group[feats[1], feats[2]] <- amat_group[feats[1], feats[2]] + 1
            #amat_group[feat2, feat1] <- amat_group[feat2, feat1] + 1 # Comment out to eliminate duplicate entries
          }
        }
        # if only 1 rule, then increment entry w/ itself as row and col
      } else {
        feat1 <- strsplit(rule[1], "_")[[1]][1]
        amat_group[feat1, feat1] <- amat_group[feat1, feat1] + 1
      }
    }
  }

  for(i in 1:nrow(amat_group)){
    label <- getLabel(rownames(amat_group)[i], featureGroups)
    colnames(amat_group)[i] <- label
    rownames(amat_group)[i] <- label
  }

  return(amat_group)
}

