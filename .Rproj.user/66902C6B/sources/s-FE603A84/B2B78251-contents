## Functions for making t-SNE plot


#' Make a t-SNE plot
#'
#' Makes a t-SNE plot with color coding for classification correctness and shading for predicted outcome
#' @param df dataframe of binary features
#' @param Y vector of binary outcome
#' @param A rule set to use for determining which variables will be used for training t-SNE (will use outcome plus all variables that appear in \code{A})
#' @param caseColors vector colors to use for classification correctness. The first element is for correct, second is for incorrect
#' @param boxColor color for shading
#' @param pointSize \code{cex} paramater for \code{plot} for size of points
#' @param textSize \code{cex} paramater for \code{legend} for size of text
#' @param symb \code{pch} paramater for \code{plot} for symbol of points
#' @param bottom_buffer amount by which to shift bottom of graph up for legend
#' @param all_buffer buffer on all sides of plot
#' @param legend_under_plot whether legend is under plot or inside the plot. Overrides \code{legend_coordinates}
#' @param legend_bg_col background color of legend
#' @param legend_offset distance to move legend (as a percentage). Negative values move legend to the left and down, positive values move legend to the right and up. If \code{legend_under_plot=F}, then this argument is analagous to the inset arguement for \code{\link[graphics]{plot}}
#' @param legend_position either x,y coordinates or keyword for legend position. Only used if \code{legend_under_plot=F}. Passed to the \code{x,y} argument for \code{\link[graphics]{legend}}
#' @param legend_position position of legend
#' @param shade_name label in legend for shaded points
#' @param jitter_factor factor by which to jitter points
#' @param jitter_amount amount by which to jitter points
#' @param max_iter maximum number of iterations to run t-SNE
#' @return t-SNE plot
#' @export
plot_tsne <- function(df, Y, A, caseColors, boxColor, pointSize=1, textSize=1, symb=19,
                      bottom_buffer=1.25, all_buffer=1,
                      legend_under_plot=T, legend_bg_col="transparent", legend_offset=c(0,0),
                      legend_position="bottomright", shade_name="Positive classification",
                      jitter_factor=1, jitter_amount=NULL,
                      max_iter=1000){

  ## Train tsne using outcome and features in A
  incl <- c() # variables to include
  for (rule in A) {
    for (cond in rule) {
      split <- strsplit(cond, "_")[[1]]
      is_neg <- as.numeric(split[length(split)]=="neg")
      incl <- c(incl, paste(split[1:(length(split)-is_neg)], collapse="_"))
    }
  }

  train <- cbind.data.frame(Y, df[, colnames(df)[which(colnames(df) %in% incl)]])


  if(F){
  for(i in 1:length(A)){
    for(j in 1:length(A[[i]])){
      incl <- c(strsplit(A[[i]][j], "_")[[1]][1], incl) # append new feature to beginning of list of features to be included
      if(is.na(match(incl[1], incl[-1]))){ # Make sure feature isn't already included
        for(k in 1:ncol(df)){
          feat <- strsplit(names(df)[k], "_")[[1]][1]
          if(feat==incl[1]){
            train <- cbind(train, df[,k])
          }
        }
        #train %>% dplyr::mutate(!!feat := df[[feat]])
      }
    }
  }
  names(train) <- 1:ncol(train) # done so that no two columns share a name
  }


  tsne <- Rtsne::Rtsne(train, dims=2, perplexity=5, verbose=F, max_iter=max_iter, check_duplicates=F)#, partial_pca=T)
  Yhat <- rep(0, length(Y)) # Classification based on A

  ## Jitter
  tsne$Y <- jitter(tsne$Y, factor=jitter_factor, amount=jitter_amount)
  ## Standardize coordinates
  rangeX <- max(tsne$Y[,1])-min(tsne$Y[,1])
  rangeY <- max(tsne$Y[,2])-min(tsne$Y[,2])
  tsne$Y[,1] <- tsne$Y[,1]/rangeX*100
  tsne$Y[,2] <- tsne$Y[,2]/rangeY*100

  ## Dimensions for boxes for shading
  boxWidth = (max(tsne$Y[,1]) - min(tsne$Y[,1]))/25
  boxHeight = (max(tsne$Y[,2]) - min(tsne$Y[,2]))/25

  ## Make shading slightly transparent
  #boxColor <- add.alpha(boxColor, alpha=.7)

  ## Plot boxes
  par(mar=c(bottom_buffer,0,0,0)+all_buffer)
  plot(2*max(tsne$Y), xlim=c(min(tsne$Y[,1]), max(tsne$Y[,1])), ylim=c(min(tsne$Y[,2]), max(tsne$Y[,2])),
       xlab="", ylab="",
       xaxt='n', yaxt='n')
  for(i in 1:length(A)){ # For each rule we want to plot
    # One rule
    if(length(A[[i]])==1){
      # Features and values
      split <- strsplit(A[[i]], "_")[[1]]
      is_neg <- as.numeric(split[length(split)]=="neg")
      feature <- paste(split[1:(length(split)-is_neg)], collapse="_")
      value <- as.numeric(!is_neg)

      if(F){
      split_rule <- strsplit(A[[i]], "_")
      feature <- paste(split_rule[[1]][1], split_rule[[1]][2], sep="_")
      value <- as.numeric(length(split_rule[[1]])==2) # if length 2, then no _neg, so positive (1) value
      }

      tsne_pos <- tsne$Y[df[,feature]==value,]
      Yhat[Yhat==0] <- as.numeric(df[,feature]==value)[Yhat==0] # If satisfies rule, then classified as positive (if doesn't, not recoded)

      # Multiple rules
    } else {
      split_rules <- strsplit(A[[i]], "_")
      features <- c()
      values <- c()
      for(j in 1:length(A[[i]])){
        split <- split_rules[[j]]
        is_neg <- as.numeric(split[length(split)]=="neg")
        features[j] <- paste(split[1:(length(split)-is_neg)], collapse="_")
        values[j] <- as.numeric(!is_neg)

        if(F) {
        features[j] <- paste(split_rules[[j]][1], split_rules[[j]][2], sep="_")
        values[j] <- as.numeric(length(split_rules[[j]])==2) # if length 2, then no _neg, so positive (1) value
        }
      }

      values_matrix <- matrix(rep(values, nrow(df)), ncol=length(values), byrow=T)
      tsne_pos <- tsne$Y[ apply(df[,features]==values_matrix, 1, all), ] # tsne coordinates of cases that satisfy the rule
      Yhat[Yhat==0] <- as.numeric(apply(df[,features]==values_matrix, 1, all))[Yhat==0]
    }
    if(!is.null(nrow(tsne_pos))){ # Check multiple rows
      for(j in 1:nrow(tsne_pos)){ # For each case that satisfies the rule
        # Draw a box around it
        rect(xleft=tsne_pos[j,1]-boxWidth/2, xright=tsne_pos[j,1]+boxWidth/2,
             ybottom=tsne_pos[j,2]-boxHeight/2, tsne_pos[j,2]+boxHeight/2,
             border=NA,
             lwd=10,
             col=boxColor)
        #density=boxDensity, lwd=.5, angle=boxAngles[i])
      }
    } else { # Vector, only 1 "row"
      rect(xleft=tsne_pos[1]-boxWidth/2, xright=tsne_pos[1]+boxWidth/2,
           ybottom=tsne_pos[2]-boxHeight/2, tsne_pos[2]+boxHeight/2,
           border=NA,
           lwd=10,
           col=boxColor)
    }
  }

  ## Plot points
  par(new=T, mar=c(bottom_buffer,0,0,0)+all_buffer)
  plot(tsne$Y, col=ifelse(Y==Yhat, yes=caseColors[1], no=caseColors[2]), # Color coded for correct prediction
       pch=symb,#pch=symbols2[Y+1], # Symbol coded for actual value
       cex = pointSize,
       xlim=c(min(tsne$Y[,1]), max(tsne$Y[,1])), ylim=c(min(tsne$Y[,2]), max(tsne$Y[,2])),
       xaxt='n', yaxt='n', xlab="", ylab="")

  ## Legend
  if(legend_under_plot){
    legend(x=min(tsne$Y[,1])+legend_offset[1], y=min(tsne$Y[,2])-1+legend_offset[2], ncol=2,
           legend=c(shade_name, NA, "Correct classification", "Classification error"),
           bg=legend_bg_col, col=c(boxColor, NA, caseColors),
           pch=c(15, NA, rep(symb, 2)), box.lty = 0, cex=textSize,
           xpd=T)
  } else { ## legend inside plot
    legend(legend_position, inset=legend_offset,
           legend=c(shade_name, "Correct classification", "Classification error"),
           bg=legend_bg_col, col=c(boxColor, caseColors),
           pch=c(15, rep(symb, 2)), box.lty = 0, cex=textSize)
  }
}
