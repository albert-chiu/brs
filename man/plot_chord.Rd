% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chord.R
\name{plot_chord}
\alias{plot_chord}
\title{Make Chord Diagram}
\usage{
plot_chord(
  fit,
  featureNames,
  featureGroups,
  linkColors,
  gridColors,
  bgLinkColor = "gray",
  maxLen,
  lengths = c(1:maxLen),
  dims = c(1, maxLen),
  minProp = 0,
  textSize = 1,
  line_arg = 1,
  side_mar = 0,
  top_mar = 0
)
}
\arguments{
\item{fit}{output from function \code{BRS}}

\item{featureNames}{names of features as they appeared in df}

\item{featureGroups}{dataframe whose first column is the features (only the prefixes, not the values) in \code{featureNames} and whose second column is the corresponding labels to display}

\item{linkColors}{vector of colors to use for chords}

\item{gridColors}{vector of colors to use for arcs}

\item{minProp}{minimum proportion of times a relationship must appear to be plotted}

\item{textSize}{\code{cex} parameter to be passed to \code{par}}
}
\value{
chord diagram
}
\description{
Makes a chord diagram representing the frequency of interactions between features
}
