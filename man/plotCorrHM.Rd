% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sajr.R
\name{plotCorrHM}
\alias{plotCorrHM}
\title{plots symmetrical correlation heatmap}
\usage{
plotCorrHM(
  data = NULL,
  cor = NULL,
  method = "sp",
  main = "",
  norm = T,
  zeroIncenter = norm,
  cols = NULL,
  reorder = NULL,
  hclust.method = "complete",
  ...
)
}
\arguments{
\item{data}{data (set of columns) to be used to plot heatmap.}

\item{cor}{correlation matrix to be used instead of data}

\item{method}{correlation method}

\item{main}{graphical parameter}

\item{norm}{specifies whether rows should be normolized}

\item{zeroIncenter}{specifies whether color cols\link{length(cols)/2} should correspons to zero correlation}

\item{cols}{colors used to plot heatmap}

\item{reorder}{heatmap parameter}

\item{hclust.method}{}

\item{...}{other parameters to be passed to heatmap function}
}
\description{
plots symmetrical correlation heatmap
}
