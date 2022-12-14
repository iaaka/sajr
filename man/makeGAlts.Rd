% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sajr.R
\name{makeGAlts}
\alias{makeGAlts}
\title{Finds all alternative within gene}
\usage{
makeGAlts(segs, ints)
}
\arguments{
\item{segs}{data.frame with all gene segments}

\item{ints}{data.frame with all gene introns}
}
\value{
list of dataframes with all gene alternatives
}
\description{
Cuts gene by constitutive exons
}
