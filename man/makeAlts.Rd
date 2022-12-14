% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sajr.R
\name{makeAlts}
\alias{makeAlts}
\title{Splits genes into splicing alternatives}
\usage{
makeAlts(seg, ann.gff, remove.exn.ext = F)
}
\arguments{
\item{seg}{segment annotation (first element of output of \code{\link{loadSAData}} function)}

\item{ann.gff}{name of gff file}

\item{remove.exn.ext}{specifies whether alternative TSS/polyA that aren't supplied by alternative splicing (i.e. just an extension of other exons) should be removed}
}
\value{
data.frame with all found alternatives.
sites is types of splicing sites: \link{a}cceptor and \link{d}onor
segs is list of alternative segments separated by ';'. Each segment is denoted by ranks of splice sites that form segment's boundaries, rank == NA correspons to alternative TSS or polyA
ints the same as segs but for introns
}
\description{
Cuts genes by constitutive segments into regions that
can be spliced alternatively
}
