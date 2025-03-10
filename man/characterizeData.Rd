% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/characterizeData.R
\name{characterizeData}
\alias{characterizeData}
\title{characterizeData()}
\usage{
characterizeData(
  bam,
  gtf,
  fasta,
  out,
  nCPU = 4,
  numReads = NULL,
  nBam = NULL,
  seed = 1234,
  stage = "geneFilt",
  phred = 33
)
}
\arguments{
\item{bam}{a single bam file, or a directory containing several bam files used to characterize data}

\item{gtf}{gtfFile containing the reference annotation}

\item{fasta}{fastaFile containing the reference genome}

\item{out}{output directory}

\item{nCPU}{number of parallel CPU processes [default = 4]}

\item{numReads}{max number of reads in bam file to be used for data characterization, by default all reads are used}

\item{nBam}{number of bam files used for bias computation, only used if multiple bam files are used, by default all bam files are used}

\item{seed}{random seed for reproducibility [default = 1234]}

\item{stage}{stage from where to start running: geneFilt, readBams, plotStats [default = geneFilt]}

\item{phred}{phred encoding [default = 33]}
}
\value{
object of class "DataCharacteristics"
}
\description{
characterize Smart-seq3 dataset
}
