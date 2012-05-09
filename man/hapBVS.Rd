\name{hapBVS}
\alias{hapBVS}
\title{Function to estimate and report a set of haplotypes given a subset of variants}
\description{
  This function takes a subset of variants and estimates a set of haplotypes.  Only haplotypes with a frequency greater than min.Hap.freq are reported. 
}

\usage{
hapBVS(G,min.Hap.freq)
}

\arguments{

  \item{G}{ 
    an \emph{(n x g)} matrix of a subset of \emph{g} SNPs of interest that are each coded as a numeric variable that corresponds to the number of copies of minor alleles (0|1|2)}

  \item{min.Hap.freq}{
    the minimum haplotype frequency of which an estimated haplotype is reported}
  
}



\value{
  This function outputs a matrix of estimated haplotypes.
}


\author{Melanie Quintana <maw27.wilson@gmail.com>}


\keyword{ haplotypes }

