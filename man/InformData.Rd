\name{InformData}
\alias{InformData}
\docType{data}
\title{PNAT Study-based Simulation: Informative Data.}

\description{
  PNAT study-based simulated data set of 122 variants as described in Quintana and Conti (submitted).  The
  first column represents the disease status of the individual, the
  remaining columns the counts of minor alleles (0|1|2) for each variant.  The simulation was created by usingthe genotype data from a systems-based candidate gene study of smoking cessation as part ofthe Pharmacogenetics of Nicotine Addiction and Treatment Consortium.  In particular, data set was formed from genotypes of 122 variants within 789 individuals. The 122 variants are from 7 unique gene regions and thus are comprised of a great deal of correlation between the markers within each gene.  In this simulation we assumed that the predictor-level covariate corresponding to the gene \emph{CHRNB2} was informative with regards to which variants are associated with smoking cessation.}

\usage{data(InformData)}

\value{
  A list of the following items:
\item{data}{A data set with 122 variants from 789 individuals.}
\item{cov}{A set of dummy variables indicating the gene of each variant.  This set of dummy variables is used as the predictor-level covariates within an informative analysis (inform=TRUE).}
\item{genes}{A vector indicating the gene of each variant in the data set.}
}

\references{ 
  Quintana M, Conti D (Submitted). \emph{Informing Variable Selection via Bayesian Model 
  Uncertainty}.
}


\keyword{datasets}
