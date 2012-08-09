\name{summaryBVS}
\alias{summaryBVS}
\title{Calculates Posterior Summaries for BVS Methods}
\description{
  This function calculates the global and marginal 
  Bayes Factors that give the strength of evidence of there being an
  association in the overall set of variants of interest, the individual genes
  of interest (if specified) and the individual variants of interest.
}

\usage{
summaryBVS(BVS.out,data=data,forced=NULL,cov=NULL,burnin=1000,regions=NULL,
           rare=FALSE,mult.regions=FALSE,inform=FALSE)
}

\arguments{

  \item{BVS.out}{Output from \code{\link{sampleBVS}} or \code{\link{enumerateBVS}} }
  
  \item{data}{ 
    an \emph{(n x (p+1))} dimensional data frame where the first column corresponds to the response variable that is presented as a factor variable corresponding to an individuals disease status
    (0|1),and the final \emph{p} columns are the SNPs of interest each coded as a numeric variable that corresponds to the number of copies of minor alleles (0|1|2)}

  \item{forced}{ 
    an optional \emph{(n x c)} matrix of \emph{c} confounding variables that one wishes to adjust the analysis for and that will be forced into every model.}

  \item{burnin}{
    an integer indicating the length of the burnin.}

  \item{regions}{
    an optional \emph{p} dimensional vector of character strings giving the names of the
    regions (example can be gene names or pathway names) for each of the variants in data set.  If a region vector is given, the function will report regional BF.}
  
  \item{inform}{
    if inform=TRUE corresponds to iBMU algorithm of Quintana and Conti (Submitted) that incorporates user specified external predictor-level covariates into the variant selection algorithm.}

  \item{cov}{ 
    an optional \emph{(p x q)} dimensional matrix of \emph{q} predictor-level covariates (needed when inform=TRUE) that the user wishes to incorporate into the estimation of the marginal inclusion probabilities using the iBMU algorithm}

  \item{rare}{
    if rare=TRUE corresponds to the Bayesian Risk index (BRI) algorithm of Quintana and Conti (2011) that constructs a risk index based on the multiple rare variants within each model.  The marginal likelihood of each model is then calculated based on the corresponding risk index.}

  \item{mult.regions}{
    when rare=TRUE if mult.regions=TRUE then we include multiple region specific risk indices in each model.  If mult.regions=FALSE a single risk index is computed for all variants in the model.}
  
}

\details{
  Global and marginal Bayes factors (BF) are computed based on calculating the
  posterior probabilities of each of the unique models that were visited
  in \code{\link{sampleBVS}} or all models that were enumerated in \code{\link{enumerateBVS}}.  The global BF tests 
  the hypothesis that there is an association in the overall set of variants. BF's are also calculated at the regional (if regions are specified) and the variant level.  At the
  regional level, BF are computed for
  the overall evidence of at least one of the variants within the region of
  interest being associated.  Posterior estimates for the coefficients are also reported.  Finally, if inform=TRUE
  posterior estimates of the effects of the posterior-level covariates on the marginal inclusion probabilities are reported. 
}

\value{
  This function outputs a list of the following values:

  \item{Global}{Global Bayes Factor giving the strength of evidence that at least one variant within the analysis is associated 
    with the outcome of interest}

  \item{MargBF}{Marginal variant specific Bayes Factors giving the strength of evidence that each one of the variants are associated
    with the outcome of interest}

  \item{Marg.RBF}{Regional level Bayes Factors giving the strength of evidence that at least one variant within the region is associated
    with the outcome of interest}

  \item{PostAlpha}{If inform=TRUE gives that posterior estimates of the effects of the posterior-level covariates on the marginal
    inclusion probabilities.}

  \item{PostCoef}{Posterior estimates for the coefficients of each variant if rare=FALSE and of the risk index if rare=TRUE}

  \item{Which}{Matrix of the unique models as well as their prior probability and posterior probability}

  \item{Which.r}{Matrix indicating which regions are included in each of the unique models given in Which}

  \item{Coef}{Matrix indicating the coefficients of the variants (or risk index) included in each unique model}
}


\references{ 
  Quintana M, Conti D (2011). \emph{Incorporating Model Uncertainty in Detecting Rare Variants: 
  The Bayesian Risk Index}. Genetic Epidemiology 35:638-649.

  Quintana M, Conti D (Submitted). \emph{Integrative Variable Selection via Bayesian Model 
  Uncertainty}.
}

\author{Melanie Quintana <maw27.wilson@gmail.com>}

\examples{
## RARE VARIANT BRI EXAMPLE
  ## Load the data for Rare variant example
  data(RareData)

  ## Load the results from running sampleBVS on rare variant data for 100K iterations
  data(RareBVS.out)

  ## Summarize output with a burn in of 1000 iterations
  ## Results from summary found in data(RareResults)
  RareResults = summaryBVS(RareBVS.out,data=RareData,burnin=1000,rare=TRUE)

## INFORMATIVE iBMU EXAMPLE
  ##Load the data for the informative example
  data(InformData)

  ## Load the results from running sampleBVS with inform=FALSE for 100K iterations
  data(InformBVS.NI.out)

  ## Summarize output
  ## Results from summary found in data(Informresults.NI)
  
  Informresults.NI = summaryBVS(InformBVS.NI.out,data=InformData$data,burnin=1000,
                                regions=InformData$genes,inform=FALSE)

  ## Load the results from running sampleBVS with inform=TRUE for 100K iterations
  data(InformBVS.I.out)

  ## Summarize output
  ## Results from summary found in data(Informresults.I)
 
  Informresults.I = summaryBVS(InformBVS.I.out,data=InformData$data,
                               cov=as.matrix(InformData$cov),burnin=1000,
                               regions=InformData$genes,inform=TRUE)

}

\keyword{ Posterior Summaries }

