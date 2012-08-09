\name{enumerateBVS}
\alias{enumerateBVS}
\title{Function to Enumerate all models for Bayesian Variant Selection Methods}
\description{
  This function enumerates and calculates summaries for all models in the model space.  Not recommended for problems where p>20.  
}

\usage{
enumerateBVS(data,forced=NULL,cov=NULL,a1=0,rare=FALSE,mult.regions=FALSE,
             regions=NULL,hap=FALSE,inform=FALSE)
}

\arguments{
  
  \item{data}{ 
    a \emph{(n x (p+1))} dimensional data frame where the first column corresponds to the response variable that is presented as a factor variable corresponding to an individuals disease status
    (0|1),and the final \emph{p} columns are the SNPs of interest each coded as a numeric variable that corresponds to the number of copies of minor alleles (0|1|2)}

  \item{forced}{ 
    an optional \emph{(n x c)} matrix of \emph{c} confounding variables that one wishes to adjust the analysis for and that will be forced into every model.}
  
  \item{inform}{
    if inform=TRUE corresponds to the iBMU algorithm of Quintana and Conti (Submitted) that incorporates user specified external predictor-level covariates into the variant selection algorithm.}

   \item{cov}{ 
    an optional \emph{(p x q)} dimensional matrix of \emph{q} predictor-level covariates that need to be specified if inform=TRUE that the user wishes to incorporate into the estimation of the marginal inclusion probabilities using the iBMU algorithm}

  \item{a1}{
    a \emph{q} dimensional vector of specified effects of each predictor-level covariate to be used when inform=TRUE.}

  \item{rare}{
    if rare=TRUE corresponds to the Bayesian Risk index (BRI) algorithm of Quintana and Conti (2011) that constructs a risk index based on the multiple rare variants within each model.  The marginal likelihood of each model is then calculated based on the corresponding risk index.}

  \item{mult.regions}{
    when rare=TRUE if mult.regions=TRUE then we include multiple region specific risk indices in each model.  If mult.regions=FALSE a single risk index is computed for all variants in the model.}

  \item{regions}{
    if mult.regions=TRUE regions is a \emph{p} dimensional character or factor vector identifying the user defined region of each variant.}

  \item{hap}{
    if hap=TRUE we estimate a set of haplotypes from the multiple variants within each model and the marginal likelihood of each model is calculated based on the set of estimated haplotypes.}}


\value{
  This function outputs a list of the following values:
 
  \item{fitness}{
    A vector of the fitness values (log(Model likelihood) - log(Model Prior)) of each 
    enumerated model. }
 
  \item{logPrM}{
    A vector of the log Model Priors of each enumerated model.}
  
  \item{which}{
    A vector identifying the character representation of each model indicator vector.}

  \item{coef}{
    If rare=FALSE we report a matrix where each row corresponds to the estimated coefficients
    for all variables within each enumerated model.  If rare=TRUE we 
    report a vector where each entry corresponds to the estimated coefficient of the risk index (or multiple risk indices if mult.regions = TRUE) 
    corresponding to each enumerated model.}

  \item{alpha}{
    If inform=FALSE that is simply a vector of 0's.  If inform=TRUE we report a matrix where each row
    corresponds to the specified effects (alpha's) of each predictor-level covariate for each enumerated
    model. }
  

}


\references{ 
  Quintana M, Conti D (2011). \emph{Incorporating Model Uncertainty in Detecting Rare Variants: 
  The Bayesian Risk Index}. Genetic Epidemiology 35:638-649.

  Quintana M, Conti D (Submitted). \emph{Integrative Variable Selection via Bayesian Model 
  Uncertainty}.
}

\author{Melanie Quintana <maw27.wilson@gmail.com>}

\examples{
## Load the data for Rare variant example
data(RareData)

## Enumerate model space for a subset of 5 variants and save output to BVS.out 
## for rare variant example.
RareBVS.out <- enumerateBVS(data=RareData[,1:6],rare=TRUE)
}

\keyword{ model enumeration }

