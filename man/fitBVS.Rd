\name{fitBVS}
\alias{fitBVS}
\title{Function to calculate fitness for each model for Bayesian Variant Selection Methods}
\description{
  This function takes one of the models and calculates the fitness/cost value of the model.  
}

\usage{
fitBVS(Z,data,forced=NULL,cov=NULL,a1=NULL,rare=FALSE,mult.regions=FALSE,
       regions=NULL,hap=FALSE,inform=FALSE,which=NULL,which.char=NULL)
}

\arguments{
  \item{Z}{ a \emph{p} dimensional vector specifying a model of interest.  In particular if the jth value of the vector is \emph{0} the jth variant is not included in the model
    and if the jth value of the vector is \emph{1} the jth variant is included in the model.}
  
  \item{data}{ 
    a \emph{(n x (p+1))} dimensional data frame where the first column corresponds to the response variable that is presented as a factor variable corresponding to an individuals disease status
    (0|1),and the final \emph{p} columns are the SNPs of interest each coded as a numeric variable that corresponds to the number of copies of minor alleles (0|1|2)}

  \item{forced}{ 
    an optional \emph{(n x c)} matrix of \emph{c} confounding variables that one wishes to adjust the analysis for and that will be forced into every model.}
  
  \item{inform}{
    if inform=TRUE corresponds to the iBMU algorithm of Quintana and Conti (Submitted) that incorporates user specified external predictor-level covariates into the variant selection algorithm.}

  \item{cov}{ 
    an optional \emph{(p x q)} dimensional matrix of \emph{q} predictor-level covariates (need when inform=TRUE) that the user wishes to incorporate into the estimation of the marginal inclusion probabilities using the iBMU algorithm}

  \item{a1}{
    a \emph{q} dimensional vector of specified (or sampled) effects of each predictor-level covariate to be used when inform=TRUE.}

  \item{rare}{
    if rare=TRUE corresponds to the Bayesian Risk index (BRI) algorithm of Quintana and Conti (2011) that constructs a risk index based on the multiple rare variants within each model.  The marginal likelihood of each model is then calculated based on the corresponding risk index.}

   \item{mult.regions}{
    when rare=TRUE if mult.regions=TRUE then we include multiple region specific risk indices in each model.  If mult.regions=FALSE a single risk index is computed for all variants in the model.}

   \item{regions}{
    if mult.regions=TRUE regions is a \emph{p} dimensional character or factor vector identifying the user defined region of each variant.}

  \item{hap}{
    if hap=TRUE we estimate a set of haplotypes from the multiple variants within each model and the marginal likelihood of each model is calculated based on the set of estimated haplotypes.}
  
  \item{which}{
    optional current which matrix of sampled models from sampleBVS that is used to see if a model has already been sampled so that that fitness does not have to be recalculated. }

  \item{which.char}{
    optional vector that identifies that current models that have been sampled from sampleBVS that is also used to determine if a model has already been sampled.}
  
}

\details{
Uses the glm function to calculate the marginal likelihood and fitness function of the model of interest.  If rare = TRUE the marginal likelihood is based on the 
risk index produced from the subset of variants within the model of interest and if hap = TRUE the marginal likelihood is based on the estimated haplotypes produced
from the subset of variants within the model of interest.}

\value{
  This function outputs a vector of the following values:

  \item{coef}{
    If rare=FALSE we report a vector where each value corresponds to the estimated coefficients
    for all variables within the model of interest. If rare=TRUE we 
    report a value corresponding to the estimated coefficient of the risk index (or risk indices if multi.regions=TRUE) 
    corresponding to each model of interest.}  
 
  \item{fitness}{
    The value of the fitness function (log(Model likelihood) - log(Model Prior)) of the 
    model of interest.} 
 
  \item{logPrM}{
    The value of the log prior on the model of interest.}
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
p = dim(RareData)[2] -1

## Fit the Null model
fit.null = fitBVS(rep(0,p),data=RareData,rare=TRUE)}

\keyword{ fitness }

