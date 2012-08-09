\name{sampleBVS}
\alias{sampleBVS}
\title{Sampling Algorithm for Bayesian Variant Selection Methods}
\description{
  This function performs a basic MH Sampling algorithm to sample models from the model space when enumeration is not possible.  For informative marginal inclusion probabilities the algorithm also performs a basic MCMC algorithm to sample the effects of the predictor-level covariates (alpha).  
}

\usage{
sampleBVS(data,forced=NULL,inform=FALSE,cov=NULL,rare=FALSE,mult.regions=FALSE,
          regions=NULL,hap=FALSE,iter=10000,save.iter=0,outfile=NULL,
          status.file=NULL,old.results=NULL)
}

\arguments{
  
  \item{data}{ 
    an \emph{(n x (p+1))} dimensional data frame where the first column corresponds to the response variable that is presented as a factor variable corresponding to an individuals disease status
    (0|1),and the final \emph{p} columns are the SNPs of interest each coded as a numeric variable that corresponds to the number of copies of minor alleles (0|1|2)}

  \item{forced}{ 
    an optional \emph{(n x c)} dimensional matrix of \emph{c} confounding variables that one wishes to adjust the analysis for and that will be forced into every model.}
  
  \item{inform}{
    if inform=TRUE corresponds to the iBMU algorithm of Quintana and Conti (Submitted) that incorporates user specified external predictor-level covariates into the variant selection algorithm.}

  \item{cov}{ 
    an optional \emph{(p x q)} dimensional matrix of \emph{q} predictor-level covariates (needed when inform=TRUE) that the user wishes to incorporate into the estimation of the marginal inclusion probabilities using the iBMU algorithm}

  \item{rare}{
    if rare=TRUE corresponds to the Bayesian Risk index (BRI) algorithm of Quintana and Conti (2011) that constructs a risk index based on the multiple rare variants within each model.  The marginal likelihood of each model is then calculated based on the corresponding risk index.}

  \item{mult.regions}{
    when rare=TRUE if mult.regions=TRUE then we include multiple region specific risk indices in each model.  If mult.regions=FALSE a single risk index is computed for all variants in the model.}

  \item{regions}{
    if mult.regions=TRUE regions is a \emph{p} dimensional character or factor vector identifying the user defined region of each variant.}

  \item{hap}{
    if hap=TRUE we estimate a set of haplotypes from the multiple variants within each model and the marginal likelihood of each model is calculated based on the set of estimated haplotypes.}
  
  \item{iter}{
    the number of iterations to run the algorithm.}
  
  \item{save.iter}{
    the number of iterations between each checkpoint. A checkpoint file is
    written every save.iter iterations.}
  
  \item{outfile}{
    character string giving the pathname of the checkpoint file 
    to save the output of the algorithm to.}
  
  \item{status.file}{
    character string giving the pathname of the file to write the status of
    the algorithm.}

  \item{old.results}{
    old output from sampleBVS that has been run for a subset of the total number of iterations that the user wanted to run.  if specified the sampling algorithm will start from the last sampled model in old.results.  To be used if sampleBVS has been interrupted for some reason.   }
  
}

\details{
  The algorithm is run for a chosen number of iterations where we randomly add and remove variants from the current model based on a basic MH algorithm.  If inform = TRUE we also incorporate a set of predictor-level covariates that are provided by the user and use a MCMC algorithm to sample the effects of the covariates on the marginal inclusion probabilities.  Convergence of the algorithm can
  be determined by running two independent runs of the algorithm with
  different starting values and examining the marginal Bayes factors for each variant under each independent run.
}

\value{
  This function outputs a list of the following values to the file
  write.out if this file is specified for every save.iter number of
  iterations:
 
  \item{fitness}{
    A vector of the fitness values (log(Model likelihood) - log(Model Prior)) of each 
    model sampled at each iteration of the algorithm. }
 
  \item{logPrM}{
    A vector of the Model Priors of each model sampled at each iteration 
    of the algorithm.}
  
  \item{which}{
    A vector identifying the character representation of each model sampled.}

  \item{coef}{
    If rare=FALSE we report a matrix where each row corresponds to the estimated coefficients
    for all variables within each model sampled at each iteration of the algorithm.  If rare=TRUE we 
    report a vector where each entry corresponds to the estimated coefficient of the risk index (or multiple risk indices if mult.regions=TRUE) 
    corresponding to each enumerated model.}

  \item{alpha}{
    If inform=FALSE that is simply a vector of 0's.  If inform=TRUE we report a matrix where each row
    corresponds to the estimated effects (alpha's) of each predictor-level covariate for each model 
    sampled at each iteration of the algorithm. }
  

}


\references{ 
  Quintana M, Conti D (2011). \emph{Incorporating Model Uncertainty in Detecting Rare Variants: 
  The Bayesian Risk Index}. Genetic Epidemiology 35:638-649.

  Quintana M, Conti D (Submitted). \emph{Informing Variable Selection via Bayesian Model 
  Uncertainty}.
}

\author{Melanie Quintana <maw27.wilson@gmail.com>}

\examples{
## Rare Variant BRI example  
  ## Load the data for Rare variant example
  data(RareData)

  ## Run algorithm for 100 iterations for rare variant example.
  ## NOTE: Results from a more realistic run with 100K
  ## iterations can be found in data(RareBVS.out).
  RareBVS.out <- sampleBVS(data=RareData,iter=100,rare=TRUE)

  ## Run algorithm for 100 iterations for multiple region rare 
  ## variant example.
  p = dim(RareData)[2]-1
  regions = c(rep("Region1",(p/2)),rep("Region2",(p/2)))
  RareBVS.out <- sampleBVS(data=RareData,iter=100,rare=TRUE,mult.regions=TRUE,regions=regions)

## Informative iBMU Example
  ##Load the data for the informative example
  data(InformData)

  ## Run algorithm for 100 iterations for informative data example.
  ## This run is the basic Bayes model uncertainty algorithm with inform=FALSE 
  ## NOTE: Results from a more realistic run with 100K
  ## iterations can be found in data(InformBVS.NI.out).
  InformBVS.NI.out = sampleBVS(InformData$data,inform=FALSE,iter=100)

  ## Run algorithm for 100 iterations for informative data example.
  ## This run corresponds to the iBMU algorithm with inform=TRUE 
  ## and dichotomous predictor-level covariates indicating the gene of each variant.  
  ## NOTE: Results from a more realistic run with 100K
  ## iterations can be found in data(InformBVS.I.out).
  InformBVS.I.out = sampleBVS(InformData$data,inform=TRUE,
                              cov=as.matrix(InformData$cov),iter=100)
}

\keyword{ model search }

