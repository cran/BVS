\name{plotBVS}
\alias{plotBVS}

\title{ Image Plots for top Variant and Gene Inclusions }
\description{
  This function allows the user to create image plots of the top variants
  and top Genes included in the top models.  Variants and Genes are ordered based on marginal
  BF which are plotted on the right axis.  The
  width of the inclusion blocks are proportional to the posterior model
  probability that the variant or Gene is included in.
}

\usage{
plotBVS(results, num.models=100, num.snps=20, plot.coef=FALSE, true.coef=NULL,
        main=NULL, genes=NULL, type="s",...)
}

\arguments{
  \item{results}{output list from \code{\link{summaryBVS}}.}
  \item{num.models}{the number of the top models to place on the x-axis.}
  \item{num.snps}{ If inc.type="s", the number of the top variants to place
    on the y-axis.}
  \item{plot.coef}{ Only to be used for rare variant analysis when rare=TRUE.  When plot.coef=TRUE, the log(OR) of the risk indices 
    specified by each of the top models are plotted on the x axis}
  \item{type}{specifies if we want to plot the variant inclusion ("s") or
    gene inclusion ("g") }
  \item{true.coef}{optional vector giving the true odds ratios of each of the variants (if results are from a simulation)}
  \item{main}{optional vector giving the title of the plot}
  \item{genes}{
    an optional vector of character strings giving the names of the
    genes for each of the variants in data set needed when plotting type is "g"
    or can be added to include the gene names of each variant on the y axis.}
  \item{...}{General parameters for plotting functions}
}


\author{ Melanie Quintana <maw27.wilson@gmail.com> }

\examples{
## RARE VARIANT BRI EXAMPLE
  ## Load the data for Rare variant example
  data(RareData)

  ## Load the results from running sampleBVS on rare variant data for 100K iterations
  data(RareBVS.out)

  ## Load summary results 
  data(RareResults)
  
  ## Plot the variant inclusions in the top 100 models for the top 10 variants 
  plotBVS(RareResults,num.models=100,num.snps=10)

  ##Include the estimated log(OR) of the risk indices for the top models
  plotBVS(RareResults,num.models=100,num.snps=10,plot.coef=TRUE)

\donttest{
## INFORMATIVE iBMU EXAMPLE
  ##Load the data for the informative example
  data(InformData)

  ## Load the results from running sampleBVS with inform=FALSE for 100K iteration
  data(InformBVS.NI.out)

  ## Load summary results
  data(Informresults.NI)

  ## Make SNP and Gene inclusion plots
  plotBVS(Informresults.NI,num.models=50,num.snps=10,genes=InformData$genes)
  plotBVS(Informresults.NI,num.models=50,num.snps=10,genes=InformData$genes,type="g")

  ## Load the results from running sampleBVS with inform=TRUE for 100K iterations
  data(InformBVS.I.out)

  ## load summary results
  data(Informresults.I)

  ## Make SNP and Gene inclusion plots
  plotBVS(Informresults.I,num.models=50,num.snps=10,genes=InformData$genes)
  plotBVS(Informresults.I,num.models=50,num.snps=10,genes=InformData$genes,type="g")}


}


\keyword{ image plot }