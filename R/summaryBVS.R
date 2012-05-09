summaryBVS = function(BVS.out,data=data,forced=NULL,cov=NULL,burnin=1000,genes=NULL,rare=FALSE,inform=FALSE){
	if(burnin>0){
	  which = BVS.out$which[-c(1:burnin),]
	  if(rare==FALSE){
	    coef = BVS.out$coef[-c(1:burnin),]}
	  if(rare==TRUE){
	  	coef = BVS.out$coef[-c(1:burnin)]}
	  fitness = BVS.out$fitness[-c(1:burnin)]
	  PrM = BVS.out$PrM[-c(1:burnin)]
	  if(is.matrix(BVS.out$alpha)==FALSE){
	  	a1 = BVS.out$alpha[-c(1:burnin)]}
	  if(is.matrix(BVS.out$alpha)==TRUE){
	  	a1 = BVS.out$alpha[-c(1:burnin),]}
	  }
	if(burnin==0){
	  which = BVS.out$which
	  coef = BVS.out$coef
	  fitness = BVS.out$fitness
	  PrM = BVS.out$PrM
	  a1 = BVS.out$alpha}
	  
	if(is.matrix(a1)==FALSE){
		a1 = as.matrix(a1)}
	 
	
	##Parameters
	p = dim(which)[2]
	snps = colnames(which)
	c = dim(cov)[2]
	if(length(cov)==0){c=1}
	a0 = qnorm((1-2^(-1/p)))
	
	##coef.ind
	if(rare==FALSE){coef.ind=1:p}
	if(rare==TRUE){coef.ind=1}
	   
	##Make sure Null model is in the results
	null.fit = fitBVS(rep(0,p),data=data,forced=NULL,cov=cov,a1=rep(0,c),rare=rare,inform=inform)
	if(sum(apply(which,1,paste,collapse="")==paste(rep(0,p),collapse=""))==0){
		which = rbind(rep(0,p),which)
		coef = rbind(rep(0,length(coef.ind)),coef)
		fitness = c(null.fit[length(coef.ind)+1],fitness)
		PrM = c(null.fit[length(coef.ind)+2],PrM)}
	rownames(which) = c(1:dim(which)[1])
	
	##Post expectation of alpha
	if(is.matrix(BVS.out$alpha)==FALSE){
	  	post.alpha = mean(a1)}
	if(is.matrix(BVS.out$alpha)==TRUE){
	  	post.alpha = apply(a1,2,mean)}
		
	##If inform==TRUE average prior inclusion probability across multiple values of alpha!
	if(inform==TRUE){
	  eta = t(a0+as.matrix(cov)%*%t(a1))
	  inc.prob = function(x){mean(pnorm(0,mean=x,lower.tail=FALSE))}
	  Prior.marg = apply(eta,2,inc.prob)
	  lprob.inc = log(Prior.marg)
	  lprob.ninc = log(1-Prior.marg)
	  NullPrior.marg = pnorm(0,mean=a0,lower.tail=FALSE)
	  ##Global Priors
	  Prior.null = exp(sum(lprob.ninc))	    
      Prior.alt = 1-Prior.null}
      
	if(inform==FALSE){
	  NullPrior.marg = (1/(p+1))
      Prior.null = .5
      Prior.alt =.5	
	}
	
	 
	##get unique models and recalc fitness with avg. prior inc probability if inform == TRUE
	u.which = unique(which)
	u.ind = as.numeric(rownames(u.which))
	if(rare==FALSE){
	   u.coef = coef[u.ind,]}
	if(rare==TRUE){
	   u.coef = coef[u.ind]}
	u.fitness = fitness[u.ind]
	new.lPrM = log(PrM[u.ind])
	if(inform==TRUE){   
	  new.lPrM = apply(u.which,1,function(x){t(x)%*%lprob.inc + t(1-x)%*%lprob.ninc})
	  u.fitness = fitness[u.ind] + log(PrM[u.ind]) - new.lPrM}
	
	##Calculate Posterior model probabilities
    PrMgivenD = exp(-u.fitness+min(u.fitness))/sum(exp(-u.fitness+min(u.fitness)))
	 
	
    ##Global Posterior Prob & BF
	Post.null <- PrMgivenD[apply(u.which,1,paste,collapse="")==paste(rep(0,p),collapse="")]
	Post.alt <- 1 - Post.null
	Global.BF = (Post.alt/Post.null)
	
	
	##Marginal Inclusion probabilities
	Post.marg = t(u.which!=0)%*%as.matrix(PrMgivenD)
	Post.marg.n = t(u.which==0)%*%as.matrix(PrMgivenD)
	Global.marg.BF <- log(Post.marg) - log(Post.marg.n) - log(NullPrior.marg) +  log(1-NullPrior.marg)
	which = cbind(u.which,exp(new.lPrM),PrMgivenD)
	colnames(which) = c(colnames(data)[-1],"Prior","PostProb")
	
	##Gene Inclusion probabilities
	which.g = NULL
	Gene.marg.BF = NULL
	Gene.post.inc = NULL
	if(length(genes)>0){
		u.genes = unique(genes)
		which.genes = matrix(0,nrow=p,ncol=length(u.genes))
		names(which.genes) = u.genes
		for(j in 1:length(u.genes)) {
            which.genes[genes == u.genes[j], j] <- 1
        }
        which.g <- u.which %*% which.genes
        names(which.g) = u.genes
        probne0.g <- as.numeric((PrMgivenD * 10000) %*% (which.g > 0)) / 10000
        probe0.g <- as.numeric((PrMgivenD * 10000) %*% (which.g == 0)) / 10000
        NullPrior.marg.g <- (1-NullPrior.marg)^apply(which.genes,2,sum)
        Gene.post.inc = probne0.g
        Gene.marg.BF <- log(probne0.g) - log(probe0.g) - log(1-NullPrior.marg.g) + log(NullPrior.marg.g)
        Gene.marg.BF = exp(Gene.marg.BF)
		}
	
	##Posterior estimates of coefficients and sd
	post.coef=NULL
	if(is.matrix(coef)==TRUE){
	  ##Posterior estimate of OR	
	  post.coef = t(u.coef)%*%as.matrix(PrMgivenD)
	  }
    

	##Save Results
	results=list(Global.BF,exp(Global.marg.BF),Gene.marg.BF,post.alpha,post.coef,which,which.g,u.coef)
	names(results)=c("Global","MargBF","Marg.GBF","PostAlpha","PostCoef","Which","Which.g","Coef")  
    return(results)
	}