sampleBVS = function(data,forced=NULL,inform=FALSE,cov=NULL,rare=FALSE,hap=FALSE,
                     iter=10000,save.iter=0,outfile=NULL,status.file=NULL){
	library(MASS)
	library(msm)
   
    	
	##Set up sampled model matrix 
	M = dim(data)[2]-1
	c = dim(cov)[2]
	if(length(cov)==0){
		c=1}
	if(rare==FALSE){
		which.ind = 1:M} 
	if(rare==TRUE){
	    which.ind = 1}
	which = matrix(NA,nrow=iter+1,ncol=M+length(which.ind)+2+c)
	if(c>0){
	  colnames(which) = c(colnames(data[,-1]),paste("Coef.",which.ind,sep=""),paste("Alpha",1:c,sep=""),"Fitness","PrM")}
	models.char = NULL
	
	
	##Initialize parameters
	Z.current = rep(0,M)
	Z.current[sample(1:M,5)]=1
	a0 = qnorm((1-2^(-1/M)))
	a1.current = rep(0,c)
	lower.bound = rep(0,M)
	lower.bound[Z.current==0] = -Inf
	upper.bound = rep(Inf,M)
	upper.bound[Z.current==0] = 0
	t.current = rtnorm(M,lower=lower.bound,upper=upper.bound)
	fit.current = fitBVS(Z.current,data=data,forced=forced,cov=cov,a1=a1.current,rare=rare,
	                     hap=hap,inform=inform,which=which,models.char=models.char)
	which[1,] = c(Z.current,fit.current[which.ind],a1.current,fit.current[-which.ind][1:2])
	models.char = paste(Z.current,collapse="")
	
	
		
	##Run mutation for iter
	for(i in 1:iter){
		if(length(status.file)>0){
		   cat("Processing iter",i,"\n",file=status.file,append=TRUE)}
		if(length(status.file)==0){
		   cat("Processing iter",i,"\n")}   
		   
		
		##Sample new a1 if we are using the informative prior.  If not keep a1=0
		if(inform==TRUE){
		  ##Sample a1 | current model and current t
	      V.hat = solve(diag(1,c,c)+t(cov)%*%cov)
	      alpha.hat = V.hat%*%(t(cov)%*%(t.current-a0))
	      a1.current = mvrnorm(1,mu=alpha.hat,Sigma=V.hat)
	    
	      ##Sample latent variable | current a1 and current model 
	      lower.bound = rep(0,M)
	      lower.bound[Z.current==0] = -Inf
	      upper.bound = rep(Inf,M)
	      upper.bound[Z.current==0] = 0
          mean = a0+as.matrix(cov)%*%as.matrix(a1.current)
	      t.current = rtnorm(M,mean=mean,lower=lower.bound,upper=upper.bound)
	      
	      ##Calculate new current fitness | current a1
	      fit.current = fitBVS(Z.current,data=data,forced=forced,cov=cov,a1=a1.current,rare=rare,
	                           hap=hap,inform=inform,which=which,models.char=models.char)}
	    
	      ##Sample Model
	      model.type=0:1
		  mutate.ind = sample(1:M,1)
		  Z.new = Z.current
		  Z.new[mutate.ind] = sample(model.type[!model.type==Z.current[mutate.ind]],1)
		  
		
		  ##Calculate acceptance ratio 
		  fit.new = fitBVS(Z.new,data=data,forced=forced,cov=cov,a1=a1.current,rare=rare,
		                   hap=hap,inform=inform,which=which,models.char=models.char)
		  a = min(1,exp(fit.current[length(which.ind)+1]-fit.current[length(which.ind)+1])/exp(as.numeric(fit.new[length(which.ind)+1])-fit.current[length(which.ind)+1]))
		  rand = runif(1)
		  ##If rand is less than or equal to a accept the new model and replace Z.current and fit.current
		  if(rand<=a){
		  	Z.current = Z.new
		  	fit.current = fit.new}	
		
		  	
		##Place fitness results of iteration in which matrix
		which[c(1:dim(which)[1])[is.na(which[,1])][1],] = c(Z.current,fit.current[which.ind],a1.current,fit.current[-which.ind][1:2])
		models.char = c(models.char,paste(Z.current,collapse=""))
		
		  
		##Save results after every save.iter iteration if save.iter>0
		  
		if(i %% save.iter ==0 & save.iter>0){
			results = list(as.numeric(which[is.na(which[,1])==FALSE,(M+1+length(which.ind)+c)]),as.numeric(which[is.na(which[,1])==FALSE,(M+2+length(which.ind)+c)]),
	                       which[is.na(which[,1])==FALSE,1:M],
	                       which[is.na(which[,1])==FALSE,(which.ind+M)],which[is.na(which[,1])==FALSE,(M+length(which.ind)+1):(M+length(which.ind)+c)])
	        names(results) = c("fitness","PrM","which","coef","alpha")
	        save(results,file=outfile)
			}}
	
	##Return final results		
	which = which[is.na(which[,1])==FALSE,]
	results = list(as.numeric(which[,(M+length(which.ind)+1+c)]),as.numeric(which[,(M+length(which.ind)+2+c)]),which[,1:M],
	               which[,(which.ind+M)],which[,(M+length(which.ind)+1):(M+length(which.ind)+c)])
	names(results) = c("fitness","PrM","which","coef","alpha")
	return(results)
	}