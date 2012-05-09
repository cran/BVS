enumerateBVS = function(data,forced=NULL,cov=NULL,a1=0,rare=FALSE,hap=FALSE,inform=FALSE){
	p = dim(data)[2]-1
	
	#### create all possible models
	model.type=0:1
	all.models <- lapply(vector("list", p), function(v) { model.type } )
	all.models <- expand.grid(all.models)
	
	##Get results for all models
	results <- apply(all.models,1,fitBVS,data=data,forced=forced,cov=cov,a1=a1,rare=rare,hap=hap,inform=inform)
	coef = results[1,]
	fitness = results[2,]
	PrM = results[3,]
	which = all.models
	alpha = rep(a1,dim(all.models)[1])
	results = list(fitness,PrM,which,coef,alpha)
	names(results) = c("fitness","PrM","which","coef","alpha")
	return(results)
	}