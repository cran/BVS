plotBVS = function(results, num.models=100, num.snps=20, plot.coef=FALSE, true.coef=NULL,main=NULL, genes=NULL, type="s",...) {
  
  which <- data.frame(results$Which[,1:(dim(results$Which)[2]-2)])
  for(i in 1:dim(which)[2]){which[,i] = as.numeric(which[,i])}
  which.g <- results$Which.g
  u.genes = unique(genes)
  snps <- colnames(results$Which)[-c((dim(which)[2]+1):(dim(which)[2]+3))]
  p = length(snps)
  num.snps = min(p,num.snps)
  num.models = min((dim(which)[1]-1),num.models)
  postprob <- as.numeric(results$Which[,dim(results$Which)[2]])
  bf <- as.numeric(results$MargBF)
  bf.g <- as.numeric(results$Marg.GBF)
  if(plot.coef==TRUE){
  	coef = results$Coef
  }

  null.ind = c(1:dim(which)[1])[apply(which,1,paste,collapse="")==paste(rep(0,p),collapse="")]
  null.post <- postprob[null.ind]
  which = as.matrix(which[-null.ind,])
  if(length(which.g)>0){
  	  which.g = as.matrix(which.g[-null.ind,])}  
  postprob = postprob[-null.ind]
  model.ord <- order(-postprob)
  which <- as.matrix(which[model.ord, ])
  if(length(which.g)>0){
      which.g <- as.matrix(which.g[model.ord, ])}  
  postprob <- postprob[model.ord]
  if(plot.coef==TRUE){
  	coef = coef[-null.ind]
  	coef = coef[model.ord]
  }	
  

  ## Graphic Parameters
  if(type=="s"){
  nmodel <- num.models
    nvar <- num.snps
    ordr <- order(-bf)
    snps = snps[ordr][1:nvar]
    genes = genes[ordr][1:nvar]
    rownms <- paste(snps,genes,sep="\n")
    clr <- c("#FFFFFF", "#A020F0", "#0000CD")
    ordr <- ordr[1:nvar]
    color.matrix <- which[1:(nmodel), ordr[1:nvar]] + 2
    if(length(true.coef)>0){
        prob.labels <- paste("Marg BF:",as.numeric(round(bf[ordr[1:nvar]], 2)),
                             " \nTrue OR:",round(true.coef[ordr[1:nvar]],2),sep="")}
    if(length(true.coef)==0){
    	prob.labels <- paste("Marg BF:",as.numeric(round(bf[ordr[1:nvar]], 2)))}
    

    
  ## matrix of colors white, purple(la), blue(dom), red(rec)    
  keep.mar <- par(mar=c(5, 6, 4, 2) + 0.1)
  par(las=1, mar=c(8, 10, 5, 10), ps=10, font=2)
  maintitle = main
  if(length(main)==0){
     maintitle=paste("SNP Inclusions of Top Models \nGlobal BF=",round(results$Global,1))}
  prob.axis <- postprob[1:(nmodel)]
  prob.axis <- prob.axis/sum(prob.axis) 
  image(c(0, cumsum(prob.axis)), 1:nvar, as.matrix(color.matrix), col=clr, 
        xlab="Models Ranked by Post. Prob.", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1),
        main=maintitle)

  xat <- (cumsum(prob.axis) + c(0, cumsum(prob.axis[-nmodel]))) / 2
  
  if(plot.coef==FALSE){
      axis(1, at=xat, labels=c(1:num.models))}
  if(plot.coef==TRUE){
  	  beta.labels <- round(coef[1:(nmodel)],1)
  	  if(nmodel>5){
         beta.labels[5:nmodel] = NA}
  	  axis(1, at=xat,labels=beta.labels)
  }    
  axis(2, at=1:nvar, labels=rownms)
  axis(4, at=1:nvar, labels=prob.labels)
  par(mar=keep.mar)}
  
  if(type=="g"){
    nmodel <- num.models
    nvar <- length(u.genes)
    ordr <- order(-bf.g)
    u.genes = u.genes[ordr]
    rownms <- u.genes
    clr <- c("#FFFFFF", "#A020F0", "#0000CD")
    ordr <- ordr[1:nvar]
    which.g[which.g>1] = 1
    color.matrix <- which.g[1:(nmodel), ordr[1:nvar]] + 2
    if(length(true.coef)>0){
        prob.labels <- paste("Marg BF:",as.numeric(round(bf.g[ordr[1:nvar]], 2)),
                             " \nTrue OR:",round(true.coef[ordr[1:nvar]],2),sep="")}
    if(length(true.coef)==0){
    	prob.labels <- paste("Marg BF:",as.numeric(round(bf.g[ordr[1:nvar]], 2)))}
    

    
  ## matrix of colors white, purple(la), blue(dom), red(rec)    
  keep.mar <- par(mar=c(5, 6, 4, 2) + 0.1)
  par(las=1, mar=c(8, 10, 5, 10), ps=10, font=2)
  maintitle = main
  if(length(main)==0){
     maintitle=paste("Gene Inclusions of Top Models \nGlobal BF=",round(results$Global,1))}
  prob.axis <- postprob[1:(nmodel)]
  prob.axis <- prob.axis/sum(prob.axis) 
  image(c(0, cumsum(prob.axis)), 1:nvar, as.matrix(color.matrix), col=clr, 
        xlab="Models Ranked by Post. Prob.", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1),
        main=maintitle)

  xat <- (cumsum(prob.axis) + c(0, cumsum(prob.axis[-nmodel]))) / 2
  if(plot.coef==FALSE){
      axis(1, at=xat, labels=c(1:num.models))}
  if(plot.coef==TRUE){
  	  beta.labels <- round(coef[1:(nmodel)],1)
  	  if(nmodel>5){
         beta.labels[5:nmodel] = NA}
  	  axis(1, at=xat,labels=beta.labels)
  }  
  axis(2, at=1:nvar, labels=rownms)
  axis(4, at=1:nvar, labels=prob.labels)
  par(mar=keep.mar)}
}


	

	

