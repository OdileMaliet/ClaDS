library(parallel)

proposalGeneratorFactoryDE_gibbs <- function(proba.gibbs,p=0.01,var=1e-6,burn=0,n.thin=0,decreasing.var=0,alpha=log(0.1)/1000,N.sigma=0,allHyp=F){
  
  returnProposal <- function(chains,x,n){
    N=n*(burn)
    npar=ncol(chains[[1]])-2
    if(length(var)==1){
      var=rep(var,npar)
    }
    if(length(decreasing.var)==1){
      decreasing.var=rep(decreasing.var,npar)
    }
    gibbs=sample(1:npar,sample(1:npar,1,prob = proba.gibbs))
    if(allHyp){
      if(1 %in% gibbs | 2 %in% gibbs | 3 %in% gibbs){
      gibbs=1:npar
    }}
    
    ind1.1=sample(1:length(chains),1)
    ind2.1=sample(1:length(chains),1)
    if(n.thin==0 | (n-N)<=n.thin){
      
      if(floor(N)==0){
        
        ind1.2=sample(1:n,1)
        ind2.2=sample(1:n,1)
      }else{
        ind1.2=sample(N:n,1)
        ind2.2=sample(N:n,1)
      }
    }else{
      ind1.2=sample(seq(floor(N),n,length.out = n.thin),1)
      ind2.2=sample(seq(floor(N),n,length.out = n.thin),1)
    }
    # }
    u=runif(1,0,1)
    
    if(u>p){
      # gamma=2.36/sqrt(2*npar)
      # gamma=runif(1,0,2*(2.36/sqrt(2*length(gibbs))))
      gamma=rexp(1,1/(2.36/sqrt(2*length(gibbs))))
    }else{
      gamma=1
    }
    
    if(any(decreasing.var[gibbs]>0)){
      e=rnorm(length(gibbs),0,var[gibbs]+decreasing.var[gibbs]*exp(alpha*n))
    }else{
      e=rnorm(length(gibbs),0,var[gibbs])
    }
    
    
    rep=x
    # print(x)
    rep[gibbs]=rep[gibbs]+e+(gamma*(chains[[ind1.1]][ind1.2,1:npar]-chains[[ind2.1]][ind2.2,1:npar]))[gibbs]
    if(n<N.sigma){rep[1]=x[1]}
    
    return(rep)
  }
  
  return(list(returnProposal=returnProposal))
}


mcmcSamplerDE_gibbs <- function(likelihood,proba.gibbs, Nchain=3, prior = NULL, startvalue, startmodel = NULL, iterations=10000, proposalGenerator = NULL, consoleupdates=1000000, thin = NULL){
  
  require(coda)
  require(compiler)
  require(MASS)
  
  ###############################################
  # Target definitions 
  
  if (is.null(prior)){
    prior <- function(x){
      return(0) 
    }
  }
  
  
  numPars = length(startvalue[[1]]) + length(startmodel)
  
  if(is.null(proposalGenerator)){
    proposalGenerator = proposalGeneratorFactory(rep(1,numPars))
  }
  
  ####### CREATE CHAIN
  
  chains=list()
  currentLPs=list()
  former=list()
  for(i in 1:Nchain){
    if(i==1 | !identical(startvalue[[i]],startvalue[[1]])){
      chain = array(dim = c(1,numPars+2))
      chain[1,1:numPars] = c(startvalue[[i]], startmodel)
      colnames(chain) = c(1:numPars, "LL", "LP")
      form=likelihood(startvalue[[i]])
      chain[1, (numPars+1):(numPars+2)] = c(form$LP,form$LP)
    }else{
      chain=chains[[1]]
      form=former[[1]]
    }
    
    former[[i]]=form
    currentLPs[[i]] = chain[1, (numPars+2)]
    chains[[i]]=chain
  }
  
  
  ##### Sampling
  
  classFields = list(
    likelihood = likelihood, 
    
    prior = prior,
    
    startvalue = startvalue, 
    numPars = numPars,
    indexLL = numPars + 1,
    indexLP = numPars + 2,
    currentLPs = currentLPs,
    chains = chains, 
    optimize = optimize, 
    consoleupdates=consoleupdates, 
    thin = thin,
    Nchain=Nchain,
    proposalGenerator = proposalGenerator,
    codaChain = NULL,
    post=likelihood,
    former=former
  )
  
  class(classFields) <- append(class(classFields),"mcmcSampler")
  return(classFields)
}

prepare_ClaDS=function(tree,sample_fraction, Nchain=3,nlambda = 1000,nt = 30,l0=0.1,s0=1,model_id="ClaDS2"){
  nedges=nrow(tree$edge)
  # the likelihood and posterior functions
  if(model_id == "ClaDS1"){
    likelihood = createLikelihood_ClaDS1(tree,nlambda = nlambda,nt = nt, conv=1e-1) # the likelihood function
  }else{ 
    likelihood = createLikelihood_ClaDS2(tree,nlambda = nlambda,nt = nt, conv=1e-1) # the likelihood function
  }
  relToAbs = likelihood$relToAbs
  likelihood = likelihood$ll
  prior=function(x){0}  # the prior function, here a flat prior
  
  alpha_effect=relToAbs(c(0,rep(1,tree$Nnode*2)))[-1]
  
  posterior=function(param,former=NULL){        # the posterior function
    param2=param
    param2[1]=exp(param2[1])
    param2[2]=exp(param2[2])
    param2[-(1:4)]=param2[4]+param[-(1:4)]+alpha_effect*param[2]
    t=try(likelihood(param2,sample_fraction,former))
    if(inherits(t,"try-error")){
      return(list(LP=-Inf))
    }else if(is.nan(t$LL) | is.infinite(t$LL)){
      return(list(LP=-Inf))
    }else{
      Pr=prior(param2)
      t$Pr=Pr
      t$LP=t$LL+Pr
      return(t)}}
  
  start=lapply(1:Nchain,function(i){c(sigma=log(s0),alpha=log(1),mu=0,lambda=c(rnorm(nedges+1,log(l0),0))) })#rel.to.abs(tree,c(log(0.001),rep(log(A[i]),nedges)))) })
  
  npar=length(start[[1]])
  g=5 #max number of parameter updated at each iteration
  G=exp(-(1:(npar))/g)
  
  testGenerator <- proposalGeneratorFactoryDE_gibbs(decreasing.var=c(rep(1e-3,2),rep(1e-1,2),rep(1e-2,nedges)),alpha=log(0.1)/200,
                                                    var=c(1e-4,1e-4,rep(1e-4,2),rep(1e-4,nedges)),
                                                    proba.gibbs=c(G))
  
  sampler=mcmcSamplerDE_gibbs(posterior,startvalue = start,proposalGenerator = testGenerator,Nchain = Nchain,consoleupdates = 1)
  sampler$alpha_effect=alpha_effect
  sampler$relToAbs=relToAbs
  return(sampler)
}

fit_ClaDS <- function(mcmcSampler, iterations, post, thin=NULL,nCPU=1,colnames=c("sigma","alpha","mu","l_0")){
  post=mcmcSampler$post
  alpha_effect=mcmcSampler$alpha_effect
  relToAbs=mcmcSampler$relToAbs
  if(is.null(mcmcSampler$thin)){
    if(is.null(thin)){
      mcmcSampler$thin=1
    }else{
      mcmcSampler$thin=thin
    }
  }else{
    if(!is.null(thin)){
      if (!mcmcSampler$thin==thin){
        warning("thinning factor has been modified")
        mcmcSampler$thin=thin
      }
    }
  }
  lastvalue = nrow(mcmcSampler$chains[[1]])
  k=lastvalue
  currentchain=lapply(1:mcmcSampler$Nchain,function(i){mcmcSampler$chains[[i]][lastvalue,]})
  former=mcmcSampler$former
  i=lastvalue+1
  
  while (i <((lastvalue+iterations/mcmcSampler$thin)+1)){
    currentchain=mclapply(1:mcmcSampler$Nchain,function(j){
      current=mcmcSampler$chains[[j]][k,]
      form=former[[j]]
      for (l in 1:mcmcSampler$thin){
        proposal = try(mcmcSampler$proposalGenerator$returnProposal(mcmcSampler$chains,current[1:mcmcSampler$numPars],k))
        proposalEval <- try(post(proposal,form),silent = T)
        if(!is.null(proposal) & !inherits(proposal,"try-error") & !inherits(proposalEval,"try-error") & inherits(proposalEval,"list")){
          if(!is.nan(proposalEval$LP) & proposalEval$LP<Inf & proposalEval$LP>-Inf){
            probab = exp(proposalEval$LP - current[mcmcSampler$numPars+2])
            
            if (is.nan(probab) | runif(1) < probab){
              current=c(proposal, proposalEval$LP,proposalEval$LP)
              form=proposalEval
            }}}
      }
      return(list(current=current, former=form))
    },mc.cores = nCPU)
    
    if(sum(!sapply(currentchain,is.null))==mcmcSampler$Nchain ) {
      mcmcSampler$chains=lapply(1:mcmcSampler$Nchain,function(j){rbind(mcmcSampler$chains[[j]],currentchain[[j]]$current)})
      mcmcSampler$former=lapply(1:mcmcSampler$Nchain,function(j){currentchain[[j]]$former})
      
      flush.console()
      k=k+1
      i=i+1
    }
    if( i %% mcmcSampler$consoleupdates == 0 ) cat("\r","MCMC in progress",(i-2)*mcmcSampler$thin,"of",iterations+mcmcSampler$thin*(lastvalue-1),"please wait!","\r")}
  
  message("Done.")
  for(i in 1:mcmcSampler$Nchain){
    colnames(mcmcSampler$chains[[i]])[1:length(colnames)]=colnames
  }
  mcmcSampler$codaChain = mcmc.list(lapply(1:mcmcSampler$Nchain,function(i){mcmc(mcmcSampler$chains[[i]])}))
  mcmcSampler$alpha_effect=alpha_effect
  mcmcSampler$relToAbs=relToAbs
  return(mcmcSampler)
}

plotChains_ClaDS = function(sampler,burn=1/2,thin=1,param=c("sigma","alpha","mu","LP")){
  npar=ncol(sampler$codaChain[[1]])-2
  plot(mcmc.list(lapply(sampler$codaChain,function(x){mcmc(x[-(1:(floor(nrow(x)*(1-burn)))),param])})))
}

getMAPS_ClaDS = function(sampler, burn=1/2, thin=1){
  nR=nrow(sampler$codaChain[[1]])
  npar=ncol(sampler$codaChain[[1]])-2
  
  chains=mcmc.list(lapply(sampler$codaChain,function(x){mcmc(x[-(1:(1*floor(nR/2))),])}))
  
  for(k in 1:length(chains)){
    for(j in 1:nrow(chains[[1]])){
      chains[[k]][j,5:npar]=chains[[k]][j,4]+chains[[k]][j,2]*sampler$alpha_effect+
        chains[[k]][j,5:npar]
    }}
  
  MAPS=sapply(1:npar, function(i){D=density(c(chains[[1]][seq(1,nrow(chains[[1]]),10),i],
                                              chains[[2]][seq(1,nrow(chains[[1]]),10),i],
                                              chains[[3]][seq(1,nrow(chains[[1]]),10),i]))
  return(D$x[which.max(D$y)])})
  
  MAPS[1:2]=exp(MAPS[1:2])
  return(MAPS)
}
