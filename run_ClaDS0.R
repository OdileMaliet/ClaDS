create_model=function(start_val,posterior,proposalKernel,tuning=1,pamhType="UnivarNormalApprox"){
  model=list()
  model$start_val=start_val
  model$posterior=posterior
  model$proposalKernel=proposalKernel
  model$tuning=tuning
  model$pamhType=pamhType
  model$proposalfunction2 <- function(param, tuning, indice){
    proposalKernel=model$proposalKernel
    new = proposal_kernel(param[indice], tuning, method=proposalKernel)
    hastings = new$hastings
    param[indice] = new$moves
    
    returnProp <- list(
      proposal=param,
      hastings=hastings
    )
    
    return(returnProp)
  } 
  
  if(pamhType=="UnivarNormalApprox"){
    model$proposalfunction1 <- function(param, mcmc,indice, nrow, nparam, meanMcmc, varMcmc){
      new = rnorm(1,meanMcmc[indice],sqrt(varMcmc[indice]))
      par=param[indice]
      hastings = ((((new-meanMcmc[indice])^2)/varMcmc[indice])-(((par-meanMcmc[indice])^2)/varMcmc[indice]))/2
      #   hastings = dnorm(par,meanMcmc[indice],sqrt(varMcmc[indice]),log=T)-dnorm(new,meanMcmc[indice],sqrt(varMcmc[indice]),log=T)
      #   print(paste(hastings,hastings1,sep=" ; "))
      param[indice]=new
      returnProp <- list(
        proposal=param,
        hastings=hastings
      )
      
      return(returnProp)
    } 
  }else if (pamhType=="MultivarNormalApprox"){
    model$proposalfunction1 <- function(param, mcmc,indice, nrow, nparam, meanMcmc, varMcmc){
      new = rnorm(nparam,meanMcmc,sqrt(varMcmc))
      hastings = sum(((new-meanMcmc)^2)/varMcmc)-sum(((param-meanMcmc)^2)/varMcmc)
      returnProp <- list(
        proposal=new,
        hastings=hastings
      )
      
      return(returnProp)
    } 
  }else if (pamhType=="UnivarEmpiricalDensity"){
    model$proposalfunction1 <- function(param, mcmc,indice, generate, dens){
      new = generate[[indice]](1)
      if(is.na(new)) new=param[indice]
      par=param[indice]
      hastings = log(dens[[indice]](par))-log(dens[[indice]](new))
      param[indice]=new
      returnProp <- list(
        proposal=param,
        hastings=hastings
      )
      
      return(returnProp)
    } 
  }
  return(model)
}

autoMetropolisGibbs <-
  function(model, startvalue, iterations, consoleupdates=1000, thin=100, autoOptimize=TRUE, filename,...){
    
    proposalKernel = model$proposalKernel
    if(missing(startvalue)){
      startvalue <- model$start_val
    }
    nparam <- length(startvalue)
    ## --- Required packages
    # require(data.table)
    require(coda)
    ## ---
    mcmc_step = matrix(ncol=nparam+3, nrow=1)
    post = model$posterior(startvalue)
    if(is.infinite(post$LP)) stop("\r","Check if the starting values matches the priors","\r")
    current = post$LP
    mcmc_step[1:(nparam+3)] = c(startvalue, current, post$LL, post$Pr)
    colnames(mcmc_step) <- c(1:nparam,"LP","LL","Pr")
    rownames(mcmc_step) <- 0
    hastings <- 0
    count <- 1
    mcmc_acceptance <- NULL
    
    if(missing(filename)){
      filename <- paste("mcmc.log.txt")
    }
    
    ## --- Indices for auto-optimization
    kernel = model$proposalKernel
    familyKernel = "uniform"
    if(kernel == "bactrian" | kernel == "Mbactrian" | kernel == "MbactrianLog" | kernel == "MbactrianTriangle" | kernel == "MbactrianLaplace"| kernel == "bactrianTriangle"| kernel == "bactrianLaplace") familyKernel <- "bactrian"
    param <- list(...)
    if(is.null(param[["acceptance"]])){
      Poptimal <- ifelse(familyKernel=="bactrian",0.3, 0.44)}else{Poptimal <- param$acceptance}
    if(is.null(param[["update"]])){update <- 1000 }else{update <- param$update}
    if(is.null(param[["adaptation"]])){adaptation <- update*4 }else{ adaptation <- param$adaptation}
    if(is.null(param[["verbose"]])){verbose <- TRUE }else{ verbose <- param$verbose}
    if(is.null(param[["random"]])){random <- FALSE }else{ random <- param$random}
    if(is.null(param[["blocks"]])){blocks <- 1 }else{ blocks <- param$blocks}
    
    tuning = model$tuning
    
    
    if(length(tuning)!=nparam) tuning = rep(tuning[1],nparam)
    variables = 1:nparam
    
    
    ## --- Prepare the text file (can be used with tracer)
    write.table(data.frame("state"=rownames(mcmc_step),mcmc_step), file=filename, col.names=TRUE, row.names=FALSE, sep = "\t" , quote=FALSE)
    
    ## --- if autoOptimization==TRUE; we optimize for nrounds during the burnin period
    
    if(autoOptimize==TRUE){
      nbrounds<-adaptation%/%update
      nlbatch<-adaptation/nbrounds
      
      # Acceptance matrix
      mcmc_acceptance <- matrix(ncol=nparam, nrow=nbrounds)
      
      message("\r","Proposal kernel: ",kernel," with optimal acceptance rate set to: ",Poptimal,"\r")
      message("\r","Optimization of the tuning values for the proposal (",nbrounds," rounds of size: ",nlbatch,") please wait!","\r")
      acc_val <- numeric(nparam)
      iterOptim = adaptation
      
      # iteration loop (here it's an adaptation phase we don't records the results!)
      for (i in 1:iterOptim){
        
        ## --- Loop over the variables (Metropolis-within-Gibbs) algorithm (we can make it random by sampling from a vector of randomized indices)
        for(j in variables){
          
          proposalValues = model$proposalfunction2(mcmc_step[1:nparam], tuning = tuning[j], indice=j)
          proposal = proposalValues$proposal
          hastings = proposalValues$hastings
          
          ## ------ Compute the ratio
          
          post_val = model$posterior(proposal)
          newlik = post_val$LP
          probval = newlik - current + hastings
          probab = exp(probval)
          
          # NOTE: there are several waste of times: generation of new proposals, computation of the likelihood for both the tree and priors
          
          ## ------ Evaluate the ratio
          
          if (probval > 0 || probab > runif(1)){
            current = newlik
            mcmc_step[1:nparam] = proposal
            mcmc_step[nparam+1] = newlik
            mcmc_step[nparam+2] = post_val$LL
            mcmc_step[nparam+3] = post_val$P
            acc_val[j] = acc_val[j] + 1
          }
        }# End loop over variables
        
        
        
        ## ------ Update the tuning parameter
        
        if(i %% update == 0){
          Pjump = acc_val/update
          mcmc_acceptance[count,variables] <- Pjump
          
          if(any(Pjump<0.001)==TRUE | any(Pjump>0.999)==TRUE){
            indMin = indMax = NULL
            
            if(any(Pjump<0.001)==TRUE){
              indMin <- which(Pjump<0.001)
              tuning[indMin] = tuning[indMin]/100
            }
            if(any(Pjump>0.999)==TRUE){
              indMax <- which(Pjump>0.999)
              tuning[indMax] = tuning[indMax]*100
            }
            
            indtot <- c(indMin,indMax)
            tuning[-indtot] = tuning[-indtot] * (tan((pi/2)*Pjump[-indtot])/tan((pi/2)*Poptimal))
          }else{
            
            # eq. 9 in Yang & Rodriguez 2013 - PNAS
            tuning = tuning * (tan((pi/2)*Pjump)/tan((pi/2)*Poptimal))
          }
          
          # reset the counter
          acc_val[]=0
          count=count+1
        }
        
        if( i %% 100 == 0 & verbose==TRUE){
          cat("\r","burnin MCMC in progress",i,"of",iterOptim,"please wait!","\r")
        }
        flush.console()
        
      }
      
      message("\r","Optimization terminated","\r")
      mcmc_acceptance = mcmc(mcmc_acceptance)  # save as a CODA object
    }
    
    
    ## Here if we have optimized the tuning parameters, we continue the mcmc with the burnin values as starting points
    # open a connection to append files
    conn <- file(filename, open="a")
    
    for (i in 1:iterations){
      
      ## --- Loop over the variables (Metropolis-within-Gibbs) algorithm
      j <- ifelse(random==FALSE, (i-1)%%(nparam)+1, sample(variables,blocks)) # we can put a function outside the loop to avoid the if statement
      
      proposalValues = model$proposalfunction2(mcmc_step[1:nparam], tuning = tuning[j], indice=j)
      proposal = proposalValues$proposal
      hastings = proposalValues$hastings
      
      ## ------ Compute the ratio
      
      post_val = model$posterior(proposal)
      newlik = post_val$LP
      probval = newlik - current + hastings
      probab = exp(probval)
      
      ## ------ Evaluate the ratio
      
      if (probval > 0 || probab > runif(1)){
        current = newlik
        mcmc_step[1:nparam] = proposal
        mcmc_step[nparam+1] = newlik
        mcmc_step[nparam+2] = post_val$LL
        mcmc_step[nparam+3] = post_val$Pr
      }
      
      
      ## ------- Save the chain to a file
      if( i %% thin == 0 ){
        # There is a faster way to save the data?
        write.table(mcmc_step, file=conn, append=TRUE, sep = "\t", col.names=FALSE, row.names=i, quote=FALSE)
      }
      
      if( i %% consoleupdates == 0 & verbose == TRUE){
        cat("\r","MCMC in progress",i,"of",iterations,"please wait!","\r")
      } 
      flush.console() 
    }
    
    # close connection
    close(conn)
    message("\r","Processing the mcmc chain, please wait!","\r")
    
    
    ## ------- Do the Metropolis-within-Gibbs for separate blocks?
    
    chain <- read.table(filename, sep="\t", header=TRUE, row.names=1)
    # To speed up the computations
    # chain <- fread(filename, sep="\t", header=T, skip=1)
    
    results <- list(chain=mcmc(chain), finetune=tuning, acceptance=mcmc_acceptance)
    return(results)
  }


run_ClaDS0=function(tree,name,pamhLocalName,proposalKernel="bactrian",
                       iteration=10000000, thin=20000,update=1000, adaptation=10000,
                       seed=NULL, nCPU=3){
  if(! is.null(seed)) set.seed(seed)
  pamhType="None"
  ntips=tree$Nnode+1
  nedges=2*tree$Nnode
  npar=nedges+3

  times=as.numeric(branching.times(tree))
  target=function(x) -tess.likelihood(times,x,0)
  mini=optimize(target,c(0,10))$minimum
  
  likelihood_relative=createLikelihood_ClaDS0(tree)$ll

  alphaP1=1
  betaP1=log(1.1)
  P1=alphaP1*log(betaP1)-log(gamma(alphaP1))
    
  start=c(1,1,1,rnorm(nedges,0,log(1.00000)))
    
  target <- function(x) {if(x[1]<0|any(x[c(1)]<0)){
    return(list(LL=-Inf,LP=-Inf,P=0))
  }else{
    LL=likelihood_relative(sigma=(x[1]),alpha=exp(x[2]),c(exp(x[3]),exp(x[2]+x[-(1:3)])) )
    if (is.nan(LL)) LL=-Inf
    Pr=P1-2*log(x[1])-betaP1/(x[1])
    return(list(LL=LL,LP=LL+Pr,Pr=Pr))}}
  
  model=create_model(start_val = start,posterior = target,proposalKernel = proposalKernel,tuning = 0.1,pamhType = pamhType)
  
  ptm <- proc.time()
  sampler=mclapply(1:3,function(j){set.seed(j); autoMetropolisGibbs(model, iterations=iteration, consoleupdates=1000000, 
                                                                    thin=thin, autoOptimize=TRUE, filename=paste(pamhLocalName,j,sep=""),
                                                                    update=update, adaptation=adaptation, verbose=F)},mc.cores = nCPU)
  
  rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),-c((npar+1):(npar+3))])}))
  gelman=try(gelman.diag(rep))
  if(! inherits(gelman,"try-error")) {gelman=max(gelman$psrf[,2])} else {gelman=2}
  print(gelman)
  
  while(gelman>1.05){
    sampler2=mclapply(1:3,function(j){
      set.seed(j)
      modelI=model
      modelI$tuning=sampler[[j]]$finetune
      modelI$start_val=sampler[[j]]$chain[nrow(sampler[[j]]$chain),1:npar]
      autoMetropolisGibbs(modelI, iterations=iteration, consoleupdates=1000000, thin=thin, autoOptimize=F, 
                          filename=paste(pamhLocalName,j,sep=""), 
                          update=1000,adaptation=10000,verbose=F)},mc.cores = nCPU)
    for(j in 1:3){sampler[[j]]$chain=mcmc(rbind(sampler[[j]]$chain,sampler2[[j]]$chain[-1,]))}
    rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),-c((npar+1):(npar+3))])}))
    gelman=try(gelman.diag(rep))
    if(! inherits(gelman,"try-error")) {gelman=max(gelman$psrf[,2])} else {gelman=2}
    print(gelman)
  }
  
  if(! is.null(name)){save(tree,sampler,file = paste(name))}
  return(sampler)
}
