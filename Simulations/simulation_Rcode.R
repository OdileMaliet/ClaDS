# choose your working directory
path="~/Documents/ClaDS/"# here put the path to the downloaded folder
setwd(path)

#### loading the needed functions ####

library(BayesianTools)
library(vioplot)
library(apTreeshape)
library(RPANDA)
library(TESS)

source(paste0(path,"sim_ClaDS.R"))

source(paste0(path,"likelihood_ClaDS0.R"))
source(paste0(path,"utils.R"))
source(paste0(path,"proposal.R"))
source(paste0(path,"run_ClaDS0.R"))

source(paste0(path,"ClaDS1_likelihood_functions.R"))
source(paste0(path,"ClaDS2_likelihood_functions.R"))
source(paste0(path,"ClaDS3_likelihood_functions.R"))
source(paste0(path,"fit_ClaDS.R"))

# this one needs to be compiled ; for this you need to open a terminal, go to your 
# working directory, and enter
# R CMD SHLIB diversif_lognormal.c
dyn.load(paste0(path,"diversif_lognormal.so"))

#### test of the subfunctions Phi and Khi in the likelihood ####

# # # # # # # # # # # # # # # # # # # # # # # # #
#           with constant extinction            #
# # # # # # # # # # # # # # # # # # # # # # # # #

check_Psi_ClaDS1=function(t,sigma,alpha,epsilon,lambda_0,sample_frac,Nsim,mlambda=0.0001,Mlambda=100,nlambda=1000,nt=100,method="Higham08.b",NtipMax=500,conv=1e-10){
  
  phi=Phi_ClaDS1(sigma = sigma,alpha = alpha,mlambda = mlambda,Mlambda = Mlambda,nlambda = nlambda,mu = epsilon,f = sample_frac,tini=0,tf=t,by=t/nt,method=method)
  khi=Khi_ClaDS1(phi = phi,s=0,t=t,func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$mu,
                 timePhi=phi$fun[,1],nt=1000,method="Magnus",conv=conv)
  simulateKhi=0
  n=0
  i=0
  n_spec=c()
  for (j in 1:Nsim){
    cat("\r",j,"\r")
    obj=sim_ClaDS(theta=1,lamb_par=lambda_0,mu_par=epsilon,lamb_shift=alpha,sigma=sigma,condition="time",new_lamb_law = "lognormal*shift",
                  new_mu_law="uniform",prune.extinct = T,time.stop = t,taxa.stop = NtipMax, mu_min = epsilon, mu_max = epsilon)
    tree=obj$tree
    n_spec=c(n_spec,length(obj$lamb))
    if(!is.null(tree)){
      if (length(tree)==2){ntip=1}else{ntip=tree$Nnode+1}
      
      if(ntip==NtipMax){
        simulateKhi=simulateKhi+0
        n=n+1
      }else{
        simulateKhi=simulateKhi+((1-sample_frac)^(ntip-1)*sample_frac*ntip)
        n=n+(1-(1-sample_frac)^ntip)
      }
      
      i=i+1
      # n=n+1
      
    }else{
      n=n+0
    }
  }
  simulatedKhi=simulateKhi/j
  simulatedPhi=1-(n/j)
  lind=which.min(abs(phi$expLambda-lambda_0))
  computedKhi=khi[lind]
  timePhi=phi$fun[,1]
  tind=which.min(abs(timePhi-t))
  computedPhi=phi$fun[tind,lind+1]
  BD_Phi=1-(sample_frac*(lambda_0-mu)/(sample_frac*lambda_0+(lambda_0*(1-sample_frac)-mu)*exp(-t*(lambda_0-mu))))
  BD_Khi=sample_frac*((lambda_0-mu)^2)*exp(-t*(lambda_0-mu))/((sample_frac*lambda_0+(lambda_0*(1-sample_frac)-mu)*exp(-t*(lambda_0-mu)))^2)
  return(list(simulatedKhi=simulatedKhi,computedKhi=computedKhi,
              simulatedPhi=simulatedPhi,computedPhi=computedPhi,
              BD_Phi=BD_Phi, BD_Khi=BD_Khi,n_spec=n_spec))
}

dataFramePsi=data.frame()
nspec=list()

for(k in 1:500){
    set.seed(k)
    sigma=runif(1,0,1)
    alpha=runif(1,0.3,1.3)
    lambda_0=runif(1,0,0.5)
    mu=runif(1,0,1.5)*lambda_0
    sample_frac=runif(1,0,1)
    t=rexp(1,0.1)
    nlambda=300
    mlambda=1e-5
    Mlambda=100
    nt=30
    
    cat(" l0=",lambda_0,", s=",sigma,", mu=", mu,", t=",t,", f=",sample_frac, ", al=", alpha, "\n",sep="")
    cP=try(check_Psi_ClaDS1(t,sigma,alpha,epsilon=mu,lambda_0=lambda_0,sample_frac=sample_frac,Nsim=10000,nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda, method = "FFT", NtipMax = 200,conv=1e-3))
    
    if(! inherits(cP,"try-error")){
        dataFramePsi=rbind(dataFramePsi,data.frame(t=t,sigma=sigma,alpha=alpha,mu=mu,lambda_0=lambda_0,f=sample_frac,
                                                   computedPhi=cP$computedPhi,simulatedPhi=cP$simulatedPhi,
                                                   computedKhi=cP$computedKhi,simulatedKhi=cP$simulatedKhi,
                                                   nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda,seed=k,
                                                   Nsim=10000,conv=1e-3,
                                                   BD_Phi=cP$BD_Phi, BD_Khi=cP$BD_Khi))
      nspec[[length(nspec)+1]]=cP$n_spec
    }
}

# saved in /check_lik/MC_ctMu.Rdata

# # # # # # # # # # # # # # # # # # # # # # # # #
#             with constant turnover           #
# # # # # # # # # # # # # # # # # # # # # # # # #

check_Psi_ClaDS2=function(t,sigma,alpha,epsilon,lambda_0,sample_frac,Nsim,mlambda=0.0001,Mlambda=100,nlambda=1000,nt=100,method="Higham08.b",NtipMax=500,conv=1e-10){
  
  phi=Phi_ClaDS2(sigma = sigma,alpha = alpha,mlambda = mlambda,Mlambda = Mlambda,nlambda = nlambda,mu = epsilon,f = sample_frac,tini=0,tf=t,by=t/nt,method=method)
  khi=Khi_ClaDS2(phi = phi,s=0,t=t,func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$epsilon,
                 timePhi=phi$fun[,1],nt=1000,method="Magnus",conv=conv)
  simulateKhi=0
  n=0
  i=0
  n_spec=c()
  for (j in 1:Nsim){
    cat("\r",j,"\r")
    obj=sim_ClaDS(theta=1,lamb_par=lambda_0,mu_par=epsilon,lamb_shift=alpha,sigma=sigma,condition="time",new_lamb_law = "lognormal*shift",
                  new_mu_law="turnover",prune.extinct = T,time.stop = t,taxa.stop = NtipMax)
    tree=obj$tree
    n_spec=c(n_spec,length(obj$lamb))
    
    if(!is.null(tree)){
      if (length(tree)==2){ntip=1}else{ntip=tree$Nnode+1}
      
      if(ntip==NtipMax){
        simulateKhi=simulateKhi+0
        n=n+1
      }else{
        simulateKhi=simulateKhi+((1-sample_frac)^(ntip-1)*sample_frac*ntip)
        n=n+(1-(1-sample_frac)^ntip)
      }
      
      i=i+1
      # n=n+1
      
    }else{
      n=n+0
    }
  }
  simulatedKhi=simulateKhi/j
  simulatedPhi=1-(n/j)
  lind=which.min(abs(phi$expLambda-lambda_0))
  computedKhi=khi[lind]
  timePhi=phi$fun[,1]
  tind=which.min(abs(timePhi-t))
  computedPhi=phi$fun[tind,lind+1]
  mu=epsilon*lambda_0
  BD_Phi=1-(sample_frac*(lambda_0-mu)/(sample_frac*lambda_0+(lambda_0*(1-sample_frac)-mu)*exp(-t*(lambda_0-mu))))
  BD_Khi=sample_frac*((lambda_0-mu)^2)*exp(-t*(lambda_0-mu))/((sample_frac*lambda_0+(lambda_0*(1-sample_frac)-mu)*exp(-t*(lambda_0-mu)))^2)
  return(list(simulatedKhi=simulatedKhi,computedKhi=computedKhi,
              simulatedPhi=simulatedPhi,computedPhi=computedPhi,
              BD_Phi=BD_Phi, BD_Khi=BD_Khi,n_spec=n_spec))
}

dataFramePsi=data.frame()
nspec=list()

for(k in 1:500){
  set.seed(k)
  sigma=runif(1,0,1)
  alpha=runif(1,0.3,1.3)
  lambda_0=runif(1,0,0.5)
  mu=runif(1,0,1.5)
  sample_frac=runif(1,0,1)
  t=rexp(1,0.1)
  nlambda=300
  mlambda=1e-5
  Mlambda=100
  nt=30
  
  cat(" l0=",lambda_0,", s=",sigma,", mu=", mu,", t=",t,", f=",sample_frac, ", al=", alpha, "\n",sep="")
  cP=try(check_Psi_ClaDS1(t,sigma,alpha,epsilon=mu,lambda_0=lambda_0,sample_frac=sample_frac,Nsim=10000,nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda, method = "FFT", NtipMax = 200,conv=1e-3))
  
  if(! inherits(cP,"try-error")){
    dataFramePsi=rbind(dataFramePsi,data.frame(t=t,sigma=sigma,alpha=alpha,mu=mu,lambda_0=lambda_0,f=sample_frac,
                                               computedPhi=cP$computedPhi,simulatedPhi=cP$simulatedPhi,
                                               computedKhi=cP$computedKhi,simulatedKhi=cP$simulatedKhi,
                                               nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda,seed=k,
                                               Nsim=10000,conv=1e-3,
                                               BD_Phi=cP$BD_Phi, BD_Khi=cP$BD_Khi))
    nspec[[length(nspec)+1]]=cP$n_spec
  }
}

# saved in /check_lik/MC_turnover.Rdata

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#     with constant speciation and lognormal extinction     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

check_Psi_ClaDS3=function(t,sigma,alpha,epsilon,lambda_0,sample_frac,Nsim,mlambda=0.0001,Mlambda=100,nlambda=1000,nt=100,method="Higham08.b",NtipMax=500,conv=1e-10){
  
  phi=Phi_ClaDS3(sigma = sigma,alpha = alpha,mlambda = mlambda,Mlambda = Mlambda,nlambda = nlambda,mu = lambda_0,f = sample_frac,tini=0,tf=t,by=t/nt,method=method)
  khi=Khi_ClaDS3(phi = phi,s=0,t=t,func="Khi",lambda1=0,lambda2=0,lambdas=phi$lambda,M=phi$M,mu=phi$epsilon,
                 timePhi=phi$fun[,1],nt=1000,method="Magnus",conv=conv)
  simulateKhi=0
  n=0
  i=0
  n_spec=c()
  for (j in 1:Nsim){
    cat("\r",j,"\r")
    obj=sim_ClaDS(theta=1,lamb_par=lambda_0,mu_par=epsilon,mu_shift=alpha,sigma_mu=sigma,condition="time",new_lamb_law = "uniform",
                  new_mu_law="lognormal*shift",prune.extinct = T,time.stop = t,taxa.stop = NtipMax,lamb_max = lambda_0, lamb_min = lambda_0)
    tree=obj$tree
    n_spec=c(n_spec,length(obj$lamb))
    
    if(!is.null(tree)){
      if (length(tree)==2){ntip=1}else{ntip=tree$Nnode+1}
      
      if(ntip==NtipMax){
        simulateKhi=simulateKhi+0
        n=n+1
      }else{
        simulateKhi=simulateKhi+((1-sample_frac)^(ntip-1)*sample_frac*ntip)
        n=n+(1-(1-sample_frac)^ntip)
      }
      
      i=i+1
      # n=n+1
      
    }else{
      n=n+0
    }
  }
  simulatedKhi=simulateKhi/j
  simulatedPhi=1-(n/j)
  lind=which.min(abs(phi$expLambda-epsilon))
  computedKhi=khi[lind]
  timePhi=phi$fun[,1]
  tind=which.min(abs(timePhi-t))
  computedPhi=phi$fun[tind,lind+1]
  mu=epsilon
  BD_Phi=1-(sample_frac*(lambda_0-mu)/(sample_frac*lambda_0+(lambda_0*(1-sample_frac)-mu)*exp(-t*(lambda_0-mu))))
  BD_Khi=sample_frac*((lambda_0-mu)^2)*exp(-t*(lambda_0-mu))/((sample_frac*lambda_0+(lambda_0*(1-sample_frac)-mu)*exp(-t*(lambda_0-mu)))^2)
  return(list(simulatedKhi=simulatedKhi,computedKhi=computedKhi,
              simulatedPhi=simulatedPhi,computedPhi=computedPhi,
              BD_Phi=BD_Phi, BD_Khi=BD_Khi,n_spec=n_spec))
}

dataFramePsi=data.frame()
nspec=list()

for(k in 1:500){
  set.seed(k)
  sigma=runif(1,0,1)
  alpha=runif(1,0.3,1.6)
  lambda_0=runif(1,0,0.5)
  mu=runif(1,0,1.5)*lambda_0
  sample_frac=runif(1,0,1)
  t=rexp(1,0.1)
  nlambda=300
  mlambda=1e-5
  Mlambda=1000
  nt=100
  
  cat(" l0=",lambda_0,", s=",sigma,", mu=", mu,", t=",t,", f=",sample_frac, ", al=", alpha, "\n",sep="")
  cP=try(check_Psi_ClaDS1(t,sigma,alpha,epsilon=mu,lambda_0=lambda_0,sample_frac=sample_frac,Nsim=10000,nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda, method = "FFT", NtipMax = 200,conv=1e-3))
  
  if(! inherits(cP,"try-error")){
    dataFramePsi=rbind(dataFramePsi,data.frame(t=t,sigma=sigma,alpha=alpha,mu=mu,lambda_0=lambda_0,f=sample_frac,
                                               computedPhi=cP$computedPhi,simulatedPhi=cP$simulatedPhi,
                                               computedKhi=cP$computedKhi,simulatedKhi=cP$simulatedKhi,
                                               nt=nt,nlambda = nlambda,mlambda = mlambda,Mlambda = Mlambda,seed=k,
                                               Nsim=10000,conv=1e-3,
                                               BD_Phi=cP$BD_Phi, BD_Khi=cP$BD_Khi))
    nspec[[length(nspec)+1]]=cP$n_spec
  }
}

# saved in /check_lik/MC_ChangeDeath.Rdata

#### ClaDS tree shapes ####

# # # # # # # # # # # # # # # # # # # # # # # # #
#           with constant extinction            #
# # # # # # # # # # # # # # # # # # # # # # # # #

dataFrame=data.frame()
for(mean in seq(0.5,1.3,0.1)){
  for(sigma in rev(c(0,0.05,0.1,0.2,0.3,0.5,0.8,1,1.5,2))){
    for(mu in c(0.00,0.01,0.05,0.09)){
      print(c(mean,sigma,mu))
      print("")
      set.seed(replicate)
      for(replicate in 1:100){
        tree=NULL
        j=1
        while(is.null(tree) & j<1000){
          cat("\r","Turn ",replicate,";",j,"\r")
          a=sim_ClaDS(theta=1,lamb_par=0.1,
                      mu_par=mu,
                      lamb_shift=mean/exp((sigma^2)/2),
                      sigma=sigma,taxa.stop = 200,
                      condition="taxa",new_mu_law = "uniform",
                      new_lamb_law = "lognormal*shift",mu_min=mu,mu_max=mu,
                      prune.extinct = T)
          if(!length(a)<=2) tree=a$tree
          j=j+1
        }
        if(! is.null(tree)){
          true.rate=a$lamb[a$rates]
          dataFrame=rbind(dataFrame,data.frame(seed=replicate,gamma=gammaStat(tree),beta=maxlik.betasplit(tree)$max_lik,
                                               maxRate=max(true.rate),minRate=min(true.rate),l0=0.1,
                                               sigma=sigma,mean=mean,mu=mu,age=max(node.depth.edgelength(tree))))
        }
      }
    }
  }
}


# saved in /tree_shapes/treeShapes_ctMu.Rdata

# # # # # # # # # # # # # # # # # # # # # # # # #
#             with constant turnover            #
# # # # # # # # # # # # # # # # # # # # # # # # #

dataFrame=data.frame()
for(mean in seq(0.5,1.3,0.1)){
  for(sigma in rev(c(0,0.05,0.1,0.2,0.3,0.5,0.8,1,1.5,2))){
    for(mu in c(0,0.1,0.5,0.9)){
      print(c(mean,sigma,mu))
      print("")
      for(replicate in 1:100){
        tree=NULL
        j=1
        while(is.null(tree) & j<1000){
          cat("\r","Turn ",replicate,";",j,"\r")
          a=sim_ClaDS(theta=1,lamb_par=0.1,
                                      mu_par=mu,
                                      lamb_shift=mean/exp((sigma^2)/2),
                                      sigma=sigma,taxa.stop = 200,
                                      condition="taxa",new_mu_law = "turnover",
                                      new_lamb_law = "lognormal*shift",
                                      prune.extinct = T)
          if(!length(a)<=2) tree=a$tree
          j=j+1
        }
        if(! is.null(tree)){
          true.rate=a$lamb[a$rates]
          dataFrame=rbind(dataFrame,data.frame(seed=replicate,gamma=gammaStat(tree),beta=maxlik.betasplit(tree)$max_lik,
                                             maxRate=max(true.rate),minRate=min(true.rate),l0=0.1,
                                             sigma=sigma,mean=mean,mu=mu,age=max(node.depth.edgelength(tree))))
        }
      }
    }
  }
}

# saved in /tree_shapes/treeShapes_turnover.Rdata

#### test of ClaDS0 ####

# # # # # # # # # # # # # # # # # # # # # # # # #
#                on ClaDS0 trees                #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=paste0(rep(paste0(rep(paste0("LNormMH_mcmcDE",c(50,100,200)),each=6),
             c("_SP_alpha12_sigma","_SP_alpha11_sigma",
               "_SP_sigma","_SP_alpha95_sigma",
               "_SP_alpha9_sigma","_SP_alpha7_sigma")),each=6),0:5)

Ntips=rep(c(50,100,200),each=36)
Sigmas=rep(rep(log(1+(0:5)*0.1),6),3)
Alphas=rep(rep(c(1.2,1.1,1,0.95,0.9,0.7),each=6),3)

id=1 # in 1:108
i=1 # in 1:20

name=Names[id]
ntip=Ntips[id]
sigma=Sigmas[id]
alpha=Alphas[id]
    
obj= sim_ClaDS( lamb_par=0.1,      # initial speciation rate
                mu_par=0,          # turnover rate (extinction/speciation)
                sigma=sigma,           # standard deviation of the new rates law
                lamb_shift=alpha,      # trend parameter alpha
                condition="taxa",    # the stoppping condition to use (can also be time)
                taxa.stop = ntip,      # the number of tips to simulate (if condition=="taxa")
                prune.extinct = T)   # should extincti taxa be pruned from the result? (default to T)
    
    
# the function returns a list with the tree and the associated rates. Here is how you get to them :
tree = obj$tree
true.rates = obj$lamb[obj$rates]
  
sampler = run_ClaDS0(tree=tree,        
                    name=paste0(name,i),    
                    nCPU=1,    
                    pamhLocalName = paste0("local_",name,"_",i),  
                    iteration=500000, 
                    thin=2000, 
                    update=1000, adaptation=5000)
    
save(tree,true.rate,sampler,file = paste(name,i))


# # # # # # # # # # # # # # # # # # # # # # # # #
#              on one shift trees               #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=c("LNMHShift_d3","LNMHShift_d2","LNMHShift_d1",
        "LNMHShift_x1","LNMHShift_x1-5","LNMHShift_x2","LNMHShift_x3","LNMHShift_x4")

lamb_range=c(0.025,0.03,0.05,0.1,0.15,0.2,0.3,0.4,1)

id=1 # in 1:8
i = 1 # in 1:100

name=Names[id]

obj= sim_ClaDS( lamb_par=0.1, 
                mu_par=0,        
                lamb_min = lamb_range[id], lamb_max = lamb_range[id+1],
                condition="taxa",  
                taxa.stop = 200,   
                prune.extinct = T)  
    
    
# the function returns a list with the tree and the associated rates. Here is how you get to them :
tree = obj$tree
true.rates = obj$lamb[obj$rates]
    
sampler = run_ClaDS0(tree=tree,        
                     name=paste0(name,i),    
                     nCPU=1,    
                     pamhLocalName = paste0("local_",name,"_",i),  
                     iteration=500000, 
                     thin=2000, 
                     update=1000, adaptation=5000)
    
save(tree,true.rate,sampler,file = paste(name,i))


# # # # # # # # # # # # # # # # # # # # # # # # #
#   on trees with logbrownian rate variations   #
# # # # # # # # # # # # # # # # # # # # # # # # #

id=1 # in 1:3
replicate=1 # in 1:10
ntips=50 # in c(50,100,200)

obj= sim_ClaDS( lamb_par=0.1,
                mu_par=0,    
                sigma=0.05*id,      
                lamb_shift=1,    
                condition="taxa",  
                taxa.stop = ntips,  
                prune.extinct = T, 
                new_lamb_law = "logbrownian",
                return.all.extinct =   F)  

tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

save.image(paste0("LBtree",ntips,"_",id,"_",replicate,".Rdata"))

sampler = run_ClaDS0(tree=tree,       
                     name=paste0("Cl0_LBtree",ntips,"_",id,"_",replicate,".Rdata"),
                     nCPU=1,               
                     pamhLocalName = paste0("local",ntips,"_",id,"_",replicate),
                     iteration=500000,
                     thin=2000,
                     update=1000, adaptation=5000)


# # # # # # # # # # # # # # # # # # # # # # # # #
#    on trees with evolving extinction rate     #
# # # # # # # # # # # # # # # # # # # # # # # # #

id=1 # in 1:3
replicate=5 # in 1:5
ntips=50 # in c(50,100)

obj= sim_ClaDS( lamb_par=0.1, lamb_max = 0.1, lamb_min = 0.1, # parameters for the uniform draw (here constant speciation rate)
                mu_par=0.05,
                condition="taxa", 
                taxa.stop = ntips,     
                new_lamb_law = "uniform",
                new_mu_law="lognormal*shift",
                mu_shift=1,sigma_mu=0.2*id,
                return.all.extinct =   F)   

tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]
save.image(paste0("CDtree",ntips,"_",id,"_",replicate,".Rdata"))

sampler = run_ClaDS0(tree=tree,       
                     name=paste0("Cl0_CDtree",ntips,"_",id,"_",replicate,".Rdata"),
                     nCPU=1,             
                     pamhLocalName = paste0("local",ntips,"_",id,"_",replicate),   
                     iteration=500000, 
                     thin=2000,    
                     update=1000, adaptation=5000) 


#### comparison with other methods ####

# # # # # # # # # # # # # # # # # # # # # # # # #
#        computation of the DR statistics       #
# # # # # # # # # # # # # # # # # # # # # # # # #

DR=function(tree){
  
  dr_stat=function(tree,tip){
    root=tree$Nnode+2
    internal=function(node){
      if(node==root){
        return(0)
      }else{
        edge=which(tree$edge[,2]==node)
        return(internal(tree$edge[edge,1])/2+tree$edge.length[edge])
      }
    }
    
    1/internal(tip)
  }
  
  sapply(1:(tree$Nnode+1),function(i){dr_stat(tree,i)})
}

# # # # # # # # # # # # # # # # # # # # # # # # #
#         preparing BAMM runs _ function        #
# # # # # # # # # # # # # # # # # # # # # # # # #

writeTemplate_withBAMMprior=function(tree,outfile,treeName){
  fileName <- 'template.txt'
  template=readChar(fileName, file.info(fileName)$size)
  prior=setBAMMpriors(tree, total.taxa = NULL, traits = NULL,
                      outfile = NULL,
                      Nmax = 1000, suppressWarning = FALSE)
  prior_char=paste0(names(prior)[1]," = ",prior[1], "\n\n",
                    names(prior)[2]," = ",prior[2], "\n\n",
                    names(prior)[3]," = ",prior[3], "\n\n",
                    names(prior)[4]," = ",prior[4], "\n\n")
  text=paste0("treefile = ",treeName,".tre \n\n outName =",treeName,"\n\n",prior_char, "\n\n", template)
  write(text, file = outfile, sep = "")
}

# then run from terminal using command    ./bamm -c outfile 


# # # # # # # # # # # # # # # # # # # # # # # # #
#               for ClaDS0 trees                #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=c(c(paste0("LNormMH_mcmcDE200_SP_alpha12_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_alpha11_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_alpha95_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_alpha9_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_alpha7_sigma",c(0:5))),
        c(paste0("LNormMH_mcmcDE100_SP_alpha12_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_alpha11_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_alpha95_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_alpha9_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_alpha7_sigma",c(0:5))),
        c(paste0("LNormMH_mcmcDE50_SP_alpha12_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_alpha11_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_alpha95_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_alpha9_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_alpha7_sigma",c(0:5))))

for(name in Names){
  for(id in 1:20){
    t=try(load(paste0(name," ",id)))
    if(! inherits(t, "try-error")){
      write.tree(tree,file=paste0(name,"_",id,".tre"))
      writeTemplate_withBAMMprior(tree, treeName = paste0(name,"_",id), 
                                  outfile = paste0("template_",name,"_",id,".txt"))
    }
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # #
#             for one shift trees               #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=(paste0("LNormMHShift_",c(paste0("d",1:3),paste0("x",c(1,"1-5",2:4)))))

for(name in Names){
  for(id in 1:50){
    t=try(load(paste0(name," ",id)))
    if(! inherits(t, "try-error")){
      write.tree(tree,file=paste0(name,"_",id,".tre"))
      writeTemplate_withBAMMprior(tree, treeName = paste0(name,"_",id), 
                                  outfile = paste0("template_",name,"_",id,".txt"))
    }
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # #
#  for trees with logbrownian rate variations   #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=c(paste0("Cl0_LBtree50_",1:3),paste0("Cl0_LBtree100_",1:3),paste0("Cl0_LBtree200_",1:3))

for(name in Names){
  for(id in 1:10){
    t=try(load(paste0(name,"_",id,".Rdata")))
    if(! inherits(t, "try-error")){
      write.tree(tree,file=paste0(name,"_",id,".tre"))
      writeTemplate_withBAMMprior(tree, treeName = paste0(name,"_",id), 
                                  outfile = paste0("template_",name,"_",id,".txt"))
    }
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # #
#    on trees with evolving extinction rate     #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=c(paste0("Cl0_CDtree50_",1:3),paste0("Cl0_CDtree100_",1:3),paste0("Cl0_CDtree200_",1:3))

for(name in Names){
  for(id in 1:10){
    t=try(load(paste0(name,"_",id,".Rdata")))
    if(! inherits(t, "try-error")){
      write.tree(tree,file=paste0(name,"_",id,".tre"))
      writeTemplate_withBAMMprior(tree, treeName = paste0(name,"_",id), 
                                  outfile = paste0("template_",name,"_",id,".txt"))
    }
  }
}

#### extract the results from the MCMC chains ####
colorPal=colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 

# # # # # # # # # # # # # # # # # # # # # # # # #
#                on ClaDS0 trees                #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=c(c(paste0("LNormMH_mcmcDE200_SP_alpha12_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_alpha11_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_alpha95_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_alpha9_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE200_SP_alpha7_sigma",c(0:5))),
        c(paste0("LNormMH_mcmcDE100_SP_alpha12_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_alpha11_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_alpha95_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_alpha9_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE100_SP_alpha7_sigma",c(0:5))),
        c(paste0("LNormMH_mcmcDE50_SP_alpha12_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_alpha11_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_alpha95_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_alpha9_sigma",c(0:5)),
          paste0("LNormMH_mcmcDE50_SP_alpha7_sigma",c(0:5))))



data_BAMM=data.frame()


if(T){
  Sigmas=rep(rep(log(1+0.1*rep(0:5,6))),3)
  Alphas=rep(rep(c(1.2,1.1,1,0.95,0.9,0.7), each=6),3)
  Ntips=rep(c(200,100,50),each=36)
  for(id_tree in 1:20){
    name_id=0
    
    for(name in Names){
      name_id=name_id+1
      
      do=F
      ntip=Ntips[name_id]
      nedges=ntip*2-2
      do=(sum(data_BAMM$name==name & data_BAMM$seed==id_tree)==0)
      if(! do) do = (data_BAMM$ngen[data_BAMM$name==name & data_BAMM$seed==id_tree]<801)
      if(do & file.exists(paste0(name,"_",id_tree,"_event_data.txt"))){
        tree=read.tree(paste0(name,"_",id_tree,".tre"))
        ntip=tree$Nnode+1
        edata=try(getEventData(tree, eventdata = paste0(name,"_",id_tree,"_event_data.txt"), burnin=0.2 ))
        if(! inherits(edata,"try-error")){
          
          load(paste0(name," ",id_tree))
          
          tips_edges=sapply(1:(tree$Nnode+1),function(i){which(tree$edge[,2]==i)})
          
          tot_time<-max(node.age(tree)$ages)
          f.lamb <-function(t,y){y[1]}
          f.mu<-function(t,y){y[1]}
          lamb_par<-c(0.1)
          mu_par<-c(0.001)
          result_cst <- fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par, f=1,cst.lamb=TRUE,cst.mu = T,fix.mu=F,dt=1e-3)
          
          BD_rates=c(abs(result_cst$lamb_par),abs(result_cst$mu_par))
          MSE_BD=mean((BD_rates[1]-true.rate)^2)
          MSE_log_BD=mean((log(BD_rates[1])-log(true.rate))^2)
          
          #### get ClaDS0 results
          rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),])}))
          rep2=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),-c((2*ntip+2):(2*ntip+4))])}))
          
          meanRate = rep(0,nedges+3)
          for(j in 1:3){
            meanRate=meanRate+colSums(cbind((rep[[j]][,1]),rep[[j]][,2],(rep[[j]][,3]),(rep[[j]][,2]+rep[[j]][,4:(nedges+3)])))
          }
          meanRate=meanRate/(3*nrow(rep[[1]]))
          rate=exp(get_rates(tree,meanRate[-c(1:2)]))
          sigma_inf=meanRate[1]
          alpha_inf=exp(meanRate[2])
          lambda_0_inf=rate[1]
          ClaDS0_rate=rate[-1]
          
          #### DR's
          DR_rates=DR(tree)
          tip_rates=true.rate[tips_edges]
          ClaDS_tip_rates=ClaDS0_rate[tips_edges]
          
          #### and BAMM's
          summary(edata)
          
          mcmcout <- read.csv(paste0(name,"_",id_tree,"_mcmc_out.txt"), header=T)
          burnstart <- floor(0.1 * nrow(mcmcout))
          postburn <- mcmcout[burnstart:nrow(mcmcout), ]
          ES_Nshifts=effectiveSize(postburn$N_shifts)
          ES_LP=effectiveSize(postburn$logLik)
          
          BAMM_rates=getMeanBranchLengthTree(edata)$phy$edge.length
          BAMM_div_rates=getMeanBranchLengthTree(edata,rate = "ndr")$phy$edge.length
          BAMM_ext_rates=getMeanBranchLengthTree(edata,rate = "extinction")$phy$edge.length
          
          lm_BAMM=try(lm(BAMM_rates~true.rate))
          
          lm_log_BAMM=try(lm(log(BAMM_rates)~log(true.rate)))
          if(inherits(lm_log_BAMM, "try-error")){
            lm_log_BAMM=list(coefficients=rep(NA,2))
          }
          
          lm_BAMM_div=try(lm(BAMM_div_rates~true.rate))
          if(inherits(lm_BAMM_div, "try-error")){
            lm_BAMM_div=list(coefficients=rep(NA,2))
          }
          lm_log_BAMM_div=try(lm(log(BAMM_div_rates)~log(true.rate)))
          if(inherits(lm_log_BAMM_div, "try-error")){
            lm_log_BAMM_div=list(coefficients=rep(NA,2))
          }
          lm_ClaDS=try(lm(ClaDS0_rate~true.rate))
          if(inherits(lm_ClaDS, "try-error")){
            lm_ClaDS=list(coefficients=rep(NA,2))
          }
          lm_log_ClaDS=try(lm(log(ClaDS0_rate)~log(true.rate)))
          if(inherits(lm_log_ClaDS, "try-error")){
            lm_log_ClaDS=list(coefficients=rep(NA,2))
          }
          
          lm_ClaDS_tips=try(lm(ClaDS_tip_rates~tip_rates))
          if(inherits(lm_ClaDS_tips, "try-error")){
            lm_ClaDS_tips=list(coefficients=rep(NA,2))
          }
          lm_log_ClaDS_tips=try(lm(log(ClaDS_tip_rates)~log(tip_rates)))
          if(inherits(lm_log_ClaDS_tips, "try-error")){
            lm_log_ClaDS_tips=list(coefficients=rep(NA,2))
          }
          lm_DR=try(lm(DR_rates~tip_rates))
          if(inherits(lm_DR, "try-error")){
            lm_DR=list(coefficients=rep(NA,2))
          }
          lm_log_DR=try(lm(log(DR_rates)~log(tip_rates)))
          if(inherits(lm_log_DR, "try-error")){
            lm_log_DR=list(coefficients=rep(NA,2))
          }
          
          if(inherits(lm_BAMM, "try-error")){
            lm_BAMM=list(coefficients=rep(NA,2))
          }
          
          MSE_BAMM=mean((true.rate-BAMM_rates)^2)
          MSE_log_BAMM=mean((log(true.rate)-log(BAMM_rates))^2)
          MSE_BAMM_div=mean((true.rate-BAMM_div_rates)^2)
          MSE_log_BAMM_div=mean((log(true.rate)-log(BAMM_div_rates))^2)
          MSE_ClaDS=mean((true.rate-ClaDS0_rate)^2)
          MSE_log_ClaDS=mean((log(true.rate)-log(ClaDS0_rate))^2)
          MSE_ClaDS_tips=mean((tip_rates-ClaDS_tip_rates)^2)
          MSE_log_ClaDS_tips=mean((log(tip_rates)-log(ClaDS_tip_rates))^2)
          MSE_DR=mean((tip_rates-DR_rates)^2)
          MSE_log_DR=mean((log(tip_rates)-log(DR_rates))^2)
          
          var_true=var(true.rate)
          var_log_true=var(log(true.rate))
          var_true_tips=var(tip_rates)
          var_log_true_tips=var(log(tip_rates))
          
          var_BAMM=var(BAMM_rates)
          var_log_BAMM=var(log(BAMM_rates))
          var_BAMM_div=var(BAMM_div_rates)
          var_log_BAMM_div=var(log(BAMM_div_rates))
          var_ClaDS=var(ClaDS0_rate)
          var_log_ClaDS=var(log(ClaDS0_rate))
          var_ClaDS_tips=var(ClaDS_tip_rates)
          var_log_ClaDS_tips=var(log(ClaDS_tip_rates))
          var_DR=var(DR_rates)
          var_log_DR=var(log(DR_rates))
          
          do=(sum(data_BAMM$name==name & data_BAMM$seed==id_tree)==0)
          line_df=data.frame(name=name,seed=id_tree,ES_Nshifts=ES_Nshifts,ES_LP=ES_LP,ntip=Ntips[name_id],ngen=length(edata$numberEvents),
                             sigma_inf=sigma_inf,alpha_inf=alpha_inf,mean_new_rate=exp(log(alpha_inf)+(sigma_inf^2)/2),
                             var_log_true_tips=var_log_true_tips,var_true_tips=var_true_tips,
                             var_true=var_true,var_log_true=var_log_true,
                             sigma=Sigmas[name_id],shift=Alphas[name_id],ntip=ntip,
                             cor_BAMM=cor(BAMM_rates,true.rate),slope_BAMM=lm_BAMM$coefficients[2],
                             cor_log_BAMM=cor(log(BAMM_rates),log(true.rate)),slope_log_BAMM=lm_log_BAMM$coefficients[2],
                             MSE_BAMM=MSE_BAMM,MSE_log_BAMM=MSE_log_BAMM,
                             rel_error_BAMM=exp(mean(log(BAMM_rates)-log(true.rate))),
                             var_BAMM=var_BAMM,var_log_BAMM=var_log_BAMM,BAMM_ext=mean(BAMM_ext_rates),
                             cor_BAMM_div=cor(BAMM_div_rates,true.rate),slope_BAMM_div=lm_BAMM_div$coefficients[2],
                             cor_log_BAMM_div=cor(log(BAMM_div_rates),log(true.rate)),slope_log_BAMM_div=lm_log_BAMM_div$coefficients[2],
                             MSE_BAMM_div=MSE_BAMM_div,MSE_log_BAMM_div=MSE_log_BAMM_div,
                             rel_error_BAMM_div=exp(mean(log(BAMM_div_rates)-log(true.rate))),
                             var_BAMM_div=var_BAMM_div,var_log_BAMM_div=var_log_BAMM_div,
                             var_log_ClaDS=var_log_ClaDS,var_ClaDS=var_ClaDS,
                             MSE_log_BD=MSE_log_BD,MSE_BD=MSE_BD,BD_birth=BD_rates[1],BD_death=BD_rates[2],rel_error_BD=exp(sum(log(BD_rates[1]/true.rate))/nedges),
                             cor_ClaDS=cor(ClaDS0_rate,true.rate),slope_ClaDS=lm_ClaDS$coefficients[2],
                             cor_log_ClaDS=cor(log(ClaDS0_rate),log(true.rate)),slope_log_ClaDS=lm_log_ClaDS$coefficients[2],
                             MSE_ClaDS=MSE_ClaDS,MSE_log_ClaDS=MSE_log_ClaDS,
                             rel_error_ClaDS=exp(mean(log(ClaDS0_rate)-log(true.rate))),
                             cor_DR=cor(DR_rates,tip_rates),slope_DR=lm_DR$coefficients[2],
                             cor_log_DR=cor(log(DR_rates),log(tip_rates)),slope_log_DR=lm_log_DR$coefficients[2],
                             MSE_DR=MSE_DR,MSE_log_DR=MSE_log_DR,
                             rel_error_DR=exp(mean(log(DR_rates)-log(tip_rates))),
                             var_DR=var_DR,var_log_DR=var_log_DR,
                             cor_ClaDS_tips=cor(ClaDS_tip_rates,tip_rates),slope_ClaDS_tips=lm_ClaDS_tips$coefficients[2],
                             cor_log_ClaDS_tips=cor(log(ClaDS_tip_rates),log(tip_rates)),slope_log_ClaDS_tips=lm_log_ClaDS_tips$coefficients[2],
                             MSE_ClaDS_tips=MSE_ClaDS_tips,MSE_log_ClaDS_tips=MSE_log_ClaDS_tips,
                             rel_error_ClaDS_tips=exp(mean(log(ClaDS_tip_rates)-log(tip_rates))),
                             var_ClaDS_tips=var_ClaDS_tips,var_log_ClaDS_tips=var_log_ClaDS_tips)
          
          if(!do){
            data_BAMM[data_BAMM$name==name & data_BAMM$seed==id_tree,]=line_df
            
          }else{
            data_BAMM=rbind(data_BAMM,line_df)
          }
          
        }
        save(data_BAMM,file="ClaDS0_BAMM.Rdata")
      }
    }
  }
}



# # # # # # # # # # # # # # # # # # # # # # # # #
#              on one shift trees               #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=(paste0("LNormMHShift_",c(paste0("d",1:3),paste0("x",c(1,"1-5",2:4)))))

data_OneShift=data.frame()


if(T){
  name_id=0
  ntip=200
  nedges=2*ntip-2
  for(id_tree in 1:50){
    for(name in Names){
      name_id=name_id+1
      do=F
      do=(sum(data_OneShift$name==name & data_OneShift$seed==id_tree)==0)
      if(! do) do = (data_OneShift$ngen[data_OneShift$name==name & data_OneShift$seed==id_tree]<801)
      if(do & file.exists(paste0(name,"_",id_tree,"_event_data.txt"))){
        tree=read.tree(paste0(name,"_",id_tree,".tre"))
        
        edata=try(getEventData(tree, eventdata = paste0(name,"_",id_tree,"_event_data.txt"), burnin=0.2 ))
        if(! inherits(edata,"try-error")){
          
          load(paste0(name," ",id_tree))
          tips_edges=sapply(1:(tree$Nnode+1),function(i){which(tree$edge[,2]==i)})
          
          tot_time<-max(node.age(tree)$ages)
          f.lamb <-function(t,y){y[1]}
          f.mu<-function(t,y){y[1]}
          lamb_par<-c(1)
          mu_par<-c(0.1)
          result_cst <- fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par, f=1,cst.lamb=TRUE,cst.mu = T,fix.mu=F,dt=1e-3)
          
          BD_rates=c(abs(result_cst$lamb_par),abs(result_cst$mu_par))
          MSE_BD=mean((BD_rates[1]-true.rate)^2)
          MSE_log_BD=mean((log(BD_rates[1])-log(true.rate))^2)
          
          #### get ClaDS0 results
          rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),])}))
          rep2=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),-c((2*ntip+2):(2*ntip+4))])}))
          
          meanRate = rep(0,nedges+3)
          for(j in 1:3){
            meanRate=meanRate+colSums(cbind((rep[[j]][,1]),rep[[j]][,2],(rep[[j]][,3]),(rep[[j]][,2]+rep[[j]][,4:(nedges+3)])))
          }
          meanRate=exp(meanRate/(3*nrow(rep[[1]])))
          rate=exp(get_rates(tree,log(meanRate[-c(1:2)])))
          lambda_0_inf=rate[1]
          ClaDS0_rate=rate[-1]
          
          #### DR's
          DR_rates=DR(tree)
          tip_rates=true.rate[tips_edges]
          ClaDS_tip_rates=ClaDS0_rate[tips_edges]
          
          #### and BAMM's
          summary(edata)

          mcmcout <- read.csv(paste0(name,"_",id_tree,"_mcmc_out.txt"), header=T)
          burnstart <- floor(0.1 * nrow(mcmcout))
          postburn <- mcmcout[burnstart:nrow(mcmcout), ]
          ES_Nshifts=effectiveSize(postburn$N_shifts)
          ES_LP=effectiveSize(postburn$logLik)
          
          BAMM_rates=getMeanBranchLengthTree(edata)$phy$edge.length
          BAMM_div_rates=getMeanBranchLengthTree(edata,rate = "ndr")$phy$edge.length
          BAMM_ext_rates=getMeanBranchLengthTree(edata,rate = "extinction")$phy$edge.length
          
          lm_BAMM=try(lm(BAMM_rates~true.rate))
          
          lm_log_BAMM=try(lm(log(BAMM_rates)~log(true.rate)))
          if(inherits(lm_log_BAMM, "try-error")){
            lm_log_BAMM=list(coefficients=rep(NA,2))
          }
          
          lm_BAMM_div=try(lm(BAMM_div_rates~true.rate))
          if(inherits(lm_BAMM_div, "try-error")){
            lm_BAMM_div=list(coefficients=rep(NA,2))
          }
          lm_log_BAMM_div=try(lm(log(BAMM_div_rates)~log(true.rate)))
          if(inherits(lm_log_BAMM_div, "try-error")){
            lm_log_BAMM_div=list(coefficients=rep(NA,2))
          }
          lm_ClaDS=try(lm(ClaDS0_rate~true.rate))
          if(inherits(lm_ClaDS, "try-error")){
            lm_ClaDS=list(coefficients=rep(NA,2))
          }
          lm_log_ClaDS=try(lm(log(ClaDS0_rate)~log(true.rate)))
          if(inherits(lm_log_ClaDS, "try-error")){
            lm_log_ClaDS=list(coefficients=rep(NA,2))
          }
          
          lm_ClaDS_tips=try(lm(ClaDS_tip_rates~tip_rates))
          if(inherits(lm_ClaDS_tips, "try-error")){
            lm_ClaDS_tips=list(coefficients=rep(NA,2))
          }
          lm_log_ClaDS_tips=try(lm(log(ClaDS_tip_rates)~log(tip_rates)))
          if(inherits(lm_log_ClaDS_tips, "try-error")){
            lm_log_ClaDS_tips=list(coefficients=rep(NA,2))
          }
          lm_DR=try(lm(DR_rates~tip_rates))
          if(inherits(lm_DR, "try-error")){
            lm_DR=list(coefficients=rep(NA,2))
          }
          lm_log_DR=try(lm(log(DR_rates)~log(tip_rates)))
          if(inherits(lm_log_DR, "try-error")){
            lm_log_DR=list(coefficients=rep(NA,2))
          }

          if(inherits(lm_BAMM, "try-error")){
            lm_BAMM=list(coefficients=rep(NA,2))
          }
          
          MSE_BAMM=mean((true.rate-BAMM_rates)^2)
          MSE_log_BAMM=mean((log(true.rate)-log(BAMM_rates))^2)
          MSE_BAMM_div=mean((true.rate-BAMM_div_rates)^2)
          MSE_log_BAMM_div=mean((log(true.rate)-log(BAMM_div_rates))^2)
          MSE_ClaDS=mean((true.rate-ClaDS0_rate)^2)
          MSE_log_ClaDS=mean((log(true.rate)-log(ClaDS0_rate))^2)
          MSE_ClaDS_tips=mean((tip_rates-ClaDS_tip_rates)^2)
          MSE_log_ClaDS_tips=mean((log(tip_rates)-log(ClaDS_tip_rates))^2)
          MSE_DR=mean((tip_rates-DR_rates)^2)
          MSE_log_DR=mean((log(tip_rates)-log(DR_rates))^2)
          
          var_true=var(true.rate)
          var_log_true=var(log(true.rate))
          var_true_tips=var(tip_rates)
          var_log_true_tips=var(log(tip_rates))
          
          var_BAMM=var(BAMM_rates)
          var_log_BAMM=var(log(BAMM_rates))
          var_BAMM_div=var(BAMM_div_rates)
          var_log_BAMM_div=var(log(BAMM_div_rates))
          var_ClaDS=var(ClaDS0_rate)
          var_log_ClaDS=var(log(ClaDS0_rate))
          var_ClaDS_tips=var(ClaDS_tip_rates)
          var_log_ClaDS_tips=var(log(ClaDS_tip_rates))
          var_DR=var(DR_rates)
          var_log_DR=var(log(DR_rates))
          
          do=(sum(data_OneShift$name==name & data_OneShift$seed==id_tree)==0)
          
          n_new_branches=sum(true.rate!=0.1)
          new_rate=true.rate[true.rate!=0.1][1]
          line_df=data.frame(name=name,seed=id_tree,ES_Nshifts=ES_Nshifts,ES_LP=ES_LP,ntip=200,ngen=length(edata$numberEvents),
                             n_new_branches=n_new_branches,new_rate=new_rate,
                             var_true=var_true,var_log_true=var_log_true,var_true_tips=var_true_tips,var_log_true_tips=var_log_true_tips,
                             cor_BAMM=cor(BAMM_rates,true.rate),slope_BAMM=lm_BAMM$coefficients[2],
                             cor_log_BAMM=cor(log(BAMM_rates),log(true.rate)),slope_log_BAMM=lm_log_BAMM$coefficients[2],
                             MSE_BAMM=MSE_BAMM,MSE_log_BAMM=MSE_log_BAMM,
                             rel_error_BAMM=exp(mean(log(BAMM_rates)-log(true.rate))),
                             var_BAMM=var_BAMM,var_log_BAMM=var_log_BAMM,BAMM_ext=mean(BAMM_ext_rates),
                             cor_BAMM_div=cor(BAMM_div_rates,true.rate),slope_BAMM_div=lm_BAMM_div$coefficients[2],
                             cor_log_BAMM_div=cor(log(BAMM_div_rates),log(true.rate)),slope_log_BAMM_div=lm_log_BAMM_div$coefficients[2],
                             MSE_BAMM_div=MSE_BAMM_div,MSE_log_BAMM_div=MSE_log_BAMM_div,
                             rel_error_BAMM_div=exp(mean(log(BAMM_div_rates)-log(true.rate))),
                             var_BAMM_div=var_BAMM_div,var_log_BAMM_div=var_log_BAMM_div,
                             var_log_ClaDS=var_log_ClaDS,var_ClaDS=var_ClaDS,
                             MSE_log_BD=MSE_log_BD,MSE_BD=MSE_BD,BD_birth=BD_rates[1],BD_death=BD_rates[2],rel_error_BD=exp(sum(log(BD_rates[1]/true.rate))/nedges),
                             cor_ClaDS=cor(ClaDS0_rate,true.rate),slope_ClaDS=lm_ClaDS$coefficients[2],
                             cor_log_ClaDS=cor(log(ClaDS0_rate),log(true.rate)),slope_log_ClaDS=lm_log_ClaDS$coefficients[2],
                             MSE_ClaDS=MSE_ClaDS,MSE_log_ClaDS=MSE_log_ClaDS,
                             rel_error_ClaDS=exp(mean(log(ClaDS0_rate)-log(true.rate))),
                             cor_DR=cor(DR_rates,tip_rates),slope_DR=lm_DR$coefficients[2],
                             cor_log_DR=cor(log(DR_rates),log(tip_rates)),slope_log_DR=lm_log_DR$coefficients[2],
                             MSE_DR=MSE_DR,MSE_log_DR=MSE_log_DR,
                             rel_error_DR=exp(mean(log(DR_rates)-log(tip_rates))),
                             var_DR=var_DR,var_log_DR=var_log_DR,
                             cor_ClaDS_tips=cor(ClaDS_tip_rates,tip_rates),slope_ClaDS_tips=lm_ClaDS_tips$coefficients[2],
                             cor_log_ClaDS_tips=cor(log(ClaDS_tip_rates),log(tip_rates)),slope_log_ClaDS_tips=lm_log_ClaDS_tips$coefficients[2],
                             MSE_ClaDS_tips=MSE_ClaDS_tips,MSE_log_ClaDS_tips=MSE_log_ClaDS_tips,
                             rel_error_ClaDS_tips=exp(mean(log(ClaDS_tip_rates)-log(tip_rates))),
                             var_ClaDS_tips=var_ClaDS_tips,var_log_ClaDS_tips=var_log_ClaDS_tips)
          if(!do){
            data_OneShift[data_OneShift$name==name & data_OneShift$seed==id_tree,]=line_df
            
          }else{
            data_OneShift=rbind(data_OneShift,line_df)
          }
          
        }}
      save(data_OneShift,file="./BAMM_OS.Rdata")
    }
  }
}


# # # # # # # # # # # # # # # # # # # # # # # # #
#   on trees with logbrownian rate variations   #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=c(paste0("LBtree50_",1:3),paste0("LBtree100_",1:3),paste0("LBtree200_",1:3))

data_LB=data.frame()


if(T){
  Sigmas=rep(0.03*(1:3),3)
  Ntips=rep(c(50,100,200),each=3)
  name_id=0
  for(name in Names){
    name_id=name_id+1
    ntip=Ntips[name_id]
    nedges=2*ntip-2
    for(id_tree in 1:10){
      do=F
      do=(sum(data_LB$name==name & data_LB$seed==id_tree)==0)
      if(! do) do = (data_LB$ngen[data_LB$name==name & data_LB$seed==id_tree]<801)
      if(do & file.exists(paste0("Cl0_",name,"_",id_tree,"_event_data.txt"))){
        load(paste0(name,"_",id_tree,".Rdata"))
        tree=read.tree(paste0("Cl0_",name,"_",id_tree,".tre"))
        
        edata=try(getEventData(tree, eventdata = paste0("Cl0_",name,"_",id_tree,"_event_data.txt"), burnin=0.2 ))
        if(! inherits(edata,"try-error")){
          
          load(paste0("Cl0_",name,"_",id_tree,".Rdata"))
          
          tips_edges=sapply(1:(tree$Nnode+1),function(i){which(tree$edge[,2]==i)})
          tot_time<-max(node.age(tree)$ages)
          f.lamb <-function(t,y){y[1]}
          f.mu<-function(t,y){y[1]}
          lamb_par<-c(1)
          mu_par<-c(0.1)
          result_cst <- fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par, f=1,cst.lamb=TRUE,cst.mu = T,fix.mu=F,dt=1e-3)
          
          BD_rates=c(abs(result_cst$lamb_par),abs(result_cst$mu_par))
          MSE_BD=mean((BD_rates[1]-speciation_rates)^2)
          MSE_log_BD=mean((log(BD_rates[1])-log(speciation_rates))^2)
          
          #### get ClaDS0 results
          rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),])}))
          rep2=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),-c((2*ntip+2):(2*ntip+4))])}))
          
          meanRate = rep(0,nedges+3)
          for(j in 1:3){
            meanRate=meanRate+colSums(cbind((rep[[j]][,1]),rep[[j]][,2],(rep[[j]][,3]),(rep[[j]][,2]+rep[[j]][,4:(nedges+3)])))
          }
          meanRate=exp(meanRate/(3*nrow(rep[[1]])))
          rate=exp(get_rates(tree,log(meanRate[-c(1:2)])))
          sigma_inf=log(meanRate[1])
          alpha_inf=meanRate[2]
          lambda_0_inf=rate[1]
          ClaDS0_rate=rate[-1]
          
          #### DR's
          DR_rates=DR(tree)
          tip_rates=speciation_rates[tips_edges]
          ClaDS_tip_rates=ClaDS0_rate[tips_edges]
          
          #### and BAMM's
          summary(edata)
          tree=read.tree(paste0("Cl0_",name,"_",id_tree,".tre"))
          
          mcmcout <- read.csv(paste0("Cl0_",name,"_",id_tree,"_mcmc_out.txt"), header=T)
          burnstart <- floor(0.1 * nrow(mcmcout))
          postburn <- mcmcout[burnstart:nrow(mcmcout), ]
          ES_Nshifts=effectiveSize(postburn$N_shifts)
          ES_LP=effectiveSize(postburn$logLik)
          
          BAMM_rates=getMeanBranchLengthTree(edata)$phy$edge.length
          BAMM_div_rates=getMeanBranchLengthTree(edata,rate = "ndr")$phy$edge.length
          BAMM_ext_rates=getMeanBranchLengthTree(edata,rate = "extinction")$phy$edge.length
          lm_BAMM=try(lm(BAMM_rates~speciation_rates))
          lm_BAMM_div=try(lm(BAMM_div_rates~speciation_rates))
          if(inherits(lm_BAMM_div, "try-error")){
            lm_BAMM_div=list(coefficients=rep(NA,2))
          }
          lm_log_BAMM_div=try(lm(log(BAMM_div_rates)~log(speciation_rates)))
          if(inherits(lm_log_BAMM_div, "try-error")){
            lm_log_BAMM_div=list(coefficients=rep(NA,2))
          }
          lm_log_BAMM=try(lm(log(BAMM_rates)~log(speciation_rates)))
          if(inherits(lm_log_BAMM, "try-error")){
            lm_log_BAMM=list(coefficients=rep(NA,2))
          }
          lm_ClaDS=try(lm(ClaDS0_rate~speciation_rates))
          if(inherits(lm_ClaDS, "try-error")){
            lm_ClaDS=list(coefficients=rep(NA,2))
          }
          lm_log_ClaDS=try(lm(log(ClaDS0_rate)~log(speciation_rates)))
          if(inherits(lm_log_ClaDS, "try-error")){
            lm_log_ClaDS=list(coefficients=rep(NA,2))
          }
          
          lm_ClaDS_tips=try(lm(ClaDS_tip_rates~tip_rates))
          if(inherits(lm_ClaDS_tips, "try-error")){
            lm_ClaDS_tips=list(coefficients=rep(NA,2))
          }
          lm_log_ClaDS_tips=try(lm(log(ClaDS_tip_rates)~log(tip_rates)))
          if(inherits(lm_log_ClaDS_tips, "try-error")){
            lm_log_ClaDS_tips=list(coefficients=rep(NA,2))
          }
          lm_DR=try(lm(DR_rates~tip_rates))
          if(inherits(lm_DR, "try-error")){
            lm_DR=list(coefficients=rep(NA,2))
          }
          lm_log_DR=try(lm(log(DR_rates)~log(tip_rates)))
          if(inherits(lm_log_DR, "try-error")){
            lm_log_DR=list(coefficients=rep(NA,2))
          }

          if(inherits(lm_BAMM, "try-error")){
            lm_BAMM=list(coefficients=rep(NA,2))
          }
          
          MSE_BAMM=mean((speciation_rates-BAMM_rates)^2)
          MSE_log_BAMM=mean((log(speciation_rates)-log(BAMM_rates))^2)
          MSE_BAMM_div=mean((speciation_rates-BAMM_div_rates)^2)
          MSE_log_BAMM_div=mean((log(speciation_rates)-log(BAMM_div_rates))^2)
          MSE_ClaDS=mean((speciation_rates-ClaDS0_rate)^2)
          MSE_log_ClaDS=mean((log(speciation_rates)-log(ClaDS0_rate))^2)
          MSE_ClaDS_tips=mean((tip_rates-ClaDS_tip_rates)^2)
          MSE_log_ClaDS_tips=mean((log(tip_rates)-log(ClaDS_tip_rates))^2)
          MSE_DR=mean((tip_rates-DR_rates)^2)
          MSE_log_DR=mean((log(tip_rates)-log(DR_rates))^2)
          
          var_true=var(speciation_rates)
          var_log_true=var(log(speciation_rates))
          var_true_tips=var(tip_rates)
          var_log_true_tips=var(log(tip_rates))
          
          var_BAMM=var(BAMM_rates)/var_true
          var_log_BAMM=var(log(BAMM_rates))/var_log_true
          var_BAMM_div=var(BAMM_div_rates)/var_true
          var_log_BAMM_div=var(log(BAMM_div_rates))/var_log_true
          var_ClaDS=var(ClaDS0_rate)/var_true
          var_log_ClaDS=var(log(ClaDS0_rate))/var_log_true
          var_ClaDS_tips=var(ClaDS_tip_rates)/var_true_tips
          var_log_ClaDS_tips=var(log(ClaDS_tip_rates))/var_log_true_tips
          var_DR=var(DR_rates)/var_true_tips
          var_log_DR=var(log(DR_rates))/var_log_true_tips
          
          do=(sum(data_LB$name==name & data_LB$seed==id_tree)==0)
          
          n_new_branches=sum(speciation_rates!=0.1)
          new_rate=speciation_rates[speciation_rates!=0.1][1]
          line_df=data.frame(name=name,seed=id_tree,ES_Nshifts=ES_Nshifts,ES_LP=ES_LP,ntips=ntip,ngen=length(edata$numberEvents),
                             sigma=Sigmas[name_id],sigma_inf=sigma_inf,alpha_inf=alpha_inf,shift=1,
                             n_new_branches=n_new_branches,new_rate=new_rate,
                             cor_BAMM=cor(BAMM_rates,speciation_rates),slope_BAMM=lm_BAMM$coefficients[2],
                             cor_log_BAMM=cor(log(BAMM_rates),log(speciation_rates)),slope_log_BAMM=lm_log_BAMM$coefficients[2],
                             MSE_BAMM=MSE_BAMM,MSE_log_BAMM=MSE_log_BAMM,
                             rel_error_BAMM=exp(mean(log(BAMM_rates)-log(speciation_rates))),
                             var_BAMM=var_BAMM,var_log_BAMM=var_log_BAMM,BAMM_ext=mean(BAMM_ext_rates),
                             cor_BAMM_div=cor(BAMM_div_rates,speciation_rates),slope_BAMM_div=lm_BAMM_div$coefficients[2],
                             cor_log_BAMM_div=cor(log(BAMM_div_rates),log(speciation_rates)),slope_log_BAMM_div=lm_log_BAMM_div$coefficients[2],
                             MSE_BAMM_div=MSE_BAMM_div,MSE_log_BAMM_div=MSE_log_BAMM_div,
                             rel_error_BAMM_div=exp(mean(log(BAMM_div_rates)-log(speciation_rates))),
                             var_BAMM_div=var_BAMM_div,var_log_BAMM_div=var_log_BAMM_div,
                             var_log_ClaDS=var_log_ClaDS,var_ClaDS=var_ClaDS,
                             MSE_log_BD=MSE_log_BD,MSE_BD=MSE_BD,BD_birth=BD_rates[1],BD_death=BD_rates[2],rel_error_BD=exp(sum(log(BD_rates[1]/speciation_rates))/nedges),
                             cor_ClaDS=cor(ClaDS0_rate,speciation_rates),slope_ClaDS=lm_ClaDS$coefficients[2],
                             cor_log_ClaDS=cor(log(ClaDS0_rate),log(speciation_rates)),slope_log_ClaDS=lm_log_ClaDS$coefficients[2],
                             MSE_ClaDS=MSE_ClaDS,MSE_log_ClaDS=MSE_log_ClaDS,
                             rel_error_ClaDS=exp(mean(log(ClaDS0_rate)-log(speciation_rates))),
                             cor_DR=cor(DR_rates,tip_rates),slope_DR=lm_DR$coefficients[2],
                             cor_log_DR=cor(log(DR_rates),log(tip_rates)),slope_log_DR=lm_log_DR$coefficients[2],
                             MSE_DR=MSE_DR,MSE_log_DR=MSE_log_DR,
                             rel_error_DR=exp(mean(log(DR_rates)-log(tip_rates))),
                             var_DR=var_DR,var_log_DR=var_log_DR,
                             cor_ClaDS_tips=cor(ClaDS_tip_rates,tip_rates),slope_ClaDS_tips=lm_ClaDS_tips$coefficients[2],
                             cor_log_ClaDS_tips=cor(log(ClaDS_tip_rates),log(tip_rates)),slope_log_ClaDS_tips=lm_log_ClaDS_tips$coefficients[2],
                             MSE_ClaDS_tips=MSE_ClaDS_tips,MSE_log_ClaDS_tips=MSE_log_ClaDS_tips,
                             rel_error_ClaDS_tips=exp(mean(log(ClaDS_tip_rates)-log(tip_rates))),
                             var_ClaDS_tips=var_ClaDS_tips,var_log_ClaDS_tips=var_log_ClaDS_tips)
          if(!do){
            data_LB[data_LB$name==name & data_LB$seed==id_tree,]=line_df
            
          }else{
            data_LB=rbind(data_LB,line_df)
          }
          
        }}
      save(data_LB,file="LBdata.Rdata")
      
    }
  }
  
  
}


# # # # # # # # # # # # # # # # # # # # # # # # #
#    on trees with evolving extinction rate     #
# # # # # # # # # # # # # # # # # # # # # # # # #


Names=c(paste0("CDtree50_",1:3),paste0("CDtree100_",1:3),paste0("CDtree200_",1:3))

data_CD=data.frame()


if(T){
  Sigmas=rep(0.03*(1:3),3)
  Ntips=rep(c(50,100,200),each=3)
  name_id=0
  for(name in Names){
    name_id=name_id+1
    ntip=Ntips[name_id]
    nedges=2*ntip-2
    for(id_tree in 1:10){
      do=F
      do=(sum(data_CD$name==name & data_CD$seed==id_tree)==0)
      if(! do) do = (data_CD$ngen[data_CD$name==name & data_CD$seed==id_tree]<801)
      if(do & file.exists(paste0("Cl0_",name,"_",id_tree,"_event_data.txt"))){
        load(paste0(name,"_",id_tree,".Rdata"))
        tree=read.tree(paste0("Cl0_",name,"_",id_tree,".tre"))
        
        edata=try(getEventData(tree, eventdata = paste0("Cl0_",name,"_",id_tree,"_event_data.txt"), burnin=0.2 ))
        if(! inherits(edata,"try-error")){
          
          load(paste0("Cl0_",name,"_",id_tree,".Rdata"))
          
          tips_edges=sapply(1:(tree$Nnode+1),function(i){which(tree$edge[,2]==i)})
          tot_time<-max(node.age(tree)$ages)
          f.lamb <-function(t,y){y[1]}
          f.mu<-function(t,y){y[1]}
          lamb_par<-c(1)
          mu_par<-c(0.1)
          result_cst <- fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par, f=1,cst.lamb=TRUE,cst.mu = T,fix.mu=F,dt=1e-3)
          
          BD_rates=c(abs(result_cst$lamb_par),abs(result_cst$mu_par))
          MSE_BD=mean((BD_rates[1]-speciation_rates)^2)
          MSE_log_BD=mean((log(BD_rates[1])-log(speciation_rates))^2)
          
          #### get ClaDS0 results
          rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),])}))
          rep2=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),-c((2*ntip+2):(2*ntip+4))])}))
          
          meanRate = rep(0,nedges+3)
          for(j in 1:3){
            meanRate=meanRate+colSums(cbind((rep[[j]][,1]),rep[[j]][,2],(rep[[j]][,3]),(rep[[j]][,2]+rep[[j]][,4:(nedges+3)])))
          }
          meanRate=exp(meanRate/(3*nrow(rep[[1]])))
          rate=exp(get_rates(tree,log(meanRate[-c(1:2)])))
          sigma_inf=log(meanRate[1])
          alpha_inf=meanRate[2]
          lambda_0_inf=rate[1]
          ClaDS0_rate=rate[-1]
          
          #### DR's
          DR_rates=DR(tree)
          tip_rates=speciation_rates[tips_edges]
          ClaDS_tip_rates=ClaDS0_rate[tips_edges]
          
          #### and BAMM's
          summary(edata)
          tree=read.tree(paste0("Cl0_",name,"_",id_tree,".tre"))
          
          mcmcout <- read.csv(paste0("Cl0_",name,"_",id_tree,"_mcmc_out.txt"), header=T)
          burnstart <- floor(0.1 * nrow(mcmcout))
          postburn <- mcmcout[burnstart:nrow(mcmcout), ]
          ES_Nshifts=effectiveSize(postburn$N_shifts)
          ES_LP=effectiveSize(postburn$logLik)
          
          BAMM_rates=getMeanBranchLengthTree(edata)$phy$edge.length
          BAMM_div_rates=getMeanBranchLengthTree(edata,rate = "ndr")$phy$edge.length
          BAMM_ext_rates=getMeanBranchLengthTree(edata,rate = "extinction")$phy$edge.length
          lm_BAMM=try(lm(BAMM_rates~speciation_rates))
          lm_BAMM_div=try(lm(BAMM_div_rates~speciation_rates))
          if(inherits(lm_BAMM_div, "try-error")){
            lm_BAMM_div=list(coefficients=rep(NA,2))
          }
          lm_log_BAMM_div=try(lm(log(BAMM_div_rates)~log(speciation_rates)))
          if(inherits(lm_log_BAMM_div, "try-error")){
            lm_log_BAMM_div=list(coefficients=rep(NA,2))
          }
          lm_log_BAMM=try(lm(log(BAMM_rates)~log(speciation_rates)))
          if(inherits(lm_log_BAMM, "try-error")){
            lm_log_BAMM=list(coefficients=rep(NA,2))
          }
          lm_ClaDS=try(lm(ClaDS0_rate~speciation_rates))
          if(inherits(lm_ClaDS, "try-error")){
            lm_ClaDS=list(coefficients=rep(NA,2))
          }
          lm_log_ClaDS=try(lm(log(ClaDS0_rate)~log(speciation_rates)))
          if(inherits(lm_log_ClaDS, "try-error")){
            lm_log_ClaDS=list(coefficients=rep(NA,2))
          }
          
          lm_ClaDS_tips=try(lm(ClaDS_tip_rates~tip_rates))
          if(inherits(lm_ClaDS_tips, "try-error")){
            lm_ClaDS_tips=list(coefficients=rep(NA,2))
          }
          lm_log_ClaDS_tips=try(lm(log(ClaDS_tip_rates)~log(tip_rates)))
          if(inherits(lm_log_ClaDS_tips, "try-error")){
            lm_log_ClaDS_tips=list(coefficients=rep(NA,2))
          }
          lm_DR=try(lm(DR_rates~tip_rates))
          if(inherits(lm_DR, "try-error")){
            lm_DR=list(coefficients=rep(NA,2))
          }
          lm_log_DR=try(lm(log(DR_rates)~log(tip_rates)))
          if(inherits(lm_log_DR, "try-error")){
            lm_log_DR=list(coefficients=rep(NA,2))
          }
          
          if(inherits(lm_BAMM, "try-error")){
            lm_BAMM=list(coefficients=rep(NA,2))
          }
          
          MSE_BAMM=mean((speciation_rates-BAMM_rates)^2)
          MSE_log_BAMM=mean((log(speciation_rates)-log(BAMM_rates))^2)
          MSE_BAMM_div=mean((speciation_rates-BAMM_div_rates)^2)
          MSE_log_BAMM_div=mean((log(speciation_rates)-log(BAMM_div_rates))^2)
          MSE_ClaDS=mean((speciation_rates-ClaDS0_rate)^2)
          MSE_log_ClaDS=mean((log(speciation_rates)-log(ClaDS0_rate))^2)
          MSE_ClaDS_tips=mean((tip_rates-ClaDS_tip_rates)^2)
          MSE_log_ClaDS_tips=mean((log(tip_rates)-log(ClaDS_tip_rates))^2)
          MSE_DR=mean((tip_rates-DR_rates)^2)
          MSE_log_DR=mean((log(tip_rates)-log(DR_rates))^2)
          
          var_true=var(speciation_rates)
          var_log_true=var(log(speciation_rates))
          var_true_tips=var(tip_rates)
          var_log_true_tips=var(log(tip_rates))
          
          var_BAMM=var(BAMM_rates)/var_true
          var_log_BAMM=var(log(BAMM_rates))/var_log_true
          var_BAMM_div=var(BAMM_div_rates)/var_true
          var_log_BAMM_div=var(log(BAMM_div_rates))/var_log_true
          var_ClaDS=var(ClaDS0_rate)/var_true
          var_log_ClaDS=var(log(ClaDS0_rate))/var_log_true
          var_ClaDS_tips=var(ClaDS_tip_rates)/var_true_tips
          var_log_ClaDS_tips=var(log(ClaDS_tip_rates))/var_log_true_tips
          var_DR=var(DR_rates)/var_true_tips
          var_log_DR=var(log(DR_rates))/var_log_true_tips
          
          do=(sum(data_CD$name==name & data_CD$seed==id_tree)==0)
          
          n_new_branches=sum(speciation_rates!=0.1)
          new_rate=speciation_rates[speciation_rates!=0.1][1]
          line_df=data.frame(name=name,seed=id_tree,ES_Nshifts=ES_Nshifts,ES_LP=ES_LP,ntips=ntip,ngen=length(edata$numberEvents),
                             sigma=Sigmas[name_id],sigma_inf=sigma_inf,alpha_inf=alpha_inf,shift=1,
                             n_new_branches=n_new_branches,new_rate=new_rate,
                             cor_BAMM=cor(BAMM_rates,speciation_rates),slope_BAMM=lm_BAMM$coefficients[2],
                             cor_log_BAMM=cor(log(BAMM_rates),log(speciation_rates)),slope_log_BAMM=lm_log_BAMM$coefficients[2],
                             MSE_BAMM=MSE_BAMM,MSE_log_BAMM=MSE_log_BAMM,
                             rel_error_BAMM=exp(mean(log(BAMM_rates)-log(speciation_rates))),
                             var_BAMM=var_BAMM,var_log_BAMM=var_log_BAMM,BAMM_ext=mean(BAMM_ext_rates),
                             cor_BAMM_div=cor(BAMM_div_rates,speciation_rates),slope_BAMM_div=lm_BAMM_div$coefficients[2],
                             cor_log_BAMM_div=cor(log(BAMM_div_rates),log(speciation_rates)),slope_log_BAMM_div=lm_log_BAMM_div$coefficients[2],
                             MSE_BAMM_div=MSE_BAMM_div,MSE_log_BAMM_div=MSE_log_BAMM_div,
                             rel_error_BAMM_div=exp(mean(log(BAMM_div_rates)-log(speciation_rates))),
                             var_BAMM_div=var_BAMM_div,var_log_BAMM_div=var_log_BAMM_div,
                             var_log_ClaDS=var_log_ClaDS,var_ClaDS=var_ClaDS,
                             MSE_log_BD=MSE_log_BD,MSE_BD=MSE_BD,BD_birth=BD_rates[1],BD_death=BD_rates[2],rel_error_BD=exp(sum(log(BD_rates[1]/speciation_rates))/nedges),
                             cor_ClaDS=cor(ClaDS0_rate,speciation_rates),slope_ClaDS=lm_ClaDS$coefficients[2],
                             cor_log_ClaDS=cor(log(ClaDS0_rate),log(speciation_rates)),slope_log_ClaDS=lm_log_ClaDS$coefficients[2],
                             MSE_ClaDS=MSE_ClaDS,MSE_log_ClaDS=MSE_log_ClaDS,
                             rel_error_ClaDS=exp(mean(log(ClaDS0_rate)-log(speciation_rates))),
                             cor_DR=cor(DR_rates,tip_rates),slope_DR=lm_DR$coefficients[2],
                             cor_log_DR=cor(log(DR_rates),log(tip_rates)),slope_log_DR=lm_log_DR$coefficients[2],
                             MSE_DR=MSE_DR,MSE_log_DR=MSE_log_DR,
                             rel_error_DR=exp(mean(log(DR_rates)-log(tip_rates))),
                             var_DR=var_DR,var_log_DR=var_log_DR,
                             cor_ClaDS_tips=cor(ClaDS_tip_rates,tip_rates),slope_ClaDS_tips=lm_ClaDS_tips$coefficients[2],
                             cor_log_ClaDS_tips=cor(log(ClaDS_tip_rates),log(tip_rates)),slope_log_ClaDS_tips=lm_log_ClaDS_tips$coefficients[2],
                             MSE_ClaDS_tips=MSE_ClaDS_tips,MSE_log_ClaDS_tips=MSE_log_ClaDS_tips,
                             rel_error_ClaDS_tips=exp(mean(log(ClaDS_tip_rates)-log(tip_rates))),
                             var_ClaDS_tips=var_ClaDS_tips,var_log_ClaDS_tips=var_log_ClaDS_tips)
          if(!do){
            data_CD[data_CD$name==name & data_CD$seed==id_tree,]=line_df
            
          }else{
            data_CD=rbind(data_CD,line_df)
          }
          
        }}
      save(data_CD,file="CDdata.Rdata")
      
    }
  }
  
  
}

#### type I and type II error for the oneShift trees ####
library(fields)
library(plotrix)
library(coda)
library(vioplot)
library(apTreeshape)
library(caroline)

data_shift=data.frame()

for(Name in c("LNormMHShift_d3","LNormMHShift_d2",
              "LNormMHShift_d1","LNormMHShift_x1",
              "LNormMHShift_x1-5","LNormMHShift_x2",
              "LNormMHShift_x3","LNormMHShift_x4")){
  continue=T
  k=1
  while(k<501 & continue){
    name=paste(Name,k)
    k=k+1
      TRY=try(load(file=name))
      if(inherits(TRY,"try-error")){
        continue=F
      }else{
        print(name)
        phylo=tree
        ntip = tree$Nnode+1
        nedges=ntip*2-2
        edge = tree$edge
        edge.length = phylo$edge.length
        
        rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),])}))
        rep2=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),-c((2*ntip+2):(2*ntip+4))])}))
        
        n=nrow(rep2[[1]])
        n_sample=1000
        N=3
        n_sample=min(n_sample,n*N)
        ind=sample(1:(N*n),n_sample,replace = F)
        ind_chain=1+floor((ind-1)/n)
        ind=ind-(ind_chain-1)*n
        sampled_chain=mcmc(t(sapply(1:n_sample,function(j){
          vect=rep2[[ind_chain[j]]][ind[j],1:(nedges+3)]
          vect[-c(1:3)]=get_rates(tree,c(vect[3],(vect[-c(1:3)]+vect[2])))[-1]
          return(vect)})))

        rep=matrix(0,nedges,nedges)
        for(i in 1:n_sample){
          cat("\r",i,"\r")
          M=sapply(1:nedges,function(j){sampled_chain[i,4:(nedges+3)]-sampled_chain[i,j+3]})
          rep=rep+(M>0)
        }
        cat("\r","","\r")
        rep=rep/n_sample
        M=sapply(1:nedges,function(j){(true.rate-true.rate[j])/true.rate})
        seuilM=0
        M=matrix(0,nedges,nedges)+1*(M>seuilM)-1*(M<(-seuilM))
        Colors= colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 )
        seuil=0.05
        rep_seuil=matrix(0,nedges,nedges)+(rep<seuil)*(1-rep)-((1-rep)<seuil)*(rep)
        composite=rep_seuil
        for(j in 1:(nedges)){
          composite[j,1:j]=M[j,1:j]
        }
        image(rep,col=Colors)
        Colors= colorRampPalette(c("steelblue2","gray90","darkorange"))( 3 )
        par(mfrow=c(1,2))
        plot.with.rate(tree,true.rate,lwd=2)
        image(composite,col=Colors)
        lines(c(0,1),c(0,1))
        
        shiftM=0
        shiftr=0
        for(j in 1:(nedges-1)){
          shiftM=shiftM+sum(M[(j+1):nedges,j]!=0)
          shiftr=shiftr+sum(rep_seuil[j,(j+1):nedges]!=0)
        }
        Error1=0
        for(j in 1:(nedges-1)){
          Error1=Error1+sum(rep_seuil[j,(j+1):nedges]!=0 & M[(j+1):nedges,j]==0)
        }
        Error1=Error1*(1/(((nedges-1)*nedges/2)-shiftM))
        Error2=0
        for(j in 1:(nedges-1)){
          Error2=Error2+sum(rep_seuil[j,(j+1):nedges]==0 & M[(j+1):nedges,j]!=0)
        }
        Error2=Error2/shiftM
        ErrorS=0
        for(j in 1:(nedges-1)){
          ErrorS=ErrorS+sum(rep_seuil[j,(j+1):nedges] * M[(j+1):nedges,j] <0)
        }
        ErrorS=ErrorS/shiftr
        if(is.nan(ErrorS)) ErrorS=0
        if(is.nan(Error2)) Error2=0
        
        print(ErrorS)
        print(Error1)
        print(Error2)

        data_shift=rbind(data_shift,data.frame(name=name,amplitude=true[true!=0.1][1],
                                               size=table(as.factor(true))[levels(as.factor(true))!="0.1"],
                                               Error1=Error1,
                                               Error2=Error2,ErrorS=ErrorS))

        data_shift[nrow(data_shift),]
        row.names(data_shift)=NULL}
  }
  save(data_shift,file="Shift_Analysis.Rdata")
}


#### test of ClaDS2 ####

# # # # # # # # # # # # # # # # # # # # # # # # #
#                on ClaDS2 trees                #
# # # # # # # # # # # # # # # # # # # # # # # # #

id=1 # in 1:8
replicate=1 # in 1:5
l0=0.1
epsilon=(rep(c(0.9,0.1),4))[id]
sig=rep(c(0.7,0,0.2,0.2),each=2)[id]
al=rep(c(0.9,1,1,0.9),each=2)[id]*exp(-(sig^2)/2)

obj= sim_ClaDS( lamb_par=l0, 
                mu_par=epsilon,
                sigma=sig,
                lamb_shift=al,
                condition="taxa", 
                taxa.stop = 100)   
tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

sampler_ClaDS0 = run_ClaDS0(tree=tree, 
                     name=paste0("noDeath_",id,"_",replicate,".Rdata"),   
                     nCPU=1,          
                     pamhLocalName = paste0("local",id,"_",replicate,"_"), 
                     iteration=10000000, thin=20000,update=1000, adaptation=10000)  

# prepare the sampler
sampler = prepare_ClaDS(tree=tree, 
                        sample_fraction=1, 
                        Nchain = 3,      
                        nlambda = 1000,   
                        nt = 30,           
                        model_id="ClaDS2",  
                        res_ClaDS0 = sampler_ClaDS0) 

for(i in 1:5000){
  sampler = fit_ClaDS(           
    mcmcSampler = sampler,      
    iterations = 1000,           
    thin = 50,  
    nCPU = 1) 
  save.image(paste0("Cl2_",id,"_",replicate,".Rdata"))
}

#### comparison with other methods ####

#### extract the results from the MCMC chains ####

# # # # # # # # # # # # # # # # # # # # # # # # #
#                ClaDS2 results                 #
# # # # # # # # # # # # # # # # # # # # # # # # #

data_Cl2=data.frame()
arg=0
for(id in c(1:8)){
  print(id)
  chainCor=list()
  for(replicate in c(1:5)){
    t=try(load(paste("Cl2_",id,"_",replicate,".Rdata",sep="")))
    do=!(inherits(t,"try-error"))
    if (do){ 
      nR=nrow(sampler$codaChain[[1]])
      chains=mcmc.list(lapply(sampler$codaChain,function(x){mcmc(x[-(1:floor(nR/2)),1:npar])}))
      gelm=try(gelman.diag(chains)$psrf[,1])
      if(inherits(gelm,"try-error") | is.nan(gelm[1])) {
        gelm=10
        ind.gelm=npar+1
      }else{
        ind.gelm=which.max(gelm)
      }
      true_div=true.rate*(1-muT)
      true_ext=true.rate*(muT)
      
      plot.chains=mcmc.list(lapply(sampler$codaChain,function(x){mcmc(x[-(1:(floor(1*nR/2))),])}))
      
      
      MAPS=sapply(1:npar, function(i){D=density(c(plot.chains[[1]][seq(1,nrow(plot.chains[[1]]),10),i],                                       plot.chains[[2]][seq(1,nrow(plot.chains[[1]]),10),i],
                                                  plot.chains[[3]][seq(1,nrow(plot.chains[[1]]),10),i]))
      return(D$x[which.max(D$y)])})
      
      MAPS[c(1,2,4)]=exp(MAPS[c(1,2,4)])
      rateMAPS=exp(MAPS[-c(1:4,(npar+1):(npar+3))])
      lm=try(lm(rateMAPS~true.rate))
      if(inherits(lm,"try-error")){
        lm=list(coefficients=rep(NA,2))
      }
      lm_log=try(lm(log(rateMAPS)~log(true.rate)))
      if(inherits(lm_log,"try-error")){
        lm_log=list(coefficients=rep(NA,2))
      }
      
      div_MAPS=rateMAPS*(1-MAPS[3])
      lmMAPS_div=try(lm(div_MAPS~true_div))
      if(inherits(lmMAPS_div,"try-error")){
        lmMAPS_div=list(coefficients=rep(NA,2))
      }
      
      ext_MAPS=rateMAPS*(MAPS[3])
      lmMAPS_ext=try(lm(ext_MAPS~true_ext))
      if(inherits(lmMAPS_ext,"try-error")){
        lmMAPS_ext=list(coefficients=rep(NA,2))
      }
      
      MSE=mean((rateMAPS-true.rate)^2)
      MSE_div=mean((div_MAPS-true.rate*(1-muT))^2)
      MSE_ext=mean((ext_MAPS-true.rate*(muT))^2)
      data_Cl2=rbind(data_Cl2,data.frame(nrow=nR,
                                         minBL=max(node.depth.edgelength(tree))/min(tree$edge.length),
                                         sigma=sig,sigma_inf=MAPS[1],
                                         alpha=al,
                                         alpha_inf=MAPS[2],
                                         ind.gelm=ind.gelm,mu=muT,nt=size_nt,
                                         mu_inf=MAPS[3],
                                         id=id,replicate=replicate,
                                         lambda_0=l0,lambda_0_inf=MAPS[4],
                                         rel_err=exp(mean(log(rateMAPS)-log(true.rate))),
                                         rel_err_diff=exp(mean(log(rateMAPS*(1-MAPS[3]))
                                                               -log(true.rate*(1-muT)))),
                                         cor=cor(true.rate,rateMAPS), MSE=MSE,MSE_div=MSE_div,MSE_ext=MSE_ext,
                                         slope=lm$coefficients[2],
                                         slope_div=lmMAPS_div$coefficients[2],
                                         inter_div=lmMAPS_div$coefficients[1],
                                         slope_ext=lmMAPS_ext$coefficients[2],
                                         inter_ext=lmMAPS_ext$coefficients[1],
                                         rel_err_ext=exp(mean(log(ext_MAPS)-log(true_ext))),
                                         cor_div=cor(true_div,div_MAPS),
                                         cor_ext=cor(true_ext,ext_MAPS),cor_ext_log=cor(log(true_ext),log(ext_MAPS)),
                                         slope_log=lm_log$coefficients[2],
                                         cor_log=cor(log(true.rate),log(rateMAPS))))
      
    }    }
  
}
save(data_Cl2,file="data_Cl2_Cl2.Rdata")


# # # # # # # # # # # # # # # # # # # # # # # # #
#              BAMM and DR results              #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=1:8
data_BAMM=data.frame()

plot_tree=T
plot_cor=T

if(T){
  name_id=0
  Mu=rep(c(0.9,0.1),4)
  Sigmas=rep(c(0.7,0,0.2,0.2),each=2)
  Alphas=c(0.9,0.9,1,1,1,1,0.8,0.8)*exp(-Sigmas^2/2)
  for(nameSc in Names){
    name_id=name_id+1
    for(id_tree in 1:5){
      do=T
      do=(sum(data_BAMM$name==nameSc & data_BAMM$seed==id_tree)==0)
      if(! do) do = (data_BAMM$ngen[data_BAMM$name==nameSc & data_BAMM$seed==id_tree]<801)
      if(do & file.exists(paste0(nameSc,"_",id_tree,"_event_data.txt"))){
        tree=read.tree(paste0(nameSc,"_",id_tree,".tre"))
        
        edata=try(getEventData(tree, eventdata = paste0(nameSc,"_",id_tree,"_event_data.txt"), burnin=0.2 ))
        if(! inherits(edata,"try-error")){
          t=try(load(paste0("~/Documents/ClaDS/Revisions/cleaned_data/test_Cl2/ClaDS2/Cl2_",nameSc,"_",id_tree,".Rdata")))
          
          summary(edata)
          if(plot_tree){
            par(mfrow=c(1,2))
            plot.with.rate(tree,true.rate,
                           log=T,        # should the rates be plotted on a log scale?
                           lwd=3)        # width of the branches
            
            plot.bammdata(edata, lwd=2,  legend=T, pal = colorPal )
            title(main=paste0(nameSc,"_",id_tree))
          }
          mcmcout <- read.csv(paste0(nameSc,"_",id_tree,"_mcmc_out.txt"), header=T)
          burnstart <- floor(0.1 * nrow(mcmcout))
          postburn <- mcmcout[burnstart:nrow(mcmcout), ]
          ES_Nshifts=effectiveSize(postburn$N_shifts)
          ES_LP=effectiveSize(postburn$logLik)
          
          BAMM_rates=getMeanBranchLengthTree(edata)$phy$edge.length/1000
          BAMM_div_rates=getMeanBranchLengthTree(edata,rate = "ndr")$phy$edge.length/1000
          BAMM_extinction=getMeanBranchLengthTree(edata,rate = "extinction")$phy$edge.length/1000
          lm=try(lm(BAMM_rates~true.rate))
          
          lm_log=try(lm(log(BAMM_rates)~log(true.rate)))
          if(inherits(lm_log, "try-error")){
            lm_log=list(coefficients=rep(NA,2))
          }
          
          true_div=(true.rate*(1-muT))
          lm_div=try(lm(BAMM_div_rates~true_div))
          if(inherits(lm_div, "try-error")){
            lm_div=list(coefficients=rep(NA,2))
          }
          
          true_ext=(true.rate*(muT))
          lm_ext=try(lm(BAMM_extinction~true_ext))
          if(inherits(lm_ext, "try-error")){
            lm_ext=list(coefficients=rep(NA,2))
          }
          
          if(plot_cor){
            par(mfrow=c(1,2))
            plot.with.rate(tree,BAMM_rates,log=T,lwd=2)
            plot(true.rate,BAMM_rates)
            title(main=cor(BAMM_rates,true.rate))
          }
          if(inherits(lm, "try-error")){
            lm=list(coefficients=rep(NA,2))
          }else{
            if(plot_cor & !is.na(lm$coefficients[2])) abline(a=lm$coefficients[1],b=lm$coefficients[2],lwd=2,col="red")
          }
          
          MSE=mean((BAMM_rates-true.rate)^2)
          MSE_div=mean((BAMM_div_rates-true.rate*(1-muT))^2)
          
          data_line=data.frame(name=nameSc,seed=id_tree,ES_Nshifts=ES_Nshifts,ES_LP=ES_LP,ntip=200,ngen=length(edata$numberEvents),
                               sigma=Sigmas[name_id],shift=Alphas[name_id],
                               MSE=MSE,MSE_div=MSE_div,mean_extinction=mean(BAMM_extinction),
                               mean_turn=mean(BAMM_extinction/BAMM_rates),geom_mean_turn=exp(mean(log(BAMM_extinction)-log(BAMM_rates))),
                               cor=cor(BAMM_rates,true.rate),slope=lm$coefficients[2],slope_div=lm_div$coefficients[2],
                               slope_ext=lm_ext$coefficients[2],
                               cor_log=cor(log(BAMM_rates),log(true.rate)),slope_log=lm_log$coefficients[2],
                               rel_error=exp(mean(log(BAMM_rates)-log(true.rate))),
                               cor_div=cor(BAMM_div_rates,true_div),cor_ext=cor(BAMM_extinction,true_ext),
                               rel_error_div=exp(mean(log(BAMM_div_rates)-log(true.rate*(1-Mu[name_id])))))
          do=(sum(data_BAMM$name==nameSc & data_BAMM$seed==id_tree)==0)
          if(!do){
            data_BAMM[data_BAMM$name==nameSc & data_BAMM$seed==id_tree,]=data_line
            
          }else{
            data_BAMM=rbind(data_BAMM,data_line)
          }
          
        }}
    }
  }
  
  save(data_BAMM,file="BAMM_Cl2.Rdata")
  
}


#### application to the bird data ####

# # # # # # # # # # # # # # # # # # # # # # # # #
#                   run ClaDS                   #
# # # # # # # # # # # # # # # # # # # # # # # # #

clade_names_concensus=dir("/BirdFamily/MCC_Clades/constraint_trees")
clade_names_MCC=dir("/BirdFamily/MCC_Clades/trees")
data.names=read.csv("/BirdFamily/MCC_Clades/clade_summary.csv")

args=1 # in 1:126
name=paste0("MCCclade_",args)
tree=read.nexus(paste0("/BirdFamily/MCC_Clades/trees/",clade_names_MCC[args]))

data.id=which(data.names$Clade==sub(".MCC.*", "", clade_names_MCC[args]))
tree <- drop.tip(tree,c(levels(data.names$X1st.Outgroup)[data.names$X1st.Outgroup[data.id]],
                        levels(data.names$X2nd.Outgroup)[data.names$X2nd.Outgroup[data.id]]))
sample_fraction = (tree$Nnode+1)/data.names[data.id,"Clade.Size"]
  
if(tree$Nnode>49){
    
  sampler = run_ClaDS0(tree=tree,       
                       name=paste0("Cl0_",name,".Rdata"),
                       nCPU=1,             
                       pamhLocalName = paste0("local",name),   
                       iteration=500000, 
                       thin=2000,    
                       update=1000, adaptation=5000) 
  
  sampler = prepare_ClaDS(tree=tree, 
                          sample_fraction=1,  
                          Nchain = 3,    
                          nlambda = 1000,   
                          nt = 30,           
                          model_id="ClaDS2",  
                          res_ClaDS0 = sampler)  
  
  for(i in 1:5000){
    sampler = fit_ClaDS(           
      mcmcSampler = sampler,      
      iterations = 1000,           
      thin = 50,  
      nCPU = 1) 
    save.image(paste0("Cl2_",name,".Rdata"))
  }
}

# # # # # # # # # # # # # # # # # # # # # # # # #
#        extract the results  for ClaDS         #
# # # # # # # # # # # # # # # # # # # # # # # # #

is.tip=function(tree){
  ntip=tree$Nnode+1
  rep=rep(F,2*(ntip-1))
  for(i in 1:(2*(ntip-1))){
    if(tree$edge[i,2]<=ntip) rep[i]=T
  }
  return(rep)
}

densities=list()
all.rates=c()
clade.rates.Full=list()
all.ratesFull=c()
dataBird=data.frame()
dataTip=data.frame()
n=0
rangeX=c()
rangeY=c()
rangeXFull=c()
rangeYFull=c()
rangeXe=c()
rangeYe=c()
for( id in (1:130)[]){
  runing=NA
  sampler=NULL
  j=0
  t=try(load(paste0("Cl2_MCCclade_",id,".Rdata")))
  do=!(is.null(sampler$chains))
  if (do){
    nR=nrow(sampler$codaChain[[1]])
    n=n+1
    plot.chains=mcmc.list(lapply(sampler$codaChain,function(x){mcmc(x[-(1:(1*floor(nR/2))),])}))
    
    for(k in 1:length(plot.chains)){
      plot.chains[[k]][,5:npar]=plot.chains[[k]][,5:npar]+plot.chains[[k]][,4]+
        sapply(alpha_effect,function(x){x*plot.chains[[k]][,2]})
      plot.chains[[k]][,2]=plot.chains[[k]][,2]/(1+plot.chains[[k]][,3])
    }
    
    MAPS=sapply(1:npar, function(i){D=density(c(plot.chains[[1]][seq(1,nrow(plot.chains[[1]]),10),i],                                       plot.chains[[2]][seq(1,nrow(plot.chains[[1]]),10),i],
                                                plot.chains[[3]][seq(1,nrow(plot.chains[[1]]),10),i]))
    return(D$x[which.max(D$y)])})
    
    MAPS[c(1,2,4)]=exp(MAPS[c(1,2,4)])
    rateMAPS=(MAPS[-c(1:4,(npar+1):(npar+3))])
    tip=is.tip(tree)
    dataBird=rbind(dataBird,data.frame(id=id,
                                       sample_fraction=sample_fraction,
                                       ntip=ntip,nrow=nR,
                                       gamma=gammaStat(tree),beta=maxlik.betasplit(tree)$max_lik,
                                       alpha=MAPS[2],
                                       mnR=MAPS[2]*exp((MAPS[1]^2)/2),
                                       mu=MAPS[3],
                                       sigma=MAPS[1],
                                       l0=MAPS[4],
                                       minRate=min(exp(rateMAPS)),maxRate=max(exp(rateMAPS)),
                                       diffRate=exp(diff(range(rateMAPS))),
                                       root_age=max(node.depth.edgelength(tree)),
                                       med_tip=median(rateMAPS[tip]),var_tip=var(rateMAPS[tip]),
                                       med_full=median(rateMAPS),var_full=var(rateMAPS)))
    
    d=density(rateMAPS[tip])
    dFull=density(rateMAPS)
    de=density(exp(rateMAPS[tip]))
    all.rates=c(all.rates,rateMAPS[tip])
    all.ratesFull=c(all.ratesFull,rateMAPS)
    rangeX=range(c(rangeX,d$x))
    rangeY=range(c(rangeY,d$y))
    rangeXFull=range(c(rangeXFull,dFull$x))
    rangeYFull=range(c(rangeYFull,dFull$y))
    rangeYe=range(c(rangeYe,de$y))
    rangeXe=range(c(rangeXe,de$x))
    densities[[n]]=list(id=id,x=d$x,y=d$y,ex=de$x,ey=de$y,xFull=dFull$x,yFull=dFull$y, ntip=ntip, all.rates=rateMAPS,rates.tips=rateMAPS[tip])
  }
}

save(all.ratesFull, all.rates,dataBird,rangeXFull,rangeYFull,densities,
     file="data_birds.Rdata")



# # # # # # # # # # # # # # # # # # # # # # # # #
#                    run BAMM                   #
# # # # # # # # # # # # # # # # # # # # # # # # #

for(name_tree in Names){
  for(id in 1:126){
    t=try(load(paste0(path,"BirdFamily/Cl2_MCCclade_",id,".Rdata")))
    if(! inherits(t, "try-error")){
      tips_edges=sapply(1:(tree$Nnode+1),function(i){which(tree$edge[,2]==i)})
      par(mfrow=c(1,2))
      plot(tree, show.tip.label = F)
      tree$edge.length[tree$edge.length==0]=max(tree$edge.length)/10000
      tree=force.ultrametric(tree, method = "extend")
      plot(tree, show.tip.label = F)
      
      print(min(tree$edge.length))
      write.tree(tree,file=paste0(path,"BirdFamily/BAMM/MCCclade_",id,".tre"))
      writeTemplate_withBAMMprior(tree, treeName = paste0("MCCclade_",id), 
                                  outfile = paste0(path,"BirdFamily/BAMM/template_MCCclade_",id,".txt"))
    }
  }
}


# # # # # # # # # # # # # # # # # # # # # # # # #
#              extract BAMM results             #
# # # # # # # # # # # # # # # # # # # # # # # # #

Names=paste0("MCCclade","_",1:126)
data_BB=data.frame()
densities=list()
all.rates=c()
all.ratesFull=c()
rangeX=c()
rangeY=c()
rangeXFull=c()
rangeYFull=c()
rangeYe=c()
rangeXe=c()

plot_tree=T
plot_cor=T

if(T){
  name_id=0
  n=0
  for(name_tree in Names){
    name_id=name_id+1
    
    if(file.exists(paste0(name_tree,"_event_data.txt"))){
      n=n+1
      tree=read.tree(paste0(name_tree,".tre"))
      
      ntip=tree$Nnode+1
      nedges=nrow(tree$edge)
      
      edata=try(getEventData(tree, eventdata = paste0(name_tree,"_event_data.txt"), burnin=0.2 ))
      if(! inherits(edata,"try-error")){
        
        tips_edges=sapply(1:(tree$Nnode+1),function(i){which(tree$edge[,2]==i)})
        tot_time<-max(node.age(tree)$ages)
        f.lamb <-function(t,y){y[1]}
        f.mu<-function(t,y){y[1]}
        lamb_par<-c(1)
        mu_par<-c(0.1)
        result_cst <- fit_bd(tree,tot_time,f.lamb,f.mu,lamb_par,mu_par, f=1,cst.lamb=TRUE,cst.mu = T,fix.mu=F,dt=1e-3)
        
        BD_rates=c(abs(result_cst$lamb_par),abs(result_cst$mu_par))
        
        #### DR's results
        DR_rates=DR(tree)
        
        #### and BAMM's
        summary(edata)

        tree=read.tree(paste0(name_tree,".tre"))
        
        mcmcout <- read.csv(paste0(name_tree,"_mcmc_out.txt"), header=T)
        burnstart <- floor(0.1 * nrow(mcmcout))
        postburn <- mcmcout[burnstart:nrow(mcmcout), ]
        ES_Nshifts=effectiveSize(postburn$N_shifts)
        ES_LP=effectiveSize(postburn$logLik)
        
        BAMM_rates=getMeanBranchLengthTree(edata)$phy$edge.length
        BAMM_div_rates=getMeanBranchLengthTree(edata,rate = "ndr")$phy$edge.length
        BAMM_ext_rates=getMeanBranchLengthTree(edata,rate = "extinction")$phy$edge.length
        
        if(plot_cor){
          par(mfrow=c(1,1))
          plot.with.rate(tree,BAMM_rates,log=T,lwd=2)
        }
        
        var_BAMM=var(BAMM_rates)
        var_log_BAMM=var(log(BAMM_rates))
        var_BAMM_div=var(BAMM_div_rates)
        var_log_BAMM_div=var(log(BAMM_div_rates))
        var_DR=var(DR_rates)
        var_log_DR=var(log(DR_rates))
        
        ### compute the rates densities
        
        d=density(log(BAMM_rates[tips_edges]))
        dFull=density(log(BAMM_rates))
        de=density((BAMM_rates[tips_edges]))
        all.rates=c(all.rates,log(BAMM_rates[tips_edges]))
        all.ratesFull=c(all.ratesFull,log(BAMM_rates))
        rangeX=range(c(rangeX,d$x))
        rangeY=range(c(rangeY,d$y))
        rangeXFull=range(c(rangeXFull,dFull$x))
        rangeYFull=range(c(rangeYFull,dFull$y))
        rangeYe=range(c(rangeYe,de$y))
        rangeXe=range(c(rangeXe,de$x))
        densities[[n]]=list(id=name_id,x=d$x,y=d$y,ex=de$x,ey=de$y,xFull=dFull$x,yFull=dFull$y, ntip=ntip, all.rates=BAMM_rates,rates.tips=BAMM_rates[tips_edges])
        
        
        ### and save all this
        
        do=(sum(data_BB$name==name_tree)==0)
        
        line_df=data.frame(name=name_tree,ES_Nshifts=ES_Nshifts,ES_LP=ES_LP,ntips=ntip,ngen=length(edata$numberEvents),
                           var_BAMM=var_BAMM,var_log_BAMM=var_log_BAMM,BAMM_ext=mean(BAMM_ext_rates),
                           var_BAMM_div=var_BAMM_div,var_log_BAMM_div=var_log_BAMM_div,
                           BD_birth=BD_rates[1],BD_death=BD_rates[2],
                           var_DR=var_DR,var_log_DR=var_log_DR)
        if(!do){
          data_BB[data_BB$name==name_tree ,]=line_df
          
        }else{
          data_BB=rbind(data_BB,line_df)
        }
        
      }}
    save(data_BB,rangeX, rangeY,rangeXFull,rangeYFull,rangeYe, rangeXe,densities,all.ratesFull,file=("~/Documents/ClaDS/Revisions/cleaned_data/BirdFamily/data_birds_BAMM.Rdata"))
    
  }
}


# # # # # # # # # # # # # # # # # # # # # # # # #
#  contribution of intraclade variance in rates #
# # # # # # # # # # # # # # # # # # # # # # # # #

# BAMM
load(paste0(path,"Simulations/data_birds_BAMM.Rdata"))
all_rates=c()
mean_rates=c()
for (i in 1:length(densities)){
  all_rates=c(all_rates,log(densities[[i]]$all.rates))
  mean_rates=c(mean_rates,rep(mean(log(densities[[i]]$all.rates),length(densities[[i]]$all.rates))))
}

var(all_rates)
var(mean_rates)
1-var(mean_rates)/var(all_rates)

# ClaDS2
load(paste0(path,"Simulations/data_birds.Rdata"))
all_rates=c()
mean_rates=c()
for (i in 1:length(densities)){
  all_rates=c(all_rates,(densities[[i]]$all.rates))
  mean_rates=c(mean_rates,rep(mean((densities[[i]]$all.rates),length(densities[[i]]$all.rates))))
}

var(all_rates)
var(mean_rates)
1-var(mean_rates)/var(all_rates)


#### check the MCMC sampler implementations ####

# # # # # # # # # # # # # # # # # # # # # # # # #
#                    ClaDS2                     #
# # # # # # # # # # # # # # # # # # # # # # # # #


prepare_ClaDS_prior=function(tree,sample_fraction, Nchain=3,nlambda = 1000,nt = 30,l0=0.1,s0=1,model_id="ClaDS2",decreasing.var=NULL, res_ClaDS0=NULL, silent=T, factor=0.1){
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
  likelihood_relative=createLikelihood_ClaDS0(tree)$ll
  alphaP1=1
  betaP1=log(1.1)
  P1=alphaP1*log(betaP1)-log(gamma(alphaP1))
  
  posterior <- function(x,former=NULL) {
    if(x[1]<0|any(x[c(1)]<0)){
      return(list(LL=-Inf,LP=-Inf,P=0))
    }else{
      LL=likelihood_relative(sigma=(x[1]),alpha=exp(x[2]),c(exp(x[3]),exp(x[2]+x[-(1:3)])) )
      if (is.nan(LL)) LL=-Inf
      Pr=P1-2*log(x[1])-betaP1/(x[1])
      if(x[2]<log(0.4) | x[2]>log(2)){Pr=-Inf}
      if(x[3]<log(1e-5) | x[3]>log(1000)){Pr=-Inf}
      return(list(LL=LL,LP=LL+Pr,Pr=Pr))}}
  
  if(is.null(res_ClaDS0)){
    start=lapply(1:Nchain,function(i){c(sigma=(s0),alpha=log(1),lambda=(c(rnorm(nedges+1,log(l0),0)))) })#rel.to.abs(tree,c(log(0.001),rep(log(A[i]),nedges)))) })
    if(is.null(decreasing.var)) decreasing.var=c(rep(1e-3,2),rep(1e-3,2),rep(1e-3,nedges))
  }else{
    start=lapply(1:Nchain,function(i){x=res_ClaDS0[[i]]$chain[nrow(res_ClaDS0[[i]]$chain),1:(nedges+3)]
    return(c(log(x[1]),x[2],runif(1),relToAbs(x[-(1:2)])))})
    decreasing.var=sapply(1:3,function(i){res_ClaDS0[[i]]$finetune})
    decreasing.var=rowMeans(decreasing.var)
    decreasing.var=0.1*c(decreasing.var[1:2],0.1,decreasing.var[-(1:2)])
  }
  
  npar=length(start[[1]])
  g=5 #mean number of parameter updated at each iteration
  G=exp(-(1:(npar))/g)
  testGenerator <- proposalGeneratorFactoryDE_gibbs(decreasing.var=decreasing.var,alpha=log(0.1)/200,
                                                    var=c(1e-4,1e-4,rep(1e-4,2),rep(1e-4,nedges)),
                                                    proba.gibbs=c(G))
  
  sampler=mcmcSamplerDE_gibbs(posterior,startvalue = start,proposalGenerator = testGenerator,Nchain = Nchain,consoleupdates = 1000)
  return(sampler)
}


replicate=0 # in 0:999
set.seed(replicate)
sigma=rinvgamma(n = 1,shape = 1, scale = log(1.1))
alpha=exp(runif(1,log(0.4),log(2)))
l0=exp(runif(1,log(1e-5),log(1000)))

obj= sim_ClaDS( lamb_par=l0,      
                  mu_par=0,
                  sigma=sigma,  
                  lamb_shift=alpha, 
                  condition="taxa", 
                  taxa.stop = 50, 
                  prune.extinct = T) 

tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]
  
sampler = prepare_ClaDS_prior(tree=tree,          # the phylogeny
                                sample_fraction=1,  # the fraction of species of the clade present in the reconstructed phylogeny
                                Nchain = 3,         # number of MCMC chains
                                nlambda = 1000,     # number of lambdaspace steps in likelihood function
                                nt = 30,            # number of time steps in likelihood function
                                model_id="ClaDS2",  # use model_id = "ClaDS1" for constant extinction rate
                                res_ClaDS0 = NULL, silent = T)  # the output of ClaDS0 to use as a startpoint. If NULL (the default)a random startpoint is used 
  
gelman=10
id_g=0
npar=nrow(tree$edge)+3
while(gelman>1.05){
    sampler = fit_ClaDS(           
      mcmcSampler = sampler,       
      iterations = 100000,        
      thin = 5000,  
      nCPU = 3)     

nr=nrow(sampler$chains[[1]])
print(nr)
rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler$chains[[j]][-(1:ceiling(nr/10)),-c((npar+1):(npar+3))])}))
gelman=try(gelman.diag(rep))
  if(! inherits(gelman,"try-error")){
     id_g=which.max(gelman$psrf[,1])
    gelman=max(gelman$psrf[,1])
  }else{
     gelman=10
  }
  print(c(gelman, id_g))
}

sampler$chains[[1]]=mcmc(sampler$chains[[1]][rows,])
sampler$chains[[2]]=mcmc(sampler$chains[[2]][rows,])
sampler$chains[[3]]=mcmc(sampler$chains[[3]][rows,])

save(sampler,file=paste0("checkSampler_50_Cl2_",seed,".Rdata"))

load("./check_sampler/ClaDS2/prior_draws_Cl2.Rdata")
param_Cl2$alpha = log(param_Cl2$alpha)
param_Cl2$l0 = log(param_Cl2$l0)
results_Cl0 = list()

files = list.files("./", pattern = "*.Rdata", full.names = T)

results_Cl2=list()
for (i in 1:(length(files)-1)){
  load(files[i])
  nr=nrow(sampler$chains[[1]])
  rows=1:nr
  results_Cl2[[i]] =  list(mcmc(sampler$chains[[1]]),
                                   mcmc(sampler$chains[[2]]),
                                   mcmc(sampler$chains[[3]]))
}

calibrationTest(posteriorList = results_Cl2, 
                priorDraws= param_Cl2[,1:3], 
                whichParameters = c(1:3),
                thin = 1)


# # # # # # # # # # # # # # # # # # # # # # # # #
#                    ClaDS0                     #
# # # # # # # # # # # # # # # # # # # # # # # # #

run_ClaDS0_prior=function(tree,name,pamhLocalName,proposalKernel="bactrian",
                          iteration=10000000, thin=20000,update=1000, adaptation=10000,
                          seed=NULL, nCPU=3, param=c()){
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
  
  start=c(1,0,1,rnorm(nedges,0,log(1.00000)))
  
  target <- function(x) {if(x[1]<0|any(x[c(1)]<0)){
    return(list(LL=-Inf,LP=-Inf,P=0))
  }else{
    LL=likelihood_relative(sigma=(x[1]),alpha=exp(x[2]),c(exp(x[3]),exp(x[2]+x[-(1:3)])) )
    if (is.nan(LL)) LL=-Inf
    Pr=P1-2*log(x[1])-betaP1/(x[1])
    if(x[2]<log(0.4) | x[2]>log(2)){Pr=-Inf}
    if(x[3]<log(1e-5) | x[3]>log(1000)){Pr=-Inf}
    return(list(LL=LL,LP=LL+Pr,Pr=Pr))}}
  
  model=create_model(start_val = start,posterior = target,proposalKernel = proposalKernel,tuning = 0.1,pamhType = pamhType)
  
  ptm <- proc.time()
  sampler=mclapply(1:3,function(j){ autoMetropolisGibbs(model, iterations=iteration, consoleupdates=1000000, 
                                                        thin=thin, autoOptimize=TRUE, filename=paste(pamhLocalName,j,sep=""),
                                                        update=update, adaptation=adaptation, verbose=F)},mc.cores = nCPU)
  
  rep=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),-c((npar+1):(npar+3))])}))
  gelman=try(gelman.diag(rep))
  if(! inherits(gelman,"try-error")) {gelman=max(gelman$psrf[,2])} else {gelman=2}
  print(gelman)
  
  while(gelman>1.05){
    sampler2=mclapply(1:3,function(j){
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
  
  if(! is.null(name)){save(tree, param, sampler,file = paste(name))}
  return(sampler)
}

replicate=0 # in 0:999
set.seed(replicate)

sigma=rinvgamma(n = 1,shape = 1, scale = log(1.1))
alpha=exp(runif(1,log(0.4),log(2)))
l0=exp(runif(1,log(1e-5),log(1000)))


obj= sim_ClaDS( lamb_par=l0,      # initial speciation rate
                mu_par=0,          # turnover rate (extinction/speciation)
                sigma=sigma,           # standard deviation of the new rates law
                lamb_shift=alpha,      # trend parameter alpha
                condition="taxa",    # the stoppping condition to use (can also be time)
                taxa.stop = 10,      # the number of tips to simulate (if condition=="taxa")
                prune.extinct = T)   # should extincti taxa be pruned from the result? (default to T)


# the function returns a list with the tree and the associated rates. Here is how you get to them :
tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

run_ClaDS0_prior(tree,name=paste0("checkSampler_Cl0_",replicate,".Rdata"),
                 pamhLocalName=paste0("Cl0_",replicate,"_"),proposalKernel="bactrian",
                 iteration=500000, thin=2000,update=1000, adaptation=5000,
                 seed=replicate, nCPU=3, param=list(sigma=sigma,alpha=alpha,
                                                    l0=l0,speciation_rates=speciation_rates,seed=replicate))

# only keep 100 iteration per chain to save memory space 
load(paste0("checkSampler_Cl0_",replicate,".Rdata"))
nr=nrow(sampler[[1]]$chain)
rows=ceiling(seq(1,nr-1,length.out = min(nr-1,100)))
sampler[[1]]$chain=mcmc(sampler[[1]]$chain[rows,])
sampler[[2]]$chain=mcmc(sampler[[2]]$chain[rows,])
sampler[[3]]$chain=mcmc(sampler[[3]]$chain[rows,])
save(tree, param, sampler,paste0("checkSampler_Cl0_",replicate,".Rdata"))
  
   
# analyse the results

load("./prior_draws_Cl2.Rdata")
param_Cl2$alpha = log(param_Cl2$alpha)
param_Cl2$l0 = log(param_Cl2$l0)
results_Cl0 = list()

files = list.files("./", pattern = "*.Rdata", full.names = T)
paramDraws_Cl0 = matrix(ncol = 3, nrow = length(files))
colnames(paramDraws_Cl0) = colnames(param_Cl2)[1:3]
for (i in 1:(length(files))){
  load(files[i])
  print(paste(i, files[i]))
  nr=nrow(sampler[[1]]$chain)
  rows=1:nr#ceiling(seq(1,nr-1,length.out = min(nr-1,50)))
  paramDraws_Cl0[i,1:3] = c(param[[1]], param[[2]], param[[3]])
  results_Cl0[[i]] = list(mcmc(sampler[[1]]$chain),
                          mcmc(sampler[[2]]$chain),
                          mcmc(sampler[[3]]$chain))
}

paramDraws_Cl0 = as.data.frame(paramDraws_Cl0)
paramDraws_Cl0$alpha = log(paramDraws_Cl0$alpha)
paramDraws_Cl0$l0 = log(paramDraws_Cl0$l0)

calibrationTest(posteriorList = results_Cl0, 
                priorDraws= paramDraws_Cl0, 
                whichParameters = c(1:3),
                thin = 1)

