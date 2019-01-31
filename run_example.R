#### 1. Loading the functions and packages ####

# chose your working directory : put the path of the downloaded Rcode folder below
setwd("~/Documents/ClaDS/")
library(TESS)

# loading the needed functions
source("likelihood_ClaDS0.R")
source("sim_ClaDS.R")
source("utils.R")
source("proposal.R")
source("run_ClaDS0.R")
source("ClaDS1_likelihood_functions.R")
source("ClaDS2_likelihood_functions.R")
source("fit_ClaDS.R")

# this one needs to be compiled ; for this you need to open a terminal, go to your 
# working directory, and enter
# R CMD SHLIB diversif_lognormal.c
dyn.load("diversif_lognormal.so")



#### 2. Simulating trees from the model #####

set.seed(1)

obj= sim_ClaDS( lamb_par=0.1,      # initial speciation rate
                mu_par=0.5,          # turnover rate (extinction/speciation)
                sigma=0.3,           # standard deviation of the new rates law
                lamb_shift=0.9,      # trend parameter alpha
                condition="taxa",    # the stoppping condition to use (can also be time)
                taxa.stop = 20,      # the number of tips to simulate (if condition=="taxa")
                prune.extinct = T)   # should extincti taxa be pruned from the result? (default to T)


# the function returns a list with the tree and the associated rates. Here is how you get to them :
tree = obj$tree
speciation_rates = obj$lamb[obj$rates]
extinction_rates = obj$mu[obj$rates]

# and the plotting function
plot.with.rate(tree,speciation_rates,
               log=T,        # should the rates be plotted on a log scale?
               lwd=3)        # width of the branches




#### 3. Infering the model parameters on a tree, for the pure birth version (ClaDS0) ####

# this function runs the mcmc until the gelman stopping criterion (that takes several chains and 
# measures the difference between the variance within and among chains) is reached. It should be below 
# 1.05, the functions prints the gelman factor every "iteration" iterations. The function reach 3 chains 
# at the same time, that can be ran in parallel (see nCPU)

# here I run it on the tree simulated above, which is very small

sampler = run_ClaDS0(tree=tree,        # the data
              name="example.Rdata",          # the name the results will be saved on (if NULL it won't be)
              nCPU=1,             # the number of CPUs to use (3 chains are run so it only makes sense to make it 1 or 3)               
              pamhLocalName = "local",   # the function is writing in a text file to make the execution quicker, this is the name of this file
              iteration=500000,   # number of iterations after which the gelman factor is computed and printed. The function stops if it is below 1.05
              thin=2000,           # number of iterations after which the chains state is recorded
              update=1000, adaptation=5000)  # options for the initial proposal adaptation phase



#### 4. Plotting the results (ClaDS0) ####

# if the exemple in 3. takes too long to run, it is saved in the attached file "example.Rdata"

load("example_ClaDS0.Rdata")
sampler_ClaDS0=sampler

ntips=tree$Nnode+1
nedges=2*tree$Nnode
npar=nedges+3


# prune the first iterations of the chains
chains=mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),])}))
plot(mcmc.list(lapply(1:3,function(j){mcmc(sampler[[j]]$chain[-(1:10),c(1:3,nedges+4)])})))

# the mcmc above is run on the log of the relative rates (rate/(parentRate*alpha)), 
# below I compute the log of the absolute rates within the chains to be then able to get
# to the MAPs (posterior maximas) of the rates
for(k in 1:length(chains)){
  for(l in 1:nrow(chains[[k]])){
    chains[[k]][l,4:npar]=get_rates(tree,c(chains[[k]][l,3],
                                            chains[[k]][l,4:npar]+chains[[k]][l,2]))[-1]
  }
}

# extract the MAPs (argmax of the posterior marginals)
MAPS=sapply(1:npar, function(i){D=density(c(chains[[1]][,i],
                                            chains[[2]][,i],
                                            chains[[3]][,i]))
return(D$x[which.max(D$y)])})
MAPS[2:npar]=exp(MAPS[2:npar])

print(paste0("sigma = ",MAPS[1]," ; alpha = ",MAPS[2]," ; lambda_0 = ", MAPS[3]))

plot.with.rate(tree, MAPS[4:npar],log=T,lwd=3)

# with the true simulated rates on the left panel :
plot.with.rate(tree, speciation_rates, MAPS[4:npar],log=T,lwd=3, same.scale = T)

# on such a small trees the inference is not very good. It gets better on larger trees.



#### 5. Infering the model parameters (constant turnover, ClaDS2) ####

# prepare the sampler
sampler = prepare_ClaDS(tree=tree,          # the phylogeny
                        sample_fraction=1,  # the fraction of species of the clade present in the reconstructed phylogeny
                        Nchain = 3,         # number of MCMC chains
                        nlambda = 1000,     # number of lambdaspace steps in likelihood function
                        nt = 30,            # number of time steps in likelihood function
                        model_id="ClaDS2",  # use model_id = "ClaDS1" for constant extinction rate
                        res_ClaDS0 = sampler_ClaDS0)  # the output of ClaDS0 to use as a startpoint. If NULL (the default)a random startpoint is used 
# and run it
sampler = fit_ClaDS(           
  mcmcSampler = sampler,       # an object as returned by prepare ClaDS or fit_ClaDS
  iterations = 1000,            # number of iterations
  thin = 50,                   # number of iterations after which the chains state is recorded
  nCPU = 3)                    # number of cores to use (between 1 and the number of chains)


# you can continue to run a fit
sampler = fit_ClaDS(mcmcSampler = sampler,iterations = 100,thin = 50, nCPU = 3)

# and plot the results
plotChains_ClaDS(sampler)
MAPS=getMAPS_ClaDS(sampler)
plot.with.rate(tree,speciation_rates,exp(MAPS[-(1:4)]),log=T,lwd=3)  # it should, of course, be run for many more iterations
