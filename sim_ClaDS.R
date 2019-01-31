library(geiger)
library(ape)
library(phytools)

sim_ClaDS <- function (lamb_par, mu_par,theta=1,
                                       f.lamb=function(x,y){y}, f.mu=function(x,y){y}, 
                                       lamb_shift=1,mu_shift=1,
                                       new_lamb_law="lognormal*shift",new_mu_law="turnover",alpha=0, 
                                       sigma=0.1,lamb_max=1,lamb_min=0,mu_min=mu_par,mu_max=mu_par, 
                                       time.stop = 0, 
                                       sigma_mu=0,taxa.stop = Inf, return.all.extinct=TRUE, 
                                       prune.extinct=TRUE, condition="time",nShiftMax=Inf,print=F,
                                       maxRate=Inf)

# theta probability to have a shift at a speciation event (independant for the two new lineages)
# lamb_par parameters of the birth rate function, i.e. initial speciation rate if you keep f.lamb the default
# mu_par parameters of the birth rate function, i.e. initial extinction rate if you keep f.mu the default
# f.lamb , f.mu could theoretically be used to have rates variing through time but I have to modify the function first, so for now it must stay the default
# condition="time" the process stop after a certain time or number of taxa (if "taxa") is reached (or before if it goes extinct and return.all.extinct=TRUE)
# time.stop, taxa.stop time or number of tips before the process is stoped
  
# new_lamb_law the probability law from which new speciation rates are drawn. Can be "uniform", "normal", "lognormal", "normal*t" (normal with standard 
# deviation proportional to the length of the branch), "lognormal*shift", "normal+shift", "normal*shift" (lognormal or 
# normal with mode the previous speciation rate * (or +) lamb_shift)
  
# new_mu_law the probability law from which new death rates are drawn. Can be "uniform", "normal", "lognormal","normal*t" (normal with standard 
# deviation proportional to the length of the branch), "diversify" (constant diversification rate) or "turnover" (constant turnover rate)
# Yet I never used it with varying extinction rate, so I am not sure there would be no problem. To have constant extinction rate 
# set new_mu_law="uniform", mu_min=mu_max=mu_par

# relative_death For new_mu_law="diversify" or "turnover", diversification or turnover rate. Not used in other cases
# sigma=0.1 standard deviation for new_lamb_law = "normal", "lognormal", "normal*t", "lognormal*t", "lognormal*shift", "normal+shift", "normal*shift"
# sigma_mu standard deviation for new_mu_law = "normal", "lognormal", "normal*t"

# lamb_max,lamb_min,mu_min,mu_max, limits of the uniform laws
# return.all.extinct if FALSE the fuction is runned until we get a tree that don't become extinct before the stopping condition is reached 
# prune.extinct are extinct taxa removed from the phylogeny ?

# nShiftMax : if this number of shifts in the phylogeny is reached, theta is set to 0

# The function returns a list list with:
# list$tree the resulting tree
# list$rates the rate category for each branch
# list$lamb and list$mu the speciation and extinctio rate of each rate category (so that speciation rate of each branch is obtained 
# with list$lamb[list$rates])
# list$times and list$nblineages speciation times and corresponding number of lineages 

{

  relative_death=mu_par
  if(new_mu_law=="turnover"){ mu_par=mu_par*lamb_par}
  
  while (1) {
    
    lamb=lamb_par
    mu=mu_par
    rates=c()
    nShift=0
    tooHigh=F
    
    nblineages<-c(1)
    times<-c(0)
    b<-lamb
    d<-mu
    dt <- rexp(1,(b + d))
    #print(dt)
    t<-dt
    
    if ((t >= time.stop| 1>taxa.stop) & condition=="time") {
      t <- time.stop
    	alive<-1
    	rates<-c(1)
    	times<-c(times,t)
    	nblineages<-c(nblineages,1)
    	break
    }
    
    if (taxa.stop==1 & condition=="taxa") {
      alive<-1
      rates<-c(1)
      break
    }
            	
    r <- runif(1)
    
    if (r>b/(b + d)){
      #print("die")
      times<-c(times,dt)
      nblineages<-c(nblineages,0)
      alive<-rep(FALSE,1)
    }else{
      u=runif(2)
      for(i in 1:2){
        if(u[i]<theta & nShift<nShiftMax){
          nShift=nShift+1
          rates=c(rates,length(lamb)+1)
          if(new_lamb_law=="uniform"){
            lamb=c(lamb,runif(1,min=lamb_min,max=lamb_max))
          }else if (new_lamb_law=="normal"){
              new_lambda=rnorm(1,mean=lamb[[1]],sd=sigma)
              while(new_lambda<0){
                new_lambda=rnorm(1,mean=lamb[[1]],sd=sigma)
              }
              lamb=c(lamb,new_lambda)
          }else if (new_lamb_law=="lognormal"){
              lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[1]]),sdlog = sigma))
          }else if (new_lamb_law=="lognormal*shift"){
            lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[1]]*lamb_shift),sdlog = sigma))
          }else if (new_lamb_law=="lognormal*t"){
            lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[1]]),sdlog = sigma*t))
          }else if (new_lamb_law=="logbrownian"){
            lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[1]]),sdlog = sigma*sqrt(t)))
          }else if (new_lamb_law=="normal*t"){
            new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma*t)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma*t)
            }
            lamb=c(lamb,new_lambda)
          }else if (new_lamb_law=="normal+shift"){
            new_lambda=rnorm(1,mean=lamb[[1]]+lamb_shift,sd=sigma)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb[[1]]+lamb_shift,sd=sigma)
            }
            lamb=c(lamb,new_lambda)
          }else if (new_lamb_law=="normal*shift"){
            new_lambda=rnorm(1,mean=lamb[[1]]*lamb_shift,sd=sigma)
            while(new_lambda<0){
              new_lambda=rnorm(1,mean=lamb[[1]]*lamb_shift,sd=sigma)
            }
            lamb=c(lamb,new_lambda)
          }
          if(new_mu_law=="uniform"){
            mu=c(mu,runif(1,min=mu_min,max=mu_max))
          }else if (new_mu_law=="normal"){
            new_mu=rnorm(1,mean=mu[[1]],sd=sigma)
            while(new_mu<0){
              new_mu=rnorm(1,mean=mu[[1]],sd=sigma)
            }
            mu=c(mu,new_mu)
          }else if (new_mu_law=="diversify"){
            new_mu=lamb[length(lamb)]-relative_death
            mu=c(mu,max(0,new_mu))
          }else if (new_mu_law=="turnover"){
            
            new_mu=lamb[length(lamb)]*relative_death
            mu=c(mu,new_mu)
          }else if (new_mu_law=="lognormal"){
            mu=c(mu,rlnorm(1, meanlog = log(mu[[1]]),sdlog = sigma_mu))
          }else if (new_mu_law=="lognormal*shift"){
            mu=c(mu,rlnorm(1, meanlog = log(mu[[1]]*mu_shift),sdlog = sigma_mu))
          }else if (new_mu_law=="normal*t"){
            new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu*t)
            while(new_mu<0){
              new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu*t)
            }
            mu=c(mu,new_mu)
          }
        }else{
          rates=c(rates,1)
        }
      }

    	edge <- rbind(c(1, 2), c(1, 3))
    	edge.length <- rep(NA, 2)
    	stem.depth <- rep(t, 2)
    	alive <- rep(TRUE, 2)
    	times<-c(times,dt)
    	nblineages<-c(nblineages,sum(alive))
    	next.node <- 4}
    	
    	repeat {
    	  if(print){
    	    cat("\r",sum(alive),"of",taxa.stop," ; ",t, "of" ,time.stop," ; ",max(lamb), "of" ,maxRate,"\r")
    	  }
    	  if (sum(alive) == 0)	 break
    	  if (sum(alive)>=taxa.stop){#} & condition=="taxa") {
    	    if(length(lamb)<=1){
    	      b=max(lamb,.Machine$double.eps)
    	      d=mu
    	      totalrate=sum(alive)*(b+d)
    	    }else{
    	      b<-sapply(1:length(lamb),function(x){max(lamb[[x]],.Machine$double.eps)})
    	      d<-sapply(1:length(mu),function(x){mu[[x]]})
    	      totalrate=sum(b[rates][alive])+sum(d[rates][alive])}
    	    
    	    dt <- rexp(1, totalrate)
    	    t <- t + dt
    	    
    	    break
    	  }else{
    		
    	  if(length(lamb)<=1){
    	    b=lamb
    	    d=mu
    	    totalrate=sum(alive)*(b+d)
    	    live_rates=rates[alive]
    	 }else{
    	   live_rates=rates[alive]
    	   b=lamb
    	   d=mu
         totalrate=sum(b[live_rates])+sum(d[live_rates])}

        dt <- rexp(1, totalrate)
    		t <- t + dt
    		
    		if (t >= time.stop & condition=="time") {
    		  t <- time.stop
    			times<-c(times,t)
    			nblineages<-c(nblineages,sum(alive))
    			break
    		}
    		
    		if (any(lamb>maxRate)) {
    		  tooHigh=T
    		  break
    		}
    		
    		r <- runif(1)
    		s=0
    		continue=T
    		if(theta<1) {
    		  Walive=unique(live_rates)
    		}else{
    		    Walive=live_rates[order(b[live_rates],decreasing = T)]
    		  }
    		k=1
    		sumBirth=sum(b[live_rates])/totalrate
    		i=Walive[k]
    		if(r<=sumBirth){
    		while(i<=length(b) & continue){
    		  if(theta<1){
    		    s=s+b[i]*sum((live_rates==i))/totalrate
    		  }else{
    		    s=s+b[i]/totalrate
    		    }
    		  if(r<=s){
    		    continue=F
    		    #print("speciation")
    		    if(theta<1){if(length(which(rates==i & alive))>1){
    		      random_lineage = sample(which(rates==i & alive),1)
    		    }else{
    		      random_lineage = which(rates==i & alive)
    		    }}else{
    		      random_lineage = which(rates==i)
    		    }
    		    parent <- edge[random_lineage, 2]
    		    alive[random_lineage] <- FALSE
    		    edge <- rbind(edge, c(parent, next.node), c(parent,next.node + 1))
    		    next.node <- next.node + 2
    		    alive <- c(alive, TRUE, TRUE)
    		    stem.depth <- c(stem.depth, t, t)
    		    #print(stem.depth)
    		    x <- which(edge[, 2] == parent)
    		    edge.length[x] <- t - stem.depth[x]
    		    edge.length <- c(edge.length, NA, NA)
    		    #print(edge.length)
    		    times<-c(times,t)
    		    nblineages<-c(nblineages,sum(alive))
    		    u=runif(2)
    		    for(j in 1:2){
    		      if(u[j]<theta & nShift<nShiftMax){
    		        nShift=nShift+1
    		        rates=c(rates,length(lamb)+1)
    		        if(new_lamb_law=="uniform"){
    		          lamb=c(lamb,runif(1,min=lamb_min,max=lamb_max))
    		        }else if (new_lamb_law=="normal"){
    		          new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma)
    		          }
    		          lamb=c(lamb,new_lambda)
    		        }else if (new_lamb_law=="lognormal"){
    		          lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[i]]),sdlog = sigma))
    		        }else if (new_lamb_law=="lognormal*shift"){
    		          #mean lamb[[i]] iff lamb_shift==exp(-sigma^2/2)
    		          lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[i]]*lamb_shift),sdlog = sigma))
    		        }else if (new_lamb_law=="lognormal*t"){
    		          lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[i]]),sdlog = sigma*edge.length[x]))
    		        }else if (new_lamb_law=="logbrownian"){
    		          lamb=c(lamb,rlnorm(1, meanlog = log(lamb[[i]]),sdlog = sigma*sqrt(edge.length[x])))
    		        }else if (new_lamb_law=="normal*t"){
    		          new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma*edge.length[x])
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb[[i]],sd=sigma*edge.length[x])
    		          }
    		          lamb=c(lamb,new_lambda)
    		        }else if (new_lamb_law=="normal+shift"){
    		          new_lambda=rnorm(1,mean=lamb[[i]]+lamb_shift,sd=sigma)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb[[i]]+lamb_shift,sd=sigma)
    		          }
    		          lamb=c(lamb,new_lambda)
    		        }else if (new_lamb_law=="normal*shift"){
    		          new_lambda=rnorm(1,mean=lamb[[i]]*lamb_shift,sd=sigma)
    		          while(new_lambda<0){
    		            new_lambda=rnorm(1,mean=lamb[[i]]*lamb_shift,sd=sigma)
    		          }
    		          lamb=c(lamb,new_lambda)
    		        }
    		        if(new_mu_law=="uniform"){
    		          mu=c(mu,runif(1,min=mu_min,max=mu_max))
    		        }else if (new_mu_law=="normal"){
    		          new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu)
    		          while(new_mu<0){
    		            new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu)
    		          }
    		          mu=c(mu,new_mu)
    		        }else if (new_mu_law=="turnover"){
    		          new_mu=lamb[length(lamb)]*relative_death
    		          mu=c(mu,new_mu)
    		        }else if (new_mu_law=="diversify"){
    		          new_mu=lamb[length(lamb)]+relative_death
    		          mu=c(mu,new_mu)
    		        }else if (new_mu_law=="lognormal"){
    		          mu=c(mu,rlnorm(1, meanlog = log(mu[[i]]),sdlog = sigma_mu))
    		        }else if (new_mu_law=="lognormal*shift"){
    		          mu=c(mu,rlnorm(1, meanlog = log(mu[[i]]*mu_shift),sdlog = sigma_mu))
    		        }
    		        else if (new_mu_law=="normal*t"){
    		          new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu*edge.length[x])
    		          while(new_mu<0){
    		            new_mu=rnorm(1,mean=mu[[i]],sd=sigma_mu*edge.length[x])
    		          }
    		          mu=c(mu,new_mu)
    		        }
    		      }else{
    		        rates=c(rates,i)
    		      }}
    		    
    		    #print(times)
    		    #print(nblineages)
    		  }
    		  k=k+1
    		  i=Walive[k]
    		}}else{
    		  s=sumBirth
    		}
    		k=1
    		i=Walive[k]
    		while(i<=length(b) & continue){
    		  if(theta<1){
    		    s=s+d[i]*sum((live_rates==i))/totalrate
    		  }else{
    		    s=s+d[i]/totalrate
    		    }
    		  
    		  if(r<=s){    							
    		    if(theta<1){if(length(which(rates==i & alive))>1){
    		      random_lineage = sample(which(rates==i & alive),1)
    		    }else{
    		      random_lineage = which(rates==i & alive)
    		    }}else{
    		      random_lineage = which(rates==i)
    		    }
    		    continue=F
    		    edge.length[random_lineage] <- t - stem.depth[random_lineage]
    		    alive[random_lineage] <- FALSE
    		    times<-c(times,t)
    		    nblineages<-c(nblineages,sum(alive))}
    		  k=k+1
    		  i=Walive[k]
    		}
    	}}
    	
    	
  if (return.all.extinct == TRUE | sum(alive) > 0) {
    #print("return.tree")
    break
    }
    }
    							
	if ((sum(alive)==0 & prune.extinct) | (length(nblineages)==2 & !prune.extinct)) {obj<-NULL; root_length=t} #
 	else if (sum(alive)==1 & prune.extinct) {obj<-list(nbTaxa=1,"maxRate"=tooHigh); root_length=t}
 	else {
 	  edge.length[alive] <- t - stem.depth[alive]
 	  n <- -1
 	  for (i in 1:max(edge)) {
 	  if (any(edge[, 1] == i)) {
 	  edge[which(edge[, 1] == i), 1] <- n
 	  edge[which(edge[, 2] == i), 2] <- n
 	  n <- n - 1
 	  }
 	}
    
 	edge[edge > 0] <- 1:sum(edge > 0)
  tip.label <- 1:sum(edge > 0)
  mode(edge) <- "character"
  mode(tip.label) <- "character"
  obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label)

  class(obj) <- "phylo"
  obj <- old2new.phylo(obj)
   rep=rigth.order(obj,rates)
   obj=rep$tree
   rates=rep$rates
  if (prune.extinct){
    rep=prune.extinct.with.rates(obj,rates)
    obj=rep$tree
    rates=rep$rates
  }
   root_length=t-max(node.depth.edgelength(obj))

 	}
  return(list("tree"=obj,"times"=times,"nblineages"=nblineages,"rates"=rates,"lamb"=lamb,"mu"=mu,"maxRate"=tooHigh,"root_length"=root_length))
  
  }

prune.extinct.with.rates=function(phy,rates,extinct=NULL)
  #used in the simulation function to remove extinct taxa withour loosing the rate information
  
{
  obj=list(tree=phy,rates=rates)
  if(is.null(extinct)) extinct=is.extinct(obj$tree,tol=max(obj$tree$edge.length)/100000)
  nodes=extinct
  if(length(extinct)>0){
    for(i in 1:length(extinct)){
      edge=which(obj$tree$edge[,2]==which(obj$tree$tip.label==extinct[i]))
      if(obj$tree$edge[edge,1]==(obj$tree$Nnode+2)){
        edge=which(obj$tree$edge[,1]==(obj$tree$Nnode+2))
      }else{
        edge=c(which(obj$tree$edge[,1]==obj$tree$edge[edge,1]))
      }
      obj$rates=obj$rates[-edge]
      obj$tree=drop.tip(obj$tree,extinct[i])
    }}
  return(obj)
}

rigth.order=function(phy,rates){
  n=phy$Nnode+1
  root=which(sapply(1:max(phy$edge),function(x){!(x %in% phy$edge[,2])}))
  next_node=c(root)
  order=rep(0,2*n-1)
  is.tip=c()
  i=1
  while(length(next_node)>0){
    order[next_node[1]]=i
    offspring=phy$edge[phy$edge[,1]==next_node[1],2]
    next_node=c(offspring,next_node[-1])
    if(length(offspring)==0) is.tip=c(is.tip,i)
    i=i+1
  }
  phy$edge[,1]=order[phy$edge[,1]]
  phy$edge[,2]=order[phy$edge[,2]]
  order=order(phy$edge[,2])
  phy$edge=phy$edge[order,]
  if(inherits(rates,"list")){
    rates=lapply(rates,function(x){x[order]})
  }else{
    rates=rates[order]
  }
  phy$edge.length=phy$edge.length[order]
  itip=1
  inode=n+1
  newnames=c()
  for(i in 1:(2*n-1)){
    if(i %in% is.tip) {
      newnames=c(newnames,itip)
      itip=itip+1
    }else{
        newnames=c(newnames,inode)
        inode=inode+1
      }
  }
  phy$edge=matrix(newnames[phy$edge],ncol=2)
  return(list(tree=phy,rates=rates))
}

mean.rate=function(tree,rate){
  nh=nodeHeights(tree)
  times=sort(unique(c(nh[,1],nh[,2])))
  meanRate=c()
  nb=0
  keep.i=c()
  for(i in 1:(length(times)-1)){
    t=times[i]
    which.edge=(nh[,1]<=t & nh[,2]>t)
    if(sum(which.edge)>nb){
      nb=sum(which.edge)
    meanRate=c(meanRate,sum(rate[which.edge])/sum(which.edge))
    keep.i=c(keep.i,i)
    }
  }
  keep.i=c(keep.i,length(times))
  return(rbind(times[keep.i],c(meanRate[1],meanRate)))
}
# 
# X=mean.rate(tree,true.rate)
# plot.default(t(X), xaxs = "r", yaxs = "r",  type = "S")

get_rates <- function(phylo, lambda, Ancestors=NULL){
  
  
  nbtips = Ntip(phylo)
  edge = phylo$edge
  nedges=(2*(nbtips - 1))
  edge.length = phylo$edge.length
  
  if(is.null(Ancestors)){
    type = rep(NA, Ntip(phylo))
    parents = rep(NA, Ntip(phylo))
    
    for (i in 1:nedges) {
      parent=which(phylo$edge[,2]==phylo$edge[i,1])
      if(length(parent)==0){
        type[i] = 1
        parents[i] = NA
      }else if (!(phylo$edge[i,2]<=nbtips)){
        type[i] = 2    
        parents[i] = parent
      }else{
        type[i] = 3       
        parents[i] = parent
      }
    }
    
    roots = which(type ==1)
    internal = which(type ==2)
    parentsInternal = parents[type == 2]
    internalAndRoots = which(type <3)
    internalAndTerminal = which(type >1)
    parentsInternalAndTerminal = parents[type >1]
    
    ancestors=list()
    for(i in 1:nedges){
      if(type[i]==1){
        ancestors[[i]]=-1
      }else{
        ancestors[[i]]=c(parents[i],ancestors[[parents[i]]][ancestors[[parents[i]]]>0])
      }
    }
  }else{ancestors=Ancestors}
  
  
  lambda2=c(lambda[1],sapply(1:(nedges),function(i){lambda[i+1]+lambda[1]+sum(lambda[ancestors[[i]][ancestors[[i]]>0]+1])}))
  return(lambda2)
  
}


plot.with.rate=function(phylo,rate1,rate2=NULL,same.scale=T,main=NULL,lwd=1,log=F){
  Colors = colorRampPalette(c("steelblue2","paleturquoise3","palegreen2","yellow2","salmon1","darkorange", "red","red4"))( 100 ) 
  if(is.null(rate2)){
    if(log) rate1=log(rate1)
    if(isTRUE(all.equal(rep(as.numeric(rate1[1]),length(rate1)),as.numeric(rate1)))){
      col=rep(1,length(rate1))
      plot(phylo, edge.color = Colors[col], show.tip.label = F,main=main,edge.width =lwd)
      if(log){
        image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T)
      }else{
        image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
      }
    }else{
      col = round( (rate1 - min(rate1)) / diff(range(rate1))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = F,main=main,edge.width =lwd)
      if(log){
        min=min(rate1)
        max=max(rate1)
        m10=floor(min/log(10))
        M10=ceiling(max/log(10))
        if((M10-m10)<4){
          ticks=c(1,2,5)
        }else{
          ticks=1
        }
        ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
        lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
        if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
        image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
      }else{
        image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
      }
    }
  }else{
    if(log){
      rate1=log(rate1)
      rate2=log(rate2)
    }
    if(same.scale){
      min=min(min(rate1),min(rate2))
      max=max(max(rate1),max(rate2))
      par(mfrow=c(1,2))
      col = round(( (rate1 - min) / (max-min))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
      col = round(( (rate2 - min) / (max-min))*99   )+1
      plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
      par(mfrow=c(1,1))
      if(log){
        m10=floor(min/log(10))
        M10=ceiling(max/log(10))
        if((M10-m10)<4){
          ticks=c(1,2,5)
        }else{
          ticks=1
        }
        ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
        lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
        if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
        # ticks=seq(min,max,length.out = 5)
        image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
      }else{
        image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T)
      }
    }else{
      par(mfrow=c(1,2))
      if(isTRUE(all.equal(rep(rate1[1],length(rate1)),rate1))){
        col=rep(1,length(rate1))
        plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        if(log){
          
          image.plot(z = c(exp(rate1[1]),2*exp(rate1[1])),col = Colors, horizontal=T,legend.only = T)
        }else{
          image.plot(z = c(rate1[1],2*rate1[1]),col = Colors, horizontal=T,legend.only = T)
        }
      }else{
        col = round(( (rate1 - min(rate1)) / (max(rate1)-min(rate1)))*99   )+1
        plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        if(log){
          min=min(rate1)
          max=max(rate1)
          m10=floor(min/log(10))
          M10=ceiling(max/log(10))
          if((M10-m10)<4){
            ticks=c(1,2,5)
          }else{
            ticks=1
          }
          ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
          lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
          if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
          image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
        }else{
          image.plot(z = as.matrix(rate1),col = Colors, horizontal=T,legend.only = T)
        }
      }
      if(isTRUE(all.equal(rep(rate2[1],length(rate2)),rate2))){
        col=rep(1,length(rate2))
        plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        if(log){
          image.plot(z = c(exp(rate2[1]),2*exp(rate2[1])),col = Colors, horizontal=T,legend.only = T)
        }else{
          image.plot(z = c(rate2[1],2*rate2[1]),col = Colors, horizontal=T,legend.only = T)
        }
      }else{
        col = round(( (rate2 - min(rate2)) / (max(rate2)-min(rate2)))*99   )+1
        plot(phylo, edge.color = Colors[col], show.tip.label = F,edge.width =lwd)
        if(log){
          min=min(rate2)
          max=max(rate2)
          m10=floor(min/log(10))
          M10=ceiling(max/log(10))
          if((M10-m10)<4){
            ticks=c(1,2,5)
          }else{
            ticks=1
          }
          ticks=as.vector(sapply(m10:M10,function(k){return(ticks*10^k)}))
          lt=length(ticks[ticks>exp(min) & ticks<exp(max)])
          if(lt<4){ticks=c(round(exp(min),max(0,-1*m10+(lt<2))),ticks,round(exp(max),max(0,-1*M10+1+(lt<2))))}
          image.plot(z = c(min,max),col = Colors, horizontal=T,legend.only = T,axis.args=list( at=log(ticks), labels=ticks))
        }else{
          image.plot(z = as.matrix(rate2),col = Colors, horizontal=T,legend.only = T)
        }
      }
    }
    par(mfrow=c(1,1))
  }
}
