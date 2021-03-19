


rm(list=ls())
  
  library(gridExtra)
  
  
  #Parameters
  
  par(mfrow=c(1,1))
  
  r <- 0.0866       # Intrinsic rate of increase of susceptible
  mu <- 0.0002      # Natural mortality rate
  B <- 5       # Transmission coefficient
  alpha <- 0.5       # Mortality rate
  K <- 200          # Carrying Capacity
  phi <- 0.011      # Resistance loss rate
  delta <- 0.21      # Exposed to infectious period
  
  #starting conditions
  
  S.start <- 199
  E.start <- 0
  I.start <- 1
  R.start <- 0
  M.start <- 0
  g.start <- 0.005  # Initial resistance trait
  sigma <- 0.0005   # Trait variance
  
  t <- 10000
  

  
  
  
  ###############Plotting many deterministic runs
  B <- 5
  
  g.values <- seq(0,0.2,0.0025)
  sigma.values <- c(0)
  
  t <- 10000
  
  

  
  

  
  #################################Plotting many stochastic
  library(tidyverse)
  
  rho <- 0.08
  
  parameters <- c(sigma=sigma,r=r,mu=mu,B=B,alpha=alpha,K=K,rho=rho,delta=delta)
  state <- c(S = S.start,E = E.start,I = I.start,R = R.start)
  t <- 10000
  
  SIR.sto <- function(g.start,sigma,t,r,mu,B,alpha,K,phi,delta,S.start,E.start,I.start,R.start,M.start){
    XT <- 10000000 #potential number of events
    WT <- 0 #start time
    S <- S.start
    E <- E.start
    I <- I.start
    R <- R.start
    M <- M.start
    
    N <- S+E+I+R
    
    g <- rnorm(N,mean=g.start,sd=sigma)
    g[g<0] <- 0
    g[g>1] <- 1
    g.avg <- mean(g)
    
    SJ <- rep(0,XT) #Make list for time steps
    
    individuals <- rep(0,N[1]); #create vector for all individual animals 
    #and their potential classes
    if(S>0){individuals[1:S] <- 1}; #Assign S individuals into S
    if(I>0){individuals[(S+1):(S+I)] <- 2};
    #Assign R inidividuals
    if(R>0){individuals[(S+I+1):(S+I+R)] <- 3};
    #Assign R inidividuals
    
    for(i in 2:XT){
      
      if(length(which(individuals==1)) >=1){
        Strait <- which(individuals==1)
        gtrait <- sum(g[Strait])
      } else{gtrait <- 0}
      
      #Potential Events
      birth <- mu*N[i-1]
      nat <- mu*N[i-1]
      infect <- B*delta*S[i-1]*I[i-1]
      disease <- alpha*I[i-1]
      recover <- B*gtrait*I[i-1]
      
      sumrates <- birth+nat+infect+disease+recover
      if(sumrates <= 0) break
      
      #Determine time to next event (Sojourn time)
      SJ[i] <- rexp(1,rate=sumrates)
      WT <- c(WT,0)
      WT[i]=WT[i-1]+SJ[i]
      
      
      #Decide event to occur and update populations
      u<-runif(1,min=0,max=1)
      
      if(u<birth/sumrates){
        #Choose a susceptible parent
        if(N[i-1]<=1){  #sample cannot handle a list of 1 to choose from
          Mummy <- 1 #pick only one if 1 available
        } else{
          Mummy <- sample(1:length(individuals),1) #random sample
        } #end else
        #Give parent a baby
        individuals <- append(individuals,1); #New susceptible individual 
        g <- append(g,0);             #New trait for individuals
        #Determine if mutation occurs and assign new individual a trait
        mut <- runif(1,min=0,max=1); #random chance generator
        if(mut<rho){ #Mutation occurs
          g[length(g)] <- g[Mummy]+rnorm(1, mean=0, sd=sigma) #Use fixed standard deviation
        } else{g[length(g)] <- g[Mummy]} #Assign parent trait to newest individual
        g[g<0]=0
        g[g>1]=1
      } else{if(u<(birth+nat)/sumrates){
        #Choose an individual to off
        if(N[i-1]<=1){  
          Goner <- 1
        } else{
          Goner <- sample(1:length(individuals),1)
        } #end else
        individuals <- individuals[-Goner] #Kill it
        g <- g[-Goner]
      } else{if(u<(birth+nat+infect)/sumrates){
        #Choose a susceptible to infect
        if(S[i-1]<=1){  
          Diseased <- which(individuals==1,arr.ind=TRUE)
        } else{
          Diseased <- sample(which(individuals==1,arr.ind=TRUE),1)
          #Diseased <- sample(which(individuals==1,arr.ind=TRUE),1,prob=(1-g[Itrait]/gtrait))
        } #end else
        individuals[Diseased] <- 2 #Infected!
      } else{if(u<(birth+nat+infect+disease)/sumrates){
        if(I[i-1]<=1){  
          Released <- which(individuals==2,arr.ind=TRUE)
        } else{
          Released <- sample(which(individuals==2,arr.ind=TRUE),1)
        } #end else
        individuals <- individuals[-Released] #Kill it
        g <- g[-Released]
      } else{
        if(S[i-1]<=1){  
          EyeOfTiger <- which(individuals==1,arr.ind=TRUE)
        } else{
          EyeOfTiger=sample(which(individuals==1,arr.ind=TRUE),1,prob=g[Strait]/gtrait)
        } #end else
        individuals[EyeOfTiger] <- 3
      } #end else 4
      } #end else 3
      } #end else 2
      }#end else 1
      
      N <- cbind(N,length(individuals))
      S <- cbind(S,length(which(individuals==1)))
      I <- cbind(I,length(which(individuals==2)))
      R <- cbind(R,length(which(individuals==3)))
      g.avg <- append(g.avg,sum(g)/N[i])
      
      if(N[i] == 0) break
      if(WT[i] >= t) break
      
    }#end forloop
    
    evolve <- (max(g.avg,na.rm=TRUE)-g.avg[1])/sigma
    survive <- 0
    if(N[length(N)]>0) survive <- 1
    
    output <- c(N[length(N)],min(N),g.start,g.avg[1],max(g.avg,na.rm=TRUE),evolve,survive)
    return(output)
    
  } #end function
  
  
  
  library(doParallel)
  detectCores()
  c <- makeCluster(3)
  registerDoParallel(c)
  getDoParWorkers()
  
  
  rep <- 500
  
  results <- foreach(i = 1:length(g.values), .combine= cbind) %dopar% {
    sigma <- 0
    g.start <- rep(g.values[i], rep)
    ind.results <- mapply(SIR.sto,g.start,sigma=sigma,t=t,r=r,mu=mu,B=B,alpha=alpha,K=K,phi=phi,delta=delta,S.start=S.start,E.start=E.start,I.start=I.start,R.start=R.start,M.start=M.start)
    return(ind.results)
  }
  
  
  ind.results <- as.data.frame(t(results))
  colnames(ind.results) <- c("Nt","min.N","g","g.1","max.g","evolve","survive")
  ind.results$survive <- as.numeric(ind.results$survive)
  
  
  
  
  group.results <- NA
  for(i in 1:(length(ind.results[,1])/rep)){
    p.survive <- mean(ind.results$survive[(1+rep*(i-1)):(rep*i)])
    max.evolve <- mean(ind.results$evolve[(1+rep*(i-1)):(rep*i)])
    g.1 <- mean(ind.results$g.1[(1+rep*(i-1)):(rep*i)])
    group.results <- cbind(group.results,c(p.survive,ind.results$g[rep*i],max.evolve,g.1))
  }
  
  
  group.results <- as.data.frame(t(group.results[,-1]))
  colnames(group.results) <- c("p.survive","g","max.evolve","g.1")
  
  write.csv(x = as.data.frame(group.results), file = "Tier1SISRcorrected.csv")
  
  
  par(mfrow=c(1,1))
  
  ggplot(group.results,aes(g.1,p.survive)) + geom_line() + theme_minimal() +
    xlab("Initial trait value") + ylab("Probability of survival")
  