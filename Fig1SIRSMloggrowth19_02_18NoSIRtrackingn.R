  ######################Plotting one run
  
  library(deSolve)
  #Parameters
  
  r <- 0.0866       # Intrinsic rate of increase of susceptible
  mu <- 0.0002      # Natural mortality rate
  B <- 5       # Transmission coefficient
  alpha <- 0.5       # Mortality rate
  K <- 200          # Carrying Capacity
  phi <- 0.011      # Resistance loss rate
  
  Br <- 0.07       #plague reservoir infection rate
  lambda <- 0.006  #reservoir decay rate
  
  #starting conditions
  
  S.start <- 199
  I.start <- 1
  R.start <- 0
  M.start <- 0
  g.start <- 0.005  # Initial resistance trait
  sigma <- 0.0005   # Trait variance
  
  t <- 10000
  
  SIR <- function(t,state,parameters){
    with(as.list(c(state, parameters)), {
      dS <- r*(S+I+R)*(1-(S+I+R)/K)-mu*S-B*S*I-Br*S*M+phi*R
      dI <- B*S*I+Br*S*M-alpha*I-mu*I-g*I
      dR <- g*I-mu*R-phi*R
      dM <- alpha*I-lambda*M
      dg <- sigma^2*(I+M/B)/(R+1)
      list(c(dS,dI,dR,dM,dg))
    })
  }
  parameters <- c(sigma=sigma,r=r,mu=mu,B=B,alpha=alpha,K=K,lambda=lambda,Br=Br)
  state <- c(S = S.start,I = I.start,R = R.start,M = M.start,g = g.start)
  times <- seq(0, t, by = 0.1)
  out <- ode(y = state, times = times, func = SIR, parms = parameters)
  plot(out)
  
  output <- as.data.frame(out)
  
  N <- output$S+output$I+output$R
  R0 <- B/(alpha+mu+output$g)
  FoI <- B*output$S*output$I+Br*output$S*output$M
  growth <- r*(output$S+output$I+output$R)*(1-(output$S+output$I+output$R)/K)+phi*output$R
  
  par(mfrow=c(2,2))
  plot(times,N,type="l",xlab="time",ylab="N",main="Population Size")
  plot(times,R0,type="l",xlab="time",ylab="R0",main="R0")
  plot(times,FoI,type="l",xlab="time",ylab="FoI",main="Force of Infection",ylim=c(0,10))
  plot(times,growth,type="l",xlab="time",ylab="growth rate",main="Rate of increase of susceptibles")
  
  
  
  ###############Plotting many deterministic runs
  
  g.values <- seq(0,0.2,0.005)
  sigma.values <- c(0,0.0005,seq(0.005,0.2,0.005))
  
  t <- 10000
  
  
  popsize <- matrix(nrow=length(g.values),ncol=length(sigma.values),NA)
  det.trait <- matrix(nrow=length(g.values),ncol=length(sigma.values),NA)
  scaled.trait <- matrix(nrow=length(g.values),ncol=length(sigma.values),NA)
  # 
  # for(i in 1:length(g.values)){
  #   g.start <- g.values[i]
  #   for(j in 1:length(sigma.values)){
  #     sigma <- sigma.values[j]
  #     parameters <- c(sigma=sigma,r=r,mu=mu,B=B,alpha=alpha,K=K,phi=phi,lambda=lambda,Br=Br)
  #     state <- c(S = S.start,I = I.start,R = R.start,M = M.start,g = g.start)
  #     times <- seq(0, t, by = 0.1)
  #     out <- ode(y = state, times = times, func = SIR, parms = parameters)
  #     output <- as.data.frame(out)
  #     popsize[i,j] <- output$S[length(output$S)]+output$I[length(output$I)]+output$R[length(output$R)]
  #     det.trait[i,j] <- output$g[length(output$g)]
  #     scaled.trait[i,j] <- (output$g[length(output$g)]-output$g[1])/sigma
  #   }
  # }

  
  
  #################################Plotting many stochastic
  library(tidyverse)
  
  rho <- 0.08
  
  parameters <- c(r=r,mu=mu,B=B,alpha=alpha,K=K,rho=rho,lambda=lambda,Br=Br)
  state <- c(S.start = S.start,I.start = I.start,R.start = R.start,M.start=M.start)
  t <- 10000
  
  SIR.sto <- function(g.start,sigma,t,r,mu,B,alpha,K,rho,lambda,Br,S.start,I.start,R.start,M.start){
    start_time <- Sys.time()
    XT <- 10000000 #potential number of events
    WT <- rep(NA,XT) #start time
    WT[1] <- 0
    
    N <- rep(NA,XT)
    N[1] <- S.start+I.start+R.start
    S <- S.start
    I <- I.start
    R <- R.start
    M <- M.start
    
    g <- rnorm(N[1],mean=g.start,sd=sigma)
    g[g<0] <- 0
    g.avg <- rep(NA,XT)
    g.avg[1] <- mean(g)
    
    SJ <- rep(0,XT) #Make list for time steps
    
    individuals <- rep(0,N[1]); #create vector for all individual animals 
    #and their potential classes
    if(S.start>0){individuals[1:S.start] <- 1}; #Assign S individuals into S
    if(I.start>0){individuals[(S.start+1):(S.start+I.start)] <- 2};
    #Assign I individuals
    if(R.start>0){individuals[(S.start+I.start+1):(S.start+I.start+R.start)] <- 3};
    #Assign R inidividuals
    
    for(i in 2:XT){
      
      Itrait <- which(individuals==2)
      gtrait <- sum(g[Itrait])
      
      #Potential Events
      birth <- r*N[i-1]*(1-N[i-1]/K)
      nat <- mu*N[i-1]
      infect <- B*S*I+Br*S*M
      disease <- alpha*I
      recover <- gtrait
      suscept <- phi*R
      reservoir <- lambda*M
      
      sumrates <- birth+nat+infect+disease+recover+suscept+reservoir
      if(sumrates <= 0) break
      
      #Determine time to next event (Sojourn time)
      SJ[i] <- rexp(1,rate=sumrates)
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
        if(S<=1){  
          Diseased <- which(individuals==1,arr.ind=TRUE)
        } else{
          Diseased <- sample(which(individuals==1,arr.ind=TRUE),1)
        } #end else
        individuals[Diseased] <- 2 #Infected!
      } else{if(u<(birth+nat+infect+disease)/sumrates){
        if(I<=1){  
          Released <- which(individuals==2,arr.ind=TRUE)
        } else{
          Released <- sample(which(individuals==2,arr.ind=TRUE),1)
        } #end else
        individuals <- individuals[-Released] #Kill it
        g <- g[-Released]
        M <- M+1
      } else{if(u<(birth+nat+infect+disease+recover)/sumrates){
        if(I<=1){  
          EyeOfTiger <- which(individuals==2,arr.ind=TRUE)
        } else{
          EyeOfTiger=sample(which(individuals==2,arr.ind=TRUE),1,prob=g[Itrait]/gtrait)
        } #end else
        individuals[EyeOfTiger] <- 3
      } else{if(u<(birth+nat+infect+disease+recover+suscept)/sumrates){
        #Choose a susceptible to infect
        if(R<=1){  
          Susceptible <- which(individuals==3,arr.ind=TRUE)
        } else{
          Susceptible <- sample(which(individuals==3,arr.ind=TRUE),1)
        } #end else
        individuals[Susceptible] <- 1 #Recovered to susceptible
      } else{
        M <- M-1
      } #end else 6
      } #end else 5
      } #end else 4
      } #end else 3
      } #end else 2
      } #end else 1
      
      N[i] <- length(individuals)
      g.avg[i] <- sum(g)/N[i]
      
      S <- length(which(individuals==1))
      I <- length(which(individuals==2))
      R <- length(which(individuals==3))
      
      if(N[i] == 0) break
      if(WT[i] >= t) break
      
    }#end forloop
    
    N <- N[which(N!="NA")]
    evolve <- (max(g.avg,na.rm=TRUE)-g.avg[1])/sigma
    survive <- 0
    if(N[length(N)]>0) survive <- 1
    
    end_time <- Sys.time()
    
    output <- c(N[length(N)],min(N),g.start,g.avg[1],max(g.avg,na.rm=TRUE),evolve,survive,sigma)
    return(output)
    
  } #end function
  
  
  
  
  
  library(doParallel)
  detectCores()
  c <- makeCluster(3)
  registerDoParallel(c)
  getDoParWorkers()
  
  
  rep <- 100
  
  ###########
  #ORIGINAL
  ###########
  
  # results <- foreach(i = 1:length(g.values), .combine= cbind) %dopar% {
  #   print(i)
  #   ind.results <- NA
  #   g.start <- g.values[i]
  #   for(j in 1:length(sigma.values)){
  #     sigma <- sigma.values[j]
  #     d.trait <- det.trait[i,j]
  #     for(k in 1:rep){
  #       x <- append(SIR.sto(state,t,parameters,g.start,sigma),sigma)
  #       x <- append(x,d.trait)
  #       ind.results <- cbind(ind.results,x)
  #     }
  #   }
  #   ind.results <- ind.results[,-1]
  #   return(ind.results)
  # }
  
  ##########
  #NEW
  ##########
  
  results <- foreach(i = 1:length(g.values), .combine= cbind) %dopar% {
    print(i)
    g.start <- rep(g.values[i],(length(sigma.values)*rep))
    sigma <- rep(sigma.values,each=rep)
    ind.results <- mapply(SIR.sto,g.start,sigma,t=t,r=r,mu=mu,B=B,alpha=alpha,K=K,rho=rho,lambda=lambda,Br=Br,S.start = S.start,I.start = I.start,R.start = R.start,M.start=M.start)
    return(ind.results)
  }
  
  ind.results <- as.data.frame(t(results))
  # ind.results <- cbind(ind.results,NA)
  colnames(ind.results) <- c("Nt","min.N","g","g.1","max.g","evolve","survive","sigma")
  # colnames(ind.results) <- c("Nt","min.N","g","g.1","max.g","evolve","survive","sigma","determ.g","sto.det.dif")
  # d.g <- ind.results$max.g-ind.results$determ.g
  # ind.results$sto.det.dif[which(d.g >= sigma)] <- as.character("greater")
  # ind.results$sto.det.dif[which(d.g < sigma)] <- as.character("equal")
  # ind.results$sto.det.dif[which(d.g <= -sigma)] <- as.character("lesser")
  ind.results$survive <- as.numeric(ind.results$survive)
  
  write.csv(x = as.data.frame(ind.results), file = "Fig1ResultsSIRSMloggrowth.csv")
  
  
  group.results <- NA
  for(i in 1:(length(ind.results[,1])/rep)){
    p.survive <- mean(ind.results$survive[(1+rep*(i-1)):(rep*i)])
    max.evolve <- mean(ind.results$evolve[(1+rep*(i-1)):(rep*i)])
    g.1 <- mean(ind.results$g.1[(1+rep*(i-1)):(rep*i)])
    variance <- ind.results$sigma[rep*i]
    group.results <- cbind(group.results,c(p.survive,ind.results$g[rep*i],max.evolve,variance,g.1))
  }
  
  group.results <- as.data.frame(t(group.results[,-1]))
  colnames(group.results) <- c("p.survive","g","max.evolve","variance","g.1")
  
  par(mfrow=c(1,1))
  
  ggplot(group.results,aes(g.1,p.survive,group=variance)) + geom_line(aes(color=variance)) + 
    geom_line(data=group.results[which(group.results$variance==0.0005),],aes(g.1,p.survive),color="red") +
    xlab("Initial trait value") + ylab("Probability of survival")
  
  ggplot(group.results,aes(g.1,max.evolve,group=variance)) + geom_line(aes(color=variance))+ 
    geom_line(data=group.results[which(group.results$variance==0.0005),],aes(g.1,max.evolve),color="red") +
    xlab("Initial trait value") + ylab("Evolution standardized by variance")
  