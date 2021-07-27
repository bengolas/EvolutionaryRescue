###############################################################################
#The following code will define the stochastic plague model system function 
#and generate data from simulations. It runs a series of *rep* iterations 
#across ranges of average initial trait size (mu_gamma) 
#and variance (sd == standard deviation of genetic variance, 
#i.e. sd^2=variance).
#Throughout the script, 'g' is shorthand for 'gamma'.
#
#Following code for simulation output, there is another section for performing
#PRCC parameter sensitivity analysis.
###############################################################################


rm(list=ls())


#Ranges to evaluate
mugamma.values <- seq(0,0.05,0.0025) #values of mu_gamma to investigate
sd.values <- c(0,0.0005,seq(0.005,0.2,0.005)) #values of variance to investigate
rho.values <- 0.08 #multiple values of mutation/recombination rate can be explored

r = .0866       # Intrinsic rate of increase of susceptible (host)
K =  200        # Carrying capacity (host)
mu = .0002      # Natural mortality rate (host)
Bc = .073       # Airborne transmission rate
Br = .073       # Reservoir transmission rate
B = 20          # Number burrows host enters
delta = .21    # Rate of becoming infectious following exposure (host)
alpha = .5      # mortality rate (host)
lambda = .006   # Reservoir decay rate
phi = .011      # Loss of resistance
p = 0           # Proportion born resistant
rho = 0.08      # Probability of mutation/recombination in resistance rate


#prairie dog values for testing purposes
# sd <- 0.0005
# mugamma <- 0.005


PLAGUE <- function(r,K,mu,Bc,Br,B,delta,alpha,lambda,phi,p,rho,sd,mugamma){
    
    ######################################
    #Initial settings
    ######################################
    
    XT = 1000000      # Number of events
    
    WT<- rep(0,XT) #Create a time vector with space for each potential incident
    N<- rep(0,XT) #Vector to keep track of population size over time
    
    S <- 199 #Starting values for population size in each class
    E <- 1
    I <- 0
    R <- 0
    M <- 0
    
    N[1] = S + E + I + R
    
    WT.initial.size <- 0          #Initial wait time
    
    #################################
    #Individuals and their traits   #
    #################################
    
    individuals <- rep(0,N[1]) #create vector for all individual animals 
    #and their potential classes
    if(S>0){individuals[1:S]=1} #Assign S individuals into S
    if(E>0){individuals[(S+1):(S+E)]=2} 
    #Assign E individuals into E
    if(I>0){individuals[(S+E+1):(S+E+I)]=3}
    #Assign I individuals
    if(R>0){individuals[(S+E+I+1):(S+E+I+R)]=4}
    #Assign R inidividuals
    
    g <- rnorm(N[1], mean=mugamma, sd=sd) #Generate gamma values for initial population
    g[g<0]=0 #Truncate at 0 to prevent negative values
    g.avg <- rep(0,XT) #Create a vector to store average population resistance rate at each time step
    g.avg[1] = mean(g)
    
    SJ <- WT #Keep track of time in each step
    
    
    for(i in 2:XT){
      #################################
      #sum(gamma) for prairie dogs in E
      #################################
      
      Etrait=which(individuals==2)
      gtrait=sum(g[Etrait])
      
      #############################
      #Potential Events
      #############################
      
      SbornS=r*S*(1-N[i-1]/K)                               #Susceptibles born by S
      SbornR=r*R*(1-p)*(1-N[i-1]/K)                         #Susceptibles born by R
      Sdie=mu*S                                             #Natural deaths
      Edie=mu*E
      Rdie=mu*R
      Infect=Bc*S*I/N[i-1]+Br*S*M/B                         #Infecting susceptibles
      ExpInf=delta*E                                        #Exposed to Infectious
      Mort=alpha*I                                          #Disease-induced mortality
      ResGain=gtrait                                        #Exposed gain resistance
      ResLoss=phi*R                                         #Resistant become susceptible
      ResBorn=r*R*p*(1-N[i-1]/K)                            #Individuals born resistant
      Decay=lambda*M                                        #Reservoir decay
      
      sumrates=SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain+ResLoss+ResBorn+Decay
      
      if (sumrates <= 0) break #Stop-check in case typos lead to negative probabilities
      
      #Determine time to next event (Sojourn time)
      SJ[i]=rexp(1,rate=sumrates)
      WT[i]=WT[i-1]+SJ[i]
      
      
      
      #Decide event to occur and update populations
      u<-runif(1,min=0,max=1)
      
      if(u<(SbornS)/sumrates){ #Susceptible born by susceptible parent
        #Choose a susceptible parent
        if(S<=1){  #sample cannot handle a list of 1 to choose from
          MummySS=which(individuals==1,arr.ind=TRUE) #pick only one if 1 available
        } else{
          MummySS=sample(which(individuals==1,arr.ind=TRUE),1) #random sample
        } #end else
        #Give parent a baby
        individuals=append(individuals,1) #New susceptible individual 
        g=append(g,0)             #New trait for individuals
        #Determine if mutation occurs and assign new individual a trait
        mut=runif(1,min=0,max=1) #random chance generator
        if(mut<rho){ #Mutation occurs
          #Assign mutated parent trait
          g[length(g)]=g[MummySS]+rnorm(1, mean=0, sd=sd) 
        } else{g[length(g)]=g[MummySS]} #Assign parent trait to newest individual
        g[g<0]=0
      } else{if(u<(SbornS+SbornR)/sumrates){ #Susceptible born by resistant parent
        #Choose a resistant parent
        if(R<=1){  
          MummyRS=which(individuals==4,arr.ind=TRUE)
        } else{
          MummyRS=sample(which(individuals==4,arr.ind=TRUE),1)
        } #end else
        #Give parent a baby
        individuals=append(individuals,1) #New susceptible individual 
        g=append(g,0)                    #New trait for individual
        #Determine if mutation occurs and assign new individual a trait
        mut=runif(1,min=0,max=1) #random chance generator
        if(mut<rho){
          #g[length(g)]=g[MummyRS]+rnorm(1, mean=0, sd=(0.05*g[MummyRS])) #Mutation occurs 
          g[length(g)]=g[MummyRS]+rnorm(1, mean=0, sd=sd) #Use fixed standard deviation
        } else{g[length(g)]=g[MummyRS]}     #No mutation
        g[g<0]=0
      } #end SbornR
        
        else{if(u<(SbornS+SbornR+Sdie)/sumrates){ #Susceptible dies naturally
          #Choose an individual to off
          if(S<=1){  
            GonerS=which(individuals==1,arr.ind=TRUE)
          } else{
            GonerS=sample(which(individuals==1,arr.ind=TRUE),1)
          } #end else
          individuals = individuals[-GonerS] #Kill it
          g = g[-GonerS] #Remove from trait pool
        } #end Sdie
          
          else{if(u<(SbornS+SbornR+Sdie+Edie)/sumrates){ #Exposed dies naturally
            #Choose an individual to off
            if(E<=1){  
              GonerE=which(individuals==2,arr.ind=TRUE)
            } else{
              GonerE=sample(which(individuals==2,arr.ind=TRUE),1)
            } #end else
            individuals = individuals[-GonerE] #Kill it
            g = g[-GonerE] #Remove from trait pool
          } #end Edie
            
            else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie)/sumrates){ #Resistant dies naturally
              #Choose an individual to off
              if(R<=1){  
                GonerR=which(individuals==4,arr.ind=TRUE)
              } else{
                GonerR=sample(which(individuals==4,arr.ind=TRUE),1)
              } #end else
              individuals = individuals[-GonerR] #Kill it
              g = g[-GonerR] #Remove from trait pool
            } #end Rdie
              
              else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect)/sumrates){ #Infect a susceptible
                #Choose a susceptible to infect
                if(S<=1){  
                  Diseased=which(individuals==1,arr.ind=TRUE)
                } else{
                  Diseased=sample(which(individuals==1,arr.ind=TRUE),1)
                } #end else
                individuals[Diseased]=2 #Infected!
              } #end Infect
                
                else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf)/sumrates){ #Exposed to infectious
                  if(E<=1){  
                    Ohshit=which(individuals==2,arr.ind=TRUE)
                  } else{
                    Ohshit=sample(which(individuals==2,arr.ind=TRUE),1)
                  } #end else
                  individuals[Ohshit]=3 #Doomed to be released from this mortal coil
                } #end ExpInf
                  
                  else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort)/sumrates){ #Infectious dies
                    if(I<=1){  
                      Released=which(individuals==3,arr.ind=TRUE)
                    } else{
                      Released=sample(which(individuals==3,arr.ind=TRUE),1)
                    } #end else
                    individuals=individuals[-Released] #Kill
                    M = M+1 #Add to reservoir
                    g = g[-Released] #Remove from trait pool
                  } #end Mort
                    
                    else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain)/sumrates){ #Exposed gains resistance
                      if(E<=1){  
                        EyeOfTiger=which(individuals==2,arr.ind=TRUE)
                      } else{
                        EyeOfTiger=sample(Etrait,1,prob=g[Etrait]/gtrait)
                      } #end else
                      
                      individuals[EyeOfTiger]=4
                    } #end ResGain
                      
                      else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain+ResLoss)/sumrates){ #Resistant becomes susceptible
                        if(R<=1){  
                          WeakestLink=which(individuals==4,arr.ind=TRUE)
                        } else{
                          WeakestLink=sample(which(individuals==4,arr.ind=TRUE),1)
                        } #end else
                        individuals[WeakestLink]=1
                      } #end ResLoss
                        
                        else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain+ResLoss+ResBorn)/sumrates){ #New resistant born from resistant parents
                          #Choose a resistant parent
                          if(R<=1){  
                            MummyRR=which(individuals==4,arr.ind=TRUE)
                          } else{
                            MummyRR=sample(which(individuals==4,arr.ind=TRUE),1)
                          } #end else
                          #Give parent a baby
                          individuals=append(individuals,4) #New susceptible individual column
                          g=append(g,0)                    #New trait column
                          #Determine if mutation occurs and assign new individual a trait
                          mut=runif(1,min=0,max=1) #random chance generator
                          if(mut<rho){
                            #g[length(g)]=g[MummyRR]+rnorm(1, mean=0, sd=(0.05*g[MummyRR])) #Mutation occurs
                            g[length(g)]=g[MummyRR]+rnorm(1, mean=0, sd=sd) #Use fixed standard deviation
                          } else{g[length(g)]=g[MummyRR]}     #No mutation
                          g[g<0]=0
                        } #end ResBorn
                          
                          else{if(u<=(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain+ResLoss+ResBorn+Decay)/sumrates){ #Reservoir decays
                            M = M-1
                          } #end Decay
                            
                          } #end else (Decay) 
                        } #end else (ResBorn)
                      } #end else (ResLoss)
                    } #end else (ResGain)
                  } #end else (Mort)
                } #end else (ExpInf)
              } #end else (Infect)
            } #end else (Rdie)
          } #end else (Edie)
        } #end else (Sdie)
      } #end first else (SbornR)
      
      
      
      
      
      #LAST STEP of for loop, save population counts
      S=length(which(individuals==1))
      E=length(which(individuals==2))
      I=length(which(individuals==3))
      R=length(which(individuals==4))
      
      N[i]=S+E+I+R
      
      g.avg[i]=mean(g)
      
      if(N[i]<=0) break #end
      if(S<0) break #end
      if(E<0) break #end
      if(I<0) break #end
      if(R<0) break #end
      if(M<0) break #end
      if(WT[i]>10000) break #end
    }#end for loop
    
    #########################################################
    #COLLECTING RESULTS
    #########################################################
    
    #Clipping vectors to remove NA
    tm=length(WT[WT>0])+1
    N=N[1:tm]
    g.avg <- g.avg[1:tm]
    if(N[tm]==0){g.avg[tm] <- g.avg[tm-1]}
    
    
    ############
    #Outputs
    ############
    survival <- 1
    if(N[tm]==0){survival <- 0}
    g.max <- max(g.avg)
    
    
    
    OUTPUT <- c(N[length(N)],min(N),mugamma,g.avg[1],g.max,survival,sd)
    return(OUTPUT)
    
  }#End function
  

##################################################
#Generate data
##################################################

reps <- 200 #number of iterations per parameterization

#Run parallel processing
library(doParallel)
detectCores()
c <- makeCluster(3)
registerDoParallel(c)
getDoParWorkers()

library(foreach)
      
load.results <- foreach(k = 1:length(mugamma.values), .combine= cbind) %dopar% {
  mugamma <- mugamma.values[k]
  MATR <- NULL
  for(kk in 1:length(sd.values)){
    sd <- sd.values[kk]
    for(kkkk in 1:length(rho.values)){
      rho <- rho.values[kkkk]
      M <- NULL
      for(j in 1:reps){
        x <- PLAGUE(r,K,mu,Bc,Br,B,delta,alpha,lambda,phi,p,rho,sd,mugamma)
        x[which(x==Inf)] <- NA
        M <- cbind(M,x)
      }
      
      y <- rowMeans(M,na.rm=TRUE)
      MATR <- cbind(MATR,y)
    }
  }
  return(MATR)
}

#Format compiled data
P <- as.data.frame(t(load.results))
#Outputs are as follows: 
#(Nt= end population size), 
#(min.N= minimum population size), 
#(trait= mu_gamma), 
#(g.1= bar(gamma)_0), 
#(g.max= maximum average trait value over time, 
#(survivor.evolve= evolution of only surviving populations), 
#(survival= 0/1 extinct/survived), 
#(sd= standard deviation of genetic variance)
colnames(P) <- c("Nt","min.N","trait","g.1","g.max","survival","sd")



write.csv(x = as.data.frame(P), file = "PlagueData.csv")







##############################################################################
#The following code is used to generate results for the parameter
#sensitivity analysis. The code begins by defining the functions for the
#plague (PLAGUE) and system (SIR.sto) models to run simulations. Then, we
#generate LHS matrices, run simulations for parameterizations of the matrices
#with and without evolution, and compile results into survival probabilities
#and potential for ER. Finally, we run partial rank correlation coefficient 
#parameter sensitivity analysis using these results.
##############################################################################


SIR.sto <- function(mugamma,sigma,r,mu,B,alpha,K,rho,lambda,Br,phi){
  XT <- 10000000 #potential number of events
  WT <- rep(NA,XT) #start time
  WT[1] <- 0
  
  t <- 10000
  
  S.start <- K-1
  I.start <- 1
  R.start <- 0
  M.start <- 0
  
  N <- rep(NA,XT)
  N[1] <- S.start+I.start+R.start
  S <- S.start
  I <- I.start
  R <- R.start
  M <- M.start
  
  g <- rnorm(N[1],mean=mugamma,sd=sigma)
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
  
  output <- c(N[length(N)],min(N),evolve,survive,g.avg[1])
  return(output)
  
} #end function







PLAGUE <- function(mugamma,sigma,r,mu,Bc,alpha,K,rho,lambda,Br,phi,B,delta){
  
  ######################################
  #Initial settings
  ######################################
  
  #Standardize variables across models
  p <- 0
  eps <- sigma
  
  XT = 1000000;      # Number of events
  
  WT<- rep(0,XT); #Create a time vector with space for each potential incident
  S<- 0  #Create population vectors to store information following each incident
  E<- 0
  I<- 0
  R<- 0
  M<- 0
  N<- rep(0,XT);
  
  S = K-1;
  E = 1;
  I = 0;
  R = 0;
  M = 0;
  
  N[1] = S + E + I + R;
  
  WT.initial.size = 0;          #Initial wait time
  
  #################################
  #Individuals and their traits   #
  #################################
  
  individuals=rep(0,N[1]); #create vector for all individual animals 
  #and their potential classes
  if(S>0){individuals[1:S]=1}; #Assign S individuals into S
  if(E>0){individuals[(S+1):(S+E)]=2}; 
  #Assign E individuals into E
  if(I>0){individuals[(S+E+1):(S+E+I)]=3};
  #Assign I individuals
  if(R>0){individuals[(S+E+I+1):(S+E+I+R)]=4};
  #Assign R inidividuals
  
  g <- rnorm(N[1], mean=mugamma, sd=eps) #Initial g values (additive)
  #g <- rnorm(N[1], mean=g.initial.size, sd=(0.05*g.initial.size)) #Initial g values (multiplicative)
  #g <- rep(g.initial.size,N[1]) #Initial g values (all individuals same starting value)
  g[g<0]=0
  g.avg <- rep(0,XT)
  g.avg[1] = mean(g)
  
  SJ <- WT #Make list for time steps
  
  
  for(i in 2:XT){
    #############################
    #sum(g) for prairie dogs in E
    #############################
    
    Etrait=which(individuals==2)
    gtrait=sum(g[Etrait])
    
    #############################
    #Potential Events
    #############################
    
    SbornS=r*S*(1-N[i-1]/K)                               #Susceptibles born by S
    SbornR=r*R*(1-p)*(1-N[i-1]/K)                         #Susceptibles born by R
    Sdie=mu*S                                             #Natural deaths
    Edie=mu*E
    Rdie=mu*R
    Infect=Bc*S*I/N[i-1]+Br*S*M/B          #Infecting susceptibles
    ExpInf=delta*E                                       #Exposed to Infectious
    Mort=alpha*I                                          #Disease-induced mortality
    ResGain=gtrait                                             #Exposed gain resistance
    ResLoss=phi*R                                         #Resistant become susceptible
    ResBorn=r*R*p*(1-N[i-1]/K)                            #Individuals born resistant
    Decay=lambda*M                                        #Reservoir decay
    
    sumrates=SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain+ResLoss+ResBorn+Decay
    
    if (sumrates <= 0) break
    
    #Determine time to next event (Sojourn time)
    SJ[i]=rexp(1,rate=sumrates)
    WT[i]=WT[i-1]+SJ[i]
    
    
    
    #Decide event to occur and update populations
    u<-runif(1,min=0,max=1)
    
    if(u<(SbornS)/sumrates){ #Susceptible born by susceptible parent
      #Choose a susceptible parent
      if(S<=1){  #sample cannot handle a list of 1 to choose from
        MummySS=which(individuals==1,arr.ind=TRUE) #pick only one if 1 available
      } else{
        MummySS=sample(which(individuals==1,arr.ind=TRUE),1) #random sample
      } #end else
      #Give parent a baby
      individuals=append(individuals,1); #New susceptible individual 
      g=append(g,0);             #New trait for individuals
      #Determine if mutation occurs and assign new individual a trait
      mut=runif(1,min=0,max=1); #random chance generator
      if(mut<rho){ #Mutation occurs
        #g[length(g)]=g[MummySS]+rnorm(1, mean=0, sd=(0.05*g[MummySS])); #Assign mutated parent trait
        g[length(g)]=g[MummySS]+rnorm(1, mean=0, sd=eps) #Use fixed standard deviation
      } else{g[length(g)]=g[MummySS]}; #Assign parent trait to newest individual
      g[g<0]=0
    } else{if(u<(SbornS+SbornR)/sumrates){ #Susceptible born by resistant parent
      #Choose a resistant parent
      if(R<=1){  
        MummyRS=which(individuals==4,arr.ind=TRUE)
      } else{
        MummyRS=sample(which(individuals==4,arr.ind=TRUE),1)
      } #end else
      #Give parent a baby
      individuals=append(individuals,1) #New susceptible individual 
      g=append(g,0)                    #New trait for individual
      #Determine if mutation occurs and assign new individual a trait
      mut=runif(1,min=0,max=1) #random chance generator
      if(mut<rho){
        #g[length(g)]=g[MummyRS]+rnorm(1, mean=0, sd=(0.05*g[MummyRS])) #Mutation occurs 
        g[length(g)]=g[MummyRS]+rnorm(1, mean=0, sd=eps) #Use fixed standard deviation
      } else{g[length(g)]=g[MummyRS]}     #No mutation
      g[g<0]=0
    } #end SbornR
      
      else{if(u<(SbornS+SbornR+Sdie)/sumrates){ #Susceptible dies naturally
        #Choose an individual to off
        if(S<=1){  
          GonerS=which(individuals==1,arr.ind=TRUE)
        } else{
          GonerS=sample(which(individuals==1,arr.ind=TRUE),1)
        } #end else
        individuals = individuals[-GonerS] #Kill it
        g = g[-GonerS] #Remove from trait pool
      } #end Sdie
        
        else{if(u<(SbornS+SbornR+Sdie+Edie)/sumrates){ #Exposed dies naturally
          #Choose an individual to off
          if(E<=1){  
            GonerE=which(individuals==2,arr.ind=TRUE)
          } else{
            GonerE=sample(which(individuals==2,arr.ind=TRUE),1)
          } #end else
          individuals = individuals[-GonerE] #Kill it
          g = g[-GonerE] #Remove from trait pool
        } #end Edie
          
          else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie)/sumrates){ #Resistant dies naturally
            #Choose an individual to off
            if(R<=1){  
              GonerR=which(individuals==4,arr.ind=TRUE)
            } else{
              GonerR=sample(which(individuals==4,arr.ind=TRUE),1)
            } #end else
            individuals = individuals[-GonerR] #Kill it
            g = g[-GonerR] #Remove from trait pool
          } #end Rdie
            
            else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect)/sumrates){ #Infect a susceptible
              #Choose a susceptible to infect
              if(S<=1){  
                Diseased=which(individuals==1,arr.ind=TRUE)
              } else{
                Diseased=sample(which(individuals==1,arr.ind=TRUE),1)
              } #end else
              individuals[Diseased]=2 #Infected!
            } #end Infect
              
              else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf)/sumrates){ #Exposed to infectious
                if(E<=1){  
                  Ohshit=which(individuals==2,arr.ind=TRUE)
                } else{
                  Ohshit=sample(which(individuals==2,arr.ind=TRUE),1)
                } #end else
                individuals[Ohshit]=3 #Doomed to be released from this mortal coil
              } #end ExpInf
                
                else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort)/sumrates){ #Infectious dies
                  if(I<=1){  
                    Released=which(individuals==3,arr.ind=TRUE)
                  } else{
                    Released=sample(which(individuals==3,arr.ind=TRUE),1)
                  } #end else
                  individuals=individuals[-Released] #Kill
                  M = M+1 #Add to reservoir
                  g = g[-Released] #Remove from trait pool
                } #end Mort
                  
                  else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain)/sumrates){ #Exposed gains resistance
                    if(E<=1){  
                      EyeOfTiger=which(individuals==2,arr.ind=TRUE)
                    } else{
                      EyeOfTiger=sample(Etrait,1,prob=g[Etrait]/gtrait)
                    } #end else
                    
                    individuals[EyeOfTiger]=4
                  } #end ResGain
                    
                    else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain+ResLoss)/sumrates){ #Resistant becomes susceptible
                      if(R<=1){  
                        WeakestLink=which(individuals==4,arr.ind=TRUE)
                      } else{
                        WeakestLink=sample(which(individuals==4,arr.ind=TRUE),1)
                      } #end else
                      individuals[WeakestLink]=1
                    } #end ResLoss
                      
                      else{if(u<(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain+ResLoss+ResBorn)/sumrates){ #New resistant born from resistant parents
                        #Choose a resistant parent
                        if(R<=1){  
                          MummyRR=which(individuals==4,arr.ind=TRUE)
                        } else{
                          MummyRR=sample(which(individuals==4,arr.ind=TRUE),1)
                        } #end else
                        #Give parent a baby
                        individuals=append(individuals,4) #New susceptible individual column
                        g=append(g,0)                    #New trait column
                        #Determine if mutation occurs and assign new individual a trait
                        mut=runif(1,min=0,max=1) #random chance generator
                        if(mut<rho){
                          #g[length(g)]=g[MummyRR]+rnorm(1, mean=0, sd=(0.05*g[MummyRR])) #Mutation occurs
                          g[length(g)]=g[MummyRR]+rnorm(1, mean=0, sd=eps) #Use fixed standard deviation
                        } else{g[length(g)]=g[MummyRR]}     #No mutation
                        g[g<0]=0
                      } #end ResBorn
                        
                        else{if(u<=(SbornS+SbornR+Sdie+Edie+Rdie+Infect+ExpInf+Mort+ResGain+ResLoss+ResBorn+Decay)/sumrates){ #Reservoir decays
                          M = M-1
                        } #end Decay
                          
                        } #end else (Decay) 
                      } #end else (ResBorn)
                    } #end else (ResLoss)
                  } #end else (ResGain)
                } #end else (Mort)
              } #end else (ExpInf)
            } #end else (Infect)
          } #end else (Rdie)
        } #end else (Edie)
      } #end else (Sdie)
    } #end first else (SbornR)
    
    
    
    
    
    #LAST STEP of for loop, save population counts
    S=length(which(individuals==1))
    E=length(which(individuals==2))
    I=length(which(individuals==3))
    R=length(which(individuals==4))
    
    N[i]=S+E+I+R
    
    g.avg[i]=mean(g)
    
    if(N[i]<=0) break; #end
    if(S<0) break; #end
    if(E<0) break; #end
    if(I<0) break; #end
    if(R<0) break; #end
    if(M<0) break; #end
    if(WT[i]>10000) break; #end
  }#end for loop
  
  #########################################################
  #COLLECTING RESULTS
  #########################################################
  
  #Clipping vectors hack
  tm=length(WT[WT>0])+1
  N=N[1:tm]
  g.avg <- g.avg[1:tm]
  if(N[tm]==0){g.avg[tm] <- g.avg[tm-1]}
  
  
  ############
  #Outputs
  ############
  survive <- 1
  if(N[tm]==0){survive <- 0}
  evolve <- (max(g.avg,na.rm=TRUE)-g.avg[1])/sigma
  survive <- 0
  if(N[length(N)]>0) survive <- 1
  
  
  
  output <- c(N[length(N)],min(N),evolve,survive,g.avg[1])
  return(output)
  
}



##########################
#Generate LHS matrices
##########################

#Get survival and g.1 from first run, then run again with same 
#parameters but sigma=0 and mugamma=mean(g.1)

library(fitur)
library(sensitivity)
library(lhs)

set.seed(28)

n <- 400 #number of parameter combinations
sims <- 200 #number of simulations per parameter combination


#Generate LHS matrices for the reduced model (LHS.SIR) and the plague system
#model (LHS.PL) and run the given functions for each of 'n' parameter sets in the
#LHS matrix for 'sims' number of times.

LHS <- randomLHS(n,11)

LHS.SIR <- data.frame(mugamma = qunif(LHS[,1],min=0,max=0.2),
                      sigma = qunif(LHS[,2],min=0,max=0.2),
                      r = qunif(LHS[,3],min=0,max=0.1),
                      mu = qunif(LHS[,4],min=0,max=0.001),
                      B = qunif(LHS[,5],min=0.05,max=0.55),
                      alpha = qunif(LHS[,6],min=0.25,max=1),
                      K = qunif(LHS[,7],min=1,max=500),
                      rho = qunif(LHS[,8],min=0,max=1),
                      lambda = qunif(LHS[,9],min=0,max=0.02),
                      Br = qunif(LHS[,10],min=0.02,max=0.12),
                      phi = qunif(LHS[,11],min=0,max=0.02))


LHS.SIR.o <- LHS.SIR

#Repeat LHS matrix 'sims' times so we can use the mapply function instead of
#for loops to generate our data
LHS.SIR <- LHS.SIR[rep(1:nrow(LHS.SIR),each=sims),] 

SAmat.ind.SIR <- mapply(SIR.sto,mugamma=LHS.SIR$mugamma,sigma=LHS.SIR$sigma,r=LHS.SIR$r,mu=LHS.SIR$mu,
                        B=LHS.SIR$B,alpha=LHS.SIR$alpha,K=LHS.SIR$K,rho=LHS.SIR$rho,lambda=LHS.SIR$lambda,
                        Br=LHS.SIR$Br,phi=LHS.SIR$phi)

write.csv(x = as.data.frame(SAmat.ind.SIR), file = "SA_SIRSMlg.csv")



LHS <- randomLHS(n,13)

LHS.PL <- data.frame(mugamma = qunif(LHS[,1],min=0,max=0.2),
                     sigma = qunif(LHS[,2],min=0,max=0.2),
                     r = qunif(LHS[,3],min=0,max=0.1),
                     mu = qunif(LHS[,4],min=0,max=0.001),
                     Bc = qunif(LHS[,5],min=0.05,max=0.55),
                     alpha = qunif(LHS[,6],min=0.25,max=1),
                     K = qunif(LHS[,7],min=1,max=500),
                     rho = qunif(LHS[,8],min=0,max=1),
                     lambda = qunif(LHS[,9],min=0,max=0.02),
                     Br = qunif(LHS[,10],min=0.02,max=0.12),
                     phi = qunif(LHS[,11],min=0,max=0.02),
                     B = qunif(LHS[,12],min=1,max=100),
                     delta = qunif(LHS[,13],min=0.1,max=0.33))


LHS.PL.o <- LHS.PL

LHS.PL <- LHS.PL[rep(1:nrow(LHS.PL),each=sims),]

SAmat.ind.P <- mapply(PLAGUE,mugamma=LHS.PL$mugamma,sigma=LHS.PL$sigma,r=LHS.PL$r,mu=LHS.PL$mu,
                      Bc=LHS.PL$Bc,alpha=LHS.PL$alpha,K=LHS.PL$K,rho=LHS.PL$rho,lambda=LHS.PL$lambda,
                      Br=LHS.PL$Br,phi=LHS.PL$phi,B=LHS.PL$B,delta=LHS.PL$delta)

write.csv(x = as.data.frame(SAmat.ind.P), file = "SA_plague.csv")


write.csv(x = as.data.frame(LHS.PL.o), file = "SA_LHS_plague.csv")
write.csv(x = as.data.frame(LHS.SIR.o), file = "SA_LHS_SIR.csv")


#Aggregate values for probability of survival (P.survive), initial average
#trait value (g), and any other output statistic of interest (we show median 
#minimum population size here).
P.survive.PL <- rep(NA,n)
N.min.PL <- rep(NA,n)
P.survive.SIR <- rep(NA,n)
N.min.SIR <- rep(NA,n)
g.PL <- rep(NA,n)
g.SIR <- rep(NA,n)
for(i in 1:n){
  N.min.PL[i] <- median(SAmat.ind.P[2,((i-1)*sims+1):(i*sims)])
  P.survive.PL[i] <- mean(SAmat.ind.P[4,((i-1)*sims+1):(i*sims)])
  N.min.SIR[i] <- median(SAmat.ind.SIR[2,((i-1)*sims+1):(i*sims)])
  P.survive.SIR[i] <- mean(SAmat.ind.SIR[4,((i-1)*sims+1):(i*sims)])
  g.PL[i] <- mean(SAmat.ind.P[5,((i-1)*sims+1):(i*sims)])
  g.SIR[i] <- mean(SAmat.ind.SIR[5,((i-1)*sims+1):(i*sims)])
}


#Run simulations for 'no evolution' using the aggregated g values just computed
SAmat.ind.SIR.0 <- mapply(SIR.sto,mugamma=rep(g.SIR,each=sims),sigma=rep(0,length(LHS.SIR$sigma)),r=LHS.SIR$r,mu=LHS.SIR$mu,
                          B=LHS.SIR$B,alpha=LHS.SIR$alpha,K=LHS.SIR$K,rho=LHS.SIR$rho,lambda=LHS.SIR$lambda,
                          Br=LHS.SIR$Br,phi=LHS.SIR$phi)
write.csv(x = as.data.frame(SAmat.ind.SIR.0), file = "SA_SIRSMlg_0.csv")


SAmat.ind.P.0 <- mapply(PLAGUE,mugamma=rep(g.PL,each=sims),sigma=rep(0,length(LHS.PL$sigma)),r=LHS.PL$r,mu=LHS.PL$mu,
                        Bc=LHS.PL$Bc,alpha=LHS.PL$alpha,K=LHS.PL$K,rho=LHS.PL$rho,lambda=LHS.PL$lambda,
                        Br=LHS.PL$Br,phi=LHS.PL$phi,B=LHS.PL$B,delta=LHS.PL$delta)
write.csv(x = as.data.frame(SAmat.ind.P.0), file = "SA_plague_0.csv")

#Aggregate values for 'no evolution' scenarios
P.survive.PL.0 <- rep(NA,n)
N.min.PL.0 <- rep(NA,n)
P.survive.SIR.0 <- rep(NA,n)
N.min.SIR.0 <- rep(NA,n)
for(i in 1:n){
  # N.min.PL.0[i] <- median(SAmat.ind.P.0[2,((i-1)*sims+1):(i*sims)])
  P.survive.PL.0[i] <- mean(SAmat.ind.P.0[4,((i-1)*sims+1):(i*sims)])
  # N.min.SIR.0[i] <- median(SAmat.ind.SIR.0[2,((i-1)*sims+1):(i*sims)])
  P.survive.SIR.0[i] <- mean(SAmat.ind.SIR.0[4,((i-1)*sims+1):(i*sims)])
}


#determine probability of Evolutionary Rescue for each parameter set
ER.PL <- P.survive.PL - P.survive.PL.0
ER.SIR <- P.survive.SIR - P.survive.SIR.0

param.PL <- 13
param.SIR <- 11

#Retrieve original non-repeated LHS matrices from above
LHS.PL.o <- LHS.PL[which(duplicated(LHS.PL)==FALSE),]
LHS.SIR.o <- LHS.SIR[which(duplicated(LHS.SIR)==FALSE),]



########################################
#Perform PRCC sensitivity analysis
########################################


#sensitivity analysis for reduced model
Y.SIR <- pcc(LHS.SIR.o,ER.SIR,rank=TRUE)
T.SIR <- Y.SIR$PRCC*sqrt((n-2-param.SIR)/(1-Y.SIR$PRCC^2))
qt(c(0.025,0.975),n-2-param.SIR) #95% confidence interval
T.SIR

#sensitivity analysis for system model
Y.PL <- pcc(LHS.PL.o,ER.PL,rank=TRUE)
T.PL <- Y.PL$PRCC*sqrt((n-2-param.PL)/(1-Y.PL$PRCC^2))
qt(c(0.025,0.975),n-2-param.PL) #95% confidence interval
T.PL

