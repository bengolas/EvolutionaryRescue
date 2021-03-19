###############################################################################
#The following code will define the stochastic plague model system function 
#and generate data from simulations. It runs a series of *rep* iterations 
#across ranges of average initial trait size (mu_gamma) 
#and variance (sd == standard deviation of genetic variance, 
#i.e. sd^2=variance).
#Throughout the script, 'g' is shorthand for 'gamma'.
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

reps <- 2 #number of iterations per parameterization

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
