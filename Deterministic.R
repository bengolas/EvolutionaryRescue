  ######################Plotting one run
  
  #Parameters
  
  r <- 0.0866       # Intrinsic rate of increase of susceptible
  mu <- 0.0002      # Natural mortality rate
  Bc <- 0.073      # Transmission coefficient
  alpha <- 0.5       # Mortality rate
  K <- 200          # Carrying Capacity
  phi <- 0.011      # Resistance loss rate
  B <- 20
  delta <- 0.21
  phi <- 0.011
  
  Br <- 0.073       #plague reservoir infection rate
  lambda <- 0.006  #reservoir decay rate
  
  #starting conditions
  
  S.start <- 199
  E.start <- 0
  I.start <- 1
  R.start <- 0
  M.start <- 0
  g.start <- 0.005  # Initial resistance trait
  sigma <- 0.0005   # Trait variance
  

  
  ###############Plotting many deterministic runs
  
  g.values <- seq(0,0.2,0.005)
  sigma.values <- c(0,0.0005,seq(0.005,0.2,0.01))
  
library(deSolve)
library(reshape)
library(viridis)
library(tidyverse)

t <- 1000
g <- 0.005

SIR <- function(t,state,parameters){
  with(as.list(c(state,parameters)), {
    dS <- r*(S+I+R)*(1-(S+I+R)/K) -mu*S -Bc*S*I*delta-Bc*S*I*g-Br*S*M*delta-Br*S*M*g +phi*R
    dI <- Bc*S*I*delta+Br*S*M*delta-mu*I -alpha*I -g*I
    dR <- Bc*S*I*g+Br*S*M*g -mu*R -phi*R
    dM <- alpha*I-lambda*M
    list(c(dS,dI,dR,dM))
  })
}

# SIR <- function(t,state,parameters){
#   with(as.list(c(state,parameters)), {
#     dS <- r*(S+R)*(1-(S+E+I+R)/K) +phi*R -mu*S -Bc*S*I/(S+E+I+R) -Br*S*M/B
#     dE <- Bc*S*I/(S+E+I+R) +Br*S*M/B -delta*E -g*E
#     dI <- delta*E -alpha*I
#     dR <- g*E -mu*R -phi*R
#     dM <- alpha*I-lambda*M
#     list(c(dS,dE,dI,dR,dM))
#   })
# }

parameters <- c(r=r,mu=mu,B=B,alpha=alpha,g=g,K=K,Br=Br,B=B,phi=phi,delta=delta,lambda=lambda)
state <- c(S=S.start,I=I.start,R=R.start,M=M.start)
times <- seq(0, t, by=1)
out <- ode(y=state, times=times, func=SIR, parms=parameters)
output <- as.data.frame(out)
output <- melt(output,id="time")


ggplot(output, aes(as.numeric(as.character(time)),as.numeric(as.character(value)))) + 
  theme_classic() + xlab("Time") + 
  ylab("Number of individuals") + ylim(0,210) + xlim(0,250) +
  geom_line(aes(group=variable, color=variable), size=1) +
  scale_color_viridis_d(option="viridis", direction=-1, name="Class")
  