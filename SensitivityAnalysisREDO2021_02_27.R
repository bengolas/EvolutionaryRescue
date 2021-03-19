#SA uses functions from file SensitivityAnalysisFunctions to generate 
#sensitivity analysis results for SIRSMlg and plague models

#Get survival and g.1 from first run, then run again with same parameters but sigma=0 and g.start=mean(g.1)

library(fitur)
library(sensitivity)

set.seed(19)

n <- 200 #number of parameter combinations
sims <- 200 #number of simulations per parameter combination

# LHS.SIR <- data.frame(g.start = runif(n,0,0.2),
#                       sigma = runif(n,0,0.2),
#                       r = runif(n,0,0.1),
#                       mu = runif(n,0,0.001),
#                       B = runif(n,0.05,5.5),
#                       alpha = runif(n,0.25,1),
#                       K = rdunif(n,1,500),
#                       rho = runif(n,0,1),
#                       lambda = runif(n,0,0.02),
#                       Br = runif(n,0.02,0.12),
#                       phi = runif(n,0,0.02))
# 
# 
# LHS.SIR.o <- LHS.SIR
# 
# LHS.SIR <- LHS.SIR[rep(1:nrow(LHS.SIR),each=sims),]
# 
# SAmat.ind.SIR <- mapply(SIR.sto,g.start=LHS.SIR$g.start,sigma=LHS.SIR$sigma,r=LHS.SIR$r,mu=LHS.SIR$mu,
#        B=LHS.SIR$B,alpha=LHS.SIR$alpha,K=LHS.SIR$K,rho=LHS.SIR$rho,lambda=LHS.SIR$lambda,
#        Br=LHS.SIR$Br,phi=LHS.SIR$phi)
# 
# write.csv(x = as.data.frame(SAmat.ind.SIR), file = "SA_SIRSMlgREDO.csv")


LHS.PL <- data.frame(g.start = runif(n,0,0.2),
                      sigma = runif(n,0,0.2),
                      r = runif(n,0,0.1),
                      mu = runif(n,0,0.001),
                      Bc = runif(n,0.05,0.55),
                      alpha = runif(n,0.25,1),
                      K = rdunif(n,1,500),
                      rho = runif(n,0,1),
                      lambda = runif(n,0,0.02),
                      Br = runif(n,0.02,0.12),
                      phi = runif(n,0,0.02),
                      B = runif(n,0,100),
                      period = runif(n,0.1,0.33))


LHS.PL.o <- LHS.PL

LHS.PL <- LHS.PL[rep(1:nrow(LHS.PL),each=sims),]

SAmat.ind.P <- mapply(PLAGUE,g.start=LHS.PL$g.start,sigma=LHS.PL$sigma,r=LHS.PL$r,mu=LHS.PL$mu,
                    Bc=LHS.PL$Bc,alpha=LHS.PL$alpha,K=LHS.PL$K,rho=LHS.PL$rho,lambda=LHS.PL$lambda,
                    Br=LHS.PL$Br,phi=LHS.PL$phi,B=LHS.PL$B,period=LHS.PL$period)

write.csv(x = as.data.frame(SAmat.ind.P), file = "SA_plagueREDO2b.csv")


write.csv(x = as.data.frame(LHS.PL.o), file = "SA_LHS_plagueREDO2b.csv")
# write.csv(x = as.data.frame(LHS.SIR.o), file = "SA_LHS_SIRREDO.csv")


#Aggregate values
P.survive.PL <- rep(NA,n)
N.min.PL <- rep(NA,n)
# P.survive.SIR <- rep(NA,n)
# N.min.SIR <- rep(NA,n)
g.PL <- rep(NA,n)
# g.SIR <- rep(NA,n)
for(i in 1:n){
  N.min.PL[i] <- median(SAmat.ind.P[2,((i-1)*sims+1):(i*sims)])
  P.survive.PL[i] <- mean(SAmat.ind.P[4,((i-1)*sims+1):(i*sims)])
  # N.min.SIR[i] <- median(SAmat.ind.SIR[2,((i-1)*sims+1):(i*sims)])
  # P.survive.SIR[i] <- mean(SAmat.ind.SIR[4,((i-1)*sims+1):(i*sims)])
  g.PL[i] <- mean(SAmat.ind.P[5,((i-1)*sims+1):(i*sims)])
  # g.SIR[i] <- mean(SAmat.ind.SIR[5,((i-1)*sims+1):(i*sims)])
}


#Run simulations for 'no evolution' g counterparts to above
# SAmat.ind.SIR.0 <- mapply(SIR.sto,g.start=rep(g.SIR,each=sims),sigma=rep(0,length(LHS.SIR$sigma)),r=LHS.SIR$r,mu=LHS.SIR$mu,
#                         B=LHS.SIR$B,alpha=LHS.SIR$alpha,K=LHS.SIR$K,rho=LHS.SIR$rho,lambda=LHS.SIR$lambda,
#                         Br=LHS.SIR$Br,phi=LHS.SIR$phi)
# write.csv(x = as.data.frame(SAmat.ind.SIR.0), file = "SA_SIRSMlg_0REDO.csv")


SAmat.ind.P.0 <- mapply(PLAGUE,g.start=rep(g.PL,each=sims),sigma=rep(0,length(LHS.PL$sigma)),r=LHS.PL$r,mu=LHS.PL$mu,
                      Bc=LHS.PL$Bc,alpha=LHS.PL$alpha,K=LHS.PL$K,rho=LHS.PL$rho,lambda=LHS.PL$lambda,
                      Br=LHS.PL$Br,phi=LHS.PL$phi,B=LHS.PL$B,period=LHS.PL$period)
write.csv(x = as.data.frame(SAmat.ind.P.0), file = "SA_plague_0REDO2b.csv")


P.survive.PL.0 <- rep(NA,n)
N.min.PL.0 <- rep(NA,n)
P.survive.SIR.0 <- rep(NA,n)
N.min.SIR.0 <- rep(NA,n)
for(i in 1:n){
  N.min.PL.0[i] <- median(SAmat.ind.P.0[2,((i-1)*sims+1):(i*sims)])
  P.survive.PL.0[i] <- mean(SAmat.ind.P.0[4,((i-1)*sims+1):(i*sims)])
  # N.min.SIR.0[i] <- median(SAmat.ind.SIR.0[2,((i-1)*sims+1):(i*sims)])
  # P.survive.SIR.0[i] <- mean(SAmat.ind.SIR.0[4,((i-1)*sims+1):(i*sims)])
}


#determine Evolutionary Rescue for each parameter set
ER.PL <- P.survive.PL - P.survive.PL.0
ER.SIR <- P.survive.SIR - P.survive.SIR.0

param.PL <- 13
param.SIR <- 11

LHS.PL.o <- LHS.PL[which(duplicated(LHS.PL)==FALSE),]
LHS.SIR.o <- LHS.SIR[which(duplicated(LHS.SIR)==FALSE),]

#sensitivity analysis SIRSMlg
Y.SIR <- pcc(LHS.SIR.o,ER.SIR,rank=TRUE)
T.SIR <- Y.SIR$PRCC*sqrt((n-2-param.SIR)/(1-Y.SIR$PRCC^2))
qt(c(0.025,0.975),n-2-param.SIR)
T.SIR

#sensitivity analysis plague
Y.PL <- pcc(LHS.PL.o,ER.PL,rank=TRUE)
T.PL <- Y.PL$PRCC*sqrt((n-2-param.PL)/(1-Y.PL$PRCC^2))
qt(c(0.025,0.975),n-2-param.PL)
T.PL


#Make a bar chart comparison

library(tidyverse)
library(reshape2)
library(latex2exp)


T.dat <- as.data.frame(rbind(T.PL,rbind(T.SIR,0,0)))
T.dat <- as.data.frame(cbind(T.dat,c(rep("Plague",13),rep("SIRSMlg",13)),rep(rownames(T.PL),2)))
colnames(T.dat) <- c("value","model","parameter")


ggplot(T.PL) + theme_classic() + xlab("") + ylab("Sensitivity") +
  geom_bar(aes(parameter,value,fill=model),stat="identity",position=position_dodge()) +
  scale_x_discrete(labels=c(expression(alpha),"B",expression(beta[c]),expression(beta[r]),expression(mu[gamma]),"K",expression(lambda),expression(mu),expression(delta),expression(phi),"r",expression(rho),expression(sigma))) +
  scale_fill_manual(name="Model",values=c("black","grey"),labels=c("Plague","SIRSMlg")) +
  geom_hline(yintercept=c(-1.973,1.973),linetype=2)


par(mfrow=c(2,2))
for(i in 1:param.PL){
  p <- plot(LHS.PL.o[,i],ER.PL,xlab=names(LHS.PL.o[i]))
  print(p)
  abline(lm(ER.PL ~ LHS.PL.o[,i]))
}

par(mfrow=c(2,2))
for(i in 1:param.SIR){
  p <- plot(LHS.SIR.o[,i],ER.SIR,xlab=names(LHS.SIR.o[i]))
  print(p)
  abline(lm(ER.SIR ~ LHS.SIR.o[,i]))
}
