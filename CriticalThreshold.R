rm(list=ls())

setwd("C:/Users/webblab/Documents/Plague")

library(tidyverse)
library(latex2exp)
library(viridis)
library(cowplot)
packs <- c("png","grid")
lapply(packs, require, character.only=TRUE)


#Make group result averages
ind.results <- rbind(read.csv("Fig1ResultsSIRSMloggrowth_store1_20.csv"),
                     read.csv("Fig1ResultsSIRSMloggrowth_store2_100.csv"),
                     read.csv("Fig1ResultsSIRSMloggrowth_store3_100.csv"))
ind.results <- arrange(ind.results,g,sigma)
rep <- 20+100+100


SIRMlg <- NA
for(i in 1:(length(ind.results[,1])/rep)){
  g.1 <- mean(ind.results$g.1[(1+rep*(i-1)):(rep*i)])
  min.N <- median(ind.results$min.N[(1+rep*(i-1)):(rep*i)])
  p.survive <- mean(ind.results$survive[(1+rep*(i-1)):(rep*i)])
  evolve <- mean(ind.results$evolve[(1+rep*(i-1)):(rep*i)])
  SIRMlg <- cbind(SIRMlg,c(p.survive,ind.results$g[rep*i],ind.results$sigma[rep*i],evolve,g.1,min.N))
}
SIRMlg <- as.data.frame(t(SIRMlg[,-1]))
SIRMlg <- cbind(SIRMlg,NA,NA)
colnames(SIRMlg) <- c("survive","trait","sd","evolve","g.1","min.N","spline","ER")

#Add spline base values and "ER" potential
spline <- with(SIRMlg[which(SIRMlg$sd==0),],smooth.spline(x=g.1,y=survive))
SIR.spl <- predict(spline,x=SIRMlg$g.1,deriv=0)
SIR.spl.y <- SIR.spl$y
SIRMlg$spline[1:length(SIR.spl.y)] <- SIR.spl.y
SIRMlg$ER <- SIRMlg$survive-as.numeric(SIR.spl$y)



p4 <- ggplot(SIRMlg,aes(g.1,min.N,group=sd)) + 
  xlab(TeX("Mean initial trait value $\\bar{\\gamma}_0$")) + ylab("Median minimum population") +
  geom_point(aes(g.1,min.N,color=survive,size=ER)) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  scale_color_viridis(option="magma",direction=-1,name="Proportion\nsurviving",limits=c(0,1),values=c(0,0.94,1)) +
  geom_point(data=SIRMlg[which(SIRMlg$sd==0.0005 & SIRMlg$trait==0.005),],aes(g.1,min.N),color="red3",size=5.5,pch=16) +
  geom_point(data=SIRMlg[which(SIRMlg$sd==0.125 & SIRMlg$trait==0.170),],aes(g.1,min.N),color="dodgerblue",size=7,pch=18) +
  geom_text(data=SIRMlg[which(SIRMlg$sd==0.0005 & SIRMlg$trait==0.005),],aes(g.1,min.N,label="P"),color="white",size=4,fontface="bold") +
  geom_text(data=SIRMlg[which(SIRMlg$sd==0.125 & SIRMlg$trait==0.170),],aes(g.1,min.N,label="G"),color="white",size=4,fontface="bold") +
  guides(size=guide_legend(order=1), col=guide_colorbar(order=2))

p4






#Averaged across sigma

crit.thresh <- as.data.frame(cbind(rep(0:200,length(unique(ind.results$sigma))),NA, rep(unique(ind.results$sigma),each=201)))
colnames(crit.thresh) <- c("N","survive","sigma")

for(i in 1:201){
  for(j in 1:length(unique(ind.results$sigma))){
    crit.thresh[(j-1)*201+i,2] <- mean(ind.results$survive[which(ind.results$sigma==unique(ind.results$sigma)[j] & ind.results$min.N<=(i-1))])
  }
}


p.new <- ggplot(crit.thresh, aes(N,survive)) +
  xlab("Population size reached") + ylab("Probability of survival") +
  geom_point(aes(color=sigma)) +
  geom_line(aes(group=sigma)) +
  theme_classic() + xlim(0,100) +
  scale_color_viridis(option="viridis",direction=-1,name="Variance")
p.new



#Averaged across mu.gamma, not sigma

crit.thresh.2 <- as.data.frame(cbind(rep(0:200,length(unique(ind.results$g))),NA, rep(unique(ind.results$g),each=201)))
colnames(crit.thresh.2) <- c("N","survive","mu.gamma")

for(i in 1:201){
  for(j in 1:length(unique(ind.results$g))){
    crit.thresh.2[(j-1)*201+i,2] <- mean(ind.results$survive[which(ind.results$g==unique(ind.results$g)[j] & ind.results$min.N<=(i-1))])
  }
}


p.new.2 <- ggplot(crit.thresh.2, aes(N,survive)) +
  xlab("Population size reached") + ylab("Probability of survival") +
  geom_point(aes(color=mu.gamma)) +
  geom_line(aes(group=mu.gamma)) +
  theme_classic() + xlim(0,100) +
  scale_color_viridis(option="viridis",direction=-1,name="mu.gamma")
p.new.2


#Averaged across mu.gamma and sigma

crit.thresh.3 <- as.data.frame(cbind(rep(0:100,length(unique(ind.results$g))*length(unique(ind.results$sigma))),
                                     NA, 
                                     rep(unique(ind.results$g),each=101*length(unique(ind.results$sigma))),
                                     rep(unique(ind.results$sigma),101*length(unique(ind.results$g)))))
colnames(crit.thresh.3) <- c("N","survive","mu.gamma","sigma")

for(i in 1:201){
  for(j in 1:length(unique(ind.results$g))){
    for(k in 1:length(unique(ind.results$sigma))){
      crit.thresh.3[which(crit.thresh.3$mu.gamma==unique(ind.results$g)[j] & crit.thresh.3$sigma==unique(ind.results$sigma)[k])[i],2] <- mean(ind.results$survive[which(ind.results$g==unique(ind.results$g)[j] & ind.results$sigma==unique(ind.results$sigma)[k] & ind.results$min.N<=(i-1))])
    }
  }
}


p.new.3 <- ggplot(crit.thresh.3, aes(N,survive)) +
  xlab("Population size reached") + ylab("Probability of survival") +
  geom_point(aes(color=sigma, size=mu.gamma)) +
  theme_classic() +
  scale_color_viridis(option="viridis",direction=-1,name="sigma")
p.new.3
