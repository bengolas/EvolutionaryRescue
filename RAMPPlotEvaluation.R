rm(list=ls())

setwd("C:/Users/webblab/Documents/Plague")

library(tidyverse)
library(latex2exp)
library(viridis)
library(cowplot)
packs <- c("png","grid")
lapply(packs, require, character.only=TRUE)


#Make group result averages
# ind.results <- rbind(read.csv("Fig1ResultsSIR_store1_5.csv"),
#                      read.csv("Fig1ResultsSIR_store2_195.csv"),
#                      read.csv("Fig1ResultsSIR_store3_300.csv"))
ind.results <- rbind(read.csv("Fig1ResultsRtoS.csv"))
ind.results <- arrange(ind.results,g,variance)
rep <- length(ind.results[,1])/(length(unique(ind.results$g))*length(unique(ind.results$variance)))

SIRMlg <- NA
for(i in 1:(length(ind.results[,1])/rep)){
  g.1 <- mean(ind.results$g.1[(1+rep*(i-1)):(rep*i)])
  min.N <- median(ind.results$min.N[(1+rep*(i-1)):(rep*i)])
  p.survive <- mean(ind.results$survive[(1+rep*(i-1)):(rep*i)])
  evolve <- mean(ind.results$evolve[(1+rep*(i-1)):(rep*i)])
  SIRMlg <- cbind(SIRMlg,c(p.survive,ind.results$g[rep*i],ind.results$variance[rep*i],evolve,g.1,min.N))
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



p <- ggplot(SIRMlg,aes(g.1,survive,group=sd)) + geom_line(aes(color=sd),size=1.25)  + 
  xlab(TeX("Mean initial trait value $\\bar{\\gamma}_0$")) + ylab(TeX("Proportion surviving populations")) +
  geom_point(aes(g.1,survive,color=sd),size=3.5) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  ylim(-0.01,1.01) +
  scale_color_viridis(option="viridis",direction=-1,name=TeX("Variance $\\sigma^2$"), labels=c("0","0.0025","0.01","0.0225","0.04")) +
  geom_line(data=SIRMlg[which(SIRMlg$sd==0),],aes(g.1,spline),size=1.5) +
  geom_point(data=SIRMlg[which(SIRMlg$sd==0.0005 & SIRMlg$trait==0.005),],aes(g.1,survive),color="red3",size=5.5,pch=16) +
  geom_point(data=SIRMlg[which(SIRMlg$sd==0.125 & SIRMlg$trait==0.170),],aes(g.1,survive),color="dodgerblue",size=7,pch=18) +
  geom_text(data=SIRMlg[which(SIRMlg$sd==0.0005 & SIRMlg$trait==0.005),],aes(g.1,survive,label="P"),color="white",size=4,fontface="bold") +
  geom_text(data=SIRMlg[which(SIRMlg$sd==0.125 & SIRMlg$trait==0.170),],aes(g.1,survive,label="G"),color="white",size=4,fontface="bold") 


setwd("C:/Users/webblab/Documents/Plague")
print(p)
