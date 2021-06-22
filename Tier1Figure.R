

rm(list=ls())
setwd("C:/Users/webblab/Documents/Plague/SAMP/")
library(latex2exp)
library(tidyverse)


####Tier 1 processes

files <- c("Tier1SIR.csv","Tier1SIRlg.csv","Tier1SISR.csv","Tier1SEIR.csv",
           "Tier1SIRS.csv","Tier1SIRM.csv")

df <- NA
for(i in 1:length(files)){
  x <- read.csv(files[i])
  x <- cbind(x[,-1],gsub(".csv","",files[i]))
  df <- rbind(df,x)
}
df <- df[-1,]
colnames(df) <- c("p.survive","g","max.evolve","g.1","Model")
df$Model <- gsub("Tier1","",df$Model)
mypalette <- c("black","#56B4E9","#0072B2","#E69F00","#D55E00","#CC79A7")
df$Model <- factor(df$Model, levels=c("SIR","SIRlg","SISR","SEIR","SIRS","SIRM"))

ggplot(df, aes(as.numeric(as.character(g.1)),as.numeric(as.character(p.survive)))) + theme_minimal() + 
  xlab(TeX("Mean initial resistance rate $\\bar{\\gamma}_0$")) +
  ylab("Proportion surviving populations") +
  geom_point(aes(color=Model, shape=Model), size=3.5) +
  geom_line(aes(group=Model, color=Model), size=1.5) +
  scale_color_manual(values=mypalette, labels=c("SIR","Logistic","Split","Exposed","Impermanent","Reservoir")) +
  scale_shape_manual(values=c(16,15,17,18,4,8), labels=c("SIR","Logistic","Split","Exposed","Impermanent","Reservoir")) +
  geom_line(data=df[which(df$Model=="SIR"),], aes(group=1), color="black", size=1.5) +
  geom_point(data=df[which(df$Model=="SIR"),], aes(group=1), color="black", size=3.5, shape=16) +
  theme(legend.background = element_rect(fill="white", linetype="solid", color="black"),
        text = element_text(size=16),
        legend.position = c(0.8,0.6))



####Tier 2 processes

files <- c("Tier1SIR.csv","Tier2SEIRM.csv","Fig1ResultsSIRSM.csv")

df2 <- NA
for(i in 1:2){
  x <- read.csv(files[i])
  x <- cbind(x[,-1],gsub(".csv","",files[i]))
  df2 <- rbind(df2,x)
}
df2 <- df2[-1,]
colnames(df2) <- c("p.survive","g","max.evolve","g.1","Model")
df2$Model <- gsub("Tier1","",df2$Model)
df2$Model <- gsub("Tier2","",df2$Model)

x <- read.csv(files[3])
x <- x[which(x$sigma==0),]
rep <- length(x[,1])/length(unique(x$g))
for(i in 1:length(unique(x$g))){
  df2 <- rbind(df2,c(mean(x$survive[((i-1)*rep+1):(i*rep)]),x$g[i*rep],x$evolve[i*rep],x$g.1[i*rep],"SIRSM"))
}
df2$Model <- factor(df2$Model, levels=c("SIR","SEIRM","SIRSM"))
mypalette <- c("black","dimgray","snow4")

ggplot(df2, aes(as.numeric(as.character(g.1)),as.numeric(as.character(p.survive)))) + theme_minimal() + 
  xlab(TeX("Mean initial resistance rate $\\bar{\\gamma}_0$")) +
  ylab("Proportion of surviving populations") +
  geom_line(aes(group=Model, color=Model), size=1.5) +
  geom_point(aes(shape=Model, color=Model), size=3.5) +
  scale_color_manual(values=mypalette, labels=c("SIR","Reservoir-Exposed","Reservoir-Impermanent")) +
  scale_shape_manual(values=c(16,15,17), labels=c("SIR","Reservoir-Exposed","Reservoir-Impermanent")) +
  geom_line(data=df2[which(df2$Model=="SIR"),], aes(group=1), color="black", size=1.5) +
  geom_point(data=df2[which(df2$Model=="SIR"),], aes(group=1), color="black", size=3.5, shape=16) +
  theme(legend.background = element_rect(fill="white", linetype="solid", color="black"),
        text = element_text(size=16),
        legend.position = c(0.8,0.6))



####Tier 3 processes

setwd("C:/Users/webblab/Documents/Plague/SAMP/")
files <- c("Tier1SIR.csv","Tier3SIMSRS.csv","Tier2SIRMlg.csv")

df3 <- NA
for(i in 1:1){
  x <- read.csv(files[i])
  x <- cbind(x[,-1],gsub(".csv","",files[i]))
  df3 <- rbind(df3,x)
}
df3 <- df3[-1,]
colnames(df3) <- c("p.survive","g","max.evolve","g.1","Model")
df3$Model <- gsub("Tier1","",df3$Model)

for(i in 2:2){
  x <- read.csv(files[i])
  x <- cbind(x[,-1],gsub(".csv","",files[i]))
  colnames(x) <- c("p.survive","g","max.evolve","g.1","Model")
  df3 <- rbind(df3,x)
}
df3$Model <- gsub("Tier3","",df3$Model)

df3$Model <- factor(df3$Model, levels=c("SIR","SIMSRS"))
mypalette <- c("black","dimgray")

ggplot(df3, aes(as.numeric(as.character(g.1)),as.numeric(as.character(p.survive)))) + theme_minimal() + 
  xlab(TeX("Mean initial trait value $\\bar{\\gamma}_0$")) +
  ylab("Proportion of surviving populations") +
  geom_line(aes(group=Model, color=Model), size=1.5) +
  geom_point(aes(shape=Model, color=Model), size=3.5) +
  scale_color_manual(values=mypalette, labels=c("SIR","Reservoir-Impermanent-Split")) +
  scale_shape_manual(values=c(16,15), labels=c("SIR","Reservoir-Impermanent-Split")) +
  geom_line(data=df3[which(df3$Model=="SIR"),], aes(group=1), color="black", size=1.5) +
  geom_point(data=df3[which(df3$Model=="SIR"),], aes(group=1), color="black", size=3.5, shape=16) +
  theme(legend.background = element_rect(fill="white", linetype="solid", color="black"),
        text = element_text(size=16),
        legend.position = c(0.8,0.6))



###Back-check

setwd("C:/Users/webblab/Documents/Plague/SAMP/")
files <- c("Tier1SIR.csv","Tier3SIMSRS.csv","Tier2SIRMlg.csv")

df3 <- NA
for(i in 1:1){
  x <- read.csv(files[i])
  x <- cbind(x[,-1],gsub(".csv","",files[i]))
  df3 <- rbind(df3,x)
}
df3 <- df3[-1,]
colnames(df3) <- c("p.survive","g","max.evolve","g.1","Model")
df3$Model <- gsub("Tier1","",df3$Model)
setwd("C:/Users/webblab/Documents/Plague/")

#Make group result averages
ind.results <- rbind(read.csv("Fig1ResultsSIRMloggrowth_store1_5.csv"),
                     read.csv("Fig1ResultsSIRMloggrowth_store2_20.csv"),
                     read.csv("Fig1ResultsSIRMloggrowth_store3_175.csv"),
                     read.csv("Fig1ResultsSIRMloggrowth_store4_400.csv"))
ind.results <- arrange(ind.results,g,sigma)
rep <- 5+20+175+400

SIRMlg <- NA
for(i in 1:(length(ind.results[,1])/rep)){
  g.1 <- mean(ind.results$g.1[(1+rep*(i-1)):(rep*i)])
  min.N <- mean(ind.results$min.N[(1+rep*(i-1)):(rep*i)])
  p.survive <- mean(ind.results$survive[(1+rep*(i-1)):(rep*i)])
  evolve <- mean(ind.results$evolve[(1+rep*(i-1)):(rep*i)])
  SIRMlg <- cbind(SIRMlg,c(p.survive,ind.results$g[rep*i],ind.results$sigma[rep*i],evolve,g.1,min.N))
}
SIRMlg <- as.data.frame(t(SIRMlg[,-1]))
SIRMlg <- cbind(SIRMlg,NA,NA)
colnames(SIRMlg) <- c("survive","trait","sd","evolve","g.1","min.N","spline","ER")
SIRMlg <- SIRMlg[which(SIRMlg$sd==0),]

for(i in 1:length(SIRMlg[,1])){
  df3 <- rbind(df3,c(SIRMlg$survive[i],SIRMlg$g[i],SIRMlg$evolve[i],SIRMlg$g.1[i],"SIRMlg"))
}

ind.results <- rbind(read.csv("Fig1ResultsSIRSloggrowth_store1_10.csv"),
                     read.csv("Fig1ResultsSIRSloggrowth_store2_100.csv"),
                     read.csv("Fig1ResultsSIRSloggrowth_store3_500.csv"),
                     read.csv("Fig1ResultsSIRSloggrowth_store4_390.csv"))
ind.results <- arrange(ind.results,g,sigma)
rep <- 1000

SIRMlg <- NA
for(i in 1:(length(ind.results[,1])/rep)){
  g.1 <- mean(ind.results$g.1[(1+rep*(i-1)):(rep*i)])
  min.N <- mean(ind.results$min.N[(1+rep*(i-1)):(rep*i)])
  p.survive <- mean(ind.results$survive[(1+rep*(i-1)):(rep*i)])
  evolve <- mean(ind.results$evolve[(1+rep*(i-1)):(rep*i)])
  SIRMlg <- cbind(SIRMlg,c(p.survive,ind.results$g[rep*i],ind.results$sigma[rep*i],evolve,g.1,min.N))
}
SIRMlg <- as.data.frame(t(SIRMlg[,-1]))
SIRMlg <- cbind(SIRMlg,NA,NA)
colnames(SIRMlg) <- c("survive","trait","sd","evolve","g.1","min.N","spline","ER")
SIRMlg <- SIRMlg[which(SIRMlg$sd==0),]

for(i in 1:length(SIRMlg[,1])){
  df3 <- rbind(df3,c(SIRMlg$survive[i],SIRMlg$g[i],SIRMlg$evolve[i],SIRMlg$g.1[i],"SIRSlg"))
}


df3$Model <- factor(df3$Model, levels=c("SIR","SIRMlg","SIRSlg"))
mypalette <- c("black","dimgray","snow4")

ggplot(df3, aes(as.numeric(as.character(g.1)),as.numeric(as.character(p.survive)))) + theme_minimal() + 
  xlab(TeX("Mean initial trait value $\\bar{\\gamma}_0$")) +
  ylab("Proportion of surviving populations") +
  geom_line(aes(group=Model, color=Model), size=1.5) +
  geom_point(aes(shape=Model, color=Model), size=3.5) +
  scale_color_manual(values=mypalette, labels=c("SIR","Reservoir-Logistic","Impermanent-Logistic")) +
  scale_shape_manual(values=c(16,17,8), labels=c("SIR","Reservoir-Logistic","Impermanent-Logistic")) +
  geom_line(data=df3[which(df3$Model=="SIR"),], aes(group=1), color="black", size=1.5) +
  geom_point(data=df3[which(df3$Model=="SIR"),], aes(group=1), color="black", size=3.5, shape=16) +
  theme(legend.background = element_rect(fill="white", linetype="solid", color="black"),
        text = element_text(size=16),
        legend.position = c(0.8,0.6))

