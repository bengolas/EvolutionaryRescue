rm(list=ls())
setwd("C:/Users/webblab/Documents/Plague/SAMP/")
library(latex2exp)
library(tidyverse)


####Tier 1 SISR comparison
setwd("C:/Users/webblab/Documents/Plague/SAMP/")
files <- c("Tier1SIR.csv","Tier1SISR.csv","Tier1SISRdelta1.csv")

df <- NA
for(i in 1:length(files)){
  x <- read.csv(files[i])
  x <- cbind(x[,-1],gsub(".csv","",files[i]))
  df <- rbind(df,x)
}
df <- df[-1,]
colnames(df) <- c("p.survive","g","max.evolve","g.1","Model")
df$Model <- gsub("Tier1","",df$Model)
mypalette <- c("black","red","blue")
df$Model <- factor(df$Model, levels=c("SIR","SISR","SISRdelta1"))

ggplot(df, aes(as.numeric(as.character(g.1)),as.numeric(as.character(p.survive)))) + theme_minimal() + 
  xlab(TeX("Mean initial trait value $\\bar{\\gamma}_0$")) +
  ylab("Proportion surviving populations") +
  geom_point(aes(color=Model, shape=Model), size=3.5) +
  geom_line(aes(group=Model, color=Model), size=1.5) +
  scale_color_manual(values=mypalette, labels=c("SIR",expression(paste("Split, ",delta,"=0.21")),expression(paste("Split, ",delta,"=1")))) +
  scale_shape_manual(values=c(16,17,15), labels=c("SIR",expression(paste("Split, ",delta,"=0.21")),expression(paste("Split, ",delta,"=1")))) +
  geom_line(data=df[which(df$Model=="SIR"),], aes(group=1), color="black", size=1.5) +
  geom_point(data=df[which(df$Model=="SIR"),], aes(group=1), color="black", size=3.5, shape=16) +
  theme(legend.background = element_rect(fill="white", linetype="solid", color="black"),
        text = element_text(size=16),
        legend.position = c(0.8,0.6))


