rm(list=ls())

setwd("C:/Users/webblab/Documents/Plague")

library(tidyverse)
library(latex2exp)
library(viridis)
library(cowplot)
library(ggnewscale)
library(ghibli)


ind.results <- rbind(read.csv("Fig1ResultsSIRSMloggrowth_store1_20.csv"),
                     read.csv("Fig1ResultsSIRSMloggrowth_store2_100.csv"),
                     read.csv("Fig1ResultsSIRSMloggrowth_store3_100.csv"))
ind.results <- arrange(ind.results,g,sigma)
rep <- length(ind.results$g)/(length(unique(ind.results$g))*length(unique(ind.results$sigma)))

n <- 201

store <- matrix(NA,ncol=3)
store <- as.data.frame(store)

for(i in 1:length(unique(ind.results$g))){
  df <- ind.results[which(ind.results$g==unique(ind.results$g)[i]),]
  for(j in 1:n){
    store[nrow(store)+1,] <- c(mean(df$survive[which(df$min.N<=(j-1))]),(j-1),unique(ind.results$g)[i])
  }
}
store <- store[-1,]
colnames(store) <- c("Survive","pop.size","sigma")


ggplot(store, aes(pop.size,Survive,group=sigma, color=sigma)) + 
  xlab("Population size reached") + ylab("Probability of survival") +
  theme_classic() +
  theme(text = element_text(size=16)) +
  geom_line() +
  geom_point() +
  scale_color_viridis_c(option="viridis", direction=-1,name=TeX("Resistance rate $\\mu_{\\gamma}$"))
