rm(list=ls())

setwd("C:/Users/bengo_000/Documents/Plague/CODING")

library(tidyverse)
library(latex2exp)
library(viridis)
library(cowplot)
packs <- c("png","grid")
lapply(packs, require, character.only=TRUE)

SGVonly <- as.data.frame(t(read.csv("PlagueGroupDataSmall10000_onlySGV.csv"))[-1,])
colnames(SGVonly) <- c("Nt","min.N","trait","g.1","evolve","survivor.evolve","survival","variance")
Plague <- as.data.frame(t(read.csv("PlagueGroupDataSmall10000.csv"))[-1,])
colnames(Plague) <- c("Nt","min.N","trait","g.1","evolve","survivor.evolve","survival","variance")

df <- as.data.frame(cbind(SGVonly$trait, SGVonly$g.1, SGVonly$variance, Plague$survival-SGVonly$survival))
colnames(df) <- c("trait","g.1","variance","survive")

ggplot(df) + 
  xlim(0,0.06) + ylim(0,0.025) +
  xlab(TeX("Mean initial trait value $\\bar{\\gamma}_0$")) + ylab(TeX("Variance $\\sigma^2$")) +
  geom_point(aes(g.1,variance^2,color=survive),size=3.5) +
  theme_classic() +
  theme(text = element_text(size=16)) +
  scale_color_viridis(option="magma",direction=1,name="Difference\nin survival",limits=c(-0.15,1.15))
