setwd("C:/Users/bengo_000/Documents/Plague/CODING")

library(tidyverse)
library(cowplot)
library(reshape2)
library(viridis)
library(latex2exp)
packs <- c("png","grid")
lapply(packs, require, character.only = TRUE) 

BTPD <- as.data.frame(read.csv("BTPD_example.csv"))

BTPD$N[which(BTPD$N>=30)] <- ((BTPD$N[which(BTPD$N>=30)]-30)/170)*5+30

BTPD <- melt(BTPD,id.vars=c("WT","g.avg"),measure.vars=c("M","R","I","E","S","N"))

CAGS <- as.data.frame(read.csv("CAGS_example3.csv"))

CAGS <- melt(CAGS,id.vars=c("WT","g.avg"),measure.vars=c("M","R","I","E","S","N"))


img <- readPNG("InsetPL.png")
g <- rasterGrob(img, interpolate=TRUE)

plasma_pal <- c("yellow3",rev(viridis::magma(n=5)[-5]))

p1 <- ggplot(BTPD) + theme_classic() + 
  theme(legend.position = "none",text = element_text(size = 16)) +
  xlab("Time (days)") + ylab("Number of individuals") +
  geom_line(data=BTPD[which(BTPD$variable!="N"),], aes(WT,value,group=variable,color=variable),size=1.5) +
  scale_color_manual(values=plasma_pal) +
  annotation_custom(g, xmin=75, xmax=250, ymin=50, ymax=225)

p2 <- ggplot(CAGS) + theme_classic() + 
  theme(legend.title = element_blank(),text = element_text(size = 16)) +
  xlab("Time (days)") + ylab("Number of individuals") +
  geom_line(data=CAGS[which(CAGS$variable!="N"),], aes(WT,value,group=variable,color=variable),size=1.5) +
  scale_color_manual(values=plasma_pal) +
  guides(color = guide_legend(reverse = TRUE))

p3 <- ggplot() + theme_classic() + 
  theme(legend.position = "none",text = element_text(size = 16)) +
  xlab("Time (days)") + ylab("Number of individuals") + 
  geom_line(data=BTPD[which(BTPD$variable=="N"),], aes(WT,value),color="red3",size=1.5) +
  geom_line(data=CAGS[which(CAGS$variable=="N"),], aes(WT,value),color="dodgerblue",size=1.5) +
  scale_y_continuous(breaks=c(0,10,20,30,35), labels=c("0","10","20","30","200"))

p4 <- ggplot() + theme_classic() + 
  theme(legend.title = element_blank(),text = element_text(size = 16),
        legend.key = element_rect(color="white",size = 5),
        legend.key.size = unit(1.5, 'lines')) +
  xlab("Time (days)") + 
  ylab(TeX("Relative change in mean resistance rate")) +
  geom_line(data=BTPD[which(BTPD$variable=="N"),], aes(WT,(g.avg-g.avg[1])/0.0005,color="red3"),size=1.5) +
  geom_line(data=CAGS[which(CAGS$variable=="N"),], aes(WT,(g.avg-g.avg[1])/0.1,color="dodgerblue"),size=1.5) +
  scale_color_manual(name="species",values=c("dodgerblue"="dodgerblue","red3"="red3"), labels=c("Ground\nsquirrel","Prairie\ndog")) +
  guides(color = guide_legend(reverse = TRUE))



plot_grid(p1,p2,p3,p4,labels="AUTO")

