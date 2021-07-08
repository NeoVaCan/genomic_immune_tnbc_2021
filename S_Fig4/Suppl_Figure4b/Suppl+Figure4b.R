rm(list=ls())

library(Hmisc)
library(ggplot2)
library(ggpubr)


Purity=read.table(file="Purity.txt", header = TRUE, sep = "\t" )

Purity<-as.data.frame(Purity[-c(1),])


pdf("Purity_boxplot.pdf", useDingbats = FALSE, height = 4, width = 4)
ggplot(Purity,aes(x="",y=Purity))+
  geom_boxplot()+
  geom_jitter(aes(x="",y=Purity), col = "red")+
  labs(x="Parallel metastases")+
  theme_classic()
dev.off()
