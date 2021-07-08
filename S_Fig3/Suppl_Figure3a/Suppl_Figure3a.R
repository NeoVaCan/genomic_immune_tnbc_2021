rm(list=ls())

library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)
library(lemon)
library(tidyverse)


IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)

IO_mets<-subset(IO,(Response %in% "Autopsy"))

#Remove unused rows
IO_mets<-IO_mets[-c(5),]

Cell.type.boxplot <- IO_mets[,c('Tissue', 'APM','IFN.Gamma','Prolif','TIS','B.Cells','T.Cells','CD45', 'CD8.T.Cells','Exhausted.CD8','Treg','Macrophages','DC' ,'NK.Cells','NK.CD56dim.Cells', 'Mast.Cells','Neutrophils', 'Th1.Cells', 'MHC2')]
Cell.type.boxplot
x.m<-melt(Cell.type.boxplot, id.var="Tissue")
x.m

pdf("dotplot_celltype_APM_TIS_Prolif.pdf", useDingbats = FALSE, height = 6, width = 16)
ggplot(x.m, aes(x = as.factor(Tissue), y = value, group = Tissue, fill=Tissue)) + ylim(0,8.5)+
  geom_boxplot(fill="white", color="grey", lwd = 0.2)+
  geom_dotplot(binaxis='y', stackdir='center', binpositions="all", dotsize = 1.0, binwidth = 0.5,color = "NA")+
  stat_compare_means(size=3, label.y = 8, label.x=2)+
  scale_x_discrete(limits=c(
    "diaphragm",
    "lung",
    "mediastinal_lymph_nodes",
    "pericardium",
    "chest wall"))+
  scale_fill_manual(values = c(
    'lung' = '#1F78B4',
    'diaphragm' = "#E31A1C",
    'mediastinal_lymph_nodes' = "#FF7F00",
    'pericardium' = "#33A02C",
    'chest wall' = "purple4")) +
  theme_classic()+  facet_rep_wrap(~factor(x.m$variable), scales = "fixed", repeat.tick.labels = TRUE, ncol= 6)+
  coord_capped_cart(bottom='both', left='both')+
  theme(text = element_text(family="Helvetica",size=13),legend.position="bottom", legend.text= element_text(size=11), legend.title=element_blank(), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), strip.text = element_text(size=11), axis.text.y=element_text( size = 10))
dev.off()