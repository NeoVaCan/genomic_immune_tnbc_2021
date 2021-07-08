rm(list=ls())

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(RColorBrewer)


setwd("~/Desktop/Figures_scripts/Figure4a/")
IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)
IO$Date <- as.Date(IO$Date, format = "%d-%b-%y")

#Remove unused rows
IO <- IO[-c(18),]


Timing <-data.frame(Timing = rep(c("Sequential", "Parallel"), c(13,19)))

# Extract and format columns for the APM heatmap
IO2<-cbind(IO[,c("Tissue","TN", "APM", "MHC2")], Timing)

Genes=read.table(file="IO_adjusted_log2_expression_data.txt", header = TRUE, sep = "\t",row.names =1 )

# Select genes related to the APM signature
Genes.APM<-Genes[,c("B2M", "NLRC5", "PSMB9", "PSMB10", "TAP1","TAP2", "STAT1", "STAT2","HLA.A","HLA.B","HLA.C", "HLA.E","HLA.F","HLA.DQA1","HLA.DQA2","HLA.DRB5" ,"HLA.DRB1", "HLA.DOA","HLA.DOB",  "HLA.DMA","HLA.DMB", "HLA.DPB1", "HLA.DQB1")]

HLA.IO<-merge(IO2,Genes.APM,by=0, all=TRUE)
rownames(HLA.IO)<- HLA.IO[,1]
HLA.IO <- HLA.IO[,-c(1)]

scores.scaled <- as.data.frame(apply(HLA.IO[,c(3,4,6:28)], 2, scale))
rownames(scores.scaled)<-rownames(HLA.IO)
my_sample_col <- as.data.frame(HLA.IO[,"Timing"])
annotation_row<-data.frame(my_sample_col,HLA.IO[, 1:2])
colnames(annotation_row)<-c("Timing", "Parallel", "TN")

rownames(annotation_row)<-rownames(scores.scaled)

#Color for Response and TILs annotation
anno_col = list(
  Timing = c(Sequential = "darkorchid4", Parallel = "darkgrey"),
  Parallel = c(
    'lung' = '#1F78B4',
    'diaphragm' = "#E31A1C",
    'mediastinal_lymph_nodes' = "#FF7F00",
    'pericardium' = "#33A02C",
    'chest wall' = "purple4"),
  TN = c(BLIA = "orange", BLIS = "deepskyblue2" , MES = "gold", "NA" = "white"))

pdf("heatmap_APM-HLA-flip.pdf", useDingbats = FALSE, height = 10, width =9)
scores.scaled<-na.omit(scores.scaled)
pheatmap(as.matrix(t(scores.scaled)),  row_order = FALSE, fontsize_col = 8,  fontsize_row = 8, breaks = seq(min(t(scores.scaled)),abs(min(t(scores.scaled))),length.out = 101),
         legend_breaks = c(-3,-2,-1,0,1,2,3), annotation_col = annotation_row, annotation_colors = anno_col, 
         reverse_colorscale = FALSE, clustering_distance_cols = "canberra", cellheight = 11, 
         color = colorRampPalette(rev(brewer.pal(n = 12, name ="RdYlBu")))(100))
dev.off()