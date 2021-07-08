rm(list=ls())

library(corrplot)
library(Hmisc)
library(ggplot2)
library(ggpubr)
library(pheatmap)


IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)

IO.purity<-subset(IO, Days %in% c("799","2033"))

# Remove unused rows and columns
IO.purity<-IO.purity[-c(1,75),7:54]


IO.purity.rcorr = rcorr(as.matrix(IO.purity),  type = "spearman")
IO.purity.coeff = IO.purity.rcorr$r
IO.purity.p = IO.purity.rcorr$P

pdf("Purity_TMB_GE_correlations_symm_matrix_.pdf", useDingbats = FALSE, height = 8, width = 8)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0", "#FFFFFF","#FDDBC7", "#F4A582","#D6604D","#B2182B", "#67001F"))
corrplot::corrplot(IO.purity.rcorr$r, sig.level = 0.05, insig = "blank", hcluster.method = "ward", order = "hclust", 
         p.mat = IO.purity.rcorr$P, type = "upper", tl.cex = 0.8, col = col2 (50))
dev.off()