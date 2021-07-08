rm(list=ls())

library(corrplot)
library(Hmisc)
library(ggplot2)
library(pheatmap)


setwd("~/Desktop/Figures_scripts/Suppl_Figure3c/")

IO=read.table(file="IO_signature_scores_TILs.txt", header = TRUE, sep = "\t", row.names = 3)
IO$Date <- as.Date(IO$Date, format = "%d-%b-%y") 

TILs<-IO[,c("TILs_Counts", "CD8.T.Cells", "Cytotoxic.Cells", "T.Cells", "Th1.Cells","Exhausted.CD8")]

TILs_T.Cells<-na.omit(TILs)

TILs_T.Cells_lm<-lm(TILs_Counts~T.Cells, data=IO)
TILs_T.Cells_lm

summary(TILs_T.Cells_lm)


TILs_T.Cells.rcorr = rcorr(as.matrix(TILs_T.Cells), type = "spearman")
TILs_T.Cells.coeff = TILs_T.Cells.rcorr$r
TILs_T.Cells.p = TILs_T.Cells.rcorr$P

col2 <- colorRampPalette(c("yellow","red", "purple4"))

pdf("CD4.CD8_Cell_symm_matrix.pdf", useDingbats = FALSE, height = 7, width = 8)
corrplot::corrplot(TILs_T.Cells.rcorr$r, sig.level = 0.05, insig = "blank", p.mat = TILs_T.Cells.rcorr$P, type = "upper", tl.cex = 0.8, col = col2 (50))
dev.off()


suppl.data<-read.delim("Suppl.Table2.txt")
cd4.cd8<- subset(suppl.data[, c(1:4,18,19,21:26)], Specimen %in% c("FFPE","Frozen_tumor/FFPE"))
cd4.cd8$CD4<-gsub("NE","0",cd4.cd8$CD4)
cd4.cd8$CD8<-gsub("NE","0",cd4.cd8$CD8)
cd4.cd8=type.convert(cd4.cd8)
cd4.cd8$Date<-as.Date(cd4.cd8$Date, format = "%d-%b-%y") 

cd4.cd8.correlations<-merge(cd4.cd8[,c(1,5,6)], TILs[,-c(1)], by.x=1, by.y=0, all=FALSE)

cd4.cd8.correlations<-na.omit(cd4.cd8.correlations[,-c(1)])

cd4.cd8.correlations.rcorr = rcorr(as.matrix(cd4.cd8.correlations), type = "spearman")
cd4.cd8.correlations.coeff = cd4.cd8.correlations.rcorr$r
cd4.cd8.correlations.p = cd4.cd8.correlations.rcorr$P

pdf("TILs_Cell_symm_matrix.pdf", useDingbats = FALSE, height = 7, width = 8)
corrplot::corrplot(cd4.cd8.correlations.rcorr$r, sig.level = 0.05, insig = "blank", p.mat = cd4.cd8.correlations$P, type = "upper", tl.cex = 0.8, col = col2 (50))
dev.off()



