rm(list=ls())

library(corrplot)
library(Hmisc)
library(ggplot2)
library(pheatmap)


IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)

#Remove unused rows
IO <- IO[-c(18),]


Genes=read.table(file="IO_adjusted_log2_expression_data.txt", header = TRUE, sep = "\t",row.names =1 )

# Select genes related to the APM signature
Genes.APM<-Genes[,c("B2M", "NLRC5", "PSMB9", "PSMB10", "TAP1","TAP2", "STAT1", "STAT2","HLA.A","HLA.B","HLA.C", "HLA.E","HLA.F","HLA.DQA1","HLA.DQA2","HLA.DRB5" ,"HLA.DRB1", "HLA.DOA","HLA.DOB",  "HLA.DMA","HLA.DMB", "HLA.DPB1", "HLA.DQB1")]

HLA.IO<-merge(IO[,c("APM", "MHC2")],Genes.APM,by=0, all=TRUE)
rownames(HLA.IO)<- HLA.IO[,1]
HLA.IO <- HLA.IO[,-c(1)]

HLAimbalance=read.table(file="HLAimbalance_summary_numbers.txt", header = TRUE, sep = "\t", row.names = 1)

# Merge HLAs gene expression with HLAimbalance
HLAimbalance.IO<-merge(HLA.IO,HLAimbalance, by=0, all.y=TRUE)
rownames(HLAimbalance.IO)<-HLAimbalance.IO[,1]

# Remove unused rows and columns
HLAimbalance.IO<-HLAimbalance.IO[-c(17),-c(1,27:32)]

HLAimbalance.IO.rcorr = rcorr(as.matrix(HLAimbalance.IO),  type = "spearman")
HLAimbalance.IO.coeff = HLAimbalance.IO.rcorr$r
HLAimbalance.IO.p = HLAimbalance.IO.rcorr$P

pdf("HLAimbalance_APM_correlations_symm_matrix_.pdf", useDingbats = FALSE, height = 8, width = 8)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0", "#FFFFFF","#FDDBC7", "#F4A582","#D6604D","#B2182B", "#67001F"))
corrplot::corrplot(HLAimbalance.IO.rcorr$r, sig.level = 0.05, insig = "blank",  p.mat = HLAimbalance.IO.rcorr$P,
         type = "upper", tl.cex = 1.0, col = col2 (50))
dev.off()