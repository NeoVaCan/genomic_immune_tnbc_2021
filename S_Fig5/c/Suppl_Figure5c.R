rm(list=ls())

library(ggplot2)
library(ggpubr)


IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE,row.names = 1)
IO_TMB<-subset(IO, Days %in% c("799","2033"))

#remove unused rows
IO_TMB<-IO_TMB[-c(1),]
IO_TMB[1,2] <-"M1"
tcga_brca<-read.table(file="brca_tcga_merged.txt", header = TRUE, sep="\t", na.strings = "NA")



#get how many tumor from each tumor type are in the pan dataset
table(tcga_brca$Type)

#Normalize mutation counts
tcga_brca$Mutation.Count.Normalized<-tcga_brca$Mutation.Count/40
means<-aggregate(tcga_brca$Mutation.Count.Normalized, by=list(tcga_brca$Type),mean, paired=FALSE)
median<-aggregate(tcga_brca$Mutation.Count.Normalized, by=list(tcga_brca$Type),median, paired=FALSE)
means
median

pdf("TMBnormalized_TCGA_Delta.pdf", useDingbats = FALSE,height = 4, width = 5)
ggplot(data = tcga_brca, mapping = aes(x = reorder(Type, Mutation.Count.Normalized,FUN = "median"), y = log10(Mutation.Count.Normalized)))+
  geom_hline(yintercept = 0.185635, colour = 'black', linetype='dotted')+
  geom_hline(yintercept = 0.4988405, colour = 'black', linetype='dotted')+
  geom_point(color = "grey")+ 
  geom_jitter(data = IO_TMB, mapping = aes(x = reorder(Tissue, TMB.norm,FUN = "median"), y = log10(TMB.norm),color = Tissue), width = 0.2, height = 0.1)+
  stat_summary(
    geom = "crossbar",
    fun = "median",
    col = "red",
    width = 0.4,
    lwd = 0.3,
  )+
  scale_x_discrete(limits=c(
    "BRCA","TNBC", "M1", "chest wall", "diaphragm","lung","mediastinal_lymph_nodes", "pericardium"), labels=c(
      "TCGA BRCA \n (primary)","TCGA TNBC \n (primary)", "M1", "chest wall", "diaphragm","lung","lymph nodes", "pericardium"))+
  scale_colour_manual(values = c(
    'BRCA' = 'grey',
    'TNBC' = 'grey',
    'lung' = '#1F78B4',
    'diaphragm' = "#E31A1C",
    'mediastinal_lymph_nodes' = "#FF7F00",
    'pericardium' = "#33A02C",
    'chest wall' = "purple4",
    'M1' = "black"))+
  scale_y_continuous(name = "Mutations per Mb", labels =c("0", "0.1", "1", "10", "100"))+ labs(x = "")+
  theme_classic()+
  theme(legend.position="none")+rotate_x_text(angle = 45)
dev.off()
