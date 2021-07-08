
rm(list=ls())

library(ggplot2)
library(pheatmap)

IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)
IO$Date <- as.Date(IO$Date, format = "%d-%b-%y")
IO
heatmap_sig<-IO[, c('Date','Tissue', 'Response', 'Prolif',
                    'Stroma',
                    'Lymphoid',
                    'Myeloid',
                    'Endothelial.Cells',
                    'APM',
                    'MHC2',
                    'IFN.Gamma',
                    'Cytotoxicity',
                    'Immunoproteasome',
                    'Apoptosis',
                    'Inflammatory.Chemokines',
                    'Hypoxia',
                    'MAGEs',
                    'Glycolytic.Activity',
                    'IFN.Downstream',
                    'Myeloid.Inflam',
                    'B.Cells',
                    'CD45',
                    'CD8.T.Cells',
                    'Cytotoxic.Cells',
                    'DC',
                    'Exhausted.CD8',
                    'Macrophages',
                    'Mast.Cells',
                    'Neutrophils',
                    'NK.CD56dim.Cells',
                    'NK.Cells',
                    'T.Cells',
                    'Th1.Cells',
                    'Treg',
                    'TIS',
                    'ARG1',
                    'NOS2',
                    'IDO1',
                    'PDL1',
                    'CTLA4',
                    'IL10',
                    'PDL2',
                    'B7.H3',
                    'TIGIT',
                    'TGF.Beta',
                    'PD1',
                    'MMR.loss',
                    'APM.loss',
                    'JAKSTAT.loss')]

heatmap_sig$Date <- format(as.Date(heatmap_sig$Date), "%d-%b-%y")
heatmap_sig
heatmap_sig[heatmap_sig < 0] <- 0

#remove unused rows
heatmap_sig <- heatmap_sig[-c(18),]

scores.scaled <- as.data.frame(apply(heatmap_sig[,4:49], 2, scale))
rownames(scores.scaled)<-rownames(heatmap_sig)
my_sample_col <- data.frame(Timing = rep(c("Sequential", "Parallel"), c(13,19)))
annotation_row<-data.frame(my_sample_col,heatmap_sig[,2])
colnames(annotation_row)<-c("Timing", "Parallel metastases" )
rownames(annotation_row)<-rownames(scores.scaled)
scores.scaled
annotation_row

anno_col = list(
  Timing = c(Sequential = "darkorchid4", Parallel = "darkgrey"),
  "Parallel metastases" = c(
    'lung' = '#1F78B4',
    'diaphragm' = "#E31A1C",
    'mediastinal_lymph_nodes' = "#FF7F00",
    'pericardium' = "#33A02C",
    'chest wall' = "purple4"))

pdf("heatmap_IOsignatures_unsupervised.pdf", useDingbats = FALSE, height = 8, width =11)
scores.scaled<-na.omit(scores.scaled)
pheatmap(as.matrix(scores.scaled),  row_order = FALSE, fontsize_col = 8,  fontsize_row = 8, 
         breaks = seq(min(t(scores.scaled)),abs(min(t(scores.scaled))),length.out = 101),
         legend_breaks = c(-3,-2,-1,0,1,2,3), annotation_row = annotation_row, 
         annotation_colors = anno_col, reverse_colorscale = FALSE, clustering_distance_rows = "canberra", 
         cellheight = 11, color = colorRampPalette(rev(brewer.pal(n = 12, name ="RdYlBu")))(100))
dev.off()

```
