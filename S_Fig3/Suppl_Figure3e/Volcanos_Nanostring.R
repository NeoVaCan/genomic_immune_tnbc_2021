# Volcano plots Nanostring figure 5

library(ggrepel)
library(ggpubr)

genes <- read.table("DE_GroupingResults.csv", header = TRUE, sep = ",", row.names = 1)
genes
genes$Genes<-rownames(genes)
genes$Significant <- ifelse(genes$FDR < 0.05, "FDR < 0.05", "NS")
a<-ggplot(genes[1:46,], aes(x = logFC, y = -log10(PValue))) +ggtitle(label="Pre-chemo vs Autopsy")+
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "black")) +
  geom_point(data = genes[33:771,],aes(x = logFC, y = -log10(PValue), color = Significant)) +
  scale_color_manual(values = c("red", "grey"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom",axis.title.x=element_blank()) +
  geom_text_repel(
    data = subset(genes[1:46,], FDR < 0.05),
    aes(label = Genes),
    size = 3,
    max.overlaps = 40,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

v1<-a + geom_text_repel(
  data = subset(genes[46:771,], FDR < 0.01),
  aes(label = Genes),
  size = 3,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines")
)

genes <- read.table("Nanostring-scores.matrix.lung_vs_chest_wall.edgeR.DE_results", header = TRUE, sep = "\t", row.names = 1)
genes
genes$Genes<-rownames(genes)
genes$Significant <- ifelse(genes$FDR < 0.05, "FDR < 0.05", "NS")
b<-ggplot(genes[1:46,], aes(x = logFC, y = -log10(PValue))) +ggtitle(label="Lung vs Chest wall")+
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("red", "black")) +
  geom_point(data = genes[46:810,],aes(x = logFC, y = -log10(PValue), color = Significant)) +
  scale_color_manual(values = c("red", "grey"))+
  theme_bw(base_size = 12) + theme(legend.position = "bottom",axis.title.x=element_blank()) +
  geom_text_repel(
    data = subset(genes[1:46,], FDR < 0.05),
    aes(label = Genes),
    size = 3,
    max.overlaps = 40,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )

v2<-b + geom_text_repel(
  data = subset(genes[33:810,], FDR < 0.01),
  aes(label = Genes),
  size = 3,
  box.padding = unit(0.35, "lines"),
  point.padding = unit(0.3, "lines")
)

figure <- ggarrange(v1,v2,
                    ncol = 1, nrow = 2, align = "hv", common.legend = TRUE, legend = "bottom") 
pdf("Volcanos-nanostring.pdf", useDingbats = FALSE, height = 9, width = 9)
figure
dev.off()
