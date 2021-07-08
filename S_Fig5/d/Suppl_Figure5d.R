rm(list=ls())

library(ggplot2)
library(ggpubr)


IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE,row.names = 1)

IO_autopsy<-subset(IO, Response %in% "Autopsy")

#remove unused rows
IO_autopsy<-IO_autopsy[-c(5),]

IO_autopsy$Immunophenotype<-IO_autopsy$APM+IO_autopsy$IFN.Gamma

pdf("Immuno_mut_synchronous.pdf", useDingbats = FALSE, height = 4, width = 7)
ggplot(IO_autopsy, aes(x = Immunophenotype , y = TMB.norm))+
  geom_point(aes(color = Tissue)) +
  geom_smooth(method = "lm")+
  scale_colour_manual(values = c(
    'lung' = '#1F78B4',
    'diaphragm' = "#E31A1C",
    'mediastinal_lymph_nodes' = "#FF7F00",
    'pericardium' = "#33A02C",
    'chest wall' = "purple4"))+
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = rownames(IO_autopsy)),size = 3)+
  labs(x="Immunophenotype",y="Mutations per Mb (TMB)",color = "Parallel")+
  theme(text = element_text(family="Helvetica", size=10), axis.title = element_text(size = 10), 
        axis.text.x = element_text(size=9), axis.text.y = element_text(size=9), 
        legend.position = "bottom", panel.grid = element_blank(),plot.margin=margin(t = 1, b = 1, unit = "pt"))
dev.off()