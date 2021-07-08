rm(list=ls())

library(ggplot2)
library(ggpubr)
library(dplyr)


setwd("~/Desktop/Figures_scripts/Figure4b")

HLAimbalance<-read.delim(file="HLAimbalance_summary.txt", header = TRUE, sep = "\t",)

pdf("lohla_imbalance_M1.pdf",useDingbats = FALSE,  height = 4, width = 10)
ggplot(HLAimbalance, aes(x=Sample_ID, y = Allele, fill=HLAimbalance, na.rm=TRUE))+
  geom_raster(aes (x=Sample_ID))+
  geom_tile(colour="black")+
  #coord_equal(ratio = 1)+
  coord_fixed(ratio=0.3)+
  labs(x = "", y = "")+    
  scale_y_discrete(expand = c(0, 0),limits=c("HLA-A*11:01:01:01","HLA-A*34:02:01","HLA-B*27:05:02:01", "HLA-B*08:01:01:01", "HLA-C*02:02:02:01", "HLA-C*07:01:01:01" )) +
  scale_fill_manual(values = c('green', 'red'))+
  scale_x_discrete(limits=c("M1", "M15", "M16", "M17", "M18", "M5A", "M5B", "M5C", "M6A", "M6B", "M7", "M8A", "M8B", "M9", "M10C", "M10D", "M10E", "M11A", "M11B", "M12A","M12B"),
                   expand = c(0, 0))+
  labs(fill = "HLA imbalance\n (P-value < 0.01)")+
  theme_void()+
  theme(text = element_text(family = "Helvetica", size = 8),axis.title.y=element_blank(),
        axis.title.x=element_blank(), axis.text.y=element_text( size = 8, face = "plain", 
        color = "black"), axis.text.x = element_text( size = 8, face = "plain", color = "black"), 
        legend.text=element_text(size = 8), legend.title=element_text( size = 8),
        axis.ticks = element_blank(),plot.margin=margin(t = 1, b= 1, unit = "pt"))+
  rotate_x_text(angle = 45)
dev.off()