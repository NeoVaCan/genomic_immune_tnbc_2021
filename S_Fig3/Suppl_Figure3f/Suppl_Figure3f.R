rm(list=ls())

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(ggrepel)

IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)


# Plot MHC II and Th1
Th1.Cells.median<-median(IO$Th1.Cells)
MHC2.median<-median(IO$MHC2)
Th1.Cells.median
MHC2.median
fit1<-lm(Th1.Cells~MHC2, data=IO)
summary(fit1)


a<-ggplot(IO, aes(x = Th1.Cells , y = MHC2))+
  geom_point(aes(color = Response)) +
  geom_smooth(method = "lm")+
  geom_hline(yintercept = 3.03, colour = 'black', linetype='dotted')+
  geom_vline(xintercept = 4.60, colour = 'black', linetype='dotted')+
  scale_color_manual(values = c("black","green", "red", "orange", "yellow")) +
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = rownames(IO)),size = 3)+
  labs(x="Th1 Cells",y="MHC Class II")+
  theme(text = element_text(family="Helvetica", size=10),axis.title = element_text(size = 10), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9), legend.position = "top", panel.grid = element_blank(),plot.margin=margin(t = 1, b = 1, unit = "pt"))


# Plot MHC I and CD8 T Cells
CD8.T.Cells.median<-median(IO$CD8.T.Cells)
APM.median<-median(IO$APM)
CD8.T.Cells.median
APM.median
fit2<-lm(CD8.T.Cells~APM, data=IO)
summary(fit2)


b<-ggplot(IO, aes(x = CD8.T.Cells , y = APM))+
  geom_point(aes(color = Response)) +
  geom_smooth(method = "lm")+
  geom_hline(yintercept = 4.26, colour = 'black', linetype='dotted')+
  geom_vline(xintercept = 4.21, colour = 'black', linetype='dotted')+
  scale_color_manual(values = c("black","green", "red", "orange", "yellow")) +
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = rownames(IO)),size = 3)+
  labs(x="CD8 T Cells",y="APM")+
  theme(text = element_text(family="Helvetica", size=10),axis.title = element_text(size = 10), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9), legend.position = "top", panel.grid = element_blank(),plot.margin=margin(t = 1, b = 1, unit = "pt"))

figure <- ggarrange(a,b,
                    ncol = 1, nrow = 2, align = "hv", common.legend = TRUE, legend = "bottom")
pdf("Immunophenotypes_APM_MHC2_CD8_th1.pdf", useDingbats = FALSE, height = 8, width = 9)
figure
dev.off()