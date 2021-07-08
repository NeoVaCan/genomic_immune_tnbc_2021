
rm(list=ls())

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(ggrepel)

IO<-read.delim(file="IO_signatures.txt", header = TRUE, sep = "\t", row.names = 1)

fit<-lm(IFN.Gamma~APM, data=IO)


summary(fit)


pdf("Immunophenotype.pdf", useDingbats = FALSE, height = 4, width = 7)
ggplot(IO, aes(x = IFN.Gamma , y = APM))+
  geom_point(aes(color = Response)) +
  geom_smooth(method = "lm")+
  geom_hline(yintercept = 4.25, colour = 'black', linetype='dotted')+
  geom_vline(xintercept = 3.70, colour = 'black', linetype='dotted')+
  scale_color_manual(values = c("black","green", "red", "orange", "yellow")) +
  theme_bw(base_size = 12)+
  geom_text_repel(aes(label = rownames(IO)), size = 3)+
  labs(x=expression(IFN~gamma~ signaling),y="APM", color = "Best response")+
  theme(text = element_text(family="Helvetica", size=10), axis.title = element_text(size = 10),
        axis.text.x = element_text(size=9), axis.text.y = element_text(size=9), legend.position = "top", 
        panel.grid = element_blank(),plot.margin=margin(t = 1, b = 1, unit = "pt"))+
  annotate("text", x = 6.2, y = 3.6, label = "paste(italic(p) < 2.2e-16)", colour="black", size = 3, parse=TRUE)
dev.off()