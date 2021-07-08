rm(list=ls())

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggthemes)
library(ggrepel)

TP53=read.table(file="ddPCR_TP53_longitudial_tumors.txt", header = TRUE, sep = "\t")

RESP= read.table(file="Response_duration.txt", header = TRUE, sep = "\t")
RESP$Start <- as.Date(RESP$Start, format = "%Y-%m-%d")
RESP$End <- as.Date(RESP$End, format = "%Y-%m-%d")

RESP1<-RESP

RESP1[1,2]<-"2012-11-15"
TP53
TP53<-TP53[-c(1,2,3,5,13,19),]
TP53$Date <- as.Date(TP53$Date, format = "%d-%b-%y") 
TP53

pdf("longitudinal_TP53ddPCR.pdf", useDingbats = FALSE, height = 4, width = 6)
ggplot(data = TP53, aes(x = Date, y = Mutation.Allele.Frequency))+ylim(0, 1)+
  scale_x_date(breaks=as.Date(c("2013-11-28","2015-01-28","2018-06-15")), date_labels = c("d373","d799","d2,033"))+
  geom_point(aes(y = Mutation.Allele.Frequency, shape = Tissue, color = Tissue))+
  stat_summary(fun = mean, geom = 'line', colour = "blue", aes(y = Mutation.Allele.Frequency, group = 1))+
  scale_shape_manual(values=c(16,17), label=c("blood","tumor"))+
  scale_color_manual(values=c("blue","blue"), guide=FALSE)+
  guides(shape=guide_legend(override.aes = list (color = "blue")))+
  labs(shape = 'TP53fs (disease burden)')+
  geom_rect(data = RESP1, aes(xmin = as.Date(Start,"%Y-%m-%d"), xmax = as.Date(End,"%Y-%m-%d"), fill = Response), 
            ymin = -Inf, ymax = Inf, alpha = 0.2,
            inherit.aes = FALSE, show.legend = FALSE)+
  scale_fill_manual(values = c("green", "blue", "grey", "red", "orange", "yellow" )) +
  labs(y = "Mutation Allele Frequency")+
  theme_classic()+
  theme(text = element_text(family="Helvetica"),legend.position="right",axis.title.x=element_blank(), 
        ,plot.margin=margin(t = 1, b = 1, unit = "pt"))
dev.off()